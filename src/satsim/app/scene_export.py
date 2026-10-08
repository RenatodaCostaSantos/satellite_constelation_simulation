"""Exportação da "cena" para o viewer: Python calcula, o viewer só desenha (ADR 0001, GUI).

``build_scene`` propaga cada satélite numa grade uniforme de tempo e devolve um ``dict`` no
formato do contrato JSON do ADR: trajetórias em ECI (km) e o subponto (lat, lon, alt) em
graus/km, estações de solo e passagens (vazio até a Semana 3).

Subponto na Semana 1 (aproximação esférica): r_ECEF = Rz(θG)·r_ECI, com
θG(t) = θG(época) + (dθG/dt)·t (IAU-82, ver ``build_scene``), lat = arcsin(z/r),
lon = arctan2(y, x) em [−180°, 180°) e alt = r − R⊕. O geodésico WGS84 entra na S2-03.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from typing import Any

import numpy as np
from numpy.typing import NDArray

from satsim.astrotime import datetime_to_jd, gmst, gmst_rate
from satsim.constants import R_EARTH
from satsim.propagators import Propagator

SCHEMA_VERSION = 1
"""Versão do contrato da cena (``meta.schema_version``)."""

ECI_DECIMALS = 6
"""Casas decimais de ``eci_km`` (1e-6 km = 1 mm). Dá folga ao teste de consistência
lla ↔ eci a 1e-6°: o erro de longitude é amplificado por 1/cos(lat) perto de 82°; com 1 cm o
pior caso fica em 7,5e-7°, com 1 mm em 5,1e-7° (custo: +13 kB no JSON)."""

ANGLE_DECIMALS = 6
"""Casas decimais de lat/lon [°] (1e-6° ≈ 0,1 m na superfície)."""

ALT_DECIMALS = 4
"""Casas decimais de ``alt_km`` (1e-4 km = 0,1 m)."""

_FRAME = "ECI simplificado (x -> equinócio vernal)"


@dataclass(frozen=True)
class SceneSatellite:
    """Satélite da cena: dados de exibição + propagador.

    Provisório: será substituído pela classe ``Satellite`` do ADR na S5-01.
    """

    id: str
    name: str
    color: str  # cor CSS, ex.: "#00e5ff"
    propagator: Propagator


@dataclass(frozen=True)
class SceneStation:
    """Estação de solo da cena (graus e metros).

    Provisório: será substituído por ``GroundStation`` do ADR na Semana 4.
    """

    id: str
    name: str
    lat_deg: float
    lon_deg: float
    alt_m: float


def eci_to_ecef(r_eci: NDArray[np.float64], theta: NDArray[np.float64]) -> NDArray[np.float64]:
    """Rotação ECI → ECEF por θG em torno de z (rotação do referencial, ADR 0001).

    x' = cos θ·x + sin θ·y, y' = −sin θ·x + cos θ·y, z' = z. ``r_eci`` (N, 3), ``theta`` (N,).
    """
    c, s = np.cos(theta), np.sin(theta)
    x, y, z = r_eci[:, 0], r_eci[:, 1], r_eci[:, 2]
    return np.stack([c * x + s * y, -s * x + c * y, z], axis=-1)


def subpoint_spherical(r_ecef: NDArray[np.float64]) -> NDArray[np.float64]:
    """Subponto esférico: (lat [°], lon [°] em [−180, 180), alt [km]) com shape (N, 3)."""
    r_norm = np.linalg.norm(r_ecef, axis=-1)
    lat = np.degrees(np.arcsin(r_ecef[:, 2] / r_norm))
    lon = np.degrees(np.arctan2(r_ecef[:, 1], r_ecef[:, 0]))
    lon = np.mod(lon + 180.0, 360.0) - 180.0
    alt_km = (r_norm - R_EARTH) / 1e3
    return np.stack([lat, lon, alt_km], axis=-1)


def build_scene(
    satellites: Sequence[SceneSatellite],
    stations: Sequence[SceneStation],
    epoch: datetime,
    duration_s: float,
    step_s: float,
) -> dict[str, Any]:
    """Monta a cena (contrato JSON do ADR 0001) amostrando todos os satélites no tempo.

    Args:
        satellites: satélites a exportar (qualquer quantidade, R7).
        stations: estações de solo.
        epoch: instante inicial da cena (tz-aware). Cada propagador pode ter época própria; o
            tempo é convertido para cada um.
        duration_s: duração da janela [s]; deve ser múltiplo de ``step_s``.
        step_s: passo da amostragem [s].

    Returns:
        ``dict`` serializável com ``meta``, ``earth``, ``satellites``, ``stations`` e
        ``passes``; cada trajetória tem N = duration_s/step_s + 1 amostras.

    Raises:
        ValueError: época naive, passo não positivo ou duração que não é múltiplo do passo.
    """
    if epoch.tzinfo is None or epoch.utcoffset() is None:
        raise ValueError("epoch deve ser um datetime tz-aware")
    if step_s <= 0 or duration_s < 0:
        raise ValueError("step_s deve ser positivo e duration_s não negativo")
    n_steps = round(duration_s / step_s)
    if not math.isclose(n_steps * step_s, duration_s, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("duration_s deve ser múltiplo de step_s")

    t = np.arange(n_steps + 1) * float(step_s)
    epoch_jd = datetime_to_jd(epoch)
    # θG(t) = θG(época) + (dθG/dt)·t, ambos da fórmula IAU-82. Numa janela de horas o termo
    # quadrático é desprezível (~1e-18 rad/s na taxa), e evita-se o ruído de ≈ 40 µs do JD em
    # float a cada amostra; é exatamente o modelo linear que o contrato entrega ao viewer.
    gmst0 = float(gmst(epoch_jd))
    rate = float(gmst_rate(epoch_jd))
    theta = gmst0 + rate * t

    sats_out = []
    for sat in satellites:
        offset_s = (epoch - sat.propagator.epoch).total_seconds()
        r_eci, _ = sat.propagator.propagate(t + offset_s)
        lla = subpoint_spherical(eci_to_ecef(r_eci, theta))
        lla_rounded = np.column_stack(
            [
                np.round(lla[:, 0], ANGLE_DECIMALS),
                np.round(lla[:, 1], ANGLE_DECIMALS),
                np.round(lla[:, 2], ALT_DECIMALS),
            ]
        )
        sats_out.append(
            {
                "id": sat.id,
                "name": sat.name,
                "color": sat.color,
                "eci_km": np.round(r_eci / 1e3, ECI_DECIMALS).tolist(),
                "lla": lla_rounded.tolist(),
            }
        )

    return {
        "meta": {
            "schema_version": SCHEMA_VERSION,
            "epoch_utc": epoch.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "duration_s": float(duration_s),
            "step_s": float(step_s),
            "frame": _FRAME,
        },
        "earth": {
            "gmst0_rad": gmst0,
            "gmst_rate_rad_s": rate,
        },
        "satellites": sats_out,
        "stations": [
            {
                "id": st.id,
                "name": st.name,
                "lat_deg": st.lat_deg,
                "lon_deg": st.lon_deg,
                "alt_m": st.alt_m,
            }
            for st in stations
        ],
        "passes": [],
    }
