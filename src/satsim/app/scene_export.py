"""Exportação da "cena" para o viewer: Python calcula, o viewer só desenha (ADR 0001, GUI).

``build_scene`` propaga cada satélite numa grade uniforme de tempo e devolve um ``dict`` no
formato do contrato JSON do ADR: trajetórias em ECI (km) e o subponto (lat, lon, alt) em
graus/km, estações de solo e passagens (vazio até a Semana 3).

Subponto geodésico WGS84 (desde a S2-07): r_ECEF = Rz(θG)·r_ECI (``satsim.frames``), com
θG(t) = θG(época) + (dθG/dt)·t (IAU-82, ver ``build_scene``), convertido em latitude geodésica,
longitude em [−180°, 180°) e altitude elipsoidal (``satsim.geodesy``). Na Semana 1 o subponto
era esférico (lat = arcsin(z/r), alt = r − R⊕); a diferença chega a ~0,2° em latitude e ~21 km
em altitude. ``meta.geodesy = "WGS84"`` identifica cenas com o subponto geodésico.

Estação, visibilidade e passagens (S3-06, campos opcionais)
-----------------------------------------------------------
Com ``el_min_deg``, a cena ganha, sem mudar ``schema_version`` (todos os campos novos são
opcionais; o viewer funciona sem eles):

- ``satellites[i].footprint_deg[]``: raio angular do círculo de visibilidade a ``el_min``, por
  amostra (``observer.footprint_half_angle``, Sem4 §9.4);
- ``stations[j].looks[sat_id] = {az_deg[], el_deg[], range_km[]}``: o satélite visto da
  estação, por amostra (``GroundStation.look_angles`` sobre o mesmo r_ECEF do subponto);
- ``passes[]``: passagens acima de ``el_min`` na janela (``passes.find_passes``), com instantes
  em UTC ISO 8601;
- ``meta.el_min_deg`` e, se um TLE for dado, ``meta.tle = {name, norad_id, epoch_utc,
  age_days_at_start, source}``.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from typing import Any

import numpy as np
from numpy.typing import NDArray

from satsim.astrotime import datetime_to_jd, gmst, gmst_rate
from satsim.frames import eci_to_ecef
from satsim.geodesy import ecef_to_geodetic
from satsim.observer import GroundStation, footprint_half_angle
from satsim.passes import Pass, find_passes
from satsim.propagators import Propagator
from satsim.tle import Tle

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

LOOK_DECIMALS = 4
"""Casas decimais de az/el [°] em ``looks`` (1e-4° ≈ 1 m a 500 km de distância)."""

RANGE_DECIMALS = 3
"""Casas decimais de ``range_km`` em ``looks`` (1 m)."""

FOOTPRINT_DECIMALS = 6
"""Casas decimais de ``footprint_deg`` (1e-6°, ≈ 0,1 m na superfície)."""

_FRAME = "ECI simplificado (x -> equinócio vernal)"
_GEODESY = "WGS84"


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

    @classmethod
    def from_ground_station(cls, station_id: str, station: GroundStation) -> "SceneStation":
        """Estação da cena a partir de uma ``observer.GroundStation`` (S3-06)."""
        return cls(station_id, station.name, station.lat_deg, station.lon_deg, station.alt)

    def ground_station(self) -> GroundStation:
        """A ``GroundStation`` equivalente (para calcular az/el e passagens)."""
        return GroundStation.from_degrees(self.name, self.lat_deg, self.lon_deg, self.alt_m)


def _round_ms(dt: datetime) -> datetime:
    """``datetime`` em UTC arredondado ao milissegundo mais próximo."""
    dt = dt.astimezone(UTC)
    return (
        dt
        + timedelta(milliseconds=round(dt.microsecond / 1000))
        - timedelta(microseconds=dt.microsecond)
    )


def _iso_utc(dt: datetime) -> str:
    """``datetime`` → ISO 8601 em UTC, arredondado ao milissegundo: 2019-12-10T00:46:08.063Z."""
    return _round_ms(dt).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def _pass_json(p: Pass, station_id: str) -> dict[str, Any]:
    return {
        "sat_id": p.sat_id,
        "station_id": station_id,
        "aos_utc": _iso_utc(p.aos),
        "tca_utc": _iso_utc(p.tca),
        "los_utc": _iso_utc(p.los),
        "max_el_deg": round(p.max_elevation_deg, 3),
        "az_aos_deg": round(p.az_aos_deg, 2),
        "az_tca_deg": round(p.az_tca_deg, 2),
        "az_los_deg": round(p.az_los_deg, 2),
        "duration_s": round(p.duration_s, 1),
    }


def subpoint_geodetic(r_ecef: NDArray[np.float64]) -> NDArray[np.float64]:
    """Subponto geodésico WGS84: (lat [°], lon [°] em [−180, 180), alt [km]) com shape (N, 3)."""
    lat, lon, alt = ecef_to_geodetic(r_ecef)
    lon_deg = np.mod(np.degrees(lon) + 180.0, 360.0) - 180.0
    return np.stack([np.degrees(lat), lon_deg, alt / 1e3], axis=-1)


def build_scene(
    satellites: Sequence[SceneSatellite],
    stations: Sequence[SceneStation],
    epoch: datetime,
    duration_s: float,
    step_s: float,
    el_min_deg: float | None = None,
    tle: Tle | None = None,
    tle_source: str | None = None,
) -> dict[str, Any]:
    """Monta a cena (contrato JSON do ADR 0001) amostrando todos os satélites no tempo.

    Args:
        satellites: satélites a exportar (qualquer quantidade, R7).
        stations: estações de solo.
        epoch: instante inicial da cena (tz-aware). Cada propagador pode ter época própria; o
            tempo é convertido para cada um.
        duration_s: duração da janela [s]; deve ser múltiplo de ``step_s``.
        step_s: passo da amostragem [s].
        el_min_deg: elevação mínima [°] (S3-06). Se dada, a cena ganha footprint, ``looks`` de
            cada estação e as passagens; se ``None``, esses campos ficam de fora (cenas antigas).
        tle: TLE do satélite real, para ``meta.tle`` (época e idade no início da cena).
        tle_source: origem do TLE (``"network"``, ``"cache"``, ``"bundled"`` ou um arquivo).

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
    # A época da cena é arredondada ao milissegundo, a resolução de meta.epoch_utc: amostras,
    # passagens e o relógio do viewer partem exatamente do mesmo instante (ex.: época de um TLE).
    epoch = _round_ms(epoch)

    t = np.arange(n_steps + 1) * float(step_s)
    epoch_jd = datetime_to_jd(epoch)
    # θG(t) = θG(época) + (dθG/dt)·t, ambos da fórmula IAU-82. Numa janela de horas o termo
    # quadrático é desprezível (~1e-18 rad/s na taxa), e evita-se o ruído de ≈ 40 µs do JD em
    # float a cada amostra; é exatamente o modelo linear que o contrato entrega ao viewer.
    gmst0 = float(gmst(epoch_jd))
    rate = float(gmst_rate(epoch_jd))
    theta = gmst0 + rate * t

    ground = [st.ground_station() for st in stations]
    looks: list[dict[str, Any]] = [{} for _ in stations]
    passes_out: list[dict[str, Any]] = []
    sats_out = []
    for sat in satellites:
        # cada propagador tem a sua época; a cena tem a sua (S3-02)
        r_eci, _ = sat.propagator.propagate_utc(epoch, t)
        r_ecef = eci_to_ecef(r_eci, theta)
        lla = subpoint_geodetic(r_ecef)
        lla_rounded = np.column_stack(
            [
                np.round(lla[:, 0], ANGLE_DECIMALS),
                np.round(lla[:, 1], ANGLE_DECIMALS),
                np.round(lla[:, 2], ALT_DECIMALS),
            ]
        )
        sat_out = {
            "id": sat.id,
            "name": sat.name,
            "color": sat.color,
            "eci_km": np.round(r_eci / 1e3, ECI_DECIMALS).tolist(),
            "lla": lla_rounded.tolist(),
        }
        if el_min_deg is not None:
            el_min = math.radians(el_min_deg)
            # a partir da altitude publicada em lla: o JSON é verificável por quem o lê
            footprint = np.degrees(footprint_half_angle(lla_rounded[:, 2] * 1e3, el_min))
            sat_out["footprint_deg"] = np.round(footprint, FOOTPRINT_DECIMALS).tolist()
            for j, (st, gs) in enumerate(zip(stations, ground, strict=True)):
                az, el, rng = gs.look_angles(r_ecef)
                looks[j][sat.id] = {
                    "az_deg": np.round(np.degrees(az), LOOK_DECIMALS).tolist(),
                    "el_deg": np.round(np.degrees(el), LOOK_DECIMALS).tolist(),
                    "range_km": np.round(rng / 1e3, RANGE_DECIMALS).tolist(),
                }
                found = find_passes(sat.propagator, gs, epoch, duration_s, el_min_deg,
                                    step_s=min(20.0, float(step_s)), sat_id=sat.id)  # fmt: skip
                passes_out += [_pass_json(p, st.id) for p in found]
        sats_out.append(sat_out)
    passes_out.sort(key=lambda p: (p["aos_utc"], p["sat_id"], p["station_id"]))

    meta: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "epoch_utc": _iso_utc(epoch) if epoch.microsecond else f"{epoch:%Y-%m-%dT%H:%M:%SZ}",
        "duration_s": float(duration_s),
        "step_s": float(step_s),
        "frame": _FRAME,
        "geodesy": _GEODESY,
    }
    if el_min_deg is not None:
        meta["el_min_deg"] = float(el_min_deg)
    if tle is not None:
        meta["tle"] = {
            "name": tle.name,
            "norad_id": tle.norad_id,
            "epoch_utc": _iso_utc(tle.epoch),
            "age_days_at_start": round((epoch - tle.epoch).total_seconds() / 86400.0, 4),
            "source": tle_source or "",
        }

    stations_out = []
    for st, st_looks in zip(stations, looks, strict=True):
        st_out: dict[str, Any] = {
            "id": st.id,
            "name": st.name,
            "lat_deg": st.lat_deg,
            "lon_deg": st.lon_deg,
            "alt_m": st.alt_m,
        }
        if el_min_deg is not None:
            st_out["looks"] = st_looks
        stations_out.append(st_out)

    return {
        "meta": meta,
        "earth": {
            "gmst0_rad": gmst0,
            "gmst_rate_rad_s": rate,
        },
        "satellites": sats_out,
        "stations": stations_out,
        "passes": passes_out,
    }
