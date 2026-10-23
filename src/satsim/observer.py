"""Estação de solo e coordenadas topocêntricas: ECEF → ENU → azimute, elevação e distância.

Para saber onde um observador vê o satélite, o vetor estação → satélite é escrito no referencial
local ENU (Leste, Norte, Zênite) da estação. As linhas da matriz ENU são os versores locais, dados
pela latitude geodésica φ e pela longitude λ da estação:

    E = (−sin λ,         cos λ,          0)
    N = (−sin φ·cos λ,  −sin φ·sin λ,   cos φ)
    U = ( cos φ·cos λ,   cos φ·sin λ,   sin φ)

    enu = M · (r_sat − r_estação)        (ambos em ECEF)

Daí, com ρ = |enu| (a distância):

    az = atan2(E, N) mod 2π    medido a partir do norte verdadeiro, no sentido horário
                               (0° = norte, 90° = leste, 180° = sul, 270° = oeste)
    el = arcsin(U / ρ)         em relação ao horizonte GEODÉSICO local (plano ⟂ à normal ao
                               elipsoide WGS84), de −90° a 90°; 90° = zênite

Sem refração: a atmosfera curva a luz e faz o satélite parecer mais alto, até ~0,5° perto do
horizonte e quase nada acima de ~10°. Para visibilidade por rádio, com elevação mínima de 10°, é
desprezível.

A estação entra como configuração: ``SAO_BENTO_DO_SAPUCAI`` abaixo é o único lugar com as
coordenadas (a S6-01 as moverá para um arquivo de configuração).
"""

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import R_EARTH
from satsim.frames import eci_to_ecef, gmst_at
from satsim.geodesy import geodetic_to_ecef
from satsim.propagators.base import Propagator, as_time_array

FloatOrArray = float | NDArray[np.float64]


def _as_scalar_if_0d(x: NDArray[np.float64]) -> FloatOrArray:
    return float(x) if x.ndim == 0 else x


@dataclass(frozen=True)
class GroundStation:
    """Estação de solo em coordenadas geodésicas WGS84.

    Args:
        name: nome para exibição.
        lat: latitude geodésica φ [rad].
        lon: longitude λ [rad] (positiva para leste).
        alt: altitude elipsoidal h [m].
    """

    name: str
    lat: float
    lon: float
    alt: float

    @classmethod
    def from_degrees(
        cls, name: str, lat_deg: float, lon_deg: float, alt_m: float
    ) -> "GroundStation":
        """Estação a partir de latitude e longitude em graus e altitude em metros."""
        return cls(name, math.radians(lat_deg), math.radians(lon_deg), float(alt_m))

    @property
    def lat_deg(self) -> float:
        return math.degrees(self.lat)

    @property
    def lon_deg(self) -> float:
        return math.degrees(self.lon)

    def ecef(self) -> NDArray[np.float64]:
        """Posição da estação em ECEF [m], shape (3,) (``geodesy.geodetic_to_ecef``)."""
        return geodetic_to_ecef(self.lat, self.lon, self.alt)

    def enu_matrix(self) -> NDArray[np.float64]:
        """Matriz (3, 3) ECEF → ENU: linhas E, N, U (versores locais; ver docstring do módulo).

        É ortonormal com determinante +1 (rotação): E × N = U.
        """
        sp, cp = math.sin(self.lat), math.cos(self.lat)
        sl, cl = math.sin(self.lon), math.cos(self.lon)
        return np.array(
            [
                [-sl, cl, 0.0],  # E: leste
                [-sp * cl, -sp * sl, cp],  # N: norte
                [cp * cl, cp * sl, sp],  # U: zênite (normal ao elipsoide)
            ]
        )

    def enu(self, r_ecef: ArrayLike) -> NDArray[np.float64]:
        """Vetor estação → alvo em ENU [m]: M·(r − r_estação).

        Args:
            r_ecef: posição do alvo em ECEF [m], shape (3,) ou (N, 3).

        Returns:
            (E, N, U) com o shape da entrada.
        """
        rel = np.asarray(r_ecef, dtype=float) - self.ecef()
        return rel @ self.enu_matrix().T

    def look_angles(self, r_ecef: ArrayLike) -> tuple[FloatOrArray, FloatOrArray, FloatOrArray]:
        """Azimute, elevação e distância do alvo vistos da estação, vetorizado.

        Args:
            r_ecef: posição do alvo em ECEF [m], shape (3,) ou (N, 3).

        Returns:
            (az [rad] em [0, 2π), el [rad] em [−π/2, π/2], distância [m]); escalares para entrada
            (3,), arrays (N,) para (N, 3).
        """
        enu = self.enu(r_ecef)
        e, n, u = enu[..., 0], enu[..., 1], enu[..., 2]
        rng = np.sqrt(e * e + n * n + u * u)
        az = np.mod(np.arctan2(e, n), 2.0 * math.pi)
        el = np.arcsin(np.clip(u / rng, -1.0, 1.0))
        return _as_scalar_if_0d(az), _as_scalar_if_0d(el), _as_scalar_if_0d(rng)


# Coordenadas APROXIMADAS, a confirmar pelo responsável do projeto (prompt da Semana 3, A.6):
# São Bento do Sapucaí (SP), ≈ −22,69°, −45,73°, ≈ 900 m.
SAO_BENTO_DO_SAPUCAI = GroundStation.from_degrees("São Bento do Sapucaí", -22.69, -45.73, 900.0)


@dataclass(frozen=True)
class Observation:
    """Satélite visto da estação: azimute [rad], elevação [rad] e distância [m], shape (N,)."""

    az: NDArray[np.float64]
    el: NDArray[np.float64]
    range: NDArray[np.float64]

    @property
    def az_deg(self) -> NDArray[np.float64]:
        return np.degrees(self.az)

    @property
    def el_deg(self) -> NDArray[np.float64]:
        return np.degrees(self.el)


def observe(propagator: Propagator, station: GroundStation, t: ArrayLike) -> Observation:
    """Azimute, elevação e distância do satélite vistos de ``station`` nos instantes t.

    Propaga (qualquer ``Propagator``), passa ao ECEF com o GMST de ``propagator.epoch + t``
    (``frames.eci_to_ecef``) e chama ``station.look_angles``.

    Args:
        propagator: Kepler, MeanJ2, Cowell, SGP4…
        station: estação de solo.
        t: segundos desde ``propagator.epoch``, escalar ou (N,).

    Returns:
        ``Observation`` com campos de shape (N,) (N = 1 para t escalar).
    """
    t = as_time_array(t)
    r_eci, _ = propagator.propagate(t)
    r_ecef = eci_to_ecef(r_eci, gmst_at(propagator.epoch, t))
    az, el, rng = station.look_angles(r_ecef)
    return Observation(np.atleast_1d(az), np.atleast_1d(el), np.atleast_1d(rng))


def footprint_half_angle(alt_m: ArrayLike, el_min: ArrayLike) -> FloatOrArray:
    """Raio angular geocêntrico do footprint λ = arccos(R⊕/(R⊕ + h)·cos ε) − ε (Sem4 §9.4).

    É o ângulo, no centro da Terra, entre o subponto do satélite e a borda da região de onde ele
    é visto acima da elevação ε (Terra esférica de raio R⊕). Ex.: h = 400 km, ε = 0 → 19,8°;
    h = 420 km, ε = 10° → ≈ 12,5°. Será usado no viewer (S3-06).

    Args:
        alt_m: altitude do satélite h [m], escalar ou array.
        el_min: elevação mínima ε [rad], compatível por broadcasting.

    Returns:
        λ [rad] com o shape da entrada (float para escalar).
    """
    h = np.asarray(alt_m, dtype=float)
    eps = np.asarray(el_min, dtype=float)
    lam = np.arccos(R_EARTH / (R_EARTH + h) * np.cos(eps)) - eps
    return _as_scalar_if_0d(lam)
