"""Subponto geodésico e traço no solo (ground track), base matemática do mapa 2D da GUI.

O subponto é o ponto da superfície diretamente abaixo do satélite: a posição ECI do propagador é
girada para ECEF pelo GMST (``frames``) e convertida para latitude, longitude e altitude
geodésicas WGS84 (``geodesy``). Como a Terra gira por baixo da órbita, o traço se desloca para
oeste a cada revolução: Δλ ≈ −(ω⊕ − Ω̇)·Tnod por volta (Sem4 §7.1).

Funciona com qualquer ``Propagator`` (Kepler, MeanJ2, Cowell), sem caso especial.
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.frames import eci_to_ecef, gmst_at
from satsim.geodesy import ecef_to_geodetic
from satsim.propagators.base import Propagator, as_time_array


def _wrap_pi(angle: NDArray[np.float64]) -> NDArray[np.float64]:
    """Reduz ângulos [rad] a [−π, π)."""
    return np.mod(angle + np.pi, 2.0 * np.pi) - np.pi


def subpoint(
    propagator: Propagator, t: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Subponto geodésico WGS84 do satélite.

    Args:
        propagator: qualquer ``Propagator``.
        t: segundos desde ``propagator.epoch``, escalar ou (N,).

    Returns:
        (lat, lon, alt): latitude geodésica [rad], longitude [rad] em (−π, π] e altitude
        elipsoidal [m], cada um com shape (N,) (N = 1 para t escalar).
    """
    t = as_time_array(t)
    r_eci, _ = propagator.propagate(t)
    r_ecef = eci_to_ecef(r_eci, gmst_at(propagator.epoch, t))
    lat, lon, alt = ecef_to_geodetic(r_ecef)
    return np.atleast_1d(lat), np.atleast_1d(lon), np.atleast_1d(alt)


def split_at_antimeridian(
    lon: ArrayLike, lat: ArrayLike
) -> list[tuple[NDArray[np.float64], NDArray[np.float64]]]:
    """Quebra o traço onde ele cruza o antimeridiano (±180°), para desenhar sem riscos no mapa.

    Um salto |Δλ| > π entre amostras consecutivas indica que o traço saiu por uma borda do mapa e
    entrou pela outra; ali começa um novo segmento. Nenhum ponto é criado nem descartado.

    Args:
        lon: longitudes [rad] em (−π, π], shape (N,).
        lat: latitudes [rad], shape (N,).

    Returns:
        Lista de segmentos (lon, lat), na ordem do tempo; dentro de cada um, |Δλ| ≤ π.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    breaks = np.flatnonzero(np.abs(np.diff(lon)) > np.pi) + 1
    return list(zip(np.split(lon, breaks), np.split(lat, breaks), strict=True))


def ascending_node_crossings(
    t: ArrayLike, lat: ArrayLike, lon: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Instantes e longitudes dos nós ascendentes (latitude cruzando 0 de sul para norte).

    Para cada par de amostras com lat[k] < 0 ≤ lat[k+1], interpola linearmente o instante em que
    lat = 0 e a longitude nesse instante (com a diferença de longitude reduzida a [−π, π), para
    não interpolar através do antimeridiano).

    Args:
        t: tempos [s], shape (N,), crescentes.
        lat: latitudes [rad], shape (N,).
        lon: longitudes [rad], shape (N,).

    Returns:
        (t_nodes, lon_nodes): instantes [s] e longitudes [rad] em [−π, π), shape (K,).
    """
    t = np.asarray(t, dtype=float)
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    k = np.flatnonzero((lat[:-1] < 0.0) & (lat[1:] >= 0.0))
    frac = -lat[k] / (lat[k + 1] - lat[k])
    t_nodes = t[k] + frac * (t[k + 1] - t[k])
    lon_nodes = _wrap_pi(lon[k] + frac * _wrap_pi(lon[k + 1] - lon[k]))
    return t_nodes, lon_nodes
