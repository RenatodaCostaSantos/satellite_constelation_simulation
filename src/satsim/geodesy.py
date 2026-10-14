"""Coordenadas geodésicas no elipsoide WGS84 ↔ ECEF, e a aproximação esférica.

Latitude **geodésica** φ é o ângulo entre o plano do equador e a *normal ao elipsoide* no ponto (é
a latitude dos mapas e do GPS); latitude **geocêntrica** é o ângulo entre o equador e a reta que
liga o ponto ao centro da Terra. As duas coincidem no equador e nos polos e diferem até ~0,19° a
45°. A altitude h é medida ao longo da normal, a partir da superfície do elipsoide.

Ida (fechada):  N = a/√(1 − e² sin²φ)
                x = (N + h) cos φ cos λ,  y = (N + h) cos φ sin λ,  z = (N(1 − e²) + h) sin φ
Volta (iterativa, Bowring 1976): latitude paramétrica β e
                φ = atan2(z + e'²·b·sin³β, p − e²·a·cos³β),  tan β = (1 − f) tan φ,
iterada até |Δφ| < 1e-12 rad. A altitude usa h = p·cos φ + z·sin φ − a·√(1 − e² sin²φ), que não
divide por cos φ nem por sin φ e vale nos polos, no equador e para h negativo.

Shapes: ECEF (3,) ↔ escalares; (N, 3) ↔ (N,). Ângulos em rad, distâncias em m.
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import R_EARTH, WGS84_A, WGS84_B, WGS84_E2, WGS84_F

FloatOrArray = float | NDArray[np.float64]

_EP2 = WGS84_E2 / (1.0 - WGS84_E2)  # e'² = (a² − b²)/b²


def _as_scalar_if_0d(x: NDArray[np.float64]) -> FloatOrArray:
    return float(x) if x.ndim == 0 else x


def geodetic_to_ecef(lat: ArrayLike, lon: ArrayLike, h: ArrayLike) -> NDArray[np.float64]:
    """Geodésico WGS84 (φ, λ, h) → ECEF [m].

    Args:
        lat: latitude geodésica φ [rad], escalar ou (N,).
        lon: longitude λ [rad], compatível por broadcasting.
        h: altitude elipsoidal [m], compatível por broadcasting.

    Returns:
        r_ecef [m], shape (3,) para entradas escalares ou (N, 3).
    """
    lat, lon, h = np.broadcast_arrays(*(np.asarray(x, dtype=float) for x in (lat, lon, h)))
    sin_lat, cos_lat = np.sin(lat), np.cos(lat)
    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * sin_lat**2)
    return np.stack(
        [
            (n + h) * cos_lat * np.cos(lon),
            (n + h) * cos_lat * np.sin(lon),
            (n * (1.0 - WGS84_E2) + h) * sin_lat,
        ],
        axis=-1,
    )


def ecef_to_geodetic(
    r_ecef: ArrayLike, tol: float = 1e-12, max_iter: int = 10
) -> tuple[FloatOrArray, FloatOrArray, FloatOrArray]:
    """ECEF → geodésico WGS84 (φ, λ, h), pelo método iterativo de Bowring.

    Args:
        r_ecef: posição ECEF [m], shape (3,) ou (N, 3). A origem (r = 0) não é suportada.
        tol: critério de parada sobre |Δφ| [rad].
        max_iter: número máximo de iterações (converge em 2–3 para pontos até GEO).

    Returns:
        (lat, lon, h): latitude geodésica [rad] em [−π/2, π/2], longitude [rad] em (−π, π]
        (0 sobre o eixo z) e altitude elipsoidal [m]; escalares para entrada (3,), (N,) para
        (N, 3).

    Raises:
        RuntimeError: se não convergir em ``max_iter`` iterações.
    """
    r = np.asarray(r_ecef, dtype=float)
    x, y, z = r[..., 0], r[..., 1], r[..., 2]
    p = np.hypot(x, y)
    lon = np.arctan2(y, x)
    # chute inicial: latitude paramétrica de um ponto na superfície
    beta = np.arctan2(z, (1.0 - WGS84_F) * p)
    lat = np.arctan2(z, p)
    for _ in range(max_iter):
        lat_new = np.arctan2(
            z + _EP2 * WGS84_B * np.sin(beta) ** 3, p - WGS84_E2 * WGS84_A * np.cos(beta) ** 3
        )
        beta = np.arctan2((1.0 - WGS84_F) * np.sin(lat_new), np.cos(lat_new))
        converged = np.all(np.abs(lat_new - lat) < tol)
        lat = lat_new
        if converged:
            break
    else:
        raise RuntimeError(f"ecef_to_geodetic não convergiu em {max_iter} iterações")
    sin_lat, cos_lat = np.sin(lat), np.cos(lat)
    h = p * cos_lat + z * sin_lat - WGS84_A * np.sqrt(1.0 - WGS84_E2 * sin_lat**2)
    return _as_scalar_if_0d(lat), _as_scalar_if_0d(lon), _as_scalar_if_0d(h)


def ecef_to_spherical(
    r_ecef: ArrayLike,
) -> tuple[FloatOrArray, FloatOrArray, FloatOrArray]:
    """ECEF → (latitude geocêntrica, longitude, r − R⊕): a aproximação esférica da Semana 1.

    Args:
        r_ecef: posição ECEF [m], shape (3,) ou (N, 3).

    Returns:
        (lat_gc, lon, alt): lat_gc = arcsin(z/r) [rad], lon [rad] em (−π, π], alt = r − R⊕ [m];
        escalares para entrada (3,), (N,) para (N, 3).
    """
    r = np.asarray(r_ecef, dtype=float)
    r_norm = np.linalg.norm(r, axis=-1)
    lat = np.arcsin(r[..., 2] / r_norm)
    lon = np.arctan2(r[..., 1], r[..., 0])
    return _as_scalar_if_0d(lat), _as_scalar_if_0d(lon), _as_scalar_if_0d(r_norm - R_EARTH)
