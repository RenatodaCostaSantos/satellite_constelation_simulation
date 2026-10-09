"""Rotação entre os referenciais ECI e ECEF pelo tempo sideral de Greenwich (GMST).

ECI aqui é o equatorial inercial **do equinócio da data, sem precessão nem nutação** (ADR 0001):
não é o GCRS/ICRF do Skyfield, e os dois diferem por ~0,37° em 2026 (precessão acumulada desde
J2000).
ECEF gira com a Terra; os dois compartilham o eixo z (sem movimento polar) e diferem por uma
rotação em torno de z pelo ângulo θ = GMST (UT1 ≈ UTC):

    r_ecef = Rz(θ)·r_eci,   Rz(θ) = [[cos θ, sin θ, 0], [−sin θ, cos θ, 0], [0, 0, 1]]

Rz(θ) é uma rotação *passiva* (dos eixos): a Terra girou θ para leste, então um ponto fixo no
espaço aparece θ mais a oeste no ECEF. Para velocidades, desconta-se o arrasto da rotação:

    v_ecef = Rz(θ)·(v_eci − ω⊕ × r_eci),   ω⊕ = (0, 0, OMEGA_EARTH)

Shapes: posições/velocidades (3,) ou (N, 3); θ escalar ou (N,) (um ângulo por linha).
"""

from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.astrotime import datetime_to_jd, gmst, jd_from_epoch
from satsim.constants import OMEGA_EARTH

FloatOrArray = float | NDArray[np.float64]


def rot_z(theta: ArrayLike) -> NDArray[np.float64]:
    """Matriz de rotação passiva Rz(θ) = [[c, s, 0], [−s, c, 0], [0, 0, 1]].

    Args:
        theta: ângulo [rad], escalar ou array (N,).

    Returns:
        Shape (3, 3) para θ escalar, (N, 3, 3) para θ (N,).
    """
    theta = np.asarray(theta, dtype=float)
    c, s = np.cos(theta), np.sin(theta)
    zero, one = np.zeros_like(theta), np.ones_like(theta)
    return np.stack(
        [
            np.stack([c, s, zero], axis=-1),
            np.stack([-s, c, zero], axis=-1),
            np.stack([zero, zero, one], axis=-1),
        ],
        axis=-2,
    )


def _rotate(vec: NDArray[np.float64], theta: NDArray[np.float64]) -> NDArray[np.float64]:
    """Aplica Rz(θ) linha a linha, sem montar as matrizes: (x, y, z) → (cx + sy, −sx + cy, z)."""
    c, s = np.cos(theta), np.sin(theta)
    x, y, z = vec[..., 0], vec[..., 1], vec[..., 2]
    return np.stack([c * x + s * y, -s * x + c * y, z], axis=-1)


def _omega_cross(r: NDArray[np.float64]) -> NDArray[np.float64]:
    """ω⊕ × r com ω⊕ = (0, 0, OMEGA_EARTH): (−ω y, ω x, 0)."""
    return np.stack(
        [-OMEGA_EARTH * r[..., 1], OMEGA_EARTH * r[..., 0], np.zeros_like(r[..., 2])], axis=-1
    )


def eci_to_ecef(
    r: ArrayLike, theta: ArrayLike, v: ArrayLike | None = None
) -> NDArray[np.float64] | tuple[NDArray[np.float64], NDArray[np.float64]]:
    """ECI → ECEF: r_ecef = Rz(θ)·r_eci; v_ecef = Rz(θ)·(v_eci − ω⊕ × r_eci).

    Args:
        r: posição ECI [m], shape (3,) ou (N, 3).
        theta: GMST θG [rad], escalar ou (N,).
        v: velocidade ECI [m/s], mesmo shape de ``r`` (opcional).

    Returns:
        r_ecef, ou (r_ecef, v_ecef) se ``v`` for dado; shapes iguais aos de entrada.
    """
    r = np.asarray(r, dtype=float)
    theta = np.asarray(theta, dtype=float)
    r_ecef = _rotate(r, theta)
    if v is None:
        return r_ecef
    v = np.asarray(v, dtype=float)
    return r_ecef, _rotate(v - _omega_cross(r), theta)


def ecef_to_eci(
    r_ecef: ArrayLike, theta: ArrayLike, v: ArrayLike | None = None
) -> NDArray[np.float64] | tuple[NDArray[np.float64], NDArray[np.float64]]:
    """ECEF → ECI: r_eci = Rz(θ)ᵀ·r_ecef; v_eci = Rz(θ)ᵀ·v_ecef + ω⊕ × r_eci.

    Args:
        r_ecef: posição ECEF [m], shape (3,) ou (N, 3).
        theta: GMST θG [rad], escalar ou (N,).
        v: velocidade ECEF [m/s], mesmo shape de ``r_ecef`` (opcional).

    Returns:
        r_eci, ou (r_eci, v_eci) se ``v`` for dado; shapes iguais aos de entrada.
    """
    r_ecef = np.asarray(r_ecef, dtype=float)
    theta = np.asarray(theta, dtype=float)
    r_eci = _rotate(r_ecef, -theta)  # Rz(θ)ᵀ = Rz(−θ)
    if v is None:
        return r_eci
    v = np.asarray(v, dtype=float)
    return r_eci, _rotate(v, -theta) + _omega_cross(r_eci)


def gmst_at(epoch: datetime, t: ArrayLike) -> FloatOrArray:
    """GMST θG [rad] em t segundos desde ``epoch`` (fórmula IAU-82 completa, via ``astrotime``).

    Args:
        epoch: ``datetime`` tz-aware (o t = 0 do propagador).
        t: segundos desde ``epoch``, escalar ou (N,).

    Returns:
        θG [rad] em [0, 2π), com o shape de ``t``.
    """
    return gmst(jd_from_epoch(datetime_to_jd(epoch), t))
