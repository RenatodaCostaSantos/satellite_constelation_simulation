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
    # zeros_like/ones_like criam arrays de 0 e de 1 com o mesmo shape e dtype de theta (escalar
    # 0-D ou (N,)). O np.stack exige entradas de shapes iguais, então os 0 e 1 constantes da
    # matriz precisam ter o shape de c e s; com um 0 literal o stack falharia para θ (N,).
    zero, one = np.zeros_like(theta), np.ones_like(theta)
    # stack interno (axis=-1) monta cada linha [.., .., ..] no último eixo; o externo (axis=-2)
    # empilha as 3 linhas no penúltimo: (3, 3) para θ escalar, (N, 3, 3) para θ (N,).
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
    # vec[..., k]: os "..." (Ellipsis) valem "todos os eixos anteriores, quantos forem", e k
    # indexa o ÚLTIMO eixo, o das coordenadas. Para vec (3,) equivale a vec[k] (um número); para
    # (N, 3), a vec[:, k] (shape (N,)); para (P, N, 3), a vec[:, :, k] (shape (P, N)).
    x, y, z = vec[..., 0], vec[..., 1], vec[..., 2]
    # Operação inversa: axis=-1 recoloca as três componentes no último eixo, de modo que a
    # saída tem sempre o shape da entrada, sem nenhum if.
    return np.stack([c * x + s * y, -s * x + c * y, z], axis=-1)


def _omega_cross(r: NDArray[np.float64]) -> NDArray[np.float64]:
    """ω⊕ × r com ω⊕ = (0, 0, OMEGA_EARTH): (−ω y, ω x, 0)."""
    # r[..., 1] = y e r[..., 0] = x de cada posição (ver _rotate); zeros_like(r[..., 2]) dá a
    # componente z nula com o shape de uma coluna de r, como o stack exige.
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
    # v_eci é a velocidade do satélite em relação ao referencial ECI. A entrada v é a velocidade
    # em relação ao ECEF (que gira com a Terra), escrita nos eixos do ECEF. Duas etapas:
    # 1) _rotate(v, -theta): reescreve o mesmo vetor nos eixos do ECI (θ = GMST, o ângulo que a
    #    Terra girou; é uma troca de eixos, que não altera o tamanho do vetor);
    # 2) + ω⊕ × r_eci: soma a velocidade linear que o próprio ECEF tem naquele ponto por estar
    #    girando (≈ 430 m/s em São Bento; zero sobre o eixo z).
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
