"""Aceleração perturbadora do harmônico zonal J2 (Sem4 §9.3).

O achatamento da Terra soma ao potencial central µ/r a função perturbadora (Sem4 §2.2)

    R(r) = −µ·J2·R⊕²·(3 sin²φ − 1)/(2 r³),   sin φ = z/r (latitude geocêntrica),

e a aceleração é o seu gradiente, a = ∇R. Em coordenadas cartesianas ECI (Sem4 §9.3):

    c = 1,5·J2·µ·R⊕²/r⁵
    a = c·[x(5z²/r² − 1), y(5z²/r² − 1), z(5z²/r² − 3)]

Ela é uma correção pequena (~1e-3) à gravidade central, que continua dominante: a aceleração
total sempre aponta para dentro. No equador a correção é radial e para dentro
(−1,5·J2·µ·R⊕²/r⁴), deixando a gravidade ali um pouco mais forte; sobre os polos é radial e para
fora (+3·J2·µ·R⊕²/r⁴), deixando-a um pouco mais fraca. Nesses pontos ela é radial e não exerce
torque. Nas latitudes intermediárias surge uma componente não central, perpendicular a r e
voltada para o plano do equador, de módulo 3·J2·µ·R⊕²·sin φ·cos φ/r⁴ (máxima em φ = 45°). O
torque r × a dessa componente é horizontal (sem componente z, por isso ℓz se conserva) e gira ℓ
em torno de z: é a precessão do nó, Ω̇ (Sem4 eq. 23).

A energia conservada sob gravidade central + J2 é ε = v²/2 − µ/r − R(r).
"""

from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import J2, MU, R_EARTH


class AccelerationModel(Protocol):
    """Aceleração perturbadora usada pelo ``CowellPropagator``.

    ``__call__(t, r, v) → a``: t em segundos desde a época do propagador; r [m] e v [m/s] em ECI
    com shape (3,); devolve a aceleração [m/s²] em ECI com shape (3,). O termo central −µr/r³
    não entra aqui: é somado pelo propagador.
    """

    def __call__(
        self, t: float, r: NDArray[np.float64], v: NDArray[np.float64]
    ) -> NDArray[np.float64]: ...


def j2_acceleration(
    r: ArrayLike, mu: float = MU, j2: float = J2, re: float = R_EARTH
) -> NDArray[np.float64]:
    """Aceleração de J2 em ECI (Sem4 §9.3), vetorizada.

    Args:
        r: posição ECI [m], shape (3,) ou (N, 3).
        mu: parâmetro gravitacional [m³/s²].
        j2: harmônico zonal J2 (não normalizado).
        re: raio equatorial de referência de J2 [m].

    Returns:
        Aceleração [m/s²] com o mesmo shape de ``r``.
    """
    r = np.asarray(r, dtype=float)
    x, y, z = r[..., 0], r[..., 1], r[..., 2]
    r2 = x * x + y * y + z * z
    c = 1.5 * j2 * mu * re**2 / (r2 * r2 * np.sqrt(r2))
    s = 5.0 * z * z / r2
    return np.stack([c * x * (s - 1.0), c * y * (s - 1.0), c * z * (s - 3.0)], axis=-1)


def j2_potential(
    r: ArrayLike, mu: float = MU, j2: float = J2, re: float = R_EARTH
) -> NDArray[np.float64]:
    """Função perturbadora de J2 R(r) = −µ·J2·R⊕²·(3 sin²φ − 1)/(2 r³) [m²/s²] (Sem4 §2.2).

    A aceleração é ∇R; a energia conservada é v²/2 − µ/r − R. Shape (3,) → (), (N, 3) → (N,).
    """
    r = np.asarray(r, dtype=float)
    r2 = np.sum(r * r, axis=-1)
    sin2_phi = r[..., 2] ** 2 / r2
    return -mu * j2 * re**2 * (3.0 * sin2_phi - 1.0) / (2.0 * r2 * np.sqrt(r2))


class J2Acceleration:
    """Modelo de aceleração de J2 (implementa ``AccelerationModel``).

    Args:
        mu: parâmetro gravitacional [m³/s²].
        j2: harmônico zonal J2 (não normalizado).
        re: raio equatorial de referência de J2 [m].
    """

    def __init__(self, mu: float = MU, j2: float = J2, re: float = R_EARTH) -> None:
        self.mu = mu
        self.j2 = j2
        self.re = re

    def __call__(
        self, t: float, r: NDArray[np.float64], v: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        return j2_acceleration(r, mu=self.mu, j2=self.j2, re=self.re)

    def potential(self, r: ArrayLike) -> NDArray[np.float64]:
        """Função perturbadora R(r) com os parâmetros deste modelo (ver ``j2_potential``)."""
        return j2_potential(r, mu=self.mu, j2=self.j2, re=self.re)
