"""Decaimento orbital por arrasto em órbita circular: taxa, perda por revolução e tempo de vida.

A cadeia causal do Sem3 (ρ → F_D → dE/dt → da/dt → h(t)) em forma fechada, para órbita circular
(r = a, v = √(µ/a)) e altitude esférica h = a − R⊕:

    da/dt = −C_D·(A/m)·ρ(h)·√(µ·a)            (Sem3 eq. 14)
    Δa_rev = −2π·C_D·(A/m)·ρ(h)·a²            (Sem3 eq. 15: µ desaparece)
    dh/dt = −C_D·(A/m)·ρ(h)·√(µ·(R⊕ + h))     (Sem3 eq. 22, com ρ(h) da tabela por partes)

O tempo de vida integra a eq. 22 até a altitude de reentrada. A estimativa de uma linha
t_vida ≈ H/|ḣ(h0)| (Sem3 eq. 24) é o "tempo para descer uma escala de altura"; a integração com
a atmosfera por partes dá, por exemplo, 4,64 anos a partir de 500 km com C_D·A/m = 0,011, contra
5,0 da eq. 24 (Sem3 §5.1).

Todas as funções recebem o produto ``cd_a_m`` = C_D·A/m [m²/kg], o único parâmetro que a
trajetória enxerga (Sem3 §2.4). Ver ``satsim.forces.drag`` para a aceleração usada no Cowell.
"""

import math

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp

from satsim.constants import MU, R_EARTH, SECONDS_PER_DAY
from satsim.forces.atmosphere import ExponentialAtmosphere

FloatOrArray = float | NDArray[np.float64]

SECONDS_PER_YEAR = 365.25 * SECONDS_PER_DAY
"""Ano juliano [s] (= 3,15576e7 s; o código do Sem3 §6.1 usa 3,156e7)."""

_MAX_LIFETIME_YEARS = 1e5  # limite da integração: acima disso o tempo de vida é tratado como ∞


def _atmosphere(atmosphere: ExponentialAtmosphere | None) -> ExponentialAtmosphere:
    return atmosphere if atmosphere is not None else ExponentialAtmosphere()


def _as_scalar_if_0d(x: NDArray[np.float64]) -> FloatOrArray:
    return float(x) if x.ndim == 0 else x


def decay_rate(
    h_m: ArrayLike, cd_a_m: float, atmosphere: ExponentialAtmosphere | None = None
) -> FloatOrArray:
    """Taxa de decaimento da altitude dh/dt = da/dt [m/s] em órbita circular (Sem3 eq. 14).

    dh/dt = −C_D·(A/m)·ρ(h)·√(µ·(R⊕ + h)), vetorizada sobre ``h_m``; negativa (a órbita encolhe).

    Args:
        h_m: altitude [m], escalar ou array.
        cd_a_m: C_D·A/m [m²/kg].
        atmosphere: modelo de densidade (padrão: ``ExponentialAtmosphere()``, Vallado).
    """
    h = np.asarray(h_m, dtype=float)
    rho = np.asarray(_atmosphere(atmosphere).density(h))
    return _as_scalar_if_0d(-cd_a_m * rho * np.sqrt(MU * (R_EARTH + h)))


def decay_per_rev(
    h_m: ArrayLike, cd_a_m: float, atmosphere: ExponentialAtmosphere | None = None
) -> FloatOrArray:
    """Perda de semieixo por revolução Δa_rev = −2π·C_D·(A/m)·ρ·a² [m] (Sem3 eq. 15).

    Supõe ρ constante ao longo da volta (órbita circular). Ex.: 400 km, C_D·A/m = 0,011 →
    ≈ −11,8 m por revolução.
    """
    h = np.asarray(h_m, dtype=float)
    rho = np.asarray(_atmosphere(atmosphere).density(h))
    return _as_scalar_if_0d(-2.0 * math.pi * cd_a_m * rho * (R_EARTH + h) ** 2)


def lifetime(
    h0_m: float,
    cd_a_m: float,
    atmosphere: ExponentialAtmosphere | None = None,
    h_reentry_m: float = 120e3,
) -> float:
    """Tempo de vida [anos julianos] de uma órbita circular, integrando a Sem3 eq. 22.

    Integra dh/dt de ``h0_m`` até ``h_reentry_m`` com ``solve_ivp`` e um evento terminal (como
    no código do Sem3 §6.1). A 120 km a reentrada é iminente (Sem3 Figura 2).

    Args:
        h0_m: altitude inicial [m].
        cd_a_m: C_D·A/m [m²/kg]; 0 dá ``math.inf``.
        atmosphere: modelo de densidade (padrão: Vallado).
        h_reentry_m: altitude que encerra a integração [m].

    Returns:
        Tempo de vida [anos]; ``math.inf`` se a reentrada não ocorrer em 1e5 anos. 0 se
        ``h0_m ≤ h_reentry_m``.
    """
    if h0_m <= h_reentry_m:
        return 0.0
    if cd_a_m <= 0.0:
        return math.inf
    atm = _atmosphere(atmosphere)

    def rhs(t: float, y: NDArray[np.float64]) -> list[float]:
        return [float(decay_rate(y[0], cd_a_m, atm))]

    def reentry(t: float, y: NDArray[np.float64]) -> float:
        return float(y[0] - h_reentry_m)

    reentry.terminal = True  # type: ignore[attr-defined]
    reentry.direction = -1  # type: ignore[attr-defined]
    t_max = _MAX_LIFETIME_YEARS * SECONDS_PER_YEAR
    sol = solve_ivp(
        rhs, (0.0, t_max), [h0_m], events=reentry, rtol=1e-8, atol=1e-2, max_step=t_max / 1e3
    )
    if sol.status != 1:
        return math.inf
    return float(sol.t_events[0][0]) / SECONDS_PER_YEAR
