"""Projeto de órbita: condição heliossíncrona (SSO) e a Tabela 1 do Sem4 (§4).

Uma órbita é heliossíncrona quando o plano orbital precessa à mesma taxa média com que o Sol se
move em ascensão reta, Ω̇ = OMEGA_SUN ≈ 0,9856°/dia (Sem4 eq. 26): o ângulo entre o plano e a
direção do Sol (o LTAN, hora local do nó ascendente) fica constante ao longo do ano. Igualando a
taxa secular de J2 (Sem4 eq. 23) a OMEGA_SUN e isolando cos i (Sem4 eq. 27):

    cos i = −2·Ω☉·a^(7/2)·(1 − e²)² / (3·√µ·J2·R⊕²)

O sinal negativo faz toda SSO ser retrógrada (i > 90°). Como |cos i| ≤ 1, existe um teto de
semieixo, atingido em i = 180°:

    a_max = (3·√µ·J2·R⊕² / (2·Ω☉))^(2/7) ≈ 12 352 km    (h_max ≈ 5 974 km)

Teoria secular de primeira ordem só em J2 (ver ``satsim.secular``).
"""

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import J2, MU, OMEGA_EARTH, OMEGA_SUN, R_EARTH, SECONDS_PER_DAY
from satsim.secular import nodal_period, secular_rates

FloatOrArray = float | NDArray[np.float64]


def _as_scalar_if_0d(x: NDArray[np.float64]) -> FloatOrArray:
    return float(x) if x.ndim == 0 else x


def sso_inclination(a: ArrayLike, e: ArrayLike = 0.0) -> FloatOrArray:
    """Inclinação heliossíncrona i [rad] para semieixo ``a`` e excentricidade ``e`` (Sem4 eq. 27).

    cos i = −2·Ω☉·a^(7/2)·(1 − e²)² / (3·√µ·J2·R⊕²), com Ω☉ = ``OMEGA_SUN``.

    Args:
        a: semieixo maior [m], escalar ou array.
        e: excentricidade, 0 ≤ e < 1, compatível com ``a`` por broadcasting.

    Returns:
        i [rad] em (π/2, π]; float para entradas escalares.

    Raises:
        ValueError: se algum |cos i| > 1 (a acima do teto ``sso_max_semi_major_axis``) ou e fora
            de [0, 1).
    """
    a, e = np.broadcast_arrays(np.asarray(a, dtype=float), np.asarray(e, dtype=float))
    if np.any(e < 0.0) or np.any(e >= 1.0):
        raise ValueError("somente órbitas elípticas são suportadas: exige 0 <= e < 1")
    cos_i = -2.0 * OMEGA_SUN * a**3.5 * (1.0 - e**2) ** 2 / (3.0 * math.sqrt(MU) * J2 * R_EARTH**2)
    if np.any(np.abs(cos_i) > 1.0):
        raise ValueError(
            "não existe órbita heliossíncrona: |cos i| > 1 (semieixo acima do teto "
            f"{sso_max_semi_major_axis() / 1e3:.0f} km)"
        )
    return _as_scalar_if_0d(np.arccos(cos_i))


def sso_max_semi_major_axis() -> float:
    """Teto de semieixo de uma SSO circular, a_max = (3·√µ·J2·R⊕²/(2·Ω☉))^(2/7) [m] (Sem4 §4.3).

    É o semieixo em que cos i = −1 (i = 180°); acima dele J2 não consegue precessar o plano
    na taxa do Sol.
    """
    return (3.0 * math.sqrt(MU) * J2 * R_EARTH**2 / (2.0 * OMEGA_SUN)) ** (2.0 / 7.0)


@dataclass(frozen=True)
class SsoRow:
    """Uma linha da Tabela 1 do Sem4 (§4.3): SSO circular na altitude h.

    h_km: altitude [km]; i_deg: inclinação heliossíncrona [°]; nodal_period_s: período nodal
    Tnod [s] (Sem4 eq. 34); rev_per_day: revoluções por dia nodal (Sem4 eq. 37).
    """

    h_km: float
    i_deg: float
    nodal_period_s: float
    rev_per_day: float


def sso_row(h: float, e: float = 0.0) -> SsoRow:
    """Linha da Tabela 1 (Sem4 §4.3) para a altitude ``h`` [m].

    a = R⊕ + h, i = ``sso_inclination(a, e)``, Tnod = ``nodal_period(a, e, i)`` e
    rev/dia = (ψ̇ + ω̇)/(ω⊕ − Ω̇) (Sem4 eq. 37): voltas nodais por dia nodal, isto é, por
    intervalo entre duas passagens do mesmo meridiano sob o plano orbital.

    Args:
        h: altitude [m] (a = R⊕ + h).
        e: excentricidade (a Tabela 1 usa órbitas circulares).

    Returns:
        SsoRow com h em km e i em graus.
    """
    a = R_EARTH + h
    i = sso_inclination(a, e)
    rates = secular_rates(a, e, i)
    rev_per_day = (rates.psi_dot + rates.argp_dot) / (OMEGA_EARTH - rates.raan_dot)
    return SsoRow(
        h_km=h / 1e3,
        i_deg=math.degrees(i),
        nodal_period_s=nodal_period(a, e, i),
        rev_per_day=rev_per_day,
    )


def ltan_drift_min_per_year(raan_dot: ArrayLike) -> FloatOrArray:
    """Deriva do LTAN [min/ano] para uma taxa de precessão Ω̇ [rad/s] (Sem4 §6.1).

    O LTAN anda (Ω̇ − Ω☉) por unidade de tempo; convertido em hora local (15°/h) e acumulado em
    um ano trópico: (Ω̇ − Ω☉)·365,2422 dias / (15°/h) · 60 min/h.
    """
    deg_per_year = np.degrees(np.asarray(raan_dot, dtype=float) - OMEGA_SUN) * (
        365.2422 * SECONDS_PER_DAY
    )
    return _as_scalar_if_0d(deg_per_year / 15.0 * 60.0)
