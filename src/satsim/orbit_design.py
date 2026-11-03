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

Órbita congelada e SSO calibrada (card S4-05)
---------------------------------------------
``frozen_eccentricity`` acha a excentricidade em que ω̇ = 0 com ω = 90° (Sem4 §6.3), agora com
J4 além de J2 e J3 (taxas do Apêndice A.5). ``sso_inclination_calibrated`` resolve a condição
heliossíncrona **no próprio Cowell**: a fórmula de 1ª ordem só em J2 (eq. 27) erra Ω̇ em ≈ −0,3%
com J2–J4 (LTAN −4,5 min/ano na órbita de projeto), e mesmo a teoria J2–J4 de 1ª ordem deixa
−0,08% (termos de J2²). A calibração numérica vale para qualquer conjunto de forças.

Órbita de repetição e faseamento da constelação (card S5-01)
------------------------------------------------------------
``repeat_orbit`` resolve juntas a condição heliossíncrona (eq. 27) e a de repetição
N·Tnod = D·Tdia,nod (eq. 37) para (a, i), só com J2 de 1ª ordem, como o solver do Sem4 §9.1:
reproduz a Tabela 3 do Sem4. ``constellation_phases`` dá as fases ψ0 dos P satélites num único
plano (Sem4 eq. 42, Tabela 4).
"""

import logging
import math
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq

from satsim.analysis import fit_rate, osculating_elements
from satsim.constants import J2, J3, J4, MU, OMEGA_EARTH, OMEGA_SUN, R_EARTH, SECONDS_PER_DAY
from satsim.elements import KeplerianElements
from satsim.ephemeris import sun_ra_dec
from satsim.forces import AccelerationModel
from satsim.mean_osc import conversion_models, mean_to_osculating
from satsim.propagators.cowell import CowellPropagator
from satsim.secular import nodal_period, secular_rates, zonal_mean_rates

logger = logging.getLogger(__name__)

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


def repeat_orbit(n_revs: int, n_days: int, e: float = 0.0) -> tuple[float, float]:
    """SSO de repetição: (a, i) com N revoluções nodais em D dias nodais (Sem4 eqs. 27 e 37).

    Resolve em a, por ``brentq`` (solver do Sem4 §9.1), a condição de repetição

        N/D = (ψ̇ + ω̇)/(ω⊕ − Ω̇)        (Sem4 eq. 37)

    com i = ``sso_inclination(a, e)`` (eq. 27) e as taxas seculares de 1ª ordem só em J2
    (``secular_rates``, eqs. 23–25). Numa SSO, Ω̇ = Ω☉ e o dia nodal é o dia solar médio (§7.2).
    Teoria só J2: o semieixo que fecha o traço no Cowell J2–J4 é ≈ 19 m maior (calibrado no
    S5-02). Com e = 0 (padrão) a correção da excentricidade congelada é de 1,4 cm em a (§7.3).

    Args:
        n_revs: N, revoluções nodais no ciclo.
        n_days: D, dias do ciclo. Se gcd(N, D) > 1 a órbita existe, mas repete em D/gcd dias
            (degenerada, Sem4 §7.4); a função não recusa, quem chama decide.
        e: excentricidade média.

    Returns:
        (a [m], i [rad]).

    Raises:
        ValueError: se N, D não forem positivos ou se não houver solução entre 150 e 2 500 km
            de altitude (faixa do Sem4 §9.1).
    """
    if n_revs <= 0 or n_days <= 0:
        raise ValueError("n_revs e n_days devem ser positivos")
    target = n_revs / n_days

    def revs_per_day_minus_target(a: float) -> float:
        i = sso_inclination(a, e)
        rates = secular_rates(a, e, i)
        return (rates.psi_dot + rates.argp_dot) / (OMEGA_EARTH - rates.raan_dot) - target

    lo, hi = R_EARTH + 150e3, R_EARTH + 2500e3
    if not revs_per_day_minus_target(lo) * revs_per_day_minus_target(hi) < 0.0:
        raise ValueError(f"sem órbita de repetição {n_revs}/{n_days} entre 150 e 2 500 km")
    a = float(brentq(revs_per_day_minus_target, lo, hi, xtol=1e-6))
    return a, float(sso_inclination(a, e))


def constellation_phases(n_sats: int) -> NDArray[np.float64]:
    """Fases iniciais ψ0 [rad] de P satélites igualmente faseados num plano (Sem4 eq. 42).

    ψ0_j = −j·2π/P (mod 2π), j = 0, …, P − 1: S1 em 0 e cada satélite seguinte 360°/P atrás
    (P = 6: 0°, 300°, 240°, 180°, 120°, 60°, Tabela 4). O satélite j cruza o nó ascendente
    j·Tnod/P depois de S1, e o seu nó cai j·D/P passos de grade a oeste (eq. 43). Em órbita
    quase circular, fasear por ψ ou por u = ω + ν dá o mesmo resultado a menos de ~e.

    Args:
        n_sats: P ≥ 1.

    Returns:
        ψ0 [rad] em [0, 2π), shape (P,).
    """
    if n_sats < 1:
        raise ValueError("n_sats deve ser ≥ 1")
    return np.mod(-2.0 * np.pi * np.arange(n_sats) / n_sats, 2.0 * np.pi)


def frozen_eccentricity(a: float, i: float, j4: bool = True) -> float:
    """Excentricidade congelada: raiz de ω̇(e) = 0 em ω = 90° (Sem4 §6.3, taxas do Apêndice A.5).

    Com ``j4=False`` usa só J2 e J3 e reproduz a Sem4 eq. 33, e = −(J3/2J2)(R⊕/a)·sin i, a menos
    dos termos em e² e do p = a(1 − e²) das taxas completas (≈ 6e-6 relativo na órbita de
    projeto). Com J4 (padrão) o termo de J4 em ω̇ desloca a raiz em ≈ +0,3%.

    Args:
        a: semieixo maior médio [m].
        i: inclinação média [rad].
        j4: incluir J4 nas taxas.

    Returns:
        e congelada [-].

    Raises:
        ValueError: se não houver raiz no intervalo de busca (ex.: i perto do equador).
    """
    j4_value = J4 if j4 else 0.0

    def argp_dot(e: float) -> float:
        return float(zonal_mean_rates(a, e, i, math.pi / 2.0, j4=j4_value).argp_dot)

    e_eq33 = abs(-J3 / (2.0 * J2) * (R_EARTH / a) * math.sin(i))
    lo, hi = 1e-3 * e_eq33, min(10.0 * e_eq33, 0.5)
    if not argp_dot(lo) * argp_dot(hi) < 0.0:
        raise ValueError("ω̇(e) = 0 sem raiz no intervalo de busca: não há órbita congelada")
    return float(brentq(argp_dot, lo, hi, xtol=1e-15, rtol=1e-14))


def raan_for_ltan(ltan_h: ArrayLike, jd: ArrayLike) -> FloatOrArray:
    """Ω [rad] em [0, 2π) que dá o LTAN pedido na data ``jd`` (inverso da Sem4 eq. 28).

    Ω = α☉ + 15°·(LTAN − 12h), com α☉ a ascensão reta do Sol (``ephemeris.sun_ra_dec``, equador e
    equinócio da data). Ex.: LTAN 22h30 em 2026-10-23 00:00 UTC → Ω = 207,5235° + 157,5° =
    5,0235° (Apêndice A.7).

    Args:
        ltan_h: hora local do nó ascendente [h].
        jd: Data Juliana (UTC usado como TT).
    """
    ra_sun, _ = sun_ra_dec(jd)
    raan = ra_sun + np.radians(15.0 * (np.asarray(ltan_h, dtype=float) - 12.0))
    return _as_scalar_if_0d(np.mod(raan, 2.0 * np.pi))


def sso_inclination_calibrated(
    a: float,
    e: float,
    epoch: datetime,
    models: Sequence[AccelerationModel],
    days: float = 10.0,
    argp: float = math.pi / 2.0,
    raan: float = 0.0,
    tol_deg_per_day: float = 1e-6,
    max_iter: int = 8,
) -> float:
    """Inclinação média heliossíncrona calibrada no Cowell(models) por secante em i.

    Para cada i candidata: elementos médios (a, e, i, Ω, ω, ψ = 0) → ``mean_to_osculating`` →
    Cowell(models) por ``days`` dias → Ω̇ = ``fit_rate`` do Ω osculador (amostras a cada 60 s).
    A secante zera Ω̇ − OMEGA_SUN partindo da fórmula de 1ª ordem só J2 (``sso_inclination``,
    eq. 27); o 2º ponto escala cos i pela razão OMEGA_SUN/Ω̇ medido (Ω̇ ∝ cos i).

    A conversão médio → osculador usa só a gravidade zonal de ``models``
    (``mean_osc.conversion_models``); arrasto, SRP e Sol/Lua entram só na propagação.

    Args:
        a: semieixo maior médio [m].
        e: excentricidade média (em geral a congelada, ``frozen_eccentricity``).
        epoch: época, ``datetime`` tz-aware.
        models: forças do Cowell (ex.: ``(ZonalAcceleration(),)``).
        days: duração da medida de Ω̇ [dias].
        argp: ω médio [rad] (90° na órbita congelada).
        raan: Ω médio na época [rad]; só importa com Sol/Lua.
        tol_deg_per_day: tolerância em |Ω̇ − OMEGA_SUN| [°/dia] (1e-6 °/d ≈ 0,0015 min/ano de
            LTAN).
        max_iter: número máximo de medidas no Cowell.

    Returns:
        i calibrada [rad]. Se ``max_iter`` se esgotar, devolve a última e registra um aviso.
    """
    # Convenção do projeto: a conversão médio → osculador usa só a gravidade zonal (J2–J4) de
    # `models`; arrasto, SRP e Sol/Lua ficam só na propagação (mean_osc.conversion_models).
    conv_models = conversion_models(models)
    # Instantes de amostragem da medida de Ω̇: de 0 a `days` dias, a cada 60 s (o +1 s garante que
    # o último instante, days·86400, entre no arange). 10 dias → 14 401 amostras.
    t = np.arange(0.0, days * SECONDS_PER_DAY + 1.0, 60.0)
    # Tolerância em Ω̇ convertida de °/dia para rad/s, a unidade de fit_rate e de OMEGA_SUN.
    tol = math.radians(tol_deg_per_day) / SECONDS_PER_DAY

    def residual(i: float) -> float:
        """f(i) = Ω̇ medido no Cowell com inclinação média i, menos a taxa do Sol [rad/s].

        A SSO calibrada é a raiz f(i) = 0. Cada chamada é uma "medida": ≈ 1,2 s de conversão
        mais ≈ 1,7 s de Cowell em 10 dias.
        """
        # 1. Elementos MÉDIOS candidatos: o projeto (a, e, ω, Ω) com a inclinação i em teste;
        #    ψ = 0 (fase na órbita) é arbitrário, Ω̇ não depende dela.
        mean = KeplerianElements(a=a, e=e, i=i, raan=raan, argp=argp, psi=0.0)
        # 2. Médio → osculador (S4-04): o estado instantâneo cuja média no Cowell é `mean`.
        #    Sem isso, o semieixo médio ficaria ≈ 9,5 km alto e Ω̇ sairia errado.
        osc = mean_to_osculating(mean, epoch, conv_models).elements
        # 3. Cowell com TODAS as forças pedidas, partindo do estado osculador.
        prop = CowellPropagator(osc, epoch, models)
        # 4. Ω osculador nos instantes t (oscila a cada volta) e ajuste de uma reta Ω(t) por
        #    mínimos quadrados: a inclinação da reta é o Ω̇ médio. As oscilações de curto período
        #    se anulam no ajuste de 10 dias (≈ 150 voltas).
        raan_dot = fit_rate(t, osculating_elements(prop, t).raan)
        logger.debug("i = %.6f°: Ω̇ = %.7f °/d", math.degrees(i), math.degrees(raan_dot) * 86400)
        # 5. Quanto falta (ou sobra) para a precessão acompanhar o Sol (0,985647 °/dia).
        return raan_dot - OMEGA_SUN

    # Ponto de partida i0: a fórmula analítica de 1ª ordem só com J2 (Sem4 eq. 27). Na órbita de
    # projeto dá 97,4396°, ≈ 0,024° abaixo da resposta com J2–J4.
    i0 = float(sso_inclination(a, e))
    # 1ª medida: com J2–J4 o Cowell precessa ≈ 0,3% devagar demais em i0 (f0 < 0).
    f0 = residual(i0)
    # Se por acaso i0 já satisfaz a tolerância (ex.: models só com J2), não há o que corrigir.
    if abs(f0) < tol:
        return i0
    # 2º ponto, para a secante ter dois pontos: Ω̇ ∝ cos i (eq. 23), então escalar cos i pela
    # razão (taxa desejada)/(taxa medida) corrige a maior parte do erro numa tacada.
    # f0 + OMEGA_SUN é o próprio Ω̇ medido em i0.
    i1 = math.acos(math.cos(i0) * OMEGA_SUN / (f0 + OMEGA_SUN))
    # Método da secante: até max_iter medidas no total (a 1ª já foi feita acima).
    for _ in range(max_iter - 1):
        # Mede o resíduo no ponto atual.
        f1 = residual(i1)
        # Convergiu: |Ω̇ − Ω☉| < tol (1e-6 °/d ≈ 0,0015 min/ano de deriva da LTAN).
        if abs(f1) < tol:
            return i1
        # Secante: a reta pelos dois últimos pontos (i0, f0) e (i1, f1) cruza zero em
        # i1 − f1·(i1 − i0)/(f1 − f0). Esse cruzamento vira o novo ponto; o atual vira o
        # anterior. Como f(i) é quase linear perto da raiz, bastam 1 ou 2 passos.
        i0, f0, i1 = i1, f1, i1 - f1 * (i1 - i0) / (f1 - f0)
    # Esgotou as medidas sem atingir a tolerância: avisa e devolve a melhor estimativa.
    logger.warning("sso_inclination_calibrated não convergiu em %d medidas", max_iter)
    return i1
