"""Taxas seculares de J2 sobre os elementos médios (Sem4 §3, eqs. 23–25 e 34).

O achatamento da Terra (J2) não muda a forma média da órbita: a, e e i médios ficam constantes,
e três ângulos derivam linearmente com o tempo (Sem4 §3.5):

    p = a(1 − e²),   k = 1,5·J2·(R⊕/p)²,   n = √(µ/a³)

    Ω̇ = −k·n·cos i                                (eq. 23, regressão do nó)
    ω̇ = ½·k·n·(5cos²i − 1)                        (eq. 24, rotação da linha dos apsides)
    ψ̇ = n·[1 + k·√(1 − e²)·(1 − 1,5·sin²i)]       (eq. 25, movimento médio perturbado)

O período nodal (de um nó ascendente ao seguinte) é Tnod = 2π/(ψ̇ + ω̇) (eq. 34).

Limitações de ``secular_rates``: teoria secular de **primeira ordem só em J2**. Não inclui J3
(que cancela ω̇ na órbita congelada), J4, nem as oscilações de curto período dos elementos
osculadores. As expressões não têm divisão por e nem por sin i, logo valem também para órbitas
circulares e equatoriais.

J2, J3 e J4 de 1ª ordem (card S4-05)
------------------------------------
``zonal_mean_rates`` soma as taxas de J3 e J4 (Apêndice A.5 do prompt da Semana 4), deduzidas
pelo método do Sem4 §3: média do potencial ⟨R⟩ numa revolução e equações planetárias de
Lagrange. Com n = √(µ/a³), η = √(1 − e²), p = a·η², q_k = (R⊕/p)^k, s = sin i, c = cos i:

    Ω̇ = −(3/2) n J2 q2 c
         − (3/8) n J3 q3 e (15s² − 4)(c/s) sin ω
         − (15/32) n J4 q4 c [(14s² − 8) + e²(7s² − 6) + e² sin²ω (28s² − 12)]

    ω̇ = (3/4) n J2 q2 (4 − 5s²)
         − (3/8) n J3 q3 [s²(5s² − 4) + e²(4 − 35s²c²)] sin ω / (e s)
         − (15/128) n J4 q4 [(64 − 224s² + 168s⁴) + e²(48 − 112s² + 63s⁴)
                              + sin²ω (56s⁴ − 48s²) + e² sin²ω (48 − 280s² + 252s⁴)]

    ψ̇ = n + (3/4) n J2 q2 η (2 − 3s²)
         + (3/8) n J3 q3 η (1 − 4e²)(5s² − 4) s sin ω / e
         − (15/128) n J4 q4 η [(28s⁴ − 24s²) + sin²ω (48s² − 56s⁴)
                                + e²(35s⁴ − 60s² + 24) + e² sin²ω (140s⁴ − 120s²)]

    ė = (3/8) n J3 q3 η² (5s² − 4) s cos ω + (15/16) n J4 q4 η² e (7s² − 6) s² sin ω cos ω
    i̇ = −(3/8) n J3 q3 e (5s² − 4) c cos ω − (15/16) n J4 q4 e² (7s² − 6) s c sin ω cos ω

J2 não muda e nem i médios; J3 e J4 sim, de forma acoplada a ω (período longo). Os termos de J3
em ω̇ e ψ̇ vão com 1/e: é o que permite a órbita congelada (o termo de J3 em ω̇ cancela o de J2
numa excentricidade finita, ``orbit_design.frozen_eccentricity``), mas também faz as taxas
divergirem em e → 0 (com J3 ≠ 0 exige-se e > 0). O termo J3 de Ω̇ tem c/s: singular em órbita
equatorial. A soma ψ̇ + ω̇ (que dá o período nodal) também tem 1/e com J3, mas os dois termos
quase se cancelam.

O que fica de fora é a 2ª ordem (termos de J2²): na SSO de projeto ela responde pelo resíduo de
≈ −0,08% entre o Ω̇ desta teoria e o do Cowell J2–J4 (Sem4 §10, item 1). A calibração numérica
de ``orbit_design.sso_inclination_calibrated`` contorna isso para qualquer conjunto de forças.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import J2, J3, J4, MU, R_EARTH

FloatOrArray = float | NDArray[np.float64]


# @dataclass gera automaticamente __init__, __repr__ e __eq__ a partir dos campos anotados abaixo;
# frozen=True torna as instâncias imutáveis (atribuir a um campo levanta FrozenInstanceError).
@dataclass(frozen=True)
class SecularRates:
    """Contêiner (só dados, sem cálculo) das taxas seculares de J2 [rad/s] (Sem4 eqs. 23–25).

    Não é construído diretamente pelo usuário: é o valor de retorno de ``secular_rates``, que faz
    as contas e empacota os resultados aqui para acesso por nome em vez de posição numa tupla::

        rates = secular_rates(a=7000e3, e=0.001, i=np.radians(98.0))
        rates.raan_dot      # Ω̇ [rad/s]

    Attributes:
        n: movimento médio kepleriano √(µ/a³).
        raan_dot: Ω̇, regressão do nó ascendente (eq. 23).
        argp_dot: ω̇, rotação da linha dos apsides (eq. 24).
        psi_dot: ψ̇, movimento médio perturbado (eq. 25).

    Escalares para entrada escalar, arrays com o shape do broadcasting das entradas.
    """

    n: FloatOrArray
    raan_dot: FloatOrArray
    argp_dot: FloatOrArray
    psi_dot: FloatOrArray


def _as_scalar_if_0d(x: NDArray[np.float64]) -> FloatOrArray:
    return float(x) if x.ndim == 0 else x


def secular_rates(
    a: ArrayLike,
    e: ArrayLike,
    i: ArrayLike,
    mu: float = MU,
    j2: float = J2,
    re: float = R_EARTH,
) -> SecularRates:
    """Taxas seculares de J2 de primeira ordem (Sem4 eqs. 23, 24 e 25).

    Função de módulo (não é método de ``SecularRates``): calcula n, Ω̇, ω̇ e ψ̇ e devolve os
    quatro valores empacotados numa instância de ``SecularRates``.

    Args:
        a: semieixo maior médio [m], escalar ou array.
        e: excentricidade média, 0 ≤ e < 1, compatível com ``a`` por broadcasting.
        i: inclinação média [rad], compatível com ``a``.
        mu: parâmetro gravitacional [m³/s²].
        j2: harmônico zonal J2 (não normalizado).
        re: raio equatorial de referência de J2 [m].

    Returns:
        SecularRates com os campos n, raan_dot, argp_dot e psi_dot [rad/s]; campos escalares se
        todas as entradas forem escalares.

    Raises:
        ValueError: se algum a ≤ 0 ou e fora de [0, 1).
    """
    a, e, i = np.broadcast_arrays(*(np.asarray(x, dtype=float) for x in (a, e, i)))
    if np.any(a <= 0.0):
        raise ValueError("semieixo maior deve ser positivo")
    if np.any(e < 0.0) or np.any(e >= 1.0):
        raise ValueError("somente órbitas elípticas são suportadas: exige 0 <= e < 1")
    n = np.sqrt(mu / a**3)
    p = a * (1.0 - e**2)
    k = 1.5 * j2 * (re / p) ** 2
    cos_i, sin_i = np.cos(i), np.sin(i)
    raan_dot = -k * n * cos_i
    argp_dot = 0.5 * k * n * (5.0 * cos_i**2 - 1.0)
    psi_dot = n * (1.0 + k * np.sqrt(1.0 - e**2) * (1.0 - 1.5 * sin_i**2))
    return SecularRates(*(_as_scalar_if_0d(x) for x in (n, raan_dot, argp_dot, psi_dot)))


def nodal_period(
    a: ArrayLike,
    e: ArrayLike,
    i: ArrayLike,
    mu: float = MU,
    j2: float = J2,
    re: float = R_EARTH,
) -> FloatOrArray:
    """Período nodal Tnod = 2π/(ψ̇ + ω̇) [s] (Sem4 eq. 34), com as taxas de ``secular_rates``.

    O argumento de latitude u = ω + ν avança em média ψ̇ + ω̇; Tnod é o tempo entre duas
    passagens consecutivas pelo nó ascendente. ω̇ é o de J2 isolado (ver limitações do módulo).
    """
    rates = secular_rates(a, e, i, mu=mu, j2=j2, re=re)
    return _as_scalar_if_0d(np.asarray(2.0 * np.pi / (rates.psi_dot + rates.argp_dot)))


@dataclass(frozen=True)
class ZonalMeanRates:
    """Taxas médias de 1ª ordem com J2, J3 e J4 [rad/s] (Apêndice A.5; ver docstring do módulo).

    Valor de retorno de ``zonal_mean_rates``. Com J3 = J4 = 0, n, raan_dot, argp_dot e psi_dot são
    os de ``secular_rates`` e e_dot = i_dot = 0.

    Attributes:
        n: movimento médio kepleriano √(µ/a³).
        raan_dot: Ω̇.
        argp_dot: ω̇.
        psi_dot: ψ̇ (anomalia média).
        e_dot: ė [1/s].
        i_dot: i̇.
    """

    n: FloatOrArray
    raan_dot: FloatOrArray
    argp_dot: FloatOrArray
    psi_dot: FloatOrArray
    e_dot: FloatOrArray
    i_dot: FloatOrArray


def zonal_mean_rates(
    a: ArrayLike,
    e: ArrayLike,
    i: ArrayLike,
    argp: ArrayLike,
    mu: float = MU,
    j2: float = J2,
    j3: float = J3,
    j4: float = J4,
    re: float = R_EARTH,
) -> ZonalMeanRates:
    """Taxas médias de 1ª ordem com J2, J3 e J4 (Apêndice A.5 do prompt da Semana 4).

    Cada harmônico entra separado (efeito de 1ª ordem, sem termos cruzados): com ``j3 = j4 = 0``
    o resultado reproduz ``secular_rates``.

    Args:
        a: semieixo maior médio [m], escalar ou array.
        e: excentricidade média, 0 ≤ e < 1 (e > 0 se ``j3 ≠ 0``: os termos de J3 vão com 1/e).
        i: inclinação média [rad] (sin i ≠ 0 se ``j3 ≠ 0``).
        argp: argumento do perigeu médio ω [rad].
        mu: parâmetro gravitacional [m³/s²].
        j2, j3, j4: harmônicos zonais (não normalizados).
        re: raio equatorial de referência [m].

    Returns:
        ZonalMeanRates com n, raan_dot, argp_dot, psi_dot, e_dot e i_dot; campos escalares se
        todas as entradas forem escalares.

    Raises:
        ValueError: se algum a ≤ 0, e fora de [0, 1) ou, com ``j3 ≠ 0``, e = 0.
    """
    a, e, i, w = np.broadcast_arrays(*(np.asarray(x, dtype=float) for x in (a, e, i, argp)))
    if np.any(a <= 0.0):
        raise ValueError("semieixo maior deve ser positivo")
    if np.any(e < 0.0) or np.any(e >= 1.0):
        raise ValueError("somente órbitas elípticas são suportadas: exige 0 <= e < 1")
    if j3 != 0.0 and np.any(e == 0.0):
        raise ValueError("com J3 as taxas médias divergem em e = 0: exige e > 0")
    n = np.sqrt(mu / a**3)
    e2 = e**2
    eta = np.sqrt(1.0 - e2)
    p = a * eta**2
    q2, q3, q4 = (re / p) ** 2, (re / p) ** 3, (re / p) ** 4
    s, c = np.sin(i), np.cos(i)
    s2 = s**2
    s4 = s2**2
    sw, cw = np.sin(w), np.cos(w)
    sw2 = sw**2

    # J2: secular puro (Sem4 eqs. 23–25)
    raan_dot = -1.5 * n * j2 * q2 * c
    argp_dot = 0.75 * n * j2 * q2 * (4.0 - 5.0 * s2)
    psi_dot = n + 0.75 * n * j2 * q2 * eta * (2.0 - 3.0 * s2)
    e_dot = np.zeros_like(n)
    i_dot = np.zeros_like(n)

    if j3 != 0.0:
        k3 = 0.375 * n * j3 * q3
        raan_dot = raan_dot - k3 * e * (15.0 * s2 - 4.0) * (c / s) * sw
        argp_dot = argp_dot - k3 * (s2 * (5.0 * s2 - 4.0) + e2 * (4.0 - 35.0 * s2 * c**2)) * sw / (
            e * s
        )
        psi_dot = psi_dot + k3 * eta * (1.0 - 4.0 * e2) * (5.0 * s2 - 4.0) * s * sw / e
        e_dot = e_dot + k3 * eta**2 * (5.0 * s2 - 4.0) * s * cw
        i_dot = i_dot - k3 * e * (5.0 * s2 - 4.0) * c * cw

    if j4 != 0.0:
        k4 = n * j4 * q4
        raan_dot = raan_dot - (15.0 / 32.0) * k4 * c * (
            (14.0 * s2 - 8.0) + e2 * (7.0 * s2 - 6.0) + e2 * sw2 * (28.0 * s2 - 12.0)
        )
        argp_dot = argp_dot - (15.0 / 128.0) * k4 * (
            (64.0 - 224.0 * s2 + 168.0 * s4)
            + e2 * (48.0 - 112.0 * s2 + 63.0 * s4)
            + sw2 * (56.0 * s4 - 48.0 * s2)
            + e2 * sw2 * (48.0 - 280.0 * s2 + 252.0 * s4)
        )
        psi_dot = psi_dot - (15.0 / 128.0) * k4 * eta * (
            (28.0 * s4 - 24.0 * s2)
            + sw2 * (48.0 * s2 - 56.0 * s4)
            + e2 * (35.0 * s4 - 60.0 * s2 + 24.0)
            + e2 * sw2 * (140.0 * s4 - 120.0 * s2)
        )
        e_dot = e_dot + (15.0 / 16.0) * k4 * eta**2 * e * (7.0 * s2 - 6.0) * s2 * sw * cw
        i_dot = i_dot - (15.0 / 16.0) * k4 * e2 * (7.0 * s2 - 6.0) * s * c * sw * cw

    return ZonalMeanRates(
        *(_as_scalar_if_0d(x) for x in (n, raan_dot, argp_dot, psi_dot, e_dot, i_dot))
    )
