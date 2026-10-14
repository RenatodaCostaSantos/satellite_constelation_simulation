"""Taxas seculares de J2 sobre os elementos médios (Sem4 §3, eqs. 23–25 e 34).

O achatamento da Terra (J2) não muda a forma média da órbita: a, e e i médios ficam constantes,
e três ângulos derivam linearmente com o tempo (Sem4 §3.5):

    p = a(1 − e²),   k = 1,5·J2·(R⊕/p)²,   n = √(µ/a³)

    Ω̇ = −k·n·cos i                                (eq. 23, regressão do nó)
    ω̇ = ½·k·n·(5cos²i − 1)                        (eq. 24, rotação da linha dos apsides)
    ψ̇ = n·[1 + k·√(1 − e²)·(1 − 1,5·sin²i)]       (eq. 25, movimento médio perturbado)

O período nodal (de um nó ascendente ao seguinte) é Tnod = 2π/(ψ̇ + ω̇) (eq. 34).

Limitações: teoria secular de **primeira ordem só em J2**. Não inclui J3 (que cancela ω̇ na
órbita congelada; tratada na Semana 5), J4, nem as oscilações de curto período dos elementos
osculadores. As expressões não têm divisão por e nem por sin i, logo valem também para órbitas
circulares e equatoriais.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import J2, MU, R_EARTH

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
