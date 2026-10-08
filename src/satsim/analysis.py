"""Análise de trajetórias: elementos osculadores, elementos médios e taxas ajustadas.

Elementos osculadores são os elementos keplerianos do estado instantâneo (r, v): a órbita de dois
corpos que o satélite seguiria se as perturbações sumissem naquele instante. Sob J2 eles oscilam
a cada revolução (curto período, ~10 km no semieixo em LEO) em torno de valores médios que só
derivam secularmente (Sem4 §3). Este módulo mede essas duas componentes numa trajetória
propagada, para comparar com a teoria secular (``satsim.secular``).
"""

import numpy as np
from numpy.typing import ArrayLike

from satsim.elements import KeplerianElements, state_to_elements
from satsim.propagators.base import Propagator, as_time_array


def osculating_elements(propagator: Propagator, t: ArrayLike) -> KeplerianElements:
    """Elementos osculadores ao longo da trajetória de ``propagator``.

    Args:
        propagator: qualquer ``Propagator``.
        t: tempos [s] desde a época do propagador, shape () ou (N,).

    Returns:
        KeplerianElements com campos de shape (N,); raan, argp e psi em [0, 2π) (use
        ``np.unwrap`` para séries contínuas).
    """
    r, v = propagator.propagate(as_time_array(t))
    return state_to_elements(r, v)


def mean_elements(
    t: ArrayLike, elements: KeplerianElements, period: float
) -> tuple[float, float, float]:
    """Médias de a, e e i sobre o maior número inteiro de revoluções contido em ``t``.

    A média sobre revoluções inteiras elimina as oscilações de curto período (período orbital e
    harmônicos), restando a parte média. A janela é [t0, t0 + K·period], com
    K = ⌊(t_N − t0)/period⌋; a média é a integral trapezoidal dividida pela duração coberta pelas
    amostras da janela.

    Args:
        t: tempos [s], shape (N,), crescentes.
        elements: elementos (osculadores) amostrados em ``t``, campos de shape (N,).
        period: período da oscilação a remover [s] (em geral o período nodal).

    Returns:
        (a, e, i) médios [m, -, rad].

    Raises:
        ValueError: se ``t`` não contiver ao menos uma revolução completa.
    """
    t = as_time_array(t)
    n_revs = int(np.floor((t[-1] - t[0]) / period))
    if n_revs < 1:
        raise ValueError("t deve cobrir ao menos uma revolução completa")
    window = t <= t[0] + n_revs * period
    tw = t[window]

    def _avg(x: ArrayLike) -> float:
        xw = np.broadcast_to(np.asarray(x, dtype=float), t.shape)[window]
        return float(np.trapezoid(xw, tw) / (tw[-1] - tw[0]))

    return _avg(elements.a), _avg(elements.e), _avg(elements.i)


def fit_rate(t: ArrayLike, angle: ArrayLike) -> float:
    """Taxa média de um ângulo por ajuste linear de mínimos quadrados [rad/s].

    O ângulo é desembrulhado (``np.unwrap``) antes do ajuste, para que os saltos de 2π de
    valores em [0, 2π) não quebrem a reta.

    Args:
        t: tempos [s], shape (N,).
        angle: ângulo [rad] em cada t, shape (N,).
    """
    t = as_time_array(t)
    slope, _ = np.polyfit(t, np.unwrap(np.asarray(angle, dtype=float)), 1)
    return float(slope)
