"""Análise de trajetórias: elementos osculadores, elementos médios e taxas ajustadas.

Elementos osculadores são os elementos keplerianos do estado instantâneo (r, v): a órbita de dois
corpos que o satélite seguiria se as perturbações sumissem naquele instante. Sob J2 eles oscilam
a cada revolução (curto período, ~10 km no semieixo em LEO) em torno de valores médios que só
derivam secularmente (Sem4 §3). Este módulo mede essas duas componentes numa trajetória
propagada, para comparar com a teoria secular (``satsim.secular``).

Também ajusta o coeficiente balístico C_D·A/m do Cowell a uma trajetória de referência (S3-03).
"""

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import minimize_scalar

from satsim.elements import KeplerianElements, state_to_elements
from satsim.forces import DragAcceleration, J2Acceleration
from satsim.forces.drag import AltitudeModel
from satsim.propagators.base import Propagator, as_time_array
from satsim.propagators.cowell import CowellPropagator


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


def cowell_j2_drag(
    reference: Propagator,
    cd_a_m: float,
    cd: float = 2.2,
    altitude_model: AltitudeModel = "geodetic",
    rtol: float = 1e-10,
    atol: float = 1e-3,
) -> CowellPropagator:
    """Cowell (J2 + arrasto com C_D·A/m = ``cd_a_m``) partindo do estado de ``reference`` em t = 0.

    O produto é repartido como C_D = ``cd`` e A/m = ``cd_a_m``/``cd``; só o produto importa.
    Com ``cd_a_m = 0`` o arrasto é nulo (J2 apenas).
    """
    r0, v0 = reference.propagate(0.0)
    drag = DragAcceleration(cd=cd, area_over_mass=cd_a_m / cd, altitude_model=altitude_model)
    return CowellPropagator.from_state(
        r0[0], v0[0], reference.epoch, (J2Acceleration(), drag), rtol=rtol, atol=atol
    )


def position_error(propagator: Propagator, reference: Propagator, t: ArrayLike) -> np.ndarray:
    """Erro de posição 3D |r − r_ref| [m] nos instantes t [s desde a época de ``reference``].

    ``propagator`` pode ter outra época: usa ``propagate_utc``.
    """
    t = as_time_array(t)
    r, _ = propagator.propagate_utc(reference.epoch, t)
    r_ref, _ = reference.propagate(t)
    return np.linalg.norm(r - r_ref, axis=1)


def fit_ballistic_coefficient(
    reference: Propagator,
    window_s: float = 72 * 3600.0,
    step_s: float = 600.0,
    bounds: tuple[float, float] = (0.0, 0.03),
    cd: float = 2.2,
    altitude_model: AltitudeModel = "geodetic",
    xatol: float = 1e-6,
) -> float:
    """C_D·A/m [m²/kg] que minimiza o RMS do erro de posição do Cowell contra ``reference``.

    O Cowell (J2 + arrasto, ``cowell_j2_drag``) parte do estado de ``reference`` em t = 0 e é
    comparado a cada ``step_s`` em [0, ``window_s``]. A busca é unidimensional e limitada
    (``scipy.optimize.minimize_scalar``, método ``bounded``) em ``bounds``. É uma versão mínima de
    determinação de órbita: em vez de radares, o SGP4 faz o papel das observações, e só o arrasto
    é ajustado (o estado inicial é o do SGP4).

    O valor ajustado absorve o que o modelo não tem (atmosfera real do dia, atitude da ISS, J3/J4
    e Sol/Lua do SGP4): é um parâmetro efetivo, não o C_D·A/m físico.

    Args:
        reference: propagador de referência (ex.: ``SGP4Propagator``).
        window_s: janela do ajuste [s] (padrão 72 h).
        step_s: espaçamento das amostras do erro [s].
        bounds: intervalo de busca de C_D·A/m [m²/kg].
        cd: C_D usado para repartir o produto (não afeta o resultado).
        altitude_model: altitude do arrasto (``"geodetic"`` ou ``"spherical"``).
        xatol: tolerância absoluta em C_D·A/m.
    """
    t = np.arange(0.0, window_s + 0.5 * step_s, step_s)

    def rms(cd_a_m: float) -> float:
        cowell = cowell_j2_drag(reference, cd_a_m, cd=cd, altitude_model=altitude_model)
        err = position_error(cowell, reference, t)
        return float(np.sqrt(np.mean(err**2)))

    result = minimize_scalar(rms, bounds=bounds, method="bounded", options={"xatol": xatol})
    return float(result.x)
