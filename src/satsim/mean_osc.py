"""Conversão médio ↔ osculador: a ponte entre o projeto e a simulação (card S4-04).

Os elementos de projeto (os da SSO 91/6) são **médios**: descrevem a órbita sem as oscilações de
curto período. O Cowell precisa de um estado **osculador**: o (r, v) exato num instante. Entregar
os médios ao Cowell como se fossem osculadores erra o semieixo médio em ≈ +9,5 km na órbita de
projeto, o perigeu não congela e o LTAN deriva ≈ −12 min/ano (Apêndice A.6 do prompt da Semana 4;
relatório da S2-06).

Este módulo resolve o problema numericamente, para **qualquer conjunto de forças**: a média de um
estado é a média dos elementos osculadores do Cowell ao longo de uma revolução nodal, e a
conversão inversa é um ponto fixo sobre essa média. Teorias analíticas (Brouwer–Lyddane, Kozai)
ficam fora do escopo; podem servir de chute inicial do ponto fixo no futuro.

Quais forças entram na conversão (convenção do projeto, S4-05)
--------------------------------------------------------------
As funções aceitam qualquer conjunto de forças, mas a convenção é converter **só com a gravidade
zonal** (``conversion_models``) e acrescentar arrasto, SRP e Sol/Lua apenas na propagação. É a
mesma divisão do SGP4, em que a relação médio ↔ osculador só tem termos gravitacionais e o
arrasto entra como deriva secular dos elementos médios. Motivos: (1) o J2 responde pelos ≈ 9,5 km
da conversão, e as demais forças mudam o resultado em metros; (2) com forças conservativas e
simétricas o elemento médio é constante (a, e, i) ou uma taxa constante (Ω, ω), o que o torna útil
para projeto; com arrasto o semieixo cai dentro da janela (≈ 3,9 m por volta a 510 km) e a média
passaria a depender de onde a janela começa.

Por que a média do vetor excentricidade
---------------------------------------
Numa órbita quase circular (e ≈ 1e-3), as oscilações de curto período de J2 no vetor
excentricidade (e·cos ω, e·sin ω) têm amplitude ~J2·(R⊕/a)² ~ 1e-3, do tamanho do próprio e: ao
longo de uma volta, a ponta do vetor descreve uma pequena curva em torno do valor médio, às vezes
passando perto da origem. O ω osculador então gira dezenas de graus por volta e o e osculador
oscila de ~0 ao dobro. Médias separadas de e e de ω ficam enviesadas (a média de um módulo é
maior que o módulo da média; a média de um ângulo que dá voltas não tem sentido). A média do
vetor é linear e bem definida mesmo para e → 0; e e ω médios saem dela:
e = |⟨vetor⟩|, ω = atan2(⟨e·sin ω⟩, ⟨e·cos ω⟩).

Por que a 1ª iteração corrige Δa ≈ 9,5 km
-----------------------------------------
O semieixo osculador oscila em torno do médio com o termo de curto período de J2 (Kozai, e ≈ 0):

    a_osc − a_médio ≈ (3/2)·J2·(R⊕²/a)·sin²i·cos 2u

Na órbita de projeto, J2·R⊕²/a = 6,39 km e sin²i = 0,983, logo a amplitude é ≈ 9,4 km. Os
elementos de projeto começam em ω = 90°, ψ = 0, isto é, em u = ω + ν = 90° (o ponto mais ao norte),
onde cos 2u = −1: ali o osculador está ≈ 9,4 km **abaixo** do médio. Usar o médio como osculador
põe o satélite 9,4 km acima do que deveria naquele ponto, e o semieixo médio resultante fica ≈ 9,5
km alto: é o +9,5 km do relatório da S2-06. A 1ª iteração mede exatamente esse desvio e o
subtrai; as seguintes só corrigem o que a não linearidade e J3/J4 deixam (74 m, 0,25 m, …).
"""

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.analysis import osculating_elements
from satsim.elements import KeplerianElements, mean_to_true, state_to_elements, true_to_mean
from satsim.forces import AccelerationModel, J2Acceleration, ZonalAcceleration
from satsim.propagators.base import Propagator
from satsim.propagators.cowell import CowellPropagator
from satsim.secular import nodal_period

logger = logging.getLogger(__name__)

FloatOrArray = float | NDArray[np.float64]

DEFAULT_TOL: tuple[float, float, float] = (1e-3, 1e-9, 1e-9)
"""Tolerâncias do ponto fixo: |Δa| [m], |Δ(vetor e)| (cada componente) [-], |Δi| [rad]."""


@dataclass(frozen=True)
class MeanElements:
    """Elementos médios por janela: a [m], e [-], i [rad], argp [rad] (escalares ou arrays).

    e e argp vêm da média do vetor excentricidade (ver docstring do módulo).
    """

    a: FloatOrArray
    e: FloatOrArray
    i: FloatOrArray
    argp: FloatOrArray

    @property
    def ex(self) -> FloatOrArray:
        """e·cos ω médio."""
        return self.e * np.cos(self.argp)

    @property
    def ey(self) -> FloatOrArray:
        """e·sin ω médio."""
        return self.e * np.sin(self.argp)


@dataclass(frozen=True)
class MeanToOscResult:
    """Resultado de ``mean_to_osculating``.

    elements: elementos osculadores na época (``KeplerianElements`` escalar).
    residuals: histórico, shape (k, 4): |Δa| [m], |Δ(e cos ω)|, |Δ(e sin ω)|, |Δi| [rad] entre o
        alvo médio e a média medida, uma linha por iteração (antes da correção daquela iteração).
    converged: tolerâncias atingidas dentro de ``iters``.
    period: período nodal usado nas janelas [s].
    """

    elements: KeplerianElements
    residuals: NDArray[np.float64]
    converged: bool
    period: float

    @property
    def iterations(self) -> int:
        """Número de propagações feitas (linhas de ``residuals``)."""
        return len(self.residuals)


def conversion_models(models: Sequence[AccelerationModel]) -> tuple[AccelerationModel, ...]:
    """Só a gravidade zonal de ``models`` (``J2Acceleration``, ``ZonalAcceleration``).

    Convenção do projeto para a conversão médio ↔ osculador (docstring do módulo): arrasto, SRP e
    Sol/Lua ficam de fora da conversão e entram só na propagação.
    """
    return tuple(m for m in models if isinstance(m, (J2Acceleration, ZonalAcceleration)))


def _scalar_if_0d(x: NDArray[np.float64]) -> FloatOrArray:
    return float(x) if np.ndim(x) == 0 else x


def window_mean_elements(
    propagator: Propagator, t0: ArrayLike, period: float, step_s: float = 10.0
) -> MeanElements:
    """Médias trapezoidais de a, i e do vetor excentricidade numa revolução nodal a partir de t0.

    Cada janela [t0, t0 + period] é amostrada com passo ≤ ``step_s`` (as extremidades caem
    exatamente nas bordas) e a média é ∫x dt / period pela regra do trapézio. A média do vetor
    (e·cos ω, e·sin ω) é bem definida para e → 0 (docstring do módulo).

    Args:
        propagator: qualquer ``Propagator`` (em geral um ``CowellPropagator``).
        t0: início da janela [s desde a época], escalar ou (K,) para K janelas.
        period: duração da janela [s] (o período nodal, ``secular.nodal_period``).
        step_s: passo máximo de amostragem [s].

    Returns:
        ``MeanElements`` com campos escalares (t0 escalar) ou (K,).
    """
    t0 = np.asarray(t0, dtype=float)
    n = int(np.ceil(period / step_s)) + 1
    tau = np.linspace(0.0, period, n)
    t = (np.atleast_1d(t0)[:, None] + tau[None, :]).ravel()
    el = osculating_elements(propagator, t)

    def avg(x: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.trapezoid(np.reshape(x, (-1, n)), tau, axis=1) / period

    a = avg(el.a)
    i = avg(el.i)
    ex = avg(el.e * np.cos(el.argp))
    ey = avg(el.e * np.sin(el.argp))
    shape = t0.shape

    def out(x: NDArray[np.float64]) -> FloatOrArray:
        return _scalar_if_0d(np.reshape(x, shape))

    return MeanElements(
        a=out(a),
        e=out(np.hypot(ex, ey)),
        i=out(i),
        argp=out(np.mod(np.arctan2(ey, ex), 2.0 * np.pi)),
    )


def osculating_to_mean(
    elements_or_state: KeplerianElements | tuple[ArrayLike, ArrayLike],
    epoch: datetime,
    models: Sequence[AccelerationModel],
    step_s: float = 10.0,
    **cowell_kwargs,
) -> MeanElements:
    """Elementos médios de um estado osculador: média de uma revolução nodal do Cowell.

    O período nodal da janela vem dos próprios elementos médios (``secular.nodal_period``),
    recalculado até estabilizar (≤ 1e-6 s, no máximo 5 passes): um período 12 s errado (o de um
    semieixo 9,5 km fora) enviesaria a média em dezenas de metros.

    Args:
        elements_or_state: elementos osculadores ou estado (r0, v0) ECI na época.
        epoch: época, ``datetime`` tz-aware.
        models: forças do Cowell; por convenção só a gravidade zonal (``conversion_models``).
        step_s: passo da média [s].
        **cowell_kwargs: repassados ao ``CowellPropagator`` (``rtol``, ``atol``...).
    """
    prop = CowellPropagator(elements_or_state, epoch, models, **cowell_kwargs)
    if isinstance(elements_or_state, KeplerianElements):
        start = elements_or_state
    else:
        r0, v0 = prop.initial_state
        start = state_to_elements(r0, v0)
    period = float(nodal_period(start.a, start.e, start.i))
    for _ in range(5):
        mean = window_mean_elements(prop, 0.0, period, step_s)
        new_period = float(nodal_period(mean.a, mean.e, mean.i))
        converged = abs(new_period - period) <= 1e-6
        period = new_period
        if converged:
            break
    return mean


def _osculating_from_vector(x: NDArray[np.float64], raan: float, u0: float) -> KeplerianElements:
    """x = (a, e·cos ω, e·sin ω, i) + Ω + argumento de latitude u0 → elementos osculadores."""
    a, ex, ey, i = x
    e = float(np.hypot(ex, ey))
    argp = float(np.mod(np.arctan2(ey, ex), 2.0 * np.pi))
    nu = u0 - argp
    psi = float(np.mod(true_to_mean(nu, e), 2.0 * np.pi))
    return KeplerianElements(a=float(a), e=e, i=float(i), raan=raan, argp=argp, psi=psi)


def mean_to_osculating(
    mean: KeplerianElements,
    epoch: datetime,
    models: Sequence[AccelerationModel],
    iters: int = 6,
    tol: tuple[float, float, float] = DEFAULT_TOL,
    step_s: float = 10.0,
    **cowell_kwargs,
) -> MeanToOscResult:
    """Estado osculador cujo movimento no Cowell(models) tem as médias pedidas (ponto fixo).

    Sobre x = (a, e·cos ω, e·sin ω, i) (card S4-04):

        x_osc ← x_médio-alvo
        repita: propague o Cowell(models) a partir de x_osc por 1 período nodal (passo step_s)
                m ← window_mean_elements(...)
                x_osc ← x_osc + (x_médio-alvo − m)
        até |Δa| < tol[0], |Δe_vet| < tol[1] (cada componente), |Δi| < tol[2] (ou iters esgotado)

    Ω vem do alvo. A fase inicial preserva o argumento de latitude do alvo, u0 = ω + ν: numa órbita
    quase circular ω e ν isolados mudam muito entre médio e osculador (o ω osculador da órbita de
    projeto sai 270°, contra 90° médio), mas a posição do satélite na órbita é u, e é u que deve
    ser preservado. O período nodal é o do alvo (``secular.nodal_period``, só J2): a diferença
    para o período real com J3/J4 (~5 ppm) é a mesma em todas as iterações e não muda o resultado.

    Args:
        mean: elementos médios-alvo (a, e, i, Ω, ω, ψ), ``KeplerianElements`` escalar.
        epoch: época, ``datetime`` tz-aware.
        models: forças do Cowell; por convenção só a gravidade zonal (``conversion_models``).
        iters: número máximo de propagações.
        tol: tolerâncias (|Δa| [m], |Δe_vet| [-], |Δi| [rad]).
        step_s: passo da média [s].
        **cowell_kwargs: repassados ao ``CowellPropagator``.

    Returns:
        ``MeanToOscResult`` com os elementos osculadores e o histórico de resíduos. Se ``iters``
        se esgotar sem atingir ``tol``, ``converged`` é ``False`` e um aviso é registrado.
    """
    target = np.array(
        [mean.a, mean.e * np.cos(mean.argp), mean.e * np.sin(mean.argp), mean.i], dtype=float
    )
    u0 = float(mean.argp + mean_to_true(mean.psi, mean.e))
    raan = float(mean.raan)
    period = float(nodal_period(mean.a, mean.e, mean.i))
    tol_vec = np.array([tol[0], tol[1], tol[1], tol[2]])

    x = target.copy()
    residuals = []
    converged = False
    for _ in range(iters):
        osc = _osculating_from_vector(x, raan, u0)
        prop = CowellPropagator(osc, epoch, models, **cowell_kwargs)
        m = window_mean_elements(prop, 0.0, period, step_s)
        delta = target - np.array([m.a, m.ex, m.ey, m.i])
        residuals.append(np.abs(delta))
        if np.all(np.abs(delta) < tol_vec):
            converged = True
            break
        x = x + delta
    if not converged:
        logger.warning(
            "mean_to_osculating não convergiu em %d iterações (último resíduo %s)",
            iters,
            residuals[-1],
        )
    return MeanToOscResult(
        elements=_osculating_from_vector(x, raan, u0),
        residuals=np.array(residuals),
        converged=converged,
        period=period,
    )
