"""Determinação de órbita mínima: ajusta o estado inicial do Cowell a efemérides (card S4-06).

É o que operadores fazem com as medidas de radar, com o SGP4 no papel dos radares: em vez de
partir do estado do SGP4 em t = 0 (método da Semana 3), procura-se o estado inicial (r0, v0) e,
opcionalmente, o C_D·A/m que fazem a trajetória do Cowell passar o mais perto possível das
posições de referência numa janela (72 h por padrão):

    min  Σ_k |r_Cowell(t_k; r0, v0, C_D·A/m) − r_ref(t_k)|²     t_k = 0, step_s, …, window_s

Mínimos quadrados não lineares (``scipy.optimize.least_squares``, Jacobiano por diferenças
finitas) sobre x = (Δr0 [km], Δv0 [m/s], C_D·A/m [1e-3 m²/kg]), a partir do estado da referência
em t = 0. As unidades deixam todos os parâmetros com ordem de grandeza 1, o que mantém o
Jacobiano bem condicionado e os passos das diferenças finitas comparáveis.

Por que o ajuste resolve o "piso" do estado em t = 0
----------------------------------------------------
O estado do SGP4 numa época tem erro de ~1 km (é uma teoria analítica truncada). Partindo dele,
o Cowell herda esse erro, e um erro de semieixo vira deriva ao longo da órbita que cresce
linearmente com o tempo: há um piso que nenhuma força a mais remove. O ajuste distribui o erro
pela janela inteira, e o que sobra é só a diferença de modelo entre o Cowell e o SGP4.

Cuidado com bounds: um parâmetro parado no limite trava o otimizador (no experimento do prompt,
o C_D·A/m preso em 0 não saía do lugar). Sem arrasto (``fit_drag=False``) ajustam-se só os 6
parâmetros do estado; com arrasto, o C_D·A/m parte de ``cd_a_m0`` > 0, longe do limite inferior.
O C_D·A/m ajustado é **efetivo**: absorve a atmosfera do dia e o que falta ao modelo de forças.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares

from satsim.forces import AccelerationModel, DragAcceleration
from satsim.forces.drag import AltitudeModel
from satsim.propagators.base import Propagator
from satsim.propagators.cowell import CowellPropagator

_KM = 1e3
_CDAM_UNIT = 1e-3  # unidade interna de C_D·A/m no vetor de parâmetros [m²/kg]


@dataclass(frozen=True)
class OrbitFitResult:
    """Resultado de ``fit_initial_state``.

    r0, v0: estado inicial ajustado em ECI [m, m/s] na época da referência, shape (3,).
    cd_a_m: C_D·A/m ajustado [m²/kg] (efetivo), ou ``None`` sem arrasto.
    rms_km: RMS do erro de posição 3D na janela, no ajuste final [km].
    rms_initial_km: RMS partindo do chute inicial (antes do ajuste) [km].
    history_km: RMS de cada avaliação da função (inclui as do Jacobiano) [km].
    nfev: número de avaliações (cada uma é uma propagação do Cowell na janela).
    success: convergência reportada por ``least_squares``.
    epoch: época (t = 0) da referência.
    models: forças usadas (já com o arrasto ajustado, se houver).
    """

    r0: NDArray[np.float64]
    v0: NDArray[np.float64]
    cd_a_m: float | None
    rms_km: float
    rms_initial_km: float
    history_km: tuple[float, ...]
    nfev: int
    success: bool
    epoch: datetime
    models: tuple[AccelerationModel, ...]

    def propagator(self, **cowell_kwargs) -> CowellPropagator:
        """Cowell com o estado e as forças ajustados (``cowell_kwargs`` vai ao construtor)."""
        return CowellPropagator((self.r0, self.v0), self.epoch, self.models, **cowell_kwargs)


def fit_initial_state(
    reference: Propagator,
    models: Sequence[AccelerationModel],
    window_s: float = 72 * 3600.0,
    step_s: float = 600.0,
    fit_drag: bool = True,
    cd: float = 2.2,
    cd_a_m0: float = 0.002,
    altitude_model: AltitudeModel = "geodetic",
    rtol: float = 1e-10,
    atol: float = 1e-3,
    initial: tuple[ArrayLike, ArrayLike] | None = None,
) -> OrbitFitResult:
    """Ajusta (r0, v0) e, opcionalmente, C_D·A/m do Cowell às posições de ``reference``.

    Args:
        reference: efemérides de referência (ex.: ``SGP4Propagator``); define a época.
        models: forças do Cowell **sem** o arrasto (ele é acrescentado se ``fit_drag``).
        window_s: janela do ajuste [s].
        step_s: espaçamento das posições de referência [s].
        fit_drag: ajustar também C_D·A/m (acrescenta um ``DragAcceleration``).
        cd: C_D para repartir o produto C_D·A/m (só o produto importa).
        cd_a_m0: chute inicial de C_D·A/m [m²/kg], > 0 (longe do limite inferior).
        altitude_model: altitude do arrasto.
        rtol, atol: tolerâncias do DOP853 em cada propagação (as do experimento do A.8).
        initial: chute inicial (r0, v0) [m, m/s]; padrão: o estado de ``reference`` em t = 0.

    Returns:
        ``OrbitFitResult``.

    Raises:
        ValueError: ``models`` já contém arrasto com ``fit_drag=True``, ou ``cd_a_m0 <= 0``.
    """
    models = tuple(models)
    if fit_drag and any(isinstance(m, DragAcceleration) for m in models):
        raise ValueError("com fit_drag=True o arrasto é acrescentado aqui: tire-o de models")
    if fit_drag and cd_a_m0 <= 0.0:
        raise ValueError("cd_a_m0 deve ser positivo (um parâmetro no limite trava o ajuste)")
    epoch = reference.epoch
    t = np.arange(0.0, window_s + 0.5 * step_s, step_s)
    r_ref, _ = reference.propagate(t)
    if initial is None:
        r_start, v_start = (x[0] for x in reference.propagate(0.0))
    else:
        r_start, v_start = (np.asarray(x, dtype=float) for x in initial)
    history: list[float] = []

    def unpack(x: NDArray[np.float64]):
        r0 = r_start + x[0:3] * _KM
        v0 = v_start + x[3:6]
        cd_a_m = float(x[6]) * _CDAM_UNIT if fit_drag else None
        return r0, v0, cd_a_m

    def build_models(cd_a_m: float | None) -> tuple[AccelerationModel, ...]:
        if cd_a_m is None:
            return models
        drag = DragAcceleration(cd=cd, area_over_mass=cd_a_m / cd, altitude_model=altitude_model)
        return (*models, drag)

    def residuals(x: NDArray[np.float64]) -> NDArray[np.float64]:
        r0, v0, cd_a_m = unpack(x)
        prop = CowellPropagator((r0, v0), epoch, build_models(cd_a_m), rtol=rtol, atol=atol)
        r, _ = prop.propagate(t)
        res = ((r - r_ref) / _KM).ravel()
        history.append(float(np.sqrt(np.mean(np.sum(((r - r_ref) / _KM) ** 2, axis=1)))))
        return res

    x0 = np.zeros(7 if fit_drag else 6)
    lower = np.full(x0.size, -np.inf)
    if fit_drag:
        x0[6] = cd_a_m0 / _CDAM_UNIT
        lower[6] = 0.0
    sol = least_squares(
        residuals,
        x0,
        bounds=(lower, np.full(x0.size, np.inf)),
        # passo relativo das diferenças finitas: 0,1 m em r0, 0,1 mm/s em v0. Bem acima do ruído
        # da integração (que muda a sequência de passos a cada estado), bem abaixo da escala em
        # que o problema deixa de ser linear. Com 1e-6 o ajuste parava a ≈ 0,3 m do estado certo.
        diff_step=1e-4,
        xtol=1e-8,
        ftol=1e-8,
        gtol=1e-8,
    )
    r0, v0, cd_a_m = unpack(sol.x)
    rms = float(np.sqrt(np.mean(np.sum(sol.fun.reshape(-1, 3) ** 2, axis=1))))
    return OrbitFitResult(
        r0=r0,
        v0=v0,
        cd_a_m=cd_a_m,
        rms_km=rms,
        rms_initial_km=history[0],
        history_km=tuple(history),
        nfev=len(history),
        success=bool(sol.success),
        epoch=epoch,
        models=build_models(cd_a_m),
    )
