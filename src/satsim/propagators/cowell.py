"""Propagador numérico pelo método de Cowell: integra r̈ = −µr/r³ + Σ a_k(t, r, v).

O estado é integrado diretamente em coordenadas cartesianas ECI com ``scipy.integrate.solve_ivp``
(DOP853, Runge-Kutta de ordem 8 com passo adaptativo e saída densa). As perturbações entram como
modelos de aceleração (``satsim.forces.AccelerationModel``): J2 nesta semana, arrasto na
Semana 3, sem mudar o propagador.

Integração preguiçosa com cache
-------------------------------
Nada é integrado no construtor. Cada chamada a ``propagate(t)`` integra apenas de onde a última
integração parou até o maior t pedido e guarda o trecho (com saída densa); chamadas seguintes
reutilizam os trechos já integrados. Assim, pedir os mesmos instantes de novo, ou fora de ordem,
devolve sempre os mesmos valores. Só tempos t ≥ 0 são aceitos (sem integração para trás).
"""

from collections.abc import Sequence
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import OdeSolution, solve_ivp

from satsim.constants import MU
from satsim.elements import KeplerianElements, elements_to_state
from satsim.forces import AccelerationModel, J2Acceleration
from satsim.propagators.base import Propagator, as_time_array

_DEFAULT_MODELS: tuple[AccelerationModel, ...] = (J2Acceleration(),)


class CowellPropagator(Propagator):
    """Integra numericamente a equação de movimento com perturbações (método de Cowell).

    **Elementos osculadores ≠ médios.** A condição inicial são elementos *osculadores*: o
    estado (r, v) exato na época. Elementos nominais de projeto costumam ser *médios* (os do
    ``MeanJ2Propagator``); usados aqui como osculadores, o semieixo médio resultante fica alguns
    km fora do nominal (ex.: +9,5 km na SSO de projeto), o movimento médio muda em ~1e-3 e a
    posição se afasta da do ``MeanJ2Propagator`` nominal em até ~1 300 km em 24 h. Não é erro de
    integração: é a diferença entre as duas descrições. A conversão média → osculador
    (Brouwer/Kozai) não é implementada; ``satsim.analysis.mean_elements`` mede os elementos
    médios de uma trajetória.

    Args:
        initial: elementos osculadores na época (``KeplerianElements`` com campos escalares) ou
            estado ``(r0, v0)`` em ECI [m, m/s], cada um com shape (3,).
        epoch: instante t = 0, ``datetime`` tz-aware (UTC recomendado).
        models: acelerações perturbadoras; ``()`` integra o problema de dois corpos puro.
        mu: parâmetro gravitacional do termo central [m³/s²].
        rtol: tolerância relativa do integrador.
        atol: tolerância absoluta do integrador (m nas posições, m/s nas velocidades). Em LEO
            ela domina sobre rtol·|y|; com 1e-3 a energia deriva ~1e-10/dia (1,1e-9 em 10 dias),
            com 1e-6 fica em ~4e-11 em 10 dias, com ~50% mais passos.

    Raises:
        ValueError: se ``epoch`` não tiver fuso horário (naive) ou o estado inicial não tiver
            shape (3,).
    """

    def __init__(
        self,
        initial: KeplerianElements | tuple[ArrayLike, ArrayLike],
        epoch: datetime,
        models: Sequence[AccelerationModel] = _DEFAULT_MODELS,
        mu: float = MU,
        rtol: float = 1e-11,
        atol: float = 1e-6,
    ) -> None:
        if epoch.tzinfo is None or epoch.utcoffset() is None:
            raise ValueError("epoch deve ser um datetime tz-aware (ex.: tzinfo=timezone.utc)")
        if isinstance(initial, KeplerianElements):
            r0, v0 = elements_to_state(initial, mu)
            r0, v0 = r0[0], v0[0]
        else:
            r0, v0 = (np.asarray(x, dtype=float) for x in initial)
            if r0.shape != (3,) or v0.shape != (3,):
                raise ValueError("o estado inicial (r0, v0) deve ter shape (3,) cada")
        self._epoch = epoch
        self._models = tuple(models)
        self._mu = mu
        self._rtol = rtol
        self._atol = atol
        self._y0 = np.concatenate([r0, v0])
        # Cache: trechos de saída densa contíguos [0, ends[0]], [ends[0], ends[1]], ...
        self._segments: list[OdeSolution] = []
        self._ends: list[float] = []
        self._y_end = self._y0.copy()

    @property
    def epoch(self) -> datetime:
        """Época (t = 0), tz-aware."""
        return self._epoch

    @property
    def models(self) -> tuple[AccelerationModel, ...]:
        """Modelos de aceleração perturbadora."""
        return self._models

    @property
    def initial_state(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estado (r0, v0) em ECI na época, cada um com shape (3,)."""
        return self._y0[:3].copy(), self._y0[3:].copy()

    def _rhs(self, t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        r, v = y[:3], y[3:]
        r_norm = np.sqrt(r @ r)
        acc = -self._mu / r_norm**3 * r
        for model in self._models:
            acc = acc + model(t, r, v)
        return np.concatenate([v, acc])

    def _extend(self, t_max: float) -> None:
        """Integra do fim do cache até ``t_max`` e guarda o trecho."""
        t_start = self._ends[-1] if self._ends else 0.0
        if t_max <= t_start:
            return
        sol = solve_ivp(
            self._rhs,
            (t_start, t_max),
            self._y_end,
            method="DOP853",
            rtol=self._rtol,
            atol=self._atol,
            dense_output=True,
        )
        if not sol.success:
            raise RuntimeError(f"integração falhou: {sol.message}")
        self._segments.append(sol.sol)
        self._ends.append(t_max)
        self._y_end = sol.y[:, -1].copy()

    def propagate(self, t: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estado ECI em t segundos desde ``epoch``.

        Args:
            t: escalar ou sequência (N,) de tempos [s], todos ≥ 0, em qualquer ordem.

        Returns:
            (r, v) em ECI [m, m/s], cada um com shape (N, 3) (N = 1 para t escalar).

        Raises:
            ValueError: se algum t < 0 (ou não finito).
        """
        t = as_time_array(t)
        if not np.all(np.isfinite(t)) or np.any(t < 0.0):
            raise ValueError("CowellPropagator só aceita t >= 0 (finito)")
        y = np.empty((t.size, 6))
        at_epoch = t == 0.0
        y[at_epoch] = self._y0
        if not np.all(at_epoch):
            self._extend(float(t.max()))
            # Trecho k cobre [ends[k-1], ends[k]]; searchsorted(left) põe t = ends[k] no trecho k.
            seg_idx = np.searchsorted(np.asarray(self._ends), t, side="left")
            for k in np.unique(seg_idx[~at_epoch]):
                sel = (seg_idx == k) & ~at_epoch
                y[sel] = self._segments[k](t[sel]).T
        return y[:, :3], y[:, 3:]
