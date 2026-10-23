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

Reentrada (S3-03)
-----------------
Com arrasto (``DragAcceleration`` entre os modelos), a integração para num evento terminal quando
a altitude cai abaixo de ``stop_altitude_m`` (padrão 100 km, base da tabela atmosférica). O
instante fica em ``reentry_time``; pedidos de instantes posteriores devolvem ``NaN``, com aviso.
"""

import logging
from collections.abc import Sequence
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import OdeSolution, solve_ivp

from satsim.constants import MU
from satsim.elements import KeplerianElements, elements_to_state
from satsim.forces import AccelerationModel, DragAcceleration, J2Acceleration, altitude
from satsim.propagators.base import Propagator, as_time_array

logger = logging.getLogger(__name__)

_DEFAULT_MODELS: tuple[AccelerationModel, ...] = (J2Acceleration(),)
DEFAULT_STOP_ALTITUDE_M = 100e3


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
        stop_altitude_m: altitude [m] do evento terminal de reentrada, ativo só quando há um
            ``DragAcceleration`` em ``models`` (medida com o ``altitude_model`` dele); ``None``
            desliga o evento.

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
        stop_altitude_m: float | None = DEFAULT_STOP_ALTITUDE_M,
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
        drag = next((m for m in self._models if isinstance(m, DragAcceleration)), None)
        self._stop_altitude_m = stop_altitude_m if drag is not None else None
        self._stop_altitude_model = drag.altitude_model if drag is not None else "geodetic"
        self._reentry_time: float | None = None

    @classmethod
    def from_state(
        cls,
        r: ArrayLike,
        v: ArrayLike,
        epoch: datetime,
        models: Sequence[AccelerationModel] = _DEFAULT_MODELS,
        **kwargs,
    ) -> "CowellPropagator":
        """Cowell a partir do estado (r, v) em ECI na época, ex.: o do SGP4 em t = 0 (S3-03).

        Equivale a ``CowellPropagator((r, v), epoch, models, **kwargs)``; ``kwargs`` vai para o
        construtor (``mu``, ``rtol``, ``atol``, ``stop_altitude_m``).
        """
        return cls((r, v), epoch, models, **kwargs)

    @property
    def reentry_time(self) -> float | None:
        """Instante [s desde a época] do evento de reentrada, se já integrado; senão ``None``."""
        return self._reentry_time

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
        if t_max <= t_start or self._reentry_time is not None:
            return
        events = None
        if self._stop_altitude_m is not None:
            stop_m, model = self._stop_altitude_m, self._stop_altitude_model

            def reentry(t: float, y: NDArray[np.float64]) -> float:
                return float(altitude(y[:3], model)) - stop_m

            reentry.terminal = True  # type: ignore[attr-defined]
            reentry.direction = -1  # type: ignore[attr-defined]
            events = [reentry]
        sol = solve_ivp(
            self._rhs,
            (t_start, t_max),
            self._y_end,
            method="DOP853",
            rtol=self._rtol,
            atol=self._atol,
            dense_output=True,
            events=events,
        )
        if not sol.success:
            raise RuntimeError(f"integração falhou: {sol.message}")
        self._segments.append(sol.sol)
        self._ends.append(float(sol.t[-1]))  # = t_max, ou o instante do evento
        self._y_end = sol.y[:, -1].copy()
        if sol.status == 1:  # evento terminal: reentrada
            self._reentry_time = float(sol.t_events[0][0])
            logger.warning(
                "reentrada: altitude < %.0f km em t = %.0f s (%.2f h) desde a época",
                self._stop_altitude_m / 1e3,
                self._reentry_time,
                self._reentry_time / 3600.0,
            )

    def propagate(self, t: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estado ECI em t segundos desde ``epoch``.

        Args:
            t: escalar ou sequência (N,) de tempos [s], todos ≥ 0, em qualquer ordem.

        Returns:
            (r, v) em ECI [m, m/s], cada um com shape (N, 3) (N = 1 para t escalar); linhas
            ``NaN`` para instantes depois da reentrada (com aviso).

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
            after = np.zeros(t.shape, dtype=bool)
            if self._reentry_time is not None:
                after = t > self._reentry_time
                y[after] = np.nan
                if np.any(after):
                    logger.warning(
                        "%d instantes depois da reentrada (t > %.0f s) devolvidos como NaN",
                        int(after.sum()),
                        self._reentry_time,
                    )
            todo = ~at_epoch & ~after
            # Trecho k cobre [ends[k-1], ends[k]]; searchsorted(left) põe t = ends[k] no trecho k.
            seg_idx = np.searchsorted(np.asarray(self._ends), t, side="left")
            for k in np.unique(seg_idx[todo]):
                sel = (seg_idx == k) & todo
                y[sel] = self._segments[k](t[sel]).T
        return y[:, :3], y[:, 3:]
