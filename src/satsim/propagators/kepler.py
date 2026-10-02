"""Propagador kepleriano puro (dois corpos, sem perturbações).

Sob gravidade central, cinco elementos ficam constantes e só a anomalia média avança
linearmente com o tempo (Sem4 §1): ψ(t) = ψ0 + n·t, com movimento médio n = √(µ/a³). O estado
em cada instante sai de ``elements_to_state``, vetorizado sobre t.
"""

import math
from dataclasses import replace
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import MU
from satsim.elements import KeplerianElements, elements_to_state
from satsim.propagators.base import Propagator, as_time_array


class KeplerPropagator(Propagator):
    """Propaga uma órbita kepleriana a partir de elementos na época.

    Args:
        elements: elementos clássicos na época (campos escalares; ``psi`` = ψ0).
        epoch: instante t = 0, ``datetime`` tz-aware (UTC recomendado).
        mu: parâmetro gravitacional [m³/s²].

    Raises:
        ValueError: se ``epoch`` não tiver fuso horário (naive), ou se a ≤ 0.
    """

    def __init__(self, elements: KeplerianElements, epoch: datetime, mu: float = MU) -> None:
        if epoch.tzinfo is None or epoch.utcoffset() is None:
            raise ValueError("epoch deve ser um datetime tz-aware (ex.: tzinfo=timezone.utc)")
        if not elements.a > 0.0:
            raise ValueError("semieixo maior deve ser positivo")
        self._elements = elements
        self._epoch = epoch
        self._mu = mu
        self._mean_motion = math.sqrt(mu / elements.a**3)

    @property
    def epoch(self) -> datetime:
        """Época (t = 0), tz-aware."""
        return self._epoch

    @property
    def elements(self) -> KeplerianElements:
        """Elementos na época."""
        return self._elements

    @property
    def mean_motion(self) -> float:
        """Movimento médio n = √(µ/a³) [rad/s]."""
        return self._mean_motion

    @property
    def period(self) -> float:
        """Período kepleriano Tkep = 2π/n [s] (Sem4, tabela de símbolos)."""
        return 2.0 * math.pi / self._mean_motion

    def propagate(self, t: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estado ECI em t segundos desde ``epoch``.

        Args:
            t: escalar ou sequência (N,) de tempos [s]; pode ser negativo.

        Returns:
            (r, v) em ECI [m, m/s], cada um com shape (N, 3) (N = 1 para t escalar).
        """
        t = as_time_array(t)
        psi = self._elements.psi + self._mean_motion * t
        return elements_to_state(replace(self._elements, psi=psi), self._mu)
