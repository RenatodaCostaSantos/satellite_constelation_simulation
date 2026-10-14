"""Propagador de elementos médios com as taxas seculares de J2 (Sem4 §3.5).

Sob J2 (primeira ordem, só efeitos seculares), a, e e i médios ficam constantes e três ângulos
avançam linearmente (Sem4 eqs. 23–25):

    Ω(t) = Ω0 + Ω̇·t,   ω(t) = ω0 + ω̇·t,   ψ(t) = ψ0 + ψ̇·t

O estado em cada instante sai de ``elements_to_state`` vetorizado sobre t. Com J2 = 0 o
resultado é idêntico ao ``KeplerPropagator``.
"""

import math
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import J2, MU, R_EARTH
from satsim.elements import KeplerianElements, elements_to_state
from satsim.propagators.base import Propagator, as_time_array
from satsim.secular import SecularRates, secular_rates


class MeanJ2Propagator(Propagator):
    """Propaga elementos **médios** sob as taxas seculares de J2.

    A saída é a órbita média: a posição difere da osculadora (a de um integrador numérico com
    J2, ``CowellPropagator``) pelas oscilações de curto período, de ~10 km em LEO. Os elementos
    de entrada devem ser médios; tratar elementos osculadores como médios gera um erro de
    movimento médio que cresce ao longo da trajetória.

    Args:
        elements: elementos médios na época (campos escalares; ``psi`` = ψ0).
        epoch: instante t = 0, ``datetime`` tz-aware (UTC recomendado).
        mu: parâmetro gravitacional [m³/s²].
        j2: harmônico zonal J2 (não normalizado).
        re: raio equatorial de referência de J2 [m].

    Raises:
        ValueError: se ``epoch`` não tiver fuso horário (naive), se a ≤ 0 ou e fora de [0, 1).
    """

    def __init__(
        self,
        elements: KeplerianElements,
        epoch: datetime,
        mu: float = MU,
        j2: float = J2,
        re: float = R_EARTH,
    ) -> None:
        if epoch.tzinfo is None or epoch.utcoffset() is None:
            raise ValueError("epoch deve ser um datetime tz-aware (ex.: tzinfo=timezone.utc)")
        self._elements = elements
        self._epoch = epoch
        self._mu = mu
        self._rates = secular_rates(elements.a, elements.e, elements.i, mu=mu, j2=j2, re=re)

    @property
    def epoch(self) -> datetime:
        """Época (t = 0), tz-aware."""
        return self._epoch

    @property
    def elements(self) -> KeplerianElements:
        """Elementos médios na época."""
        return self._elements

    @property
    def rates(self) -> SecularRates:
        """Taxas seculares Ω̇, ω̇, ψ̇ e o movimento médio n [rad/s] (Sem4 eqs. 23–25)."""
        return self._rates

    @property
    def period(self) -> float:
        """Período nodal Tnod = 2π/(ψ̇ + ω̇) [s] (Sem4 eq. 34)."""
        return 2.0 * math.pi / (self._rates.psi_dot + self._rates.argp_dot)

    def elements_at(self, t: ArrayLike) -> KeplerianElements:
        """Elementos médios em t segundos desde ``epoch``.

        Args:
            t: escalar ou sequência (N,) de tempos [s]; pode ser negativo.

        Returns:
            KeplerianElements com todos os campos de shape (N,). Ω, ω e ψ não são reduzidos a
            [0, 2π): crescem continuamente com t.
        """
        t = as_time_array(t)
        el, rates = self._elements, self._rates
        return KeplerianElements(
            a=np.full_like(t, el.a),
            e=np.full_like(t, el.e),
            i=np.full_like(t, el.i),
            raan=el.raan + rates.raan_dot * t,
            argp=el.argp + rates.argp_dot * t,
            psi=el.psi + rates.psi_dot * t,
        )

    def propagate(self, t: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estado ECI em t segundos desde ``epoch``, a partir dos elementos médios.

        Args:
            t: escalar ou sequência (N,) de tempos [s]; pode ser negativo.

        Returns:
            (r, v) em ECI [m, m/s], cada um com shape (N, 3) (N = 1 para t escalar).
        """
        return elements_to_state(self.elements_at(t), self._mu)
