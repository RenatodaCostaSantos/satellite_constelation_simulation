"""Interface base de todos os propagadores (ADR 0001, seção Interfaces).

Todo propagador implementa ``propagate(t) → (r, v)`` vetorizado: ``t`` em segundos desde a
``epoch`` do propagador, estado em ECI e SI, sempre com shape (N, 3). Kepler, J2 médio, Cowell e
SGP4 são implementações intercambiáveis desta ABC.
"""

from abc import ABC, abstractmethod
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray


class Propagator(ABC):
    """Propagador de órbita: estado ECI em função do tempo desde a época."""

    @property
    @abstractmethod
    def epoch(self) -> datetime:
        """Época do propagador (UTC, tz-aware): o instante t = 0."""

    @abstractmethod
    def propagate(self, t: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """t: segundos desde ``epoch``, shape () ou (N,).

        Retorna (r, v) em ECI, SI [m, m/s], sempre com shape (N, 3).
        """

    def propagate_utc(
        self, epoch: datetime, t: ArrayLike
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estado ECI nos instantes ``epoch + t``, para qualquer ``epoch`` (S3-02).

        Converte para segundos desde ``self.epoch`` e chama ``propagate``. Útil quando vários
        propagadores, cada um com a sua época (ex.: o TLE de cada satélite), são amostrados numa
        grade de tempo comum. Propagadores que exigem t ≥ 0 (``CowellPropagator``) levantam
        ``ValueError`` se algum instante cair antes de ``self.epoch``.

        Args:
            epoch: referência dos tempos ``t``, ``datetime`` tz-aware.
            t: segundos desde ``epoch``, escalar ou (N,).

        Returns:
            (r, v) em ECI [m, m/s], cada um com shape (N, 3).

        Raises:
            ValueError: se ``epoch`` for naive.
        """
        if epoch.tzinfo is None or epoch.utcoffset() is None:
            raise ValueError("epoch deve ser um datetime tz-aware (ex.: tzinfo=timezone.utc)")
        offset_s = (epoch - self.epoch).total_seconds()
        return self.propagate(as_time_array(t) + offset_s)


def as_time_array(t: ArrayLike) -> NDArray[np.float64]:
    """Converte escalar ou sequência de tempos [s] em ndarray float 1D de shape (N,).

    Raises:
        ValueError: se ``t`` tiver mais de uma dimensão.
    """
    arr = np.atleast_1d(np.asarray(t, dtype=float))
    if arr.ndim != 1:
        raise ValueError(f"t deve ser escalar ou 1D, recebido shape {arr.shape}")
    return arr
