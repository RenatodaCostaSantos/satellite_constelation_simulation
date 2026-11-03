"""Satélite = dados + propagador (ADR 0001, decisão 2 e seção Interfaces; card S5-01).

Um ``Satellite`` junta a identidade e os dados de exibição de um satélite (id, nome, cor) e,
opcionalmente, os seus parâmetros físicos (C_D, A/m) a um ``Propagator`` qualquer (Kepler,
MeanJ2, MeanZonal, Cowell, SGP4). Toda a física fica no propagador: ``propagate`` e
``propagate_utc`` apenas delegam a ele, com o mesmo contrato de shapes (t: (N,) → r, v: (N, 3)).

Adicionar satélites é instanciar mais ``Satellite``; ``satsim.constellation.Constellation``
empilha vários em arrays (P, N, 3). Substitui o ``SceneSatellite`` provisório da cena
(``satsim.app.scene_export``), que continua existindo por compatibilidade.
"""

from dataclasses import dataclass
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.propagators import Propagator

DEFAULT_COLOR = "#00e5ff"
"""Cor CSS padrão de um satélite na cena (a da SSO-1 das semanas anteriores)."""


@dataclass(frozen=True)
class Satellite:
    """Satélite: identidade, dados de exibição e o propagador que o move (ADR, Interfaces).

    Attributes:
        id: identificador único e estável (ex.: ``"s1"``), usado como chave na cena JSON.
        name: nome de exibição (ex.: ``"S1"``).
        propagator: qualquer ``Propagator``; a época do satélite é a do propagador.
        color: cor CSS do satélite na cena (ex.: ``"#00e5ff"``).
        cd: coeficiente de arrasto C_D [-], só informativo (o arrasto, se houver, já está nas
            forças do propagador).
        area_over_mass: razão área/massa A/m [m²/kg], idem.
    """

    id: str
    name: str
    propagator: Propagator
    color: str = DEFAULT_COLOR
    cd: float | None = None
    area_over_mass: float | None = None

    @property
    def epoch(self) -> datetime:
        """Época do satélite: a do propagador (t = 0 de ``propagate``)."""
        return self.propagator.epoch

    def propagate(self, t: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estado ECI [m, m/s] em ``t`` segundos desde ``self.epoch``; cada um com shape (N, 3)."""
        return self.propagator.propagate(t)

    def propagate_utc(
        self, epoch: datetime, t: ArrayLike
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estado ECI nos instantes ``epoch + t`` (``Propagator.propagate_utc``, S3-02)."""
        return self.propagator.propagate_utc(epoch, t)
