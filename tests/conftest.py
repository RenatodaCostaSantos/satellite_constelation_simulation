"""Fixtures compartilhadas entre módulos de teste (Semana 5).

A constelação de referência do Apêndice A.1 do prompt da Semana 5 (6 satélites, Cowell completo,
i calibrada dada) é construída **uma vez por sessão**: o ``CowellPropagator`` guarda em cache a
saída densa já integrada, então as janelas de 24 h usadas em vários módulos (grade, cobertura,
contatos) integram o Cowell uma única vez. Os testes não devem modificar os objetos.

Semana 6 (S6-01): a constelação de 18 satélites do Apêndice A.1 do prompt da Semana 6 (mesma
órbita, ``DESIGN_91_6_P18``) e os seus nós ascendentes de 1 dia (≈ 20 s, a integração de 1 dia
dos 18 no Cowell completo) também são fixtures de sessão: o S6-02 as reaproveita.
"""

import math
import time
from datetime import datetime

import pytest
from golden.reference_values import CONSTELLATION_A1

from satsim.constants import SECONDS_PER_DAY
from satsim.constellation import DESIGN_91_6_P18, Constellation, sso_constellation
from satsim.forces import (
    DragAcceleration,
    SolarRadiationPressure,
    ThirdBodyAcceleration,
    ZonalAcceleration,
)
from satsim.groundtrack import node_crossings

A1_EPOCH = datetime.fromisoformat(CONSTELLATION_A1.epoch_utc)


def _full_models(epoch: datetime) -> tuple:
    a1 = CONSTELLATION_A1
    return (
        ZonalAcceleration(),
        ThirdBodyAcceleration("sun", epoch),
        ThirdBodyAcceleration("moon", epoch),
        SolarRadiationPressure(epoch, cr=a1.srp_cr, area_over_mass=a1.srp_area_over_mass),
        DragAcceleration(cd=2.2, area_over_mass=a1.drag_cd_a_m / 2.2),
    )


@pytest.fixture(scope="session")
def a1_models() -> tuple:
    """Cowell completo do A.1 (A.7 da Semana 4): J2–J4 + Sol/Lua + SRP + arrasto geodésico."""
    return _full_models(A1_EPOCH)


@pytest.fixture(scope="session")
def reference_constellation_built(a1_models) -> tuple[Constellation, float]:
    """A constelação de referência do A.1 (a de projeto, i dada) e o tempo de construção [s]."""
    t0 = time.perf_counter()
    constellation = sso_constellation(
        A1_EPOCH, a1_models, i=math.radians(CONSTELLATION_A1.i_calibrated_deg)
    )
    return constellation, time.perf_counter() - t0


@pytest.fixture(scope="session")
def reference_constellation(reference_constellation_built) -> Constellation:
    return reference_constellation_built[0]


@pytest.fixture(scope="session")
def constellation18_built(a1_models) -> tuple[Constellation, float]:
    """A constelação de 18 do A.1 da Semana 6 (a de projeto, i dada) e o tempo de construção."""
    t0 = time.perf_counter()
    constellation = DESIGN_91_6_P18.build(
        A1_EPOCH, a1_models, i=math.radians(CONSTELLATION_A1.i_calibrated_deg)
    )
    return constellation, time.perf_counter() - t0


@pytest.fixture(scope="session")
def constellation18(constellation18_built) -> Constellation:
    return constellation18_built[0]


@pytest.fixture(scope="session")
def constellation18_nodes(constellation18) -> list[tuple]:
    """Nós ascendentes de 1 dia de cada um dos 18, ``node_crossings`` a 30 s: [(t, lon), …] em
    ordem de satélite. Integra o Cowell completo de 1 dia dos 18 (cache para a sessão)."""
    return [node_crossings(s, A1_EPOCH, SECONDS_PER_DAY, 30.0) for s in constellation18]
