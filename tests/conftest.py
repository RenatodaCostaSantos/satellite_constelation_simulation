"""Fixtures compartilhadas entre módulos de teste (Semana 5).

A constelação de referência do Apêndice A.1 do prompt da Semana 5 (6 satélites, Cowell completo,
i calibrada dada) é construída **uma vez por sessão**: o ``CowellPropagator`` guarda em cache a
saída densa já integrada, então as janelas de 24 h usadas em vários módulos (grade, cobertura,
contatos) integram o Cowell uma única vez. Os testes não devem modificar os objetos.
"""

import math
import time
from datetime import datetime

import pytest
from golden.reference_values import CONSTELLATION_A1

from satsim.constellation import Constellation, sso_constellation
from satsim.forces import (
    DragAcceleration,
    SolarRadiationPressure,
    ThirdBodyAcceleration,
    ZonalAcceleration,
)

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
