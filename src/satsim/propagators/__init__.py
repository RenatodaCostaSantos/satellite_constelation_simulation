"""Propagadores de órbita: todos implementam ``Propagator.propagate(t) → (r, v)``."""

from satsim.propagators.base import Propagator, as_time_array
from satsim.propagators.cowell import CowellPropagator
from satsim.propagators.kepler import KeplerPropagator
from satsim.propagators.mean_j2 import MeanJ2Propagator

__all__ = [
    "CowellPropagator",
    "KeplerPropagator",
    "MeanJ2Propagator",
    "Propagator",
    "as_time_array",
]
