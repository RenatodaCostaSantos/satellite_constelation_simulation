"""Propagadores de órbita: todos implementam ``Propagator.propagate(t) → (r, v)``."""

from satsim.propagators.base import Propagator, as_time_array
from satsim.propagators.kepler import KeplerPropagator

__all__ = ["KeplerPropagator", "Propagator", "as_time_array"]
