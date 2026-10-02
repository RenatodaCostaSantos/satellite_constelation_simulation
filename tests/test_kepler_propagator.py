"""Testes do propagador kepleriano (S1-05)."""

import math
from datetime import UTC, datetime

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from satsim.constants import MU
from satsim.elements import KeplerianElements, elements_to_state
from satsim.propagators import KeplerPropagator, Propagator, as_time_array

DEG = math.pi / 180.0
EPOCH = datetime(2026, 9, 29, 12, 0, 0, tzinfo=UTC)

# Órbita de projeto 91/6 (Sem4, Apêndice A.2 do prompt): a = 6888,089 km, e = 1,074e-3,
# i = 97,4396°, ω = 90° (congelada); raan e psi como no card S1-04.
DESIGN = KeplerianElements(
    a=6888.089e3, e=1.074e-3, i=97.4396 * DEG, raan=30.0 * DEG, argp=90.0 * DEG, psi=45.0 * DEG
)


@pytest.fixture
def prop() -> KeplerPropagator:
    return KeplerPropagator(DESIGN, EPOCH)


def test_is_propagator(prop: KeplerPropagator) -> None:
    assert isinstance(prop, Propagator)
    assert prop.epoch == EPOCH


def test_period_design_orbit(prop: KeplerPropagator) -> None:
    # Sem4, Apêndice A.2 do prompt: Tkep = 5689,30 s para a = 6888,089 km
    assert prop.period == pytest.approx(5689.30, abs=0.05)
    assert prop.mean_motion == pytest.approx(2 * math.pi / prop.period, rel=1e-15)


def test_conservation_30_days(prop: KeplerPropagator) -> None:
    t = np.arange(0.0, 30 * 86400.0 + 1.0, 60.0)
    r, v = prop.propagate(t)
    energy = 0.5 * np.sum(v * v, axis=1) - MU / np.linalg.norm(r, axis=1)
    ell = np.cross(r, v)
    ell_norm = np.linalg.norm(ell, axis=1)
    ell_hat = ell / ell_norm[:, None]

    np.testing.assert_allclose(energy, energy[0], rtol=1e-10)
    np.testing.assert_allclose(ell_norm, ell_norm[0], rtol=1e-10)
    # direção de ℓ constante: ângulo entre ℓ̂(t) e ℓ̂(0), via |ℓ̂(t) × ℓ̂(0)| = sin(ângulo).
    # (arccos do produto escalar não serve: arccos(1 − ε) ≈ 1,5e-8 já no arredondamento.)
    sin_angle = np.linalg.norm(np.cross(ell_hat, ell_hat[0]), axis=1)
    assert np.max(sin_angle) < 1e-10


def test_returns_after_one_period(prop: KeplerPropagator) -> None:
    r, v = prop.propagate([0.0, prop.period])
    assert np.linalg.norm(r[1] - r[0]) < 1e-4  # m
    assert np.linalg.norm(v[1] - v[0]) < 1e-7  # m/s


def test_circular_velocity() -> None:
    # Sem4, Apêndice A.2 do prompt: v (circular) = 7607 m/s para a = 6888,089 km
    el = KeplerianElements(6888.089e3, 0.0, 97.4396 * DEG, 0.0, 0.0, 0.0)
    _, v = KeplerPropagator(el, EPOCH).propagate(np.linspace(0.0, 6000.0, 50))
    np.testing.assert_allclose(np.linalg.norm(v, axis=1), 7607.0, atol=1.0)


def test_propagate_zero_matches_elements_to_state(prop: KeplerPropagator) -> None:
    r, v = prop.propagate([0.0])
    r0, v0 = elements_to_state(DESIGN)
    np.testing.assert_array_equal(r, r0)
    np.testing.assert_array_equal(v, v0)


def test_max_latitude_one_orbit(prop: KeplerPropagator) -> None:
    # Sem4 eq. 12 / Apêndice A.2 do prompt: latitude máxima = 180° − i = 82,56°
    r, _ = prop.propagate(np.linspace(0.0, prop.period, 20001))
    lat = np.degrees(np.arcsin(r[:, 2] / np.linalg.norm(r, axis=1)))
    assert np.max(lat) == pytest.approx(82.56, abs=0.01)


def test_cross_validation_solve_ivp(prop: KeplerPropagator) -> None:
    def two_body(_t: float, s: np.ndarray) -> np.ndarray:
        r = s[:3]
        return np.concatenate([s[3:], -MU * r / np.linalg.norm(r) ** 3])

    r0, v0 = prop.propagate(0.0)
    t_eval = np.linspace(0.0, 86400.0, 145)  # a cada 10 min, por 1 dia
    sol = solve_ivp(
        two_body,
        (0.0, 86400.0),
        np.concatenate([r0[0], v0[0]]),
        method="DOP853",
        rtol=1e-13,
        atol=1e-6,
        t_eval=t_eval,
    )
    assert sol.success
    r_kep, _ = prop.propagate(t_eval)
    err = np.linalg.norm(sol.y[:3].T - r_kep, axis=1)
    assert np.max(err) < 1.0  # m


def test_naive_epoch_rejected() -> None:
    with pytest.raises(ValueError):
        KeplerPropagator(DESIGN, datetime(2026, 9, 29, 12, 0, 0))


@pytest.mark.parametrize("t", [0.0, 1234.5, np.float64(10.0), [5.0], (1.0, 2.0, 3.0)])
def test_shapes_scalar_and_array(prop: KeplerPropagator, t) -> None:
    n = len(np.atleast_1d(t))
    r, v = prop.propagate(t)
    assert r.shape == (n, 3) and v.shape == (n, 3)


def test_array_matches_pointwise(prop: KeplerPropagator) -> None:
    t = np.linspace(-3000.0, 50000.0, 23)
    r, v = prop.propagate(t)
    for k, tk in enumerate(t):
        rk, vk = prop.propagate(float(tk))
        np.testing.assert_array_equal(r[k], rk[0])
        np.testing.assert_array_equal(v[k], vk[0])


def test_as_time_array() -> None:
    assert as_time_array(3.0).shape == (1,)
    assert as_time_array([1, 2, 3]).dtype == np.float64
    with pytest.raises(ValueError):
        as_time_array([[1.0, 2.0]])
