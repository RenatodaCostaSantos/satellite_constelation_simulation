"""Testes do propagador de elementos médios com J2 (S2-01).

Referências: tests/golden/reference_values.py (DESIGN_ORBIT, Sem4 §7.5), Sem4 eqs. 23–25 e 34.
"""

import math
import time
from datetime import UTC, datetime

import numpy as np
import pytest
from golden.reference_values import DESIGN_ORBIT

from satsim.elements import KeplerianElements
from satsim.propagators import KeplerPropagator, MeanJ2Propagator, Propagator
from satsim.secular import nodal_period, secular_rates

DEG = math.pi / 180.0
EPOCH = datetime(2026, 9, 29, 12, 0, 0, tzinfo=UTC)

# Órbita de projeto 91/6 (DESIGN_ORBIT, Sem4 §7.5) com raan0 = 30°, argp0 = 90° e ψ0 = 45°.
DESIGN = KeplerianElements(
    a=DESIGN_ORBIT.a_km * 1e3,
    e=DESIGN_ORBIT.e,
    i=DESIGN_ORBIT.i_deg * DEG,
    raan=30.0 * DEG,
    argp=DESIGN_ORBIT.argp_deg * DEG,
    psi=45.0 * DEG,
)


@pytest.fixture
def prop() -> MeanJ2Propagator:
    return MeanJ2Propagator(DESIGN, EPOCH)


def test_is_propagator(prop: MeanJ2Propagator) -> None:
    assert isinstance(prop, Propagator)
    assert prop.epoch == EPOCH
    assert prop.elements == DESIGN


def test_rejects_naive_epoch() -> None:
    with pytest.raises(ValueError):
        MeanJ2Propagator(DESIGN, datetime(2026, 9, 29, 12, 0, 0))


def test_rates_and_period(prop: MeanJ2Propagator) -> None:
    expected = secular_rates(DESIGN.a, DESIGN.e, DESIGN.i)
    assert prop.rates == expected
    assert prop.period == nodal_period(DESIGN.a, DESIGN.e, DESIGN.i)
    # DESIGN_ORBIT (Sem4 §7.5): Tnod = 5696,70 s
    assert prop.period == pytest.approx(DESIGN_ORBIT.tnod_s, abs=0.05)


def test_shapes(prop: MeanJ2Propagator) -> None:
    r, v = prop.propagate(100.0)
    assert r.shape == (1, 3) and v.shape == (1, 3)
    r, v = prop.propagate(np.arange(5.0))
    assert r.shape == (5, 3) and v.shape == (5, 3)
    el = prop.elements_at(np.arange(5.0))
    for x in (el.a, el.e, el.i, el.raan, el.argp, el.psi):
        assert np.shape(x) == (5,)


def test_scalar_and_array_t_agree(prop: MeanJ2Propagator) -> None:
    t = np.array([0.0, 1234.5, 86400.0])
    r, v = prop.propagate(t)
    for k, tk in enumerate(t):
        rk, vk = prop.propagate(float(tk))
        np.testing.assert_array_equal(r[k], rk[0])
        np.testing.assert_array_equal(v[k], vk[0])


def test_t0_matches_kepler(prop: MeanJ2Propagator) -> None:
    r, v = prop.propagate(0.0)
    rk, vk = KeplerPropagator(DESIGN, EPOCH).propagate(0.0)
    np.testing.assert_array_equal(r, rk)
    np.testing.assert_array_equal(v, vk)


def test_zero_j2_matches_kepler() -> None:
    t = np.linspace(0.0, 86400.0, 101)
    r, v = MeanJ2Propagator(DESIGN, EPOCH, j2=0.0).propagate(t)
    rk, vk = KeplerPropagator(DESIGN, EPOCH).propagate(t)
    np.testing.assert_allclose(r, rk, rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(v, vk, rtol=0.0, atol=1e-9)


def test_radius_within_mean_ellipse(prop: MeanJ2Propagator) -> None:
    r, _ = prop.propagate(np.arange(0.0, 86400.0 + 1.0, 10.0))
    r_norm = np.linalg.norm(r, axis=1)
    a, e = DESIGN.a, DESIGN.e
    assert np.all(r_norm >= a * (1.0 - e) - 1e-6)
    assert np.all(r_norm <= a * (1.0 + e) + 1e-6)


def test_elements_after_one_day(prop: MeanJ2Propagator) -> None:
    el = prop.elements_at(86400.0)
    rates = prop.rates
    assert el.raan[0] == pytest.approx(DESIGN.raan + rates.raan_dot * 86400.0, rel=1e-15)
    assert el.argp[0] == pytest.approx(DESIGN.argp + rates.argp_dot * 86400.0, rel=1e-15)
    assert el.psi[0] == pytest.approx(DESIGN.psi + rates.psi_dot * 86400.0, rel=1e-15)
    assert el.a[0] == DESIGN.a and el.e[0] == DESIGN.e and el.i[0] == DESIGN.i
    # SSO: Ω avança ≈ 0,9856°/dia (Sem4 eq. 26)
    assert math.degrees(el.raan[0] - DESIGN.raan) == pytest.approx(0.9856, abs=1e-4)


def test_state_consistent_with_elements_at(prop: MeanJ2Propagator) -> None:
    # O plano orbital do estado propagado gira com Ω(t): ℓ̂ = (sin i sin Ω, −sin i cos Ω, cos i)
    t = np.array([0.0, 43200.0, 86400.0])
    r, v = prop.propagate(t)
    el = prop.elements_at(t)
    ell = np.cross(r, v)
    ell_hat = ell / np.linalg.norm(ell, axis=1, keepdims=True)
    expected = np.stack(
        [
            np.sin(el.i) * np.sin(el.raan),
            -np.sin(el.i) * np.cos(el.raan),
            np.cos(el.i),
        ],
        axis=1,
    )
    np.testing.assert_allclose(ell_hat, expected, rtol=0.0, atol=1e-12)


def test_one_day_at_10s_is_fast(prop: MeanJ2Propagator) -> None:
    t = np.arange(0.0, 86400.0 + 1.0, 10.0)
    assert t.shape == (8641,)
    start = time.perf_counter()
    r, v = prop.propagate(t)
    elapsed = time.perf_counter() - start
    assert r.shape == (8641, 3)
    assert elapsed < 1.0
