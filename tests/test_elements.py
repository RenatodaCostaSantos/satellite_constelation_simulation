"""Testes de elementos orbitais ↔ estado cartesiano (S1-04)."""

import math

import numpy as np
import pytest

from satsim.constants import MU
from satsim.elements import (
    KeplerianElements,
    eccentric_to_true,
    elements_to_state,
    mean_to_true,
    solve_kepler,
    state_to_elements,
    true_to_eccentric,
    true_to_mean,
)

DEG = math.pi / 180.0

# Órbita de projeto 91/6 (Sem4, Apêndice A.2 do prompt): a, e e i da SSO congelada;
# raan, argp e psi são os valores pedidos pelo card S1-04.
DESIGN = KeplerianElements(
    a=6888.089e3, e=1.074e-3, i=97.4396 * DEG, raan=30.0 * DEG, argp=90.0 * DEG, psi=45.0 * DEG
)


def _angle_diff(x: float, y: float) -> float:
    """Diferença angular com wrap em 2π, em [0, π]."""
    d = (x - y) % (2.0 * math.pi)
    return min(d, 2.0 * math.pi - d)


def _assert_states_close(r1, v1, r2, v2) -> None:
    np.testing.assert_allclose(r1, r2, rtol=0.0, atol=1e-6)  # m
    np.testing.assert_allclose(v1, v2, rtol=0.0, atol=1e-9)  # m/s


# --- Equação de Kepler -------------------------------------------------------------------------


@pytest.mark.parametrize("e", [0.0, 1e-3, 0.1, 0.7, 0.95])
def test_kepler_residual(e: float) -> None:
    psi = np.linspace(-math.pi, math.pi, 1000)
    big_e = solve_kepler(psi, e)
    assert np.max(np.abs(big_e - e * np.sin(big_e) - psi)) < 1e-12


def test_kepler_preserves_revolution() -> None:
    psi = np.array([-7.0, 3.5, 12.0])
    big_e = solve_kepler(psi, 0.3)
    np.testing.assert_allclose(big_e - 0.3 * np.sin(big_e), psi, rtol=0.0, atol=1e-12)


def test_kepler_scalar_returns_float() -> None:
    assert isinstance(solve_kepler(1.0, 0.1), float)


@pytest.mark.parametrize("e", [-0.1, 1.0, 1.5])
def test_kepler_rejects_non_elliptic(e: float) -> None:
    with pytest.raises(ValueError):
        solve_kepler(0.5, e)


def test_anomaly_conversions_roundtrip() -> None:
    nu = np.linspace(-math.pi + 1e-3, math.pi, 500)
    for e in (0.0, 0.01, 0.5, 0.9):
        np.testing.assert_allclose(eccentric_to_true(true_to_eccentric(nu, e), e), nu, atol=1e-12)
        np.testing.assert_allclose(mean_to_true(true_to_mean(nu, e), e), nu, atol=1e-11)


# --- Órbita de projeto: ida e volta ------------------------------------------------------------


def test_design_elements_state_elements() -> None:
    r, v = elements_to_state(DESIGN)
    el = state_to_elements(r[0], v[0])
    assert abs(el.a - DESIGN.a) < 1e-6
    assert abs(el.e - DESIGN.e) < 1e-9
    for name in ("i", "raan", "argp", "psi"):
        assert _angle_diff(getattr(el, name), getattr(DESIGN, name)) < 1e-9, name


def test_design_state_elements_state() -> None:
    r0, v0 = elements_to_state(DESIGN)
    r1, v1 = elements_to_state(state_to_elements(r0[0], v0[0]))
    _assert_states_close(r0, v0, r1, v1)


# --- Casos degenerados (Sem4 §3.3) -------------------------------------------------------------


@pytest.mark.parametrize(
    "el",
    [
        # circular inclinada
        KeplerianElements(6888.089e3, 0.0, 97.4396 * DEG, 30 * DEG, 0.0, 45 * DEG),
        # circular inclinada com argp não nulo na entrada (vai para u)
        KeplerianElements(6888.089e3, 0.0, 51.64 * DEG, 200 * DEG, 70 * DEG, 10 * DEG),
        # equatorial excêntrica
        KeplerianElements(7000e3, 0.1, 0.0, 0.0, 60 * DEG, 120 * DEG),
        # circular equatorial
        KeplerianElements(7000e3, 0.0, 0.0, 0.0, 0.0, 250 * DEG),
        # retrógrada equatorial excêntrica
        KeplerianElements(7000e3, 0.1, math.pi, 0.0, 60 * DEG, 120 * DEG),
        # retrógrada equatorial circular
        KeplerianElements(7000e3, 0.0, math.pi, 0.0, 0.0, 300 * DEG),
    ],
    ids=["circ-incl", "circ-incl-argp", "equat-exc", "circ-equat", "retro-exc", "retro-circ"],
)
def test_degenerate_roundtrip(el: KeplerianElements) -> None:
    r0, v0 = elements_to_state(el)
    back = state_to_elements(r0[0], v0[0])
    assert np.isfinite([back.a, back.e, back.i, back.raan, back.argp, back.psi]).all()
    r1, v1 = elements_to_state(back)
    _assert_states_close(r0, v0, r1, v1)


def test_degenerate_conventions() -> None:
    # circular: argp = 0 e todo o ângulo em u = argp + nu
    el = KeplerianElements(6888.089e3, 0.0, 51.64 * DEG, 200 * DEG, 70 * DEG, 10 * DEG)
    back = state_to_elements(*(x[0] for x in elements_to_state(el)))
    assert back.argp == 0.0
    assert _angle_diff(back.psi, 80 * DEG) < 1e-9
    assert _angle_diff(back.raan, 200 * DEG) < 1e-9
    # equatorial: raan = 0, linha dos nós no eixo x
    el = KeplerianElements(7000e3, 0.1, 0.0, 0.0, 60 * DEG, 120 * DEG)
    back = state_to_elements(*(x[0] for x in elements_to_state(el)))
    assert back.raan == 0.0
    assert _angle_diff(back.argp, 60 * DEG) < 1e-9
    # circular equatorial: raan = argp = 0, longitude verdadeira em psi
    el = KeplerianElements(7000e3, 0.0, 0.0, 0.0, 0.0, 250 * DEG)
    back = state_to_elements(*(x[0] for x in elements_to_state(el)))
    assert back.raan == 0.0 and back.argp == 0.0
    assert _angle_diff(back.psi, 250 * DEG) < 1e-9


def test_state_to_elements_rejects_hyperbolic() -> None:
    r = np.array([7000e3, 0.0, 0.0])
    v = np.array([0.0, 1.5 * math.sqrt(2 * MU / 7000e3), 0.0])
    with pytest.raises(ValueError):
        state_to_elements(r, v)


# --- Invariantes -------------------------------------------------------------------------------


def _design_orbit_samples(n: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    psi = np.linspace(0.0, 2 * math.pi, n)
    el = KeplerianElements(DESIGN.a, DESIGN.e, DESIGN.i, DESIGN.raan, DESIGN.argp, psi)
    return elements_to_state(el)


def test_invariants() -> None:
    a, e = DESIGN.a, DESIGN.e
    r, v = _design_orbit_samples()
    r_norm = np.linalg.norm(r, axis=1)
    assert np.all(r_norm >= a * (1 - e) * (1 - 1e-14))
    assert np.all(r_norm <= a * (1 + e) * (1 + 1e-14))

    energy = 0.5 * np.sum(v * v, axis=1) - MU / r_norm
    np.testing.assert_allclose(energy, -MU / (2 * a), rtol=1e-12)

    ell = np.cross(r, v)
    ell_norm = np.linalg.norm(ell, axis=1)
    np.testing.assert_allclose(ell_norm, math.sqrt(MU * a * (1 - e**2)), rtol=1e-12)
    np.testing.assert_allclose(np.sum(ell * r, axis=1) / (ell_norm * r_norm), 0.0, atol=1e-14)

    # i > 90° ⇒ órbita retrógrada ⇒ ℓz < 0
    assert DESIGN.i > math.pi / 2
    assert np.all(ell[:, 2] < 0)


def test_sem4_eq12_latitude() -> None:
    # Sem4 eq. 12: z/r = sin i · sin u, com u = argp + nu
    psi = np.linspace(0.0, 2 * math.pi, 20001)
    el = KeplerianElements(DESIGN.a, DESIGN.e, DESIGN.i, DESIGN.raan, DESIGN.argp, psi)
    r, _ = elements_to_state(el)
    r_norm = np.linalg.norm(r, axis=1)
    u = DESIGN.argp + np.asarray(mean_to_true(psi, DESIGN.e))
    np.testing.assert_allclose(r[:, 2] / r_norm, math.sin(DESIGN.i) * np.sin(u), atol=1e-12)

    # Latitude máxima = 180° − i = 82,56° (Sem4, Apêndice A.2 do prompt)
    lat_max_deg = math.degrees(np.max(np.arcsin(r[:, 2] / r_norm)))
    assert lat_max_deg == pytest.approx(82.56, abs=0.01)


# --- Vetorização -------------------------------------------------------------------------------


def test_vectorization_matches_scalar() -> None:
    psi = np.linspace(-math.pi, 3 * math.pi, 37)
    el = KeplerianElements(DESIGN.a, DESIGN.e, DESIGN.i, DESIGN.raan, DESIGN.argp, psi)
    r, v = elements_to_state(el)
    assert r.shape == (37, 3) and v.shape == (37, 3)
    for k, p in enumerate(psi):
        rk, vk = elements_to_state(
            KeplerianElements(DESIGN.a, DESIGN.e, DESIGN.i, DESIGN.raan, DESIGN.argp, float(p))
        )
        assert rk.shape == (1, 3)
        np.testing.assert_array_equal(r[k], rk[0])
        np.testing.assert_array_equal(v[k], vk[0])


def test_state_to_elements_vectorized() -> None:
    r, v = _design_orbit_samples(50)
    el = state_to_elements(r, v)
    assert el.a.shape == (50,)
    for k in (0, 17, 49):
        single = state_to_elements(r[k], v[k])
        assert el.a[k] == pytest.approx(single.a, rel=1e-15)
        assert el.psi[k] == pytest.approx(single.psi, abs=1e-15)
