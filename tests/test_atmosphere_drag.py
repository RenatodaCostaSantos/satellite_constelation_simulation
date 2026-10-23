"""Testes da atmosfera, do arrasto, do decaimento e do Cowell com arrasto × SGP4 (S3-03).

Referências: tests/golden/reference_values.py — ATMOSPHERE_A3, DECAY_TABLE3, DRAG_DECAY_A4,
LIFETIMES_A4, COWELL_VS_SGP4_A5 e BALLISTIC_FIT_A5 (Apêndices A.3–A.5 do prompt da Semana 3;
Sem3 eqs. 5–6, 14, 15, 22 e Tabela 3).
"""

import logging
import math
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pytest
from golden.reference_values import (
    ATMOSPHERE_A3,
    BALLISTIC_FIT_A5,
    COWELL_VS_SGP4_A5,
    DECAY_TABLE3,
    DRAG_DECAY_A4,
    LIFETIMES_A4,
)

from satsim.analysis import cowell_j2_drag, fit_ballistic_coefficient, position_error
from satsim.constants import MU, OMEGA_EARTH, R_EARTH, WGS84_A, WGS84_B
from satsim.decay import decay_per_rev, decay_rate, lifetime
from satsim.elements import state_to_elements
from satsim.forces import DragAcceleration, ExponentialAtmosphere, J2Acceleration, altitude
from satsim.propagators import CowellPropagator, SGP4Propagator

EPOCH = datetime(2026, 10, 16, 0, 0, 0, tzinfo=UTC)
SAMPLE_TEXT = (Path(__file__).parent / "data" / "iss_sample.tle").read_text()
D = DRAG_DECAY_A4


def circular_state(h_m: float, polar: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Estado circular no equador (x = R⊕ + h), velocidade em y (equatorial) ou z (polar)."""
    r = R_EARTH + h_m
    v = math.sqrt(MU / r)
    return np.array([r, 0.0, 0.0]), np.array([0.0, 0.0, v] if polar else [0.0, v, 0.0])


# --- Atmosfera (A.3, Sem3 Tabela 3) ------------------------------------------------------------


def test_density_at_table_points_a3():
    atm = ExponentialAtmosphere()
    for h0_km, rho0, _ in ATMOSPHERE_A3:
        assert atm.density(h0_km * 1e3) == pytest.approx(rho0, rel=1e-12)


@pytest.mark.parametrize("row", DECAY_TABLE3, ids=lambda r: f"{r.h_km:.0f}km")
def test_density_matches_sem3_table3_within_5pct(row):
    assert ExponentialAtmosphere().density(row.h_km * 1e3) == pytest.approx(
        row.rho_table3_kg_m3, rel=0.05
    )


def test_density_continuous_at_band_boundaries():
    atm = ExponentialAtmosphere()
    for h0_km, _, _ in ATMOSPHERE_A3[1:]:
        h0 = h0_km * 1e3
        below, above = atm.density(h0 - 1e-6), atm.density(h0)
        assert abs(above / below - 1.0) < 1e-3


def test_density_vectorized_and_monotonic():
    atm = ExponentialAtmosphere()
    h = np.linspace(100e3, 1200e3, 500)
    rho = atm.density(h)
    assert rho.shape == h.shape
    np.testing.assert_array_equal(rho, [atm.density(x) for x in h])
    assert np.all(np.diff(rho) < 0.0)
    assert isinstance(atm.density(400e3), float)


def test_density_below_100km_uses_first_band_with_warning(caplog):
    atm = ExponentialAtmosphere()
    with caplog.at_level(logging.WARNING, logger="satsim.forces.atmosphere"):
        rho = atm.density(np.array([80e3, 90e3]))
        atm.density(50e3)  # um aviso por instância
    np.testing.assert_array_equal(rho, ATMOSPHERE_A3[0][1])
    assert len(caplog.records) == 1 and "abaixo de 100 km" in caplog.records[0].message


def test_density_above_1000km_extrapolates_last_band():
    h0_km, rho0, scale_km = ATMOSPHERE_A3[-1]
    rho = ExponentialAtmosphere().density(1200e3)
    assert rho == pytest.approx(rho0 * math.exp(-200.0 / scale_km), rel=1e-12)


def test_atmosphere_rejects_bad_table():
    with pytest.raises(ValueError):
        ExponentialAtmosphere([(200.0, 1e-10, 37.0), (100.0, 1e-7, 6.0)])


# --- Aceleração de arrasto (Sem3 §2.4, eqs. 5–6) ----------------------------------------------


def test_drag_acceleration_magnitude_400km_sem3():
    drag = DragAcceleration(D.cd, D.area_over_mass_m2_kg, altitude_model="spherical",
                            corotating=False)  # fmt: skip
    r = np.array([R_EARTH + 400e3, 0.0, 0.0])
    v = np.array([0.0, D.v_400_m_s, 0.0])
    a = np.linalg.norm(drag(0.0, r, v))
    assert a == pytest.approx(D.drag_accel_400_m_s2, rel=0.03)
    g = MU / np.linalg.norm(r) ** 2
    assert g == pytest.approx(D.g_400_m_s2, rel=1e-3)
    assert a / g == pytest.approx(D.drag_over_gravity_400, rel=0.03)


def _cos(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


def test_drag_direction_antiparallel_to_relative_velocity():
    r = np.array([5.0e6, 3.0e6, 2.5e6])
    v = np.array([-3.0e3, 4.0e3, 5.0e3])
    drag = DragAcceleration()
    assert _cos(drag(0.0, r, v), drag.relative_velocity(r, v)) == pytest.approx(-1.0, abs=1e-12)
    no_rot = DragAcceleration(corotating=False)
    assert _cos(no_rot(0.0, r, v), v) == pytest.approx(-1.0, abs=1e-12)


def test_relative_velocity_equatorial_point_polar_orbit():
    r, v = circular_state(400e3, polar=True)
    v_rel = DragAcceleration().relative_velocity(r, v)
    expected = math.hypot(np.linalg.norm(v), OMEGA_EARTH * np.linalg.norm(r))
    assert np.linalg.norm(v_rel) == pytest.approx(expected, rel=1e-14)


def test_drag_vectorized_matches_scalar():
    rng = np.random.default_rng(1)
    r = rng.normal(size=(5, 3))
    r = (R_EARTH + 450e3) * r / np.linalg.norm(r, axis=1, keepdims=True)
    v = rng.normal(size=(5, 3)) * 4e3
    drag = DragAcceleration()
    a = drag(0.0, r, v)
    assert a.shape == (5, 3)
    for k in range(5):
        np.testing.assert_allclose(a[k], drag(0.0, r[k], v[k]), rtol=1e-14)


def test_geodetic_vs_spherical_altitude():
    # no equador as duas coincidem; sobre o polo a geodésica é 21,4 km maior (raio polar menor)
    eq = np.array([WGS84_A + 400e3, 0.0, 0.0])
    pole = np.array([0.0, 0.0, WGS84_B + 400e3])
    assert altitude(eq, "geodetic") == pytest.approx(400e3, abs=1e-3)
    assert altitude(eq, "spherical") == pytest.approx(400e3, abs=1e-3)
    assert altitude(pole, "geodetic") == pytest.approx(400e3, abs=1e-3)
    assert altitude(pole, "geodetic") - altitude(pole, "spherical") == pytest.approx(
        WGS84_A - WGS84_B, abs=1e-3
    )


def test_drag_rejects_bad_parameters():
    with pytest.raises(ValueError):
        DragAcceleration(cd=-1.0)
    with pytest.raises(ValueError):
        DragAcceleration(altitude_model="plana")


# --- Decaimento (Sem3 eqs. 14, 15, 22; Tabela 3; A.4) ------------------------------------------


def test_decay_per_rev_eq15():
    assert decay_per_rev(400e3, D.cd_a_m_m2_kg) == pytest.approx(D.da_per_rev_400_m, rel=0.03)


@pytest.mark.parametrize("row", DECAY_TABLE3, ids=lambda r: f"{r.h_km:.0f}km")
def test_decay_rate_sem3_table3(row):
    dhdt_m_day = -decay_rate(row.h_km * 1e3, D.cd_a_m_m2_kg) * 86400.0
    assert dhdt_m_day == pytest.approx(row.dhdt_table3_m_day, rel=0.05)
    # arredondado na precisão impressa no A.4, o valor calculado é exatamente o do A.4
    assert round(dhdt_m_day, row.dhdt_a4_decimals) == row.dhdt_a4_m_day


def test_decay_rate_vectorized():
    h = np.array([row.h_km for row in DECAY_TABLE3]) * 1e3
    rates = decay_rate(h, D.cd_a_m_m2_kg)
    assert rates.shape == h.shape and np.all(rates < 0.0)
    np.testing.assert_allclose(rates, [decay_rate(x, D.cd_a_m_m2_kg) for x in h], rtol=1e-15)


@pytest.mark.parametrize("case", LIFETIMES_A4, ids=lambda c: f"{c.h0_km:.0f}km-{c.cd_a_m_m2_kg}")
def test_lifetime_a4(case):
    years = lifetime(case.h0_km * 1e3, case.cd_a_m_m2_kg, h_reentry_m=D.h_reentry_km * 1e3)
    assert years == pytest.approx(case.lifetime_years, rel=0.05)


def test_lifetime_edge_cases():
    assert lifetime(100e3, 0.011) == 0.0
    assert lifetime(500e3, 0.0) == math.inf


def test_cowell_drag_only_decay_matches_eq15_and_day():
    """Cowell só com arrasto (sem J2, altitude esférica, sem co-rotação) a 400 km."""
    r0, v0 = circular_state(400e3)
    drag = DragAcceleration(D.cd, D.area_over_mass_m2_kg, altitude_model="spherical",
                            corotating=False)  # fmt: skip
    cowell = CowellPropagator.from_state(r0, v0, EPOCH, (drag,))
    a0 = R_EARTH + 400e3
    period = 2.0 * math.pi * math.sqrt(a0**3 / MU)
    r, v = cowell.propagate([0.0, period, 86400.0])
    a = np.asarray(state_to_elements(r, v).a)
    assert a[1] - a[0] == pytest.approx(D.da_per_rev_400_m, rel=0.05)
    assert a[0] - a[2] == pytest.approx(D.da_per_day_400_m, rel=0.05)


def test_negligible_drag_reproduces_j2_only():
    r0, v0 = circular_state(420e3)
    t = np.linspace(0.0, 86400.0, 25)
    j2_only = CowellPropagator.from_state(r0, v0, EPOCH, (J2Acceleration(),))
    no_drag = CowellPropagator.from_state(
        r0, v0, EPOCH, (J2Acceleration(), DragAcceleration(area_over_mass=0.0))
    )
    np.testing.assert_allclose(no_drag.propagate(t)[0], j2_only.propagate(t)[0], rtol=0, atol=1e-3)


def test_reentry_event_and_nan_after(caplog):
    r0, v0 = circular_state(150e3)
    cowell = CowellPropagator.from_state(
        r0, v0, EPOCH, (J2Acceleration(), DragAcceleration(area_over_mass=0.02))
    )
    with caplog.at_level(logging.WARNING, logger="satsim.propagators.cowell"):
        r, v = cowell.propagate([0.0, 600.0, 86400.0])
    assert cowell.reentry_time is not None and cowell.reentry_time < 86400.0
    assert np.all(np.isfinite(r[:2])) and np.all(np.isnan(r[2])) and np.all(np.isnan(v[2]))
    assert any("reentrada" in rec.message for rec in caplog.records)
    assert any("NaN" in rec.message for rec in caplog.records)
    # o estado no evento está na altitude de parada (100 km)
    r_evt, _ = cowell.propagate(cowell.reentry_time)
    assert float(altitude(r_evt[0])) == pytest.approx(100e3, abs=1.0)


def test_from_state_equals_constructor():
    r0, v0 = circular_state(500e3)
    a = CowellPropagator.from_state(r0, v0, EPOCH)
    b = CowellPropagator((r0, v0), EPOCH)
    np.testing.assert_array_equal(a.propagate(3600.0)[0], b.propagate(3600.0)[0])


# --- Cowell (J2 + arrasto) × SGP4, TLE de exemplo (A.5) ----------------------------------------


@pytest.fixture(scope="module")
def iss() -> SGP4Propagator:
    return SGP4Propagator.from_tle_text(SAMPLE_TEXT)


@pytest.fixture(scope="module")
def fitted_cd_a_m(iss) -> float:
    return fit_ballistic_coefficient(iss)  # janela de 72 h, a cada 600 s (≈ 12 s)


T_CHECK = np.array([24.0, 72.0, 168.0]) * 3600.0


@pytest.mark.parametrize("row", COWELL_VS_SGP4_A5[:3], ids=lambda r: r.label)
def test_cowell_vs_sgp4_reproduces_a5(iss, row):
    cowell = cowell_j2_drag(iss, row.cd_a_m_m2_kg, altitude_model=row.altitude_model)
    err_km = position_error(cowell, iss, T_CHECK) / 1e3
    expected = [row.err_24h_km, row.err_72h_km, row.err_7d_km]
    np.testing.assert_allclose(err_km, expected, rtol=0.06)


def test_fitted_ballistic_coefficient_a5(iss, fitted_cd_a_m):
    lo, hi = BALLISTIC_FIT_A5.fitted_range_m2_kg
    assert lo <= fitted_cd_a_m <= hi
    assert fitted_cd_a_m == pytest.approx(BALLISTIC_FIT_A5.fitted_cd_a_m_m2_kg, rel=0.02)
    err_km = position_error(cowell_j2_drag(iss, fitted_cd_a_m), iss, T_CHECK) / 1e3
    assert np.all(err_km <= np.array(BALLISTIC_FIT_A5.max_err_km))  # critério: ≤ 2 / 8 / 40 km
    row = COWELL_VS_SGP4_A5[3]
    np.testing.assert_allclose(err_km, [row.err_24h_km, row.err_72h_km, row.err_7d_km],
                               rtol=0.1)  # fmt: skip
