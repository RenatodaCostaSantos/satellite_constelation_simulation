"""Testes do modelo médio J2–J4, da órbita congelada e da SSO calibrada (S4-05).

Referências: tests/golden/reference_values.py (MEAN_ZONAL_A5 e SSO_CALIBRATED_A7, Apêndices A.5 e
A.7 do prompt da Semana 4; a diferença esperada no período nodal medido está documentada no
golden).
"""

import math
from datetime import datetime

import numpy as np
import pytest
from golden.reference_values import MEAN_ZONAL_A5, SSO_CALIBRATED_A7

from satsim.analysis import fit_rate, osculating_elements
from satsim.astrotime import datetime_to_jd
from satsim.constants import J3, J4, MU, R_EARTH
from satsim.elements import KeplerianElements, mean_to_true
from satsim.forces import (
    DragAcceleration,
    SolarRadiationPressure,
    ThirdBodyAcceleration,
    ZonalAcceleration,
)
from satsim.mean_osc import conversion_models, mean_to_osculating, window_mean_elements
from satsim.orbit_design import (
    frozen_eccentricity,
    ltan_drift_min_per_year,
    raan_for_ltan,
    sso_inclination,
    sso_inclination_calibrated,
)
from satsim.propagators import CowellPropagator, MeanJ2Propagator, MeanZonalPropagator
from satsim.secular import secular_rates, zonal_mean_rates

A5 = MEAN_ZONAL_A5
A7 = SSO_CALIBRATED_A7
EPOCH = datetime.fromisoformat(A7.epoch_utc)
MODELS = (ZonalAcceleration(),)
DEG_DAY = math.degrees(1.0) * 86400.0  # rad/s → °/dia
A_DESIGN = A5.a_km * 1e3
I_DESIGN = math.radians(A5.i_deg)


@pytest.fixture(scope="module")
def i_calibrated() -> float:
    return sso_inclination_calibrated(A_DESIGN, A7.e_frozen, EPOCH, MODELS)


def _calibrated_mean(i: float) -> KeplerianElements:
    raan = math.radians(A7.raan0_mean_vs_cowell_deg)
    return KeplerianElements(a=A_DESIGN, e=A7.e_frozen, i=i, raan=raan, argp=math.pi / 2, psi=0.0)


@pytest.fixture(scope="module")
def cowell_30d(i_calibrated):
    """SSO calibrada, estado convertido, Cowell J2–J4, 30 dias a 60 s: Ω̇ e u̇ ajustados."""
    mean = _calibrated_mean(i_calibrated)
    osc = mean_to_osculating(mean, EPOCH, MODELS).elements
    prop = CowellPropagator(osc, EPOCH, MODELS)
    t = np.arange(0.0, 30 * 86400.0 + 1.0, 60.0)
    el = osculating_elements(prop, t)
    u = el.argp + mean_to_true(el.psi, el.e)
    return mean, t, fit_rate(t, el.raan), fit_rate(t, u)


# --- Taxas de 1ª ordem (A.5) -----------------------------------------------------------------


def test_zonal_rates_reduce_to_j2() -> None:
    a = np.array([6700e3, 6888.089e3, 7200e3, 12000e3])
    e = np.array([0.0, 1.074e-3, 0.02, 0.3])
    i = np.radians([28.5, 97.4396, 63.4, 140.0])
    z = zonal_mean_rates(a, e, i, argp=1.1, j3=0.0, j4=0.0)
    s = secular_rates(a, e, i)
    for got, ref in ((z.n, s.n), (z.raan_dot, s.raan_dot), (z.argp_dot, s.argp_dot)):
        np.testing.assert_allclose(got, ref, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(z.psi_dot, s.psi_dot, rtol=1e-12, atol=0.0)
    assert np.all(z.e_dot == 0.0) and np.all(z.i_dot == 0.0)


@pytest.mark.parametrize(
    ("j3", "j4", "raan_ref", "argp_ref", "u_ref"),
    [
        (0.0, 0.0, A5.raan_dot_j2_deg_day, A5.argp_dot_j2_deg_day, A5.u_dot_j2_deg_day),
        (0.0, J4, A5.raan_dot_j2_j4_deg_day, A5.argp_dot_j2_j4_deg_day, A5.u_dot_j2_j4_deg_day),
        (J3, J4, A5.raan_dot_j2_j3_j4_deg_day, A5.argp_dot_j2_j3_j4_deg_day,
         A5.u_dot_j2_j3_j4_deg_day),
    ],
)  # fmt: skip
def test_zonal_rates_match_a5(j3, j4, raan_ref, argp_ref, u_ref) -> None:
    r = zonal_mean_rates(A_DESIGN, A5.e, I_DESIGN, math.radians(A5.argp_deg), j3=j3, j4=j4)
    assert r.raan_dot * DEG_DAY == pytest.approx(raan_ref, rel=1e-6)
    # o A.5 imprime ω̇ com 5 decimais (4 algarismos no caso congelado): tolerância de arredondamento
    assert r.argp_dot * DEG_DAY == pytest.approx(argp_ref, abs=5e-6)
    assert (r.psi_dot + r.argp_dot) * DEG_DAY == pytest.approx(u_ref, abs=5e-5)


def test_zonal_rates_frozen_orbit_e_and_i_constant() -> None:
    # Em ω = 90° os termos de período longo em ė e i̇ (cos ω) se anulam: é o ponto fixo
    r = zonal_mean_rates(A_DESIGN, A5.e, I_DESIGN, math.pi / 2)
    assert abs(r.e_dot) < 1e-25 and abs(r.i_dot) < 1e-25
    r = zonal_mean_rates(A_DESIGN, A5.e, I_DESIGN, 0.0)
    assert r.e_dot != 0.0


def test_zonal_rates_reject_circular_with_j3() -> None:
    with pytest.raises(ValueError, match="e > 0"):
        zonal_mean_rates(A_DESIGN, 0.0, I_DESIGN, 0.0)
    zonal_mean_rates(A_DESIGN, 0.0, I_DESIGN, 0.0, j3=0.0)  # sem J3 não há 1/e


# --- Órbita congelada -------------------------------------------------------------------------


def test_frozen_eccentricity() -> None:
    e_eq33 = -J3 / (2.0 * 1.08262668e-3) * (R_EARTH / A_DESIGN) * math.sin(I_DESIGN)
    assert e_eq33 == pytest.approx(A5.frozen_e_eq33, rel=1e-5)
    e_j3 = frozen_eccentricity(A_DESIGN, I_DESIGN, j4=False)
    assert e_j3 == pytest.approx(A5.frozen_e_eq33, rel=1e-5)
    assert e_j3 == pytest.approx(A5.frozen_e_j2_j3, rel=1e-6)
    e_j4 = frozen_eccentricity(A_DESIGN, I_DESIGN)
    assert e_j4 == pytest.approx(A5.frozen_e_j2_j4, rel=1e-4)
    rates = zonal_mean_rates(A_DESIGN, e_j4, I_DESIGN, math.pi / 2)
    assert abs(rates.argp_dot * DEG_DAY) < 1e-9


# --- Propagador médio -------------------------------------------------------------------------


def test_mean_zonal_propagator_reduces_to_mean_j2() -> None:
    el = KeplerianElements(7000e3, 0.01, math.radians(98.0), 0.3, 1.2, 0.4)
    t = np.linspace(-86400.0, 5 * 86400.0, 41)
    mz = MeanZonalPropagator(el, EPOCH, j3=0.0, j4=0.0)
    mj2 = MeanJ2Propagator(el, EPOCH)
    got, ref = mz.elements_at(t), mj2.elements_at(t)
    for field in ("a", "e", "i", "raan", "argp", "psi"):
        np.testing.assert_allclose(getattr(got, field), getattr(ref, field), rtol=1e-10, atol=1e-9)
    r, _ = mz.propagate(t)
    r_ref, _ = mj2.propagate(t)
    assert np.max(np.linalg.norm(r - r_ref, axis=1)) < 1e-2
    assert mz.period == pytest.approx(mj2.period, rel=1e-12)


def test_mean_zonal_propagator_frozen_orbit_stays_frozen() -> None:
    e_f = frozen_eccentricity(A_DESIGN, I_DESIGN)
    el = KeplerianElements(A_DESIGN, e_f, I_DESIGN, 0.0, math.pi / 2, 0.0)
    mz = MeanZonalPropagator(el, EPOCH)
    end = mz.elements_at(365 * 86400.0)
    assert end.e[0] == pytest.approx(e_f, rel=1e-9)
    assert end.argp[0] == pytest.approx(math.pi / 2, abs=1e-9)
    assert end.a[0] == A_DESIGN  # sem arrasto, a médio não muda


def test_mean_zonal_vs_cowell(cowell_30d) -> None:
    # Modelo médio J2–J4 de 1ª ordem nos elementos de projeto × Cowell J2–J4 (A.7)
    mean, t, raan_dot_cowell, u_dot_cowell = cowell_30d
    mz = MeanZonalPropagator(mean, EPOCH)
    el = mz.elements_at(t)
    raan_dot_mean = fit_rate(t, el.raan)
    ratio = raan_dot_mean / raan_dot_cowell - 1.0
    assert abs(ratio) <= A7.raan_dot_ratio_tol
    assert ratio == pytest.approx(A7.raan_dot_ratio_minus_1, abs=5e-5)
    tnod_cowell = 2.0 * math.pi / u_dot_cowell
    assert mz.period == pytest.approx(A7.tnod_theory_s, abs=1e-3)
    assert abs(mz.period / tnod_cowell - 1.0) * 1e6 <= A7.tnod_ppm_tol
    assert tnod_cowell == pytest.approx(A7.tnod_cowell_satsim_s, abs=2e-3)  # ver golden
    # O resíduo de 2ª ordem vira dezenas de km ao longo da órbita: modelo de projeto, não de posição
    u_dot_mean = fit_rate(t, el.argp + mean_to_true(el.psi, el.e))
    along_km = A_DESIGN * abs(u_dot_cowell - u_dot_mean) * t[-1] / 1e3
    assert 20.0 < along_km < 2 * A7.along_track_30d_km


@pytest.mark.parametrize(
    ("altitude_model", "cowell_ref"),
    [
        ("geodetic", A7.a_dot_cowell_geodetic_m_day),
        ("spherical", A7.a_dot_cowell_spherical_m_day),
    ],
)
def test_mean_drag_vs_cowell(i_calibrated, altitude_model, cowell_ref) -> None:
    drag = DragAcceleration(
        cd=2.2, area_over_mass=A7.cd_a_m / 2.2, altitude_model=altitude_model, corotating=False
    )
    mean = _calibrated_mean(i_calibrated)
    mz = MeanZonalPropagator(mean, EPOCH, drag=drag)
    t = np.linspace(0.0, 10 * 86400.0, 101)
    a_dot_mean = np.polyfit(t, mz.elements_at(t).a, 1)[0] * 86400.0
    # Cowell J2–J4 + arrasto a partir do estado convertido (conversão só com gravidade zonal)
    models = (*MODELS, drag)
    osc = mean_to_osculating(mean, EPOCH, conversion_models(models)).elements
    prop = CowellPropagator(osc, EPOCH, models)
    period = mz.period
    starts = np.arange(int(10 * 86400.0 // period)) * period
    a_mean = window_mean_elements(prop, starts, period).a
    a_dot_cowell = np.polyfit(starts + period / 2, a_mean, 1)[0] * 86400.0
    assert a_dot_cowell == pytest.approx(cowell_ref, rel=0.01)
    assert abs(a_dot_mean / a_dot_cowell - 1.0) <= A7.a_dot_ratio_tol
    # ⟨ρ⟩ ao longo da órbita × densidade em h = a − R⊕ (A.7: −29,69 m/d esférica)
    if altitude_model == "spherical":
        rho_h = drag.atmosphere.density(A_DESIGN - R_EARTH)
        a_dot_h = -A7.cd_a_m * rho_h * math.sqrt(MU * A_DESIGN) * 86400.0
        assert a_dot_h == pytest.approx(A7.a_dot_eq14_at_mean_altitude_m_day, rel=2e-3)


def test_mean_zonal_rejects_naive_epoch() -> None:
    el = KeplerianElements(A_DESIGN, A5.e, I_DESIGN, 0.0, math.pi / 2, 0.0)
    with pytest.raises(ValueError, match="tz-aware"):
        MeanZonalPropagator(el, datetime(2026, 10, 23))


# --- SSO calibrada e LTAN ---------------------------------------------------------------------


def test_raan_for_ltan_matches_a7() -> None:
    jd = datetime_to_jd(EPOCH)
    assert math.degrees(raan_for_ltan(A7.ltan_h, jd)) == pytest.approx(A7.raan_demo_deg, abs=1e-3)
    # LTAN 12h: o plano contém o Sol, Ω = α☉
    assert math.degrees(raan_for_ltan(12.0, jd)) == pytest.approx(A7.ra_sun_deg, abs=1e-3)
    raan = raan_for_ltan(np.array([0.0, 6.0, 18.0]), jd)
    assert np.shape(raan) == (3,) and np.all((raan >= 0.0) & (raan < 2 * math.pi))


def test_sso_inclination_calibrated(i_calibrated, cowell_30d) -> None:
    assert math.degrees(i_calibrated) == pytest.approx(
        A7.i_calibrated_deg, abs=A7.i_calibrated_tol_deg
    )
    # A fórmula de 1ª ordem só J2 fica ≈ 0,024° abaixo
    assert math.degrees(i_calibrated - sso_inclination(A_DESIGN, A7.e_frozen)) > 0.02
    _, _, raan_dot_cowell, _ = cowell_30d
    assert abs(ltan_drift_min_per_year(raan_dot_cowell)) < A7.ltan_drift_tol_min_year


def test_sso_inclination_calibrated_sentinel_2a() -> None:
    i_cal = sso_inclination_calibrated(
        A7.sentinel_a_m, A7.sentinel_e, EPOCH, MODELS, argp=math.radians(A7.sentinel_argp_deg)
    )
    i_j2 = sso_inclination(A7.sentinel_a_m, A7.sentinel_e)
    err_cal = math.degrees(i_cal) - A7.sentinel_i_real_deg
    err_j2 = math.degrees(i_j2) - A7.sentinel_i_real_deg
    assert abs(err_cal) < A7.sentinel_i_tol_deg
    assert abs(err_cal) < abs(err_j2)
    assert math.degrees(i_j2) == pytest.approx(A7.sentinel_i_first_order_j2_deg, abs=1e-4)
    assert math.degrees(i_cal) == pytest.approx(A7.sentinel_i_calibrated_deg, abs=1e-4)


def test_conversion_models_keeps_only_zonal_gravity() -> None:
    zonal = ZonalAcceleration()
    models = (
        zonal,
        DragAcceleration(),
        SolarRadiationPressure(EPOCH),
        ThirdBodyAcceleration("sun", EPOCH),
    )
    assert conversion_models(models) == (zonal,)
