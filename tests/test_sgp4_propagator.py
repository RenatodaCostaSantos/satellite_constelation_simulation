"""Testes do SGP4Propagator (S3-02) e de ``Propagator.propagate_utc``.

Referências: tests/golden/reference_values.py (ISS_TLE_A1, SGP4_REGRESSION_A1 e SGP4_SANITY_A1,
Apêndice A.1 do prompt da Semana 3: experimento independente com sgp4 2.27 e Skyfield 1.55).
"""

import logging
import math
from datetime import UTC, datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pytest
from golden.reference_values import ISS_TLE_A1, SGP4_REGRESSION_A1, SGP4_SANITY_A1

from satsim.app.scene_export import SceneSatellite, SceneStation, build_scene
from satsim.constants import R_EARTH, WGS84_F
from satsim.elements import KeplerianElements, state_to_elements
from satsim.frames import eci_to_ecef, gmst_at
from satsim.geodesy import ecef_to_spherical
from satsim.groundtrack import subpoint
from satsim.propagators import (
    CowellPropagator,
    KeplerPropagator,
    MeanJ2Propagator,
    Propagator,
    SGP4Propagator,
)
from satsim.tle import Tle

SAMPLE_TEXT = (Path(__file__).parent / "data" / "iss_sample.tle").read_text()


@pytest.fixture
def iss() -> SGP4Propagator:
    return SGP4Propagator.from_tle_text(SAMPLE_TEXT)


# --- Regressão (A.1) e sanidade física --------------------------------------------------------


def test_regression_vectors_a1(iss):
    t = np.array([row.t_h for row in SGP4_REGRESSION_A1]) * 3600.0
    r, v = iss.propagate(t)
    for k, row in enumerate(SGP4_REGRESSION_A1):
        np.testing.assert_allclose(r[k], np.array(row.r_km) * 1e3, rtol=0, atol=1.0)
        np.testing.assert_allclose(v[k], np.array(row.v_km_s) * 1e3, rtol=0, atol=1e-3)


def test_physical_sanity_at_epoch(iss):
    r, v = iss.propagate(0.0)
    assert np.linalg.norm(r[0]) == pytest.approx(SGP4_SANITY_A1.r0_km * 1e3, abs=5e3)
    assert np.linalg.norm(v[0]) == pytest.approx(SGP4_SANITY_A1.v0_km_s * 1e3, abs=10.0)
    osc = state_to_elements(r[0], v[0])
    assert math.degrees(osc.i) == pytest.approx(SGP4_SANITY_A1.osc_inclination_deg, abs=0.05)
    assert iss.period == pytest.approx(SGP4_SANITY_A1.period_s, abs=1.0)


def test_epoch_and_constructors(iss):
    assert iss.epoch == Tle.from_text(SAMPLE_TEXT).epoch
    assert abs((iss.epoch - datetime.fromisoformat(ISS_TLE_A1.epoch_utc_iso)).total_seconds()) < 1
    assert iss.tle.norad_id == ISS_TLE_A1.norad_id
    direct = SGP4Propagator(Tle.from_text(SAMPLE_TEXT))
    np.testing.assert_array_equal(direct.propagate(3600.0)[0], iss.propagate(3600.0)[0])


# --- Vetorização ------------------------------------------------------------------------------


def test_vectorized_equals_scalar_calls_any_order_and_negative(iss):
    t = np.array([86400.0, -3600.0, 0.0, 12345.6, -86400.0, 600.0])  # fora de ordem, negativos
    r, v = iss.propagate(t)
    assert r.shape == v.shape == (t.size, 3)
    for k, tk in enumerate(t):
        rk, vk = iss.propagate(tk)
        assert rk.shape == (1, 3)
        np.testing.assert_array_equal(r[k], rk[0])
        np.testing.assert_array_equal(v[k], vk[0])
    assert np.all(np.isfinite(r))
    assert iss.last_error_codes.tolist() == [0]  # da última chamada (escalar)


# --- propagate_utc (base.py) ------------------------------------------------------------------


def test_propagate_utc_matches_propagate(iss):
    r_utc, v_utc = iss.propagate_utc(iss.epoch + timedelta(hours=1), [0.0])
    r, v = iss.propagate([3600.0])
    np.testing.assert_allclose(r_utc, r, rtol=0, atol=1e-6)
    np.testing.assert_allclose(v_utc, v, rtol=0, atol=1e-9)


def test_propagate_utc_shifted_epoch_and_t_sum(iss):
    # epoch deslocado + t somam o mesmo instante, inclusive antes da época do TLE
    t = np.array([0.0, 600.0, 7200.0])
    for shift_s in (-5000.0, 0.0, 3600.0, 86400.0):
        r_a, _ = iss.propagate_utc(iss.epoch + timedelta(seconds=shift_s), t)
        r_b, _ = iss.propagate(t + shift_s)
        np.testing.assert_allclose(r_a, r_b, rtol=0, atol=1e-6)


def test_propagate_utc_other_timezone(iss):
    brt = datetime(2019, 12, 9, 13, 38, 29, tzinfo=timezone(timedelta(hours=-3)))
    r_a, _ = iss.propagate_utc(brt, [0.0])
    r_b, _ = iss.propagate_utc(brt.astimezone(UTC), [0.0])
    np.testing.assert_array_equal(r_a, r_b)


def test_propagate_utc_rejects_naive(iss):
    with pytest.raises(ValueError, match="tz-aware"):
        iss.propagate_utc(datetime(2019, 12, 10), [0.0])


def test_propagate_utc_cowell_before_epoch_raises(iss):
    r0, v0 = iss.propagate(0.0)
    cowell = CowellPropagator((r0[0], v0[0]), iss.epoch)
    cowell.propagate_utc(iss.epoch + timedelta(minutes=10), [0.0, 60.0])  # ok
    with pytest.raises(ValueError, match="t >= 0"):
        cowell.propagate_utc(iss.epoch - timedelta(seconds=1), [0.0])


# --- Erro do SGP4 simulado --------------------------------------------------------------------


class _FailingSatrec:
    """Envolve o Satrec real e marca com o código 6 (decaimento) os instantes de fr ≥ limite."""

    def __init__(self, real, fr_limit: float):
        self._real, self._fr_limit = real, fr_limit

    def __getattr__(self, name):
        return getattr(self._real, name)

    def sgp4_array(self, jd, fr):
        errors, r, v = self._real.sgp4_array(jd, fr)
        errors = np.where(fr >= self._fr_limit, 6, errors).astype(np.uint8)
        return errors, r, v


def test_sgp4_error_gives_nan_rows_and_single_warning(iss, monkeypatch, caplog):
    real = iss._satrec
    fr_limit = real.jdsatepochF + 2.0  # falha a partir de t = 2 dias
    monkeypatch.setattr(iss, "_satrec", _FailingSatrec(real, fr_limit))
    t = np.arange(0.0, 4 * 86400.0 + 1, 86400.0)  # 0, 1, 2, 3, 4 dias
    with caplog.at_level(logging.WARNING, logger="satsim.propagators.sgp4_ref"):
        r, v = iss.propagate(t)  # não levanta exceção
    failed = np.array([False, False, True, True, True])
    assert np.all(np.isnan(r[failed])) and np.all(np.isnan(v[failed]))
    assert np.all(np.isfinite(r[~failed])) and np.all(np.isfinite(v[~failed]))
    assert iss.last_error_codes.tolist() == [0, 0, 6, 6, 6]
    warnings = [rec for rec in caplog.records if rec.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "3 de 5" in warnings[0].message and "código 6" in warnings[0].message


# --- Skyfield (oráculo independente) ----------------------------------------------------------


def test_ecef_matches_skyfield_itrs():
    skyfield_api = pytest.importorskip("skyfield.api")
    framelib = pytest.importorskip("skyfield.framelib")
    ts = skyfield_api.load.timescale(builtin=True)
    sat = skyfield_api.EarthSatellite(ISS_TLE_A1.line1, ISS_TLE_A1.line2, ISS_TLE_A1.name, ts)
    iss = SGP4Propagator.from_tle_text(SAMPLE_TEXT)

    t = np.arange(0.0, 3 * 86400.0 + 1, 600.0)
    r_eci, _ = iss.propagate(t)
    r_ecef = eci_to_ecef(r_eci, gmst_at(iss.epoch, t))
    times = ts.from_datetimes([iss.epoch + timedelta(seconds=float(tk)) for tk in t])
    r_sky = sat.at(times).frame_xyz(framelib.itrs).m.T

    diff = np.linalg.norm(r_ecef - r_sky, axis=1)
    # critério do card: < 150 m; referência A.1: máx. 86 m, média 70 m (UT1 ≈ UTC, sem polo)
    assert diff.max() < 150.0
    assert diff.max() == pytest.approx(SGP4_SANITY_A1.skyfield_ecef_max_diff_m, abs=5.0)
    assert diff.mean() == pytest.approx(SGP4_SANITY_A1.skyfield_ecef_mean_diff_m, abs=5.0)


# --- Intercambialidade ------------------------------------------------------------------------


def test_subpoint_latitude_bounded_by_inclination(iss):
    t = np.arange(0.0, 86400.0 + 1, 60.0)
    lat, lon, alt = subpoint(iss, t)
    r_eci, _ = iss.propagate(t)
    lat_gc, _, _ = ecef_to_spherical(eci_to_ecef(r_eci, gmst_at(iss.epoch, t)))
    lat_gc_max = np.max(np.abs(np.degrees(lat_gc)))
    lat_gd_max = np.max(np.abs(np.degrees(lat)))
    # geocêntrica ≤ i (51,7°); geodésica = geocêntrica + f·sin 2φ·R⊕/(R⊕ + h) (ver golden)
    assert 51.0 < lat_gc_max <= SGP4_SANITY_A1.lat_max_geocentric_deg
    phi = math.radians(lat_gc_max)
    dphi = math.degrees(WGS84_F * math.sin(2 * phi) * R_EARTH / (R_EARTH + float(np.mean(alt))))
    assert lat_gd_max == pytest.approx(lat_gc_max + dphi, abs=0.01)
    assert np.all((alt > 380e3) & (alt < 460e3))  # faixa de altitude da ISS


def _iss_like_elements() -> KeplerianElements:
    return KeplerianElements(
        a=R_EARTH + 420e3, e=0.0007, i=math.radians(51.64), raan=1.0, argp=0.3, psi=0.0
    )


def _all_propagators(iss: SGP4Propagator) -> list[Propagator]:
    el = _iss_like_elements()
    return [
        KeplerPropagator(el, iss.epoch),
        MeanJ2Propagator(el, iss.epoch),
        CowellPropagator(el, iss.epoch),
        iss,
    ]


def test_same_code_runs_with_all_propagators(iss):
    """subpoint e build_scene rodam sem alteração com Kepler, MeanJ2, Cowell e SGP4."""
    props = _all_propagators(iss)
    t = np.arange(0.0, 3 * 3600.0 + 1, 60.0)
    for prop in props:
        lat, lon, alt = subpoint(prop, t)
        assert lat.shape == lon.shape == alt.shape == t.shape
        assert np.max(np.abs(np.degrees(lat))) <= 51.9
    sats = [SceneSatellite(f"s{k}", type(p).__name__, "#ffffff", p) for k, p in enumerate(props)]
    station = SceneStation("sbs", "São Bento do Sapucaí", -22.69, -45.73, 900.0)
    scene = build_scene(sats, [station], iss.epoch + timedelta(minutes=5), 3600.0, 60.0)
    assert len(scene["satellites"]) == 4
    for sat in scene["satellites"]:
        assert len(sat["lla"]) == 61
        assert np.all(np.isfinite(np.asarray(sat["eci_km"])))
