"""Testes da estação de solo e das coordenadas topocêntricas (S3-04).

Referências: tests/golden/reference_values.py (STATION_A6, FOOTPRINT_A6 e ISS_TLE_A1, Apêndices
A.1, A.2 e A.6 do prompt da Semana 3; Sem4 §9.4 para o footprint).
"""

import math
from datetime import timedelta
from pathlib import Path

import numpy as np
import pytest
from golden.reference_values import FOOTPRINT_A6, ISS_TLE_A1, STATION_A6

from satsim.constants import R_EARTH, WGS84_A
from satsim.elements import KeplerianElements
from satsim.geodesy import geodetic_to_ecef
from satsim.observer import (
    SAO_BENTO_DO_SAPUCAI,
    GroundStation,
    Observation,
    footprint_half_angle,
    observe,
)
from satsim.propagators import (
    CowellPropagator,
    KeplerPropagator,
    MeanJ2Propagator,
    SGP4Propagator,
)

SAMPLE_TEXT = (Path(__file__).parent / "data" / "iss_sample.tle").read_text()
ORIGIN = GroundStation.from_degrees("origem", 0.0, 0.0, 0.0)  # lat 0, lon 0, h 0
ABS_DEG = 1e-6


def _deg(x):
    return math.degrees(float(x))


# --- Geometria exata (estação em lat 0°, lon 0°, h = 0) ----------------------------------------
# Nessa estação, U = +x, E = +y e N = +z em ECEF.


def test_zenith_satellite():
    az, el, rng = ORIGIN.look_angles([WGS84_A + 500e3, 0.0, 0.0])
    assert _deg(el) == pytest.approx(90.0, abs=ABS_DEG)
    assert rng == pytest.approx(500e3, abs=1e-6)


@pytest.mark.parametrize(
    ("offset", "az_deg"),
    [((0, 0, 1), 0.0), ((0, 1, 0), 90.0), ((0, 0, -1), 180.0), ((0, -1, 0), 270.0)],
    ids=["norte", "leste", "sul", "oeste"],
)
def test_azimuth_convention(offset, az_deg):
    # satélite 500 km acima, deslocado 200 km para norte, leste, sul ou oeste
    r = np.array([WGS84_A + 500e3, 0.0, 0.0]) + 200e3 * np.array(offset, dtype=float)
    az, el, _ = ORIGIN.look_angles(r)
    assert _deg(az) == pytest.approx(az_deg, abs=ABS_DEG)
    assert 0.0 < _deg(el) < 90.0


def test_horizon_and_footprint_half_angle():
    h = 400e3
    lam0 = math.acos(R_EARTH / (R_EARTH + h))
    assert math.degrees(lam0) == pytest.approx(FOOTPRINT_A6.h400_el0_deg, abs=1e-3)
    assert footprint_half_angle(h, 0.0) == pytest.approx(lam0, abs=1e-15)

    def sat_at(lon):  # satélite no equador, a R⊕ + h, deslocado em longitude
        return (R_EARTH + h) * np.array([math.cos(lon), math.sin(lon), 0.0])

    _, el, _ = ORIGIN.look_angles(sat_at(lam0))
    assert _deg(el) == pytest.approx(0.0, abs=ABS_DEG)
    eps = math.radians(10.0)
    _, el, _ = ORIGIN.look_angles(sat_at(footprint_half_angle(h, eps)))
    assert _deg(el) == pytest.approx(10.0, abs=ABS_DEG)
    _, el, _ = ORIGIN.look_angles(sat_at(math.pi))  # do outro lado da Terra
    assert el < 0.0


def test_footprint_iss_typical_value():
    lam = footprint_half_angle(420e3, math.radians(10.0))
    assert math.degrees(lam) == pytest.approx(FOOTPRINT_A6.h420_el10_deg, abs=0.1)
    lam_arr = footprint_half_angle(np.array([400e3, 420e3]), math.radians(10.0))
    assert lam_arr.shape == (2,) and isinstance(lam, float)


# --- Matriz ENU e estação padrão ---------------------------------------------------------------


@pytest.mark.parametrize("station", [ORIGIN, SAO_BENTO_DO_SAPUCAI,
                                     GroundStation.from_degrees("polo", 89.9, 120.0, 3000.0)],
                         ids=["origem", "sbs", "polo"])  # fmt: skip
def test_enu_matrix_is_rotation(station):
    m = station.enu_matrix()
    np.testing.assert_allclose(m @ m.T, np.eye(3), atol=1e-15)
    assert np.linalg.det(m) == pytest.approx(1.0, abs=1e-15)
    np.testing.assert_allclose(np.cross(m[0], m[1]), m[2], atol=1e-15)  # E × N = U


def test_sao_bento_station_ecef_a2():
    s = SAO_BENTO_DO_SAPUCAI
    assert (s.lat_deg, s.lon_deg, s.alt) == pytest.approx(
        (STATION_A6.lat_deg, STATION_A6.lon_deg, STATION_A6.alt_m), abs=1e-12
    )
    np.testing.assert_allclose(s.ecef(), geodetic_to_ecef(s.lat, s.lon, s.alt), rtol=0, atol=1e-6)
    np.testing.assert_allclose(s.ecef(), STATION_A6.ecef_m, rtol=0, atol=0.1)


def test_station_itself_is_origin_of_enu():
    s = SAO_BENTO_DO_SAPUCAI
    np.testing.assert_allclose(s.enu(s.ecef()), 0.0, atol=1e-9)
    # um ponto 1 km acima ao longo da normal está no zênite
    _, el, rng = s.look_angles(s.ecef() + 1e3 * s.enu_matrix()[2])
    assert _deg(el) == pytest.approx(90.0, abs=1e-6) and rng == pytest.approx(1e3, abs=1e-9)


# --- Vetorização -------------------------------------------------------------------------------


def test_look_angles_vectorized():
    rng_gen = np.random.default_rng(3)
    r = rng_gen.normal(size=(20, 3)) * 7e6
    az, el, d = SAO_BENTO_DO_SAPUCAI.look_angles(r)
    assert az.shape == el.shape == d.shape == (20,)
    assert np.all((az >= 0.0) & (az < 2 * math.pi))
    for k in range(20):
        az_k, el_k, d_k = SAO_BENTO_DO_SAPUCAI.look_angles(r[k])
        assert isinstance(az_k, float)
        assert (az_k, el_k, d_k) == pytest.approx((az[k], el[k], d[k]), rel=1e-14)


# --- observe com qualquer propagador -----------------------------------------------------------


@pytest.fixture(scope="module")
def iss() -> SGP4Propagator:
    return SGP4Propagator.from_tle_text(SAMPLE_TEXT)


def test_observe_runs_with_all_propagators(iss):
    el0 = KeplerianElements(a=R_EARTH + 420e3, e=0.0007, i=math.radians(51.64), raan=1.0,
                            argp=0.3, psi=0.0)  # fmt: skip
    epoch = iss.epoch
    t = np.arange(0.0, 6 * 3600.0 + 1, 60.0)
    for prop in (KeplerPropagator(el0, epoch), MeanJ2Propagator(el0, epoch),
                 CowellPropagator(el0, epoch), iss):  # fmt: skip
        obs = observe(prop, SAO_BENTO_DO_SAPUCAI, t)
        assert isinstance(obs, Observation)
        assert obs.az.shape == obs.el.shape == obs.range.shape == t.shape
        assert np.all(obs.range > 300e3)
    scalar = observe(iss, SAO_BENTO_DO_SAPUCAI, 600.0)
    assert scalar.el.shape == (1,)


def test_observe_matches_skyfield_altaz(iss):
    skyfield_api = pytest.importorskip("skyfield.api")
    ts = skyfield_api.load.timescale(builtin=True)
    sat = skyfield_api.EarthSatellite(ISS_TLE_A1.line1, ISS_TLE_A1.line2, ISS_TLE_A1.name, ts)
    topos = skyfield_api.wgs84.latlon(STATION_A6.lat_deg, STATION_A6.lon_deg,
                                      elevation_m=STATION_A6.alt_m)  # fmt: skip

    t = np.arange(0.0, 3 * 86400.0 + 1, 30.0)
    obs = observe(iss, SAO_BENTO_DO_SAPUCAI, t)
    visible = obs.el > 0.0
    assert visible.sum() > 100
    times = ts.from_datetimes([iss.epoch + timedelta(seconds=float(x)) for x in t[visible]])
    alt, az, dist = (sat - topos).at(times).altaz()

    daz = (obs.az_deg[visible] - az.degrees + 180.0) % 360.0 - 180.0
    del_ = obs.el_deg[visible] - alt.degrees
    drange = obs.range[visible] - dist.m
    # critérios do card: 0,05°, 0,05° e 150 m; referência A.2: 0,026°, 0,009° e 74 m
    assert np.max(np.abs(daz)) < 0.05
    assert np.max(np.abs(del_)) < 0.05
    assert np.max(np.abs(drange)) < 150.0
    assert np.max(np.abs(daz)) == pytest.approx(STATION_A6.skyfield_max_daz_deg, abs=0.003)
    assert np.max(np.abs(del_)) == pytest.approx(STATION_A6.skyfield_max_del_deg, abs=0.002)
    assert np.max(np.abs(drange)) == pytest.approx(STATION_A6.skyfield_max_drange_m, abs=5.0)
