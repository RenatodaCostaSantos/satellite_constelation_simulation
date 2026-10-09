"""Testes dos quadros ECI ↔ ECEF e da geodésia WGS84 (S2-03).

Referências: tests/golden/reference_values.py (WGS84_REF, Apêndice A.3 do prompt da Semana 2;
NIMA TR8350.2). Skyfield (oráculo) só na parte ECEF ↔ geodésico: o ECI do satsim é o equatorial
da data sem precessão/nutação e não é comparável ao GCRS do Skyfield.
"""

import math
from datetime import UTC, datetime

import numpy as np
import pytest
from golden.reference_values import WGS84_REF

from satsim.astrotime import datetime_to_jd, gmst
from satsim.constants import OMEGA_EARTH, WGS84_A, WGS84_B
from satsim.frames import ecef_to_eci, eci_to_ecef, gmst_at, rot_z
from satsim.geodesy import ecef_to_geodetic, ecef_to_spherical, geodetic_to_ecef

EPOCH = datetime(2026, 9, 29, 12, 0, 0, tzinfo=UTC)
RNG = np.random.default_rng(3)

# --- ECI ↔ ECEF --------------------------------------------------------------------------------


def test_rotation_sign() -> None:
    # A.3: r_eci = (1, 0, 0) com θ = 90° → longitude −90° no ECEF (a Terra girou 90° para leste)
    r_ecef = eci_to_ecef([1.0, 0.0, 0.0], math.pi / 2)
    np.testing.assert_allclose(r_ecef, [0.0, -1.0, 0.0], atol=1e-15)
    assert math.degrees(math.atan2(r_ecef[1], r_ecef[0])) == pytest.approx(-90.0, abs=1e-12)


def test_rot_z_matches_function_and_is_orthogonal() -> None:
    theta = RNG.uniform(0.0, 2 * np.pi, 7)
    mats = rot_z(theta)
    assert mats.shape == (7, 3, 3) and rot_z(0.3).shape == (3, 3)
    r = RNG.normal(size=(7, 3)) * 7e6
    np.testing.assert_allclose(np.einsum("nij,nj->ni", mats, r), eci_to_ecef(r, theta), atol=1e-8)
    identity = np.broadcast_to(np.eye(3), mats.shape)
    np.testing.assert_allclose(mats @ mats.transpose(0, 2, 1), identity, atol=1e-15)
    np.testing.assert_allclose(rot_z(math.pi / 2), [[0, 1, 0], [-1, 0, 0], [0, 0, 1]], atol=1e-16)


def test_eci_ecef_roundtrip() -> None:
    # Critério S2-03 "< 1e-9 m": abaixo da resolução do float64 em LEO (ulp de 7e6 m = 0,93 nm);
    # testado como erro relativo de máquina, |Δr| < 1e-15·|r| (≈ 7 nm em 7000 km).
    n = 10_000
    u = RNG.normal(size=(n, 3))
    r = u / np.linalg.norm(u, axis=1, keepdims=True) * RNG.uniform(6.6e6, 4.2e7, n)[:, None]
    v = RNG.normal(size=(n, 3)) * 7.5e3
    theta = RNG.uniform(0.0, 2 * np.pi, n)
    r_ecef, v_ecef = eci_to_ecef(r, theta, v)
    r_back, v_back = ecef_to_eci(r_ecef, theta, v_ecef)
    r_norm = np.linalg.norm(r, axis=1)
    assert np.all(np.linalg.norm(r_back - r, axis=1) < 1e-15 * r_norm)
    assert np.max(np.linalg.norm(v_back - v, axis=1)) < 1e-9


def test_ground_point_velocity() -> None:
    # Ponto fixo no solo: ECEF constante, velocidade ECI = ω⊕ × r; de volta, v_ecef = 0
    r_ecef = geodetic_to_ecef(
        math.radians(WGS84_REF.sbs_lat_deg),
        math.radians(WGS84_REF.sbs_lon_deg),
        WGS84_REF.sbs_alt_m,
    )
    theta = 1.234
    r_eci, v_eci = ecef_to_eci(r_ecef, theta, np.zeros(3))
    np.testing.assert_allclose(v_eci, np.cross([0.0, 0.0, OMEGA_EARTH], r_eci), atol=1e-12)
    _, v_ecef = eci_to_ecef(r_eci, theta, v_eci)
    np.testing.assert_allclose(v_ecef, np.zeros(3), atol=1e-9)


def test_gmst_at_matches_astrotime() -> None:
    t = np.array([0.0, 3600.0, 86400.0])
    expected = gmst(datetime_to_jd(EPOCH) + t / 86400.0)
    np.testing.assert_allclose(gmst_at(EPOCH, t), expected, rtol=0.0, atol=1e-12)
    assert isinstance(gmst_at(EPOCH, 0.0), float)


def test_frames_shapes() -> None:
    assert eci_to_ecef(np.ones(3), 0.5).shape == (3,)
    r_ecef, v_ecef = eci_to_ecef(np.ones((4, 3)), np.zeros(4), np.ones((4, 3)))
    assert r_ecef.shape == (4, 3) and v_ecef.shape == (4, 3)
    assert ecef_to_eci(np.ones((4, 3)), 0.5).shape == (4, 3)  # θ escalar para todas as linhas


# --- Geodésia ----------------------------------------------------------------------------------


def test_known_values() -> None:
    # A.3: (a, 0, 0) → (0°, 0°, 0); (0, 0, b) → (90°, ·, |h| < 1e-6 m)
    lat, lon, h = ecef_to_geodetic([WGS84_A, 0.0, 0.0])
    assert (lat, lon) == (0.0, 0.0) and abs(h) < 1e-6
    lat, _, h = ecef_to_geodetic([0.0, 0.0, WGS84_B])
    assert lat == pytest.approx(math.pi / 2, abs=1e-15) and abs(h) < 1e-6
    lat, _, h = ecef_to_geodetic([0.0, 0.0, -WGS84_B - 500e3])
    assert lat == pytest.approx(-math.pi / 2, abs=1e-15)
    assert h == pytest.approx(500e3, abs=1e-6)
    np.testing.assert_allclose(geodetic_to_ecef(math.pi / 2, 0.0, 0.0), [0, 0, WGS84_B], atol=1e-6)
    assert WGS84_B == pytest.approx(WGS84_REF.b_m, abs=1e-6)
    assert (WGS84_A - WGS84_B) / 1e3 == pytest.approx(WGS84_REF.equator_minus_pole_km, abs=1e-3)


def test_geocentric_minus_geodetic_at_45() -> None:
    # A.3: latitude geocêntrica − geodésica a 45° (h = 0) = −0,1924°
    lat_gc, _, _ = ecef_to_spherical(geodetic_to_ecef(math.radians(45.0), 0.0, 0.0))
    assert math.degrees(lat_gc) - 45.0 == pytest.approx(
        WGS84_REF.geocentric_minus_geodetic_45_deg, abs=1e-3
    )


def _sample_geodetic(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """n pontos uniformes na esfera (lat, lon), h ∈ [−10 km, 40 000 km], mais casos de contorno."""
    lat = np.arcsin(RNG.uniform(-1.0, 1.0, n))
    lon = RNG.uniform(-np.pi, np.pi, n)
    h = RNG.uniform(-10e3, 40_000e3, n)
    edge_lat = [np.pi / 2, -np.pi / 2, 0.0, 0.0, np.pi / 2, -np.pi / 2, 0.0, np.pi / 2]
    edge_lon = [0.0, 1.0, 0.0, 2.0, -3.0, 0.5, np.pi, 0.0]
    edge_h = [0.0, 0.0, 0.0, -10e3, 40_000e3, -10e3, 0.0, 400e3]
    return np.r_[lat, edge_lat], np.r_[lon, edge_lon], np.r_[h, edge_h]


def test_geodetic_roundtrip() -> None:
    # Critério S2-03: ecef → geo → ecef < 1 mm; geo → ecef → geo |Δh| < 1 mm, |Δlat| < 1e-10 rad,
    # inclusive polos e equador (10 000 pontos + contorno)
    lat, lon, h = _sample_geodetic(10_000)
    r_ecef = geodetic_to_ecef(lat, lon, h)
    lat2, lon2, h2 = ecef_to_geodetic(r_ecef)
    assert np.max(np.abs(lat2 - lat)) < 1e-10
    assert np.max(np.abs(h2 - h)) < 1e-3
    off_pole = np.abs(lat) < np.pi / 2 - 1e-6  # nos polos a longitude é indefinida
    dlon = np.angle(np.exp(1j * (lon2 - lon)))[off_pole]
    assert np.max(np.abs(dlon)) < 1e-12
    r_back = geodetic_to_ecef(lat2, lon2, h2)
    assert np.max(np.linalg.norm(r_back - r_ecef, axis=1)) < 1e-3


def test_geodesy_shapes() -> None:
    r = geodetic_to_ecef(0.1, 0.2, 300.0)
    assert r.shape == (3,)
    out = ecef_to_geodetic(r)
    assert all(isinstance(x, float) for x in out)
    assert all(isinstance(x, float) for x in ecef_to_spherical(r))
    lat, lon, h = (np.full(5, 0.1), np.linspace(-1, 1, 5), np.zeros(5))
    r = geodetic_to_ecef(lat, lon, h)
    assert r.shape == (5, 3)
    assert all(x.shape == (5,) for x in ecef_to_geodetic(r))
    assert all(x.shape == (5,) for x in ecef_to_spherical(r))
    assert geodetic_to_ecef(lat, lon, 0.0).shape == (5, 3)  # broadcasting de h escalar


def test_spherical_approximation() -> None:
    r_ecef = np.array([WGS84_A + 500e3, 0.0, 0.0])
    lat, lon, alt = ecef_to_spherical(r_ecef)
    assert (lat, lon) == (0.0, 0.0) and alt == pytest.approx(500e3, abs=1e-6)


# --- Skyfield (oráculo, só ECEF ↔ geodésico) ---------------------------------------------------


def _skyfield_points() -> list[tuple[float, float, float]]:
    """20 pontos (lat°, lon°, elevação m), incluindo São Bento do Sapucaí (A.3)."""
    pts = [(WGS84_REF.sbs_lat_deg, WGS84_REF.sbs_lon_deg, WGS84_REF.sbs_alt_m)]
    pts += [(0.0, 0.0, 0.0), (89.9, 10.0, 0.0), (-89.9, -170.0, 1000.0), (45.0, 179.9, 500e3)]
    lat = np.degrees(np.arcsin(RNG.uniform(-1.0, 1.0, 15)))
    lon = RNG.uniform(-180.0, 180.0, 15)
    h = RNG.uniform(-100.0, 1000e3, 15)
    pts += list(zip(lat.tolist(), lon.tolist(), h.tolist(), strict=True))
    return pts


def test_geodetic_to_ecef_vs_skyfield() -> None:
    sf = pytest.importorskip("skyfield.api")
    for lat_deg, lon_deg, h_m in _skyfield_points():
        expected = sf.wgs84.latlon(lat_deg, lon_deg, elevation_m=h_m).itrs_xyz.m
        got = geodetic_to_ecef(math.radians(lat_deg), math.radians(lon_deg), h_m)
        assert np.linalg.norm(got - expected) < 10.0


def test_ecef_to_geodetic_vs_skyfield() -> None:
    sf = pytest.importorskip("skyfield.api")
    toposlib = pytest.importorskip("skyfield.toposlib")
    units = pytest.importorskip("skyfield.units")
    # ITRS → ITRS não envolve rotação; o instante só é exigido pela API (.at(t))
    t = sf.load.timescale(builtin=True).utc(2026, 9, 29, 12)
    for lat_deg, lon_deg, h_m in _skyfield_points():
        r = geodetic_to_ecef(math.radians(lat_deg), math.radians(lon_deg), h_m)
        pos = toposlib.ITRSPosition(units.Distance(m=r)).at(t)
        sf_lat, sf_lon = sf.wgs84.latlon_of(pos)
        sf_h = sf.wgs84.height_of(pos).m
        lat, lon, h = ecef_to_geodetic(r)
        back = geodetic_to_ecef(sf_lat.radians, sf_lon.radians, sf_h)
        assert np.linalg.norm(back - r) < 10.0  # posição equivalente à do Skyfield < 10 m
        assert abs(h - sf_h) < 10.0
        assert abs(math.degrees(lat) - sf_lat.degrees) < 1e-4
