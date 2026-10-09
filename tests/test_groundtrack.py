"""Testes do subponto e do traço no solo (S2-04).

Referências: tests/golden/reference_values.py (DESIGN_ORBIT, Sem4 §7.1 e §7.5; GROUND_TRACK_S2,
card S2-04 do prompt da Semana 2).
"""

import math
from dataclasses import replace
from datetime import UTC, datetime

import numpy as np
import pytest
from golden.reference_values import DESIGN_ORBIT, GROUND_TRACK_S2

from satsim.constants import R_EARTH
from satsim.elements import KeplerianElements
from satsim.frames import eci_to_ecef, gmst_at
from satsim.geodesy import ecef_to_geodetic
from satsim.groundtrack import ascending_node_crossings, split_at_antimeridian, subpoint
from satsim.propagators import CowellPropagator, KeplerPropagator, MeanJ2Propagator

EPOCH = datetime(2026, 9, 29, 12, 0, 0, tzinfo=UTC)
T_1D = np.arange(0.0, 86400.0 + 1.0, 10.0)  # 24 h a 10 s

# SSO de projeto (DESIGN_ORBIT) e ISS-like (a = R⊕ + 420 km, i = 51,64°); Ω0 = 30° e ψ0 = 0.
SSO = KeplerianElements(
    a=DESIGN_ORBIT.a_km * 1e3,
    e=DESIGN_ORBIT.e,
    i=math.radians(DESIGN_ORBIT.i_deg),
    raan=math.radians(30.0),
    argp=math.radians(DESIGN_ORBIT.argp_deg),
    psi=0.0,
)
ISS_LIKE = KeplerianElements(
    a=R_EARTH + 420e3,
    e=0.0,
    i=math.radians(GROUND_TRACK_S2.iss_like_i_deg),
    raan=math.radians(30.0),
    argp=0.0,
    psi=0.0,
)


@pytest.fixture(scope="module")
def sso_track() -> tuple[np.ndarray, ...]:
    prop = MeanJ2Propagator(SSO, EPOCH)
    lat, lon, alt = subpoint(prop, T_1D)
    r, _ = prop.propagate(T_1D)
    return lat, lon, alt, r


# --- Subponto ----------------------------------------------------------------------------------


def test_subpoint_matches_frames_and_geodesy() -> None:
    prop = MeanJ2Propagator(SSO, EPOCH)
    t = np.array([0.0, 1234.0, 50000.0])
    lat, lon, alt = subpoint(prop, t)
    r, _ = prop.propagate(t)
    expected = ecef_to_geodetic(eci_to_ecef(r, gmst_at(EPOCH, t)))
    for got, exp in zip((lat, lon, alt), expected, strict=True):
        np.testing.assert_array_equal(got, exp)


def test_subpoint_shapes() -> None:
    prop = KeplerPropagator(SSO, EPOCH)
    assert all(x.shape == (1,) for x in subpoint(prop, 0.0))
    assert all(x.shape == (5,) for x in subpoint(prop, np.arange(5.0)))


@pytest.mark.parametrize(
    "make", [KeplerPropagator, MeanJ2Propagator, CowellPropagator], ids=lambda c: c.__name__
)
def test_subpoint_any_propagator(make) -> None:
    # Critério S2-04: mesma interface para Kepler, MeanJ2 e Cowell, sem alteração de código
    lat, lon, alt = subpoint(make(SSO, EPOCH), np.arange(0.0, 6000.0, 60.0))
    assert lat.shape == lon.shape == alt.shape == (100,)
    assert np.all(np.abs(lat) <= math.pi / 2) and np.all(np.abs(lon) <= math.pi)
    assert np.all((alt > 400e3) & (alt < 600e3))


# --- Latitude máxima ---------------------------------------------------------------------------


def test_sso_max_latitude(sso_track) -> None:
    # DESIGN_ORBIT (Sem4 eq. 12): max |arcsin(z/r)| = 180° − i = 82,56°;
    # GROUND_TRACK_S2: max |lat geodésica| ≈ 82,61°
    lat, _, _, r = sso_track
    lat_gc = np.arcsin(r[:, 2] / np.linalg.norm(r, axis=1))
    assert math.degrees(np.max(np.abs(lat_gc))) == pytest.approx(DESIGN_ORBIT.lat_max_deg, abs=0.01)
    assert math.degrees(np.max(np.abs(lat))) == pytest.approx(
        GROUND_TRACK_S2.sso_lat_max_geodetic_deg, abs=0.1
    )


def test_iss_like_max_latitude() -> None:
    # GROUND_TRACK_S2: latitude geocêntrica máxima = i = 51,64°
    r, _ = MeanJ2Propagator(ISS_LIKE, EPOCH).propagate(T_1D)
    lat_gc = np.arcsin(r[:, 2] / np.linalg.norm(r, axis=1))
    assert math.degrees(np.max(np.abs(lat_gc))) == pytest.approx(
        GROUND_TRACK_S2.iss_like_i_deg, abs=0.01
    )


# --- Nós ascendentes ---------------------------------------------------------------------------


def test_sso_ascending_nodes(sso_track) -> None:
    # DESIGN_ORBIT (Sem4 §7.1, §7.5): Δλ por revolução = −23,736° e Tnod = 5696,70 s
    lat, lon, _, _ = sso_track
    t_nodes, lon_nodes = ascending_node_crossings(T_1D, lat, lon)
    assert t_nodes.size == 15
    dlon = np.degrees(np.diff(lon_nodes))
    dlon = (dlon + 180.0) % 360.0 - 180.0  # módulo 360°
    np.testing.assert_allclose(dlon, DESIGN_ORBIT.dlambda_per_rev_deg, rtol=0.0, atol=0.01)
    np.testing.assert_allclose(np.diff(t_nodes), DESIGN_ORBIT.tnod_s, rtol=0.0, atol=0.05)


def test_91_6_cycle_closure() -> None:
    # Critério S2-04: após 6 dias (91 revoluções nodais), o nó k+91 coincide com o nó k a menos
    # de 0,01° (GROUND_TRACK_S2: ≈ 1,3e-4°). a é o valor do Sem4; o solver de repetição é da S5.
    assert DESIGN_ORBIT.n_revs == 91 and DESIGN_ORBIT.d_days == 6
    t = np.arange(0.0, 6 * 86400.0 + 2 * DESIGN_ORBIT.tnod_s, 10.0)
    lat, lon, _ = subpoint(MeanJ2Propagator(SSO, EPOCH), t)
    _, lon_nodes = ascending_node_crossings(t, lat, lon)
    assert lon_nodes.size > 91
    gap = math.degrees(lon_nodes[91] - lon_nodes[0])
    gap = (gap + 180.0) % 360.0 - 180.0
    assert abs(gap) < 0.01


def test_ascending_nodes_synthetic() -> None:
    # lat = sin(2πt/100): sobe por 0 em t = 0, 100, 200; lon cruza o antimeridiano em t = 100
    t = np.arange(-5.0, 250.0, 7.0)
    lat = np.sin(2 * np.pi * t / 100.0)
    lon = np.angle(np.exp(1j * (np.pi - 0.5 + 0.01 * t)))  # passa por +π e volta a −π
    t_nodes, lon_nodes = ascending_node_crossings(t, lat, lon)
    np.testing.assert_allclose(t_nodes, [0.0, 100.0, 200.0], atol=0.3)
    expected = np.angle(np.exp(1j * (np.pi - 0.5 + 0.01 * t_nodes)))
    np.testing.assert_allclose(lon_nodes, expected, atol=1e-12)


# --- Antimeridiano -----------------------------------------------------------------------------


def test_split_at_antimeridian(sso_track) -> None:
    # Critério S2-04: nenhum segmento com |Δλ| > 180°; nenhum ponto criado nem descartado
    lat, lon, _, _ = sso_track
    segments = split_at_antimeridian(lon, lat)
    for seg_lon, seg_lat in segments:
        assert seg_lon.shape == seg_lat.shape
        assert np.all(np.abs(np.diff(seg_lon)) <= math.pi)
    assert sum(seg_lon.size for seg_lon, _ in segments) == lon.size
    np.testing.assert_array_equal(np.concatenate([s for s, _ in segments]), lon)


def _sso_crossings(raan_deg: float, duration_s: float, step_s: float) -> int:
    el = replace(SSO, raan=math.radians(raan_deg))
    t = np.arange(0.0, duration_s + step_s / 2, step_s)
    lat, lon, _ = subpoint(MeanJ2Propagator(el, EPOCH), t)
    return len(split_at_antimeridian(lon, lat)) - 1


@pytest.mark.parametrize("raan_deg", [0.0, 30.0, 100.0, 150.0, 250.0])
def test_sso_antimeridian_crossings_24h(raan_deg: float) -> None:
    # GROUND_TRACK_S2 (corrigido na revisão da S2-04): SSO retrógrada → 16,17 cruzamentos/dia,
    # 16 ou 17 em 24 h conforme Ω0 (o card dizia {15, 16})
    n = _sso_crossings(raan_deg, 86400.0, 10.0)
    assert n in GROUND_TRACK_S2.sso_antimeridian_crossings_24h


def test_sso_antimeridian_crossings_91_6_cycle() -> None:
    # GROUND_TRACK_S2: no ciclo de 6 dias, (360° + 23,736°)/360° × 91 = 97,0 cruzamentos
    rate = (360.0 - DESIGN_ORBIT.dlambda_per_rev_deg) / 360.0 * 86400.0 / DESIGN_ORBIT.tnod_s
    assert rate == pytest.approx(GROUND_TRACK_S2.sso_antimeridian_crossings_per_day, abs=0.005)
    n = _sso_crossings(30.0, DESIGN_ORBIT.d_days * 86400.0, 30.0)
    assert n == GROUND_TRACK_S2.sso_antimeridian_crossings_91_6_cycle


def test_split_without_crossing_is_single_segment() -> None:
    lon = np.linspace(-1.0, 1.0, 20)
    segments = split_at_antimeridian(lon, np.zeros(20))
    assert len(segments) == 1
    np.testing.assert_array_equal(segments[0][0], lon)
