"""Testes da exportação da cena para o viewer (S1-08)."""

import json
import math
from datetime import UTC, datetime, timedelta

import numpy as np
import pytest
from golden.reference_values import DESIGN_ORBIT

from satsim.app.scene_export import (
    SCHEMA_VERSION,
    SceneSatellite,
    SceneStation,
    build_scene,
    eci_to_ecef,
    subpoint_spherical,
)
from satsim.constants import R_EARTH
from satsim.elements import KeplerianElements
from satsim.propagators import KeplerPropagator

EPOCH = datetime(2026, 9, 29, 12, 0, 0, tzinfo=UTC)
DURATION_S = 6 * 3600
STEP_S = 10
STATION = SceneStation("sbs", "São Bento do Sapucaí", -22.69, -45.73, 900.0)

# Órbita de projeto 91/6 (Sem4 §7.5, golden DESIGN_ORBIT)
D = DESIGN_ORBIT
SSO_ELEMENTS = KeplerianElements(
    a=D.a_km * 1e3,
    e=D.e,
    i=math.radians(D.i_deg),
    raan=math.radians(30.0),
    argp=math.radians(D.argp_deg),
    psi=0.0,
)


def _sso(epoch: datetime = EPOCH, sat_id: str = "sso1") -> SceneSatellite:
    return SceneSatellite(sat_id, "SSO-1", "#00e5ff", KeplerPropagator(SSO_ELEMENTS, epoch))


def _example(i_deg: float = 51.64, sat_id: str = "kep1") -> SceneSatellite:
    el = KeplerianElements(6798e3, 0.0, math.radians(i_deg), math.radians(120.0), 0.0, 0.0)
    return SceneSatellite(sat_id, "Exemplo kepleriano", "#ffb020", KeplerPropagator(el, EPOCH))


@pytest.fixture(scope="module")
def scene() -> dict:
    return build_scene([_sso(), _example()], [STATION], EPOCH, DURATION_S, STEP_S)


# --- Estrutura do JSON (contrato do ADR 0001) --------------------------------------------------


def test_structure(scene: dict) -> None:
    assert set(scene) == {"meta", "earth", "satellites", "stations", "passes"}
    assert scene["meta"]["schema_version"] == SCHEMA_VERSION
    assert scene["meta"]["epoch_utc"] == "2026-09-29T12:00:00Z"
    assert scene["meta"]["duration_s"] == DURATION_S and scene["meta"]["step_s"] == STEP_S
    assert set(scene["earth"]) == {"gmst0_rad", "gmst_rate_rad_s"}
    assert scene["passes"] == []
    assert scene["stations"] == [
        {
            "id": "sbs",
            "name": "São Bento do Sapucaí",
            "lat_deg": -22.69,
            "lon_deg": -45.73,
            "alt_m": 900.0,
        }
    ]
    n = DURATION_S // STEP_S + 1
    for sat in scene["satellites"]:
        assert set(sat) == {"id", "name", "color", "eci_km", "lla"}
        assert len(sat["eci_km"]) == n and len(sat["lla"]) == n
        assert all(len(p) == 3 for p in sat["eci_km"]) and all(len(p) == 3 for p in sat["lla"])


def test_json_serializable(scene: dict) -> None:
    text = json.dumps(scene)
    assert json.loads(text) == scene


# --- SSO de projeto ----------------------------------------------------------------------------


def test_sso_latitude_and_longitude_range(scene: dict) -> None:
    lla = np.array(scene["satellites"][0]["lla"])
    # Sem4 eq. 12 / golden: latitude máxima = 180° − i = 82,56°
    assert np.max(np.abs(lla[:, 0])) <= D.lat_max_deg + 0.05
    assert np.max(lla[:, 0]) == pytest.approx(D.lat_max_deg, abs=0.1)
    assert np.all((lla[:, 1] >= -180.0) & (lla[:, 1] < 180.0))


def test_sso_radius_and_altitude(scene: dict) -> None:
    eci = np.array(scene["satellites"][0]["eci_km"])
    lla = np.array(scene["satellites"][0]["lla"])
    r_km = np.linalg.norm(eci, axis=1)
    a_km, e = D.a_km, D.e
    tol = 1e-5  # arredondamento de eci_km (1e-6 km por componente)
    assert np.all(r_km >= a_km * (1 - e) - tol) and np.all(r_km <= a_km * (1 + e) + tol)
    np.testing.assert_allclose(lla[:, 2], r_km - R_EARTH / 1e3, atol=2e-4)


# --- Consistência: lla recalculado a partir de eci_km e do GMST da cena -------------------------


def test_lla_consistent_with_eci_and_gmst(scene: dict) -> None:
    earth = scene["earth"]
    t = np.arange(DURATION_S // STEP_S + 1) * STEP_S
    theta = earth["gmst0_rad"] + earth["gmst_rate_rad_s"] * t  # modelo linear do contrato
    for sat in scene["satellites"]:
        eci_km = np.array(sat["eci_km"])
        lla = np.array(sat["lla"])
        recomputed = subpoint_spherical(eci_to_ecef(eci_km * 1e3, theta))
        np.testing.assert_allclose(recomputed[:, 0], lla[:, 0], atol=1e-6)
        dlon = (recomputed[:, 1] - lla[:, 1] + 180.0) % 360.0 - 180.0
        assert np.max(np.abs(dlon)) < 1e-6


def test_eci_to_ecef_rotation() -> None:
    # θ = 90°: o eixo x do ECI aparece em −y no ECEF (a Terra girou 90° para leste)
    r = eci_to_ecef(np.array([[1.0, 0.0, 0.0]]), np.array([math.pi / 2]))
    np.testing.assert_allclose(r, [[0.0, -1.0, 0.0]], atol=1e-15)


# --- Multi-satélite (R7) -----------------------------------------------------------------------


@pytest.mark.parametrize("n_sats", [1, 2, 3])
def test_multi_satellite(n_sats: int) -> None:
    sats = [_sso(), _example(), _example(i_deg=30.0, sat_id="kep2")][:n_sats]
    scene = build_scene(sats, [STATION], EPOCH, 3600, STEP_S)
    assert [s["id"] for s in scene["satellites"]] == [s.id for s in sats]
    assert all(len(s["lla"]) == 361 for s in scene["satellites"])


def test_propagator_epoch_offset() -> None:
    # Propagador com época 1 h antes da cena: o instante 0 da cena é t = 3600 s do propagador
    early = EPOCH - timedelta(hours=1)
    scene_a = build_scene([_sso(epoch=early)], [], EPOCH, 600, STEP_S)
    r_expected, _ = KeplerPropagator(SSO_ELEMENTS, early).propagate(3600.0)
    np.testing.assert_allclose(
        scene_a["satellites"][0]["eci_km"][0], r_expected[0] / 1e3, atol=1e-5
    )


# --- Continuidade da longitude -----------------------------------------------------------------


def test_longitude_continuity(scene: dict) -> None:
    for sat in scene["satellites"]:
        lon = np.array(sat["lla"])[:, 1]
        jumps = np.diff(lon)
        wraps = np.abs(jumps) > 180.0  # saltos de ±360° no antimeridiano
        unwrapped = jumps - 360.0 * np.sign(jumps) * wraps
        assert np.max(np.abs(unwrapped)) < 5.0
        assert np.all(np.abs(np.abs(jumps[wraps]) - 360.0) < 5.0)
        # 6 h ≈ 3,8 voltas: o traço cruza o antimeridiano algumas vezes, nunca dezenas
        assert 1 <= np.count_nonzero(wraps) <= 10


# --- Validação de entradas ---------------------------------------------------------------------


def test_rejects_bad_inputs() -> None:
    with pytest.raises(ValueError):
        build_scene([_sso()], [], datetime(2026, 9, 29, 12), 600, STEP_S)
    with pytest.raises(ValueError):
        build_scene([_sso()], [], EPOCH, 605, STEP_S)
    with pytest.raises(ValueError):
        build_scene([_sso()], [], EPOCH, 600, 0)
