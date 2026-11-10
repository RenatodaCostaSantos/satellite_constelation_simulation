"""Testes da cena (S1-08; geodésico na S2-07; ISS real na S3-06; Satellite na S5-01;
constelação na S5-06)."""

import importlib.util
import json
import math
import urllib.error
import urllib.request
from datetime import UTC, datetime, timedelta
from pathlib import Path

import numpy as np
import pytest
from golden.reference_values import (
    CONSTELLATION18_A1_S6,
    CONSTELLATION_A1,
    COVERAGE_A5,
    DEMO_SCENE_A7,
    DEMO_SCENE_S5,
    DEMO_SCENE_S6,
    DESIGN_ORBIT,
    FOOTPRINT_A6,
    PASSES_A2,
    PASSES_A2_META,
    SCENE_S2,
)

from satsim.app.scene_export import (
    COMPACT_PRECISION,
    FULL_PRECISION,
    MAINTENANCE_FIELDS,
    SCHEMA_VERSION,
    SceneConstellation,
    SceneSatellite,
    SceneStation,
    build_scene,
    subpoint_geodetic,
)
from satsim.astrotime import datetime_to_jd
from satsim.constants import R_EARTH
from satsim.constellation import (
    CONSTELLATION_DESIGNS,
    DESIGN_91_6_P18,
    Constellation,
    constellation_palette,
)
from satsim.coverage import roll_for_half_grid, swath_half_angle
from satsim.elements import KeplerianElements, state_to_elements
from satsim.ephemeris import beta_angle, ltan_hours
from satsim.frames import eci_to_ecef, gmst_at
from satsim.geodesy import ecef_to_geodetic, ecef_to_spherical, geodetic_to_ecef
from satsim.groundtrack import node_crossings
from satsim.maintenance import maintenance_status
from satsim.observer import SAO_BENTO_DO_SAPUCAI, footprint_half_angle, observe
from satsim.orbit_design import constellation_phases
from satsim.passes import find_passes
from satsim.propagators import (
    CowellPropagator,
    KeplerPropagator,
    MeanJ2Propagator,
    SGP4Propagator,
)
from satsim.satellite import Satellite
from satsim.tle import Tle

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
    # S2-07: campo novo e compatível; o resto do contrato fica igual
    assert set(scene["meta"]) == {
        "schema_version",
        "epoch_utc",
        "duration_s",
        "step_s",
        "frame",
        "geodesy",
    }
    assert scene["meta"]["geodesy"] == "WGS84"
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
    eci = np.array(scene["satellites"][0]["eci_km"])
    # Sem4 eq. 12 / golden: latitude geocêntrica máxima = 180° − i = 82,56°
    lat_gc = np.degrees(np.arcsin(eci[:, 2] / np.linalg.norm(eci, axis=1)))
    assert np.max(np.abs(lat_gc)) == pytest.approx(D.lat_max_deg, abs=0.01)
    # SCENE_S2 (S2-07): latitude geodésica máxima ≈ 82,61°
    assert np.max(np.abs(lla[:, 0])) == pytest.approx(SCENE_S2.sso_lat_max_geodetic_deg, abs=0.1)
    assert np.all((lla[:, 1] >= -180.0) & (lla[:, 1] < 180.0))


def test_sso_radius_and_altitude(scene: dict) -> None:
    eci = np.array(scene["satellites"][0]["eci_km"])
    lla = np.array(scene["satellites"][0]["lla"])
    r_km = np.linalg.norm(eci, axis=1)
    a_km, e = D.a_km, D.e
    tol = 1e-5  # arredondamento de eci_km (1e-6 km por componente)
    assert np.all(r_km >= a_km * (1 - e) - tol) and np.all(r_km <= a_km * (1 + e) + tol)
    # altitude elipsoidal: entre a esférica (r − R⊕) e ela + (a − b) = 21,4 km (S2-07)
    dalt = lla[:, 2] - (r_km - R_EARTH / 1e3)
    assert np.all(dalt >= -2e-4) and np.all(dalt <= 21.4)


# --- Consistência: lla recalculado a partir de eci_km e do GMST da cena -------------------------


def test_lla_consistent_with_eci_and_gmst(scene: dict) -> None:
    earth = scene["earth"]
    t = np.arange(DURATION_S // STEP_S + 1) * STEP_S
    theta = earth["gmst0_rad"] + earth["gmst_rate_rad_s"] * t  # modelo linear do contrato
    for sat in scene["satellites"]:
        eci_km = np.array(sat["eci_km"])
        lla = np.array(sat["lla"])
        # S2-07: lla = ecef_to_geodetic(eci_to_ecef(eci)) a 1e-6°
        lat, lon, alt = ecef_to_geodetic(eci_to_ecef(eci_km * 1e3, theta))
        np.testing.assert_allclose(np.degrees(lat), lla[:, 0], atol=1e-6)
        dlon = (np.degrees(lon) - lla[:, 1] + 180.0) % 360.0 - 180.0
        assert np.max(np.abs(dlon)) < 1e-6
        np.testing.assert_allclose(alt / 1e3, lla[:, 2], atol=2e-4)


def test_geodetic_vs_spherical_limits(scene: dict) -> None:
    # SCENE_S2 (S2-07): |lat geodésica − geocêntrica| ≤ 0,2° e alt − (r − R⊕) em [−25, +25] km
    earth = scene["earth"]
    t = np.arange(DURATION_S // STEP_S + 1) * STEP_S
    theta = earth["gmst0_rad"] + earth["gmst_rate_rad_s"] * t
    lo, hi = SCENE_S2.alt_vs_spherical_km
    for sat in scene["satellites"]:
        lla = np.array(sat["lla"])
        lat_gc, _, alt_sph = ecef_to_spherical(eci_to_ecef(np.array(sat["eci_km"]) * 1e3, theta))
        dlat = lla[:, 0] - np.degrees(lat_gc)
        assert np.max(np.abs(dlat)) <= SCENE_S2.max_abs_dlat_vs_spherical_deg
        dalt = lla[:, 2] - alt_sph / 1e3
        assert np.all((dalt >= lo) & (dalt <= hi))


def test_eci_to_ecef_rotation() -> None:
    # θ = 90°: o eixo x do ECI aparece em −y no ECEF (a Terra girou 90° para leste). Desde a
    # S2-07 a cena usa satsim.frames.eci_to_ecef (a função local da S1-08 foi removida).
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


# --- S2-07: propagadores com J2 e a ISS-like ---------------------------------------------------

ISS_LIKE = KeplerianElements(
    R_EARTH + 420e3, 0.0, math.radians(SCENE_S2.iss_like_lat_max_geocentric_deg),
    math.radians(120.0), 0.0, 0.0,
)  # fmt: skip


def test_any_propagator_same_function() -> None:
    # S2-07: Kepler, MeanJ2 e Cowell exportam com a mesma função, sem código por tipo
    sats = [
        SceneSatellite("sso", "SSO", "#00e5ff", MeanJ2Propagator(SSO_ELEMENTS, EPOCH)),
        SceneSatellite("iss_j2", "ISS-like J2", "#ffb020", MeanJ2Propagator(ISS_LIKE, EPOCH)),
        SceneSatellite("iss_kep", "ISS-like Kepler", "#9aa4b2", KeplerPropagator(ISS_LIKE, EPOCH)),
        SceneSatellite("iss_cow", "ISS-like Cowell", "#ff4d8d", CowellPropagator(ISS_LIKE, EPOCH)),
    ]
    scene = build_scene(sats, [STATION], EPOCH, 3600, 60)
    assert [s["id"] for s in scene["satellites"]] == ["sso", "iss_j2", "iss_kep", "iss_cow"]
    assert all(len(s["lla"]) == 61 for s in scene["satellites"])
    # em t = 0 os três ISS-like coincidem (MeanJ2 e Kepler partem dos mesmos elementos)
    np.testing.assert_allclose(scene["satellites"][1]["lla"][0], scene["satellites"][2]["lla"][0])


def test_iss_like_max_latitude() -> None:
    # SCENE_S2 (S2-07): ISS-like 51,64° geocêntrica e 51,82° geodésica (o card dizia ≈ 51,7°;
    # corrigido no golden file com a dedução)
    sat = SceneSatellite("iss", "ISS-like", "#ffb020", MeanJ2Propagator(ISS_LIKE, EPOCH))
    scene = build_scene([sat], [], EPOCH, 6 * 3600, 20)
    eci = np.array(scene["satellites"][0]["eci_km"])
    lla = np.array(scene["satellites"][0]["lla"])
    lat_gc = np.degrees(np.arcsin(eci[:, 2] / np.linalg.norm(eci, axis=1)))
    assert np.max(np.abs(lat_gc)) == pytest.approx(
        SCENE_S2.iss_like_lat_max_geocentric_deg, abs=0.01
    )
    assert np.max(np.abs(lla[:, 0])) == pytest.approx(
        SCENE_S2.iss_like_lat_max_geodetic_deg, abs=0.01
    )


def test_iss_like_j2_vs_kepler_after_48h() -> None:
    # SCENE_S2: em 48 h a regressão do nó da ISS-like com J2 vale ≈ −9,9° (Ω̇ = −4,947 °/dia)
    mean_j2 = MeanJ2Propagator(ISS_LIKE, EPOCH)
    d_raan = math.degrees(mean_j2.elements_at(2 * 86400.0).raan[0] - ISS_LIKE.raan)
    assert d_raan == pytest.approx(SCENE_S2.iss_raan_regression_48h_deg, abs=0.05)


# --- S3-06: ISS real, estação, footprint e passagens -------------------------------------------
# Referências: PASSES_A2 e PASSES_A2_META (Apêndice A.2 do prompt da Semana 3), FOOTPRINT_A6
# (A.6). Cena determinística: TLE de teste A.1, início na época do TLE, 3 dias a 20 s.

S3_STEP_S = 20
S3_DURATION_S = 3 * 86400
S3_EL_MIN = 10.0
SBS = SceneStation.from_ground_station("sbs", SAO_BENTO_DO_SAPUCAI)


@pytest.fixture(scope="module")
def iss_tle() -> Tle:
    return Tle.from_text((Path(__file__).parent / "data" / "iss_sample.tle").read_text())


@pytest.fixture(scope="module")
def iss_scene(iss_tle) -> dict:
    iss = SceneSatellite("iss", "ISS", "#ffb020", SGP4Propagator(iss_tle))
    return build_scene([iss], [SBS], iss_tle.epoch, S3_DURATION_S, S3_STEP_S,
                       el_min_deg=S3_EL_MIN, tle=iss_tle, tle_source="arquivo")  # fmt: skip


def _scene_seconds(iso: str, scene: dict) -> float:
    start = datetime.fromisoformat(scene["meta"]["epoch_utc"].replace("Z", "+00:00"))
    return (datetime.fromisoformat(iso.replace("Z", "+00:00")) - start).total_seconds()


def test_s3_scene_passes_match_a2(iss_scene, iss_tle) -> None:
    passes = iss_scene["passes"]
    assert len(passes) == PASSES_A2_META.n_passes_10deg
    # A.2 conta desde a época do TLE; a cena, desde meta.epoch_utc (a mesma, ao milissegundo)
    start = datetime.fromisoformat(iss_scene["meta"]["epoch_utc"].replace("Z", "+00:00"))
    assert iss_scene["meta"]["epoch_utc"] == "2019-12-09T16:38:29.363Z"
    offset = (start - iss_tle.epoch).total_seconds()
    for p, ref in zip(passes, PASSES_A2, strict=True):
        assert p["sat_id"] == "iss" and p["station_id"] == "sbs"
        for key, ref_s in (("aos_utc", ref.aos_s), ("tca_utc", ref.tca_s), ("los_utc", ref.los_s)):
            got = _scene_seconds(p[key], iss_scene) + offset
            assert got == pytest.approx(ref_s, abs=PASSES_A2_META.tol_time_s)
        assert p["max_el_deg"] == pytest.approx(ref.max_el_deg, abs=PASSES_A2_META.tol_el_deg)
        assert p["duration_s"] == pytest.approx(ref.duration_s, abs=1.0)


def test_s3_scene_meta_tle(iss_scene, iss_tle) -> None:
    meta = iss_scene["meta"]
    assert meta["schema_version"] == SCHEMA_VERSION == 1
    assert meta["el_min_deg"] == S3_EL_MIN
    tle = meta["tle"]
    assert set(tle) == {"name", "norad_id", "epoch_utc", "age_days_at_start", "source"}
    assert tle["norad_id"] == 25544 and tle["source"] == "arquivo"
    assert tle["epoch_utc"] == "2019-12-09T16:38:29.363Z"
    assert tle["age_days_at_start"] == pytest.approx(0.0, abs=1e-4)


def test_s3_looks_lengths_ranges_and_crossings(iss_scene) -> None:
    n = S3_DURATION_S // S3_STEP_S + 1
    looks = iss_scene["stations"][0]["looks"]["iss"]
    assert set(looks) == {"az_deg", "el_deg", "range_km"}
    assert all(len(v) == n for v in looks.values())
    az = np.array(looks["az_deg"])
    el = np.array(looks["el_deg"])
    assert np.all((az >= 0.0) & (az < 360.0)) and np.all((el >= -90.0) & (el <= 90.0))
    assert np.all(np.array(looks["range_km"]) > 300.0)
    t = np.arange(n) * S3_STEP_S
    # el cruza el_min a menos de 1 passo (±20 s) de cada AOS (subindo) e LOS (descendo)
    above = el >= S3_EL_MIN
    rises = t[1:][~above[:-1] & above[1:]]
    sets = t[1:][above[:-1] & ~above[1:]]
    aos = np.array([_scene_seconds(p["aos_utc"], iss_scene) for p in iss_scene["passes"]])
    los = np.array([_scene_seconds(p["los_utc"], iss_scene) for p in iss_scene["passes"]])
    assert len(rises) == len(aos) and len(sets) == len(los)
    assert np.max(np.abs(rises - aos)) <= S3_STEP_S
    assert np.max(np.abs(sets - los)) <= S3_STEP_S


def test_s3_looks_consistent_with_observe(iss_scene, iss_tle) -> None:
    looks = iss_scene["stations"][0]["looks"]["iss"]
    k = np.array([0, 1000, 5000, 12000])
    t = k * S3_STEP_S
    start = datetime.fromisoformat(iss_scene["meta"]["epoch_utc"].replace("Z", "+00:00"))
    offset = (start - iss_tle.epoch).total_seconds()  # < 1 ms (época arredondada ao ms)
    obs = observe(SGP4Propagator(iss_tle), SAO_BENTO_DO_SAPUCAI, t + offset)
    np.testing.assert_allclose(np.array(looks["el_deg"])[k], obs.el_deg, atol=1e-3)
    np.testing.assert_allclose(np.array(looks["range_km"])[k], obs.range / 1e3, atol=2e-3)


def test_s3_footprint_matches_formula(iss_scene) -> None:
    sat = iss_scene["satellites"][0]
    fp = np.array(sat["footprint_deg"])
    alt_km = np.array(sat["lla"])[:, 2]  # a altitude publicada na própria cena
    expected = np.degrees(footprint_half_angle(alt_km * 1e3, math.radians(S3_EL_MIN)))
    np.testing.assert_allclose(fp, expected, rtol=0, atol=1e-6)
    # valor típico: ≈ 12,5° para a ISS a ~420 km com el_min = 10° (A.6)
    assert np.median(fp) == pytest.approx(FOOTPRINT_A6.h420_el10_deg, abs=0.3)


def test_s3_new_fields_are_optional(scene: dict) -> None:
    # sem el_min_deg, a cena é a mesma das semanas anteriores (o viewer funciona sem os campos)
    assert "el_min_deg" not in scene["meta"] and "tle" not in scene["meta"]
    assert all("footprint_deg" not in s for s in scene["satellites"])
    assert all("looks" not in st for st in scene["stations"])
    assert scene["passes"] == []
    # S4-06: sem sun_geometry, nada de LTAN, β e eclipse
    for field in ("ltan_h", "beta_deg", "eclipse"):
        assert all(field not in s for s in scene["satellites"])


def test_s3_multi_satellite_iss_sgp4_and_sso_j2(iss_tle) -> None:
    iss = SceneSatellite("iss", "ISS", "#ffb020", SGP4Propagator(iss_tle))
    sso = SceneSatellite("sso", "SSO", "#00e5ff", MeanJ2Propagator(SSO_ELEMENTS, EPOCH))
    start = iss_tle.epoch.replace(microsecond=0)
    sc = build_scene([iss, sso], [SBS], start, 6 * 3600, 20, el_min_deg=10.0, tle=iss_tle)
    assert [s["id"] for s in sc["satellites"]] == ["iss", "sso"]
    assert set(sc["stations"][0]["looks"]) == {"iss", "sso"}
    assert all(len(s["footprint_deg"]) == 1081 for s in sc["satellites"])
    # cada satélite com a sua época: em t = 0 da cena, a SSO está onde o MeanJ2 a põe
    r_sso, _ = MeanJ2Propagator(SSO_ELEMENTS, EPOCH).propagate_utc(start, [0.0])
    np.testing.assert_allclose(sc["satellites"][1]["eci_km"][0], r_sso[0] / 1e3, atol=1e-5)
    assert [p["aos_utc"] for p in sc["passes"]] == sorted(p["aos_utc"] for p in sc["passes"])


def test_s3_scene_passes_coherent_with_find_passes(iss_scene, iss_tle) -> None:
    # o viewer mostra as mesmas passagens que scripts/list_passes.py (mesma função)
    start = datetime.fromisoformat(iss_scene["meta"]["epoch_utc"].replace("Z", "+00:00"))
    direct = find_passes(SGP4Propagator(iss_tle), SAO_BENTO_DO_SAPUCAI, start, S3_DURATION_S)
    assert len(direct) == len(iss_scene["passes"])
    for p, q in zip(iss_scene["passes"], direct, strict=True):
        # os instantes do JSON são arredondados ao milissegundo
        dt = _scene_seconds(p["aos_utc"], iss_scene) - (q.aos - start).total_seconds()
        assert abs(dt) <= 0.5e-3
        assert p["max_el_deg"] == pytest.approx(q.max_elevation_deg, abs=1e-3)


def test_s3_scene_station_from_ground_station() -> None:
    assert SBS.ground_station() == SAO_BENTO_DO_SAPUCAI
    assert (SBS.lat_deg, SBS.lon_deg, SBS.alt_m) == (-22.69, -45.73, 900.0)


# --- S3-06: build_demo.py offline ---------------------------------------------------------------


def _load_build_demo():
    path = Path(__file__).resolve().parents[1] / "scripts" / "build_demo.py"
    spec = importlib.util.spec_from_file_location("build_demo", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _scene_from_html(html: str) -> dict:
    start = html.index("window.SATSIM_SCENE = ") + len("window.SATSIM_SCENE = ")
    end = html.index(";\nwindow.SATSIM_TEXTURE", start)
    return json.loads(html[start:end])


def test_build_demo_offline_without_network(tmp_path, monkeypatch, reference_constellation):
    def no_network(*args, **kwargs):
        raise urllib.error.URLError("sem rede (teste)")

    monkeypatch.setattr(urllib.request, "urlopen", no_network)
    out = tmp_path / "index.html"
    demo = _load_build_demo()
    # A constelação de referência do conftest (mesma época, Cowell completo já integrado na
    # sessão) no lugar da do script, que difere só pelo a calibrado (+19 m): poupa ≈ 7 s
    monkeypatch.setattr(
        demo, "build_constellation", lambda start, n_sats=6: reference_constellation
    )
    # Início fixo na data do A.7: os dois TLEs embarcados (ISS e Sentinel-2A) servem
    demo.main(["--offline", "--cache-dir", str(tmp_path / "cache"), "-o", str(out),
               "--start", DEMO_SCENE_A7.start_utc])  # fmt: skip
    html = out.read_text(encoding="utf-8")
    scene = _scene_from_html(html)
    assert scene["meta"]["tle"]["source"] == "bundled"
    ids = [s["id"] for s in scene["satellites"]]
    assert ids == ["iss", *DEMO_SCENE_S5.member_ids, "s2a"]
    assert scene["meta"]["duration_s"] == 24 * 3600 and scene["meta"]["step_s"] == 20
    assert "looks" in scene["stations"][0] and "footprint_deg" in scene["satellites"][0]
    assert 'src="vendor' not in html and "<script src=" not in html  # tudo embutido
    # S4-06: LTAN, β e eclipse em todos os satélites, N amostras cada
    n = 24 * 3600 // 20 + 1
    for sat in scene["satellites"]:
        assert len(sat["ltan_h"]) == len(sat["beta_deg"]) == len(sat["eclipse"]) == n
        assert all(isinstance(x, bool) for x in sat["eclipse"])
    # S1 é a SSO de projeto do A.7 (S4-06): LTAN constante em 22h30 e eclipse de ≈ 36%
    sso = scene["satellites"][1]
    ltan = np.array(sso["ltan_h"])
    assert np.all(np.abs(ltan - DEMO_SCENE_A7.ltan_h) <= DEMO_SCENE_A7.ltan_tol_h)
    lo, hi = DEMO_SCENE_A7.ltan_range_h
    assert lo - 2e-4 <= ltan.min() and ltan.max() <= hi + 2e-4
    fraction = float(np.mean(sso["eclipse"]))
    assert fraction == pytest.approx(DEMO_SCENE_A7.eclipse_fraction,
                                     abs=DEMO_SCENE_A7.eclipse_fraction_tol)  # fmt: skip
    assert fraction == pytest.approx(DEMO_SCENE_A7.eclipse_fraction_a4, abs=1e-3)
    beta = np.array(sso["beta_deg"])
    assert DEMO_SCENE_A7.beta_range_deg[0] - 0.01 <= beta.min()
    assert beta.max() <= DEMO_SCENE_A7.beta_range_deg[1] + 0.01
    _check_demo_constellation(scene)


def _check_demo_constellation(scene: dict) -> None:
    """Critérios de cena do S5-06 sobre a demo de 1 dia com a constelação no Cowell completo."""
    ref = DEMO_SCENE_S5
    meta = scene["meta"]["constellation"]
    assert meta["member_ids"] == list(ref.member_ids) and meta["n_sats"] == 6
    assert (meta["n_revs"], meta["n_days"]) == (91, 6)
    assert meta["grid_deg"] == pytest.approx(360.0 / 91, abs=1e-6)
    assert meta["roll_deg"] == pytest.approx(COVERAGE_A5.roll_half_grid_deg, abs=1e-4)
    members = [s for s in scene["satellites"] if s["id"] in ref.member_ids]
    others = [s for s in scene["satellites"] if s["id"] not in ref.member_ids]
    # 91 nós somados em 1 dia (A.3), em ordem de tempo e dentro da janela
    assert sum(len(s["nodes"]["t_s"]) for s in members) == ref.nodes_1d
    for s in members:
        t_nodes = np.array(s["nodes"]["t_s"])
        assert len(t_nodes) == len(s["nodes"]["lon_deg"]) and np.all(np.diff(t_nodes) > 0)
        assert 0.0 <= t_nodes[0] and t_nodes[-1] <= 24 * 3600
    # LTAN dos seis constante em 24 h (22,500 ± 0,01 h) e iguais entre si
    ltan = np.array([s["ltan_h"] for s in members])
    assert np.all(np.abs(ltan - ref.ltan_h) <= ref.ltan_tol_h)
    assert np.max(np.ptp(ltan, axis=0)) <= ref.ltan_spread_tol_h
    # faixa de meia grade e manutenção do S5-05 em cada membro; os demais sem os campos
    for s in members:
        assert s["swath_deg"] == pytest.approx(180.0 / 91, rel=0.02)
        m = s["maintenance"]
        assert set(m) == set(MAINTENANCE_FIELDS)
        assert m["days_to_track_correction"] == pytest.approx(6.95, rel=0.02)  # A.7
    for s in others:
        assert not {"nodes", "swath_deg", "maintenance"} & set(s)


def test_s4_sun_geometry_fields_match_ephemeris() -> None:
    sc = build_scene([_sso()], [STATION], EPOCH, 3 * 3600, 60, sun_geometry=True)
    sat = sc["satellites"][0]
    n = 3 * 3600 // 60 + 1
    assert len(sat["ltan_h"]) == len(sat["beta_deg"]) == len(sat["eclipse"]) == n
    t = np.arange(n) * 60.0
    r, v = KeplerPropagator(SSO_ELEMENTS, EPOCH).propagate(t)
    jd = datetime_to_jd(EPOCH) + t / 86400.0
    ltan = ltan_hours(state_to_elements(r, v).raan, jd)
    np.testing.assert_allclose(sat["ltan_h"], ltan, atol=1e-4)
    np.testing.assert_allclose(sat["beta_deg"], np.degrees(beta_angle(r, v, jd)), atol=1e-4)
    # numa volta de ~95 min o satélite entra e sai da sombra (β ≈ 0 em Ω = 30° nesta data)
    assert any(sat["eclipse"]) and not all(sat["eclipse"])
    assert json.loads(json.dumps(sc)) == sc


def _fast_constellation(start: datetime, n_sats: int = 6) -> Constellation:
    """Constelação 91/6 barata (J2 médio): para testes que não verificam a dinâmica (TLE,
    tamanho do arquivo). O tamanho da cena depende do número de amostras, não do propagador.
    Com ``n_sats`` (S6-02), a de 18 com a paleta por família."""
    palette = constellation_palette(CONSTELLATION_DESIGNS[n_sats])
    sats = [
        Satellite(
            f"s{j + 1}",
            f"S{j + 1}",
            MeanJ2Propagator(
                KeplerianElements(
                    SSO_ELEMENTS.a,
                    SSO_ELEMENTS.e,
                    SSO_ELEMENTS.i,
                    SSO_ELEMENTS.raan,
                    SSO_ELEMENTS.argp,
                    float(psi),
                ),
                start,
            ),
            color=palette[j],
        )
        for j, psi in enumerate(constellation_phases(n_sats))
    ]
    return Constellation(start, sats)


def test_build_demo_network_failure_falls_back(tmp_path, monkeypatch) -> None:
    # sem --offline, a rede falha e get_tle cai no TLE embarcado: a demo sai igual
    def no_network(*args, **kwargs):
        raise urllib.error.URLError("sem rede (teste)")

    monkeypatch.setattr(urllib.request, "urlopen", no_network)
    out = tmp_path / "index.html"
    demo = _load_build_demo()
    monkeypatch.setattr(demo, "build_constellation", _fast_constellation)
    demo.main(["--cache-dir", str(tmp_path / "cache"), "-o", str(out)])
    assert _scene_from_html(out.read_text(encoding="utf-8"))["meta"]["tle"]["source"] == "bundled"


def test_build_demo_deterministic_tle_epoch(tmp_path, monkeypatch) -> None:
    out = tmp_path / "index.html"
    tle_file = Path(__file__).parent / "data" / "iss_sample.tle"
    demo = _load_build_demo()
    monkeypatch.setattr(demo, "build_constellation", _fast_constellation)
    demo.main(["--tle-file", str(tle_file), "--start", "tle-epoch", "-o", str(out)])
    scene = _scene_from_html(out.read_text(encoding="utf-8"))
    # TLE da ISS de 2019: o do Sentinel-2A (2026) fica a anos da cena e é omitido
    assert [s["id"] for s in scene["satellites"]] == ["iss", *DEMO_SCENE_S5.member_ids]
    iss_passes = [p for p in scene["passes"] if p["sat_id"] == "iss"]
    # 24 h após a época: as passagens 1 e 2 do A.2 (AOS em 29 258,7 s e 64 638,7 s)
    assert len(iss_passes) == 2
    assert iss_passes[0]["max_el_deg"] == pytest.approx(PASSES_A2[0].max_el_deg, abs=0.05)


# --- S5-01: Satellite e Constellation na cena ---------------------------------------------------


def _as_satellite(sat: SceneSatellite) -> Satellite:
    return Satellite(sat.id, sat.name, sat.propagator, color=sat.color)


def test_s5_scene_accepts_satellite_same_contract(scene: dict) -> None:
    # A mesma cena do fixture, montada com a classe Satellite do ADR: JSON idêntico
    sats = [_as_satellite(s) for s in (_sso(), _example())]
    got = build_scene(sats, [STATION], EPOCH, DURATION_S, STEP_S)
    assert got == scene


def test_s5_scene_accepts_constellation_and_mixed_list() -> None:
    # Constellation (época própria, diferente da cena) e lista mista Satellite + SceneSatellite
    sso, example = _sso(), _example()
    constellation = Constellation(EPOCH - timedelta(hours=2), [_as_satellite(sso)])
    from_constellation = build_scene(constellation, [STATION], EPOCH, 3600, 60)
    from_scene_sat = build_scene([sso], [STATION], EPOCH, 3600, 60)
    assert from_constellation == from_scene_sat
    mixed = build_scene([_as_satellite(sso), example], [STATION], EPOCH, 3600, 60)
    assert [s["id"] for s in mixed["satellites"]] == ["sso1", "kep1"]
    assert [s["color"] for s in mixed["satellites"]] == ["#00e5ff", "#ffb020"]
    assert mixed["meta"]["schema_version"] == SCHEMA_VERSION


# --- S5-06: constelação na cena -------------------------------------------------------------------


def _fast_info(constellation: Constellation, **kwargs) -> SceneConstellation:
    return SceneConstellation("Constelação 91/6", 91, 6, tuple(constellation.ids), **kwargs)


def test_s5_constellation_fields_and_meta() -> None:
    c = _fast_constellation(EPOCH)
    roll = math.radians(23.1957)
    status = maintenance_status("s1", -23.71 / 86400, SSO_ELEMENTS.a, SSO_ELEMENTS.i)
    info = _fast_info(c, roll_rad=roll, maintenance={"s1": status})
    sc = build_scene([_example(), *c], [STATION], EPOCH, 6 * 3600, 60, constellation_info=info)
    assert sc["meta"]["schema_version"] == SCHEMA_VERSION  # campos opcionais, mesma versão
    assert sc["meta"]["constellation"] == {
        "name": "Constelação 91/6", "n_revs": 91, "n_days": 6, "n_sats": 6,
        "grid_deg": round(360 / 91, 6), "member_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
        "roll_deg": 23.1957,
    }  # fmt: skip
    kep, s1, s2 = sc["satellites"][:3]
    assert not {"nodes", "swath_deg", "maintenance"} & set(kep)  # fora da constelação
    # nós: os mesmos de node_crossings, em segundos desde o início e graus
    t_ref, lon_ref = node_crossings(c[0].propagator, EPOCH, 6 * 3600)
    np.testing.assert_allclose(s1["nodes"]["t_s"], t_ref, atol=1e-3)
    np.testing.assert_allclose(s1["nodes"]["lon_deg"], np.degrees(lon_ref), atol=1e-6)
    assert len(s1["nodes"]["t_s"]) in (3, 4)  # ≈ 3,8 voltas em 6 h
    # faixa: eq. 45 no raio médio da órbita nas amostras da cena
    r, _ = c[0].propagator.propagate_utc(EPOCH, np.arange(0.0, 6 * 3600 + 1, 60.0))
    lam = math.degrees(swath_half_angle(float(np.mean(np.linalg.norm(r, axis=1))), roll))
    assert s1["swath_deg"] == pytest.approx(lam, abs=1e-6)
    # manutenção só de quem tem indicador, com os campos do resumo
    assert s1["maintenance"]["days_to_track_correction"] == pytest.approx(6.95, rel=0.02)
    assert set(s1["maintenance"]) == set(MAINTENANCE_FIELDS) and "maintenance" not in s2
    assert json.loads(json.dumps(sc, allow_nan=False)) == sc


def test_s5_constellation_infinite_maintenance_and_no_roll() -> None:
    # sem decaimento a próxima correção é infinita: vira null (JSON válido)
    c = _fast_constellation(EPOCH)
    status = maintenance_status("s1", 0.0, SSO_ELEMENTS.a, SSO_ELEMENTS.i)
    sc = build_scene(c, [STATION], EPOCH, 3600, 60,
                     constellation_info=_fast_info(c, maintenance={"s1": status}))  # fmt: skip
    m = sc["satellites"][0]["maintenance"]
    assert m["days_to_track_correction"] is None and m["track_delta_v_per_year_m_s"] == 0.0
    assert sc["meta"]["constellation"]["roll_deg"] is None
    assert "swath_deg" not in sc["satellites"][0]
    json.dumps(sc, allow_nan=False)


def test_s5_constellation_rejects_missing_member() -> None:
    c = _fast_constellation(EPOCH)
    info = SceneConstellation("x", 91, 6, ("s1", "s9"))
    with pytest.raises(ValueError, match="s9"):
        build_scene(c, [STATION], EPOCH, 3600, 60, constellation_info=info)


def test_s5_old_scenes_have_no_constellation_fields(scene: dict) -> None:
    # compatibilidade: sem constellation_info, nada novo na cena (o viewer funciona sem eles)
    assert "constellation" not in scene["meta"]
    for sat in scene["satellites"]:
        assert not {"nodes", "swath_deg", "maintenance"} & set(sat)


def test_build_demo_six_days_under_size_limit(tmp_path, monkeypatch) -> None:
    """--days 6, offline: o ciclo inteiro abaixo de 15 MB (a constelação barata tem o mesmo
    número de amostras que a do Cowell; a demo real de 6 dias mede 11,2 MB)."""

    def no_network(*args, **kwargs):
        raise urllib.error.URLError("sem rede (teste)")

    monkeypatch.setattr(urllib.request, "urlopen", no_network)
    out = tmp_path / "index.html"
    demo = _load_build_demo()
    monkeypatch.setattr(demo, "build_constellation", _fast_constellation)
    demo.main(["--offline", "--cache-dir", str(tmp_path / "cache"), "-o", str(out),
               "--start", DEMO_SCENE_S5.start_utc, "--days", "6"])  # fmt: skip
    size_mb = out.stat().st_size / 1e6
    assert size_mb < DEMO_SCENE_S5.max_file_mb
    scene = _scene_from_html(out.read_text(encoding="utf-8"))
    assert scene["meta"]["duration_s"] == 6 * 86400
    assert scene["meta"]["step_s"] == DEMO_SCENE_S5.step_6d_s
    ids = [s["id"] for s in scene["satellites"]]
    assert ids == ["iss", *DEMO_SCENE_S5.member_ids, "s2a"]
    # 6 dias: cada satélite dá as 91 voltas do ciclo (91 nós, ± 1 pela borda da janela)
    for s in scene["satellites"][1:7]:
        assert abs(len(s["nodes"]["t_s"]) - 91) <= 1
    with pytest.raises(SystemExit):
        demo.main(["--offline", "--days", "7"])


def test_build_demo_constellation_recipe() -> None:
    # a receita do script, sem propagar: 6 membros no Cowell completo, época = início da cena
    demo = _load_build_demo()
    start = datetime.fromisoformat(DEMO_SCENE_S5.start_utc)
    c = demo.build_constellation(start)
    assert c.ids == list(DEMO_SCENE_S5.member_ids) and c.epoch == start
    assert all(len(s.propagator.models) == 5 for s in c)
    assert demo.A_CALIBRATED_M == pytest.approx(6888108.374, abs=0.01)  # S5-02
    assert set(demo.STEP_BY_DAYS) == set(range(1, 7))
    for days, step in demo.STEP_BY_DAYS.items():
        assert (days * 86400) % step == 0


# --- S6-02: cenas de 6 e 18 satélites, variantes e interpolação do viewer -------------------------

S6 = DEMO_SCENE_S6
A1_EPOCH = datetime.fromisoformat(CONSTELLATION_A1.epoch_utc)
A_DESIGN = CONSTELLATION_A1.a_design_m
NEW_META_FIELDS = {"tracks_per_day", "daily_grid_deg", "families"}


def _info18(constellation: Constellation) -> SceneConstellation:
    roll = DESIGN_91_6_P18.roll_for_daily_coverage(A_DESIGN)
    return SceneConstellation.from_design(DESIGN_91_6_P18, "Constelação 91/6",
                                          constellation.ids, roll_rad=roll)  # fmt: skip


def test_s6_scene_18_members_families_and_nodes(constellation18) -> None:
    """Cena de 18 (Cowell da fixture): 18 membros, 273 trilhas/dia, famílias 0/1/2 com 6 cada,
    273 nós somados em 1 dia e faixa coerente com a rolagem de 8,1835°."""
    sc = build_scene(constellation18, [STATION], A1_EPOCH, 86400, 60,
                     constellation_info=_info18(constellation18))  # fmt: skip
    meta = sc["meta"]["constellation"]
    assert meta["n_sats"] == 18 and meta["member_ids"] == list(DESIGN_91_6_P18.ids)
    assert meta["tracks_per_day"] == 273 and meta["families"] == 3
    assert meta["daily_grid_deg"] == pytest.approx(360 / 273, abs=1e-6)
    assert meta["roll_deg"] == pytest.approx(S6.roll18_deg, abs=1e-4)
    families = [s["family"] for s in sc["satellites"]]
    assert families == [0, 1, 2] * 6 and all(families.count(f) == 6 for f in range(3))
    assert sum(len(s["nodes"]["t_s"]) for s in sc["satellites"]) == CONSTELLATION18_A1_S6.p18_nodes
    roll = DESIGN_91_6_P18.roll_for_daily_coverage(A_DESIGN)  # a da cena (8,18348°)
    for sat, s in zip(constellation18, sc["satellites"], strict=True):
        r, _ = sat.propagate_utc(A1_EPOCH, np.arange(0.0, 86400 + 1, 60.0))
        lam = math.degrees(swath_half_angle(float(np.mean(np.linalg.norm(r, axis=1))), roll))
        assert s["swath_deg"] == pytest.approx(lam, abs=1e-6)
        assert s["swath_deg"] == pytest.approx(180 / 273, rel=0.02)  # meia grade diária
    assert json.loads(json.dumps(sc, allow_nan=False)) == sc


def _strip_s6_fields(scene: dict) -> dict:
    out = json.loads(json.dumps(scene))
    for key in NEW_META_FIELDS:
        out["meta"]["constellation"].pop(key, None)
    for sat in out["satellites"]:
        sat.pop("family", None)
    return out


def test_s6_scene_6_identical_except_new_optional_fields() -> None:
    """A cena de 6 da demo é a da S5-06 (mesma rolagem, passo e precisão) mais os campos
    opcionais novos."""
    demo = _load_build_demo()
    c = _fast_constellation(EPOCH)
    info_new = demo.constellation_info(c)
    info_s5 = SceneConstellation("Constelação 91/6", 91, 6, tuple(c.ids),
                                 roll_rad=roll_for_half_grid(A_DESIGN, 91),
                                 maintenance=info_new.maintenance)  # fmt: skip
    sats = [_example(), *c]
    kwargs = {"el_min_deg": 10.0, "sun_geometry": True}
    step = demo.STEP_BY_SATS_DAYS[(6, 1)]
    new = build_scene(sats, [STATION], EPOCH, 3 * 3600, step, constellation_info=info_new,
                      precision=demo.PRECISION_BY_SATS[6], **kwargs)  # fmt: skip
    s5 = build_scene(sats, [STATION], EPOCH, 3 * 3600, step, constellation_info=info_s5, **kwargs)
    assert NEW_META_FIELDS <= set(new["meta"]["constellation"])
    assert new["meta"]["constellation"]["tracks_per_day"] == 91
    assert new["meta"]["constellation"]["families"] == 1
    assert _strip_s6_fields(new) == s5
    assert demo.PRECISION_BY_SATS[6] == FULL_PRECISION
    assert {d: demo.STEP_BY_SATS_DAYS[(6, d)] for d in range(1, 7)} == demo.STEP_BY_DAYS


def test_s6_compact_precision_and_variants_contract() -> None:
    c = _fast_constellation(EPOCH, n_sats=18)
    variants = [{"id": "p6", "label": "6 satélites", "href": "index.html"},
                {"id": "p18", "label": "18 satélites", "href": "constelacao18.html"}]  # fmt: skip
    sc = build_scene(c, [STATION], EPOCH, 3600, 60, el_min_deg=10.0, constellation_info=_info18(c),
                     variants=variants, precision=COMPACT_PRECISION)  # fmt: skip
    assert sc["meta"]["variants"] == variants and sc["meta"]["schema_version"] == SCHEMA_VERSION
    full = build_scene(c, [STATION], EPOCH, 3600, 60, el_min_deg=10.0)
    # mesma posição a menos do arredondamento: 1 m em ECI, 1e-5° em lat/lon
    for a, b in zip(sc["satellites"], full["satellites"], strict=True):
        np.testing.assert_allclose(a["eci_km"], b["eci_km"], atol=5e-4 + 1e-12)
        np.testing.assert_allclose(a["lla"], b["lla"], atol=5e-4)
    assert len(json.dumps(sc)) < 0.85 * len(json.dumps(full))
    for bad in ({"id": "x", "label": "x", "href": "/demo/x.html"},
                {"id": "x", "label": "x", "href": "https://example.org/x.html"},
                {"id": "x", "href": "x.html"}):  # fmt: skip
        with pytest.raises(ValueError):
            build_scene(c, [STATION], EPOCH, 60, 60, variants=[bad])
    with pytest.raises(ValueError, match="18 satélites"):
        SceneConstellation.from_design(DESIGN_91_6_P18, "x", ("s1", "s2"))
    # cenas antigas: sem variants nem families (o viewer abre como antes)
    old = build_scene(c, [STATION], EPOCH, 60, 60,
                      constellation_info=SceneConstellation("x", 91, 6, tuple(c.ids)))  # fmt: skip
    assert "variants" not in old["meta"]
    assert not NEW_META_FIELDS & set(old["meta"]["constellation"])
    assert all("family" not in s for s in old["satellites"])


def test_s6_build_demo_all_variants(tmp_path, monkeypatch) -> None:
    """--all --offline: os quatro arquivos, com meta.variants cruzadas e relativas, < 15 MB cada
    (constelação barata: o tamanho depende do número de amostras, não do propagador)."""

    def no_network(*args, **kwargs):
        raise urllib.error.URLError("sem rede (teste)")

    monkeypatch.setattr(urllib.request, "urlopen", no_network)
    demo = _load_build_demo()
    monkeypatch.setattr(demo, "build_constellation", _fast_constellation)
    demo.main(["--all", "--offline", "--cache-dir", str(tmp_path / "cache"),
               "--out-dir", str(tmp_path), "--start", DEMO_SCENE_S5.start_utc])  # fmt: skip
    names = {(6, 1): "index.html", (18, 1): "constelacao18.html",
             (6, 6): "ciclo_6dias.html", (18, 6): "ciclo_6dias_18.html"}  # fmt: skip
    assert sorted(p.name for p in tmp_path.glob("*.html")) == sorted(names.values())
    for (n_sats, days), name in names.items():
        path = tmp_path / name
        assert path.stat().st_size / 1e6 < S6.max_file_mb, name
        sc = _scene_from_html(path.read_text(encoding="utf-8"))
        meta = sc["meta"]
        assert meta["constellation"]["n_sats"] == n_sats and meta["duration_s"] == days * 86400
        assert meta["step_s"] == demo.STEP_BY_SATS_DAYS[(n_sats, days)]
        assert meta["epoch_utc"] == "2026-10-23T00:00:00Z"  # mesma época: o instante se preserva
        variants = meta["variants"]
        assert [v["id"] for v in variants] == ["p6", "p18"]
        assert [v["label"] for v in variants] == ["6 satélites", "18 satélites"]
        for v in variants:
            href = v["href"]
            assert "/" not in href and ":" not in href  # relativo, na mesma pasta
            assert (tmp_path / href).exists()
            sibling = _scene_from_html((tmp_path / href).read_text(encoding="utf-8"))
            assert sibling["meta"]["duration_s"] == meta["duration_s"]  # irmão de mesma duração
        assert meta["constellation"]["roll_deg"] == pytest.approx(
            S6.roll18_deg if n_sats == 18 else COVERAGE_A5.roll_half_grid_deg, abs=1e-4
        )
    with pytest.raises(SystemExit):
        demo.main(["--offline", "--sats", "7"])
    with pytest.raises(SystemExit):
        demo.main(["--all", "--offline", "-o", str(tmp_path / "x.html")])


def test_s6_default_output_names() -> None:
    demo = _load_build_demo()
    assert demo.default_output(6, 1).name == "index.html"
    assert demo.default_output(18, 1).name == "constelacao18.html"
    assert demo.default_output(6, 6).name == "ciclo_6dias.html"
    assert demo.default_output(18, 6).name == "ciclo_6dias_18.html"
    assert demo.default_output(18, 3).name == "ciclo_3dias_18.html"
    assert set(demo.STEP_BY_SATS_DAYS) == {(p, d) for p in (6, 18) for d in range(1, 7)}
    for (_, days), step in demo.STEP_BY_SATS_DAYS.items():
        assert (days * 86400) % step == 0


# --- Interpolação do viewer (reimplementada em Python) ---------------------------------------


def _viewer_sample(lla: np.ndarray, step: float, t: np.ndarray) -> np.ndarray:
    """O ``sampleAt`` do viewer (S6-02): Hermite cúbica (Catmull-Rom) no vetor unitário de
    (lat, lon) e na altitude, tangentes por diferença central (unilateral nas pontas)."""
    n = len(lla)
    f = np.clip(t / step, 0.0, n - 1 - 1e-12)
    k = np.floor(f).astype(int)
    u = (f - k)[:, None]
    i0, i3 = np.maximum(k - 1, 0), np.minimum(k + 2, n - 1)
    d1, d2 = (k + 1 - i0)[:, None], (i3 - k)[:, None]
    lat, lon = np.radians(lla[:, 0]), np.radians(lla[:, 1])
    unit = np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    vals = np.column_stack([unit, lla[:, 2]])
    p0, p1, p2, p3 = vals[i0], vals[k], vals[k + 1], vals[i3]
    h00, h10 = 2 * u**3 - 3 * u**2 + 1, u**3 - 2 * u**2 + u
    h01, h11 = -2 * u**3 + 3 * u**2, u**3 - u**2
    v = h00 * p1 + h10 * (p2 - p0) / d1 + h01 * p2 + h11 * (p3 - p1) / d2
    xyz, alt = v[:, :3], v[:, 3]
    r = np.linalg.norm(xyz, axis=1)
    return np.column_stack([np.degrees(np.arcsin(xyz[:, 2] / r)),
                            np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0])), alt])  # fmt: skip


def _linear_lla(lla: np.ndarray, step: float, t: np.ndarray) -> np.ndarray:
    """O ``sampleAt`` linear em lat/lon das Semanas 3–5 (para comparação)."""
    k = np.minimum(np.floor(t / step).astype(int), len(lla) - 2)
    u = (t / step - k)[:, None]
    a, b = lla[k], lla[k + 1]
    d = b - a
    d[:, 1] = np.mod(d[:, 1] + 180.0, 360.0) - 180.0
    return a + u * d


def _position_error_km(sat, step: float, method, decimals) -> float:
    """Maior erro de posição [km] do método de interpolação entre amostras da cena (1 dia)."""
    t_samples = np.arange(0.0, 86400.0 + 1.0, step)
    r, _ = sat.propagate_utc(A1_EPOCH, t_samples)
    lla = subpoint_geodetic(eci_to_ecef(r, gmst_at(A1_EPOCH, t_samples)))
    lla = np.column_stack([np.round(lla[:, :2], decimals.angle), np.round(lla[:, 2], decimals.alt)])
    t = np.concatenate([t_samples[:-1] + q * step for q in (0.1, 0.25, 0.5, 0.75, 0.9)])
    r_true, _ = sat.propagate_utc(A1_EPOCH, t)
    truth = eci_to_ecef(r_true, gmst_at(A1_EPOCH, t))
    est = method(lla, step, t)
    pos = geodetic_to_ecef(np.radians(est[:, 0]), np.radians(est[:, 1]), est[:, 2] * 1e3)
    return float(np.max(np.linalg.norm(pos - truth, axis=1)) / 1e3)


def test_s6_viewer_interpolation_error_per_step(reference_constellation) -> None:
    """Com o passo de cada (P, dias), o erro de posição do método do viewer fica ≤ 5 km."""
    demo = _load_build_demo()
    sat = reference_constellation[0]  # Cowell completo, 1 dia já integrado na sessão
    errors = {}
    for (n_sats, _), step in sorted(demo.STEP_BY_SATS_DAYS.items(), key=lambda kv: kv[1]):
        if step in errors:
            continue
        errors[step] = _position_error_km(sat, step, _viewer_sample,
                                          demo.PRECISION_BY_SATS[n_sats])  # fmt: skip
    assert max(errors.values()) <= S6.interp_limit_km, errors
    assert errors[120] == pytest.approx(S6.interp_error_120s_km, abs=0.05)
    # o método antigo (linear em lat/lon) passa de 5 km já com 60 s, perto da latitude máxima
    assert _position_error_km(sat, 60.0, _linear_lla, FULL_PRECISION) > S6.interp_limit_km


def test_s6_viewer_implements_the_tested_interpolation() -> None:
    """O viewer usa o método reimplementado acima (vetor unitário + Hermite) e o seletor."""
    static = Path(__file__).resolve().parents[1] / "src" / "satsim" / "app" / "static"
    js = (static / "viewer.js").read_text(encoding="utf-8")
    for snippet in ("function unitsOf(lla)", "var h00 = 2 * u3 - 3 * u2 + 1;",
                    "((p2 - p0) / d1)", "((p3 - p1) / d2)", "Math.asin(v[2] / r)",
                    "meta.variants", "function readHash()", "tracks_per_day"):  # fmt: skip
        assert snippet in js, snippet
    html = (static / "viewer.html").read_text(encoding="utf-8")
    assert 'id="variants"' in html and "info-family" in html
