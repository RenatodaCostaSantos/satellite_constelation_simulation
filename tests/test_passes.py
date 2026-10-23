"""Testes do detector de passagens (S3-05).

Referências: tests/golden/reference_values.py — PASSES_A2 e PASSES_A2_META (Apêndice A.2 do
prompt da Semana 3; AOS da passagem 9 corrigido), ISS_TLE_A1 (A.1). Nenhum teste acessa a rede.
"""

import importlib.util
import math
import urllib.request
from datetime import UTC, datetime, timedelta
from pathlib import Path

import numpy as np
import pytest
from golden.reference_values import ISS_TLE_A1, PASSES_A2, PASSES_A2_META

from satsim.constants import R_EARTH
from satsim.elements import KeplerianElements
from satsim.observer import SAO_BENTO_DO_SAPUCAI, observe
from satsim.passes import Pass, find_passes, find_passes_many, format_passes
from satsim.propagators import MeanJ2Propagator, SGP4Propagator
from satsim.tle import Tle, tle_checksum

ROOT = Path(__file__).resolve().parents[1]
SAMPLE_TEXT = (Path(__file__).parent / "data" / "iss_sample.tle").read_text()
M = PASSES_A2_META
WINDOW = timedelta(days=M.window_days)


@pytest.fixture(scope="module")
def iss() -> SGP4Propagator:
    return SGP4Propagator.from_tle_text(SAMPLE_TEXT)


@pytest.fixture(scope="module")
def passes10(iss) -> list[Pass]:
    return find_passes(iss, SAO_BENTO_DO_SAPUCAI, iss.epoch, WINDOW, el_min_deg=10.0, sat_id="ISS")


@pytest.fixture(scope="module")
def passes5(iss) -> list[Pass]:
    return find_passes(iss, SAO_BENTO_DO_SAPUCAI, iss.epoch, WINDOW, el_min_deg=5.0, sat_id="ISS")


def seconds(p_time: datetime, epoch: datetime) -> float:
    return (p_time - epoch).total_seconds()


def _az_diff(a: float, b: float) -> float:
    return (a - b + 180.0) % 360.0 - 180.0


# --- Tabela de referência (A.2) ----------------------------------------------------------------


def test_reference_table_a2(iss, passes10):
    assert len(passes10) == M.n_passes_10deg == len(PASSES_A2)
    for p, ref in zip(passes10, PASSES_A2, strict=True):
        assert seconds(p.aos, iss.epoch) == pytest.approx(ref.aos_s, abs=M.tol_time_s)
        assert seconds(p.tca, iss.epoch) == pytest.approx(ref.tca_s, abs=M.tol_time_s)
        assert seconds(p.los, iss.epoch) == pytest.approx(ref.los_s, abs=M.tol_time_s)
        assert p.duration_s == pytest.approx(ref.duration_s, abs=M.tol_time_s)
        assert p.max_elevation_deg == pytest.approx(ref.max_el_deg, abs=M.tol_el_deg)
        for az, az_ref in ((p.az_aos_deg, ref.az_aos_deg), (p.az_tca_deg, ref.az_tca_deg),
                           (p.az_los_deg, ref.az_los_deg)):  # fmt: skip
            assert abs(_az_diff(az, az_ref)) <= M.tol_az_deg
        assert p.sat_id == "ISS" and not p.starts_inside and not p.ends_inside


def test_times_are_utc_aware_microseconds(passes10):
    for p in passes10:
        for dt in (p.aos, p.tca, p.los):
            assert dt.tzinfo is not None and dt.utcoffset() == timedelta(0)
        assert p.duration_s == pytest.approx((p.los - p.aos).total_seconds(), abs=1e-6)


# --- Contagem e contenção ----------------------------------------------------------------------


def test_count_5_and_10_degrees_and_containment(passes5, passes10):
    assert len(passes5) == M.n_passes_5deg
    assert len(passes10) == M.n_passes_10deg
    for p in passes10:  # cada passagem a 10° está dentro de uma a 5°
        assert any(q.aos <= p.aos and p.los <= q.los for q in passes5)


def test_sorted_no_overlap_positive_duration(passes5):
    for p in passes5:
        assert p.duration_s > 0.0 and p.aos < p.tca < p.los
    for a, b in zip(passes5, passes5[1:], strict=False):
        assert a.los < b.aos


# --- Consistência ------------------------------------------------------------------------------


def test_elevation_at_aos_los_and_maximum_at_tca(iss, passes10):
    for p in passes10:
        t = np.array([seconds(x, iss.epoch) for x in (p.aos, p.los)])
        el = observe(iss, SAO_BENTO_DO_SAPUCAI, t).el_deg
        np.testing.assert_allclose(el, 10.0, rtol=0, atol=1e-3)
        t_tca = seconds(p.tca, iss.epoch)
        el_tca, el_before, el_after = observe(
            iss, SAO_BENTO_DO_SAPUCAI, [t_tca, t_tca - 2.0, t_tca + 2.0]
        ).el_deg
        assert el_before < el_tca and el_after < el_tca
        # o TCA guardado é arredondado ao microssegundo: ~1e-9° de diferença ao reavaliar
        assert p.max_elevation_deg == pytest.approx(el_tca, abs=1e-6)


def test_step_60_gives_same_passes(iss, passes10):
    coarse = find_passes(iss, SAO_BENTO_DO_SAPUCAI, iss.epoch, WINDOW, step_s=60.0, sat_id="ISS")
    assert len(coarse) == len(passes10)
    for a, b in zip(coarse, passes10, strict=True):
        assert abs((a.aos - b.aos).total_seconds()) < 1e-2
        assert abs((a.tca - b.tca).total_seconds()) < 1e-1
        assert abs((a.los - b.los).total_seconds()) < 1e-2


# --- Bordas da janela --------------------------------------------------------------------------


def test_window_starting_inside_a_pass(iss, passes10):
    p1 = passes10[0]
    start = p1.tca
    got = find_passes(iss, SAO_BENTO_DO_SAPUCAI, start, timedelta(hours=2))
    assert got[0].starts_inside and got[0].aos == start
    assert abs((got[0].los - p1.los).total_seconds()) < 1e-2
    assert not got[0].ends_inside


def test_window_ending_inside_a_pass(iss, passes10):
    p2 = passes10[1]
    end = p2.aos + (p2.los - p2.aos) / 2
    got = find_passes(iss, SAO_BENTO_DO_SAPUCAI, iss.epoch, end - iss.epoch)
    assert len(got) == 2
    assert got[-1].ends_inside and got[-1].los == end
    assert not got[-1].starts_inside


def test_window_entirely_inside_a_pass(iss, passes10):
    p1 = passes10[0]
    got = find_passes(iss, SAO_BENTO_DO_SAPUCAI, p1.aos + timedelta(seconds=30), 60.0)
    assert len(got) == 1 and got[0].starts_inside and got[0].ends_inside


def test_rejects_naive_start_and_bad_step(iss):
    with pytest.raises(ValueError, match="tz-aware"):
        find_passes(iss, SAO_BENTO_DO_SAPUCAI, datetime(2019, 12, 10), 3600.0)
    with pytest.raises(ValueError):
        find_passes(iss, SAO_BENTO_DO_SAPUCAI, iss.epoch, 3600.0, step_s=0.0)


# --- Skyfield (oráculo independente) -----------------------------------------------------------


@pytest.mark.parametrize("el_min", [10.0, 5.0])
def test_matches_skyfield_find_events(iss, el_min):
    skyfield_api = pytest.importorskip("skyfield.api")
    ts = skyfield_api.load.timescale(builtin=True)
    sat = skyfield_api.EarthSatellite(ISS_TLE_A1.line1, ISS_TLE_A1.line2, ISS_TLE_A1.name, ts)
    s = SAO_BENTO_DO_SAPUCAI
    topos = skyfield_api.wgs84.latlon(s.lat_deg, s.lon_deg, elevation_m=s.alt)
    t0, t1 = ts.from_datetime(iss.epoch), ts.from_datetime(iss.epoch + WINDOW)
    t, events = sat.find_events(topos, t0, t1, altitude_degrees=el_min)
    assert list(events) == [0, 1, 2] * (len(events) // 3)  # AOS, TCA, LOS completos
    sky = np.array([seconds(x, iss.epoch) for x in t.utc_datetime()]).reshape(-1, 3)
    mine = find_passes(iss, SAO_BENTO_DO_SAPUCAI, iss.epoch, WINDOW, el_min_deg=el_min)
    assert len(mine) == len(sky)
    ours = np.array([[seconds(x, iss.epoch) for x in (p.aos, p.tca, p.los)] for p in mine])
    # critério do card: < 1 s; referência A.2: ≤ 0,19 s
    assert np.max(np.abs(ours - sky)) < 1.0
    assert np.max(np.abs(ours - sky)) <= M.skyfield_max_dt_s + 0.01


# --- Vários satélites e outros propagadores ----------------------------------------------------


def _shifted_tle(days: float) -> Tle:
    """Cópia do TLE de teste com a época deslocada de ``days`` (checksum recalculado)."""
    tle = Tle.from_text(SAMPLE_TEXT)
    day = float(tle.line1[20:32]) + days
    line1 = tle.line1[:20] + f"{day:012.8f}" + tle.line1[32:]
    return Tle("ISS deslocada", line1[:68] + str(tle_checksum(line1)), tle.line2)


def test_find_passes_many_sorted_with_sat_id(iss):
    shifted = SGP4Propagator(_shifted_tle(0.25))
    many = find_passes_many({"iss": iss, "iss-6h": shifted}, SAO_BENTO_DO_SAPUCAI, iss.epoch,
                            WINDOW)  # fmt: skip
    assert [p.aos for p in many] == sorted(p.aos for p in many)
    for sat_id, prop in (("iss", iss), ("iss-6h", shifted)):
        alone = find_passes(prop, SAO_BENTO_DO_SAPUCAI, iss.epoch, WINDOW, sat_id=sat_id)
        assert [p for p in many if p.sat_id == sat_id] == alone
        assert alone
    assert {p.sat_id for p in many} == {"iss", "iss-6h"}


def test_find_passes_with_mean_j2_propagator(iss):
    el = KeplerianElements(a=R_EARTH + 420e3, e=0.0007, i=math.radians(51.64), raan=1.0,
                           argp=0.3, psi=0.0)  # fmt: skip
    prop = MeanJ2Propagator(el, datetime(2026, 10, 23, tzinfo=UTC))
    passes = find_passes(prop, SAO_BENTO_DO_SAPUCAI, prop.epoch, timedelta(days=2))
    assert passes
    for p in passes:
        assert 60.0 <= p.duration_s <= 12 * 60.0


# --- Formatação --------------------------------------------------------------------------------


def test_format_passes_utc_and_brt(passes10):
    text = format_passes(passes10)
    first = text.splitlines()[2]
    # 1ª passagem do A.2: AOS 10/12 00:46:08 UTC → 09/12 21:46:08 BRT (UTC − 3 h)
    assert "2019-12-10 00:46:08" in first and "21:46:08" in first
    assert "BRT" in text.splitlines()[0]
    assert len(text.splitlines()) == 2 + len(passes10)
    assert "UTC+1" in format_passes(passes10, tz_offset_hours=1.0).splitlines()[0]
    assert "nenhuma passagem" in format_passes([])


# --- Script (offline) --------------------------------------------------------------------------


def _load_script():
    path = ROOT / "scripts" / "list_passes.py"
    spec = importlib.util.spec_from_file_location("list_passes", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_list_passes_script_offline(tmp_path, monkeypatch, capsys):
    def no_network(*args, **kwargs):
        raise AssertionError("o script tentou acessar a rede")

    monkeypatch.setattr(urllib.request, "urlopen", no_network)
    script = _load_script()
    out_csv = tmp_path / "p.csv"
    script.main(["--offline", "--cache-dir", str(tmp_path), "--days", "1", "--csv", str(out_csv)])
    out = capsys.readouterr().out
    assert "fonte: bundled" in out and "idade:" in out and "AOS UTC" in out
    assert out_csv.read_text().startswith("sat_id,aos_utc")


def test_list_passes_script_reproduces_a2(capsys):
    script = _load_script()
    script.main(["--tle-file", str(Path(__file__).parent / "data" / "iss_sample.tle"),
                 "--start", "tle-epoch"])  # fmt: skip
    out = capsys.readouterr().out
    assert "idade:    0.00 dias" in out and "AVISO" not in out
    assert out.count("ISS (ZARYA)  2019-12") == M.n_passes_10deg
