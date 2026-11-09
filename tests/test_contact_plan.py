"""Plano de contatos da constelação sobre São Bento do Sapucaí (card S5-04).

Referência: Apêndice A.6 do prompt da Semana 5 (CONTACT_PASSES_A6, CONTACT_PLAN_A6 em
tests/golden/reference_values.py). A constelação é a de referência do conftest (Cowell completo,
integrado uma vez por sessão): o plano de 24 h reaproveita o cache.
"""

import importlib.util
import urllib.request
from datetime import UTC, datetime, timedelta
from pathlib import Path

import numpy as np
import pytest
from golden.reference_values import CONTACT_PASSES_A6, CONTACT_PLAN_A6

from satsim.astrotime import datetime_to_jd
from satsim.constellation import Constellation
from satsim.ephemeris import sun_position_eci, sun_ra_dec
from satsim.frames import eci_to_ecef, gmst_at
from satsim.observer import SAO_BENTO_DO_SAPUCAI
from satsim.passes import (
    contact_plan,
    find_passes,
    format_contact_plan,
    sun_elevation_deg,
)
from satsim.satellite import Satellite

ROOT = Path(__file__).resolve().parents[1]
A6 = CONTACT_PLAN_A6
START = datetime.fromisoformat(A6.start_utc)
STATION = SAO_BENTO_DO_SAPUCAI
BRT = timedelta(hours=-3)


@pytest.fixture(scope="module")
def plan(reference_constellation):
    return contact_plan(reference_constellation, STATION, START, A6.duration_s, A6.el_min_deg)


def _utc(hhmmss: str) -> datetime:
    h, m, *s = (int(x) for x in hhmmss.split(":"))
    return START.replace(hour=h, minute=m, second=s[0] if s else 0)


# --- Critérios do card (A.6) ---------------------------------------------------------------------


def test_plan_counts_contact_time_and_overlaps(plan) -> None:
    assert len(plan.contacts) == pytest.approx(A6.n_passes, abs=A6.n_passes_tol)
    assert plan.total_contact_s / 60 == pytest.approx(
        A6.total_contact_min, abs=A6.total_contact_tol_min
    )
    assert len(plan.overlaps) == A6.n_overlaps
    per_sat = plan.passes_per_satellite
    assert [per_sat.get(f"s{j}", 0) for j in range(1, 7)] == list(A6.passes_per_sat)
    assert max(p.max_elevation_deg for p in plan.passes) == pytest.approx(
        A6.max_elevation_deg, abs=0.01
    )


def test_two_trains_night_and_day(plan) -> None:
    """Dois trens: noturno ≈ 21h30–23h30 BRT e diurno ≈ 9h–11h BRT (A.6)."""
    assert len(plan.trains) == A6.n_trains
    night, day = plan.trains
    assert night.lighting == "noite" and day.lighting == "dia"
    for train, (t0, t1) in ((night, A6.night_train_utc), (day, A6.day_train_utc)):
        assert abs((train.start - _utc(t0)).total_seconds()) < 60
        assert abs((train.end - _utc(t1)).total_seconds()) < 60
    # em BRT: 21h32–23h25 e 9h05–10h58
    assert (night.start + BRT).hour == 21 and (night.end + BRT).hour == 23
    assert (day.start + BRT).hour == 9 and (day.end + BRT).hour == 10


def test_aos_spacing_inside_trains(plan) -> None:
    """AOS consecutivos de um trem separados por 14–17 min (≈ Tnod/6 = 15,8 min + geometria)."""
    lo, hi = A6.train_aos_spacing_min
    for train in plan.trains:
        aos = [plan.contacts[k].pass_.aos for k in train.contact_indices]
        spacing = np.diff([(t - START).total_seconds() for t in aos]) / 60
        assert np.all((spacing >= lo) & (spacing <= hi)), spacing
        # o trem passa em fila: cada satélite segue o anterior (s_j → s_{j+1}, com volta de s6 a s1)
        ids = [int(plan.contacts[k].pass_.sat_id[1:]) for k in train.contact_indices]
        assert all((b - a) % 6 == 1 for a, b in zip(ids, ids[1:], strict=False))


def test_passes_match_a6_table(plan) -> None:
    # Regressão contra a tabela do A.6 (AOS no segundo, el máx com 2 casas, duração com 1)
    assert len(plan.contacts) == len(CONTACT_PASSES_A6)
    for c, ref in zip(plan.contacts, CONTACT_PASSES_A6, strict=True):
        p = c.pass_
        assert p.sat_id == ref.sat.lower()
        assert abs((p.aos - _utc(ref.aos_utc)).total_seconds()) < 1.0
        assert p.max_elevation_deg == pytest.approx(ref.max_el_deg, abs=0.01)
        assert p.duration_s == pytest.approx(ref.duration_s, abs=0.1)


def test_max_gap_without_contact(plan) -> None:
    # Bordas incluídas: do último LOS (13:58 UTC) ao fim da janela, ≈ 10,03 h; entre os trens 9,66 h
    night, day = plan.trains
    between = (day.start - night.end).total_seconds()
    tail = (START + timedelta(seconds=A6.duration_s) - day.end).total_seconds()
    assert plan.max_gap_s == pytest.approx(max(between, tail, (night.start - START).seconds))
    assert plan.max_gap_s / 3600 == pytest.approx(10.03, abs=0.01)


# --- Coerência ----------------------------------------------------------------------------------


def test_plan_passes_equal_find_passes_per_satellite(plan, reference_constellation) -> None:
    for sat in reference_constellation:
        alone = find_passes(sat.propagator, STATION, START, A6.duration_s, A6.el_min_deg,
                            sat_id=sat.id)  # fmt: skip
        assert [p for p in plan.passes if p.sat_id == sat.id] == alone


def test_day_night_matches_ephemeris(plan) -> None:
    # Elevação do Sol independente: posição do Sol (ECI) → ECEF → look_angles da estação
    for c in plan.contacts:
        tca = c.pass_.tca
        sun_eci = sun_position_eci(datetime_to_jd(tca))
        _, el, _ = STATION.look_angles(eci_to_ecef(sun_eci, gmst_at(tca, 0.0)))
        assert c.sun_elevation_deg == pytest.approx(np.degrees(float(el)), abs=0.01)
        assert c.lighting == ("dia" if c.sun_elevation_deg >= -6.0 else "noite")
    # noturno às ≈ 22h locais (Sol bem abaixo), diurno às ≈ 10h (Sol alto)
    night, day = plan.trains
    assert all(plan.contacts[k].sun_elevation_deg < -40 for k in night.contact_indices)
    assert all(plan.contacts[k].sun_elevation_deg > 45 for k in day.contact_indices)


def test_sun_elevation_noon_and_midnight() -> None:
    # Meio-dia solar em São Bento (λ = −45,73°): 12h + 3h03 − equação do tempo (≈ 16 min em
    # outubro) ≈ 14h47 UTC. Elevação máxima 90° − |φ − δ|; à meia-noite, −(90° − |φ + δ|).
    noon = datetime(2026, 10, 23, 14, 47, tzinfo=UTC)
    el = sun_elevation_deg(STATION, [noon, noon + timedelta(hours=12)])
    _, dec = sun_ra_dec(datetime_to_jd(noon))
    phi, delta = STATION.lat_deg, np.degrees(float(dec))
    assert el[0] == pytest.approx(90.0 - abs(phi - delta), abs=0.2)
    assert el[1] == pytest.approx(-(90.0 - abs(phi + delta)), abs=0.3)
    assert sun_elevation_deg(STATION, []).size == 0
    with pytest.raises(ValueError):
        sun_elevation_deg(STATION, [datetime(2026, 10, 23)])


# --- Sobreposições e casos de borda ---------------------------------------------------------------


def test_overlap_detected_for_coincident_satellites(reference_constellation) -> None:
    # Dois "satélites" com o mesmo propagador: cada passagem de um sobrepõe a do outro
    s1 = reference_constellation[0]
    twins = Constellation(START, [s1, Satellite("s1b", "S1b", s1.propagator)])
    plan = contact_plan(twins, STATION, START, 6 * 3600.0)  # S1 passa uma vez, às 01:02 UTC
    assert len(plan.contacts) == 2 and len(plan.overlaps) == 1
    ov = plan.overlaps[0]
    p = plan.contacts[0].pass_
    assert {ov.sat_a, ov.sat_b} == {"s1", "s1b"} and (ov.start, ov.end) == (p.aos, p.los)
    assert len(plan.trains) == 1  # as 2 passagens com o mesmo AOS


def test_empty_plan(reference_constellation) -> None:
    # 20 min sem passagem (00:00–00:20 UTC): plano vazio, intervalo = janela inteira
    plan = contact_plan(reference_constellation, STATION, START, 1200.0)
    assert plan.contacts == () and plan.trains == () and plan.overlaps == ()
    assert plan.total_contact_s == 0.0 and plan.max_gap_s == 1200.0
    assert "nenhuma passagem" in format_contact_plan(plan)


def test_format_contact_plan(plan) -> None:
    text = format_contact_plan(plan)
    assert "AOS BRT" in text and "trem" in text
    assert text.count("noite") >= 8 and text.count(" dia") >= 8
    assert "trem 1: 8 passagens, 21:32–23:25 BRT (noite)" in text
    assert "trem 2: 8 passagens, 09:05–10:58 BRT (dia)" in text
    assert "contato total: 96.7 min em 16 passagens" in text
    assert "sobreposições (duas antenas): 0" in text


# --- scripts/list_passes.py --constellation ------------------------------------------------------


def _load_script():
    path = ROOT / "scripts" / "list_passes.py"
    spec = importlib.util.spec_from_file_location("list_passes_s5", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_list_passes_constellation(monkeypatch, capsys, tmp_path, reference_constellation):
    def no_network(*args, **kwargs):
        raise AssertionError("o script tentou acessar a rede")

    monkeypatch.setattr(urllib.request, "urlopen", no_network)
    script = _load_script()
    # a constelação do script é a mesma da referência (mesma receita); reaproveita o cache
    monkeypatch.setattr(script, "build_constellation", lambda start: reference_constellation)
    out_csv = tmp_path / "c.csv"
    script.main(["--constellation", "--start", "2026-10-23T00:00:00Z", "--days", "1",
                 "--csv", str(out_csv)])  # fmt: skip
    out = capsys.readouterr().out
    assert "constelação: 91/6" in out
    assert "trem 2: 8 passagens, 09:05–10:58 BRT (dia)" in out
    assert out_csv.read_text().count("\n") == 1 + A6.n_passes
    with pytest.raises(SystemExit):
        script.main(["--constellation", "--start", "tle-epoch"])


def test_build_constellation_recipe() -> None:
    script = _load_script()
    c = script.build_constellation(START)
    assert c.ids == ["s1", "s2", "s3", "s4", "s5", "s6"] and c.epoch == START
    assert len(c[0].propagator.models) == 5
