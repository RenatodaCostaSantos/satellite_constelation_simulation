"""Constelação de 18 satélites: P = 3·D no Cowell completo (card S6-01).

Referências: Apêndice A do prompt da Semana 6 (DESIGN_A2_S6, CONSTELLATION18_A1_S6,
COVERAGE18_A3_S6, CONTACTS18_A4_S6, MAINTENANCE18_S6), em tests/golden/reference_values.py.

Custo: a constelação de 18 no Cowell completo e os seus nós de 1 dia são fixtures de sessão do
conftest (≈ 20 s, uma vez); cobertura de 1 dia, contatos e manutenção reaproveitam a integração.
As janelas de 3 dias usam a mesma constelação no Cowell só J2–J4 (≈ 10 s), como no S5-03, em
grade de 2°: os valores de 3 dias com P = 18 são idênticos aos do Cowell completo (A.3).
"""

import math
import re
from datetime import datetime, timedelta

import numpy as np
import pytest
from golden.reference_values import (
    CONSTELLATION18_A1_S6,
    CONSTELLATION_A1,
    CONTACTS18_A4_S6,
    COVERAGE18_A3_S6,
    DESIGN_A2_S6,
    MAINTENANCE18_S6,
    NODE_GRID_A3,
    ROLL_DAILY_COVERAGE_TOL_DEG_S6,
)

from satsim.constants import SECONDS_PER_DAY
from satsim.constellation import (
    CONSTELLATION_DESIGNS,
    CONSTELLATION_PALETTE,
    DESIGN_91_6_P6,
    DESIGN_91_6_P18,
    Constellation,
    ConstellationDesign,
    constellation_palette,
    sso_constellation,
)
from satsim.coverage import Access, TargetGrid, access_times, coverage_stats
from satsim.elements import KeplerianElements
from satsim.forces import ZonalAcceleration
from satsim.groundtrack import grid_metrics, node_crossings
from satsim.maintenance import maintenance_status, measure_decay_rate
from satsim.observer import SAO_BENTO_DO_SAPUCAI
from satsim.orbit_design import constellation_phases, repeat_orbit
from satsim.passes import ContactPlan, contact_plan
from satsim.propagators import MeanJ2Propagator

A1, C3, K4, M18 = CONSTELLATION18_A1_S6, COVERAGE18_A3_S6, CONTACTS18_A4_S6, MAINTENANCE18_S6
EPOCH = datetime.fromisoformat(CONSTELLATION_A1.epoch_utc)
DAY = SECONDS_PER_DAY
A_DESIGN = CONSTELLATION_A1.a_design_m
I_A1 = math.radians(CONSTELLATION_A1.i_calibrated_deg)
ROLL_6 = math.radians(C3.roll_p6_deg)  # 23,2°: a rolagem de projeto da Semana 5
ROLL_18 = DESIGN_91_6_P18.roll_for_daily_coverage(A_DESIGN)  # 8,1835°
DESIGN_ROWS = {row.n_sats: row for row in DESIGN_A2_S6}


def _deg(x: float) -> float:
    return math.degrees(x)


def _family(sats, f: int, m: int = 3) -> Constellation:
    """A família f (s_{f+1}, s_{f+1+m}, …) como constelação."""
    members = [s for j, s in enumerate(sats) if j % m == f]
    return Constellation(EPOCH, members)


# --- Desenho (teoria: Sem4 §8.2–8.3, A.2) --------------------------------------------------------


def test_design_p18() -> None:
    d, row = DESIGN_91_6_P18, DESIGN_ROWS[18]
    assert (d.n_revs, d.n_days, d.n_sats) == (91, 6, 18)
    assert d.m == 3.0 and d.n_families == 3 and d.is_uniform is row.uniform
    assert d.tracks_per_day == row.tracks_per_day
    assert d.daily_grid_deg == pytest.approx(row.daily_grid_deg, abs=5e-6)
    assert d.daily_grid_deg == pytest.approx(d.grid_deg / 3)
    assert math.radians(d.daily_grid_deg) * 6378.137 == pytest.approx(row.daily_grid_km, abs=0.5)
    assert d.phase_step_deg == row.delta_u_deg
    assert d.full_grid_days == 1
    assert d.ids == tuple(f"s{j}" for j in range(1, 19))
    assert [d.family(j) for j in range(18)] == [0, 1, 2] * 6
    assert d.node_delay_s(5696.7) == pytest.approx(5696.7 / 18)
    with pytest.raises(ValueError):
        d.family(18)


def test_design_p6_is_the_week5_constellation() -> None:
    d, row = DESIGN_91_6_P6, DESIGN_ROWS[6]
    assert d.m == 1.0 and d.n_families == 1 and d.is_uniform
    assert d.tracks_per_day == row.tracks_per_day == 91
    assert d.daily_grid_deg == pytest.approx(row.daily_grid_deg, abs=5e-6) == d.grid_deg
    assert d.phase_step_deg == row.delta_u_deg
    assert [d.family(j) for j in range(6)] == [0] * 6


def test_design_p7_not_uniform() -> None:
    d = ConstellationDesign(91, 6, 7, "P = 7")
    assert d.is_uniform is DESIGN_ROWS[7].uniform is False
    assert d.tracks_per_day is None and d.daily_grid_deg is None and d.full_grid_days is None
    assert d.phase_step_deg == pytest.approx(DESIGN_ROWS[7].delta_u_deg, abs=0.005)
    assert d.family(6) == 0
    with pytest.raises(ValueError, match="D não divide P"):
        d.roll_for_daily_coverage(A_DESIGN)


def test_design_p3_uniform_by_p_divides_d() -> None:
    # P = 3: D | P falha, P | D vale → uniforme, grade completa a cada D/P = 2 dias
    d = ConstellationDesign(91, 6, 3, "P = 3")
    assert d.is_uniform and d.m == 0.5 and d.n_families == 1
    assert d.tracks_per_day is None and d.full_grid_days == 2


def test_design_rejects_invalid() -> None:
    with pytest.raises(ValueError):
        ConstellationDesign(91, 6, 0, "vazia")


@pytest.mark.parametrize("n_sats", [6, 18])
def test_roll_for_daily_coverage_a2(n_sats: int) -> None:
    """23,1957° com P = 6 e 8,1835° ± 0,01° com P = 18 (eq. 45, λ = π/(N·P/D))."""
    eta = CONSTELLATION_DESIGNS[n_sats].roll_for_daily_coverage(A_DESIGN)
    ref = DESIGN_ROWS[n_sats].roll_daily_coverage_deg
    assert _deg(eta) == pytest.approx(ref, abs=ROLL_DAILY_COVERAGE_TOL_DEG_S6)


def test_registry_is_the_single_source() -> None:
    assert CONSTELLATION_DESIGNS == {6: DESIGN_91_6_P6, 18: DESIGN_91_6_P18}
    assert all(d.n_sats == p and (d.n_revs, d.n_days) == (91, 6)
               for p, d in CONSTELLATION_DESIGNS.items())  # fmt: skip


# --- Paleta por família --------------------------------------------------------------------------


def test_palette_p6_identical_to_week5() -> None:
    assert constellation_palette(DESIGN_91_6_P6) == CONSTELLATION_PALETTE


def test_palette_p18_by_family() -> None:
    colors = constellation_palette(DESIGN_91_6_P18)
    assert len(colors) == 18 and len(set(colors)) == 18
    assert all(re.fullmatch(r"#[0-9a-f]{6}", c) for c in colors)
    # família 0 (s1, s4, …, s16) com as cores de S1 … S6 da Semana 5
    assert colors[0::3] == CONSTELLATION_PALETTE
    # mesma posição k = j // 3, mesmo matiz: a família 1 escurece e a 2 clareia a cor base
    for k, base in enumerate(CONSTELLATION_PALETTE):
        rgb = [
            np.array([int(c[i : i + 2], 16) for i in (1, 3, 5)]) for c in colors[3 * k : 3 * k + 3]
        ]
        assert np.all(rgb[1] <= rgb[0]) and np.all(rgb[2] >= rgb[0])
        assert colors[3 * k] == base


def test_palette_not_uniform_cycles_base_colors() -> None:
    colors = constellation_palette(ConstellationDesign(91, 6, 7, "P = 7"))
    assert colors == (*CONSTELLATION_PALETTE, CONSTELLATION_PALETTE[0])


def test_build_p6_identical_to_sso_constellation() -> None:
    """Compatibilidade: com P = 6, ids, nomes, cores e estados iniciais idênticos aos da S5."""
    zonal = (ZonalAcceleration(),)
    built = DESIGN_91_6_P6.build(EPOCH, zonal, i=I_A1)
    ref = sso_constellation(EPOCH, zonal, i=I_A1)
    assert built.ids == ref.ids == ["s1", "s2", "s3", "s4", "s5", "s6"]
    assert [s.color for s in built] == [s.color for s in ref] == list(CONSTELLATION_PALETTE)
    assert [s.name for s in built] == [s.name for s in ref]
    for a, b in zip(built, ref, strict=True):
        assert a.propagator.propagate(0.0)[0] == pytest.approx(b.propagator.propagate(0.0)[0])


def test_build_p18_recipe(constellation18_built) -> None:
    c, _ = constellation18_built
    assert c.ids == list(DESIGN_91_6_P18.ids) and c.epoch == EPOCH
    assert [s.color for s in c] == list(constellation_palette(DESIGN_91_6_P18))
    assert len(c[0].propagator.models) == 5  # Cowell completo
    assert c[0].cd == 2.2  # C_D e A/m copiados do arrasto, como na S5


def test_family0_starts_on_the_week5_constellation(constellation18, reference_constellation):
    """s1, s4, …, s16 têm os estados iniciais de S1 … S6 da Semana 5 (ψ0 = −j·20° = −k·60°)."""
    r18 = [s.propagator.propagate(0.0)[0] for s in constellation18]
    r6 = [s.propagator.propagate(0.0)[0] for s in reference_constellation]
    for k in range(6):
        assert r18[3 * k] == pytest.approx(r6[k], abs=1e-6)


# --- Modelo médio: P = 18 em 1 dia --------------------------------------------------------------


def test_mean_model_grid_p18() -> None:
    """MeanJ2, P = 18, 1 dia: 273 nós com passo uniforme de 1,31868° (Sem4 §8.3)."""
    a, i = repeat_orbit(91, 6)
    sats = [
        MeanJ2Propagator(KeplerianElements(a, 0.0, i, 0.3, 0.0, float(psi)), EPOCH)
        for psi in constellation_phases(18)
    ]
    start = EPOCH - timedelta(seconds=1)  # o nó de s1 em t = 0 entra
    lon = np.concatenate([node_crossings(s, start, DAY)[1] for s in sats])
    count, gap_min, gap_max, _ = grid_metrics(lon, 273)
    ideal = DESIGN_ROWS[18].daily_grid_deg
    assert count == A1.p18_nodes
    assert _deg(gap_min) == pytest.approx(ideal, abs=A1.mean_grid_tol_deg)
    assert _deg(gap_max) == pytest.approx(ideal, abs=A1.mean_grid_tol_deg)


# --- Cowell completo, 1 dia (A.1) ----------------------------------------------------------------


def test_cowell_grid_273_tracks(constellation18_nodes) -> None:
    """273 nós uniformes em 1 dia (âncora no 1º nó de s1): 1,31556 … 1,32196° ± 0,002°."""
    lon = np.concatenate([lo for _, lo in constellation18_nodes])
    count, gap_min, gap_max, offset = grid_metrics(lon, DESIGN_91_6_P18.tracks_per_day)
    assert count == A1.p18_nodes
    assert _deg(gap_min) == pytest.approx(A1.p18_gap_min_deg, abs=A1.gap_tol_deg)
    assert _deg(gap_max) == pytest.approx(A1.p18_gap_max_deg, abs=A1.gap_tol_deg)
    assert _deg(offset) <= A1.p18_max_offset_limit_deg
    assert _deg(offset) == pytest.approx(A1.p18_max_offset_deg, abs=5e-4)


def test_neighbour_node_delay(constellation18_nodes) -> None:
    """Cada satélite cruza o nó Tnod/18 ≈ 316,5 s ± 1 s depois do vizinho da frente."""
    delays = []
    for j in range(17):
        t_a, t_b = constellation18_nodes[j][0], constellation18_nodes[j + 1][0]
        for t in t_a:
            later = t_b[(t_b > t) & (t_b < t + 2000.0)]
            if later.size:
                delays.append(later[0] - t)
    delays = np.array(delays)
    assert delays.size > 200
    assert np.all(np.abs(delays - A1.node_delay_s2_s1_s) <= A1.node_delay_tol_s)
    s1, s2 = constellation18_nodes[0][0][0], constellation18_nodes[1][0][0]
    assert s2 - s1 == pytest.approx(A1.node_delay_s2_s1_s, abs=0.01)


def test_first_nodes_regression(constellation18_nodes) -> None:
    for (t, _), ref in zip(constellation18_nodes[:4], A1.first_node_s, strict=True):
        assert t[0] == pytest.approx(ref, abs=0.006)
    # s4 cruza o nó no mesmo instante que o S2 da constelação de 6 (A.1 da Semana 5)
    assert constellation18_nodes[3][0][0] == pytest.approx(
        CONSTELLATION_A1.first_node_s[1], abs=0.006
    )


def test_family0_reproduces_week5_grid(constellation18_nodes) -> None:
    """Família 0 (s1, s4, …, s16): 91 nós, 3,95173 … 3,96432° — o A.3 da Semana 5."""
    lon = np.concatenate([constellation18_nodes[j][1] for j in range(0, 18, 3)])
    count, gap_min, gap_max, offset = grid_metrics(lon, 91)
    assert count == A1.family0_nodes == NODE_GRID_A3.p6_nodes
    assert _deg(gap_min) == pytest.approx(A1.family0_gap_min_deg, abs=A1.gap_tol_deg)
    assert _deg(gap_max) == pytest.approx(A1.family0_gap_max_deg, abs=A1.gap_tol_deg)
    assert _deg(offset) <= NODE_GRID_A3.p6_max_offset_limit_deg


def test_families_are_shifted_by_a_third_of_the_grid(constellation18_nodes) -> None:
    """Cada família é uma grade de 91; as três se intercalam a 1/3 de passo (eq. 43)."""
    grid = math.radians(DESIGN_91_6_P18.grid_deg)
    s1_lon = constellation18_nodes[0][1][0]
    for f in range(3):
        lon = np.concatenate([constellation18_nodes[j][1] for j in range(f, 18, 3)])
        assert grid_metrics(lon, 91)[0] == 91
        # todos os nós da família f ficam f/3 de passo a oeste da grade da família 0
        rel = np.mod(lon - s1_lon, grid) / grid
        expected = np.mod(-f / 3, 1.0)
        dist = np.abs(np.mod(rel - expected + 0.5, 1.0) - 0.5)
        assert dist.max() < 0.01


# --- Cobertura (grade 2°, passo 10 s, A.3) -------------------------------------------------------


@pytest.fixture(scope="module")
def grid2() -> TargetGrid:
    return TargetGrid(C3.grid_deg, C3.grid_deg)


def test_daily_coverage_p18_with_8_degrees(constellation18, grid2) -> None:
    """P = 18 com ±8,18°: 100% dos alvos em 1 dia (a grade de 273 trilhas fecha)."""
    access = access_times(constellation18, grid2, EPOCH, DAY, ROLL_18, step_s=C3.step_s)
    frac, _ = coverage_stats(access, DAY)
    assert frac == C3.fraction_p18_roll8_1d


@pytest.fixture(scope="module")
def access_p18_roll23_1d(constellation18, grid2) -> Access:
    return access_times(constellation18, grid2, EPOCH, DAY, ROLL_6, step_s=C3.step_s)


def test_redundancy_17_satellites(constellation18, grid2) -> None:
    """Sem s1, os 17 restantes com ±23,2° ainda cobrem 100% em 1 dia."""
    c17 = Constellation(EPOCH, list(constellation18)[1:])
    frac, _ = coverage_stats(access_times(c17, grid2, EPOCH, DAY, ROLL_6, step_s=C3.step_s), DAY)
    assert len(c17) == 17 and frac == C3.fraction_p17_roll23_1d


def _passes_per_target(access: Access, gap_merge_s: float) -> np.ndarray:
    """Passagens distintas por alvo (acessos a menos de ``gap_merge_s`` contam como uma)."""
    idx, t = access.target_index, access.t_s
    new = np.r_[True, idx[1:] != idx[:-1]] | (np.r_[np.inf, np.diff(t)] >= gap_merge_s)
    return np.bincount(idx[new], minlength=access.n_targets)


def test_passes_per_target_p18(access_p18_roll23_1d, grid2) -> None:
    """Com P = 18 e 23,2°, cada alvo em |lat| ≤ 60° tem 2 a 3 passagens por dia, não 3 × 2:
    os satélites vizinhos do trem chegam minutos depois e se fundem na mesma passagem."""
    n = _passes_per_target(access_p18_roll23_1d, C3.gap_merge_s)
    lat60 = np.abs(grid2.lat_deg) <= 60.0
    assert n[lat60].min() == C3.passes_per_target_1d_min
    assert np.median(n[lat60]) == C3.passes_per_target_1d_median
    assert n[lat60].max() == C3.passes_per_target_1d_max


@pytest.fixture(scope="module")
def zonal18() -> Constellation:
    """A constelação de 18 no Cowell só J2–J4: janelas de 3 dias a baixo custo."""
    return DESIGN_91_6_P18.build(EPOCH, (ZonalAcceleration(),), i=I_A1)


def test_revisit_p18_barely_changes(zonal18, grid2) -> None:
    """Revisita máxima em |lat| ≤ 60° (3 dias, 23,2°): 13,79 h ± 0,5 h com 18 satélites, a menos
    de 0,5 h da de 6 (13,97 h): os satélites extras chegam no mesmo trem."""
    lat60 = np.abs(grid2.lat_deg) <= 60.0
    revisit = {}
    for label, cons in (("p18", zonal18), ("p6", _family(zonal18, 0))):
        access = access_times(cons, grid2, EPOCH, 3 * DAY, ROLL_6, step_s=C3.step_s)
        frac, gaps = coverage_stats(access, 3 * DAY, gap_merge_s=C3.gap_merge_s)
        assert frac == 1.0
        revisit[label] = gaps[lat60].max() / 3600
    assert revisit["p18"] == pytest.approx(C3.revisit_p18_roll23_3d_lat60_h, abs=C3.revisit_tol_h)
    assert revisit["p6"] == pytest.approx(C3.revisit_p6_roll23_3d_lat60_h, abs=C3.revisit_tol_h)
    assert abs(revisit["p18"] - revisit["p6"]) < C3.revisit_tol_h
    assert revisit["p18"] <= revisit["p6"]


def test_p6_with_8_degrees_leaves_gaps(zonal18, grid2) -> None:
    """P = 6 com ±8,18° em 3 dias: 74,35% ± 2% (a faixa estreita não fecha a grade de 91)."""
    access = access_times(_family(zonal18, 0), grid2, EPOCH, 3 * DAY, ROLL_18, step_s=C3.step_s)
    frac, _ = coverage_stats(access, 3 * DAY)
    assert frac == pytest.approx(C3.fraction_p6_roll8_3d, abs=C3.fraction_p6_roll8_3d_tol)


# --- Contatos com São Bento do Sapucaí (A.4) -----------------------------------------------------


@pytest.fixture(scope="module")
def plan18(constellation18) -> ContactPlan:
    return contact_plan(constellation18, SAO_BENTO_DO_SAPUCAI, EPOCH, DAY, 10.0)


def _max_simultaneous(plan: ContactPlan) -> int:
    """Maior número de satélites visíveis ao mesmo tempo (varredura de AOS/LOS)."""
    events = sorted([(p.aos, 1) for p in plan.passes] + [(p.los, -1) for p in plan.passes],
                    key=lambda e: (e[0], e[1]))  # fmt: skip
    current = peak = 0
    for _, step in events:
        current += step
        peak = max(peak, current)
    return peak


def test_contacts_18_counts(plan18) -> None:
    assert len(plan18.contacts) == pytest.approx(K4.n_passes, abs=K4.n_passes_tol)
    assert plan18.total_contact_s / 60 == pytest.approx(
        K4.total_contact_min, abs=K4.total_contact_tol_min
    )
    assert len(plan18.overlaps) == pytest.approx(K4.n_overlaps, abs=K4.n_overlaps_tol)
    assert _max_simultaneous(plan18) == K4.max_simultaneous
    lo, hi = K4.passes_per_sat_range
    assert all(lo <= n <= hi for n in plan18.passes_per_satellite.values())
    assert len(plan18.passes_per_satellite) == 18
    assert plan18.max_gap_s / 3600 == pytest.approx(K4.max_gap_h, abs=0.05)


def test_contacts_18_trains(plan18) -> None:
    """Dois trens (noite e dia), AOS consecutivos a 4–6,5 min (≈ Tnod/18 = 5,3 min)."""
    assert len(plan18.trains) == K4.n_trains
    night, day = plan18.trains
    assert (night.lighting, day.lighting) == ("noite", "dia")
    assert len(night.contact_indices) == K4.night_train_passes
    assert len(day.contact_indices) == K4.day_train_passes
    for train, (t0, t1) in ((night, K4.night_train_utc), (day, K4.day_train_utc)):
        for got, ref in ((train.start, t0), (train.end, t1)):
            h, m = (int(x) for x in ref.split(":"))
            assert abs((got - EPOCH.replace(hour=h, minute=m)).total_seconds()) < 60
    lo, hi = K4.train_aos_spacing_min
    for train in plan18.trains:
        aos = [plan18.contacts[k].pass_.aos for k in train.contact_indices]
        spacing = np.diff([(t - EPOCH).total_seconds() for t in aos]) / 60
        assert np.all((spacing >= lo) & (spacing <= hi)), spacing


# --- Manutenção ---------------------------------------------------------------------------------


def test_maintenance_18_equals_6(constellation18, reference_constellation) -> None:
    """Mesmo a, mesmo arrasto: ȧ e ``MaintenanceStatus`` dos 18 iguais aos dos 6 (1%, 1 dia)."""
    rates18 = np.array([measure_decay_rate(s, 1.0) for s in constellation18])
    rates6 = np.array([measure_decay_rate(s, 1.0) for s in reference_constellation])
    ref = rates6.mean()
    assert np.all(np.abs(rates18 / ref - 1.0) <= M18.rel_tol)
    assert (rates18 * DAY).min() == pytest.approx(M18.a_dot_min_m_day, abs=0.01)
    assert (rates18 * DAY).max() == pytest.approx(M18.a_dot_max_m_day, abs=0.01)
    status6 = maintenance_status("s1", ref, A_DESIGN, I_A1)
    for sat, rate in zip(constellation18, rates18, strict=True):
        st = maintenance_status(sat.id, rate, A_DESIGN, I_A1)
        for field in ("days_to_track_correction", "track_delta_v_per_year_m_s",
                      "days_to_ltan_correction", "inclination_delta_v_per_year_m_s"):  # fmt: skip
            assert getattr(st, field) == pytest.approx(getattr(status6, field), rel=M18.rel_tol)
