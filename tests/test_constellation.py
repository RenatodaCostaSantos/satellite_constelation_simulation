"""Satellite, Constellation, repeat_orbit e a constelação 91/6 (card S5-01).

Referências: Apêndice A.1 e A.2 do prompt da Semana 5 (tests/golden/reference_values.py) e
Sem4 Tabelas 3 e 4. A constelação no Cowell completo é a de referência do conftest (uma vez por
sessão, i dada, ≈ 0,3 s); aqui ela é propagada por duas voltas (≈ 1 s).
"""

import dataclasses
import math
from datetime import datetime, timedelta
from types import SimpleNamespace

import numpy as np
import pytest
from golden.reference_values import (
    CONSTELLATION_A1,
    REPEAT_ORBIT_TOL_A2,
    REPEAT_ORBITS_A2,
    SSO_REPEAT_TABLE3,
)

import satsim.constellation as constellation_module
from satsim.analysis import osculating_elements
from satsim.astrotime import datetime_to_jd
from satsim.constants import R_EARTH
from satsim.constellation import CONSTELLATION_PALETTE, Constellation, sso_constellation
from satsim.elements import KeplerianElements, mean_to_true, state_to_elements
from satsim.forces import (
    J2Acceleration,
    ZonalAcceleration,
)
from satsim.mean_osc import window_mean_elements
from satsim.orbit_design import (
    constellation_phases,
    frozen_eccentricity,
    raan_for_ltan,
    repeat_orbit,
)
from satsim.propagators import CowellPropagator, KeplerPropagator
from satsim.satellite import DEFAULT_COLOR, Satellite
from satsim.secular import nodal_period, zonal_mean_rates

A1 = CONSTELLATION_A1
EPOCH = datetime.fromisoformat(A1.epoch_utc)
I_A1 = math.radians(A1.i_calibrated_deg)
LEO = KeplerianElements(7000e3, 0.001, math.radians(98.0), 0.3, 0.5, 0.0)


@pytest.fixture(scope="module")
def built(reference_constellation_built) -> tuple[Constellation, float]:
    """A constelação de referência (A.1, i dada; conftest, por sessão) e o tempo de construção."""
    return reference_constellation_built


@pytest.fixture(scope="module")
def constellation(reference_constellation) -> Constellation:
    return reference_constellation


@pytest.fixture(scope="module")
def first_nodes(constellation) -> list[np.ndarray]:
    """Instantes dos nós ascendentes (z_ECI = 0, vz > 0) de cada satélite em 2 voltas [s]."""
    t = np.arange(0.0, 11400.0 + 1.0, 10.0)
    r, _ = constellation.propagate(t)
    nodes = []
    for z in r[:, :, 2]:
        k = np.flatnonzero((z[:-1] < 0.0) & (z[1:] >= 0.0))
        nodes.append(t[k] - z[k] * (t[k + 1] - t[k]) / (z[k + 1] - z[k]))
    return nodes


# --- repeat_orbit: Sem4 Tabela 3 e A.2 ----------------------------------------------------------


@pytest.mark.parametrize("row", REPEAT_ORBITS_A2, ids=lambda r: f"{r.n_revs}/{r.d_days}")
def test_repeat_orbit_reproduces_a2(row) -> None:
    tol = REPEAT_ORBIT_TOL_A2
    a, i = repeat_orbit(row.n_revs, row.d_days)
    assert (a - R_EARTH) / 1e3 == pytest.approx(row.h_km, abs=tol.h_km)
    assert math.degrees(i) == pytest.approx(row.i_deg, abs=tol.i_deg)
    assert nodal_period(a, 0.0, i) == pytest.approx(row.tnod_s, abs=tol.tnod_s)
    # a grade só depende de N (Sem4 eq. 39)
    assert 360.0 / row.n_revs == pytest.approx(row.grid_deg, abs=5e-4)
    assert 2.0 * math.pi * R_EARTH / 1e3 / row.n_revs == pytest.approx(row.grid_km, abs=0.05)


@pytest.mark.parametrize("row", SSO_REPEAT_TABLE3, ids=lambda r: f"{r.n_revs}/{r.d_days}")
def test_repeat_orbit_reproduces_sem4_table3(row) -> None:
    # Inclui as linhas degeneradas (gcd > 1): a órbita existe, só repete antes
    a, i = repeat_orbit(row.n_revs, row.d_days)
    assert (a - R_EARTH) / 1e3 == pytest.approx(row.h_km, abs=REPEAT_ORBIT_TOL_A2.h_km)
    assert math.degrees(i) == pytest.approx(row.i_deg, abs=REPEAT_ORBIT_TOL_A2.i_deg)


def test_repeat_orbit_design_matches_constellation_default() -> None:
    # O a de projeto da fábrica é o da teoria só J2 para 91/6 (Sem4 §7.5)
    a, _ = repeat_orbit(91, 6)
    assert a == pytest.approx(constellation_module.DESIGN_SEMI_MAJOR_AXIS_M, abs=1.0)


def test_repeat_orbit_eccentricity_effect_is_tiny() -> None:
    # Sem4 §7.3: com a e congelada, Δa ≈ 1,4 cm e Δi ≈ 0,06″
    a0, i0 = repeat_orbit(91, 6)
    a1, i1 = repeat_orbit(91, 6, e=1.074e-3)
    assert abs(a1 - a0) < 0.05
    assert abs(math.degrees(i1 - i0)) * 3600.0 < 0.2


@pytest.mark.parametrize(("n", "d"), [(0, 6), (91, 0), (1, 1)])
def test_repeat_orbit_rejects_invalid(n: int, d: int) -> None:
    with pytest.raises(ValueError):
        repeat_orbit(n, d)


# --- constellation_phases: Sem4 eq. 42 e Tabela 4 -------------------------------------------------


def test_constellation_phases_table4() -> None:
    assert np.degrees(constellation_phases(6)) == pytest.approx(A1.psi0_deg, abs=1e-12)


def test_constellation_phases_general() -> None:
    assert constellation_phases(1) == pytest.approx([0.0])
    psi = constellation_phases(7)
    assert psi.shape == (7,)
    assert np.all((psi >= 0.0) & (psi < 2.0 * np.pi))
    # S_{j+1} fica 360°/P atrás de S_j
    assert np.mod(psi[:-1] - psi[1:], 2.0 * np.pi) == pytest.approx(np.full(6, 2 * np.pi / 7))
    with pytest.raises(ValueError):
        constellation_phases(0)


# --- Satellite ----------------------------------------------------------------------------------


def test_satellite_delegates_to_propagator() -> None:
    prop = KeplerPropagator(LEO, EPOCH)
    sat = Satellite("x1", "X-1", prop)
    assert sat.color == DEFAULT_COLOR == "#00e5ff"
    assert sat.cd is None and sat.area_over_mass is None
    assert sat.epoch == EPOCH
    t = np.array([0.0, 60.0, 3600.0])
    for got, want in zip(sat.propagate(t), prop.propagate(t), strict=True):
        assert np.array_equal(got, want)
    later = EPOCH + timedelta(hours=1)
    for got, want in zip(sat.propagate_utc(later, t), prop.propagate(t + 3600.0), strict=True):
        assert np.allclose(got, want, rtol=0.0, atol=1e-6)


def test_satellite_is_frozen() -> None:
    sat = Satellite("x1", "X-1", KeplerPropagator(LEO, EPOCH), cd=2.2, area_over_mass=0.005)
    with pytest.raises(dataclasses.FrozenInstanceError):
        sat.name = "outro"  # type: ignore[misc]


# --- Constellation ------------------------------------------------------------------------------


def _kepler_constellation(epoch_offsets_s: tuple[float, ...]) -> Constellation:
    sats = [
        Satellite(f"k{j}", f"K{j}", KeplerPropagator(LEO, EPOCH + timedelta(seconds=dt)))
        for j, dt in enumerate(epoch_offsets_s)
    ]
    return Constellation(EPOCH, sats)


def test_constellation_shapes_and_container() -> None:
    c = _kepler_constellation((0.0, 0.0, 0.0))
    r, v = c.propagate(np.arange(0.0, 600.0, 60.0))
    assert r.shape == v.shape == (3, 10, 3)
    r1, v1 = c.propagate(120.0)
    assert r1.shape == v1.shape == (3, 1, 3)
    assert len(c) == 3
    assert c.ids == ["k0", "k1", "k2"]
    assert [s.id for s in c] == c.ids
    assert c[1].id == "k1"
    assert c.epoch == EPOCH


def test_constellation_samples_each_satellite_at_its_own_epoch() -> None:
    # k1 tem época 1 h depois da constelação; k2, 30 min antes (ADR: t_i = t + (epoch − sat.epoch))
    c = _kepler_constellation((0.0, 3600.0, -1800.0))
    t = np.array([3600.0, 7200.0])
    r, v = c.propagate(t)
    for j, offset in enumerate((0.0, -3600.0, 1800.0)):
        r_j, v_j = c[j].propagator.propagate(t + offset)
        assert np.allclose(r[j], r_j, rtol=0.0, atol=1e-6)
        assert np.allclose(v[j], v_j, rtol=0.0, atol=1e-9)
    # propagate_utc com outra referência dá os mesmos instantes absolutos
    r_utc, _ = c.propagate_utc(EPOCH + timedelta(hours=1), t - 3600.0)
    assert np.allclose(r_utc, r, rtol=0.0, atol=1e-6)


def test_constellation_rejects_invalid() -> None:
    sat = Satellite("a", "A", KeplerPropagator(LEO, EPOCH))
    with pytest.raises(ValueError, match="tz-aware"):
        Constellation(datetime(2026, 10, 23), [sat])
    with pytest.raises(ValueError, match="ao menos um"):
        Constellation(EPOCH, [])
    with pytest.raises(ValueError, match="repetidos"):
        Constellation(EPOCH, [sat, sat])


# --- sso_constellation: a constelação de referência (A.1) ---------------------------------------


def test_sso_constellation_identity(constellation) -> None:
    assert len(constellation) == 6
    assert constellation.ids == ["s1", "s2", "s3", "s4", "s5", "s6"]
    assert [s.name for s in constellation] == ["S1", "S2", "S3", "S4", "S5", "S6"]
    colors = [s.color for s in constellation]
    assert colors == list(CONSTELLATION_PALETTE) and len(set(colors)) == 6
    assert constellation.epoch == EPOCH
    for sat in constellation:
        assert isinstance(sat.propagator, CowellPropagator)
        assert sat.epoch == EPOCH
        assert len(sat.propagator.models) == 5  # propagação com todas as forças
        assert sat.cd * sat.area_over_mass == pytest.approx(A1.drag_cd_a_m)


def test_sso_constellation_initial_osculating_plane(constellation) -> None:
    # A conversão preserva Ω e u da época; todos no mesmo plano, Ω do LTAN 22h30 (A.1)
    r, v = constellation.propagate(0.0)
    el = state_to_elements(r[:, 0], v[:, 0])
    assert np.degrees(el.raan) == pytest.approx(np.full(6, A1.raan_deg), abs=5e-4)
    u = np.degrees(np.mod(el.argp + _true_anomaly(el), 2.0 * np.pi))
    # u0 = ω + ν ≈ 90° + ψ0 (a menos de ~2e, órbita quase circular)
    du = np.mod(u - (90.0 + np.array(A1.psi0_deg)) + 180.0, 360.0) - 180.0
    assert np.abs(du).max() < 0.5


def _true_anomaly(el: KeplerianElements) -> np.ndarray:
    return np.asarray(mean_to_true(el.psi, el.e))


def test_sso_constellation_same_mean_elements(constellation) -> None:
    """a, e, i, ω e Ω médios (janela nodal a partir de t = 0) iguais para os seis.

    Medido (S5-01): a 6 888 087,6 … 6 888 089,2 m; e 1,0769e-3 … 1,0774e-3; i 97,463252° nos
    seis; ω 89,991 … 90,008°; Ω médio na janela 5,0515 … 5,0605°. A dispersão de Ω (±0,0045°) é o
    termo de curto período de J2 em Ω, que a conversão S4-04 não corrige (ela ajusta só a, o vetor
    e e i, e mantém Ω e u da época): amplitude ≈ ¾·J2·(R⊕/a)²·|cos i| ≈ 0,005°.
    """
    e_target = frozen_eccentricity(6888.089e3, I_A1)
    period = float(nodal_period(6888.089e3, e_target, I_A1))
    means = [window_mean_elements(s.propagator, 0.0, period) for s in constellation]
    a = np.array([m.a for m in means])
    e = np.array([m.e for m in means])
    i = np.degrees([m.i for m in means])
    argp = np.degrees([m.argp for m in means])
    # alvo: ψ0 diferente, mesmos a, e, i, ω. Tolerâncias: a conversão (S4-04) mais o arrasto, o
    # Sol e a Lua, que só entram na propagação (≈ 1 m em a numa volta).
    assert a == pytest.approx(np.full(6, A1.a_design_m), abs=3.0)
    assert e == pytest.approx(np.full(6, A1.e_frozen), abs=2e-6)
    assert i == pytest.approx(np.full(6, A1.i_calibrated_deg), abs=2e-5)
    assert argp == pytest.approx(np.full(6, A1.argp_deg), abs=0.05)
    # Ω médio: a média da janela [0, Tnod] está meio período à frente da época (Ω̇·Tnod/2).
    tau = np.linspace(0.0, period, 1201)
    raan_mean = np.array(
        [np.mean(np.unwrap(osculating_elements(s.propagator, tau).raan)) for s in constellation]
    )
    raan_dot = float(zonal_mean_rates(6888.089e3, e_target, I_A1, math.pi / 2).raan_dot)
    raan_mean_at_epoch = np.degrees(raan_mean - raan_dot * period / 2.0)
    assert raan_mean_at_epoch == pytest.approx(np.full(6, A1.raan_deg), abs=0.01)
    assert np.ptp(raan_mean_at_epoch) < 0.01


def test_sso_constellation_default_eccentricity_is_frozen() -> None:
    # e não dada → frozen_eccentricity(a, i) (A.1: 1,07721e-3); ψ0 da Tabela 4 nos elementos
    seen: list[KeplerianElements] = []

    def fake_mean_to_osc(mean, epoch, models, **kwargs):
        seen.append(mean)
        return SimpleNamespace(elements=mean)

    mp = pytest.MonkeyPatch()
    mp.setattr(constellation_module, "mean_to_osculating", fake_mean_to_osc)
    try:
        sso_constellation(EPOCH, (ZonalAcceleration(),), i=I_A1)
    finally:
        mp.undo()
    assert [m.e for m in seen] == pytest.approx([A1.e_frozen] * 6, abs=5e-10)
    assert np.degrees([m.psi for m in seen]) == pytest.approx(A1.psi0_deg, abs=1e-9)
    assert np.degrees([m.argp for m in seen]) == pytest.approx([90.0] * 6)
    raan = float(raan_for_ltan(A1.ltan_h, datetime_to_jd(EPOCH)))
    assert [m.raan for m in seen] == pytest.approx([raan] * 6, abs=1e-15)


def test_sso_constellation_calibrates_once_and_converts_with_zonal_only(
    monkeypatch, a1_models
) -> None:
    """Receita da Semana 4: calibração uma vez para os seis; conversão só com gravidade zonal."""
    calls: dict[str, list] = {"calib": [], "conv": []}

    def fake_calibrated(a, e, epoch, models, **kwargs):
        calls["calib"].append(tuple(models))
        return I_A1

    def fake_mean_to_osc(mean, epoch, models, **kwargs):
        calls["conv"].append(tuple(models))
        return SimpleNamespace(elements=mean)

    monkeypatch.setattr(constellation_module, "sso_inclination_calibrated", fake_calibrated)
    monkeypatch.setattr(constellation_module, "mean_to_osculating", fake_mean_to_osc)
    models = a1_models
    c = sso_constellation(EPOCH, models)
    assert len(c) == 6
    assert len(calls["calib"]) == 1
    assert len(calls["conv"]) == 6
    for used in calls["calib"] + calls["conv"]:
        assert used == (models[0],)  # só a ZonalAcceleration
    for sat in c:
        assert sat.propagator.models == models  # a propagação usa todas


def test_sso_constellation_options() -> None:
    # P = 3, paleta curta (repete em ciclo), sem arrasto (cd e A/m ficam None)
    c = sso_constellation(
        EPOCH, (J2Acceleration(),), n_sats=3, i=I_A1, e=1e-3, color_palette=("#111111", "#222222")
    )
    assert c.ids == ["s1", "s2", "s3"]
    assert [s.color for s in c] == ["#111111", "#222222", "#111111"]
    assert all(s.cd is None and s.area_over_mass is None for s in c)
    with pytest.raises(ValueError):
        sso_constellation(EPOCH, (J2Acceleration(),), n_sats=0, i=I_A1)
    with pytest.raises(ValueError):
        sso_constellation(EPOCH, (J2Acceleration(),), i=I_A1, color_palette=())


def test_sso_constellation_build_cost(built) -> None:
    # A.1: 6 conversões com i dada em 0,33 s; critério do card: < 3 s
    _, build_s = built
    assert build_s < A1.build_time_limit_s


def test_sso_constellation_propagate_shape(constellation) -> None:
    r, v = constellation.propagate(np.arange(0.0, 601.0, 60.0))
    assert r.shape == v.shape == (6, 11, 3)
    assert np.all(np.linalg.norm(r, axis=-1) > R_EARTH)


# --- Atraso no nó (Sem4 Tabela 4; A.1) -----------------------------------------------------------


def test_node_delay_between_neighbours(first_nodes) -> None:
    """No Cowell completo, o 1º nó de S_{j+1} depois do 1º nó de S_j vem Tnod/6 ≈ 949,5 s depois."""
    for j in range(5):
        t_j = first_nodes[j][0]
        t_next = first_nodes[j + 1][first_nodes[j + 1] > t_j][0]
        assert t_next - t_j == pytest.approx(A1.node_delay_s[1], abs=A1.node_delay_tol_s)


def test_node_delay_cumulative_table4(first_nodes) -> None:
    # Atraso de S_j em relação a S1, módulo Tnod (Tabela 4: 0; 949,5; …; 4 747,3 s)
    tnod = first_nodes[0][1] - first_nodes[0][0]
    t1 = first_nodes[0][0]
    delays = [np.mod(nodes[0] - t1, tnod) for nodes in first_nodes]
    assert delays == pytest.approx(A1.node_delay_s, abs=A1.node_delay_tol_s)


def test_first_nodes_regression_a1(first_nodes) -> None:
    # Regressão contra o experimento do A.1 (2 casas): satsim 4 274,463 … 3 324,934 s
    got = [nodes[0] for nodes in first_nodes]
    assert got == pytest.approx(A1.first_node_s, abs=0.02)
    assert first_nodes[1][0] - first_nodes[0][0] == pytest.approx(A1.delay_s2_s1_s, abs=0.02)
