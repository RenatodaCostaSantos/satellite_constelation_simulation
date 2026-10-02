"""Consistência interna dos valores de referência (S1-07).

Os testes recalculam números de uma tabela a partir de outros números dos PDFs; uma falha
indica erro de transcrição em tests/golden/reference_values.py (ou no próprio documento).
"""

import math
from datetime import UTC, datetime

import pytest
from golden.reference_values import (
    CONSTANTS,
    CONSTELLATION_ORBIT,
    CONSTELLATION_TABLE4,
    DESIGN_ORBIT,
    GRID_CHECK,
    J2_SECULAR_RATES,
    SSO_CIRCULAR_TABLE1,
    SSO_REPEAT_TABLE3,
)

from satsim import constants
from satsim.elements import KeplerianElements
from satsim.propagators import KeplerPropagator

D = DESIGN_ORBIT

# --- Com o código existente --------------------------------------------------------------------


def test_constants_match_golden() -> None:
    # Sem4 código §9.1 e tabela de símbolos
    assert constants.J2 == CONSTANTS.j2
    assert constants.J3 == CONSTANTS.j3
    assert constants.J4 == CONSTANTS.j4
    assert constants.MU == CONSTANTS.mu_m3_s2
    assert constants.R_EARTH == pytest.approx(CONSTANTS.r_earth_km * 1e3, rel=1e-15)
    assert constants.OMEGA_EARTH == CONSTANTS.omega_earth_rad_s
    assert constants.T_TROP_DAYS == CONSTANTS.t_trop_days


def test_kepler_period_matches_golden() -> None:
    # Sem4 §7.5: Tkep = 5689,30 s para a = 6888,089 km
    el = KeplerianElements(D.a_km * 1e3, D.e, math.radians(D.i_deg), 0.0, 0.0, 0.0)
    prop = KeplerPropagator(el, datetime(2026, 1, 1, tzinfo=UTC))
    assert prop.period == pytest.approx(D.tkep_s, abs=0.05)


# --- Tabela 1 (Sem4 §4.3) ----------------------------------------------------------------------


@pytest.mark.parametrize("row", SSO_CIRCULAR_TABLE1, ids=lambda r: f"h{r.h_km:g}")
def test_table1_rev_per_day(row) -> None:
    assert 86400.0 / row.tnod_s == pytest.approx(row.rev_per_day, abs=1e-3)


def test_table1_inclination_increases_with_altitude() -> None:
    i_values = [row.i_deg for row in SSO_CIRCULAR_TABLE1]
    assert i_values == sorted(i_values)


# --- Tabela 3 (Sem4 §7.4) ----------------------------------------------------------------------

REGULAR = [r for r in SSO_REPEAT_TABLE3 if r.gcd == 1]
DEGENERATE = [r for r in SSO_REPEAT_TABLE3 if r.gcd != 1]


@pytest.mark.parametrize("row", SSO_REPEAT_TABLE3, ids=lambda r: f"{r.n_revs}/{r.d_days}")
def test_table3_gcd_column(row) -> None:
    assert math.gcd(row.n_revs, row.d_days) == row.gcd


@pytest.mark.parametrize("row", REGULAR, ids=lambda r: f"{r.n_revs}/{r.d_days}")
def test_table3_grid(row) -> None:
    assert 360.0 / row.n_revs == pytest.approx(row.grid_deg, abs=1e-3)
    assert 2 * math.pi * CONSTANTS.r_earth_km / row.n_revs == pytest.approx(row.grid_km, abs=0.1)


def test_table3_degenerate_rows() -> None:
    # 93/6 e 90/6: gcd > 1, sem grade definida
    assert {(r.n_revs, r.d_days) for r in DEGENERATE} == {(93, 6), (90, 6)}
    assert all(r.gcd > 1 and r.grid_deg is None and r.grid_km is None for r in DEGENERATE)


# --- Órbita de projeto 91/6 (Sem4 §7.5) --------------------------------------------------------


def test_design_orbit_repeat_geometry() -> None:
    assert math.gcd(D.n_revs, D.d_days) == D.gcd == 1
    assert 86400.0 * D.d_days / D.n_revs == pytest.approx(D.tnod_s, abs=0.01)
    assert D.n_revs / D.d_days == pytest.approx(D.rev_per_day, abs=1e-3)
    assert -360.0 * D.d_days / D.n_revs == pytest.approx(D.dlambda_per_rev_deg, abs=1e-3)
    assert 360.0 / D.n_revs == pytest.approx(D.grid_deg, abs=1e-3)
    assert D.tnod_s / 60.0 == pytest.approx(D.tnod_min, abs=1e-3)


def test_design_orbit_altitude() -> None:
    assert D.a_km - CONSTANTS.r_earth_km == pytest.approx(D.h_km, abs=1e-3)
    # Sem4 §1 e §6.3: excursão de altitude 2ae = 14,8 km
    assert 2 * D.a_km * D.e == pytest.approx(D.h_excursion_km, abs=0.05)


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Discrepância no Sem4 §1: pelos elementos, a(1 − e) − R⊕ = 502,554 km e "
        "a(1 + e) − R⊕ = 517,350 km (arredondam para 502,6 e 517,3); o PDF traz 502,5 e 517,4, "
        "cuja diferença 14,9 km não bate com a excursão 2ae = 14,8 km do mesmo parágrafo."
    ),
)
def test_design_orbit_instantaneous_altitude_sem4_s1() -> None:
    # tolerância = meia unidade da última casa publicada (0,1 km)
    h_min = D.a_km * (1 - D.e) - CONSTANTS.r_earth_km
    h_max = D.a_km * (1 + D.e) - CONSTANTS.r_earth_km
    assert h_min == pytest.approx(D.h_inst_min_km, abs=0.05)
    assert h_max == pytest.approx(D.h_inst_max_km, abs=0.05)
    assert D.h_inst_max_km - D.h_inst_min_km == pytest.approx(D.h_excursion_km, abs=0.05)


# --- Consistência entre grupos -----------------------------------------------------------------


def test_cross_group_consistency() -> None:
    # Tnod − Tkep (Sem4 §7.1)
    assert D.tnod_s - D.tkep_s == pytest.approx(D.tnod_minus_tkep_s, abs=0.05)
    assert 100 * (D.tnod_s - D.tkep_s) / D.tkep_s == pytest.approx(D.tnod_minus_tkep_pct, abs=0.01)
    # Taxa heliossíncrona (Sem4 eq. 26) = Ω-ponto da órbita de projeto = constante do código
    omega_sun_deg_day = math.degrees(constants.OMEGA_SUN) * 86400.0
    assert omega_sun_deg_day == pytest.approx(CONSTANTS.omega_sun_deg_day, abs=1e-6)
    assert D.raan_dot_deg_day == CONSTANTS.omega_sun_deg_day
    assert math.radians(CONSTANTS.omega_sun_deg_day) / 86400 == pytest.approx(
        CONSTANTS.omega_sun_rad_s, rel=1e-6
    )
    # Latitude máxima = 180° − i (Sem4 eq. 12)
    assert 180.0 - D.i_deg == pytest.approx(D.lat_max_deg, abs=0.01)
    # Armadilha de normalização (Sem4 §2.1): C20 = N20·C̄20 e J2 = −C20
    assert CONSTANTS.n20 == pytest.approx(math.sqrt(5.0), abs=1e-4)
    assert CONSTANTS.n20 * CONSTANTS.cbar20 == pytest.approx(CONSTANTS.c20, abs=1e-7)
    assert -CONSTANTS.c20 == pytest.approx(CONSTANTS.j2, rel=1e-4)


def test_same_quantity_same_value_across_groups() -> None:
    """Um valor que aparece em várias tabelas deve ser o mesmo (a menos do arredondamento)."""
    t1 = next(r for r in SSO_CIRCULAR_TABLE1 if r.h_km == D.h_nominal_label_km)
    t3 = next(r for r in SSO_REPEAT_TABLE3 if (r.n_revs, r.d_days) == (D.n_revs, D.d_days))
    sso_rate = next(r for r in J2_SECULAR_RATES if r.orbit == "SSO de projeto")

    # h de projeto: 509,952 (A.2) ≈ 509,95 (Tabela 3); rótulo 510 (Tabela 1)
    assert t3.h_km == pytest.approx(D.h_km, abs=0.005)
    # i: 97,4396 (A.2) ≈ 97,440 (Tabelas 1 e 3) ≈ 97,44 (A.5)
    for i_deg in (t1.i_deg, t3.i_deg, sso_rate.i_deg):
        assert i_deg == pytest.approx(D.i_deg, abs=5e-4)
    # Tnod e rev/dia: A.2 vs Tabela 1
    assert t1.tnod_s == pytest.approx(D.tnod_s, abs=0.05)
    assert t1.rev_per_day == D.rev_per_day
    # grade: A.2 vs Tabela 3
    assert t3.grid_deg == D.grid_deg
    assert t3.grid_km == D.grid_km
    # taxas seculares: A.5 vs A.2
    assert sso_rate.raan_dot_deg_day == pytest.approx(D.raan_dot_deg_day, abs=1e-3)
    assert sso_rate.argp_dot_deg_day == D.argp_dot_j2_deg_day
    # constelação (Tabela 4) compartilha a, e, i com a órbita de projeto
    assert CONSTELLATION_ORBIT.a_km == D.a_km
    assert CONSTELLATION_ORBIT.e == D.e
    assert CONSTELLATION_ORBIT.i_deg == D.i_deg
    # verificação da grade (Sem4 §9.2) com P = 6 em 1 dia
    assert GRID_CHECK.nodes_p6_1day == D.n_revs
    assert GRID_CHECK.step_p6_deg == pytest.approx(D.grid_deg, abs=1e-3)


def test_constellation_phasing() -> None:
    # Sem4 Tabela 4: ψ0 espaçados de 60° e atraso no nó = k·Tnod/6
    n_sats = len(CONSTELLATION_TABLE4)
    assert n_sats == 6
    for k, sat in enumerate(CONSTELLATION_TABLE4):
        assert sat.psi0_deg == pytest.approx((-60.0 * k) % 360.0)
        assert sat.node_delay_s == pytest.approx(k * D.tnod_s / n_sats, abs=0.06)
