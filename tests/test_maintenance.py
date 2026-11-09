"""Manutenção da constelação: deriva do traço e do LTAN (card S5-05).

Referências: Apêndice A.7 do prompt da Semana 5 (TRACK_MAINTENANCE_A7, MAINTENANCE_A7), Sem4
§6.1 e §7.6 (MISC_SEM4) e a referência própria do teste Cowell × fórmula (MAINTENANCE_TEST_S5),
em tests/golden/reference_values.py.

Custo: o ȧ da constelação é medido em 1 dia sobre a constelação de referência do conftest
(Cowell completo, já integrado na sessão); o teste Cowell × fórmula integra 6 dias de um satélite
no Cowell só J2–J4 + arrasto (≈ 3,5 s).
"""

import json
import math
from datetime import datetime

import numpy as np
import pytest
from golden.reference_values import (
    CONSTELLATION_A1,
    MAINTENANCE_A7,
    MAINTENANCE_TEST_S5,
    MISC_SEM4,
    TRACK_MAINTENANCE_A7,
)

from satsim.constants import R_EARTH, SECONDS_PER_DAY
from satsim.constellation import sso_constellation
from satsim.forces import DragAcceleration, ZonalAcceleration
from satsim.groundtrack import node_crossings, repeat_closure
from satsim.maintenance import (
    SUN_INCLINATION_DRIFT_DEG_DAY,
    YEAR_S,
    MaintenanceStatus,
    circular_velocity,
    delta_v_for_semi_major_axis,
    inclination_delta_v_per_year,
    ltan_drift_from_inclination,
    ltan_drift_from_semi_major_axis,
    ltan_offset,
    maintenance_status,
    measure_decay_rate,
    time_to_deadband,
    time_to_ltan_band,
    track_delta_v_per_year,
    track_drift,
)

M7, REF = MAINTENANCE_A7, MAINTENANCE_TEST_S5
DAY = SECONDS_PER_DAY
A = M7.a_m
INC = math.radians(M7.i_deg)
DI_DT_SUN = math.radians(M7.di_dt_sun_deg_day) / DAY
EPOCH = datetime.fromisoformat(CONSTELLATION_A1.epoch_utc)


# --- Sem4 §6.1 e §7.6 -----------------------------------------------------------------------------


def test_sem4_track_cadence() -> None:
    # ȧ = −35 m/dia: ±5 km em 5,72 dias (A.7); ±1 / ±5 / ±20 km → 2,6 / 5,7 / 11,4 dias (Sem4 §7.6)
    a_dot = MISC_SEM4.maneuver_adot_m_day / DAY
    assert time_to_deadband(5.0, a_dot, A) / DAY == pytest.approx(5.72, abs=0.005)
    for band_km, days in MISC_SEM4.maneuver_cadence_days:
        assert time_to_deadband(band_km, a_dot, A) / DAY == pytest.approx(days, abs=0.05)


def test_sem4_ltan_sensitivities() -> None:
    # Sem4 eq. 32: δh = +1 km → −0,73 min/ano; δi = +0,01° → +1,92 min/ano (Sem4: +1,9)
    assert ltan_drift_from_semi_major_axis(1000.0, A) == pytest.approx(
        M7.ltan_per_dh_1km_min_year, abs=0.005
    )
    assert ltan_drift_from_inclination(math.radians(0.01), INC) == pytest.approx(
        M7.ltan_per_di_001deg_min_year, abs=0.005
    )
    assert ltan_drift_from_inclination(math.radians(0.01), INC) == pytest.approx(
        MISC_SEM4.ltan_drift_min_per_year_per_di_001deg, abs=0.05
    )
    arr = ltan_drift_from_inclination(np.radians([0.01, -0.01]), INC)
    assert arr.shape == (2,) and arr[0] == pytest.approx(-arr[1])


def test_track_drift_formula_and_inverse() -> None:
    a_dot = -23.71 / DAY
    lam, km = track_drift(6 * DAY, a_dot, A)
    assert km == pytest.approx(lam * R_EARTH / 1e3)
    assert km == pytest.approx(M7.closure_formula_6d_km, abs=0.005)  # A.7: 3,72 km
    # λ ∝ t²: o dobro do tempo, o quádruplo da deriva; inverso de time_to_deadband
    assert track_drift(12 * DAY, a_dot, A)[1] == pytest.approx(4 * km)
    t5 = time_to_deadband(5.0, a_dot, A)
    assert track_drift(t5, a_dot, A)[1] == pytest.approx(5.0)
    _, arr = track_drift(np.array([0.0, DAY, 2 * DAY]), a_dot, A)
    assert arr.shape == (3,) and arr[0] == 0.0
    assert time_to_deadband(5.0, 0.0, A) == math.inf


def test_circular_velocity() -> None:
    assert circular_velocity(A) == pytest.approx(M7.v_m_s, abs=0.05)


# --- A.7: tabela do traço e LTAN ------------------------------------------------------------------


@pytest.mark.parametrize("row", TRACK_MAINTENANCE_A7, ids=lambda r: f"{r.a_dot_m_day:g}")
def test_track_maintenance_table_a7(row) -> None:
    status = maintenance_status("x", row.a_dot_m_day / DAY, A, INC)
    rel = M7.rel_tol
    assert status.days_to_track_correction == pytest.approx(row.days_to_5km, rel=rel)
    assert status.track_delta_a_m == pytest.approx(row.delta_a_m, rel=rel)
    assert status.track_delta_v_per_maneuver_m_s == pytest.approx(
        row.delta_v_per_maneuver_m_s, rel=rel
    )
    assert status.track_maneuvers_per_year == pytest.approx(row.maneuvers_per_year, rel=rel)
    assert status.track_delta_v_per_year_m_s == pytest.approx(row.delta_v_per_year_m_s, rel=rel)
    # coerência entre as funções: Δv/ano = Δv por manobra × manobras/ano
    assert status.track_delta_v_per_year_m_s == pytest.approx(
        status.track_delta_v_per_maneuver_m_s * status.track_maneuvers_per_year
    )
    assert delta_v_for_semi_major_axis(row.delta_a_m, A) == pytest.approx(
        row.delta_v_per_maneuver_m_s, rel=rel
    )
    assert track_delta_v_per_year(row.a_dot_m_day / DAY, A) == pytest.approx(
        row.delta_v_per_year_m_s, rel=rel
    )


def test_ltan_drift_from_sun_a7() -> None:
    """di/dt = −9,8e-5 °/dia (S4-02): LTAN −6,87 min/ano² (−3,43·t²), ±5 min em 1,21 ano e
    4,75 m/s/ano de inclinação."""
    assert SUN_INCLINATION_DRIFT_DEG_DAY == M7.di_dt_sun_deg_day
    assert math.degrees(DI_DT_SUN) * YEAR_S == pytest.approx(M7.di_dt_sun_deg_year, rel=M7.rel_tol)
    one_year = ltan_offset(YEAR_S, DI_DT_SUN, INC)
    assert one_year == pytest.approx(M7.ltan_quadratic_coeff_min, rel=M7.rel_tol)
    assert 2 * one_year == pytest.approx(M7.ltan_accel_min_year2, rel=M7.rel_tol)
    assert ltan_offset(2 * YEAR_S, DI_DT_SUN, INC) == pytest.approx(4 * one_year)  # quadrático
    years = time_to_ltan_band(M7.ltan_band_min, DI_DT_SUN, INC) / YEAR_S
    assert years == pytest.approx(M7.years_to_ltan_band, rel=M7.rel_tol)
    assert abs(ltan_offset(years * YEAR_S, DI_DT_SUN, INC)) == pytest.approx(M7.ltan_band_min)
    assert inclination_delta_v_per_year(DI_DT_SUN, A) == pytest.approx(
        M7.inclination_delta_v_per_year_m_s, rel=M7.rel_tol
    )
    assert time_to_ltan_band(5.0, 0.0, INC) == math.inf
    arr = ltan_offset(np.array([0.0, YEAR_S]), DI_DT_SUN, INC)
    assert arr.shape == (2,) and arr[0] == 0.0


def test_track_and_ltan_budgets_are_similar_but_track_dominates_cadence() -> None:
    # Sem4 §6.2: os dois Δv/ano são parecidos (≈ 4,8 m/s), mas o traço pede uma manobra por
    # semana e o LTAN uma por ano
    status = maintenance_status("s1", TRACK_MAINTENANCE_A7[0].a_dot_m_day / DAY, A, INC)
    assert status.track_delta_v_per_year_m_s == pytest.approx(
        status.inclination_delta_v_per_year_m_s, rel=0.02
    )
    assert status.days_to_track_correction < 7.0 < 365.0 < status.days_to_ltan_correction


# --- Cowell: ȧ medido e a fórmula do traço --------------------------------------------------------


def test_measured_decay_rate_reference_constellation(reference_constellation) -> None:
    """ȧ ≈ −23,7 m/dia (A.7, rel 2%) nos seis satélites (medido −23,95 m/dia em 1 dia)."""
    rates = [measure_decay_rate(s, 1.0) * DAY for s in reference_constellation]
    ref = TRACK_MAINTENANCE_A7[0].a_dot_m_day
    for rate in rates:
        assert rate == pytest.approx(ref, rel=M7.rel_tol)
    assert np.ptp(rates) < 0.01 * abs(ref)  # mesmo a, mesmo arrasto: mesmo ȧ
    with pytest.raises(ValueError, match="2 períodos"):
        measure_decay_rate(reference_constellation[0], 0.05)


@pytest.fixture(scope="module")
def drag_s1():
    """S1 com o a calibrado (S5-02), Cowell só J2–J4 + arrasto (C_D·A/m = 0,011)."""
    models = (ZonalAcceleration(), DragAcceleration(cd=2.2, area_over_mass=0.005))
    return sso_constellation(EPOCH, models, n_sats=1, a=REF.a_calibrated_m, i=INC)[0]


def test_track_drift_formula_vs_cowell(drag_s1) -> None:
    """Fechamento de 6 dias com arrasto e o a calibrado: a fórmula com o ȧ medido fica a ≤ 10%
    do Cowell (A.4/A.7: 3,68 km medidos, 3,72 previstos; aqui, J2–J4 + arrasto: 3,84 e 3,77)."""
    _, lon = node_crossings(drag_s1, EPOCH, 6 * DAY + 2 * 5697.0)
    closure_km = repeat_closure(lon, 91) * R_EARTH / 1e3
    a_dot = measure_decay_rate(drag_s1, 6.0)
    _, formula_km = track_drift(6 * DAY, a_dot, REF.a_calibrated_m)
    assert closure_km > 0.0  # o arrasto adianta o traço para leste
    assert abs(formula_km / closure_km - 1.0) <= M7.closure_formula_tol
    assert closure_km == pytest.approx(REF.closure_6d_km, abs=REF.regression_tol_km)
    assert formula_km == pytest.approx(REF.formula_6d_km, abs=REF.regression_tol_km)
    assert a_dot * DAY == pytest.approx(REF.a_dot_6d_m_day, abs=0.01)


# --- MaintenanceStatus ----------------------------------------------------------------------------


def test_maintenance_status_per_satellite(reference_constellation) -> None:
    statuses = [
        maintenance_status(sat.id, measure_decay_rate(sat, 1.0), A, INC)
        for sat in reference_constellation
    ]
    assert [s.sat_id for s in statuses] == reference_constellation.ids
    for st in statuses:
        assert isinstance(st, MaintenanceStatus)
        assert st.days_to_track_correction == pytest.approx(6.95, rel=M7.rel_tol)
        assert st.ltan_offset_1y_min == pytest.approx(M7.ltan_quadratic_coeff_min, rel=M7.rel_tol)
        assert st.di_dt_deg_day == pytest.approx(M7.di_dt_sun_deg_day)
        assert st.track_deadband_km == 5.0 and st.ltan_band_min == 5.0
    # pronto para a cena: dicionário serializável em JSON
    data = json.loads(json.dumps(statuses[0].as_dict()))
    assert data["sat_id"] == "s1" and "track_delta_v_per_year_m_s" in data
