"""Testes do projeto de órbita heliossíncrona (S2-05).

Referências: tests/golden/reference_values.py (SSO_CIRCULAR_TABLE1, Sem4 §4.3, Tabela 1;
MISC_SEM4, Sem4 §4.3 e §6.1; DESIGN_ORBIT, Sem4 §7.5; SSO_DESIGN_S2, Apêndices A.2 e A.5 do
prompt da Semana 2), Sem4 eqs. 26, 27, 34 e 37.
"""

import math

import numpy as np
import pytest
from golden.reference_values import (
    DESIGN_ORBIT,
    MISC_SEM4,
    SSO_CIRCULAR_TABLE1,
    SSO_DESIGN_S2,
)

from satsim.constants import OMEGA_EARTH, OMEGA_SUN, R_EARTH
from satsim.orbit_design import (
    SsoRow,
    ltan_drift_min_per_year,
    sso_inclination,
    sso_max_semi_major_axis,
    sso_row,
)
from satsim.secular import secular_rates

# --- Tabela 1 ----------------------------------------------------------------------------------


@pytest.mark.parametrize("ref", SSO_CIRCULAR_TABLE1, ids=lambda r: f"h{r.h_km:g}")
def test_table1(ref) -> None:
    # SSO_CIRCULAR_TABLE1 (Sem4 Tabela 1): i (abs 1e-3°), Tnod (abs 0,1 s), rev/dia (abs 0,002)
    row = sso_row(ref.h_km * 1e3)
    assert isinstance(row, SsoRow)
    assert row.h_km == pytest.approx(ref.h_km, abs=1e-9)
    assert row.i_deg == pytest.approx(ref.i_deg, abs=1e-3)
    assert row.nodal_period_s == pytest.approx(ref.tnod_s, abs=0.1)
    assert row.rev_per_day == pytest.approx(ref.rev_per_day, abs=0.002)


def test_table1_rev_per_day_is_eq37() -> None:
    # Sem4 eq. 37: rev/dia = (ψ̇ + ω̇)/(ω⊕ − Ω̇) (voltas nodais por dia nodal)
    a = R_EARTH + 510e3
    i = sso_inclination(a)
    r = secular_rates(a, 0.0, i)
    expected = (r.psi_dot + r.argp_dot) / (OMEGA_EARTH - r.raan_dot)
    assert sso_row(510e3).rev_per_day == pytest.approx(expected, rel=1e-15)


# --- Teto de altitude --------------------------------------------------------------------------


def test_altitude_ceiling() -> None:
    # MISC_SEM4 (Sem4 §4.3): a_max = 12 352 km e h_max = 5 974 km (abs 1 km)
    a_max = sso_max_semi_major_axis()
    assert a_max / 1e3 == pytest.approx(MISC_SEM4.sso_a_max_km, abs=1.0)
    assert (a_max - R_EARTH) / 1e3 == pytest.approx(MISC_SEM4.sso_h_max_km, abs=1.0)
    # em a_max, cos i = −1 (i = 180°); acima dele não há SSO
    assert math.cos(sso_inclination(a_max)) == pytest.approx(-1.0, abs=1e-12)
    assert math.degrees(sso_inclination(a_max)) == pytest.approx(180.0, abs=1e-3)
    with pytest.raises(ValueError):
        sso_inclination(a_max * (1.0 + 1e-6))


@pytest.mark.parametrize("e", [-0.1, 1.0])
def test_rejects_invalid_eccentricity(e: float) -> None:
    with pytest.raises(ValueError):
        sso_inclination(R_EARTH + 500e3, e)


# --- Propriedades ------------------------------------------------------------------------------


def test_always_retrograde() -> None:
    # Sem4 eq. 27: cos i < 0 → toda SSO é retrógrada (i > 90°), h de 200 a 1500 km
    h = np.linspace(200e3, 1500e3, 131)
    i = sso_inclination(R_EARTH + h)
    assert i.shape == (131,)
    assert np.all(i > math.pi / 2)


def test_self_consistency() -> None:
    # Sem4 eq. 26: com i = sso_inclination(a), Ω̇ (eq. 23) = OMEGA_SUN (rel 1e-9)
    for h in (300e3, 510e3, 800e3, 1200e3, 3000e3):
        a = R_EARTH + h
        i = sso_inclination(a)
        assert secular_rates(a, 0.0, i).raan_dot == pytest.approx(OMEGA_SUN, rel=1e-9)
    # também com e ≠ 0 (órbita de projeto)
    a, e = DESIGN_ORBIT.a_km * 1e3, DESIGN_ORBIT.e
    assert secular_rates(a, e, sso_inclination(a, e)).raan_dot == pytest.approx(OMEGA_SUN, rel=1e-9)


def test_design_orbit_inclination() -> None:
    # DESIGN_ORBIT (Sem4 §7.5): a = 6888,089 km, e = 1,074e-3 → i = 97,4396°
    i = sso_inclination(DESIGN_ORBIT.a_km * 1e3, DESIGN_ORBIT.e)
    assert math.degrees(i) == pytest.approx(DESIGN_ORBIT.i_deg, abs=1e-4)


def test_eccentricity_effect() -> None:
    # SSO_DESIGN_S2 (A.2): cos i ∝ (1 − e²)²; e = 0,01 → variação relativa −2,0e-4 (abs 1e-6)
    a = R_EARTH + 510e3
    rel = math.cos(sso_inclination(a, 0.01)) / math.cos(sso_inclination(a)) - 1.0
    assert rel == pytest.approx(SSO_DESIGN_S2.cos_i_rel_change_e001, abs=1e-6)


def test_inclination_vs_altitude() -> None:
    # SSO_DESIGN_S2 (A.2): di/dh em 510 km ≈ 0,0038 °/km (abs 2e-4), por diferenças centrais
    a = R_EARTH + 510e3
    di_dh = (sso_inclination(a + 1e3) - sso_inclination(a - 1e3)) / 2.0
    assert math.degrees(di_dh) == pytest.approx(SSO_DESIGN_S2.di_dh_510_deg_per_km, abs=2e-4)
    # MISC_SEM4 (Sem4 §6.1): tan i na órbita de projeto ≈ −7,66 (abs 0,01)
    assert math.tan(math.radians(DESIGN_ORBIT.i_deg)) == pytest.approx(
        MISC_SEM4.tan_i_design, abs=0.01
    )


def test_vectorized_matches_scalar() -> None:
    a = R_EARTH + np.array([300e3, 510e3, 1200e3])
    i = sso_inclination(a, 1e-3)
    for k, ak in enumerate(a):
        assert i[k] == sso_inclination(float(ak), 1e-3)
    assert isinstance(sso_inclination(7000e3), float)


# --- Sensibilidade da LTAN (Sem4 §6.1) ---------------------------------------------------------


def test_ltan_drift_sensitivity() -> None:
    # MISC_SEM4 (Sem4 §6.1, eq. 32), na órbita de projeto, deriva de LTAN em min/ano:
    # δh = +1 km com i fixo → −0,73 (abs 0,02); δi = +0,01° com h fixo → +1,9 (abs 0,05)
    a, e, i = DESIGN_ORBIT.a_km * 1e3, DESIGN_ORBIT.e, math.radians(DESIGN_ORBIT.i_deg)
    base = ltan_drift_min_per_year(secular_rates(a, e, i).raan_dot)
    dh = ltan_drift_min_per_year(secular_rates(a + 1e3, e, i).raan_dot) - base
    di = ltan_drift_min_per_year(secular_rates(a, e, i + math.radians(0.01)).raan_dot) - base
    assert dh == pytest.approx(MISC_SEM4.ltan_drift_min_per_year_per_dh_1km, abs=0.02)
    assert di == pytest.approx(MISC_SEM4.ltan_drift_min_per_year_per_di_001deg, abs=0.05)
    # a órbita de projeto é heliossíncrona: deriva própria desprezível (< 0,01 min/ano)
    assert abs(base) < 0.01


def test_ltan_drift_zero_at_sun_rate() -> None:
    assert ltan_drift_min_per_year(OMEGA_SUN) == 0.0
