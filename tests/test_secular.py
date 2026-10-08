"""Testes das taxas seculares de J2 e do período nodal (S2-01).

Referências: tests/golden/reference_values.py (J2_SECULAR_CASES_S2 e J2_SECULAR_MISC_S2, do
Apêndice A.1 do prompt da Semana 2; DESIGN_ORBIT, Sem4 §7.5), Sem4 eqs. 23–25 e 34.
"""

import math

import numpy as np
import pytest
from golden.reference_values import DESIGN_ORBIT, J2_SECULAR_CASES_S2, J2_SECULAR_MISC_S2

from satsim.constants import MU, R_EARTH
from satsim.secular import SecularRates, nodal_period, secular_rates

DEG_DAY = math.degrees(1.0) * 86400.0  # rad/s → °/dia
I_CRIT = math.acos(1.0 / math.sqrt(5.0))  # inclinação crítica, 63,4349°

# Tolerâncias [°/dia] (raan_dot, argp_dot) de cada caso, do Apêndice A.1 do prompt.
TOLERANCES = {
    "ISS": (0.002, 0.005),
    "SSO de projeto": (1e-4, 0.002),
    "Polar exata": (None, 0.002),  # Ω̇ testado em rad/s abaixo
    "Molniya": (0.002, None),  # ω̇ testado em rad/s abaixo
}

CASES = {c.orbit: c for c in J2_SECULAR_CASES_S2}


def _rates(orbit: str) -> SecularRates:
    c = CASES[orbit]
    return secular_rates(c.a_km * 1e3, c.e, math.radians(c.i_deg))


@pytest.mark.parametrize("case", J2_SECULAR_CASES_S2, ids=lambda c: c.orbit)
def test_secular_rates_appendix_a1(case) -> None:
    rates = _rates(case.orbit)
    tol_raan, tol_argp = TOLERANCES[case.orbit]
    if tol_raan is not None:
        assert rates.raan_dot * DEG_DAY == pytest.approx(case.raan_dot_deg_day, abs=tol_raan)
    if tol_argp is not None:
        assert rates.argp_dot * DEG_DAY == pytest.approx(case.argp_dot_deg_day, abs=tol_argp)


def test_iss_plane_precession_period() -> None:
    # Apêndice A.1: 360°/|Ω̇| ≈ 72,8 dias, entre 72 e 73 dias
    days = 360.0 / abs(_rates("ISS").raan_dot * DEG_DAY)
    assert 72.0 < days < 73.0
    assert days == pytest.approx(J2_SECULAR_MISC_S2.iss_plane_precession_days, abs=0.05)


def test_sso_nodal_and_kepler_period() -> None:
    # DESIGN_ORBIT (Sem4 §7.5): Tnod = 5696,70 s, Tkep = 5689,30 s
    c = CASES["SSO de projeto"]
    a, e, i = c.a_km * 1e3, c.e, math.radians(c.i_deg)
    tnod = nodal_period(a, e, i)
    tkep = 2.0 * math.pi / secular_rates(a, e, i).n
    assert tnod == pytest.approx(DESIGN_ORBIT.tnod_s, abs=0.05)
    assert tkep == pytest.approx(DESIGN_ORBIT.tkep_s, abs=0.05)
    assert tnod - tkep == pytest.approx(DESIGN_ORBIT.tnod_minus_tkep_s, abs=0.05)


def test_sso_psi_dot_correction() -> None:
    # Apêndice A.1: ψ̇/n − 1 = −6,6e-4 (o "+0,0065%" de Sem4 §3.5 é erro de digitação)
    rates = _rates("SSO de projeto")
    assert rates.psi_dot / rates.n - 1.0 == pytest.approx(
        J2_SECULAR_MISC_S2.sso_psi_dot_over_n_minus_1, abs=2e-5
    )


def test_polar_raan_dot_is_zero() -> None:
    # Ω̇ = −k·n·cos 90° = 0 (Sem4 eq. 23)
    assert _rates("Polar exata").raan_dot == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("a_km", [6800.0, 7200.0, 26561.75])
@pytest.mark.parametrize("e", [0.0, 0.74])
def test_critical_inclination_freezes_argp(a_km: float, e: float) -> None:
    # 5cos²i − 1 = 0 em i = arccos(1/√5), para qualquer a e e (Sem4 eq. 24)
    assert abs(secular_rates(a_km * 1e3, e, I_CRIT).argp_dot) < 1e-12


def test_molniya_argp_dot_is_zero() -> None:
    # Molniya a 63,4349° (A.1): |ω̇| < 1e-12 rad/s
    assert abs(_rates("Molniya").argp_dot) < 1e-12


def test_iss_nodal_period() -> None:
    # Apêndice A.1: Tnod ≈ 5573,9 s (≈ 15,50 rev/dia)
    c = CASES["ISS"]
    tnod = nodal_period(c.a_km * 1e3, c.e, math.radians(c.i_deg))
    assert tnod == pytest.approx(J2_SECULAR_MISC_S2.iss_tnod_s, abs=0.1)
    assert 86400.0 / tnod == pytest.approx(J2_SECULAR_MISC_S2.iss_rev_per_day, abs=0.005)


def test_zero_j2_gives_kepler() -> None:
    rates = secular_rates(7000e3, 0.01, 1.0, j2=0.0)
    assert rates.raan_dot == 0.0 and rates.argp_dot == 0.0
    assert rates.psi_dot == rates.n == pytest.approx(math.sqrt(MU / 7000e3**3), rel=1e-15)


def test_scalar_inputs_return_floats() -> None:
    rates = secular_rates(7000e3, 0.0, 1.0)
    assert all(isinstance(x, float) for x in (rates.n, rates.raan_dot, rates.argp_dot))
    assert isinstance(rates.psi_dot, float)
    assert isinstance(nodal_period(7000e3, 0.0, 1.0), float)


def test_vectorized_over_altitudes_matches_scalar() -> None:
    h = np.linspace(300e3, 1200e3, 19)
    a = R_EARTH + h
    i = math.radians(97.4396)
    rates = secular_rates(a, 1e-3, i)
    tnod = nodal_period(a, 1e-3, i)
    assert rates.raan_dot.shape == (19,) and tnod.shape == (19,)
    for k, ak in enumerate(a):
        single = secular_rates(float(ak), 1e-3, i)
        assert rates.n[k] == single.n
        assert rates.raan_dot[k] == single.raan_dot
        assert rates.argp_dot[k] == single.argp_dot
        assert rates.psi_dot[k] == single.psi_dot
        assert tnod[k] == nodal_period(float(ak), 1e-3, i)


@pytest.mark.parametrize(("a", "e"), [(-7000e3, 0.0), (0.0, 0.0), (7000e3, -0.1), (7000e3, 1.0)])
def test_rejects_invalid_elements(a: float, e: float) -> None:
    with pytest.raises(ValueError):
        secular_rates(a, e, 1.0)
