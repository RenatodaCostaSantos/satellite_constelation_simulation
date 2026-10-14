"""Testes do propagador numérico (Cowell) com J2 e das ferramentas de análise (S2-02).

Referências: tests/golden/reference_values.py (COWELL_REFERENCE_S2 e afins, do Apêndice A.4 do
prompt da Semana 2; DESIGN_ORBIT, Sem4 §7.5), Sem4 eqs. 23–25 e 34.
Os critérios (faixas e tolerâncias) são os do card S2-02; os valores do A.4 são referências de
comportamento citadas nos comentários.
"""

import math
import time
from dataclasses import replace
from datetime import UTC, datetime

import numpy as np
import pytest
from golden.reference_values import (
    COWELL_ARGP_DOT_FIT_E002_DEG_DAY,
    COWELL_PSI0_DEG,
    COWELL_RAAN0_DEG,
    COWELL_REFERENCE_S2,
    CowellReferenceCase,
)

from satsim.analysis import fit_rate, mean_elements, osculating_elements
from satsim.constants import MU
from satsim.elements import KeplerianElements, elements_to_state
from satsim.forces import j2_potential
from satsim.propagators import CowellPropagator, KeplerPropagator, MeanJ2Propagator, Propagator
from satsim.secular import nodal_period, secular_rates

EPOCH = datetime(2026, 9, 29, 12, 0, 0, tzinfo=UTC)
DEG_DAY = math.degrees(1.0) * 86400.0  # rad/s → °/dia
T_10D = np.arange(0.0, 10 * 86400.0 + 1.0, 60.0)  # 10 dias a 60 s (Apêndice A.4)
T_1D = np.arange(0.0, 86400.0 + 1.0, 60.0)

REF = {c.orbit: c for c in COWELL_REFERENCE_S2}


def _elements(case: CowellReferenceCase) -> KeplerianElements:
    return KeplerianElements(
        a=case.a_km * 1e3,
        e=case.e,
        i=math.radians(case.i_deg),
        raan=math.radians(COWELL_RAAN0_DEG),
        argp=math.radians(case.argp0_deg),
        psi=math.radians(COWELL_PSI0_DEG),
    )


class Run:
    """Uma propagação de 10 dias e o que se mede nela (calculada uma vez por órbita)."""

    def __init__(self, case: CowellReferenceCase) -> None:
        self.case = case
        self.nominal = _elements(case)
        self.prop = CowellPropagator(self.nominal, EPOCH)
        self.osc = osculating_elements(self.prop, T_10D)
        self.period = nodal_period(self.nominal.a, self.nominal.e, self.nominal.i)
        self.a_mean, self.e_mean, self.i_mean = mean_elements(T_10D, self.osc, self.period)


@pytest.fixture(scope="module", params=[c.orbit for c in COWELL_REFERENCE_S2])
def run(request) -> Run:
    return Run(REF[request.param])


@pytest.fixture(scope="module")
def sso() -> Run:
    return Run(REF["SSO"])


# --- Interface ---------------------------------------------------------------------------------


def test_is_propagator(sso: Run) -> None:
    assert isinstance(sso.prop, Propagator)
    assert sso.prop.epoch == EPOCH


def test_rejects_naive_epoch() -> None:
    with pytest.raises(ValueError):
        CowellPropagator(_elements(REF["SSO"]), datetime(2026, 9, 29, 12, 0, 0))


def test_shapes_and_t0_is_initial_state() -> None:
    el = _elements(REF["SSO"])
    prop = CowellPropagator(el, EPOCH)
    r, v = prop.propagate(0.0)
    assert r.shape == (1, 3) and v.shape == (1, 3)
    r0, v0 = elements_to_state(el)
    np.testing.assert_array_equal(r, r0)
    np.testing.assert_array_equal(v, v0)
    r, v = prop.propagate(np.arange(0.0, 600.0, 60.0))
    assert r.shape == (10, 3) and v.shape == (10, 3)


def test_state_input_equals_elements_input() -> None:
    el = _elements(REF["SSO"])
    r0, v0 = elements_to_state(el)
    t = np.array([0.0, 3000.0, 6000.0])
    r_el, v_el = CowellPropagator(el, EPOCH).propagate(t)
    r_st, v_st = CowellPropagator((r0[0], v0[0]), EPOCH).propagate(t)
    np.testing.assert_array_equal(r_el, r_st)
    np.testing.assert_array_equal(v_el, v_st)


def test_rejects_bad_state_shape() -> None:
    with pytest.raises(ValueError):
        CowellPropagator((np.zeros((1, 3)), np.zeros((1, 3))), EPOCH)


def test_negative_t_raises() -> None:
    prop = CowellPropagator(_elements(REF["SSO"]), EPOCH)
    with pytest.raises(ValueError):
        prop.propagate(-1.0)
    with pytest.raises(ValueError):
        prop.propagate([0.0, 10.0, -5.0])


def test_cache_out_of_order_and_repeated_calls_are_consistent() -> None:
    el = _elements(REF["SSO"])
    t = np.array([50000.0, 100.0, 86400.0, 0.0, 30000.0])
    prop = CowellPropagator(el, EPOCH)
    r_first, v_first = prop.propagate(t)
    # repetida, fora de ordem e um a um: os mesmos valores
    r_again, v_again = prop.propagate(t[::-1])
    np.testing.assert_array_equal(r_again, r_first[::-1])
    np.testing.assert_array_equal(v_again, v_first[::-1])
    for k, tk in enumerate(t):
        rk, vk = prop.propagate(tk)
        np.testing.assert_array_equal(rk[0], r_first[k])
        np.testing.assert_array_equal(vk[0], v_first[k])
    # cache crescendo em trechos (t pedido aos poucos): mesma trajetória dentro da tolerância
    grown = CowellPropagator(el, EPOCH)
    for t_max in (100.0, 30000.0, 50000.0, 86400.0):
        grown.propagate(t_max)
    r_grown, _ = grown.propagate(t)
    np.testing.assert_allclose(r_grown, r_first, rtol=0.0, atol=0.1)


# --- Física ------------------------------------------------------------------------------------


def test_without_j2_matches_kepler() -> None:
    # Critério S2-02: models=() coincide com o KeplerPropagator em < 0,1 m após 1 dia
    el = _elements(REF["SSO"])
    t = np.arange(0.0, 86400.0 + 1.0, 600.0)
    r, _ = CowellPropagator(el, EPOCH, models=()).propagate(t)
    rk, _ = KeplerPropagator(el, EPOCH).propagate(t)
    assert np.max(np.linalg.norm(r - rk, axis=1)) < 0.1


def test_energy_and_angular_momentum_z_conserved(sso: Run) -> None:
    # Critério S2-02: ε = v²/2 − µ/r − R(r) e ℓz constantes a 1e-9 relativo em 10 dias (60 s).
    # Com R = função perturbadora de Sem4 §2.2 (a = ∇R), a energia potencial é −µ/r − R.
    r, v = sso.prop.propagate(T_10D)
    r_norm = np.linalg.norm(r, axis=1)
    energy = 0.5 * np.sum(v * v, axis=1) - MU / r_norm - j2_potential(r)
    ell_z = r[:, 0] * v[:, 1] - r[:, 1] * v[:, 0]
    assert np.max(np.abs(energy / energy[0] - 1.0)) < 1e-9
    assert np.max(np.abs(ell_z / ell_z[0] - 1.0)) < 1e-9


def test_nodal_rate_matches_theory_at_measured_mean_a(run: Run) -> None:
    # Critério S2-02: Ω̇ ajustado vs secular_rates nos elementos médios medidos, erro < 0,5%
    # (A.4: −0,09% SSO, +0,08% ISS-like, +0,12% e = 0,02)
    fitted = fit_rate(T_10D, run.osc.raan)
    theory = secular_rates(run.a_mean, run.e_mean, run.i_mean).raan_dot
    assert abs(fitted / theory - 1.0) < 5e-3
    # A.4 (referência de comportamento): Ω̇ ajustado
    assert fitted * DEG_DAY == pytest.approx(run.case.raan_dot_fit_deg_day, rel=0.3)


def test_mean_semi_major_axis_offset(run: Run) -> None:
    # A.4 (referência de comportamento, ±30%): a_médio − a_nominal = +9,5 / −6,0 / −2,5 km
    offset_km = (run.a_mean - run.nominal.a) / 1e3
    assert offset_km == pytest.approx(run.case.a_mean_minus_nominal_km, rel=0.3)


def test_argp_rate_on_e002_orbit() -> None:
    # Critério S2-02: ω̇ só na órbita e = 0,02, erro < 1% (A.4: +0,16%). Em e ≈ 1e-3 o ω
    # osculador é mal condicionado e o ajuste não tem sentido (não testado).
    r = Run(REF["e = 0,02"])
    fitted = fit_rate(T_10D, r.osc.argp)
    theory = secular_rates(r.a_mean, r.e_mean, r.i_mean).argp_dot
    assert abs(fitted / theory - 1.0) < 1e-2
    assert fitted * DEG_DAY == pytest.approx(COWELL_ARGP_DOT_FIT_E002_DEG_DAY, rel=0.3)


@pytest.mark.parametrize("orbit", ["SSO", "ISS-like"])
def test_short_period_oscillation_of_a(orbit: str) -> None:
    # Critério S2-02: pico a pico do semieixo osculador numa revolução entre 5 e 40 km
    # (A.4: ≈ 19 km SSO, ≈ 12 km ISS-like)
    case = REF[orbit]
    el = _elements(case)
    period = nodal_period(el.a, el.e, el.i)
    t = np.arange(0.0, period, 30.0)
    osc = osculating_elements(CowellPropagator(el, EPOCH), t)
    p2p_km = (np.max(osc.a) - np.min(osc.a)) / 1e3
    assert 5.0 < p2p_km < 40.0
    assert p2p_km == pytest.approx(case.a_osc_peak_to_peak_km, rel=0.3)


def test_cowell_vs_mean_j2_with_measured_mean_elements(run: Run) -> None:
    # Critério S2-02: MeanJ2 com a, e, i médios medidos (ângulos nominais) vs Cowell,
    # diferença máxima em 24 h < 20 km (A.4: 6,8 / 14,0 / 6,6 km)
    mean = replace(run.nominal, a=run.a_mean, e=run.e_mean, i=run.i_mean)
    r_mean, _ = MeanJ2Propagator(mean, EPOCH).propagate(T_1D)
    r_cowell, _ = run.prop.propagate(T_1D)
    diff_km = np.max(np.linalg.norm(r_mean - r_cowell, axis=1)) / 1e3
    assert diff_km < 20.0


def test_nominal_elements_as_osculating_drift_apart(sso: Run) -> None:
    # Armadilha osculador ≠ médio (Notas S2-02): MeanJ2 nominal vs Cowell nominal na SSO,
    # A.4: 79 / 363 / 688 / 1333 km em 1 / 6 / 12 / 24 h (±30%)
    hours = np.array([1.0, 6.0, 12.0, 24.0]) * 3600.0
    r_nom, _ = MeanJ2Propagator(sso.nominal, EPOCH).propagate(hours)
    r_cowell, _ = sso.prop.propagate(hours)
    err_km = np.linalg.norm(r_nom - r_cowell, axis=1) / 1e3
    np.testing.assert_allclose(err_km, sso.case.nominal_error_1_6_12_24h_km, rtol=0.3)


def test_one_day_at_10s_is_fast() -> None:
    # Critério S2-02: 1 dia a 10 s de passo em < 5 s (a integração conta, prop. novo)
    t = np.arange(0.0, 86400.0 + 1.0, 10.0)
    prop = CowellPropagator(_elements(REF["SSO"]), EPOCH)
    start = time.perf_counter()
    r, _ = prop.propagate(t)
    assert time.perf_counter() - start < 5.0
    assert r.shape == (8641, 3)


# --- Análise -----------------------------------------------------------------------------------


def test_osculating_elements_of_kepler_are_constant() -> None:
    el = _elements(REF["e = 0,02"])
    osc = osculating_elements(KeplerPropagator(el, EPOCH), T_1D[:50])
    assert osc.a.shape == (50,)
    np.testing.assert_allclose(osc.a, el.a, rtol=1e-12)
    np.testing.assert_allclose(osc.e, el.e, rtol=1e-9)
    np.testing.assert_allclose(osc.i, el.i, rtol=1e-12)


def test_mean_elements_removes_periodic_term() -> None:
    period = 5000.0
    t = np.arange(0.0, 3.4 * period, 10.0)  # 3 revoluções inteiras + sobra
    a = 7000e3 + 9e3 * np.sin(2 * np.pi * t / period) + 4e3 * np.cos(4 * np.pi * t / period)
    el = KeplerianElements(a, np.full_like(t, 1e-3), np.full_like(t, 1.0), 0.0, 0.0, 0.0)
    a_mean, e_mean, i_mean = mean_elements(t, el, period)
    assert a_mean == pytest.approx(7000e3, abs=1.0)
    assert e_mean == pytest.approx(1e-3, rel=1e-12)
    assert i_mean == pytest.approx(1.0, rel=1e-12)


def test_mean_elements_needs_one_revolution() -> None:
    t = np.arange(0.0, 100.0, 10.0)
    el = KeplerianElements(np.full_like(t, 7e6), 0.0, 0.0, 0.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        mean_elements(t, el, 1000.0)


def test_fit_rate_unwraps() -> None:
    t = np.linspace(0.0, 1e5, 500)
    rate = 3e-4  # rad/s: dá várias voltas
    angle = np.mod(0.5 + rate * t, 2 * np.pi)
    assert fit_rate(t, angle) == pytest.approx(rate, rel=1e-12)
