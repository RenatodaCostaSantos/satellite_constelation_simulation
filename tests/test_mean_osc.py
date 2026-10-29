"""Testes da conversão médio ↔ osculador (S4-04).

Referências: tests/golden/reference_values.py (MEAN_OSC_A6, Apêndice A.6 do prompt da Semana 4,
com a discrepância da janela documentada no golden; DESIGN_ORBIT, Sem4 §7.5).
"""

import math
from datetime import datetime

import numpy as np
import pytest
from golden.reference_values import DESIGN_ORBIT, MEAN_OSC_A6

from satsim.elements import KeplerianElements, mean_to_true
from satsim.forces import J2Acceleration, ZonalAcceleration
from satsim.mean_osc import (
    MeanElements,
    mean_to_osculating,
    osculating_to_mean,
    window_mean_elements,
)
from satsim.propagators import CowellPropagator, KeplerPropagator
from satsim.secular import nodal_period

A6 = MEAN_OSC_A6
EPOCH = datetime.fromisoformat(A6.epoch_utc)
MODELS = (ZonalAcceleration(),)
DESIGN_MEAN = KeplerianElements(
    a=DESIGN_ORBIT.a_km * 1e3,
    e=DESIGN_ORBIT.e,
    i=math.radians(DESIGN_ORBIT.i_deg),
    raan=math.radians(A6.raan0_deg),
    argp=math.radians(DESIGN_ORBIT.argp_deg),
    psi=math.radians(A6.psi0_deg),
)
PERIOD = float(nodal_period(DESIGN_MEAN.a, DESIGN_MEAN.e, DESIGN_MEAN.i))


@pytest.fixture(scope="module")
def converted():
    return mean_to_osculating(DESIGN_MEAN, EPOCH, MODELS)


def _window_series(initial: KeplerianElements, days: float) -> MeanElements:
    prop = CowellPropagator(initial, EPOCH, MODELS)
    n_windows = int(days * 86400.0 // PERIOD)
    return window_mean_elements(prop, np.arange(n_windows) * PERIOD, PERIOD)


def test_window_mean_of_kepler_orbit_is_exact() -> None:
    # Sem perturbações os elementos não oscilam: a média é o próprio elemento, também para t0 (K,)
    el = KeplerianElements(7.0e6, 0.01, 1.0, 0.5, 0.7, 0.2)
    prop = KeplerPropagator(el, EPOCH)
    m = window_mean_elements(prop, 0.0, 5000.0)
    assert m.a == pytest.approx(el.a, rel=1e-12)
    assert m.e == pytest.approx(el.e, rel=1e-9)
    assert m.i == pytest.approx(el.i, abs=1e-12)
    assert m.argp == pytest.approx(el.argp, abs=1e-9)
    series = window_mean_elements(prop, [0.0, 5000.0, 10000.0], 5000.0)
    assert np.shape(series.a) == (3,)
    np.testing.assert_allclose(series.a, el.a, rtol=1e-12)


def test_convergence_a6(converted) -> None:
    # A.6: resíduos decrescentes, Δa da 1ª iteração ≈ 9,5 km (abs 0,1 km), tolerâncias em ≤ 6
    # iterações (referência: Δa = 2e-6 m na 5ª; aqui convergiu na 4ª)
    res = converted.residuals
    assert converted.converged and converted.iterations <= 6
    assert res[0, 0] / 1e3 == pytest.approx(
        A6.residual_a_m[0] / 1e3, abs=A6.first_residual_a_tol_km
    )
    assert np.all(np.diff(res[:, 0]) < 0.0)  # |Δa| decresce a cada iteração
    assert np.all(res[-1] < [1e-3, 1e-9, 1e-9, 1e-9])
    # demais colunas do A.6, na ordem de grandeza (a janela difere: ver golden)
    assert res[0, 2] == pytest.approx(A6.residual_ey[0], rel=0.05)
    assert math.degrees(res[0, 3]) == pytest.approx(A6.residual_i_deg[0], rel=0.05)
    assert res[1, 0] == pytest.approx(A6.residual_a_m[1], rel=0.01)
    assert res[2, 0] == pytest.approx(A6.residual_a_m[2], rel=0.01)


def test_converted_osculating_elements(converted) -> None:
    # A.6 (janela kepleriana): a 6 878 668,405 m, e 7,270362e-4, i 97,444712°, ω 270°.
    # Com a janela nodal o semieixo fica ≈ 11 m acima (golden: 6 878 679,4 m); fase u0 preservada
    el = converted.elements
    assert el.a == pytest.approx(A6.nodal_osc_a_m, abs=1.0)
    assert el.a == pytest.approx(A6.osc_a_m, abs=15.0)
    assert el.e == pytest.approx(A6.osc_e, rel=0.01)
    assert math.degrees(el.i) == pytest.approx(A6.osc_i_deg, abs=1e-5)
    assert math.degrees(el.argp) == pytest.approx(A6.osc_argp_deg, abs=0.05)
    assert el.raan == DESIGN_MEAN.raan
    # u0 = ω + ν: alvo em u = 90° (ω 90°, ψ 0); osculador em ω 270°, ν ≈ 180°
    u0 = (el.argp + mean_to_true(el.psi, el.e)) % (2 * math.pi)
    assert math.degrees(u0) == pytest.approx(90.0, abs=1e-9)


def test_round_trip(converted) -> None:
    # Ida e volta: osculating_to_mean(mean_to_osculating(m)) ≈ m (Δa ≤ 1e-2 m, Δe ≤ 1e-8,
    # Δi ≤ 1e-8 rad)
    back = osculating_to_mean(converted.elements, EPOCH, MODELS)
    assert abs(back.a - DESIGN_MEAN.a) <= 1e-2
    assert abs(back.e - DESIGN_MEAN.e) <= 1e-8
    assert abs(back.i - DESIGN_MEAN.i) <= 1e-8
    assert math.degrees(back.argp) == pytest.approx(DESIGN_ORBIT.argp_deg, abs=1e-4)
    # também a partir do estado (r, v)
    r0, v0 = CowellPropagator(converted.elements, EPOCH, MODELS).initial_state
    back_state = osculating_to_mean((r0, v0), EPOCH, MODELS)
    assert abs(back_state.a - DESIGN_MEAN.a) <= 1e-2


def test_frozen_orbit_after_conversion(converted) -> None:
    # A.6: estado convertido, Cowell J2–J4, 10 dias, médias por janela nodal:
    # ω ∈ [89,5°; 90,5°], e ∈ [1,06e-3; 1,10e-3], a − alvo ∈ [−50; +50] m
    # (referência A.6: ω 90,00 … 90,16°, e 1,0740e-3 … 1,0781e-3, a − alvo −0,0 … +16,6 m)
    w = _window_series(converted.elements, 10.0)
    argp = np.degrees(np.unwrap(w.argp))
    da = w.a - DESIGN_MEAN.a
    assert A6.crit_argp_deg[0] <= argp.min() and argp.max() <= A6.crit_argp_deg[1]
    assert A6.crit_e[0] <= w.e.min() and w.e.max() <= A6.crit_e[1]
    assert A6.crit_a_minus_target_m[0] <= da.min() and da.max() <= A6.crit_a_minus_target_m[1]
    assert len(w.a) == 151  # 10 dias / Tnod (A.6: 151 janelas)


def test_mean_used_as_osculating_does_not_freeze() -> None:
    # Contraste (A.6): com o médio usado como osculador, ω varia > 15° em 10 dias (referência
    # 21,5°), e ≈ 2,8e-3 e a − alvo ≈ +9,5 km
    w = _window_series(DESIGN_MEAN, 10.0)
    argp = np.degrees(np.unwrap(w.argp))
    assert np.ptp(argp) > 15.0
    assert np.ptp(argp) == pytest.approx(A6.naive_10d_argp_range_deg, abs=0.5)
    assert np.mean(w.e) == pytest.approx(2.8e-3, rel=0.05)
    assert np.mean(w.a - DESIGN_MEAN.a) / 1e3 == pytest.approx(9.5, abs=0.1)


def test_j2_only_conversion_converges() -> None:
    # Com models=(J2Acceleration(),) a conversão também converge (caso da S2-06)
    res = mean_to_osculating(DESIGN_MEAN, EPOCH, (J2Acceleration(),))
    assert res.converged and res.iterations <= 6
    assert res.residuals[0, 0] / 1e3 == pytest.approx(9.5, abs=0.1)


def test_not_converged_warns(caplog: pytest.LogCaptureFixture) -> None:
    # iters esgotado → converged False e aviso
    with caplog.at_level("WARNING", logger="satsim.mean_osc"):
        res = mean_to_osculating(DESIGN_MEAN, EPOCH, MODELS, iters=2)
    assert not res.converged and res.iterations == 2
    assert "não convergiu" in caplog.text
