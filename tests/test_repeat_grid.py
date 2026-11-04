"""Nós ascendentes, grade de 91 trilhas e repetição calibrada (card S5-02).

Referências: Apêndice A.2 (teoria, Sem4 §9.2), A.3 (grade e precisão do nó no Cowell completo) e
A.4 (fechamento do ciclo e semieixo calibrado) do prompt da Semana 5, em
tests/golden/reference_values.py. A constelação P = 6 no Cowell completo é a de referência do
conftest (integrada uma vez por sessão); a calibração do semieixo (J2–J4, ≈ 3,6 s) e as
propagações de 6 dias usam só a gravidade zonal, como no A.4.
"""

import math
from datetime import datetime, timedelta

import numpy as np
import pytest
from golden.reference_values import (
    CONSTELLATION_A1,
    MEAN_GRID_A2,
    NODE_GRID_A3,
    REPEAT_CLOSURE_A4,
)
from scipy.optimize import brentq

from satsim.astrotime import datetime_to_jd, gmst_rate
from satsim.constants import OMEGA_EARTH, OMEGA_SUN, R_EARTH, SECONDS_PER_DAY
from satsim.constellation import Constellation, sso_constellation
from satsim.elements import KeplerianElements
from satsim.forces import ZonalAcceleration
from satsim.groundtrack import grid_metrics, node_crossings, repeat_closure, subpoint
from satsim.mean_osc import mean_to_osculating
from satsim.orbit_design import (
    constellation_phases,
    frozen_eccentricity,
    ltan_drift_min_per_year,
    repeat_orbit,
    repeat_semi_major_axis_calibrated,
    sso_inclination,
)
from satsim.propagators import CowellPropagator, KeplerPropagator, MeanJ2Propagator

A1, A2, A3, A4 = CONSTELLATION_A1, MEAN_GRID_A2, NODE_GRID_A3, REPEAT_CLOSURE_A4
A1_EPOCH = datetime.fromisoformat(A1.epoch_utc)
DAY = SECONDS_PER_DAY
I_A1 = math.radians(A1.i_calibrated_deg)
E_A1 = frozen_eccentricity(A1.a_design_m, float(sso_inclination(A1.a_design_m)))
RAAN_A1 = math.radians(A1.raan_deg)
ZONAL = (ZonalAcceleration(),)


def _deg(x: float) -> float:
    return math.degrees(x)


def _constellation_nodes(sats, start, duration_s: float, step_s: float = 30.0) -> np.ndarray:
    """Longitudes dos nós de todos os satélites, satélite a satélite (S1 primeiro: a âncora)."""
    return np.concatenate([node_crossings(s, start, duration_s, step_s)[1] for s in sats])


# --- grid_metrics e node_crossings: casos sintéticos --------------------------------------------


def test_grid_metrics_uniform_grid_with_wraparound() -> None:
    step = 2.0 * np.pi / 91
    lon = np.mod(0.3 - step * 6 * np.arange(91) + np.pi, 2 * np.pi) - np.pi  # passo de D = 6
    count, gap_min, gap_max, offset = grid_metrics(lon, 91)
    assert count == 91
    assert gap_min == pytest.approx(step, abs=1e-12)
    assert gap_max == pytest.approx(step, abs=1e-12)
    assert offset == pytest.approx(0.0, abs=1e-12)


def test_grid_metrics_offset_and_gaps() -> None:
    step = 2.0 * np.pi / 4
    lon = np.array([0.0, step + 0.1, 2 * step, 3 * step - 0.05])
    count, gap_min, gap_max, offset = grid_metrics(lon, 4)
    assert count == 4
    assert gap_min == pytest.approx(step - 0.1)
    assert gap_max == pytest.approx(step + 0.1)
    assert offset == pytest.approx(0.1)
    # um nó só: o espaçamento é a volta inteira
    assert grid_metrics([1.0], 91)[1:3] == pytest.approx((2 * np.pi, 2 * np.pi))
    with pytest.raises(ValueError):
        grid_metrics([], 91)
    with pytest.raises(ValueError):
        grid_metrics([0.0], 0)


def test_repeat_closure_and_errors() -> None:
    lon = np.array([3.1, 0.0, -3.1])
    assert repeat_closure(lon, 2) == pytest.approx(2 * np.pi - 6.2)  # −3,1 − 3,1 + 2π
    with pytest.raises(ValueError, match="3 nós"):
        repeat_closure(lon[:2], 2)


def test_node_crossings_window_and_epoch() -> None:
    # Kepler quase polar; a janela começa 1 h depois da época do propagador
    el = KeplerianElements(7000e3, 0.0, math.radians(98.0), 0.2, 0.0, 0.0)
    prop = KeplerPropagator(el, A1_EPOCH)
    start = A1_EPOCH + timedelta(hours=1)
    t, lon = node_crossings(prop, start, 3 * 6000.0)
    t_ref, lon_ref = node_crossings(prop, A1_EPOCH, 3600.0 + 3 * 6000.0, step_s=30.0)
    later = t_ref > 3600.0
    # mesma física, eixos de tempo deslocados (a amostragem é outra: interpolação ~1e-5)
    assert t == pytest.approx(t_ref[later] - 3600.0, abs=0.02)
    assert np.all(np.diff(t) > 0) and np.all((lon >= -np.pi) & (lon < np.pi))
    assert len(lon) == len(lon_ref[later])
    with pytest.raises(ValueError):
        node_crossings(prop, start, 100.0, step_s=0.0)
    with pytest.raises(ValueError):
        node_crossings(prop, start, -1.0)


# --- Teoria, modelo médio (Sem4 §9.2, A.2) ------------------------------------------------------


def _mean_model_nodes(n_sats: int, days: float) -> tuple[int, float, float, float]:
    """Grade do MeanJ2Propagator na órbita 91/6 de projeto (e = 0, como o Sem4 §9.2).

    A janela começa 1 s antes da época, para que o nó de S1 em t = 0 entre (o Sem4 conta
    t0 + k·Tnod ≤ D dias).
    """
    a, i = repeat_orbit(91, 6)
    sats = [
        MeanJ2Propagator(KeplerianElements(a, 0.0, i, 0.3, 0.0, float(psi)), A1_EPOCH)
        for psi in constellation_phases(n_sats)
    ]
    start = A1_EPOCH - timedelta(seconds=1)
    return grid_metrics(_constellation_nodes(sats, start, days * DAY), 91)


def _gmst_drift_deg(days: float) -> float:
    """Deriva da grade por a teoria usar ω⊕ = 7,292115e-5 e o ECEF girar com o GMST real."""
    rate = float(gmst_rate(datetime_to_jd(A1_EPOCH)))
    return math.degrees(abs(rate - OMEGA_EARTH) * days * DAY)


@pytest.mark.parametrize(("n_sats", "days"), [(1, 6.0), (6, 1.0)])
def test_mean_model_grid_uniform(n_sats: int, days: float) -> None:
    """P = 1 em 6 dias e P = 6 em 1 dia: 91 nós a 3,9560° (teorema P = D, Sem4 §8.3).

    Tolerância: arredondamento do A.2 (5e-5°) mais a deriva do GMST real contra o ω⊕ da teoria
    (2,5e-4° em 6 dias, 4,2e-5° em 1 dia). Medido: 3,955985 … 3,956270° (6 dias) e
    3,956000 … 3,956059° (1 dia).
    """
    count, gap_min, gap_max, _ = _mean_model_nodes(n_sats, days)
    tol = 5e-5 + _gmst_drift_deg(days)
    assert count == A2.nodes_p1_6days == A2.nodes_p6_1day
    assert _deg(gap_min) == pytest.approx(A2.grid_ideal_deg, abs=tol)
    assert _deg(gap_max) == pytest.approx(A2.grid_ideal_deg, abs=tol)


def test_mean_model_grid_p7_irregular() -> None:
    # P = 7 na grade D = 6: 107 nós, passo de 1/7 a 6/7 da grade (Sem4 §8.4)
    count, gap_min, gap_max, _ = _mean_model_nodes(7, 1.0)
    assert count == A2.nodes_p7_1day
    assert _deg(gap_min) == pytest.approx(A2.step_p7_min_deg, abs=5e-5)
    assert _deg(gap_max) == pytest.approx(A2.step_p7_max_deg, abs=5e-5)


# --- Grade no Cowell completo (A.3) --------------------------------------------------------------


@pytest.fixture(scope="module")
def p7_constellation(a1_models) -> Constellation:
    return sso_constellation(A1_EPOCH, a1_models, n_sats=7, i=I_A1)


def test_cowell_grid_p6_one_day(reference_constellation) -> None:
    """Teorema P = D no Cowell completo: 91 trilhas uniformes em 1 dia com 6 satélites."""
    lon = _constellation_nodes(reference_constellation, A1_EPOCH, DAY)
    count, gap_min, gap_max, offset = grid_metrics(lon, 91)
    assert count == A3.p6_nodes
    assert _deg(gap_min) == pytest.approx(A3.p6_gap_min_deg, abs=A3.p6_gap_tol_deg)
    assert _deg(gap_max) == pytest.approx(A3.p6_gap_max_deg, abs=A3.p6_gap_tol_deg)
    assert _deg(offset) <= A3.p6_max_offset_limit_deg  # medido 0,006447° (A.3: 0,0064°)


def test_cowell_grid_p7_breaks(p7_constellation) -> None:
    # P = 7 na mesma órbita: a grade quebra (A.3: 106 nós, afastamento 1,694°)
    lon = _constellation_nodes(p7_constellation, A1_EPOCH, DAY)
    count, gap_min, gap_max, offset = grid_metrics(lon, 91)
    assert count == A3.p7_nodes
    assert _deg(offset) > A3.p7_max_offset_limit_deg
    assert _deg(offset) == pytest.approx(A3.p7_max_offset_deg, abs=0.01)
    assert _deg(gap_min) == pytest.approx(A3.p7_gap_min_deg, abs=A3.p6_gap_tol_deg)
    assert _deg(gap_max) == pytest.approx(A3.p7_gap_max_deg, abs=A3.p6_gap_tol_deg)


# --- Precisão do nó (A.3) -----------------------------------------------------------------------


@pytest.fixture(scope="module")
def refined_nodes(reference_constellation) -> tuple:
    """S1 em 24 h: nós refinados por brentq em z_ECI = 0 (instante) e a longitude neles."""
    s1 = reference_constellation[0]
    t_guess, _ = node_crossings(s1, A1_EPOCH, DAY, step_s=30.0)

    def z(t: float) -> float:
        return float(s1.propagate(t)[0][0, 2])

    t_ref = np.array([brentq(z, t - 60.0, t + 60.0, xtol=1e-9) for t in t_guess])
    _, lon_ref, _ = subpoint(s1.propagator, t_ref)
    return s1, t_ref, lon_ref


def _node_errors(refined, step_s: float) -> tuple[float, float]:
    s1, t_ref, lon_ref = refined
    t, lon = node_crossings(s1, A1_EPOCH, DAY, step_s=step_s)
    assert len(t) == len(t_ref) == 15
    dlon = np.abs(np.mod(lon - lon_ref + np.pi, 2 * np.pi) - np.pi)
    return math.degrees(float(dlon.max())), float(np.abs(t - t_ref).max())


def test_node_precision_30s(refined_nodes) -> None:
    # A.3: 9,3e-6° (≈ 1 m) e 3e-4 s; critério ≤ 2e-5° e ≤ 1e-3 s
    err_deg, err_s = _node_errors(refined_nodes, 30.0)
    assert err_deg <= A3.node_err_30s_limit_deg
    assert err_s <= A3.node_err_30s_limit_s


def test_node_precision_falls_with_step(refined_nodes) -> None:
    # A.3: 60 s → 7,3e-5°; 30 s → 9,3e-6°; 10 s → 4,8e-7° (erro de interpolação ∝ passo²)
    err = {step: _node_errors(refined_nodes, step)[0] for step in (60.0, 30.0, 10.0)}
    assert err[10.0] < err[30.0] < err[60.0]
    assert err[60.0] == pytest.approx(A3.node_err_60s_deg, rel=0.2)
    assert err[10.0] == pytest.approx(A3.node_err_10s_deg, rel=0.2)


# --- Repetição: fechamento do ciclo e semieixo calibrado (A.4) -----------------------------------


def _s1_zonal(a: float) -> CowellPropagator:
    mean = KeplerianElements(a, E_A1, I_A1, RAAN_A1, math.pi / 2, 0.0)
    return CowellPropagator(mean_to_osculating(mean, A1_EPOCH, ZONAL).elements, A1_EPOCH, ZONAL)


def _closure_km(prop) -> float:
    _, lon = node_crossings(prop, A1_EPOCH, 6 * DAY + 2 * 5697.0)
    return repeat_closure(lon, 91) * R_EARTH / 1e3


@pytest.fixture(scope="module")
def calibrated_a() -> float:
    return repeat_semi_major_axis_calibrated(91, 6, E_A1, I_A1, A1_EPOCH, ZONAL, raan=RAAN_A1)


def test_design_a_does_not_close() -> None:
    # A.4: com o a de projeto (teoria só J2) o ciclo J2–J4 fica aberto de +1,00 km
    closure = _closure_km(_s1_zonal(A4.a_design_m))
    assert closure == pytest.approx(A4.closure_design_j2j4_km, abs=A4.closure_design_tol_km)


def test_calibrated_semi_major_axis_closes_cycle(calibrated_a) -> None:
    # A.4: Δa = +19,37 m e o ciclo fecha a −0,022 km (critério < 0,1 km)
    assert calibrated_a - A4.a_design_m == pytest.approx(A4.delta_a_m, abs=A4.delta_a_tol_m)
    closure = _closure_km(_s1_zonal(calibrated_a))
    assert abs(closure) < A4.closure_calibrated_limit_km


def test_delta_a_effect_on_ltan_is_negligible(calibrated_a) -> None:
    # Sem4 eq. 32 / §6.1: Ω̇ ∝ a^(−7/2): +19 m → −0,014 min/ano de LTAN
    delta_a = calibrated_a - A4.a_design_m
    raan_dot = OMEGA_SUN * (1.0 - 3.5 * delta_a / A4.a_design_m)
    assert ltan_drift_min_per_year(raan_dot) == pytest.approx(
        A4.ltan_drift_delta_a_min_year, abs=0.001
    )
