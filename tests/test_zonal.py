"""Testes da gravidade zonal J2–J4 (S4-01).

Referências: tests/golden/reference_values.py (ZONAL_POINTS_A2_S4 e ZONAL_A2_S4, Apêndice A.2 do
prompt da Semana 4; DESIGN_ORBIT, Sem4 §7.5), Sem4 eq. 6 (potencial zonal) e §9.3 (J2).
"""

import math
from datetime import datetime

import numpy as np
import pytest
from golden.reference_values import DESIGN_ORBIT, ZONAL_A2_S4, ZONAL_POINTS_A2_S4

from satsim.analysis import fit_rate, osculating_elements
from satsim.constants import MU, R_EARTH
from satsim.elements import KeplerianElements
from satsim.forces import ZonalAcceleration, j2_acceleration, zonal_acceleration, zonal_potential
from satsim.propagators import CowellPropagator

DEG_DAY = 86400.0 * 180.0 / math.pi
EPOCH = datetime.fromisoformat(ZONAL_A2_S4.epoch_utc)
# Elementos de projeto (Sem4 §7.5) tomados como osculadores, Ω0 = 30°, ψ0 = 0 (A.2)
DESIGN_AS_OSC = KeplerianElements(
    a=DESIGN_ORBIT.a_km * 1e3,
    e=DESIGN_ORBIT.e,
    i=math.radians(DESIGN_ORBIT.i_deg),
    raan=math.radians(ZONAL_A2_S4.raan0_deg),
    argp=math.radians(DESIGN_ORBIT.argp_deg),
    psi=math.radians(ZONAL_A2_S4.psi0_deg),
)


def _random_leo_points(n: int, seed: int = 4) -> np.ndarray:
    """n pontos com direção uniforme na esfera e altitude uniforme em 200–2000 km."""
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    return u * (R_EARTH + rng.uniform(200e3, 2000e3, size=n))[:, None]


def _rel(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """|a − b|/|b| por linha (vetorial)."""
    return np.linalg.norm(a - b, axis=-1) / np.linalg.norm(b, axis=-1)


def test_degree_2_reproduces_j2_acceleration() -> None:
    # A.2: ZonalAcceleration((2,)) = j2_acceleration (Sem4 §9.3) em 200 pontos de LEO,
    # rel ≤ 1e-12 (referência 1,3e-15)
    pts = _random_leo_points(200)
    model = ZonalAcceleration((2,))
    acc = np.array([model(0.0, r, np.zeros(3)) for r in pts])
    assert np.max(_rel(acc, j2_acceleration(pts))) <= 1e-12


def test_acceleration_is_gradient_of_potential() -> None:
    # A.2: diferenças centrais de zonal_potential (h = 1 m) × aceleração, J2–J4, rel ≤ 1e-7
    # (referência 1,4e-9)
    step = 1.0
    pts = _random_leo_points(50, seed=5)
    grad = np.stack(
        [
            (zonal_potential(pts + step * e_k) - zonal_potential(pts - step * e_k)) / (2.0 * step)
            for e_k in np.eye(3)
        ],
        axis=-1,
    )
    assert np.max(_rel(zonal_acceleration(pts), grad)) <= 1e-7


@pytest.mark.parametrize("point", ZONAL_POINTS_A2_S4, ids=lambda p: p.label)
@pytest.mark.parametrize("degree", [2, 3, 4])
def test_reference_points_a2(point, degree: int) -> None:
    # ZONAL_POINTS_A2_S4 (A.2): aJ2, aJ3 e aJ4 no equador, a 45° N e no polo, r = 7000 km,
    # rel 1e-9 sobre o vetor (componentes nulas exigem |a_k| ≤ 1e-9·|a|)
    r = np.array(point.r_km) * 1e3
    expected = np.array({2: point.a_j2, 3: point.a_j3, 4: point.a_j4}[degree])
    acc = zonal_acceleration(r, (degree,))
    assert np.linalg.norm(acc - expected) <= 1e-9 * np.linalg.norm(expected)


@pytest.mark.parametrize("degree", [2, 3, 4])
def test_parity_under_z_reflection(degree: int) -> None:
    # Sem4 eq. 6: P_n(−s) = (−1)^n P_n(s). Refletindo z → −z, os termos pares (J2, J4) refletem
    # a componente z; o ímpar (J3) inverte x e y e mantém z.
    pts = _random_leo_points(30, seed=6)
    mirrored = pts * np.array([1.0, 1.0, -1.0])
    acc = zonal_acceleration(pts, (degree,))
    acc_m = zonal_acceleration(mirrored, (degree,))
    sign = np.array([1.0, 1.0, -1.0]) if degree % 2 == 0 else np.array([-1.0, -1.0, 1.0])
    np.testing.assert_allclose(acc_m, sign * acc, rtol=1e-12, atol=1e-20)
    # o potencial é par em z para n par e ímpar para n ímpar
    pot, pot_m = zonal_potential(pts, (degree,)), zonal_potential(mirrored, (degree,))
    np.testing.assert_allclose(pot_m, (-1.0) ** degree * pot, rtol=1e-12)


def test_vectorized_matches_single() -> None:
    # (N, 3) igual a chamadas (3,), para a aceleração e o potencial
    pts = _random_leo_points(10, seed=7)
    acc, pot = zonal_acceleration(pts), zonal_potential(pts)
    assert acc.shape == (10, 3) and pot.shape == (10,)
    for k, r in enumerate(pts):
        single = zonal_acceleration(r)
        assert single.shape == (3,) and np.shape(zonal_potential(r)) == ()
        np.testing.assert_array_equal(acc[k], single)
        assert pot[k] == zonal_potential(r)


def test_model_matches_functions_and_validates_degrees() -> None:
    model = ZonalAcceleration((2, 4))
    r = _random_leo_points(1)[0]
    np.testing.assert_array_equal(model(0.0, r, np.zeros(3)), zonal_acceleration(r, (2, 4)))
    assert model.potential(r) == zonal_potential(r, (2, 4))
    # superposição: a soma dos graus separados é o modelo completo
    total = sum(zonal_acceleration(r, (n,)) for n in (2, 3, 4))
    np.testing.assert_allclose(zonal_acceleration(r), total, rtol=1e-14)
    with pytest.raises(ValueError):
        ZonalAcceleration((5,))  # sem coeficiente J5
    with pytest.raises(ValueError):
        ZonalAcceleration((1,))
    # coeficiente próprio para outro grau (a recorrência de Legendre é geral)
    j5 = ZonalAcceleration((5,), coefficients={5: -2.27e-7})
    assert np.all(np.isfinite(j5(0.0, r, np.zeros(3))))


def test_energy_conserved_in_cowell_j2_j4() -> None:
    # A.2: Cowell J2–J4, 10 dias, rtol = 1e-11: ε = v²/2 − µ/r − R(r) conservada,
    # |Δε/ε| ≤ 1e-9 (referência 3,4e-11)
    prop = CowellPropagator(DESIGN_AS_OSC, EPOCH, (ZonalAcceleration(),), rtol=1e-11)
    t = np.arange(0.0, 10 * 86400.0 + 1.0, 60.0)
    r, v = prop.propagate(t)
    energy = 0.5 * np.sum(v * v, axis=1) - MU / np.linalg.norm(r, axis=1) - zonal_potential(r)
    assert np.max(np.abs(energy / energy[0] - 1.0)) <= 1e-9


@pytest.fixture(scope="module")
def raan_rates() -> dict[tuple[int, ...], float]:
    """Ω̇ [°/dia] por fit_rate, 30 dias a 60 s, para J2, J2+J3 e J2–J4 (A.2). ≈ 15 s no total."""
    t = np.arange(0.0, 30 * 86400.0 + 1.0, 60.0)
    rates = {}
    for degrees in [(2,), (2, 3), (2, 3, 4)]:
        prop = CowellPropagator(DESIGN_AS_OSC, EPOCH, (ZonalAcceleration(degrees),))
        rates[degrees] = fit_rate(t, osculating_elements(prop, t).raan) * DEG_DAY
    return rates


def test_j4_reduces_nodal_rate(raan_rates: dict[tuple[int, ...], float]) -> None:
    # ZONAL_A2_S4 (A.2): Ω̇(J2–J4)/Ω̇(J2) = 0,997703 (abs 2e-5); valores absolutos 0,979419 e
    # 0,977169 °/dia (conferidos a 2e-6 °/dia, a precisão impressa)
    ratio = raan_rates[(2, 3, 4)] / raan_rates[(2,)]
    assert ratio == pytest.approx(ZONAL_A2_S4.raan_dot_ratio_j2j4_j2, abs=2e-5)
    assert raan_rates[(2,)] == pytest.approx(ZONAL_A2_S4.raan_dot_j2_deg_day, abs=2e-6)
    assert raan_rates[(2, 3, 4)] == pytest.approx(ZONAL_A2_S4.raan_dot_j2_j4_deg_day, abs=2e-6)


def test_j3_alone_does_not_change_nodal_rate(raan_rates: dict[tuple[int, ...], float]) -> None:
    # A.2: J2+J3 × J2, rel ≤ 1e-5 (referência 0,979418 × 0,979419)
    assert raan_rates[(2, 3)] == pytest.approx(raan_rates[(2,)], rel=1e-5)
    assert raan_rates[(2, 3)] == pytest.approx(ZONAL_A2_S4.raan_dot_j2_j3_deg_day, abs=2e-6)
