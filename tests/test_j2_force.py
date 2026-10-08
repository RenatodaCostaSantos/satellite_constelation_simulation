"""Testes da aceleração de J2 (S2-02).

Referências: Sem4 §2.2 (função perturbadora R) e §9.3 (aceleração cartesiana).
"""

import numpy as np
import pytest

from satsim.constants import J2, MU, R_EARTH
from satsim.forces import J2Acceleration, j2_acceleration, j2_potential


def _random_points(n: int, seed: int = 2) -> np.ndarray:
    """n pontos com direção uniforme na esfera e altitude uniforme em 400–2000 km."""
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    return u * (R_EARTH + rng.uniform(400e3, 2000e3, size=n))[:, None]


def test_acceleration_is_gradient_of_potential() -> None:
    # a = ∇R, R = −µ·J2·R⊕²·(3 sin²φ − 1)/(2 r³) (Sem4 §2.2); diferenças centrais com passo 1 m
    step = 1.0
    for r in _random_points(20):
        grad = np.array(
            [
                (j2_potential(r + step * e_k) - j2_potential(r - step * e_k)) / (2.0 * step)
                for e_k in np.eye(3)
            ]
        )
        acc = j2_acceleration(r)
        # rel = 1e-6 sobre o vetor: |a − ∇R| ≤ 1e-6·|a|
        assert np.linalg.norm(acc - grad) <= 1e-6 * np.linalg.norm(acc)


def test_equator_is_radial_and_attractive() -> None:
    # Sem4 §9.3: em r = (r, 0, 0), a = −1,5·J2·µ·R⊕²/r⁴ ao longo de r̂
    r = R_EARTH + 500e3
    acc = j2_acceleration([r, 0.0, 0.0])
    expected = -1.5 * J2 * MU * R_EARTH**2 / r**4
    assert acc[0] == pytest.approx(expected, rel=1e-12)
    assert acc[1] == 0.0 and acc[2] == 0.0


def test_pole_is_along_plus_z_and_repulsive() -> None:
    # Sem4 §9.3: em r = (0, 0, r), a = +3·J2·µ·R⊕²/r⁴ ao longo de +z
    r = R_EARTH + 500e3
    acc = j2_acceleration([0.0, 0.0, r])
    expected = 3.0 * J2 * MU * R_EARTH**2 / r**4
    assert acc[2] == pytest.approx(expected, rel=1e-12)
    assert acc[0] == 0.0 and acc[1] == 0.0


def test_vectorized_matches_single() -> None:
    pts = _random_points(10)
    acc = j2_acceleration(pts)
    pot = j2_potential(pts)
    assert acc.shape == (10, 3) and pot.shape == (10,)
    for k, r in enumerate(pts):
        single = j2_acceleration(r)
        assert single.shape == (3,)
        np.testing.assert_array_equal(acc[k], single)
        assert pot[k] == j2_potential(r)


def test_model_matches_function() -> None:
    model = J2Acceleration()
    r = _random_points(1)[0]
    np.testing.assert_array_equal(model(0.0, r, np.zeros(3)), j2_acceleration(r))
    assert model.potential(r) == j2_potential(r)
    # parâmetros próprios: j2 = 0 anula a aceleração
    np.testing.assert_array_equal(J2Acceleration(j2=0.0)(0.0, r, np.zeros(3)), np.zeros(3))
