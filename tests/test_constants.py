"""Testes das constantes físicas (S1-03)."""

import math

import pytest

from satsim import constants as c


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        # Sem4, tabela de símbolos
        ("MU", 3.986004418e14),
        ("R_EARTH", 6378.137e3),
        ("OMEGA_EARTH", 7.292115e-5),
        # Sem4 código §9.1 (J2) e §2.2 (J3, J4)
        ("J2", 1.08262668e-3),
        ("J3", -2.533e-6),
        ("J4", -1.620e-6),
        # Sem4, tabela de símbolos
        ("T_TROP_DAYS", 365.2422),
    ],
)
def test_constant_matches_table(name: str, expected: float) -> None:
    assert getattr(c, name) == pytest.approx(expected, rel=1e-12)


def test_helpers() -> None:
    assert c.SECONDS_PER_DAY == 86400.0
    assert c.DEG2RAD == pytest.approx(math.pi / 180.0, rel=1e-15)
    assert c.RAD2DEG * c.DEG2RAD == pytest.approx(1.0, rel=1e-15)


def test_omega_sun_deg_per_day() -> None:
    # Sem4 eq. 26: taxa heliossíncrona 0,985647 °/dia
    assert c.OMEGA_SUN * c.RAD2DEG * c.SECONDS_PER_DAY == pytest.approx(0.985647, abs=1e-6)


def test_omega_sun_rad_per_s() -> None:
    # Sem4 eq. 26: 1,991064e-7 rad/s
    assert c.OMEGA_SUN == pytest.approx(1.991064e-7, rel=1e-6)


def test_earth_rotation_over_sun_rate() -> None:
    # A precessão heliossíncrona é ~366 vezes mais lenta que a rotação da Terra
    # (razão derivada de ω⊕ e Sem4 eq. 26; valor esperado do prompt, card S1-03).
    assert c.OMEGA_EARTH / c.OMEGA_SUN == pytest.approx(366.2, rel=1e-3)


def test_j2_normalization_trap() -> None:
    # Sem4 §2.1: C̄20 = −4,8417e-4 (EGM2008, normalizado); N20 = √5; J2 = −N20·C̄20
    cbar20 = -4.8417e-4
    assert -math.sqrt(5.0) * cbar20 == pytest.approx(c.J2, rel=1e-4)


def test_j3_j4_negative() -> None:
    # Sem4 §2.2: J3 = −2,533e-6, J4 = −1,620e-6
    assert c.J3 < 0
    assert c.J4 < 0
