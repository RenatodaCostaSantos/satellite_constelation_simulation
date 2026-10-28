"""Testes da pressão de radiação solar, da sombra cilíndrica e da fração de eclipse (S4-03).

Referências: tests/golden/reference_values.py (SRP_ECLIPSE_A4, Apêndice A.4 do prompt da Semana
4), Sem4 §5.3 (sombra cilíndrica) e eq. 30 (fração de eclipse), Montenbruck & Gill (2000) §3.4.
"""

import math
from datetime import datetime

import numpy as np
import pytest
from golden.reference_values import SRP_ECLIPSE_A4

from satsim import constants
from satsim.astrotime import datetime_to_jd
from satsim.constants import AU, R_EARTH
from satsim.elements import KeplerianElements
from satsim.ephemeris import (
    beta_angle,
    eclipse_fraction_circular,
    in_shadow_cylindrical,
    sun_direction,
)
from satsim.forces import SolarRadiationPressure, ZonalAcceleration, srp_acceleration
from satsim.propagators import CowellPropagator

A4 = SRP_ECLIPSE_A4
EPOCH = datetime.fromisoformat(A4.epoch_utc)
JD = datetime_to_jd(EPOCH)
SUN_HAT = sun_direction(JD)


def _perpendicular(s_hat: np.ndarray) -> np.ndarray:
    p = np.cross(s_hat, [0.0, 0.0, 1.0])
    return p / np.linalg.norm(p)


def _angle_deg(a: np.ndarray, b: np.ndarray) -> float:
    cos = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return math.degrees(math.acos(min(1.0, max(-1.0, cos))))


def test_p_sun_constant() -> None:
    # A.4 / M&G §3.4: P☉ = 4,56e-6 N/m² a 1 UA
    assert constants.P_SUN == A4.p_sun_n_m2


def test_magnitude_and_direction() -> None:
    # A.4: C_R = 1,3, A/m = 0,01, iluminado: 5,928e-8·(AU/d)² m/s²; em 2026-10-23 00:00 UTC
    # (d = 0,995080 UA) → 5,9868e-8 (rel 1e-4); direção oposta ao Sol
    model = SolarRadiationPressure(EPOCH, cr=A4.cr, area_over_mass=A4.area_over_mass_m2_kg)
    r = (R_EARTH + 500e3) * _perpendicular(SUN_HAT)  # iluminado, d ≈ distância Terra–Sol
    acc = model(0.0, r, np.zeros(3))
    assert np.linalg.norm(acc) == pytest.approx(A4.a_epoch_m_s2, rel=1e-4)
    d_au = np.linalg.norm(r - model.sun_position(0.0)) / AU
    assert np.linalg.norm(acc) == pytest.approx(A4.a_1au_m_s2 / d_au**2, rel=1e-12)
    assert d_au == pytest.approx(A4.sun_distance_au, abs=1e-6)
    # empurra para longe do Sol: exatamente ao longo de r − s; difere de −ŝ por |r|/d ≈ 0,0026°
    sun_to_sat = r - model.sun_position(0.0)
    unit = acc / np.linalg.norm(acc)
    assert np.linalg.norm(np.cross(unit, sun_to_sat / np.linalg.norm(sun_to_sat))) < 1e-12
    assert np.dot(acc, sun_to_sat) > 0.0
    assert _angle_deg(acc, -SUN_HAT) < 0.01
    # escala linear em C_R e A/m
    double = SolarRadiationPressure(EPOCH, cr=2.6, area_over_mass=A4.area_over_mass_m2_kg)
    np.testing.assert_allclose(double(0.0, r, np.zeros(3)), 2.0 * acc, rtol=1e-14)


def test_shadow_geometry() -> None:
    # Sem4 §5.3: sombra cilíndrica. Atrás da Terra no eixo → sombra; deslocado 1,01·R⊕ do eixo →
    # iluminado; lado do Sol → iluminado
    perp = _perpendicular(SUN_HAT)
    behind = -(R_EARTH + 500e3) * SUN_HAT
    assert in_shadow_cylindrical(behind, SUN_HAT)
    assert in_shadow_cylindrical(behind + 0.99 * R_EARTH * perp, SUN_HAT)
    assert not in_shadow_cylindrical(-2.0 * R_EARTH * SUN_HAT + 1.01 * R_EARTH * perp, SUN_HAT)
    assert not in_shadow_cylindrical((R_EARTH + 500e3) * SUN_HAT, SUN_HAT)
    assert not in_shadow_cylindrical((R_EARTH + 500e3) * perp, SUN_HAT)  # terminador
    # vetorizado
    pts = np.array([behind, (R_EARTH + 500e3) * SUN_HAT])
    np.testing.assert_array_equal(in_shadow_cylindrical(pts, SUN_HAT), [True, False])


def test_zero_acceleration_in_shadow() -> None:
    model = SolarRadiationPressure(EPOCH)
    behind = -(R_EARTH + 500e3) * SUN_HAT
    np.testing.assert_array_equal(model(0.0, behind, np.zeros(3)), np.zeros(3))
    # sem sombra, a mesma posição é empurrada (para longe do Sol)
    no_shadow = SolarRadiationPressure(EPOCH, shadow="none")
    assert np.linalg.norm(no_shadow(0.0, behind, np.zeros(3))) > 5e-8
    with pytest.raises(ValueError):
        SolarRadiationPressure(EPOCH, shadow="conical")  # type: ignore[arg-type]


def test_vectorized_matches_single() -> None:
    model = SolarRadiationPressure(EPOCH)
    t = np.array([0.0, 1800.0, 3600.0])
    r = np.array([-(R_EARTH + 5e5) * SUN_HAT, [7e6, 0.0, 0.0], [0.0, 7e6, 0.0]])
    acc = srp_acceleration(r, model.sun_position(t), model.cr, model.area_over_mass)
    assert acc.shape == (3, 3)
    for k in range(3):
        np.testing.assert_array_equal(acc[k], model(t[k], r[k], np.zeros(3)))


def test_eclipse_fraction_formula() -> None:
    # Sem4 eq. 30 e A.4: caso meia-noite (β −0,94°) → 0,3767; β = 0 → β*/π; |β| > β* → 0
    a = A4.demo_a_km * 1e3
    f = eclipse_fraction_circular(a, math.radians(A4.midnight_beta_deg))
    assert f == pytest.approx(A4.midnight_eclipse_fraction_eq30, abs=1e-4)
    beta_star = math.asin(R_EARTH / a)
    assert eclipse_fraction_circular(a, 0.0) == pytest.approx(beta_star / math.pi, rel=1e-12)
    assert eclipse_fraction_circular(a, beta_star + 1e-3) == 0.0
    assert eclipse_fraction_circular(a, -(beta_star + 1e-3)) == 0.0
    betas = np.radians([0.0, 23.42, 60.0, 80.0])
    assert eclipse_fraction_circular(a, betas).shape == (4,)


def test_eclipse_fraction_demo_sso() -> None:
    # A.4: SSO de demonstração (LTAN 22h30), Cowell J2–J4, 24 h a 20 s, sombra cilíndrica:
    # fração 0,360 (abs 0,01; A.4: 0,3603); β 23,41…23,43°; eq. 30 com β médio 0,3650 (abs 0,01),
    # ≈ 34,7 min por volta (Sem4 §5.3: "~35 min")
    elements = KeplerianElements(
        a=A4.demo_a_km * 1e3,
        e=A4.demo_e,
        i=math.radians(A4.demo_i_deg),
        raan=math.radians(A4.demo_raan_deg),
        argp=math.radians(A4.demo_argp_deg),
        psi=0.0,
    )
    t = np.arange(0.0, 86400.0 + 1.0, 20.0)
    r, v = CowellPropagator(elements, EPOCH, (ZonalAcceleration(),)).propagate(t)
    jd = JD + t / 86400.0
    fraction = float(np.mean(in_shadow_cylindrical(r, sun_direction(jd))))
    beta = np.degrees(beta_angle(r, v, jd))
    assert fraction == pytest.approx(A4.demo_eclipse_fraction, abs=0.01)
    assert beta.min() == pytest.approx(A4.demo_beta_min_deg, abs=0.02)
    assert beta.max() == pytest.approx(A4.demo_beta_max_deg, abs=0.02)
    eq30 = eclipse_fraction_circular(elements.a, math.radians(beta.mean()))
    assert eq30 == pytest.approx(A4.demo_eclipse_fraction_eq30, abs=0.01)
    assert fraction == pytest.approx(eq30, abs=0.01)
    assert eq30 * A4.demo_period_s / 60.0 == pytest.approx(A4.demo_eclipse_min_per_rev, abs=0.2)
