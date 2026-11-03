"""Testes das efemérides do Sol e da Lua, LTAN, β e do terceiro corpo (S4-02).

Referências: tests/golden/reference_values.py (EPHEMERIS_A3 e THIRD_BODY_A3, Apêndice A.3 do
prompt da Semana 4: JPL DE421 via Skyfield 1.55, quadro da data), Sem4 eqs. 28 (LTAN) e 29 (β),
Montenbruck & Gill (2000) eq. 3.37 (terceiro corpo).
"""

import builtins
import math
import socket
from datetime import datetime

import numpy as np
import pytest
from golden.reference_values import EPHEMERIS_A3, THIRD_BODY_A3

from satsim import constants
from satsim.astrotime import datetime_to_jd
from satsim.constants import AU, R_EARTH
from satsim.ephemeris import (
    beta_angle,
    ltan_hours,
    moon_position_eci,
    sun_direction,
    sun_position_eci,
    sun_ra_dec,
)
from satsim.forces import ThirdBodyAcceleration, third_body_acceleration

T = THIRD_BODY_A3
EPOCH = datetime.fromisoformat(T.epoch_utc)
JD = datetime_to_jd(EPOCH)


def _angle_deg(a: np.ndarray, b: np.ndarray) -> float:
    cos = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return math.degrees(math.acos(min(1.0, max(-1.0, cos))))


def _jd(point) -> float:
    return datetime_to_jd(datetime.fromisoformat(point.date_utc))


def test_constants_match_golden() -> None:
    # A.3: MU_SUN (IERS 2010), MU_MOON (DE421), AU (IAU 2012)
    assert constants.MU_SUN == T.mu_sun_m3_s2
    assert constants.MU_MOON == T.mu_moon_m3_s2
    assert constants.AU == T.au_m


@pytest.mark.parametrize("point", EPHEMERIS_A3, ids=lambda p: p.date_utc[:10])
def test_sun_against_de421(point) -> None:
    # A.3: direção ≤ 0,02° e |Δr|/r ≤ 1e-4 (diferença de distância; ver nota no golden)
    sun = sun_position_eci(_jd(point))
    ref = np.array(point.sun_km) * 1e3
    assert _angle_deg(sun, ref) <= T.sun_dir_tol_deg
    assert abs(np.linalg.norm(sun) / np.linalg.norm(ref) - 1.0) <= T.sun_dist_rel_tol


@pytest.mark.parametrize("point", EPHEMERIS_A3, ids=lambda p: p.date_utc[:10])
def test_moon_against_de421(point) -> None:
    # A.3: direção ≤ 0,1° e |Δr| ≤ 500 km (distância, como no A.3; o vetor também cabe)
    moon = moon_position_eci(_jd(point))
    ref = np.array(point.moon_km) * 1e3
    assert _angle_deg(moon, ref) <= T.moon_dir_tol_deg
    assert abs(np.linalg.norm(moon) - np.linalg.norm(ref)) <= T.moon_dist_tol_km * 1e3
    assert np.linalg.norm(moon - ref) <= T.moon_dist_tol_km * 1e3


def test_sun_right_ascension_and_distance() -> None:
    # A.3: α☉ = 207,5235° (abs 0,001°; DE421 207,5350°) e distância 0,995080 UA em 2026-10-23
    ra, dec = sun_ra_dec(JD)
    assert math.degrees(ra) == pytest.approx(T.sun_ra_aa_deg, abs=1e-3)
    assert math.degrees(ra) == pytest.approx(T.sun_ra_de421_deg, abs=0.02)
    assert -23.5 < math.degrees(dec) < -10.0  # fim de outubro: Sol ao sul do equador
    assert np.linalg.norm(sun_position_eci(JD)) / AU == pytest.approx(T.sun_distance_au, abs=1e-6)
    np.testing.assert_allclose(np.linalg.norm(sun_direction(JD)), 1.0, rtol=1e-15)


def test_ltan() -> None:
    # Sem4 eq. 28 e A.3: Ω = 5,0235° em 2026-10-23 00:00 UTC → 22,500 h (abs 0,005 h); Ω = α☉ → 12 h
    assert ltan_hours(math.radians(T.ltan_raan_deg), JD) == pytest.approx(T.ltan_h, abs=5e-3)
    ra, _ = sun_ra_dec(JD)
    assert ltan_hours(ra, JD) == pytest.approx(12.0, abs=1e-12)
    # resultado em [0, 24) e 1 h a cada 15° de Ω
    raan = np.radians(np.arange(0.0, 360.0, 15.0))
    hours = ltan_hours(raan, JD)
    assert np.all((hours >= 0.0) & (hours < 24.0))
    np.testing.assert_allclose(np.mod(np.diff(hours), 24.0), 1.0, atol=1e-12)


def test_beta_angle() -> None:
    # Sem4 eq. 29: sin β = ŝ·ℓ̂
    s = sun_direction(JD)
    a = R_EARTH + 500e3
    perp = np.cross(s, [0.0, 0.0, 1.0])
    perp /= np.linalg.norm(perp)
    # plano contendo o Sol (r ∥ ŝ): β = 0
    assert beta_angle(a * s, 7.6e3 * perp, JD) == pytest.approx(0.0, abs=1e-12)
    # normal ao plano apontando para o Sol (ℓ ∥ ŝ): β = +90°; invertendo v, −90°
    v = np.cross(s, perp) * 7.6e3
    assert math.degrees(beta_angle(a * perp, v, JD)) == pytest.approx(90.0, abs=1e-6)
    assert math.degrees(beta_angle(a * perp, -v, JD)) == pytest.approx(-90.0, abs=1e-6)


def test_beta_matches_element_formula() -> None:
    # Órbita circular: sin β = cos δ☉·sin i·sin(Ω − α☉) + sin δ☉·cos i (Sem4 §5.2, equivalente à
    # eq. 29 com ℓ̂ = (sin i sin Ω, −sin i cos Ω, cos i))
    from satsim.elements import KeplerianElements, elements_to_state

    ra, dec = sun_ra_dec(JD)
    for i_deg, raan_deg in [(97.46, 5.0), (51.6, 120.0), (97.46, 200.0)]:
        i, raan = math.radians(i_deg), math.radians(raan_deg)
        r, v = elements_to_state(KeplerianElements(7e6, 0.0, i, raan, 0.3, 1.0))
        expected = math.asin(
            math.cos(dec) * math.sin(i) * math.sin(raan - ra) + math.sin(dec) * math.cos(i)
        )
        assert beta_angle(r[0], v[0], JD) == pytest.approx(expected, abs=1e-12)


def test_vectorized_over_jd() -> None:
    # (N,) de JD → (N, 3), igual a chamadas escalares
    jds = JD + np.linspace(0.0, 30.0, 7)
    for func in (sun_position_eci, moon_position_eci, sun_direction):
        out = func(jds)
        assert out.shape == (7, 3)
        for k, jd in enumerate(jds):
            np.testing.assert_array_equal(out[k], func(jd))
    ra, dec = sun_ra_dec(jds)
    assert ra.shape == dec.shape == (7,)
    assert ltan_hours(np.zeros(7), jds).shape == (7,)


@pytest.mark.parametrize(("body", "expected"), [("sun", T.a_sun_m_s2), ("moon", T.a_moon_m_s2)])
def test_third_body_magnitude(body: str, expected: float) -> None:
    # A.3: r = (R⊕ + 500 km, 0, 0) em 2026-10-23 00:00 UTC: |a☉| = 5,0023e-7 e
    # |a☾| = 1,2005e-6 m/s² (rel 1%)
    model = ThirdBodyAcceleration(body, EPOCH)
    r = np.array([R_EARTH + T.r_alt_km * 1e3, 0.0, 0.0])
    assert np.linalg.norm(model(0.0, r, np.zeros(3))) == pytest.approx(expected, rel=0.01)


@pytest.mark.parametrize("body", ["sun", "moon"])
def test_third_body_direction(body: str) -> None:
    # Maré: com o satélite entre a Terra e o corpo, a aceleração aponta para o corpo; do lado
    # oposto, aponta para longe dele (também "para fora" da Terra); na perpendicular, para dentro.
    model = ThirdBodyAcceleration(body, EPOCH)
    s_hat = model.body_position(0.0) / np.linalg.norm(model.body_position(0.0))
    r_mag = R_EARTH + 500e3
    a_near = model(0.0, r_mag * s_hat, np.zeros(3))
    a_far = model(0.0, -r_mag * s_hat, np.zeros(3))
    assert _angle_deg(a_near, s_hat) < 0.01
    assert _angle_deg(a_far, -s_hat) < 0.01
    perp = np.cross(s_hat, [0.0, 0.0, 1.0])
    perp /= np.linalg.norm(perp)
    a_perp = model(0.0, r_mag * perp, np.zeros(3))
    assert np.dot(a_perp, perp) < 0.0
    # maré: |a| no eixo ≈ 2× |a| na perpendicular
    assert np.linalg.norm(a_near) / np.linalg.norm(a_perp) == pytest.approx(2.0, rel=0.05)


def test_third_body_vectorized_and_validated() -> None:
    model = ThirdBodyAcceleration("moon", EPOCH)
    t = np.array([0.0, 3600.0, 86400.0])
    r = np.array([[7e6, 0.0, 0.0], [0.0, 7e6, 0.0], [0.0, 0.0, 7e6]])
    acc = third_body_acceleration(r, model.body_position(t), model.mu)
    assert acc.shape == (3, 3)
    for k in range(3):
        np.testing.assert_array_equal(acc[k], model(t[k], r[k], np.zeros(3)))
    with pytest.raises(ValueError):
        ThirdBodyAcceleration("jupiter", EPOCH)  # type: ignore[arg-type]


def test_no_network_or_files(monkeypatch: pytest.MonkeyPatch) -> None:
    # Efemérides analíticas: nenhuma rede nem arquivo (sem .bsp)
    def _blocked(*args, **kwargs):
        raise AssertionError("acesso a rede/arquivo não permitido")

    monkeypatch.setattr(socket, "socket", _blocked)
    monkeypatch.setattr(builtins, "open", _blocked)
    jds = JD + np.arange(3.0)
    assert np.all(np.isfinite(sun_position_eci(jds)))
    assert np.all(np.isfinite(moon_position_eci(jds)))
    model = ThirdBodyAcceleration("sun", EPOCH)
    assert np.all(np.isfinite(model(10.0, np.array([7e6, 0.0, 0.0]), np.zeros(3))))
