"""Testes do Cowell com as premissas do SGP4 (``satsim.sgp4_equivalent``).

Referências: a biblioteca ``sgp4`` (constantes WGS-72, regra de s e, na versão em Python puro,
o coeficiente secular de arrasto ``cc1``); TLEs congelados em tests/data. A ablação completa de
7 dias fica em scripts/compare_sgp4_equivalent.py.
"""

from pathlib import Path

import numpy as np
import pytest
from sgp4.api import WGS72
from sgp4.earth_gravity import wgs72
from sgp4.model import Satrec as PythonSatrec

from satsim.analysis import position_error
from satsim.constants import R_EARTH
from satsim.orbit_determination import fit_initial_state
from satsim.propagators import SGP4Propagator
from satsim.sgp4_equivalent import (
    BSTAR_RHO0,
    WGS72_MU,
    WGS72_RE,
    WGS72_ZONAL,
    SGP4Atmosphere,
    cd_a_m_from_bstar,
    perigee_altitude_from_tle,
    sgp4_drag,
    sgp4_equivalent_models,
    sgp4_s_altitude,
)
from satsim.tle import Tle

DATA = Path(__file__).parent / "data"
TLE_FILES = ("iss_sample.tle", "sentinel2a_sample.tle")


def _tle(name: str) -> Tle:
    return Tle.from_text((DATA / name).read_text(encoding="utf-8"))


def test_wgs72_constants_match_sgp4_library() -> None:
    assert WGS72_MU == pytest.approx(398600.8e9, rel=1e-15)
    assert WGS72_RE == pytest.approx(6378.135e3, rel=1e-15)
    assert WGS72_ZONAL == {2: wgs72.j2, 3: wgs72.j3, 4: wgs72.j4}


@pytest.mark.parametrize(
    ("perigee_km", "s_km"),
    [(None, 78.0), (400.0, 78.0), (156.0, 78.0), (130.0, 52.0), (98.0, 20.0), (90.0, 20.0)],
)
def test_s_altitude_follows_sgp4init(perigee_km, s_km) -> None:
    perigee = None if perigee_km is None else perigee_km * 1e3
    assert sgp4_s_altitude(perigee) == pytest.approx(s_km * 1e3)


def test_power_law_density() -> None:
    atm = SGP4Atmosphere()
    # em r = q0 a lei vale exatamente ρ0/R⊕ (h é a altitude esférica sobre o R⊕ do projeto)
    h_q0 = WGS72_RE + 120e3 - R_EARTH
    assert atm.density(h_q0) == pytest.approx(BSTAR_RHO0 / WGS72_RE, rel=1e-14)
    h = np.array([300e3, 600e3])
    rho = atm.density(h)
    r_minus_s = h + R_EARTH - (WGS72_RE + 78e3)
    assert rho[0] / rho[1] == pytest.approx((r_minus_s[1] / r_minus_s[0]) ** 4, rel=1e-14)


def test_drag_is_bstar_power_law_without_corotation() -> None:
    bstar = 3.8792e-5
    drag = sgp4_drag(bstar)
    r = np.array([6800e3, 1200e3, -900e3])
    v = np.array([-1500.0, 6000.0, 4000.0])
    s = WGS72_RE + 78e3
    q0 = WGS72_RE + 120e3
    expected = -(bstar / WGS72_RE) * ((q0 - s) / (np.linalg.norm(r) - s)) ** 4
    expected = expected * np.linalg.norm(v) * v  # v_rel = v: sem co-rotação
    np.testing.assert_allclose(drag(0.0, r, v), expected, rtol=1e-12)
    assert drag.cd_a_m == pytest.approx(cd_a_m_from_bstar(bstar), rel=1e-15)


@pytest.mark.parametrize("tle_file", TLE_FILES)
def test_semimajor_axis_decay_matches_sgp4(tle_file: str) -> None:
    # O SGP4 aplica a queda secular a(t) = a0·(1 − cc1·t)², isto é, da/dt = −2·a0·cc1 em t = 0
    # (ER/min). A mesma taxa, pelo arrasto equivalente: média em uma volta de
    # da/dt = (2a²/µ)·(v·a_D) ao longo das posições do próprio SGP4.
    tle = _tle(tle_file)
    py_sat = PythonSatrec.twoline2rv(tle.line1, tle.line2, WGS72)
    dadt_sgp4 = -2.0 * py_sat.a * py_sat.cc1 * WGS72_RE / 60.0  # m/s
    sgp4 = SGP4Propagator(tle)
    t = np.linspace(0.0, sgp4.period, 2000, endpoint=False)
    r, v = sgp4.propagate(t)
    acc = sgp4_drag(tle.bstar)(0.0, r, v)
    a = py_sat.a * WGS72_RE
    dadt = float(np.mean(2.0 * a**2 / WGS72_MU * np.sum(v * acc, axis=1)))
    # medido: razão 1,0029 (ISS) e 0,9993 (Sentinel-2A); a diferença é o termo de J2 de cc2
    assert dadt == pytest.approx(dadt_sgp4, rel=0.01)


def test_perigee_from_tle() -> None:
    tle = _tle("iss_sample.tle")
    perigee = perigee_altitude_from_tle(tle.mean_motion_rev_day, tle.eccentricity)
    assert 400e3 < perigee < 420e3  # ISS em dez/2019: ≈ 410 km


def test_cowell_with_sgp4_premises_recovers_sgp4() -> None:
    # Com a mesma física do SGP4 (e µ WGS-72 no termo central, inclusive no propagador do
    # resultado), o ajuste do estado inicial em 12 h fica a metros do SGP4. Medido: RMS 3,8 m e
    # 4,3 m em t = 12 h (≈ 3 s). Se o µ do WGS-72 não chegasse ao propagador do resultado, o
    # mesmo estado erraria 570 m em 12 h. Em 12 h o Sentinel-2A ainda não separa os modelos (o
    # J2–J4 do projeto sem arrasto também fica em ≈ 6 m de RMS); a separação aparece em 7 dias,
    # na ablação do script.
    tle = _tle("sentinel2a_sample.tle")
    sgp4 = SGP4Propagator(tle)
    window = 12 * 3600.0
    perigee = perigee_altitude_from_tle(tle.mean_motion_rev_day, tle.eccentricity)
    fit = fit_initial_state(sgp4, sgp4_equivalent_models(tle.bstar, perigee), window_s=window,
                            fit_drag=False, mu=WGS72_MU)  # fmt: skip
    assert fit.mu == WGS72_MU
    assert fit.rms_km < 0.010
    err = position_error(fit.propagator(rtol=1e-10, atol=1e-3), sgp4, [window])[0]
    assert err < 15.0
