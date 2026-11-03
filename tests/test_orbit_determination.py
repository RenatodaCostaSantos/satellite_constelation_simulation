"""Testes da determinação de órbita mínima (S4-06).

Referências: tests/golden/reference_values.py (HIGH_FIDELITY_A8, Apêndice A.8 do prompt da
Semana 4); TLE do Sentinel-2A congelado em tests/data/sentinel2a_sample.tle (A.1). O caso de 72 h
com todas as forças fica em scripts/compare_high_fidelity.py.
"""

from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pytest
from golden.reference_values import HIGH_FIDELITY_A8

from satsim.analysis import position_error
from satsim.forces import DragAcceleration, J2Acceleration, ZonalAcceleration
from satsim.orbit_determination import fit_initial_state
from satsim.propagators import CowellPropagator, SGP4Propagator
from satsim.tle import Tle

A8 = HIGH_FIDELITY_A8
EPOCH = datetime(2026, 10, 23, tzinfo=UTC)
SENTINEL_TLE = Path(__file__).parent / "data" / "sentinel2a_sample.tle"
R0 = np.array([6878e3, 0.0, 0.0])
V0 = np.array([0.0, -1000.0, 7547.0])


def test_fit_recovers_synthetic_state() -> None:
    # Referência = Cowell J2 a partir de um estado conhecido; chute errado em 1 km e 1 m/s
    reference = CowellPropagator((R0, V0), EPOCH, (J2Acceleration(),))
    guess = (R0 + [1e3, 0.0, 0.0], V0 + [0.0, 1.0, 0.0])
    fit = fit_initial_state(reference, (J2Acceleration(),), window_s=6 * 3600.0, fit_drag=False,
                            initial=guess, rtol=1e-11, atol=1e-6)  # fmt: skip
    assert np.linalg.norm(fit.r0 - R0) <= A8.synthetic_pos_tol_m
    assert np.linalg.norm(fit.v0 - V0) <= A8.synthetic_vel_tol_m_s
    assert fit.rms_initial_km > 10.0  # o chute erra dezenas de km ao longo da janela
    # ~0: a referência é o mesmo modelo. O piso é o ruído da integração (o DOP853 adaptativo muda
    # a sequência de passos a cada estado testado): de 0,5 a 8 mm em 6 h conforme o chute e a
    # máquina (CI: 1,4 cm). 10 cm fica acima do piso e 370 000× abaixo do RMS inicial (37 km).
    assert fit.rms_km < 1e-4
    assert fit.success and fit.cd_a_m is None
    assert fit.nfev == len(fit.history_km) and fit.history_km[0] == fit.rms_initial_km
    r, _ = fit.propagator().propagate(3 * 3600.0)
    r_ref, _ = reference.propagate(3 * 3600.0)
    assert np.linalg.norm(r - r_ref) < 1.0


def test_fit_sentinel_2a_against_sgp4_24h() -> None:
    # Teste rápido do card (A.8): J2–J4, sem arrasto, janela de 24 h a 600 s
    sgp4 = SGP4Propagator(Tle.from_text(SENTINEL_TLE.read_text(encoding="utf-8")))
    fit = fit_initial_state(sgp4, (ZonalAcceleration(),), window_s=A8.quick_window_s,
                            fit_drag=False)  # fmt: skip
    assert fit.rms_km <= A8.quick_rms_crit_km
    assert fit.rms_km < fit.rms_initial_km  # o ajuste melhora o estado do SGP4 em t = 0
    err_24h_km = position_error(fit.propagator(rtol=1e-10, atol=1e-3), sgp4, [86400.0])[0] / 1e3
    assert err_24h_km <= A8.quick_err_24h_crit_km
    assert fit.epoch == sgp4.epoch


def test_fit_rejects_drag_in_models_and_bad_guess() -> None:
    reference = CowellPropagator((R0, V0), EPOCH, (J2Acceleration(),))
    with pytest.raises(ValueError, match="arrasto"):
        fit_initial_state(reference, (J2Acceleration(), DragAcceleration()))
    with pytest.raises(ValueError, match="positivo"):
        fit_initial_state(reference, (J2Acceleration(),), cd_a_m0=0.0)
