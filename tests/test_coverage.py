"""Cobertura diária e revisita (card S5-03).

Referências: Apêndice A.5 do prompt da Semana 5 (Sem4 eq. 45 e experimento de cobertura) e a
referência própria dos testes (COVERAGE_TEST_S5), em tests/golden/reference_values.py.

Custo: a cobertura de 1 dia usa a constelação de referência do conftest (Cowell completo,
integrado uma vez por sessão) na grade de 1° do A.5. As janelas de 3 e 6 dias usam a mesma
constelação no Cowell só J2–J4 (≈ 3,6 s de integração), em grade de 2°: em 3 dias o arrasto, o
Sol e a Lua mudam a cobertura em menos de 0,01 h.
"""

import math
from datetime import datetime

import numpy as np
import pytest
from golden.reference_values import CONSTELLATION_A1, COVERAGE_A5, COVERAGE_TEST_S5

from satsim.constants import R_EARTH, SECONDS_PER_DAY
from satsim.constellation import Constellation, sso_constellation
from satsim.coverage import (
    Access,
    TargetGrid,
    access_times,
    coverage_stats,
    latitude_profile,
    roll_for_half_grid,
    swath_half_angle,
)
from satsim.forces import ZonalAcceleration

A5, REF = COVERAGE_A5, COVERAGE_TEST_S5
EPOCH = datetime.fromisoformat(CONSTELLATION_A1.epoch_utc)
DAY = SECONDS_PER_DAY
ROLL = math.radians(A5.roll_coverage_deg)


@pytest.fixture(scope="module")
def zonal_constellation() -> Constellation:
    """A constelação do A.1 no Cowell só J2–J4: janelas longas a baixo custo."""
    return sso_constellation(
        EPOCH, (ZonalAcceleration(),), i=math.radians(CONSTELLATION_A1.i_calibrated_deg)
    )


@pytest.fixture(scope="module")
def grid1() -> TargetGrid:
    return TargetGrid(1.0, 1.0)


@pytest.fixture(scope="module")
def grid2() -> TargetGrid:
    return TargetGrid(REF.grid_deg, REF.grid_deg)


# --- Eq. 45 e a rolagem de meia grade -----------------------------------------------------------


@pytest.mark.parametrize("k", range(4))
def test_swath_eq45(k: int) -> None:
    lam = swath_half_angle(A5.a_m, math.radians(A5.roll_deg[k]))
    assert math.degrees(lam) == pytest.approx(A5.swath_half_angle_deg[k], abs=1e-4)
    assert lam * R_EARTH / 1e3 == pytest.approx(A5.swath_km_sem4[k], abs=A5.swath_km_tol)
    assert lam * R_EARTH / 1e3 == pytest.approx(A5.swath_km[k], abs=0.05)


def test_swath_vectorized_symmetric_and_errors() -> None:
    lam = swath_half_angle(np.array([A5.a_m, A5.a_m]), np.radians([20.0, -20.0]))
    assert lam.shape == (2,) and lam[0] == pytest.approx(lam[1])
    assert swath_half_angle(A5.a_m, 0.0) == 0.0
    with pytest.raises(ValueError, match="horizonte"):
        swath_half_angle(A5.a_m, math.radians(80.0))
    with pytest.raises(ValueError, match="superfície"):
        swath_half_angle(R_EARTH, 0.1)


def test_roll_for_half_grid() -> None:
    eta = roll_for_half_grid(A5.a_m, 91)
    assert math.degrees(eta) == pytest.approx(A5.roll_half_grid_deg, abs=A5.roll_half_grid_tol_deg)
    assert swath_half_angle(A5.a_m, eta) == pytest.approx(math.pi / 91, abs=1e-12)
    with pytest.raises(ValueError):
        roll_for_half_grid(A5.a_m, 0)
    with pytest.raises(ValueError, match="horizonte"):
        roll_for_half_grid(A5.a_m, 2)  # meia grade de 90°: impossível


# --- TargetGrid e coverage_stats (casos sintéticos) ---------------------------------------------


def test_target_grid(grid1) -> None:
    assert grid1.band_lat_deg.size == 161 and grid1.band_lon_deg.size == 360
    assert grid1.n_targets == 161 * 360
    assert grid1.band_lat_deg[0] == -80.0 and grid1.band_lat_deg[-1] == 80.0
    assert grid1.band_lon_deg[0] == -180.0 and grid1.band_lon_deg[-1] == 179.0
    assert np.allclose(np.linalg.norm(grid1.unit, axis=1), 1.0)
    # ordenado por faixa: as 360 primeiras são a latitude −80°
    assert np.all(grid1.lat_deg[:360] == -80.0)
    assert np.degrees(np.arcsin(grid1.unit[400, 2])) == pytest.approx(grid1.lat_deg[400])
    with pytest.raises(ValueError):
        TargetGrid(0.0, 1.0)
    with pytest.raises(ValueError):
        TargetGrid(1.0, 1.0, lat_max_deg=95.0)


def test_coverage_stats_synthetic() -> None:
    # alvo 0: passagem de 0,5 h às 2 h (amostras a cada 10 min) e outra às 20 h; alvo 1: um acesso
    # às 10 h; alvo 2: nunca visto. Janela de 24 h.
    h = 3600.0
    t0 = [2 * h, 2 * h + 600, 2 * h + 1200, 2 * h + 1800 - 1, 20 * h]
    access = Access(target_index=np.array([0] * 5 + [1]), t_s=np.array([*t0, 10 * h]), n_targets=3)
    frac, revisit = coverage_stats(access, 24 * h)
    assert frac == pytest.approx(2 / 3)
    assert revisit / h == pytest.approx([18.0, 14.0, 24.0])  # 2 h → 20 h; 10 h → 24 h; nunca
    _, blackout = coverage_stats(access, 24 * h, interval="blackout")
    assert blackout[0] / h == pytest.approx(20.0 - (2.5 - 1 / 3600))  # fim da 1ª → 20 h
    assert blackout[1:] / h == pytest.approx([14.0, 24.0])
    # gap_merge pequeno: as amostras de 10 min viram passagens separadas
    _, revisit_split = coverage_stats(access, 24 * h, gap_merge_s=300.0)
    assert revisit_split[0] / h == pytest.approx(20.0 - (2.5 - 1 / 3600))
    assert access.times(0) == pytest.approx(t0)
    with pytest.raises(ValueError):
        coverage_stats(access, 24 * h, interval="outro")
    empty = Access(np.zeros(0, dtype=np.int64), np.zeros(0), 4)
    assert coverage_stats(empty, 100.0)[0] == 0.0
    assert np.all(coverage_stats(empty, 100.0)[1] == 100.0)


def test_latitude_profile(grid2) -> None:
    values = np.repeat(np.arange(grid2.band_lat_deg.size, dtype=float), grid2.band_lon_deg.size)
    lat, mean = latitude_profile(grid2, values)
    assert lat == pytest.approx(grid2.band_lat_deg) and mean == pytest.approx(np.arange(lat.size))
    assert latitude_profile(grid2, values, "max")[1] == pytest.approx(mean)
    with pytest.raises(ValueError):
        latitude_profile(grid2, values, "min")


# --- Cobertura diária (Cowell completo, grade 1°, A.5) ------------------------------------------


def test_daily_coverage_six_satellites(reference_constellation, grid1) -> None:
    """P = 6, η = 23,2°, 1 dia: 100% dos alvos em |lat| ≤ 80° vistos ao menos uma vez."""
    access = access_times(reference_constellation, grid1, EPOCH, DAY, ROLL)
    frac, revisit = coverage_stats(access, DAY)
    assert frac == A5.fraction_p6_1d
    assert np.all(revisit < DAY)


def test_daily_coverage_one_satellite(reference_constellation, grid1) -> None:
    # P = 1 no mesmo dia: A.5 48,48% ± 2% (medido 48,82%)
    s1 = Constellation(EPOCH, [reference_constellation[0]])
    frac, _ = coverage_stats(access_times(s1, grid1, EPOCH, DAY, ROLL), DAY)
    assert frac == pytest.approx(A5.fraction_p1_1d, abs=A5.fraction_p1_1d_tol)


def test_one_satellite_needs_six_days(zonal_constellation, grid2) -> None:
    s1 = Constellation(EPOCH, [zonal_constellation[0]])
    frac_6d, _ = coverage_stats(access_times(s1, grid2, EPOCH, 6 * DAY, ROLL), 6 * DAY)
    frac_5d, _ = coverage_stats(access_times(s1, grid2, EPOCH, 5 * DAY, ROLL), 5 * DAY)
    assert frac_6d == A5.fraction_p1_6d
    assert frac_5d < 1.0  # o ciclo de 6 dias é necessário


# --- Revisita em 3 dias (Cowell J2–J4, grade 2°) ------------------------------------------------


@pytest.fixture(scope="module")
def access_3d(zonal_constellation, grid2) -> tuple[Access, Access]:
    dark = access_times(zonal_constellation, grid2, EPOCH, 3 * DAY, ROLL)
    lit = access_times(
        zonal_constellation, grid2, EPOCH, 3 * DAY, ROLL,
        min_sun_elevation_deg=A5.min_sun_elevation_deg,
    )  # fmt: skip
    return dark, lit


def test_revisit_three_days(access_3d, grid2) -> None:
    """P = 6: revisita ≤ 15 h em |lat| ≤ 60° (A.5: 13,97 h; uma passagem ascendente e uma
    descendente por dia)."""
    dark, _ = access_3d
    frac, revisit = coverage_stats(dark, 3 * DAY, gap_merge_s=A5.gap_merge_s)
    lat60 = np.abs(grid2.lat_deg) <= 60.0
    assert frac == 1.0
    assert revisit[lat60].max() / 3600 <= A5.revisit_p6_3d_lat60_limit_h
    assert revisit[lat60].max() / 3600 == pytest.approx(
        REF.revisit_p6_3d_lat60_h, abs=REF.regression_tol_h
    )
    lat80_h = revisit.max() / 3600
    assert lat80_h == pytest.approx(REF.revisit_p6_3d_lat80_h, abs=REF.regression_tol_h)
    _, blackout = coverage_stats(dark, 3 * DAY, interval="blackout")
    assert blackout[lat60].max() / 3600 == pytest.approx(
        REF.blackout_p6_3d_lat60_h, abs=REF.regression_tol_h
    )
    assert np.all(blackout <= revisit)


def test_revisit_three_days_daylight_only(access_3d, grid2) -> None:
    """Só com Sol ≥ 10° no alvo: ≤ 25 h em |lat| ≤ 60° (A.5: 24,26 h); acima de ~70° N, em
    outubro, o Sol fica baixo demais e a fração vista cai para ≈ 92,5%."""
    dark, lit = access_3d
    frac, revisit = coverage_stats(lit, 3 * DAY)
    lat60 = np.abs(grid2.lat_deg) <= 60.0
    assert revisit[lat60].max() / 3600 <= A5.revisit_p6_3d_sun10_lat60_limit_h
    assert revisit[lat60].max() / 3600 == pytest.approx(
        REF.revisit_p6_3d_sun10_lat60_h, abs=REF.regression_tol_h
    )
    assert frac == pytest.approx(REF.fraction_p6_3d_sun10, abs=1e-4)
    assert frac == pytest.approx(A5.fraction_p6_3d_sun10, abs=0.005)  # A.5 em grade de 1°
    # os alvos nunca vistos com luz são os do Ártico (noite polar chegando)
    unseen = revisit == 3 * DAY
    assert np.all(grid2.lat_deg[unseen] > 60.0)
    # a restrição só remove acessos: com luz, o tempo sem acesso nunca diminui (a revisita, de
    # início a início, não é monótona: tirar um acesso pode partir uma passagem em duas)
    blackout_lit = coverage_stats(lit, 3 * DAY, interval="blackout")[1]
    blackout_dark = coverage_stats(dark, 3 * DAY, interval="blackout")[1]
    assert np.all(blackout_lit >= blackout_dark)


# --- Sanidade -----------------------------------------------------------------------------------


def test_zero_roll_gives_no_access(reference_constellation) -> None:
    grid = TargetGrid(4.0, 4.0)
    access = access_times(reference_constellation, grid, EPOCH, 3600.0, 0.0)
    assert access.t_s.size == 0
    assert coverage_stats(access, 3600.0)[0] == 0.0


def test_doubling_roll_never_reduces_coverage(reference_constellation) -> None:
    grid = TargetGrid(4.0, 4.0)
    s1 = Constellation(EPOCH, [reference_constellation[0]])
    small = access_times(s1, grid, EPOCH, DAY, math.radians(10.0))
    large = access_times(s1, grid, EPOCH, DAY, math.radians(20.0))
    f_small, r_small = coverage_stats(small, DAY, interval="blackout")
    f_large, r_large = coverage_stats(large, DAY, interval="blackout")
    assert f_large >= f_small
    seen_small = set(zip(small.target_index.tolist(), small.t_s.tolist(), strict=True))
    seen_large = set(zip(large.target_index.tolist(), large.t_s.tolist(), strict=True))
    assert seen_small <= seen_large  # todo acesso com 10° continua com 20°
    assert np.all(r_large <= r_small)  # tempo sem acesso: monótono em acessos


def test_access_times_rejects_bad_inputs(reference_constellation) -> None:
    grid = TargetGrid(10.0, 10.0)
    with pytest.raises(ValueError):
        access_times(reference_constellation, grid, EPOCH, DAY, ROLL, step_s=0.0)
    with pytest.raises(ValueError):
        access_times(reference_constellation, grid, EPOCH, -1.0, ROLL)
