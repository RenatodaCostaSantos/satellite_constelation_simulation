"""Testes do módulo de tempo (S1-06)."""

import math
from datetime import UTC, datetime, timedelta, timezone

import numpy as np
import pytest

from satsim.astrotime import datetime_to_jd, gmst, jd_from_epoch, jd_to_datetime

TWO_PI = 2.0 * math.pi


def _wrap(x: np.ndarray) -> np.ndarray:
    """Reduz ângulos a (−π, π]."""
    return (np.asarray(x) + math.pi) % TWO_PI - math.pi


# --- Data Juliana ------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("dt", "jd"),
    [
        # Definições: J2000.0, época Unix e meia-noite de 01/01/2000 (Vallado, cap. 3)
        (datetime(2000, 1, 1, 12, 0, tzinfo=UTC), 2451545.0),
        (datetime(1970, 1, 1, 0, 0, tzinfo=UTC), 2440587.5),
        (datetime(2000, 1, 1, 0, 0, tzinfo=UTC), 2451544.5),
    ],
)
def test_jd_reference_dates(dt: datetime, jd: float) -> None:
    assert datetime_to_jd(dt) == pytest.approx(jd, abs=1e-9)


def test_jd_other_timezone() -> None:
    brt = timezone(timedelta(hours=-3))
    local = datetime(2000, 1, 1, 9, 0, tzinfo=brt)  # = 12:00 UTC
    assert datetime_to_jd(local) == pytest.approx(2451545.0, abs=1e-9)
    assert jd_to_datetime(datetime_to_jd(local)) == local  # mesmo instante


def test_naive_datetime_rejected() -> None:
    with pytest.raises(ValueError):
        datetime_to_jd(datetime(2000, 1, 1, 12, 0))


def test_jd_to_datetime_is_utc_aware() -> None:
    dt = jd_to_datetime(2451545.0)
    assert dt.tzinfo is not None and dt.utcoffset() == timedelta(0)
    assert dt == datetime(2000, 1, 1, 12, 0, tzinfo=UTC)


def test_roundtrip_datetime_jd() -> None:
    rng = np.random.default_rng(42)
    start = datetime(1990, 1, 1, tzinfo=UTC)
    span_us = int((datetime(2040, 1, 1, tzinfo=UTC) - start) / timedelta(microseconds=1))
    for us in rng.integers(0, span_us, 500):
        dt = start + timedelta(microseconds=int(us))
        back = jd_to_datetime(datetime_to_jd(dt))
        assert abs((back - dt).total_seconds()) < 100e-6


def test_jd_from_epoch_vectorized() -> None:
    jd0 = 2451545.0
    t = np.array([0.0, 43200.0, 86400.0])
    np.testing.assert_allclose(jd_from_epoch(jd0, t), [jd0, jd0 + 0.5, jd0 + 1.0], atol=1e-12)
    assert isinstance(jd_from_epoch(jd0, 3600.0), float)


# --- GMST --------------------------------------------------------------------------------------


def test_gmst_j2000() -> None:
    # Vallado eq. 3-47 em T = 0: 67310,54841 s → 280,46061837°
    assert math.degrees(gmst(2451545.0)) == pytest.approx(280.46061837, abs=1e-6)


def test_gmst_vallado_example_3_5() -> None:
    # Vallado, Exemplo 3-5: 20/08/1992 12:14:00 UT1 → JD = 2448855,009722; GMST ≈ 152,5788°
    jd = datetime_to_jd(datetime(1992, 8, 20, 12, 14, 0, tzinfo=UTC))
    assert jd == pytest.approx(2448855.009722, abs=1e-6)
    assert math.degrees(gmst(jd)) == pytest.approx(152.5788, abs=1e-3)


def test_gmst_daily_rate() -> None:
    # Excesso diário sobre 360° = movimento aparente do Sol: 0,98564736629°/dia (prompt S1-06)
    for jd in (2448855.0, 2451545.0, 2461313.0):
        step = (gmst(jd + 1.0) - gmst(jd)) % TWO_PI
        assert math.degrees(step) == pytest.approx(0.98564736629, abs=1e-6)


def test_gmst_range_and_vectorization() -> None:
    jd = np.linspace(2447892.5, 2466154.5, 1000).reshape(10, 100)
    theta = gmst(jd)
    assert isinstance(theta, np.ndarray) and theta.shape == jd.shape
    assert np.all((theta >= 0.0) & (theta < TWO_PI))
    assert isinstance(gmst(2451545.0), float)


def test_gmst_matches_skyfield() -> None:
    skyfield_api = pytest.importorskip("skyfield.api")
    ts = skyfield_api.load.timescale(builtin=True)
    rng = np.random.default_rng(2026)
    # 200 datas aleatórias entre 1990-01-01 e 2040-01-01; UTC tratado como UT1 nos dois lados
    jd_1990 = datetime_to_jd(datetime(1990, 1, 1, tzinfo=UTC))
    jd_2040 = datetime_to_jd(datetime(2040, 1, 1, tzinfo=UTC))
    jd = rng.uniform(jd_1990, jd_2040, 200)
    theta_sky = ts.ut1_jd(jd).gmst * (TWO_PI / 24.0)  # horas → rad
    diff_s = np.abs(_wrap(gmst(jd) - theta_sky)) * 86400.0 / TWO_PI  # segundos de tempo
    assert np.max(diff_s) < 0.1
