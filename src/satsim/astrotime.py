"""Tempo astronômico: datetime UTC ↔ Data Juliana (JD) ↔ tempo sideral médio de Greenwich (θG).

O nome ``astrotime`` (e não ``time``) evita sombrear o módulo ``time`` da biblioteca padrão
(ADR 0001).

O GMST θG é o ângulo entre o equinócio vernal (eixo x do ECI) e o meridiano de Greenwich; é a
rotação em torno de z que leva ECI → ECEF (semana 2).

Simplificações documentadas
---------------------------
(a) **UTC é tratado como UT1.** A fórmula do GMST pede UT1 (tempo da rotação da Terra); usamos
    UTC direto. Como |UT1 − UTC| < 0,9 s, o erro máximo é ≈ 0,9 s de rotação, ≈ 0,4 km no
    equador.
(b) **Segundos intercalares são ignorados**: ``datetime`` não os representa, e
    ``timestamp()`` conta 86 400 s por dia. Intervalos que atravessam um segundo intercalar
    ficam 1 s mais curtos que o real.
(c) **Taxa do GMST ≠ OMEGA_EARTH.** A derivada da fórmula abaixo vale ≈ 7,2921158553e-5 rad/s,
    ligeiramente diferente da constante ``OMEGA_EARTH`` = 7,292115e-5 rad/s usada nos
    documentos (≈ 0,15"/dia, ≈ 5 m/dia no equador). Para converter ECI → ECEF vale a fórmula
    completa do GMST, não ``OMEGA_EARTH · t``.
"""

import math
from datetime import UTC, datetime, timedelta

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import SECONDS_PER_DAY

FloatOrArray = float | NDArray[np.float64]

JD_UNIX_EPOCH: float = 2440587.5
"""Data Juliana de 1970-01-01 00:00:00 UTC (época do Unix timestamp)."""

JD_J2000: float = 2451545.0
"""Data Juliana de J2000.0 = 2000-01-01 12:00:00 (TT na definição; aqui tratado como UT1)."""

DAYS_PER_JULIAN_CENTURY: float = 36525.0

_UNIX_EPOCH = datetime(1970, 1, 1, tzinfo=UTC)


def _as_scalar_if_0d(x: NDArray[np.float64]) -> FloatOrArray:
    return float(x) if x.ndim == 0 else x


def datetime_to_jd(dt: datetime) -> float:
    """``datetime`` tz-aware → Data Juliana (UTC).

    JD = 2440587,5 + timestamp/86400, exato para qualquer fuso: ``timestamp()`` já desconta o
    offset, de modo que 09:00−03:00 e 12:00 UTC dão o mesmo JD.

    Raises:
        ValueError: se ``dt`` for naive (sem fuso horário).
    """
    if dt.tzinfo is None or dt.utcoffset() is None:
        raise ValueError("datetime deve ser tz-aware (ex.: tzinfo=timezone.utc)")
    return JD_UNIX_EPOCH + dt.timestamp() / SECONDS_PER_DAY


def jd_to_datetime(jd: float) -> datetime:
    """Data Juliana → ``datetime`` UTC tz-aware, arredondado a microssegundos.

    A precisão é limitada pelo float do JD: ≈ 2,5e6 dias com ε ≈ 2,2e-16 dá ≈ 40 µs.
    """
    seconds = (float(jd) - JD_UNIX_EPOCH) * SECONDS_PER_DAY
    return _UNIX_EPOCH + timedelta(microseconds=round(seconds * 1e6))


def jd_from_epoch(epoch_jd: float, t: ArrayLike) -> FloatOrArray:
    """JD dos instantes t [s] desde a época: epoch_jd + t/86400, vetorizado sobre t."""
    return _as_scalar_if_0d(epoch_jd + np.asarray(t, dtype=float) / SECONDS_PER_DAY)


def gmst(jd_ut1: ArrayLike) -> FloatOrArray:
    """Tempo sideral médio de Greenwich θG [rad] em [0, 2π), vetorizado.

    Fórmula IAU-82 (Vallado, eq. 3-47), com T em séculos julianos desde J2000:

        T = (JD − 2451545,0) / 36525
        θ [s] = 67310,54841 + (876600·3600 + 8640184,812866)·T + 0,093104·T² − 6,2e-6·T³

    e θ [rad] = (θ mod 86400) · 2π/86400. Recebe JD em UT1; ver nota (a) do módulo.

    Args:
        jd_ut1: Data Juliana (UT1 ≈ UTC), escalar ou array.

    Returns:
        θG [rad] com o shape da entrada (float para escalar).
    """
    t_cent = (np.asarray(jd_ut1, dtype=float) - JD_J2000) / DAYS_PER_JULIAN_CENTURY
    theta_s = (
        67310.54841
        + (876600.0 * 3600.0 + 8640184.812866) * t_cent
        + 0.093104 * t_cent**2
        - 6.2e-6 * t_cent**3
    )
    return _as_scalar_if_0d(np.mod(theta_s, SECONDS_PER_DAY) * (2.0 * math.pi / SECONDS_PER_DAY))
