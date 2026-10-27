"""Efemérides analíticas do Sol e da Lua no ECI da data e grandezas de iluminação da SSO.

Sem arquivos externos: fórmulas de baixa precisão, vetorizadas sobre a Data Juliana. O quadro de
saída é o **equador e equinócio da data**, o mesmo "ECI" do projeto (= TEME, ADR 0001). Usar
posições referidas ao equinócio J2000 erraria a direção do Sol em ≈ 0,37° em 2026 (precessão
de 1,3972°/século; Apêndice A.3 do prompt da Semana 4).

Escala de tempo: a JD de entrada é UTC e é usada como TT (diferença de ~69 s): erro de ~0,01° na
Lua e de ~0,0008° no Sol, documentado e absorvido nas tolerâncias.

Precisão em 2026 (a cada 6 h, 1460 instantes, contra o JPL DE421; Apêndice A.3):

- Sol (Astronomical Almanac, baixa precisão): direção máx. 0,0140°, |Δr|/r ≤ 8,5e-5. A fórmula do
  Sol de Montenbruck & Gill §3.3.2 foi rejeitada: ela fixa o periélio (0,32°/século) e erra
  0,094° em 2026.
- Lua (Montenbruck & Gill 2000, §3.3.2, com +1,3972°·T na longitude para o equinócio da data):
  direção máx. 0,085°, |Δr| ≤ 445 km.

Iluminação (Sem4 §5): o LTAN (hora solar local do nó ascendente, eq. 28) e o ângulo β entre o
plano da órbita e a direção do Sol (eq. 29).
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.astrotime import DAYS_PER_JULIAN_CENTURY, JD_J2000
from satsim.constants import AU, DEG2RAD

_ARCSEC = DEG2RAD / 3600.0
PRECESSION_DEG_PER_CENTURY = 1.3972
"""Precessão geral em longitude [°/século] (Montenbruck & Gill 2000, §3.3.2)."""


def _stack_xyz(x: NDArray, y: NDArray, z: NDArray) -> NDArray[np.float64]:
    return np.stack([x, y, z], axis=-1)


def sun_position_eci(jd: ArrayLike) -> NDArray[np.float64]:
    """Posição geocêntrica do Sol no ECI da data [m] (Astronomical Almanac, baixa precisão).

    Fórmulas do Apêndice A.3 do prompt da Semana 4, com n = JD − 2451545,0:

        L = 280,460° + 0,9856474°·n          g = 357,528° + 0,9856003°·n
        λ = L + 1,915°·sin g + 0,020°·sin 2g    (latitude eclíptica β = 0)
        ε = 23,439° − 0,0000004°·n
        R = (1,00014 − 0,01671·cos g − 0,00014·cos 2g) UA
        r☉ = R·(cos λ, cos ε·sin λ, sin ε·sin λ)

    λ e ε já são da data (o periélio avança com g), por isso não há correção de precessão.

    Args:
        jd: Data Juliana (UTC usado como TT), escalar ou shape (N,).

    Returns:
        Posição [m], shape (3,) ou (N, 3).
    """
    n = np.asarray(jd, dtype=float) - JD_J2000
    mean_lon = (280.460 + 0.9856474 * n) * DEG2RAD
    g = (357.528 + 0.9856003 * n) * DEG2RAD
    lam = mean_lon + (1.915 * np.sin(g) + 0.020 * np.sin(2.0 * g)) * DEG2RAD
    eps = (23.439 - 0.0000004 * n) * DEG2RAD
    dist = (1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2.0 * g)) * AU
    return dist[..., None] * _stack_xyz(
        np.cos(lam), np.cos(eps) * np.sin(lam), np.sin(eps) * np.sin(lam)
    )


def moon_position_eci(jd: ArrayLike) -> NDArray[np.float64]:
    """Posição geocêntrica da Lua no ECI da data [m] (Montenbruck & Gill 2000, §3.3.2).

    Série do Apêndice A.3 do prompt da Semana 4, com T = (JD − 2451545,0)/36525. As longitudes
    da série de M&G referem-se ao equinócio J2000 (L0 inclui −1,3972°·T); para o equinócio da
    data soma-se de volta a precessão, λ ← λ + 1,3972°·T, e usa-se a obliquidade da data
    ε = 23,43929111° − 0,0130042°·T. Sem essa correção o erro em 2026 seria 0,46°.

    Args:
        jd: Data Juliana (UTC usado como TT), escalar ou shape (N,).

    Returns:
        Posição [m], shape (3,) ou (N, 3).
    """
    t = (np.asarray(jd, dtype=float) - JD_J2000) / DAYS_PER_JULIAN_CENTURY
    l0 = (218.31617 + 481267.88088 * t - PRECESSION_DEG_PER_CENTURY * t) * DEG2RAD
    l = (134.96292 + 477198.86753 * t) * DEG2RAD  # anomalia média da Lua  # noqa: E741
    lp = (357.52543 + 35999.04944 * t) * DEG2RAD  # anomalia média do Sol
    f = (93.27283 + 483202.01873 * t) * DEG2RAD  # distância média ao nó
    d = (297.85027 + 445267.11135 * t) * DEG2RAD  # elongação média
    sin = np.sin
    dlam = (
        22640 * sin(l)
        + 769 * sin(2 * l)
        - 4586 * sin(l - 2 * d)
        + 2370 * sin(2 * d)
        - 668 * sin(lp)
        - 412 * sin(2 * f)
        - 212 * sin(2 * l - 2 * d)
        - 206 * sin(l + lp - 2 * d)
        + 192 * sin(l + 2 * d)
        - 165 * sin(lp - 2 * d)
        + 148 * sin(l - lp)
        - 125 * sin(d)
        - 110 * sin(l + lp)
        - 55 * sin(2 * f - 2 * d)
    ) * _ARCSEC  # λ − L0
    lam = l0 + dlam
    beta = (
        18520 * sin(f + dlam + (412 * sin(2 * f) + 541 * sin(lp)) * _ARCSEC)
        - 526 * sin(f - 2 * d)
        + 44 * sin(l + f - 2 * d)
        - 31 * sin(-l + f - 2 * d)
        - 25 * sin(-2 * l + f)
        - 23 * sin(lp + f - 2 * d)
        + 21 * sin(-l + f)
        + 11 * sin(-lp + f - 2 * d)
    ) * _ARCSEC
    cos = np.cos
    dist_km = (
        385000
        - 20905 * cos(l)
        - 3699 * cos(2 * d - l)
        - 2956 * cos(2 * d)
        - 570 * cos(2 * l)
        + 246 * cos(2 * l - 2 * d)
        - 205 * cos(lp - 2 * d)
        - 171 * cos(l + 2 * d)
        - 152 * cos(l + lp - 2 * d)
    )
    lam = lam + PRECESSION_DEG_PER_CENTURY * t * DEG2RAD  # equinócio da data
    eps = (23.43929111 - 0.0130042 * t) * DEG2RAD
    cb, sb, cl, sl, ce, se = cos(beta), sin(beta), cos(lam), sin(lam), cos(eps), sin(eps)
    return (dist_km * 1e3)[..., None] * _stack_xyz(
        cb * cl, ce * cb * sl - se * sb, se * cb * sl + ce * sb
    )


def sun_direction(jd: ArrayLike) -> NDArray[np.float64]:
    """Versor ŝ do Sol no ECI da data, shape (3,) ou (N, 3)."""
    s = sun_position_eci(jd)
    return s / np.linalg.norm(s, axis=-1, keepdims=True)


def sun_ra_dec(jd: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Ascensão reta α☉ em [0, 2π) e declinação δ☉ do Sol [rad], na data, escalar ou (N,)."""
    s = sun_position_eci(jd)
    ra = np.mod(np.arctan2(s[..., 1], s[..., 0]), 2.0 * np.pi)
    dec = np.arcsin(s[..., 2] / np.linalg.norm(s, axis=-1))
    return ra, dec


def ltan_hours(raan: ArrayLike, jd: ArrayLike) -> NDArray[np.float64]:
    """Hora solar local do nó ascendente [h] em [0, 24) (Sem4 eq. 28).

    LTAN = 12 h + (Ω − α☉)/15°: o nó está ao meio-dia quando Ω = α☉ (o plano contém a direção do
    Sol). Usa o tempo solar *médio* aparente pela ascensão reta verdadeira, como na eq. 28.

    Args:
        raan: ascensão reta do nó ascendente Ω [rad], escalar ou (N,).
        jd: Data Juliana, escalar ou (N,) (broadcasting com ``raan``).
    """
    ra_sun, _ = sun_ra_dec(jd)
    hours = 12.0 + (np.asarray(raan, dtype=float) - ra_sun) / (15.0 * DEG2RAD)
    return np.mod(hours, 24.0)


def beta_angle(r: ArrayLike, v: ArrayLike, jd: ArrayLike) -> NDArray[np.float64]:
    """Ângulo β entre o plano da órbita e a direção do Sol [rad] (Sem4 eq. 29).

    sin β = ŝ·ℓ̂, ℓ = r × v (normal ao plano). β > 0 quando o Sol está do lado de +ℓ.

    Args:
        r, v: estado ECI [m, m/s], shape (3,) ou (N, 3).
        jd: Data Juliana, escalar ou (N,).

    Returns:
        β em [−π/2, π/2], shape () ou (N,).
    """
    ell = np.cross(np.asarray(r, dtype=float), np.asarray(v, dtype=float))
    ell_hat = ell / np.linalg.norm(ell, axis=-1, keepdims=True)
    sin_beta = np.sum(sun_direction(jd) * ell_hat, axis=-1)
    return np.arcsin(np.clip(sin_beta, -1.0, 1.0))
