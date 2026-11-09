"""Cobertura e revisita da constelação com apontamento fora do nadir (card S5-03).

A grade de 91 trilhas da órbita 91/6 deixa 440 km entre trilhas vizinhas no equador. Para que
todo ponto seja alcançado todo dia, cada satélite precisa enxergar até meia grade (220 km) para
cada lado do traço, girando a câmera (rolagem η) fora do nadir. O alcance transversal é o ângulo
geocêntrico (Sem4 §8.6, eq. 45)

    λ = arcsin((a/R⊕)·sin η) − η,

entre o subponto e o ponto mais distante visível com rolagem η, num satélite à distância a do
centro da Terra (esfera de raio R⊕). Com η = 23,2° e a = 6 888 km, λ = π/91: meia grade.

Fluxo:

1. ``TargetGrid``: pontos-alvo numa grade regular de latitude e longitude (|lat| ≤ lat_max);
2. ``access_times``: para cada alvo, os instantes amostrados em que o subponto de algum satélite
   está a ângulo geocêntrico ≤ λ do alvo (opcionalmente só com o alvo iluminado: câmera óptica);
3. ``coverage_stats``: fração de alvos vistos ao menos uma vez e, por alvo, o maior intervalo
   sem visita (incluindo as bordas da janela).

Geometria: esfera de raio R⊕ (como no Sem4); o subponto é a direção geocêntrica de r_ECEF, e λ
vem da distância instantânea |r| do satélite (eq. 45 com a = |r|). O Sol é visto do alvo como
direção (sem paralaxe): elevação solar = arcsin(û_alvo · ŝ_ECEF).
"""

import math
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime
from functools import cached_property

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq

from satsim.astrotime import datetime_to_jd, jd_from_epoch
from satsim.constants import R_EARTH
from satsim.ephemeris import sun_direction
from satsim.frames import eci_to_ecef, gmst_at
from satsim.satellite import Satellite

FloatOrArray = float | NDArray[np.float64]


def swath_half_angle(a: ArrayLike, roll: ArrayLike) -> FloatOrArray:
    """Alcance transversal λ [rad]: ângulo geocêntrico visível com rolagem η (Sem4 eq. 45).

    λ = arcsin((a/R⊕)·sin η) − η. Vem do triângulo centro da Terra – satélite – alvo: a linha de
    visada faz o ângulo η com o nadir, e pela lei dos senos o ângulo no alvo, entre a vertical
    local e a visada, é arcsin((a/R⊕)·sin η); λ é o que sobra para fechar os 180°. A rolagem é
    tomada em módulo (a faixa é simétrica em torno do traço).

    Args:
        a: distância do satélite ao centro da Terra [m] (o semieixo, numa órbita circular).
        roll: rolagem η [rad], |η| abaixo do horizonte, arcsin(R⊕/a).

    Returns:
        λ [rad]; float para entradas escalares.

    Raises:
        ValueError: se a visada passar do horizonte ((a/R⊕)·sin|η| > 1) ou a ≤ R⊕.
    """
    a, eta = np.broadcast_arrays(np.asarray(a, dtype=float), np.abs(np.asarray(roll, dtype=float)))
    if np.any(a <= R_EARTH):
        raise ValueError("o satélite deve estar acima da superfície: a > R⊕")
    k = a / R_EARTH * np.sin(eta)
    if np.any(k > 1.0):
        raise ValueError("rolagem além do horizonte: (a/R⊕)·sin η > 1")
    lam = np.arcsin(k) - eta
    return float(lam) if lam.ndim == 0 else lam


def roll_for_half_grid(a: float, n_tracks: int) -> float:
    """Rolagem η [rad] cujo alcance transversal é meia grade, λ = π/N (Sem4 §8.6).

    Com ela, a faixa de um satélite toca a da trilha vizinha: o ponto mais distante de qualquer
    trilha (no meio do vão de 2π/N) fica exatamente na borda. Na órbita 91/6: 23,1957°.

    Args:
        a: distância do satélite ao centro da Terra [m].
        n_tracks: N, número de trilhas da grade.

    Raises:
        ValueError: se N < 1 ou se nem a rolagem até o horizonte alcançar π/N.
    """
    if n_tracks < 1:
        raise ValueError("n_tracks deve ser ≥ 1")
    target = math.pi / n_tracks
    eta_horizon = math.asin(R_EARTH / a)
    if float(swath_half_angle(a, eta_horizon)) < target:
        raise ValueError("nem a rolagem até o horizonte alcança meia grade")
    return float(brentq(lambda eta: float(swath_half_angle(a, eta)) - target, 0.0, eta_horizon))


@dataclass(frozen=True)
class TargetGrid:
    """Pontos-alvo numa grade regular de latitude e longitude, na esfera unitária (ECEF).

    Latitudes de −lat_max a +lat_max (inclusive) a cada ``lat_step_deg``; longitudes de −180° a
    180° (exclusive) a cada ``lon_step_deg``. Os alvos ficam ordenados por faixa de latitude
    (todas as longitudes de uma latitude são contíguas), o que permite vetorizar por faixa.

    Attributes:
        lat_step_deg: passo em latitude [°].
        lon_step_deg: passo em longitude [°].
        lat_max_deg: maior |latitude| [°] (80°: além disso a SSO, que sobe até 82,6°, não cobre
            todas as longitudes).
    """

    lat_step_deg: float
    lon_step_deg: float
    lat_max_deg: float = 80.0

    def __post_init__(self) -> None:
        if self.lat_step_deg <= 0.0 or self.lon_step_deg <= 0.0:
            raise ValueError("os passos da grade devem ser positivos")
        if not 0.0 <= self.lat_max_deg <= 90.0:
            raise ValueError("lat_max_deg deve estar em [0, 90]")

    @cached_property
    def band_lat_deg(self) -> NDArray[np.float64]:
        """Latitudes das faixas [°], shape (B,)."""
        n = math.floor(2.0 * self.lat_max_deg / self.lat_step_deg + 1e-9)
        return -self.lat_max_deg + self.lat_step_deg * np.arange(n + 1)

    @cached_property
    def band_lon_deg(self) -> NDArray[np.float64]:
        """Longitudes de cada faixa [°] em [−180, 180), shape (L,)."""
        n = math.floor(360.0 / self.lon_step_deg - 1e-9) + 1
        return -180.0 + self.lon_step_deg * np.arange(n)

    @cached_property
    def lat_deg(self) -> NDArray[np.float64]:
        """Latitude de cada alvo [°], shape (M,) com M = B·L."""
        return np.repeat(self.band_lat_deg, self.band_lon_deg.size)

    @cached_property
    def lon_deg(self) -> NDArray[np.float64]:
        """Longitude de cada alvo [°], shape (M,)."""
        return np.tile(self.band_lon_deg, self.band_lat_deg.size)

    @cached_property
    def unit(self) -> NDArray[np.float64]:
        """Versores ECEF dos alvos na esfera, shape (M, 3)."""
        lat, lon = np.radians(self.lat_deg), np.radians(self.lon_deg)
        return np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])

    @property
    def n_targets(self) -> int:
        """M, número de alvos."""
        return int(self.lat_deg.size)


@dataclass(frozen=True)
class Access:
    """Acessos amostrados: pares (alvo, instante) em que um satélite alcança o alvo.

    Attributes:
        target_index: índice do alvo em ``TargetGrid`` (K,), ordenado por alvo e depois por t.
        t_s: instante [s desde o início da janela] (K,).
        n_targets: M, número de alvos da grade (inclui os nunca vistos).
    """

    target_index: NDArray[np.int64]
    t_s: NDArray[np.float64]
    n_targets: int

    def times(self, target: int) -> NDArray[np.float64]:
        """Instantes [s] dos acessos do alvo ``target``, em ordem."""
        lo, hi = np.searchsorted(self.target_index, [target, target + 1])
        return self.t_s[lo:hi]


def access_times(
    constellation: Iterable[Satellite],
    grid: TargetGrid,
    start: datetime,
    duration_s: float,
    roll: float,
    step_s: float = 10.0,
    min_sun_elevation_deg: float | None = None,
) -> Access:
    """Instantes em que algum satélite alcança cada alvo, com rolagem até ``roll`` (S5-03).

    Um alvo é alcançado numa amostra se o ângulo geocêntrico entre o alvo e o subponto do
    satélite for ≤ λ = ``swath_half_angle(|r|, roll)`` (Sem4 eq. 45). Com
    ``min_sun_elevation_deg``, só conta se o Sol estiver a essa elevação ou mais no alvo
    (imageamento óptico).

    Vetorizado por faixa de latitude: para cada faixa, só entram as amostras cujo subponto está a
    menos de λ em latitude (o ângulo geocêntrico nunca é menor que a diferença de latitudes), e o
    teste é um produto escalar (alvos da faixa × amostras candidatas). Nunca se forma o produto
    completo alvos × amostras.

    Args:
        constellation: os satélites (``Constellation`` ou qualquer iterável de ``Satellite``).
        grid: os alvos.
        start: início da janela (tz-aware).
        duration_s: duração [s]; amostras em 0, step_s, …, até duration_s (inclusive).
        roll: rolagem máxima η [rad] (0 → λ = 0: nenhum acesso).
        step_s: passo da amostragem [s].
        min_sun_elevation_deg: elevação solar mínima no alvo [°], ou ``None`` (sem restrição).

    Returns:
        ``Access`` com os pares (alvo, instante), ordenados por alvo e instante.
    """
    if step_s <= 0.0 or duration_s < 0.0:
        raise ValueError("step_s deve ser positivo e duration_s não negativo")
    t = np.arange(0.0, duration_s + 0.5 * step_s, step_s)
    theta = gmst_at(start, t)
    units, lats, lams, times = [], [], [], []
    for sat in constellation:
        r_eci, _ = sat.propagate_utc(start, t)
        r_ecef = eci_to_ecef(r_eci, theta)
        radius = np.linalg.norm(r_ecef, axis=-1)
        u = r_ecef / radius[:, None]
        units.append(u)
        lats.append(np.arcsin(np.clip(u[:, 2], -1.0, 1.0)))
        lams.append(np.asarray(swath_half_angle(radius, roll)))
        times.append(t)
    u_all = np.concatenate(units)
    lat_all = np.concatenate(lats)
    lam_all = np.concatenate(lams)
    t_all = np.concatenate(times)
    cos_lam_all = np.cos(lam_all)
    sun_all = None
    if min_sun_elevation_deg is not None:
        sun_ecef = eci_to_ecef(sun_direction(jd_from_epoch(datetime_to_jd(start), t)), theta)
        sun_all = np.tile(sun_ecef, (len(units), 1))
        sin_el_min = math.sin(math.radians(min_sun_elevation_deg))

    n_lon = grid.band_lon_deg.size
    idx_out: list[NDArray[np.int64]] = []
    t_out: list[NDArray[np.float64]] = []
    for b, lat_b in enumerate(np.radians(grid.band_lat_deg)):
        cand = np.flatnonzero((np.abs(lat_all - lat_b) <= lam_all) & (lam_all > 0.0))
        if cand.size == 0:
            continue
        targets = grid.unit[b * n_lon : (b + 1) * n_lon]
        hit = targets @ u_all[cand].T >= cos_lam_all[cand]
        if sun_all is not None:
            hit &= targets @ sun_all[cand].T >= sin_el_min
        ti, ki = np.nonzero(hit)
        idx_out.append(b * n_lon + ti)
        t_out.append(t_all[cand[ki]])

    if idx_out:
        idx = np.concatenate(idx_out).astype(np.int64)
        tt = np.concatenate(t_out)
    else:
        idx, tt = np.zeros(0, dtype=np.int64), np.zeros(0)
    order = np.lexsort((tt, idx))
    return Access(target_index=idx[order], t_s=tt[order], n_targets=grid.n_targets)


def coverage_stats(
    access: Access, duration_s: float, gap_merge_s: float = 1800.0, interval: str = "revisit"
) -> tuple[float, NDArray[np.float64]]:
    """Fração de alvos vistos e o maior intervalo entre visitas de cada alvo (S5-03).

    Acessos de um alvo separados por menos de ``gap_merge_s`` contam como uma mesma passagem
    (amostras seguidas de um sobrevoo, ou satélites vizinhos do trem, 15,8 min um atrás do
    outro). O intervalo máximo inclui as bordas da janela e depende de ``interval``:

    - ``"revisit"`` (padrão): tempo de revisita, entre os **inícios** de passagens sucessivas;
      bordas: do início da janela à 1ª passagem e da última passagem ao fim. É a definição do
      Apêndice A.5 do prompt da Semana 5 (13,97 h em |lat| ≤ 60°);
    - ``"blackout"``: tempo sem acesso, do **último** acesso de uma passagem ao **primeiro** da
      seguinte; bordas: até o 1º acesso e do último acesso ao fim (13,69 h no mesmo caso).

    A diferença é a duração da passagem: onde vários satélites do trem veem o mesmo alvo em
    sequência (latitudes altas), uma passagem dura até horas.

    Um alvo nunca visto tem intervalo máximo igual a ``duration_s``.

    Args:
        access: saída de ``access_times``.
        duration_s: duração da janela [s].
        gap_merge_s: separação mínima entre passagens distintas [s].
        interval: ``"revisit"`` ou ``"blackout"``.

    Returns:
        (fração de alvos vistos ao menos uma vez, intervalo máximo por alvo [s] com shape (M,)).
    """
    if interval not in ("revisit", "blackout"):
        raise ValueError("interval deve ser 'revisit' ou 'blackout'")
    idx = np.asarray(access.target_index)
    t = np.asarray(access.t_s, dtype=float)
    max_gap = np.full(access.n_targets, float(duration_s))
    if idx.size == 0:
        return 0.0, max_gap
    # Ordena por alvo e instante (garantido por access_times, refeito por segurança).
    order = np.lexsort((t, idx))
    idx, t = idx[order], t[order]
    # Separação para o acesso anterior do mesmo alvo (inf no 1º acesso de cada alvo).
    first = np.r_[True, idx[1:] != idx[:-1]]
    since_prev = np.where(first, np.inf, np.r_[0.0, np.diff(t)])
    # Início de passagem: 1º acesso do alvo ou separado do anterior por ≥ gap_merge_s.
    new_pass = since_prev >= gap_merge_s
    if interval == "revisit":
        # Só os inícios de passagem: intervalo de cada um ao início anterior (ou ao começo da
        # janela) e, no fim, da última passagem ao fim da janela.
        idx, t = idx[new_pass], t[new_pass]
        first = np.r_[True, idx[1:] != idx[:-1]]
        before = np.where(first, t, np.r_[0.0, np.diff(t)])
    else:
        # Todos os acessos: só os que abrem passagem trazem um intervalo (desde o acesso anterior
        # ou desde o começo da janela); dentro da passagem, zero.
        before = np.where(first, t, np.where(new_pass, since_prev, 0.0))
    starts = np.flatnonzero(first)
    ends = np.r_[starts[1:] - 1, idx.size - 1]
    gap = np.maximum(np.maximum.reduceat(before, starts), duration_s - t[ends])
    targets = idx[starts]
    max_gap[targets] = gap
    return float(targets.size / access.n_targets), max_gap


def latitude_profile(
    grid: TargetGrid, values: ArrayLike, reduce: str = "mean"
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Valores por alvo reduzidos por faixa de latitude ("mean" ou "max"): (lat [°], valor)."""
    v = np.reshape(np.asarray(values, dtype=float), (grid.band_lat_deg.size, -1))
    if reduce == "mean":
        return grid.band_lat_deg, v.mean(axis=1)
    if reduce == "max":
        return grid.band_lat_deg, v.max(axis=1)
    raise ValueError("reduce deve ser 'mean' ou 'max'")
