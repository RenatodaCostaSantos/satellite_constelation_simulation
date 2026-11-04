"""Subponto geodésico e traço no solo (ground track), base matemática do mapa 2D da GUI.

O subponto é o ponto da superfície diretamente abaixo do satélite: a posição ECI do propagador é
girada para ECEF pelo GMST (``frames``) e convertida para latitude, longitude e altitude
geodésicas WGS84 (``geodesy``). Como a Terra gira por baixo da órbita, o traço se desloca para
oeste a cada revolução: Δλ ≈ −(ω⊕ − Ω̇)·Tnod por volta (Sem4 §7.1).

Funciona com qualquer ``Propagator`` (Kepler, MeanJ2, Cowell), sem caso especial.

Nós ascendentes e a grade de trilhas (card S5-02)
-------------------------------------------------
``node_crossings`` acha os nós ascendentes de qualquer propagador numa janela que começa em
qualquer instante UTC (``propagate_utc``). ``grid_metrics`` mede a grade que os nós depositam no
equador (Sem4 §7.4 e §8.3): espaçamentos entre vizinhos e o afastamento à grade ideal de
360°/N. ``repeat_closure`` dá quanto o traço deixa de fechar depois de N revoluções.
"""

from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.frames import eci_to_ecef, gmst_at
from satsim.geodesy import ecef_to_geodetic
from satsim.propagators.base import Propagator, as_time_array


def _wrap_pi(angle: NDArray[np.float64]) -> NDArray[np.float64]:
    """Reduz ângulos [rad] a [−π, π)."""
    return np.mod(angle + np.pi, 2.0 * np.pi) - np.pi


def subpoint(
    propagator: Propagator, t: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Subponto geodésico WGS84 do satélite.

    Args:
        propagator: qualquer ``Propagator``.
        t: segundos desde ``propagator.epoch``, escalar ou (N,).

    Returns:
        (lat, lon, alt): latitude geodésica [rad], longitude [rad] em (−π, π] e altitude
        elipsoidal [m], cada um com shape (N,) (N = 1 para t escalar).
    """
    t = as_time_array(t)
    r_eci, _ = propagator.propagate(t)
    r_ecef = eci_to_ecef(r_eci, gmst_at(propagator.epoch, t))
    lat, lon, alt = ecef_to_geodetic(r_ecef)
    return np.atleast_1d(lat), np.atleast_1d(lon), np.atleast_1d(alt)


def split_at_antimeridian(
    lon: ArrayLike, lat: ArrayLike
) -> list[tuple[NDArray[np.float64], NDArray[np.float64]]]:
    """Quebra o traço onde ele cruza o antimeridiano (±180°), para desenhar sem riscos no mapa.

    Um salto |Δλ| > π entre amostras consecutivas indica que o traço saiu por uma borda do mapa e
    entrou pela outra; ali começa um novo segmento. Nenhum ponto é criado nem descartado.

    Args:
        lon: longitudes [rad] em (−π, π], shape (N,).
        lat: latitudes [rad], shape (N,).

    Returns:
        Lista de segmentos (lon, lat), na ordem do tempo; dentro de cada um, |Δλ| ≤ π.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    breaks = np.flatnonzero(np.abs(np.diff(lon)) > np.pi) + 1
    return list(zip(np.split(lon, breaks), np.split(lat, breaks), strict=True))


def ascending_node_crossings(
    t: ArrayLike, lat: ArrayLike, lon: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Instantes e longitudes dos nós ascendentes (latitude cruzando 0 de sul para norte).

    Para cada par de amostras com lat[k] < 0 ≤ lat[k+1], interpola linearmente o instante em que
    lat = 0 e a longitude nesse instante (com a diferença de longitude reduzida a [−π, π), para
    não interpolar através do antimeridiano).

    Args:
        t: tempos [s], shape (N,), crescentes.
        lat: latitudes [rad], shape (N,).
        lon: longitudes [rad], shape (N,).

    Returns:
        (t_nodes, lon_nodes): instantes [s] e longitudes [rad] em [−π, π), shape (K,).
    """
    t = np.asarray(t, dtype=float)
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    k = np.flatnonzero((lat[:-1] < 0.0) & (lat[1:] >= 0.0))
    frac = -lat[k] / (lat[k + 1] - lat[k])
    t_nodes = t[k] + frac * (t[k + 1] - t[k])
    lon_nodes = _wrap_pi(lon[k] + frac * _wrap_pi(lon[k + 1] - lon[k]))
    return t_nodes, lon_nodes


def node_crossings(
    propagator: Propagator, start: datetime, duration_s: float, step_s: float = 30.0
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Nós ascendentes de ``propagator`` na janela [start, start + duration_s] (S5-02).

    Amostra o subponto geodésico a cada ``step_s`` (``propagate_utc`` a partir de ``start``, GMST
    da fórmula IAU-82 completa) e interpola linearmente o cruzamento do equador de sul para norte
    (``ascending_node_crossings``). Latitude geodésica nula ⇔ z_ECEF = 0 ⇔ z_ECI = 0, então é o
    nó da convenção do projeto (z_ECI = 0 com vz > 0). Com passo de 30 s o erro contra o nó
    refinado por ``brentq`` é ≈ 1e-5° em longitude e ≈ 3e-4 s no instante (Apêndice A.3 do
    prompt da Semana 5).

    Args:
        propagator: qualquer ``Propagator`` (ou ``Satellite``), com qualquer época. O Cowell exige
            que a janela comece depois da sua época.
        start: início da janela, ``datetime`` tz-aware.
        duration_s: duração da janela [s] (≥ 0); o último instante é ``start + duration_s``.
        step_s: passo da amostragem [s].

    Returns:
        (t_nodes, lon_nodes): instantes [s desde ``start``] e longitudes [rad] em [−π, π), em
        ordem de tempo, shape (K,).

    Raises:
        ValueError: ``step_s`` não positivo ou ``duration_s`` negativa.
    """
    if step_s <= 0.0 or duration_s < 0.0:
        raise ValueError("step_s deve ser positivo e duration_s não negativo")
    t = np.arange(0.0, duration_s, step_s)
    t = np.append(t, float(duration_s)) if t.size == 0 or t[-1] < duration_s else t
    r_eci, _ = propagator.propagate_utc(start, t)
    lat, lon, _ = ecef_to_geodetic(eci_to_ecef(r_eci, gmst_at(start, t)))
    return ascending_node_crossings(t, np.atleast_1d(lat), np.atleast_1d(lon))


def grid_metrics(lon_nodes: ArrayLike, n_tracks: int) -> tuple[int, float, float, float]:
    """Métricas da grade de nós no equador (Sem4 §7.4, eq. 39; §8.4, Figura 12) (S5-02).

    Os nós são ordenados em longitude (mod 2π) e os espaçamentos entre vizinhos incluem a volta de
    360° (do último ao primeiro). O afastamento de cada nó à grade ideal de passo 2π/N, ancorada
    no primeiro nó de ``lon_nodes``, é a distância à linha de grade mais próxima. Para uma
    constelação, passe os nós satélite a satélite com S1 primeiro (âncora = 1º nó de S1, como no
    Apêndice A.3 do prompt da Semana 5); a âncora muda o afastamento máximo, não os espaçamentos.

    Args:
        lon_nodes: longitudes dos nós [rad], shape (K,), K ≥ 1.
        n_tracks: N, número de trilhas da grade ideal (91 na órbita 91/6).

    Returns:
        (count, gap_min, gap_max, max_offset): número de nós, menor e maior espaçamento entre
        vizinhos [rad] e maior afastamento à grade ideal [rad] (≤ π/N).

    Raises:
        ValueError: sem nós ou ``n_tracks`` < 1.
    """
    lon = np.asarray(lon_nodes, dtype=float).ravel()
    if lon.size == 0:
        raise ValueError("grid_metrics precisa de ao menos um nó")
    if n_tracks < 1:
        raise ValueError("n_tracks deve ser ≥ 1")
    ordered = np.sort(np.mod(lon, 2.0 * np.pi))
    gaps = np.diff(np.append(ordered, ordered[0] + 2.0 * np.pi))
    step = 2.0 * np.pi / n_tracks
    rel = _wrap_pi(lon - lon[0])
    offset = np.abs(rel - np.round(rel / step) * step)
    return int(lon.size), float(gaps.min()), float(gaps.max()), float(offset.max())


def repeat_closure(lon_nodes: ArrayLike, n_revs: int) -> float:
    """Fechamento do traço após N revoluções: λ_N − λ_0 em [−π, π) [rad] (Sem4 §7.3, §7.6).

    Numa órbita de repetição N/D exata, o nó N cai sobre o nó 0 (D dias depois). Positivo = o nó
    N cai a leste do nó 0 (o ciclo de N voltas terminou antes de D dias nodais). Multiplique
    por R⊕ para km no equador.

    Args:
        lon_nodes: longitudes [rad] dos nós de um satélite, em ordem de tempo, K ≥ N + 1.
        n_revs: N.

    Raises:
        ValueError: se houver menos de N + 1 nós.
    """
    lon = np.asarray(lon_nodes, dtype=float)
    if lon.size < n_revs + 1:
        raise ValueError(f"repeat_closure precisa de {n_revs + 1} nós, recebeu {lon.size}")
    return float(_wrap_pi(np.asarray(lon[n_revs] - lon[0])))
