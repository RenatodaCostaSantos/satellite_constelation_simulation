"""Passagens de um satélite sobre uma estação: AOS, TCA e LOS por busca de raízes na elevação.

Uma *passagem* é o intervalo em que o satélite fica acima da elevação mínima ε (``el_min``):

    AOS (Acquisition Of Signal)   instante em que a elevação sobe e cruza ε
    TCA (Time of Closest Approach) instante de elevação máxima (o satélite mais alto no céu)
    LOS (Loss Of Signal)          instante em que a elevação desce e cruza ε

Aqui "passagem" é visibilidade por **rádio** (satélite acima de ε). A regra óptica, satélite
iluminado com o Sol abaixo do horizonte, entra na Semana 4.

Método (``find_passes``)
------------------------
1. **Varredura vetorizada**: a elevação é calculada com ``observe`` numa grade de passo ``step_s``
   (em blocos de 6 h, para limitar a memória).
2. **Refinamento**: cada troca de sinal de el − ε entre duas amostras vizinhas é refinada por
   ``scipy.optimize.brentq`` (``xtol`` = 1e-3 s): subida → AOS, descida → LOS.
3. **TCA**: ``scipy.optimize.minimize_scalar`` (limitado entre AOS e LOS) sobre −el.
4. **Azimutes** em AOS, TCA e LOS.

**Limite do método**: uma passagem mais curta que ``step_s``, que sobe e desce acima de ε entre
duas amostras, escapa da varredura. Com ε = 10° e a ISS, a passagem mais curta do Apêndice A.2
dura 132 s; o passo padrão de 20 s tem folga. Passagens muito rasantes (máximo pouco acima de ε)
são as mais curtas: para garanti-las, reduza ``step_s``.

**Bordas da janela**: se o satélite já está acima de ε no início, a 1ª passagem começa no início
da janela (``starts_inside=True``, ``aos = start``); se ainda está acima no fim, a última termina
no fim (``ends_inside=True``, ``los = start + duration``).
"""

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta, timezone

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from satsim.observer import GroundStation, observe
from satsim.propagators.base import Propagator

BLOCK_S = 6 * 3600.0
"""Tamanho do bloco da varredura vetorizada [s]."""

XTOL_S = 1e-3
"""Tolerância das raízes (AOS, LOS) e do máximo (TCA) [s]."""


@dataclass(frozen=True)
class Pass:
    """Uma passagem do satélite ``sat_id`` sobre a estação.

    ``aos``, ``tca`` e ``los`` são ``datetime`` UTC tz-aware (resolução de microssegundos);
    elevação e azimutes em graus (azimute a partir do norte, sentido horário); ``duration_s`` =
    LOS − AOS. ``starts_inside``/``ends_inside`` marcam passagens cortadas pelas bordas da janela.
    """

    sat_id: str
    aos: datetime
    tca: datetime
    los: datetime
    max_elevation_deg: float
    az_aos_deg: float
    az_tca_deg: float
    az_los_deg: float
    duration_s: float
    starts_inside: bool = False
    ends_inside: bool = False


def _check_aware(dt: datetime, name: str) -> None:
    if dt.tzinfo is None or dt.utcoffset() is None:
        raise ValueError(f"{name} deve ser um datetime tz-aware (ex.: tzinfo=timezone.utc)")


def _as_seconds(duration: float | timedelta) -> float:
    return duration.total_seconds() if isinstance(duration, timedelta) else float(duration)


def find_passes(
    propagator: Propagator,
    station: GroundStation,
    start: datetime,
    duration: float | timedelta,
    el_min_deg: float = 10.0,
    step_s: float = 20.0,
    sat_id: str = "",
) -> list[Pass]:
    """Passagens com elevação acima de ``el_min_deg`` em [start, start + duration].

    Funciona com qualquer ``Propagator`` (só usa ``observe``) e qualquer estação. Ver o método e
    os limites na docstring do módulo.

    Args:
        propagator: Kepler, MeanJ2, Cowell, SGP4…
        station: estação de solo.
        start: início da janela, ``datetime`` tz-aware.
        duration: duração da janela [s] ou ``timedelta``.
        el_min_deg: elevação mínima ε [°].
        step_s: passo da varredura [s]; passagens mais curtas que ele podem escapar.
        sat_id: identificador copiado em cada ``Pass``.

    Returns:
        Passagens em ordem cronológica, cada uma com ``max_elevation_deg ≥ el_min_deg``.

    Raises:
        ValueError: ``start`` naive, duração negativa ou passo não positivo.
    """
    _check_aware(start, "start")
    duration_s = _as_seconds(duration)
    if duration_s < 0.0 or step_s <= 0.0:
        raise ValueError("duration deve ser ≥ 0 e step_s > 0")
    start_utc = start.astimezone(UTC)
    offset_s = (start_utc - propagator.epoch).total_seconds()
    el_min = math.radians(el_min_deg)

    def look(t_rel: float) -> tuple[float, float]:
        obs = observe(propagator, station, offset_s + t_rel)
        return float(obs.el[0]), float(obs.az[0])

    def g(t_rel: float) -> float:  # el − ε: positiva dentro da passagem
        return look(t_rel)[0] - el_min

    # 1. varredura vetorizada, em blocos
    n_steps = max(1, math.ceil(duration_s / step_s))
    t_grid = np.minimum(np.arange(n_steps + 1) * step_s, duration_s)
    above = np.empty(t_grid.size, dtype=bool)
    per_block = max(1, int(BLOCK_S // step_s))
    for k0 in range(0, t_grid.size, per_block):
        sl = slice(k0, k0 + per_block)
        above[sl] = observe(propagator, station, offset_s + t_grid[sl]).el >= el_min

    # 2. trocas de sinal → raízes por brentq
    rises, sets = [], []
    for k in np.flatnonzero(above[:-1] != above[1:]):
        root = brentq(g, t_grid[k], t_grid[k + 1], xtol=XTOL_S)
        (rises if above[k + 1] else sets).append(root)
    starts_inside, ends_inside = bool(above[0]), bool(above[-1])
    if starts_inside:
        rises.insert(0, 0.0)
    if ends_inside:
        sets.append(duration_s)

    passes = []
    for t_aos, t_los in zip(rises, sets, strict=True):
        # 3. TCA: máximo da elevação entre AOS e LOS
        if t_los > t_aos:
            res = minimize_scalar(lambda t: -look(t)[0], bounds=(t_aos, t_los), method="bounded",
                                  options={"xatol": XTOL_S})  # fmt: skip
            t_tca = float(res.x)
        else:
            t_tca = t_aos
        el_tca, az_tca = look(t_tca)
        if el_tca < el_min:  # só passagens que de fato superam ε
            continue
        # 4. azimutes nas bordas
        _, az_aos = look(t_aos)
        _, az_los = look(t_los)
        passes.append(
            Pass(
                sat_id=sat_id,
                aos=start_utc + timedelta(seconds=t_aos),
                tca=start_utc + timedelta(seconds=t_tca),
                los=start_utc + timedelta(seconds=t_los),
                max_elevation_deg=math.degrees(el_tca),
                az_aos_deg=math.degrees(az_aos),
                az_tca_deg=math.degrees(az_tca),
                az_los_deg=math.degrees(az_los),
                duration_s=t_los - t_aos,
                starts_inside=starts_inside and t_aos == 0.0,
                ends_inside=ends_inside and t_los == duration_s,
            )
        )
    return passes


def find_passes_many(
    propagators: Mapping[str, Propagator],
    station: GroundStation,
    start: datetime,
    duration: float | timedelta,
    el_min_deg: float = 10.0,
    step_s: float = 20.0,
) -> list[Pass]:
    """Passagens de vários satélites, ``{sat_id: propagador}``, ordenadas por AOS.

    Cada satélite é buscado com ``find_passes`` (cada propagador com a sua época); o resultado
    junta todos, em ordem de AOS (empates por ``sat_id``).
    """
    passes: list[Pass] = []
    for sat_id, prop in propagators.items():
        passes += find_passes(prop, station, start, duration, el_min_deg, step_s, sat_id)
    return sorted(passes, key=lambda p: (p.aos, p.sat_id))


def _tz(offset_hours: float) -> timezone:
    name = "BRT" if offset_hours == -3.0 else f"UTC{offset_hours:+g}"
    return timezone(timedelta(hours=offset_hours), name)


def format_passes(passes: Sequence[Pass], tz_offset_hours: float = -3.0) -> str:
    """Tabela de texto das passagens, com horários em UTC e no fuso local de offset fixo.

    O padrão é o horário de Brasília, BRT = UTC−03:00. O Brasil não tem horário de verão desde
    2019, então um offset fixo basta e não é preciso depender da base ``zoneinfo``/``tzdata``.
    Passagens cortadas pela janela são marcadas com ``*`` (já acima no início) ou ``+`` (ainda
    acima no fim).
    """
    tz = _tz(tz_offset_hours)
    label = tz.tzname(None)
    w = max([len("satélite"), *(len(p.sat_id) for p in passes)])
    header = (
        f"{'#':>2}  {'satélite':<{w}}  {'AOS UTC':<19}  {'AOS ' + label:<8}  {'TCA ' + label:<8}  "
        f"{'LOS ' + label:<8}  {'el máx':>6}  {'az AOS/TCA/LOS [°]':<19}  {'dur':>5}"
    )
    lines = [header, "-" * len(header)]
    for k, p in enumerate(passes, start=1):
        flag = ("*" if p.starts_inside else "") + ("+" if p.ends_inside else "")
        az = f"{p.az_aos_deg:5.1f}/{p.az_tca_deg:5.1f}/{p.az_los_deg:5.1f}"
        lines.append(
            f"{k:>2}  {p.sat_id:<{w}}  {p.aos:%Y-%m-%d %H:%M:%S}  "
            f"{p.aos.astimezone(tz):%H:%M:%S}  {p.tca.astimezone(tz):%H:%M:%S}  "
            f"{p.los.astimezone(tz):%H:%M:%S}  {p.max_elevation_deg:5.1f}°  {az:<19}  "
            f"{p.duration_s / 60:4.1f}m{flag}"
        )
    if not passes:
        lines.append("(nenhuma passagem)")
    if any(p.starts_inside or p.ends_inside for p in passes):
        lines.append("* já acima da elevação mínima no início da janela; + ainda acima no fim")
    return "\n".join(lines)
