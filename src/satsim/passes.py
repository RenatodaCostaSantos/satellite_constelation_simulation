"""Passagens de um satélite sobre uma estação: AOS, TCA e LOS por busca de raízes na elevação.

O que é uma passagem
--------------------
A estação é um ponto fixo no solo (fixo no ECEF). Dela, a **elevação** de um satélite é o ângulo
acima do horizonte em que ele aparece: 0° no horizonte, 90° no zênite. Atenção: elevação é um
ângulo, não a altura do satélite. Fixada uma elevação mínima ε (``el_min``), por exemplo 10°
para evitar montanhas, prédios e a atmosfera espessa perto do horizonte, o satélite conta como
"visível" quando a elevação é ≥ ε, em **qualquer azimute**.

Essa condição define uma região em duas dimensões, que pode ser vista de dois jeitos
equivalentes:

- **no céu**: um cone em torno da vertical da estação (o zênite), de abertura 90° − ε. Como vale
  para qualquer azimute, a região é circular: uma calota do céu centrada no zênite;
- **no chão**: um círculo em torno da estação, onde precisa estar o *subponto* do satélite (o
  ponto do solo logo abaixo dele). O raio angular, medido no centro da Terra, é o do footprint
  (``observer.footprint_half_angle``, Sem4 §9.4)::

      λ = arccos(R⊕/(R⊕ + h)·cos ε) − ε    ISS a 420 km: ε = 10° → 12,5° ≈ 1 390 km
                                                         ε = 0°  → 20,2° ≈ 2 250 km

  O círculo cresce com a altura h do satélite: quanto mais alto, de mais longe ele é visto. É
  um círculo, e não outra forma, porque a elevação só depende da distância entre o subponto e a
  estação, não da direção.

Uma *passagem* é o **intervalo de tempo** em que o satélite está dentro dessa região, ou seja,
em que o traço no solo atravessa o círculo:

    AOS (Acquisition Of Signal)   instante em que a elevação sobe e cruza ε (entra no círculo)
    TCA (Time of Closest Approach) instante de elevação máxima (o ponto do traço mais próximo
                                  da estação: o satélite mais alto no céu)
    LOS (Loss Of Signal)          instante em que a elevação desce e cruza ε (sai do círculo)

Um traço que corta o círculo pela borda dá uma passagem curta e baixa (A.2, passagem 9: 2,2 min,
11°); um que passa perto do centro dá uma passagem longa e alta (passagem 5: 6,7 min, 71,5°).

O código **não constrói o círculo**: a cada instante calcula a elevação do satélite vista da
estação (``observer.observe``) e testa el ≥ ε, que é a mesma condição. O teste no céu é exato
com o elipsoide WGS84; o círculo de raio λ é uma aproximação esférica, boa para desenhar (o
footprint do viewer, S3-06). ε é o mesmo em todos os azimutes, como se o horizonte fosse plano;
relevo real (São Bento fica na Serra da Mantiqueira) pediria uma máscara de horizonte ε(az), fora
do escopo.

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
    # --- Preparação --------------------------------------------------------------------------
    # Rejeita datetime sem fuso: "10:00" sem fuso é ambíguo (UTC? Brasília?).
    _check_aware(start, "start")
    # A duração pode vir como timedelta ou como segundos; daqui em diante, sempre segundos.
    duration_s = _as_seconds(duration)
    if duration_s < 0.0 or step_s <= 0.0:
        raise ValueError("duration deve ser ≥ 0 e step_s > 0")
    # Todos os instantes de saída (AOS, TCA, LOS) ficam em UTC, mesmo que start venha em BRT.
    start_utc = start.astimezone(UTC)
    # Internamente o tempo é contado em segundos desde start (t_rel = 0 no início da janela),
    # mas o propagador conta desde a SUA época (ex.: a época do TLE). offset_s converte:
    #   t_propagador = offset_s + t_rel
    offset_s = (start_utc - propagator.epoch).total_seconds()
    # O usuário passa ε em graus; observe devolve a elevação em radianos.
    el_min = math.radians(el_min_deg)

    def look(t_rel: float) -> tuple[float, float]:
        """(elevação, azimute) [rad] num único instante t_rel [s desde start]."""
        obs = observe(propagator, station, offset_s + t_rel)
        return float(obs.el[0]), float(obs.az[0])

    def g(t_rel: float) -> float:  # el − ε: positiva dentro da passagem
        """Função cuja raiz é o cruzamento de ε: negativa abaixo, positiva acima."""
        return look(t_rel)[0] - el_min

    # --- 1. Varredura vetorizada: acima ou abaixo de ε a cada step_s ---------------------------
    # Número de passos para cobrir a janela; ceil garante que a grade chega até o fim.
    n_steps = max(1, math.ceil(duration_s / step_s))
    # Grade 0, step_s, 2·step_s, …; o np.minimum faz o último ponto cair exatamente no fim da
    # janela (duration_s), mesmo que a duração não seja múltipla do passo.
    t_grid = np.minimum(np.arange(n_steps + 1) * step_s, duration_s)
    # above[k] = True se o satélite está acima de ε no instante t_grid[k].
    above = np.empty(t_grid.size, dtype=bool)
    # Quantas amostras cabem num bloco de 6 h: calcular dias inteiros de uma vez criaria arrays
    # grandes demais para janelas longas; em blocos, a memória fica limitada.
    per_block = max(1, int(BLOCK_S // step_s))
    for k0 in range(0, t_grid.size, per_block):
        sl = slice(k0, k0 + per_block)  # as amostras k0, k0+1, …, k0+per_block−1
        # Uma única chamada vetorizada de observe para todo o bloco (sem laço por instante).
        above[sl] = observe(propagator, station, offset_s + t_grid[sl]).el >= el_min

    # --- 2. Onde above troca de valor, a elevação cruzou ε: refina o instante exato ------------
    rises, sets = [], []  # rises: instantes de subida (AOS); sets: de descida (LOS) [s]
    # above[:-1] != above[1:] compara cada amostra com a seguinte; flatnonzero devolve os k em
    # que houve troca, ou seja, o cruzamento ocorreu entre t_grid[k] e t_grid[k+1].
    for k in np.flatnonzero(above[:-1] != above[1:]):
        # g troca de sinal nesse intervalo de step_s segundos; brentq estreita o intervalo até
        # achar a raiz g = 0 (elevação exatamente ε) com precisão de XTOL_S (1 ms).
        root = brentq(g, t_grid[k], t_grid[k + 1], xtol=XTOL_S)
        # Se depois da troca está acima (F → T), o satélite subiu: AOS. Senão (T → F): LOS.
        (rises if above[k + 1] else sets).append(root)
    # Bordas: se já está acima no 1º instante, a subida aconteceu antes da janela e não aparece
    # como troca; idem para a descida no fim. Sem correção, sobraria um LOS sem AOS (ou o
    # contrário) e o emparelhamento abaixo falharia.
    starts_inside, ends_inside = bool(above[0]), bool(above[-1])
    if starts_inside:
        rises.insert(0, 0.0)  # subida artificial no início da janela
    if ends_inside:
        sets.append(duration_s)  # descida artificial no fim da janela

    # --- 3 e 4. Emparelha subidas e descidas em passagens; acha TCA e azimutes -----------------
    passes = []
    # Subidas e descidas se alternam no tempo (sobe, desce, sobe, desce…), então a k-ésima
    # subida casa com a k-ésima descida. strict=True acusa erro se as listas não emparelharem.
    for t_aos, t_los in zip(rises, sets, strict=True):
        # 3. TCA: o instante de elevação máxima entre AOS e LOS. minimize_scalar minimiza, então
        #    minimizamos −el; method="bounded" restringe a busca ao intervalo [AOS, LOS].
        if t_los > t_aos:
            res = minimize_scalar(lambda t: -look(t)[0], bounds=(t_aos, t_los), method="bounded",
                                  options={"xatol": XTOL_S})  # fmt: skip
            t_tca = float(res.x)
        else:  # intervalo de duração nula (janela de duração 0): não há o que maximizar
            t_tca = t_aos
        el_tca, az_tca = look(t_tca)
        # Salvaguarda: por construção el ≥ ε entre AOS e LOS, mas arredondamentos nas raízes de
        # uma passagem rasante podem deixar o máximo abaixo de ε; essas não contam.
        if el_tca < el_min:  # só passagens que de fato superam ε
            continue
        # 4. Azimutes nas bordas: de onde o satélite surge (AOS) e onde some (LOS).
        _, az_aos = look(t_aos)
        _, az_los = look(t_los)
        passes.append(
            Pass(
                sat_id=sat_id,
                # Segundos desde start → datetime UTC (timedelta arredonda ao microssegundo).
                aos=start_utc + timedelta(seconds=t_aos),
                tca=start_utc + timedelta(seconds=t_tca),
                los=start_utc + timedelta(seconds=t_los),
                # Saída em graus, como no nome dos campos.
                max_elevation_deg=math.degrees(el_tca),
                az_aos_deg=math.degrees(az_aos),
                az_tca_deg=math.degrees(az_tca),
                az_los_deg=math.degrees(az_los),
                duration_s=t_los - t_aos,
                # Só a passagem que usou a subida/descida artificial é marcada como cortada.
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
