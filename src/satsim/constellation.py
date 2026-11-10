"""Constelação: P satélites propagados juntos em arrays (P, N, 3) (ADR 0001; card S5-01).

``Constellation`` empilha os estados de vários ``Satellite``, cada um com a sua época, numa
grade de tempo comum: ``t`` em segundos desde ``Constellation.epoch``, convertido para cada
satélite como no ADR (``t_i = t + (epoch − sat.epoch)``, via ``propagate_utc``).

``sso_constellation`` é a fábrica da constelação Epic of Sun (Sem4 §8, Tabela 4): P satélites
num único plano heliossíncrono, que compartilham a, e, i, Ω e ω médios e diferem só na fase
ψ0 (Sem4 eq. 42). Segue a receita da Semana 4 (docs/semana4_fechamento.md, seção 6):

1. projeto médio: a de projeto (Sem4 §7.5), e congelada (``frozen_eccentricity``), ω = 90°;
2. calibração: i heliossíncrona calibrada no Cowell só com a gravidade zonal
   (``sso_inclination_calibrated``), **uma vez para todos** (mesmos a e e);
3. Ω para o LTAN pedido na época (``raan_for_ltan``);
4. conversão médio → osculador de cada satélite só com a gravidade zonal
   (``mean_osc.conversion_models``) e propagação no Cowell com todas as forças.

É a única forma de criar a constelação nos scripts e na demo.

Desenho da constelação (card S6-01)
-----------------------------------
``ConstellationDesign`` guarda (N, D, P) e as grandezas do Sem4 §8.2–8.3 (m = P/D, uniformidade,
trilhas por dia, grade diária, atraso entre nós de vizinhos, famílias, rolagem para cobertura
diária). ``CONSTELLATION_DESIGNS`` registra os presets da órbita 91/6 com 6 e 18 satélites: é a
única fonte de P para scripts e demo; ``ConstellationDesign.build`` cria a constelação com
``sso_constellation`` e a paleta por família (``constellation_palette``).
"""

import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.astrotime import datetime_to_jd
from satsim.coverage import roll_for_half_grid
from satsim.elements import KeplerianElements
from satsim.forces import AccelerationModel, DragAcceleration
from satsim.mean_osc import conversion_models, mean_to_osculating
from satsim.orbit_design import (
    constellation_phases,
    frozen_eccentricity,
    raan_for_ltan,
    sso_inclination,
    sso_inclination_calibrated,
)
from satsim.propagators import CowellPropagator
from satsim.propagators.base import as_time_array
from satsim.satellite import Satellite

DESIGN_SEMI_MAJOR_AXIS_M = 6888.089e3
"""a médio de projeto da órbita 91/6 [m] (Sem4 §7.5, teoria só J2; calibrado no S5-02)."""

DESIGN_LTAN_H = 22.5
"""LTAN de projeto [h]: nó ascendente às 22h30, descendente às 10h30 (Sem4 §7.5)."""

CONSTELLATION_PALETTE: tuple[str, ...] = (
    "#00e5ff",  # ciano (a cor da SSO-1 das semanas anteriores)
    "#ff4d8d",  # rosa
    "#7cff6b",  # verde-limão
    "#ffd23f",  # amarelo
    "#8f7bff",  # violeta
    "#ff8a3d",  # laranja
)
"""Seis cores vivas e distintas entre si e das da ISS (#ffb020) e do Sentinel-2A (#c084fc)."""


class Constellation:
    """Conjunto de satélites propagados numa grade de tempo comum (ADR, Interfaces).

    Args:
        epoch: época comum (tz-aware): o t = 0 de ``propagate``.
        satellites: os satélites, cada um com o seu propagador e a sua época; ids únicos.

    Raises:
        ValueError: época naive, nenhum satélite ou ids repetidos.
    """

    def __init__(self, epoch: datetime, satellites: Sequence[Satellite]) -> None:
        if epoch.tzinfo is None or epoch.utcoffset() is None:
            raise ValueError("epoch deve ser um datetime tz-aware (ex.: tzinfo=timezone.utc)")
        satellites = tuple(satellites)
        if not satellites:
            raise ValueError("a constelação precisa de ao menos um satélite")
        ids = [sat.id for sat in satellites]
        if len(set(ids)) != len(ids):
            raise ValueError(f"ids de satélite repetidos: {ids}")
        self._epoch = epoch
        self._satellites = satellites

    @property
    def epoch(self) -> datetime:
        """Época comum (t = 0 de ``propagate``)."""
        return self._epoch

    @property
    def satellites(self) -> tuple[Satellite, ...]:
        """Os satélites, na ordem do eixo P."""
        return self._satellites

    @property
    def ids(self) -> list[str]:
        """Ids dos satélites, na ordem do eixo P."""
        return [sat.id for sat in self._satellites]

    def __len__(self) -> int:
        return len(self._satellites)

    def __iter__(self) -> Iterator[Satellite]:
        return iter(self._satellites)

    def __getitem__(self, index: int) -> Satellite:
        return self._satellites[index]

    def propagate(self, t: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estados ECI de todos os satélites em ``t`` segundos desde ``self.epoch``.

        Cada satélite é amostrado no mesmo instante absoluto, ``epoch + t``, convertido para a
        sua própria época (ADR: t_i = t + (epoch − sat.epoch)).

        Returns:
            (r, v) em ECI [m, m/s], cada um com shape (P, N, 3).
        """
        return self.propagate_utc(self._epoch, t)

    def propagate_utc(
        self, epoch: datetime, t: ArrayLike
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estados ECI de todos os satélites nos instantes ``epoch + t``; shape (P, N, 3) cada."""
        t = as_time_array(t)
        states = [sat.propagate_utc(epoch, t) for sat in self._satellites]
        return np.stack([r for r, _ in states]), np.stack([v for _, v in states])


def sso_constellation(
    epoch: datetime,
    models: Sequence[AccelerationModel],
    n_sats: int = 6,
    a: float = DESIGN_SEMI_MAJOR_AXIS_M,
    ltan_h: float = DESIGN_LTAN_H,
    i: float | None = None,
    e: float | None = None,
    color_palette: Sequence[str] = CONSTELLATION_PALETTE,
) -> Constellation:
    """Constelação heliossíncrona de ``n_sats`` satélites num único plano (Sem4 §8, Tabela 4).

    Elementos médios comuns: a; e = ``frozen_eccentricity(a, sso_inclination(a))`` se não for
    dada, com a inclinação de 1ª ordem só J2 (eq. 27), como na receita da Semana 4 e no A.1
    (1,07721e-3), seja a ``i`` dada ou calibrada: assim e não depende de como i foi obtida (com a
    i calibrada no lugar, e mudaria 5,6e-5 relativo, 1,077149e-3);
    i = ``sso_inclination_calibrated(a, e, epoch, conversion_models(models))`` se não for dada
    (uma calibração para todos); Ω = ``raan_for_ltan(ltan_h, jd(epoch))``; ω = 90°.
    Fases ψ0_j = ``constellation_phases(n_sats)`` (eq. 42): S1 em 0 e cada seguinte 360°/P atrás.

    Cada satélite: médio → osculador com ``conversion_models(models)`` (só gravidade zonal) e
    ``CowellPropagator`` com ``models`` (todas as forças), época ``epoch``. Ids ``s1`` … ``sP``,
    nomes ``S1`` … ``SP``; cores de ``color_palette`` (repetidas em ciclo se P for maior). C_D e
    A/m são copiados do ``DragAcceleration`` de ``models``, se houver.

    Args:
        epoch: época dos elementos médios e da constelação (tz-aware).
        models: forças do Cowell (ex.: J2–J4 + Sol/Lua + SRP + arrasto).
        n_sats: P, número de satélites.
        a: semieixo maior médio [m] (padrão: o de projeto, 6 888,089 km).
        ltan_h: LTAN [h] na época.
        i: inclinação média [rad]; ``None`` calibra no Cowell (≈ 5 s).
        e: excentricidade média; ``None`` usa a congelada J2–J4.
        color_palette: cores CSS dos satélites.

    Returns:
        ``Constellation`` com época ``epoch``.
    """
    if n_sats < 1:
        raise ValueError("n_sats deve ser ≥ 1")
    if not color_palette:
        raise ValueError("color_palette não pode ser vazia")
    conv_models = conversion_models(models)
    if e is None:
        e = frozen_eccentricity(a, float(sso_inclination(a)))
    if i is None:
        i = sso_inclination_calibrated(a, e, epoch, conv_models)
    raan = float(raan_for_ltan(ltan_h, datetime_to_jd(epoch)))
    # C_D e A/m informativos, do modelo de arrasto (se houver): o arrasto já está em `models`.
    drag = next((m for m in models if isinstance(m, DragAcceleration)), None)
    cd = drag.cd if drag is not None else None
    area_over_mass = drag.area_over_mass if drag is not None else None

    satellites = []
    for j, psi0 in enumerate(constellation_phases(n_sats)):
        mean = KeplerianElements(a=a, e=e, i=i, raan=raan, argp=math.pi / 2.0, psi=float(psi0))
        osc = mean_to_osculating(mean, epoch, conv_models).elements
        satellites.append(
            Satellite(
                id=f"s{j + 1}",
                name=f"S{j + 1}",
                propagator=CowellPropagator(osc, epoch, models),
                color=color_palette[j % len(color_palette)],
                cd=cd,
                area_over_mass=area_over_mass,
            )
        )
    return Constellation(epoch, satellites)


# --- Desenho da constelação (S6-01) ---------------------------------------------------------------


@dataclass(frozen=True)
class ConstellationDesign:
    """Constelação de P satélites igualmente faseados num plano, na órbita de repetição N/D.

    Grandezas do Sem4 §8.2–8.3 (eqs. 42 e 43): o satélite j (j = 0 … P − 1, id ``s{j+1}``)
    tem fase ψ0_j = −j·360°/P, cruza o nó Tnod/P depois do satélite j − 1 e deposita o seu nó
    D/P passos de grade a oeste. O padrão de nós é uniforme quando P | D ou D | P.

    Attributes:
        n_revs: N, revoluções do ciclo de repetição.
        n_days: D, dias do ciclo.
        n_sats: P, número de satélites.
        name: nome legível (ex.: "91/6 · 18 satélites").
    """

    n_revs: int
    n_days: int
    n_sats: int
    name: str

    def __post_init__(self) -> None:
        if min(self.n_revs, self.n_days, self.n_sats) < 1:
            raise ValueError("n_revs, n_days e n_sats devem ser ≥ 1")

    @property
    def m(self) -> float:
        """m = P/D, satélites por dia do ciclo (Sem4 §8.3)."""
        return self.n_sats / self.n_days

    @property
    def is_uniform(self) -> bool:
        """Padrão de nós uniforme: P | D ou D | P (Sem4 §8.3)."""
        return self.n_days % self.n_sats == 0 or self.n_sats % self.n_days == 0

    @property
    def n_families(self) -> int:
        """Número de famílias: m quando D | P (cada família é uma constelação N/D com P = D);
        1 nos outros casos."""
        return self.n_sats // self.n_days if self.n_sats % self.n_days == 0 else 1

    @property
    def tracks_per_day(self) -> int | None:
        """Trilhas distintas por dia, N·P/D, quando D | P (Sem4 §8.3); ``None`` se não."""
        if self.n_sats % self.n_days:
            return None
        return self.n_revs * self.n_sats // self.n_days

    @property
    def grid_deg(self) -> float:
        """Grade do ciclo, 360°/N."""
        return 360.0 / self.n_revs

    @property
    def daily_grid_deg(self) -> float | None:
        """Grade diária, 360°·D/(N·P), quando D | P (Sem4 eq. 43); ``None`` se não."""
        tracks = self.tracks_per_day
        return None if tracks is None else 360.0 / tracks

    @property
    def full_grid_days(self) -> int | None:
        """Dias para os P satélites depositarem a grade completa de N trilhas: 1 se D | P,
        D/P se P | D (revisita da grade); ``None`` se o padrão não for uniforme."""
        if self.n_sats % self.n_days == 0:
            return 1
        if self.n_days % self.n_sats == 0:
            return self.n_days // self.n_sats
        return None

    @property
    def phase_step_deg(self) -> float:
        """Δu entre satélites vizinhos, 360°/P (Sem4 eq. 42)."""
        return 360.0 / self.n_sats

    @property
    def ids(self) -> tuple[str, ...]:
        """Ids ``s1`` … ``sP``, em ordem de fase."""
        return tuple(f"s{j + 1}" for j in range(self.n_sats))

    def node_delay_s(self, t_nod: float) -> float:
        """Atraso entre os nós de satélites vizinhos, Tnod/P [s] (Sem4 eq. 42)."""
        return t_nod / self.n_sats

    def family(self, j: int) -> int:
        """Família do satélite j (0-based): j mod m quando D | P; 0 nos outros casos.

        A família f reúne s_{f+1}, s_{f+1+m}, …: uma constelação N/D completa com P = D,
        deslocada f/m de passo de grade (Apêndice A.1 do prompt da Semana 6).
        """
        if not 0 <= j < self.n_sats:
            raise ValueError(f"j fora de 0 … {self.n_sats - 1}: {j}")
        return j % self.n_families

    def roll_for_daily_coverage(self, a: float) -> float:
        """Rolagem η [rad] que fecha a grade diária: meia grade de N·P/D trilhas (Sem4 eq. 45).

        91/6: 23,1957° com P = 6 e 8,1835° com P = 18 (a = 6 888,089 km; Apêndice A.2 do prompt
        da Semana 6).

        Raises:
            ValueError: se D não dividir P (sem grade diária).
        """
        tracks = self.tracks_per_day
        if tracks is None:
            raise ValueError("sem grade diária: D não divide P")
        return roll_for_half_grid(a, tracks)

    def build(
        self,
        epoch: datetime,
        models: Sequence[AccelerationModel],
        a: float = DESIGN_SEMI_MAJOR_AXIS_M,
        ltan_h: float = DESIGN_LTAN_H,
        i: float | None = None,
        e: float | None = None,
    ) -> Constellation:
        """A constelação deste desenho por ``sso_constellation`` (receita das Semanas 4–5), com
        ``n_sats = P`` e a paleta ``constellation_palette(self)``.

        O ``a`` padrão é o de projeto da órbita 91/6; para outro N/D, passe o a da órbita.
        """
        return sso_constellation(
            epoch, models, n_sats=self.n_sats, a=a, ltan_h=ltan_h, i=i, e=e,
            color_palette=constellation_palette(self),
        )  # fmt: skip


DESIGN_91_6_P6 = ConstellationDesign(91, 6, 6, "91/6 · 6 satélites")
"""Constelação Epic of Sun da Semana 5 (Sem4 Tabela 4): P = D, Δu = 60°, 91 trilhas por dia."""

DESIGN_91_6_P18 = ConstellationDesign(91, 6, 18, "91/6 · 18 satélites")
"""P = 3·D (Sem4 §8.3): Δu = 20°, 273 trilhas por dia; s1, s4, …, s16 são a constelação de 6."""

CONSTELLATION_DESIGNS: dict[int, ConstellationDesign] = {6: DESIGN_91_6_P6, 18: DESIGN_91_6_P18}
"""Registro dos desenhos por P: a única fonte para scripts e demo."""

_TONE_MIX = 0.42
"""Fração de mistura com branco/preto que diferencia as famílias 1, 2, … na paleta."""


def _mix(color: str, target: float, frac: float) -> str:
    """Mistura a cor CSS ``#rrggbb`` com cinza ``target`` (0 = preto, 255 = branco)."""
    rgb = [int(color[k : k + 2], 16) for k in (1, 3, 5)]
    out = [round(c + (target - c) * frac) for c in rgb]
    return "#" + "".join(f"{c:02x}" for c in out)


def constellation_palette(design: ConstellationDesign) -> tuple[str, ...]:
    """Cores dos P satélites do desenho, agrupadas por família (S6-01).

    Quando D | P, o matiz vem da posição k = j // m no ciclo (``CONSTELLATION_PALETTE[k]``,
    ciclando se D > 6) e a tonalidade da família f = j mod m: f = 0 usa a cor base (s1, s4, …,
    s16 têm as cores de S1 … S6 da Semana 5); f ímpar escurece e f par clareia, cada vez mais.
    Com P = 6 devolve exatamente ``CONSTELLATION_PALETTE``. Fora de D | P, as cores base ciclam.
    """
    m = design.n_families
    base = CONSTELLATION_PALETTE
    colors = []
    for j in range(design.n_sats):
        k, f = (j // m, j % m) if m > 1 else (j, 0)
        color = base[k % len(base)]
        if f:
            frac = min(_TONE_MIX * ((f + 1) // 2), 0.8)
            color = _mix(color, 0.0 if f % 2 else 255.0, frac)
        colors.append(color)
    return tuple(colors)
