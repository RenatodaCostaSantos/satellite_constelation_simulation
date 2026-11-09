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
"""

import math
from collections.abc import Iterator, Sequence
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.astrotime import datetime_to_jd
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
