"""Pressão de radiação solar (SRP) com a sombra da Terra, para o Cowell (card S4-03).

Modelo de "bola de canhão" (cannonball): o satélite é tratado como uma superfície de área A
sempre voltada para o Sol, com coeficiente de refletividade C_R (1 = absorve tudo, 2 = espelho
perfeito; 1,3 é típico). A força empurra para longe do Sol e cai com o quadrado da distância
(Montenbruck & Gill 2000, eq. 3.75; Apêndice A.4 do prompt da Semana 4):

    a = ν · P☉ · C_R · (A/m) · (AU/|r − s|)² · (r − s)/|r − s|

    ν = 0 na sombra, 1 fora      P☉ = 4,56e-6 N/m² a 1 UA      s = posição do Sol (ECI da data)

Ordem de grandeza: com C_R = 1,3 e A/m = 0,01 m²/kg, |a| = 5,928e-8 m/s² a 1 UA, ~1e-8 da
gravidade central e ~1/10 da maré solar em LEO. O efeito secular em a é nulo fora do eclipse
(força conservativa em média); com a sombra, a parte não iluminada da volta quebra a simetria e
surgem efeitos de período longo em e e ω.

Descontinuidade na sombra. Com a sombra cilíndrica, ν salta de 1 para 0 no instante em que o
satélite cruza a borda do cilindro: a aceleração tem um degrau de ≈ 6e-8 m/s² duas vezes por
volta. Isso é aceitável para o DOP853 com passo adaptativo: o controle de erro detecta o degrau
como um erro local maior e reduz o passo em torno dele, e mesmo um passo inteiro de 60 s
"errado" deslocaria a posição em ≈ ½·Δa·h² ≈ 1e-4 m. Integradores de ordem alta perdem ordem
no degrau, mas não precisão relevante. O tratamento rigoroso (penumbra contínua, ou eventos de
entrada e saída da sombra para reiniciar a integração) fica fora do escopo.

Fora do escopo: penumbra (sombra cônica), albedo e radiação infravermelha da Terra, atitude e
área variável.
"""

from datetime import datetime
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.astrotime import datetime_to_jd, jd_from_epoch
from satsim.constants import AU, P_SUN
from satsim.ephemeris import in_shadow_cylindrical, sun_position_eci

ShadowModel = Literal["cylindrical", "none"]


def srp_acceleration(
    r: ArrayLike,
    sun_position: ArrayLike,
    cr: float,
    area_over_mass: float,
    shadow: ShadowModel = "cylindrical",
    p_sun: float = P_SUN,
) -> NDArray[np.float64]:
    """Aceleração de pressão de radiação solar (Apêndice A.4 / M&G eq. 3.75), vetorizada.

    Args:
        r: posição ECI do satélite [m], shape (3,) ou (N, 3).
        sun_position: posição do Sol no mesmo quadro [m], shape (3,) ou (N, 3).
        cr: coeficiente de refletividade C_R [-].
        area_over_mass: A/m [m²/kg].
        shadow: ``"cylindrical"`` (sombra da Terra, ν ∈ {0, 1}) ou ``"none"`` (sempre iluminado).
        p_sun: pressão de radiação a 1 UA [N/m²].

    Returns:
        Aceleração [m/s²] com o shape de broadcast das entradas.
    """
    r = np.asarray(r, dtype=float)
    s = np.asarray(sun_position, dtype=float)
    d = r - s  # do Sol para o satélite
    d_norm = np.linalg.norm(d, axis=-1, keepdims=True)
    acc = p_sun * cr * area_over_mass * (AU / d_norm) ** 2 * d / d_norm
    if shadow == "cylindrical":
        s_hat = s / np.linalg.norm(s, axis=-1, keepdims=True)
        lit = ~in_shadow_cylindrical(r, s_hat)
        acc = acc * lit[..., None]
    elif shadow != "none":
        raise ValueError(f"modelo de sombra desconhecido: {shadow!r}")
    return acc


class SolarRadiationPressure:
    """Pressão de radiação solar com sombra (implementa ``AccelerationModel``).

    O tempo do Cowell (t segundos desde a época) vira JD por ``epoch + t``; a posição do Sol vem
    de ``satsim.ephemeris.sun_position_eci`` (ECI da data).

    Args:
        epoch: instante t = 0, ``datetime`` tz-aware (o mesmo do propagador).
        cr: coeficiente de refletividade C_R [-] (padrão 1,3).
        area_over_mass: A/m [m²/kg] (padrão 0,01).
        shadow: ``"cylindrical"`` (padrão) ou ``"none"``.

    Raises:
        ValueError: modelo de sombra desconhecido ou ``epoch`` naive.
    """

    def __init__(
        self,
        epoch: datetime,
        cr: float = 1.3,
        area_over_mass: float = 0.01,
        shadow: ShadowModel = "cylindrical",
    ) -> None:
        if shadow not in ("cylindrical", "none"):
            raise ValueError(f"modelo de sombra desconhecido: {shadow!r}")
        self.epoch = epoch
        self.cr = cr
        self.area_over_mass = area_over_mass
        self.shadow = shadow
        self._epoch_jd = datetime_to_jd(epoch)

    def sun_position(self, t: ArrayLike) -> NDArray[np.float64]:
        """Posição do Sol no ECI da data [m] em t [s] desde a época, shape (3,) ou (N, 3)."""
        return sun_position_eci(jd_from_epoch(self._epoch_jd, t))

    def __call__(
        self, t: float, r: NDArray[np.float64], v: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        return srp_acceleration(r, self.sun_position(t), self.cr, self.area_over_mass, self.shadow)
