"""Perturbação de terceiro corpo (Sol e Lua) para o Cowell (card S4-02).

No referencial geocêntrico, que não é inercial (a Terra também é acelerada pelo corpo), a
aceleração perturbadora do satélite é a diferença entre a atração do corpo sobre o satélite e a
atração do corpo sobre a Terra (Montenbruck & Gill 2000, eq. 3.37):

    a = µ_b · [ (s − r)/|s − r|³ − s/|s|³ ]      s = posição do corpo (ECI da data), r = satélite

O segundo termo é a aceleração da própria Terra. O que sobra é uma **maré**: para r ≪ s,
a ≈ (µ_b/s³)·[3(ŝ·r)ŝ − r], que estica a órbita ao longo da linha Terra–corpo e a comprime na
direção perpendicular. Com o satélite entre a Terra e o corpo, a aceleração aponta para o corpo.

Ordem de grandeza (LEO, r ≈ 6 900 km): |a| ~ 2µ_b·r/s³, ≈ 5e-7 m/s² para o Sol e ≈ 1,2e-6 m/s²
para a Lua (Apêndice A.3: 5,0023e-7 e 1,2005e-6 m/s² em r = R⊕ + 500 km). É ~1e-4 de J2 e ~1e-7 da
gravidade central.

Por que o efeito em Ω̇ da SSO é pequeno: a precessão nodal média causada por um corpo distante é
da ordem de Ω̇ ~ −(3/4)·(µ_b/s³)/n · cos i · (fator da inclinação do corpo), isto é, (n_b²/n) com
n_b² = µ_b/s³ ≈ 4e-14 s⁻² (Sol) e ≈ 9e-14 s⁻² (Lua) contra n ≈ 1,1e-3 s⁻¹: ~1e-10 rad/s antes
do cos i. Na SSO (i ≈ 97,5°, cos i ≈ −0,13) sobra ~1e-11 rad/s ~ 5e-5 °/dia. O experimento do
Apêndice A.3 mede −2,45e-5 °/dia (LTAN −0,035 min/ano), cerca de 40× abaixo da estimativa
"~1e-3 °/dia" do Sem4 §10.

O efeito indireto é maior: o Sol inclina a órbita. A média do torque de maré sobre uma órbita
circular dá di/dt ~ (3/2)·(n☉²/n)·sin β·(…), ≈ −1e-4 °/dia (≈ −0,036°/ano) para a SSO de 22h30
(β ≈ 23°), a conhecida deriva de inclinação das SSO, corrigida por manobras nas missões reais.
Via J2 (dΩ̇/di = −Ω̇·tan i ≈ +7,5 °/dia por rad na SSO), essa deriva faz ΔΩ crescer com t²
(medido no S4-02: −0,0065° em 30 dias só com o Sol). Por isso um ajuste linear de Ω em 30 dias
dá uma "taxa" lunissolar maior (≈ −3,5e-4 °/dia) e dependente da janela; a Lua soma termos de
≈ 14 dias.

O SGP4 de órbitas próximas (período < 225 min) não tem termos lunissolares; esta força é
validada contra as efemérides do JPL e por magnitudes, não contra o SGP4.
"""

from datetime import datetime
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.astrotime import datetime_to_jd, jd_from_epoch
from satsim.constants import MU_MOON, MU_SUN
from satsim.ephemeris import moon_position_eci, sun_position_eci

Body = Literal["sun", "moon"]

_BODIES = {
    "sun": (MU_SUN, sun_position_eci),
    "moon": (MU_MOON, moon_position_eci),
}


def third_body_acceleration(
    r: ArrayLike, body_position: ArrayLike, mu_body: float
) -> NDArray[np.float64]:
    """a = µ_b·[(s − r)/|s − r|³ − s/|s|³] (Montenbruck & Gill 2000, eq. 3.37), vetorizada.

    Args:
        r: posição do satélite no ECI [m], shape (3,) ou (N, 3).
        body_position: posição do corpo s no mesmo quadro [m], shape (3,) ou (N, 3).
        mu_body: parâmetro gravitacional do corpo [m³/s²].

    Returns:
        Aceleração [m/s²] com o shape de broadcast de ``r`` e ``body_position``.

    A diferença dos dois termos perde ~4 algarismos (|r|/|s| ~ 1e-4 a 2e-2), o que em dupla
    precisão ainda deixa ~1e-12 relativo: dispensa a forma de Battin.
    """
    r = np.asarray(r, dtype=float)
    s = np.asarray(body_position, dtype=float)
    d = s - r
    d_norm = np.linalg.norm(d, axis=-1, keepdims=True)
    s_norm = np.linalg.norm(s, axis=-1, keepdims=True)
    return mu_body * (d / d_norm**3 - s / s_norm**3)


class ThirdBodyAcceleration:
    """Perturbação do Sol ou da Lua (implementa ``AccelerationModel``).

    O tempo do Cowell (t segundos desde a época) vira JD por ``epoch + t``; a posição do corpo
    vem de ``satsim.ephemeris`` (ECI da data, UTC usado como TT).

    Args:
        body: ``"sun"`` ou ``"moon"``.
        epoch: instante t = 0, ``datetime`` tz-aware (o mesmo do propagador).

    Raises:
        ValueError: corpo desconhecido ou ``epoch`` naive.
    """

    def __init__(self, body: Body, epoch: datetime) -> None:
        if body not in _BODIES:
            raise ValueError(f"corpo desconhecido: {body!r} (use 'sun' ou 'moon')")
        self.body = body
        self.epoch = epoch
        self.mu, self._position = _BODIES[body]
        self._epoch_jd = datetime_to_jd(epoch)

    def body_position(self, t: ArrayLike) -> NDArray[np.float64]:
        """Posição do corpo no ECI da data [m] em t [s] desde a época, shape (3,) ou (N, 3)."""
        return self._position(jd_from_epoch(self._epoch_jd, t))

    def __call__(
        self, t: float, r: NDArray[np.float64], v: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        return third_body_acceleration(r, self.body_position(t), self.mu)
