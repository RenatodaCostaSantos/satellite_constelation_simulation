"""Forças do Cowell com as mesmas premissas físicas do SGP4 (ablação "Cowell equivalente").

O SGP4 é uma teoria analítica: resolve *aproximadamente* um problema físico bem definido. Este
módulo monta, para o Cowell, *esse mesmo problema*, para que a comparação Cowell × SGP4 meça só
a diferença entre resolver as equações exatamente (integração numérica) e por série truncada:

- **Gravidade:** zonais J2, J3 e J4 com as constantes **WGS-72** com que o SGP4 trabalha e com que
  os TLEs são ajustados (µ = 398 600,8 km³/s², R⊕ = 6378,135 km; ``sgp4.earth_gravity.wgs72``).
  O Cowell precisa receber ``mu=WGS72_MU`` também no termo central.
- **Arrasto:** a atmosfera de lei de potência do SGP4 (Hoots & Roehrich 1980, Spacetrack Report
  #3; Vallado et al. 2006),

      ρ(r) = ρ0 · [(q0 − s)/(r − s)]⁴,    q0 = R⊕ + 120 km,   s = R⊕ + 78 km,

  esférica, **sem co-rotação** (v_rel = v), com a intensidade dada pelo B* do TLE.
- **Sem** Sol, Lua nem SRP: o SGP4 de órbitas próximas (período < 225 min) não os tem.

B* e o "C_D·A/m" equivalente
-----------------------------
O SGP4 define B* = ½·(C_D·A/m)·ρ0, com ρ0 = 2,461e-5 kg/m²/ER, e B* em 1/ER (ER = raio da Terra).
A aceleração resultante, em SI, é

    a_D = −½·(C_D·A/m)·ρ(r)·|v|·v = −(B*/R⊕)·[(q0 − s)/(r − s)]⁴·|v|·v.

Para usar o ``DragAcceleration`` do projeto, que calcula ½·(C_D·A/m)·ρ, a densidade é
ρ(r) = (ρ0/R⊕)·[...]⁴ [kg/m³] e o produto é C_D·A/m = 2·B*/ρ0 [m²/kg]. Esses dois números **não
são físicos** separadamente: a lei de potência dá ≈ 1e-15 kg/m³ a 400 km (a densidade real é
≈ 4e-12) e o C_D·A/m equivalente sai em m²/kg em vez de centésimos. Só o produto, que é o que a
trajetória enxerga, tem sentido: ele reproduz a queda secular do semieixo que o SGP4 aplica
(da/dt = −2·a·cc1) a 0,3% na ISS (``tests/test_sgp4_equivalent.py``).

Diferenças que ficam (e que a ablação mede)
-------------------------------------------
O SGP4 trunca a teoria: J2 entra nos termos de curto período só em 1ª ordem, J3 só nos de
longo período, J4 só no secular, e o arrasto entra só como deriva secular de a, e e da anomalia
média. O Cowell aplica as mesmas forças por inteiro, em cada ponto. Também fica de fora o ajuste
do SGP4 para perigeus abaixo de 220 km (``isimp``), que só simplifica os termos de arrasto.
"""

from typing import Final

import numpy as np
from numpy.typing import ArrayLike
from sgp4.earth_gravity import wgs72

from satsim.constants import R_EARTH
from satsim.forces.atmosphere import FloatOrArray
from satsim.forces.drag import DragAcceleration
from satsim.forces.zonal import ZonalAcceleration

WGS72_MU: Final[float] = wgs72.mu * 1e9
"""µ do WGS-72 usado pelo SGP4 [m³/s²] (398 600,8 km³/s²)."""

WGS72_RE: Final[float] = wgs72.radiusearthkm * 1e3
"""R⊕ do WGS-72 usado pelo SGP4 [m] (6378,135 km)."""

WGS72_ZONAL: Final[dict[int, float]] = {2: wgs72.j2, 3: wgs72.j3, 4: wgs72.j4}
"""J2, J3 e J4 do WGS-72 usados pelo SGP4 (não normalizados)."""

SGP4_Q0_ALTITUDE: Final[float] = 120e3
"""Altitude de referência q0 da atmosfera do SGP4 [m] (q0 = R⊕ + 120 km)."""

SGP4_S_ALTITUDE: Final[float] = 78e3
"""Altitude do parâmetro s da atmosfera do SGP4 [m] (s = R⊕ + 78 km), para perigeus ≥ 156 km."""

BSTAR_RHO0: Final[float] = 2.461e-5
"""ρ0 da definição B* = ½·(C_D·A/m)·ρ0 [kg/m²/ER] (Hoots & Roehrich 1980)."""


def sgp4_s_altitude(perigee_altitude_m: float | None) -> float:
    """Altitude de s [m] que o SGP4 usa para um perigeu dado (``sgp4init``).

    Perigeu ≥ 156 km: 78 km. Entre 98 e 156 km: perigeu − 78 km. Abaixo de 98 km: 20 km.
    ``None`` assume perigeu alto (78 km).
    """
    if perigee_altitude_m is None or perigee_altitude_m >= 156e3:
        return SGP4_S_ALTITUDE
    if perigee_altitude_m < 98e3:
        return 20e3
    return perigee_altitude_m - 78e3


class SGP4Atmosphere:
    """Atmosfera de lei de potência do SGP4: ρ(r) = (ρ0/R⊕)·[(q0 − s)/(r − s)]⁴ [kg/m³].

    A densidade só tem sentido junto com ``cd_a_m_from_bstar`` (ver o docstring do módulo).

    Args:
        perigee_altitude_m: altitude do perigeu [m]; muda s abaixo de 156 km como no SGP4
            (``sgp4_s_altitude``). ``None``: s = R⊕ + 78 km.
        earth_radius: R⊕ de q0 e s [m] (padrão: WGS-72, o do SGP4).

    ``density(h_m)`` recebe a altitude **esférica** h = |r| − R_EARTH do projeto (é o que o
    ``DragAcceleration`` com ``altitude_model="spherical"`` calcula) e reconstrói |r| antes de
    aplicar a lei; assim a diferença de 2 m entre os raios WGS-84 e WGS-72 não entra.
    """

    def __init__(
        self, perigee_altitude_m: float | None = None, earth_radius: float = WGS72_RE
    ) -> None:
        self.s = earth_radius + sgp4_s_altitude(perigee_altitude_m)
        self.q0 = earth_radius + SGP4_Q0_ALTITUDE
        self.rho_ref = BSTAR_RHO0 / earth_radius

    def density(self, h_m: ArrayLike) -> FloatOrArray:
        """ρ [kg/m³] na altitude esférica ``h_m`` [m] (sobre ``satsim.constants.R_EARTH``)."""
        r = np.asarray(h_m, dtype=float) + R_EARTH
        rho = self.rho_ref * ((self.q0 - self.s) / (r - self.s)) ** 4
        return float(rho) if rho.ndim == 0 else rho


def cd_a_m_from_bstar(bstar: float) -> float:
    """C_D·A/m [m²/kg] que, com ``SGP4Atmosphere``, reproduz o arrasto do B* [1/ER]: 2·B*/ρ0."""
    return 2.0 * bstar / BSTAR_RHO0


def sgp4_drag(bstar: float, perigee_altitude_m: float | None = None) -> DragAcceleration:
    """Arrasto com as premissas do SGP4: lei de potência, altitude esférica, sem co-rotação."""
    return DragAcceleration(
        cd=1.0,  # só o produto C_D·A/m importa
        area_over_mass=cd_a_m_from_bstar(bstar),
        atmosphere=SGP4Atmosphere(perigee_altitude_m),
        altitude_model="spherical",
        corotating=False,
    )


def sgp4_gravity() -> ZonalAcceleration:
    """Zonais J2–J4 com as constantes WGS-72 do SGP4 (use também ``mu=WGS72_MU`` no Cowell)."""
    return ZonalAcceleration(mu=WGS72_MU, re=WGS72_RE, coefficients=WGS72_ZONAL)


def perigee_altitude_from_tle(mean_motion_rev_day: float, eccentricity: float) -> float:
    """Altitude do perigeu [m] a partir de n (Kozai) e e do TLE, com µ e R⊕ do WGS-72.

    Aproximação (o SGP4 usa o semieixo de Brouwer, ~0,1% diferente): serve só para decidir s,
    que muda abaixo de 156 km.
    """
    n = mean_motion_rev_day * 2.0 * np.pi / 86400.0
    a = (WGS72_MU / n**2) ** (1.0 / 3.0)
    return float(a * (1.0 - eccentricity) - WGS72_RE)


def sgp4_equivalent_models(
    bstar: float, perigee_altitude_m: float | None = None
) -> tuple[ZonalAcceleration, DragAcceleration]:
    """Forças do Cowell com as premissas do SGP4: zonais WGS-72 + arrasto do B*.

    Use com ``CowellPropagator(..., mu=WGS72_MU)``. Sem Sol, Lua e SRP (o SGP4 de órbitas
    próximas não os tem).
    """
    return sgp4_gravity(), sgp4_drag(bstar, perigee_altitude_m)
