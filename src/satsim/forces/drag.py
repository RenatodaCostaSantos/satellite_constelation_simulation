"""Aceleração de arrasto atmosférico (Sem3 §2, eqs. 5–6).

A força de arrasto é F_D = ½·C_D·ρ·A·v_rel² (Sem3 eq. 4), oposta à velocidade em relação ao ar.
Dividindo pela massa (Sem3 eq. 6):

    a_D = −½ · C_D · (A/m) · ρ(h) · |v_rel| · v_rel,

com a forma |v_rel|·v_rel, numericamente estável, em vez do versor v_rel/|v_rel|. A termosfera
gira com a Terra, então a velocidade relativa ao ar é (Sem3 eq. 5)

    v_rel = v − ω⊕ × r,   ω⊕ = OMEGA_EARTH·ẑ.

No equador |ω⊕ × r| ≈ 465 m/s contra v ≈ 7670 m/s: até ~6% na velocidade, ~12% na força.

O parâmetro que governa o efeito é o produto C_D·A/m (o inverso do coeficiente balístico, Sem3
eq. 7): dois objetos com o mesmo C_D·A/m decaem de forma idêntica.

Altitude: geodésica × esférica
------------------------------
``altitude_model="geodetic"`` usa a altitude sobre o elipsoide WGS84,

    h = p·cos φ + z·sin φ − a·√(1 − e²·sin²φ),   p = √(x² + y²),

com φ a latitude geodésica. Ela só depende de p e z, que a rotação em torno de z não altera,
então pode ser calculada direto do r em ECI, sem GMST nem t. ``"spherical"`` usa h = |r| − R⊕,
como os PDFs. A geodésica é sempre ≥ a esférica (o elipsoide fica dentro da esfera de raio R⊕):
a diferença vai de 0 no equador a 21,4 km nos polos (≈ 21,4 km·sin²φ), o que muda ρ em até ~40%
numa órbita polar (a 400 km, H ≈ 58 km: e^(21/58) ≈ 1,44). Para a ISS (|φ| ≤ 51,6°) chega a
≈ 13 km, ≈ 25% em ρ.
"""

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import OMEGA_EARTH, R_EARTH
from satsim.forces.atmosphere import ExponentialAtmosphere
from satsim.geodesy import ecef_to_geodetic

AltitudeModel = Literal["geodetic", "spherical"]


def altitude(r: ArrayLike, model: AltitudeModel = "geodetic") -> NDArray[np.float64]:
    """Altitude [m] a partir de r em ECI (ou ECEF): geodésica WGS84 ou esférica |r| − R⊕.

    Args:
        r: posição [m], shape (3,) ou (N, 3).
        model: ``"geodetic"`` ou ``"spherical"``.

    Returns:
        Altitude com shape () ou (N,).

    Raises:
        ValueError: modelo desconhecido.
    """
    r = np.asarray(r, dtype=float)
    if model == "geodetic":
        return np.asarray(ecef_to_geodetic(r)[2])
    if model == "spherical":
        return np.linalg.norm(r, axis=-1) - R_EARTH
    raise ValueError(f"altitude_model deve ser 'geodetic' ou 'spherical', recebido {model!r}")


class DragAcceleration:
    """Arrasto atmosférico (implementa ``AccelerationModel``), Sem3 eq. 6.

    Args:
        cd: coeficiente de arrasto C_D (≈ 2,2 em LEO, Sem3 §2.2).
        area_over_mass: A/m [m²/kg] (ISS ≈ 0,004, Sem3 §2.4).
        atmosphere: modelo de densidade com ``density(h_m)``; padrão ``ExponentialAtmosphere()``.
        altitude_model: ``"geodetic"`` (WGS84) ou ``"spherical"`` (|r| − R⊕, como nos PDFs).
        corotating: ``True`` usa v_rel = v − ω⊕ × r (Sem3 eq. 5); ``False`` usa v_rel = v.

    Raises:
        ValueError: C_D ou A/m negativos, ou ``altitude_model`` desconhecido.
    """

    def __init__(
        self,
        cd: float = 2.2,
        area_over_mass: float = 0.004,
        atmosphere: ExponentialAtmosphere | None = None,
        altitude_model: AltitudeModel = "geodetic",
        corotating: bool = True,
    ) -> None:
        if cd < 0.0 or area_over_mass < 0.0:
            raise ValueError("cd e area_over_mass não podem ser negativos")
        if altitude_model not in ("geodetic", "spherical"):
            raise ValueError(
                f"altitude_model deve ser 'geodetic' ou 'spherical', recebido {altitude_model!r}"
            )
        self.cd = cd
        self.area_over_mass = area_over_mass
        self.atmosphere = atmosphere if atmosphere is not None else ExponentialAtmosphere()
        self.altitude_model: AltitudeModel = altitude_model
        self.corotating = corotating

    @property
    def cd_a_m(self) -> float:
        """C_D·A/m [m²/kg]: o único parâmetro que a trajetória enxerga (Sem3 §2.4)."""
        return self.cd * self.area_over_mass

    def relative_velocity(self, r: ArrayLike, v: ArrayLike) -> NDArray[np.float64]:
        """v_rel = v − ω⊕ × r (Sem3 eq. 5), ou v se ``corotating=False``. Shape de ``v``."""
        v = np.asarray(v, dtype=float)
        if not self.corotating:
            return v
        r = np.asarray(r, dtype=float)
        # ω⊕ × r com ω⊕ = (0, 0, OMEGA_EARTH): (−ω·y, ω·x, 0)
        omega_cross_r = np.stack(
            [-OMEGA_EARTH * r[..., 1], OMEGA_EARTH * r[..., 0], np.zeros_like(r[..., 2])], axis=-1
        )
        return v - omega_cross_r

    def __call__(
        self, t: float, r: NDArray[np.float64], v: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """a_D = −½·C_D·(A/m)·ρ(h)·|v_rel|·v_rel [m/s²] (Sem3 eq. 6); aceita (3,) ou (N, 3)."""
        v_rel = self.relative_velocity(r, v)
        rho = self.atmosphere.density(altitude(r, self.altitude_model))
        speed = np.linalg.norm(v_rel, axis=-1)
        return (-0.5 * self.cd_a_m * np.asarray(rho) * speed)[..., np.newaxis] * v_rel
