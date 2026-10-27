"""Gravidade zonal J2–J4 para o Cowell: uma única implementação para todos os harmônicos zonais.

Potencial perturbador dos zonais (m = 0, J_n = −C_n0), Sem4 eq. 6:

    R(r) = Σ_n −(µ/r)·J_n·(R⊕/r)^n·P_n(s),     s = z/r = sin φ (latitude geocêntrica)

e a aceleração é o seu gradiente (Apêndice A.2 / card S4-01 do prompt da Semana 4). Usando
∇r = r̂ e ∇s = (ẑ − s·r̂)/r:

    a = ∇R = Σ_n −µ·J_n·R⊕^n·r^−(n+2) · [ P_n'(s)·(ẑ − s·r̂) − (n+1)·P_n(s)·r̂ ]

com os polinômios de Legendre

    P2 = (3s² − 1)/2           P2' = 3s
    P3 = (5s³ − 3s)/2          P3' = (15s² − 3)/2
    P4 = (35s⁴ − 30s² + 3)/8   P4' = (140s³ − 60s)/8

calculados pela recorrência de Bonnet, (n+1)·P_{n+1} = (2n+1)·s·P_n − n·P_{n−1}, e
P'_{n+1} = P'_{n−1} + (2n+1)·P_n (válidas para qualquer grau).

Para n = 2 a fórmula reproduz exatamente ``j2_acceleration`` (Sem4 §9.3); o teste de
consistência confere em 200 pontos (rel ~1e-15).

O papel de cada termo na SSO de projeto:

- **J2** (par): o achatamento; domina a precessão do nó Ω̇ (Sem4 eq. 23) e do perigeu ω̇.
- **J3** (ímpar): a assimetria norte–sul, a Terra "em forma de pera". Por ser ímpar em s, não
  muda Ω̇ em média (a contribuição é proporcional a e·sin ω, ~1e-3 da de J2); seu efeito
  secular principal é um ė ∝ cos ω que, somado ao ω̇ de J2, tem um ponto fixo em ω = 90° e
  e = e_cong: é o termo que **congela o perigeu** (Sem4 §6.3, card S4-05).
- **J4** (par): reduz Ω̇ da SSO de projeto em 0,23% (Apêndice A.2: Ω̇(J2–J4)/Ω̇(J2) = 0,997703),
  o bastante para o LTAN derivar minutos por ano se a inclinação vier da teoria só J2.

A energia conservada sob gravidade central + zonais é ε = v²/2 − µ/r − R(r).

Tesserais (J22 e seguintes) ficam fora: a média sobre a rotação da Terra os torna pequenos em LEO.
"""

from collections.abc import Iterable, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import J2, J3, J4, MU, R_EARTH

ZONAL_COEFFICIENTS: Mapping[int, float] = {2: J2, 3: J3, 4: J4}
"""Harmônicos zonais não normalizados disponíveis por padrão (``satsim.constants``)."""


def _check_degrees(degrees: Iterable[int], coefficients: Mapping[int, float]) -> tuple[int, ...]:
    degrees = tuple(int(n) for n in degrees)
    for n in degrees:
        if n < 2:
            raise ValueError(f"grau zonal deve ser ≥ 2 (recebido {n})")
        if n not in coefficients:
            raise ValueError(f"sem coeficiente J{n} (disponíveis: {sorted(coefficients)})")
    return degrees


def _legendre(s: NDArray[np.float64], n_max: int) -> tuple[list, list]:
    """P_n(s) e P_n'(s) para n = 0 … n_max pela recorrência de Bonnet."""
    # Cada item da lista é um nparray. O append adiciona um novo np.array nas listas p e dp.
    p = [np.ones_like(s), s]
    dp = [np.zeros_like(s), np.ones_like(s)]
    for n in range(1, n_max):
        p.append(((2 * n + 1) * s * p[n] - n * p[n - 1]) / (n + 1))
        dp.append(dp[n - 1] + (2 * n + 1) * p[n])
    return p, dp


def zonal_potential(
    r: ArrayLike,
    degrees: Iterable[int] = (2, 3, 4),
    mu: float = MU,
    re: float = R_EARTH,
    coefficients: Mapping[int, float] = ZONAL_COEFFICIENTS,
) -> NDArray[np.float64]:
    """Função perturbadora zonal R(r) = Σ_n −(µ/r)·J_n·(R⊕/r)^n·P_n(z/r) [m²/s²] (Sem4 eq. 6).

    Args:
        r: posição ECI [m], shape (3,) ou (N, 3).
        degrees: graus zonais a somar (cada um com coeficiente em ``coefficients``).
        mu: parâmetro gravitacional [m³/s²].
        re: raio equatorial de referência dos harmônicos [m].
        coefficients: mapa grau → J_n (não normalizado).

    Returns:
        R com shape () para (3,) ou (N,) para (N, 3). A aceleração é ∇R; a energia conservada é
        v²/2 − µ/r − R.
    """
    degrees = _check_degrees(degrees, coefficients)
    r = np.asarray(r, dtype=float)
    r_norm = np.sqrt(np.sum(r * r, axis=-1))
    s = r[..., 2] / r_norm
    p, _ = _legendre(s, max(degrees, default=1))
    total = np.zeros_like(r_norm)
    for n in degrees:
        total = total - mu / r_norm * coefficients[n] * (re / r_norm) ** n * p[n]
    return total


def zonal_acceleration(
    r: ArrayLike,
    degrees: Iterable[int] = (2, 3, 4),
    mu: float = MU,
    re: float = R_EARTH,
    coefficients: Mapping[int, float] = ZONAL_COEFFICIENTS,
) -> NDArray[np.float64]:
    """Aceleração zonal a = ∇R em ECI (Apêndice A.2 / card S4-01 do prompt da Semana 4).

    a = Σ_n −µ·J_n·R⊕^n·r^−(n+2)·[P_n'(s)·(ẑ − s·r̂) − (n+1)·P_n(s)·r̂],  s = z/r.

    Args:
        r: posição ECI [m], shape (3,) ou (N, 3).
        degrees, mu, re, coefficients: como em ``zonal_potential``.

    Returns:
        Aceleração [m/s²] com o mesmo shape de ``r``.
    """
    degrees = _check_degrees(degrees, coefficients)
    r = np.asarray(r, dtype=float)
    r_norm = np.sqrt(np.sum(r * r, axis=-1))
    s = r[..., 2] / r_norm
    p, dp = _legendre(s, max(degrees, default=1))
    # a = radial·r̂ + axial·ẑ, com ẑ − s·r̂ = r·∇s:
    #   radial = Σ c_n·[−s·P_n' − (n+1)·P_n],   axial = Σ c_n·P_n',   c_n = −µ·J_n·R⊕^n·r^−(n+2)
    radial = np.zeros_like(r_norm)
    axial = np.zeros_like(r_norm)
    for n in degrees:
        c = -mu * coefficients[n] * re**n * r_norm ** (-(n + 2))
        radial = radial - c * (s * dp[n] + (n + 1) * p[n])
        axial = axial + c * dp[n]
    acc = (radial / r_norm)[..., None] * r
    acc[..., 2] += axial
    return acc


class ZonalAcceleration:
    """Gravidade zonal J2–J4 para o Cowell (implementa ``AccelerationModel``).

    Args:
        degrees: graus zonais incluídos, ex.: ``(2,)``, ``(2, 3)``, ``(2, 3, 4)`` (padrão).
        mu: parâmetro gravitacional [m³/s²].
        re: raio equatorial de referência dos harmônicos [m].
        coefficients: mapa grau → J_n (padrão: J2, J3, J4 de ``satsim.constants``).

    Raises:
        ValueError: se algum grau for < 2 ou não tiver coeficiente.
    """

    def __init__(
        self,
        degrees: Iterable[int] = (2, 3, 4),
        mu: float = MU,
        re: float = R_EARTH,
        coefficients: Mapping[int, float] = ZONAL_COEFFICIENTS,
    ) -> None:
        self.degrees = _check_degrees(degrees, coefficients)
        self.mu = mu
        self.re = re
        self.coefficients = dict(coefficients)

    def __call__(
        self, t: float, r: NDArray[np.float64], v: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        return zonal_acceleration(r, self.degrees, self.mu, self.re, self.coefficients)

    def potential(self, r: ArrayLike) -> NDArray[np.float64]:
        """Função perturbadora R(r) com os parâmetros deste modelo (ver ``zonal_potential``)."""
        return zonal_potential(r, self.degrees, self.mu, self.re, self.coefficients)
