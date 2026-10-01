"""Elementos orbitais clássicos ↔ estado cartesiano ECI, e a equação de Kepler.

ECI (Earth-Centered Inertial, "centrado na Terra, inercial") é o referencial cartesiano em que
posição r e velocidade v são expressas neste módulo:

- origem no centro da Terra;
- eixo z ao longo do eixo de rotação da Terra, apontando para o polo norte;
- eixo x no plano do equador, apontando para o equinócio vernal (a direção do Sol no início da
  primavera do hemisfério norte, um ponto fixo entre as estrelas);
- eixo y completa o triedro destro (y = z × x).

"Inercial" quer dizer que os eixos não giram com a Terra: um satélite sob gravidade central
descreve neste referencial uma elipse fixa, enquanto a Terra gira por baixo dela. O referencial
que gira junto com a Terra (ECEF) é obtido do ECI por uma rotação em torno de z pelo tempo
sideral θG (semana 2). Simplificação do projeto (ADR 0001): o ECI aqui ignora precessão e
nutação do eixo terrestre.

Notação dos documentos (ADR 0001): ``psi`` = anomalia média (o "M" da literatura), ``ell`` =
momento angular específico r × v (o "h" da literatura), ``nu`` = anomalia verdadeira,
``argp`` = ω, ``raan`` = Ω, ``u`` = argp + nu (argumento de latitude). Unidades SI, ângulos em
radianos. Somente órbitas elípticas (0 ≤ e < 1).
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import MU

FloatOrArray = float | NDArray[np.float64]

ECC_TOL: float = 1e-11
"""Abaixo desta excentricidade a órbita é tratada como circular (ω indefinido)."""

INC_TOL: float = 1e-11
"""Inclinação [rad] abaixo de INC_TOL (ou acima de π − INC_TOL) é tratada como equatorial."""

_TWO_PI = 2.0 * np.pi


@dataclass(frozen=True)
class KeplerianElements:
    """Elementos clássicos (a, e, i, Ω, ω, ψ).

    a: semieixo maior [m]; e: excentricidade [-]; i: inclinação [rad]; raan: ascensão reta do
    nó ascendente Ω [rad]; argp: argumento do perigeu ω [rad]; psi: anomalia média ψ [rad].
    Os campos podem ser escalares ou arrays com shapes compatíveis por broadcasting.
    """

    a: FloatOrArray
    e: FloatOrArray
    i: FloatOrArray
    raan: FloatOrArray
    argp: FloatOrArray
    psi: FloatOrArray


def _check_eccentricity(e: NDArray[np.float64]) -> None:
    if np.any(e < 0.0) or np.any(e >= 1.0):
        raise ValueError("somente órbitas elípticas são suportadas: exige 0 <= e < 1")


def _as_scalar_if_0d(x: NDArray[np.float64]) -> FloatOrArray:
    return float(x) if x.ndim == 0 else x


def solve_kepler(
    psi: ArrayLike, e: ArrayLike, tol: float = 1e-14, max_iter: int = 50
) -> FloatOrArray:
    """Resolve a equação de Kepler E − e·sin E = ψ por Newton-Raphson vetorizado.

    ψ é reduzido a [−π, π) para a iteração. Chute inicial E0 = ψ + e·sin ψ, ou ±π (com o sinal
    de ψ) quando e > 0,8, chute que converge sempre. O resultado é devolvido na mesma volta do ψ
    de entrada, isto é, E − e·sin E = ψ vale para o ψ original e não só para o reduzido.

    Args:
        psi: anomalia média [rad], escalar ou array.
        e: excentricidade, escalar ou array compatível com ``psi``.
        tol: critério de parada sobre |ΔE| [rad].
        max_iter: número máximo de iterações.

    Returns:
        Anomalia excêntrica E [rad], com o shape de ``psi`` e ``e`` após broadcasting.

    Raises:
        ValueError: se algum e < 0 ou e ≥ 1.
        RuntimeError: se não convergir em ``max_iter`` iterações.
    """
    psi = np.asarray(psi, dtype=float)
    e = np.asarray(e, dtype=float)
    _check_eccentricity(e)
    psi, e = np.broadcast_arrays(psi, e)
    psi_red = np.mod(psi + np.pi, _TWO_PI) - np.pi
    ecc_start = np.where(psi_red >= 0.0, np.pi, -np.pi)
    big_e = psi_red + e * np.sin(psi_red)
    big_e = np.where(e > 0.8, ecc_start, big_e)
    for _ in range(max_iter):
        delta = (big_e - e * np.sin(big_e) - psi_red) / (1.0 - e * np.cos(big_e))
        big_e = big_e - delta
        if np.all(np.abs(delta) < tol):
            break
    else:
        raise RuntimeError(f"equação de Kepler não convergiu em {max_iter} iterações")
    return _as_scalar_if_0d(big_e + (psi - psi_red))


def eccentric_to_true(big_e: ArrayLike, e: ArrayLike) -> FloatOrArray:
    """Anomalia excêntrica E → anomalia verdadeira ν [rad], em (−π, π]."""
    big_e = np.asarray(big_e, dtype=float)
    e = np.asarray(e, dtype=float)
    _check_eccentricity(e)
    nu = np.arctan2(np.sqrt(1.0 - e**2) * np.sin(big_e), np.cos(big_e) - e)
    return _as_scalar_if_0d(nu)


def true_to_eccentric(nu: ArrayLike, e: ArrayLike) -> FloatOrArray:
    """Anomalia verdadeira ν → anomalia excêntrica E [rad], em (−π, π]."""
    nu = np.asarray(nu, dtype=float)
    e = np.asarray(e, dtype=float)
    _check_eccentricity(e)
    big_e = np.arctan2(np.sqrt(1.0 - e**2) * np.sin(nu), e + np.cos(nu))
    return _as_scalar_if_0d(big_e)


def mean_to_true(psi: ArrayLike, e: ArrayLike) -> FloatOrArray:
    """Anomalia média ψ → anomalia verdadeira ν [rad], em (−π, π]."""
    return eccentric_to_true(solve_kepler(psi, e), e)


def true_to_mean(nu: ArrayLike, e: ArrayLike) -> FloatOrArray:
    """Anomalia verdadeira ν → anomalia média ψ = E − e·sin E [rad], em (−π, π]."""
    big_e = np.asarray(true_to_eccentric(nu, e))
    return _as_scalar_if_0d(big_e - np.asarray(e, dtype=float) * np.sin(big_e))


def elements_to_state(
    el: KeplerianElements, mu: float = MU
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Elementos clássicos → estado (r, v) em ECI.

    No plano perifocal: r_pf = (a(cos E − e), a√(1−e²) sin E, 0) e
    v_pf = (√(µa)/r)·(−sin E, √(1−e²) cos E, 0), com r = a(1 − e cos E). A rotação para ECI é a
    sequência 3-1-3 Rz(Ω)·Rx(i)·Rz(ω), aplicada pelas colunas P e Q da matriz.

    Args:
        el: elementos; ``psi`` (e os demais campos) podem ser arrays de shape (N,).
        mu: parâmetro gravitacional [m³/s²].

    Returns:
        (r, v) em ECI [m, m/s], cada um com shape (N, 3); N = 1 para entrada escalar.
    """
    a, e, inc, raan, argp, psi = np.broadcast_arrays(
        *(
            np.atleast_1d(np.asarray(x, dtype=float))
            for x in (el.a, el.e, el.i, el.raan, el.argp, el.psi)
        )
    )
    big_e = np.asarray(solve_kepler(psi, e))
    cos_e, sin_e = np.cos(big_e), np.sin(big_e)
    sqrt_1me2 = np.sqrt(1.0 - e**2)
    r_norm = a * (1.0 - e * cos_e)

    x_pf = a * (cos_e - e)
    y_pf = a * sqrt_1me2 * sin_e
    v_scale = np.sqrt(mu * a) / r_norm
    vx_pf = -v_scale * sin_e
    vy_pf = v_scale * sqrt_1me2 * cos_e

    co, so = np.cos(raan), np.sin(raan)
    cw, sw = np.cos(argp), np.sin(argp)
    ci, si = np.cos(inc), np.sin(inc)
    p_hat = np.stack([co * cw - so * sw * ci, so * cw + co * sw * ci, sw * si], axis=-1)
    q_hat = np.stack([-co * sw - so * cw * ci, -so * sw + co * cw * ci, cw * si], axis=-1)

    r = x_pf[:, None] * p_hat + y_pf[:, None] * q_hat
    v = vx_pf[:, None] * p_hat + vy_pf[:, None] * q_hat
    return r, v


def _unit(x: NDArray[np.float64]) -> NDArray[np.float64]:
    norm = np.linalg.norm(x, axis=-1, keepdims=True)
    return x / np.where(norm > 0.0, norm, 1.0)


def _signed_angle(
    a_hat: NDArray[np.float64], b_hat: NDArray[np.float64], axis_hat: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Ângulo de a → b em torno de ``axis_hat`` (sentido do movimento), em [0, 2π)."""
    sin_part = np.sum(axis_hat * np.cross(a_hat, b_hat), axis=-1)
    cos_part = np.sum(a_hat * b_hat, axis=-1)
    return np.mod(np.arctan2(sin_part, cos_part), _TWO_PI)


def state_to_elements(r: ArrayLike, v: ArrayLike, mu: float = MU) -> KeplerianElements:
    """Estado (r, v) em ECI → elementos clássicos.

    Usa ℓ = r × v, o vetor dos nós n = ẑ × ℓ e o vetor excentricidade
    e_vec = ((v² − µ/r) r − (r·v) v)/µ; a pela energia ε = v²/2 − µ/r = −µ/(2a). Os ângulos são
    obtidos por arctan2, medidos no sentido do movimento (em torno de ℓ̂): ω < π ⇔ e_vec_z > 0 e
    ν < π ⇔ r·v > 0. A inclinação usa arctan2(√(ℓx² + ℓy²), ℓz), equivalente a arccos(ℓz/|ℓ|)
    porém bem-condicionada perto de i = 0 e i = π.

    Casos degenerados (Sem4 §3.3: a singularidade é da carta de coordenadas, não da física):

    - circular (e < ECC_TOL): sem perigeu; argp = 0 e ν passa a ser o argumento de latitude u,
      medido a partir do nó ascendente;
    - equatorial (i < INC_TOL ou i > π − INC_TOL): sem linha dos nós; raan = 0 e a linha dos
      nós é tomada no eixo x, de modo que argp vira a longitude do perigeu;
    - circular e equatorial: raan = argp = 0 e ν é a longitude verdadeira (a partir do eixo x).

    Nesses casos ψ = ν (e ≈ 0) e a ida e volta pelo estado é preservada.

    Args:
        r: posição ECI [m], shape (3,) ou (N, 3).
        v: velocidade ECI [m/s], mesmo shape de ``r``.
        mu: parâmetro gravitacional [m³/s²].

    Returns:
        KeplerianElements; campos escalares para entrada (3,), arrays (N,) para (N, 3).
        raan, argp e psi em [0, 2π).

    Raises:
        ValueError: se a órbita não for elíptica (e ≥ 1 ou energia ≥ 0).
    """
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    # Um único estado (3,) é promovido a (1, 3) para que todo o cálculo abaixo seja o mesmo
    # código vetorizado; ``single`` lembra de devolver escalares no final.
    single = r.ndim == 1
    r = np.atleast_2d(r)
    v = np.atleast_2d(v)

    # Grandezas linha a linha (axis=-1): |r|, v² e o produto escalar r·v de cada estado.
    r_norm = np.linalg.norm(r, axis=-1)
    v2 = np.sum(v * v, axis=-1)
    rv = np.sum(r * v, axis=-1)
    # ℓ = r × v é perpendicular ao plano orbital (Sem4 §1): sua direção ℓ̂ define o plano e o
    # sentido do movimento; é o eixo em torno do qual medimos Ω→ω→ν mais abaixo.
    ell = np.cross(r, v)
    ell_hat = _unit(ell)

    # Energia específica ε = v²/2 − µ/r, constante sob gravidade central (Sem3 §1).
    energy = 0.5 * v2 - mu / r_norm
    # Vetor excentricidade (vetor de Laplace-Runge-Lenz / µ): aponta do centro da Terra para o
    # perigeu e seu módulo é e. [:, None] transforma os escalares (N,) em (N, 1) para que
    # multipliquem cada vetor (N, 3) linha a linha.
    e_vec = ((v2 - mu / r_norm)[:, None] * r - rv[:, None] * v) / mu
    e = np.linalg.norm(e_vec, axis=-1)
    # ε ≥ 0 ou e ≥ 1 significam trajetória parabólica/hiperbólica (escape), fora do escopo.
    if np.any(energy >= 0.0) or np.any(e >= 1.0):
        raise ValueError("somente órbitas elípticas são suportadas: exige 0 <= e < 1")
    # Da relação ε = −µ/(2a), que só depende de a (Sem3 §3.1).
    a = -mu / (2.0 * energy)
    # i é o ângulo entre ℓ e o eixo z. Mesma coisa que arccos(ℓz/|ℓ|), mas o arccos perde
    # precisão perto de 0 e π (daria i ≈ 1e-8 rad numa órbita equatorial exata).
    inc = np.arctan2(np.hypot(ell[:, 0], ell[:, 1]), ell[:, 2])

    # Máscaras booleanas (N,) dos casos degenerados (Sem4 §3.3), em que ω ou Ω não existem.
    circular = e < ECC_TOL
    equatorial = (inc < INC_TOL) | (inc > np.pi - INC_TOL)

    # Direção de referência no plano orbital a partir da qual ω é medido: a linha dos nós
    # n = ẑ × ℓ = (−ℓy, ℓx, 0), que aponta para o nó ascendente. Em órbita equatorial n = 0 (o
    # plano não cruza o equador), então usamos o eixo x por convenção, o que faz Ω = 0.
    x_hat = np.broadcast_to(np.array([1.0, 0.0, 0.0]), r.shape)
    node = np.stack([-ell[:, 1], ell[:, 0], np.zeros_like(r_norm)], axis=-1)
    node_hat = np.where(equatorial[:, None], x_hat, _unit(node))
    # Direção do perigeu, a partir da qual ν é medido. Em órbita circular não há perigeu
    # (e_vec ≈ 0, direção aleatória); tomamos o próprio nó, o que faz ω = 0 e ν = u.
    peri_hat = np.where(circular[:, None], node_hat, _unit(e_vec))

    # Ω = azimute do nó no plano equatorial, medido a partir do eixo x.
    raan = np.mod(np.arctan2(node_hat[:, 1], node_hat[:, 0]), _TWO_PI)
    # ω e ν são ângulos com sinal medidos em torno de ℓ̂ (sentido do movimento): o seno vem de
    # ℓ̂·(a × b) e o cosseno de a·b, e arctan2 resolve o quadrante. Isso equivale às regras
    # "ω < π se e_z > 0" e "ν < π se r·v > 0", sem casos especiais.
    argp = _signed_angle(node_hat, peri_hat, ell_hat)
    nu = _signed_angle(peri_hat, _unit(r), ell_hat)
    # ν → ψ pela equação de Kepler (via anomalia excêntrica), reduzido a [0, 2π).
    psi = np.mod(np.asarray(true_to_mean(nu, e)), _TWO_PI)

    fields = (a, e, inc, raan, argp, psi)
    if single:
        return KeplerianElements(*(float(x[0]) for x in fields))
    return KeplerianElements(*fields)
