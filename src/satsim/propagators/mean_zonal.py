"""Propagador de elementos médios com J2, J3 e J4 de 1ª ordem e arrasto médio (card S4-05).

Integra as EDOs dos elementos médios (a, e, i, Ω, ω, ψ) com ``solve_ivp``:

    ė, i̇, Ω̇, ω̇, ψ̇ = ``secular.zonal_mean_rates`` (Apêndice A.5)
    ȧ = −C_D·(A/m)·⟨ρ⟩·√(µa)          (Sem3 eq. 14; 0 sem arrasto)

Diferente do ``MeanJ2Propagator``, as taxas não são constantes: com J3 e J4, e, i e ω evoluem
acoplados (o termo de J3 em ω̇ vai com 1/e, e ė depende de cos ω), e o arrasto reduz a, o que
acelera todas as taxas. Por isso a integração numérica em vez de fórmulas fechadas; as EDOs são
suaves e lentas, e o integrador dá passos de horas.

⟨ρ⟩ é a densidade média ao longo da órbita média: amostrada em ``DENSITY_SAMPLES`` pontos
uniformes em anomalia média (média temporal) na altitude do modelo de arrasto (geodésica por
padrão). Numa SSO, a altitude geodésica varia ≈ 21 km entre o equador e os polos, e ⟨ρ⟩ ao longo
da órbita difere bastante da densidade em h = a − R⊕ (Apêndice A.7). A co-rotação da atmosfera não
entra na eq. 14: com ``corotating=True`` no Cowell o arrasto real difere em ~1% numa SSO.

O modelo médio é uma ferramenta de **projeto**, não de posição: o resíduo de 2ª ordem (J2²) em
Ω̇ e no período nodal, de poucos ppm, acumula ≈ 92 km ao longo da órbita em 30 dias na SSO de
projeto (Apêndice A.7). As taxas médias e o período nodal, esses sim, ficam a ≤ 0,1% e ≤ 10 ppm
do Cowell.

Limitação: com J3 ≠ 0 as taxas divergem em e → 0 (ver ``satsim.secular``). Órbitas de projeto
congeladas (e ≈ 1e-3) ficam longe disso; uma órbita que parta de e = 0 não é aceita.
"""

import math
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import OdeSolution, solve_ivp

from satsim.constants import J2, J3, J4, MU, R_EARTH
from satsim.elements import KeplerianElements, elements_to_state
from satsim.forces.drag import DragAcceleration, altitude
from satsim.propagators.base import Propagator, as_time_array
from satsim.secular import ZonalMeanRates, zonal_mean_rates

DENSITY_SAMPLES = 72
"""Pontos por revolução na média ⟨ρ⟩ ao longo da órbita (passo de 5° em anomalia média)."""


class MeanZonalPropagator(Propagator):
    """Propaga elementos **médios** sob J2, J3, J4 (1ª ordem) e arrasto médio opcional.

    Args:
        mean_elements: elementos médios na época (campos escalares; ``psi`` = ψ0).
        epoch: instante t = 0, ``datetime`` tz-aware.
        j2, j3, j4: harmônicos zonais; ``j3 = j4 = 0`` reproduz o ``MeanJ2Propagator``.
        drag: modelo de arrasto do Cowell (usa C_D·A/m, a atmosfera e o modelo de altitude dele)
            ou ``None`` (sem arrasto).
        mu: parâmetro gravitacional [m³/s²].
        re: raio equatorial de referência [m].
        rtol: tolerância relativa do integrador.

    Raises:
        ValueError: ``epoch`` naive, a ≤ 0, e fora de [0, 1) ou e = 0 com ``j3 ≠ 0``.
    """

    def __init__(
        self,
        mean_elements: KeplerianElements,
        epoch: datetime,
        j2: float = J2,
        j3: float = J3,
        j4: float = J4,
        drag: DragAcceleration | None = None,
        mu: float = MU,
        re: float = R_EARTH,
        rtol: float = 1e-12,
    ) -> None:
        if epoch.tzinfo is None or epoch.utcoffset() is None:
            raise ValueError("epoch deve ser um datetime tz-aware (ex.: tzinfo=timezone.utc)")
        self._elements = mean_elements
        self._epoch = epoch
        self._j = (j2, j3, j4)
        self._drag = drag
        self._mu = mu
        self._re = re
        self._rtol = rtol
        el = mean_elements
        self._y0 = np.array([el.a, el.e, el.i, el.raan, el.argp, el.psi], dtype=float)
        self._rates0 = self.rates(el.a, el.e, el.i, el.argp)  # valida a, e, i, ω
        # Saída densa para t ≥ 0 e para t ≤ 0, cada uma de 0 até o extremo já integrado.
        self._forward: tuple[float, OdeSolution] | None = None
        self._backward: tuple[float, OdeSolution] | None = None
        self._psi_samples = np.linspace(0.0, 2.0 * np.pi, DENSITY_SAMPLES, endpoint=False)

    @property
    def epoch(self) -> datetime:
        """Época (t = 0), tz-aware."""
        return self._epoch

    @property
    def elements(self) -> KeplerianElements:
        """Elementos médios na época."""
        return self._elements

    @property
    def period(self) -> float:
        """Período nodal na época, Tnod = 2π/(ψ̇ + ω̇) [s] com as taxas J2–J4 (Sem4 eq. 34)."""
        return 2.0 * math.pi / (self._rates0.psi_dot + self._rates0.argp_dot)

    def rates(self, a: float, e: float, i: float, argp: float) -> ZonalMeanRates:
        """Taxas médias zonais (``zonal_mean_rates``) com os J e o µ deste propagador."""
        j2, j3, j4 = self._j
        return zonal_mean_rates(a, e, i, argp, mu=self._mu, j2=j2, j3=j3, j4=j4, re=self._re)

    def mean_density(self, a: float, e: float, i: float, argp: float) -> float:
        """⟨ρ⟩ [kg/m³] ao longo da órbita média, uniforme em anomalia média (média temporal).

        Ω não entra: a altitude (geodésica ou esférica) não muda com a rotação em torno de z.
        """
        if self._drag is None:
            return 0.0
        k = len(self._psi_samples)
        el = KeplerianElements(
            a=np.full(k, a),
            e=np.full(k, e),
            i=np.full(k, i),
            raan=np.zeros(k),
            argp=np.full(k, argp),
            psi=self._psi_samples,
        )
        r, _ = elements_to_state(el, self._mu)
        h = altitude(r, self._drag.altitude_model)
        return float(np.mean(self._drag.atmosphere.density(h)))

    def semi_major_axis_rate(self, a: float, e: float, i: float, argp: float) -> float:
        """ȧ = −C_D·(A/m)·⟨ρ⟩·√(µa) [m/s] (Sem3 eq. 14); 0 sem arrasto."""
        if self._drag is None:
            return 0.0
        return -self._drag.cd_a_m * self.mean_density(a, e, i, argp) * math.sqrt(self._mu * a)

    def _rhs(self, t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        a, e, i, _, w, _ = y
        r = self.rates(a, e, i, w)
        a_dot = self.semi_major_axis_rate(a, e, i, w)
        return np.array([a_dot, r.e_dot, r.i_dot, r.raan_dot, r.argp_dot, r.psi_dot])

    def _solution(self, t_end: float) -> OdeSolution:
        """Saída densa cobrindo [0, t_end] (ou [t_end, 0]); reintegra se o trecho não basta."""
        attr = "_forward" if t_end >= 0.0 else "_backward"
        cached = getattr(self, attr)
        if cached is not None and abs(cached[0]) >= abs(t_end):
            return cached[1]
        sol = solve_ivp(
            self._rhs,
            (0.0, t_end),
            self._y0,
            method="DOP853",
            rtol=self._rtol,
            atol=1e-12 * np.maximum(np.abs(self._y0), 1.0),
            dense_output=True,
        )
        if not sol.success:
            raise RuntimeError(f"integração dos elementos médios falhou: {sol.message}")
        setattr(self, attr, (t_end, sol.sol))
        return sol.sol

    def elements_at(self, t: ArrayLike) -> KeplerianElements:
        """Elementos médios em t segundos desde ``epoch`` (t pode ser negativo).

        Returns:
            KeplerianElements com campos de shape (N,); Ω, ω e ψ crescem continuamente.
        """
        t = as_time_array(t)
        y = np.empty((6, t.size))
        for mask, t_end in ((t > 0.0, t.max()), (t < 0.0, t.min())):
            if np.any(mask):
                y[:, mask] = self._solution(float(t_end))(t[mask])
        y[:, t == 0.0] = self._y0[:, None]
        return KeplerianElements(*y)

    def propagate(self, t: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estado ECI da trajetória **média** em t segundos desde ``epoch``.

        Returns:
            (r, v) em ECI [m, m/s], cada um com shape (N, 3).
        """
        return elements_to_state(self.elements_at(t), self._mu)
