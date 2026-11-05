"""Manutenção da constelação: deriva do traço e do LTAN, Δv e próxima manobra (card S5-05).

Uma SSO de repetição sai do lugar por dois motivos independentes, que custam Δv parecido
(≈ 4,8 m/s por ano na órbita 91/6) por razões diferentes (Sem4 §6 e §7.6):

1. **Traço (arrasto).** O arrasto baixa o semieixo a uma taxa ȧ < 0; o período encurta e o traço
   se adianta para leste. Partindo da órbita nominal, o desvio no equador cresce com t²
   (Sem4 §7.6, a partir da eq. 40):

       λ(t) = 3π·|ȧ|·t² / (2·a·T_dia)          [rad; ȧ em m/s, t em s, T_dia = 86 400 s]

   (no Sem4, com ȧ em m/dia e t em dias, λ = 3π|ȧ|t²/(2a)). A banda morta ±L é atingida em
   t = √(2·a·T_dia·L_rad/(3π|ȧ|)) (eq. 41); a manobra devolve o δa = |ȧ|·t perdido, com
   Δv = v·δa/(2a) (Gauss, órbita circular). Por ano: Δv = v·|ȧ|·ano/(2a). É a manutenção
   dominante: uma manobra por semana.
2. **LTAN (Sol).** O Sol inclina o plano (S4-02: di/dt ≈ −9,8e-5 °/dia na SSO de projeto). Como
   Ω̇ ∝ cos i, um erro de inclinação δi muda a precessão por δΩ̇ = −Ω̇·tan i·δi (Sem4 eq. 32), e o
   LTAN deriva a 4 min por grau de Ω. Com δi(t) = (di/dt)·t, a deriva cresce linearmente e o LTAN
   se afasta de forma quadrática: ΔLTAN(t) = ½·k·(di/dt)·t², k = dLTAN/dt por δi. Corrigir a
   inclinação custa Δv = v·|δi| (manobra fora do plano); por ano, v·|di/dt|·ano. Uma manobra por
   ano.

``measure_decay_rate`` mede ȧ no Cowell por médias em janelas nodais (``window_mean_elements``),
e ``maintenance_status`` junta tudo num ``MaintenanceStatus`` por satélite, pronto para a cena.
As manobras não são simuladas: o indicador prevê quando e quanto.
"""

import math
from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from satsim.constants import MU, OMEGA_SUN, R_EARTH, SECONDS_PER_DAY, T_TROP_DAYS
from satsim.groundtrack import node_crossings
from satsim.mean_osc import window_mean_elements
from satsim.satellite import Satellite

FloatOrArray = float | NDArray[np.float64]

YEAR_S = T_TROP_DAYS * SECONDS_PER_DAY
"""Ano trópico [s] (365,2422 dias), o mesmo de ``orbit_design.ltan_drift_min_per_year``."""

SUN_INCLINATION_DRIFT_DEG_DAY = -9.8e-5
"""di/dt causado pelo Sol na SSO de projeto [°/dia] (S4-02; Apêndice A.7 do prompt da Semana 5)."""

MINUTES_PER_RAD_OF_RAAN = 24.0 * 60.0 / (2.0 * math.pi)
"""O LTAN anda 24 h por volta completa de Ω relativa ao Sol: 4 min por grau (Sem4 eq. 28)."""


def _scalar(x: NDArray[np.float64]) -> FloatOrArray:
    return float(x) if np.ndim(x) == 0 else x


def circular_velocity(a: float) -> float:
    """Velocidade de órbita circular v = √(µ/a) [m/s] (7 607,1 m/s na órbita de projeto)."""
    return math.sqrt(MU / a)


# --- Traço: arrasto -------------------------------------------------------------------------------


def track_drift(t: ArrayLike, a_dot: float, a: float) -> tuple[FloatOrArray, FloatOrArray]:
    """Deriva do traço no equador partindo da órbita nominal: λ(t) = 3π|ȧ|t²/(2a·T_dia).

    Sem4 §7.6, integrando a eq. 40 com δa(t) = ȧ·t: a taxa de deriva é dλ/dt = 3π|δa|/(a·T_dia)
    (o fechamento de um ciclo de D dias, 3πD·δa/a, repartido pelos dias), e cresce com o erro
    acumulado em a. Subestima em ≈ 0,5% (o efeito de ȧ em Ω̇ fica de fora, Sem4 §7.6).

    Args:
        t: tempo desde a última manobra [s], escalar ou array.
        a_dot: taxa de variação do semieixo [m/s] (negativa com arrasto; usa-se |ȧ|).
        a: semieixo [m].

    Returns:
        (λ [rad], λ·R⊕ [km]): o desvio no equador, com o shape de ``t``.
    """
    t = np.asarray(t, dtype=float)
    lam = 3.0 * math.pi * abs(a_dot) * t**2 / (2.0 * a * SECONDS_PER_DAY)
    return _scalar(lam), _scalar(lam * R_EARTH / 1e3)


def time_to_deadband(deadband_km: float, a_dot: float, a: float) -> float:
    """Tempo [s] até o traço sair da banda morta ±L no equador (Sem4 eq. 41).

    t = √(2·a·T_dia·L_rad/(3π|ȧ|)), com L_rad = L/R⊕. Sem4 §7.6: ȧ = −35 m/dia, ±5 km → 5,7 dias.
    """
    if a_dot == 0.0:
        return math.inf
    lam = deadband_km * 1e3 / R_EARTH
    return math.sqrt(2.0 * a * SECONDS_PER_DAY * lam / (3.0 * math.pi * abs(a_dot)))


def delta_v_for_semi_major_axis(delta_a: float, a: float) -> float:
    """Δv [m/s] para mudar o semieixo de δa numa órbita circular: Δv = v·|δa|/(2a) (Gauss)."""
    return circular_velocity(a) * abs(delta_a) / (2.0 * a)


def track_delta_v_per_year(a_dot: float, a: float) -> float:
    """Δv anual de manutenção do traço [m/s]: v·|ȧ|·ano/(2a) (repõe todo o a perdido no ano)."""
    return delta_v_for_semi_major_axis(abs(a_dot) * YEAR_S, a)


# --- LTAN: inclinação -----------------------------------------------------------------------------


def ltan_drift_from_inclination(delta_i: ArrayLike, i: float) -> FloatOrArray:
    """Deriva do LTAN [min/ano] causada por um erro de inclinação δi [rad] (Sem4 eq. 32).

    Ω̇ ∝ cos i ⇒ δΩ̇ = −Ω☉·tan i·δi numa SSO; convertido em hora local (4 min por grau) e
    acumulado num ano trópico. Na SSO de projeto (tan i = −7,66): δi = +0,01° → +1,92 min/ano.
    """
    d_raan_dot = -OMEGA_SUN * math.tan(i) * np.asarray(delta_i, dtype=float)
    return _scalar(d_raan_dot * YEAR_S * MINUTES_PER_RAD_OF_RAAN)


def ltan_drift_from_semi_major_axis(delta_a: ArrayLike, a: float) -> FloatOrArray:
    """Deriva do LTAN [min/ano] causada por um erro de semieixo δa [m] (Sem4 eq. 32).

    Ω̇ ∝ a^(−7/2) ⇒ δΩ̇ = −3,5·Ω☉·δa/a. Na órbita de projeto: δh = +1 km → −0,73 min/ano.
    """
    d_raan_dot = -3.5 * OMEGA_SUN * np.asarray(delta_a, dtype=float) / a
    return _scalar(d_raan_dot * YEAR_S * MINUTES_PER_RAD_OF_RAAN)


def ltan_offset(t: ArrayLike, di_dt: float, i: float) -> FloatOrArray:
    """Afastamento do LTAN [min] após t [s] com inclinação derivando a di/dt [rad/s] constante.

    δi(t) = (di/dt)·t ⇒ a deriva do LTAN cresce linearmente, k·δi(t), e o LTAN se afasta
    quadraticamente: ΔLTAN(t) = ½·k·(di/dt)·t², com k = ``ltan_drift_from_inclination`` por
    radiano (min/ano por rad, aqui em min/s). Na SSO de projeto com di/dt = −9,8e-5 °/dia:
    −6,87 min/ano² de aceleração, ΔLTAN = −3,43·t² min (t em anos).
    """
    k_min_per_s = float(ltan_drift_from_inclination(1.0, i)) / YEAR_S  # (min/s) por rad
    t = np.asarray(t, dtype=float)
    return _scalar(0.5 * k_min_per_s * di_dt * t**2)


def time_to_ltan_band(band_min: float, di_dt: float, i: float) -> float:
    """Tempo [s] até o LTAN sair da banda ±band_min (inverso de ``ltan_offset``).

    Na SSO de projeto, ±5 min com o di/dt do Sol: 1,21 ano.
    """
    accel = abs(float(ltan_offset(1.0, di_dt, i)))  # min por s², ½·k·|di/dt|
    return math.inf if accel == 0.0 else math.sqrt(band_min / accel)


def inclination_delta_v_per_year(di_dt: float, a: float) -> float:
    """Δv anual para anular a deriva de inclinação [m/s]: v·|di/dt|·ano (manobra fora do plano)."""
    return circular_velocity(a) * abs(di_dt) * YEAR_S


# --- Medida no Cowell e indicador por satélite ----------------------------------------------------


def measure_decay_rate(satellite: Satellite, days: float, step_s: float = 10.0) -> float:
    """ȧ médio [m/s] medido no propagador do satélite por médias em janelas nodais.

    O semieixo osculador oscila ≈ 9,5 km a cada meia volta (J2); a média numa janela de exatamente
    um período nodal (``window_mean_elements``) o elimina. As janelas começam nos nós ascendentes
    medidos (``node_crossings``) e duram o período nodal medido (espaçamento médio dos nós): uma
    janela mais curta ou mais longa que a volta deixaria um resto do termo de 2u, cuja fase
    escorrega de janela para janela e vira uma tendência falsa em ⟨a⟩ (com o período calculado
    do semieixo *osculador*, ≈ 12 s curto, o ȧ medido em 6 dias sairia −30 m/dia em vez de
    −24). ȧ é a inclinação da reta ajustada às médias contra o centro de cada janela (S4-05).

    Args:
        satellite: o satélite (``Satellite`` ou propagador); a medida começa na sua época.
        days: duração da medida [dias]; precisa caber ao menos 2 janelas.
        step_s: passo da amostragem dentro de cada janela [s].

    Raises:
        ValueError: se couber menos de 2 janelas.
    """
    t_nodes, _ = node_crossings(satellite, satellite.epoch, days * SECONDS_PER_DAY)
    if t_nodes.size < 3:
        raise ValueError("a medida de ȧ precisa de ao menos 2 períodos nodais")
    period = float(np.mean(np.diff(t_nodes)))
    starts = t_nodes[:-1]
    propagator = getattr(satellite, "propagator", satellite)
    a_mean = window_mean_elements(propagator, starts, period, step_s).a
    return float(np.polyfit(starts + period / 2.0, a_mean, 1)[0])


@dataclass(frozen=True)
class MaintenanceStatus:
    """Indicador de manutenção de um satélite (S5-05), pronto para a cena (S5-06).

    Attributes:
        sat_id: id do satélite.
        a_dot_m_day: decaimento do semieixo [m/dia] (negativo).
        track_deadband_km: banda morta do traço no equador [±km].
        days_to_track_correction: dias até o traço sair da banda (próxima manobra de a).
        track_delta_a_m: semieixo reposto por manobra [m].
        track_delta_v_per_maneuver_m_s: Δv de cada manobra de traço [m/s].
        track_maneuvers_per_year: manobras de traço por ano.
        track_delta_v_per_year_m_s: Δv anual de traço [m/s].
        di_dt_deg_day: deriva de inclinação adotada [°/dia].
        ltan_offset_1y_min: afastamento do LTAN previsto em 1 ano [min].
        ltan_band_min: banda do LTAN [±min].
        days_to_ltan_correction: dias até o LTAN sair da banda.
        inclination_delta_v_per_year_m_s: Δv anual de inclinação [m/s].
    """

    sat_id: str
    a_dot_m_day: float
    track_deadband_km: float
    days_to_track_correction: float
    track_delta_a_m: float
    track_delta_v_per_maneuver_m_s: float
    track_maneuvers_per_year: float
    track_delta_v_per_year_m_s: float
    di_dt_deg_day: float
    ltan_offset_1y_min: float
    ltan_band_min: float
    days_to_ltan_correction: float
    inclination_delta_v_per_year_m_s: float

    def as_dict(self) -> dict:
        """Dicionário serializável (JSON)."""
        return asdict(self)


def maintenance_status(
    sat_id: str,
    a_dot: float,
    a: float,
    i: float,
    di_dt: float = math.radians(SUN_INCLINATION_DRIFT_DEG_DAY) / SECONDS_PER_DAY,
    deadband_km: float = 5.0,
    ltan_band_min: float = 5.0,
) -> MaintenanceStatus:
    """Indicador de manutenção a partir de ȧ (medido) e di/dt (do Sol, S4-02).

    Pressupõe o satélite na órbita nominal agora (logo após uma manobra): a próxima correção do
    traço vem em ``time_to_deadband`` e a do LTAN em ``time_to_ltan_band``.

    Args:
        sat_id: id do satélite.
        a_dot: ȧ [m/s] (``measure_decay_rate``).
        a: semieixo médio [m].
        i: inclinação média [rad].
        di_dt: deriva de inclinação [rad/s] (padrão: −9,8e-5 °/dia, S4-02).
        deadband_km: banda morta do traço [±km].
        ltan_band_min: banda do LTAN [±min].
    """
    t_track = time_to_deadband(deadband_km, a_dot, a)
    delta_a = abs(a_dot) * t_track
    t_ltan = time_to_ltan_band(ltan_band_min, di_dt, i)
    return MaintenanceStatus(
        sat_id=sat_id,
        a_dot_m_day=a_dot * SECONDS_PER_DAY,
        track_deadband_km=deadband_km,
        days_to_track_correction=t_track / SECONDS_PER_DAY,
        track_delta_a_m=delta_a,
        track_delta_v_per_maneuver_m_s=delta_v_for_semi_major_axis(delta_a, a),
        track_maneuvers_per_year=YEAR_S / t_track,
        track_delta_v_per_year_m_s=track_delta_v_per_year(a_dot, a),
        di_dt_deg_day=math.degrees(di_dt) * SECONDS_PER_DAY,
        ltan_offset_1y_min=float(ltan_offset(YEAR_S, di_dt, i)),
        ltan_band_min=ltan_band_min,
        days_to_ltan_correction=t_ltan / SECONDS_PER_DAY,
        inclination_delta_v_per_year_m_s=inclination_delta_v_per_year(di_dt, a),
    )
