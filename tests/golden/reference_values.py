"""Valores de referência dos PDFs teóricos, para reuso nos testes de todas as semanas.

Fontes (docs/referencias/):
- Sem4: "Como se comporta um satélite em órbita heliossíncrona?" (orbita_heliossincrona.pdf)
- Sem3: "Por que satélites em LEO caem?" (decaimento_orbital.pdf)
Transcritos do Apêndice A do prompt da Semana 1 e conferidos contra os PDFs.

Regras:
- valores exatamente como aparecem nos PDFs (mesmo número de casas), em SI ou graus;
- a unidade vai no sufixo do nome (``_km``, ``_deg``, ``_s``, ``_deg_day``...);
- cada grupo cita a fonte; um número que aparece em mais de um grupo deve ser coerente
  (verificado em tests/test_golden_consistency.py).

Os valores de Sem4 vêm de teoria secular de primeira ordem em J2: valores de projeto preliminar.
"""

from dataclasses import dataclass

# ===============================================================================================
# A.1 Constantes — Sem4, tabela de símbolos; §2.1 (normalização), §2.2 (J_k, Clairaut);
#     eq. 26 (taxa heliossíncrona); §7.2 (dia nodal); código §9.1 (J2 com mais casas)
# ===============================================================================================


@dataclass(frozen=True)
class Constants:
    mu_m3_s2: float = 3.986004418e14  # Sem4, tabela de símbolos
    r_earth_km: float = 6378.137  # Sem4, tabela de símbolos
    omega_earth_rad_s: float = 7.292115e-5  # Sem4, tabela de símbolos
    j2: float = 1.08262668e-3  # Sem4 código §9.1 (texto §2.2: 1,082627e-3)
    j3: float = -2.533e-6  # Sem4 §2.2
    j4: float = -1.620e-6  # Sem4 §2.2
    t_trop_days: float = 365.2422  # Sem4, tabela de símbolos
    omega_sun_deg_day: float = 0.985647  # Sem4 eq. 26
    omega_sun_rad_s: float = 1.991064e-7  # Sem4 eq. 26
    cbar20: float = -4.8417e-4  # Sem4 §2.1 (EGM2008, totalmente normalizado)
    n20: float = 2.2361  # Sem4 §2.1 (= √5)
    c20: float = -1.0826e-3  # Sem4 §2.1 (não normalizado; J2 = −C20)
    j2_hydrostatic: float = 1.0814e-3  # Sem4 §2.2 (Clairaut)
    clairaut_m: float = 3.4614e-3  # Sem4 §2.2 (razão centrífuga/gravidade no equador)
    t_day_nod_sso_s: float = 86400.01  # Sem4 §7.2 (dia nodal em SSO)


CONSTANTS = Constants()

# ===============================================================================================
# A.2 Órbita de projeto 91/6 (h nominal 510 km) — Sem4 §7.5 (tabela da órbita escolhida),
#     §1 (excursão de altitude), §5.2–5.4 (β, eclipse), §6.3 (órbita congelada)
# ===============================================================================================


@dataclass(frozen=True)
class DesignOrbit:
    n_revs: int = 91  # Sem4 §7.5: N revoluções...
    d_days: int = 6  # ...em D dias
    gcd: int = 1
    h_nominal_label_km: float = 510.0  # rótulo nominal
    h_km: float = 509.952  # a − R⊕ (Sem4 §7.5 mostra 509,95)
    a_km: float = 6888.089
    e: float = 1.074e-3  # congelada (Sem4 §6.3)
    i_deg: float = 97.4396
    argp_deg: float = 90.0  # congelada (Sem4 §6.3)
    tkep_s: float = 5689.30
    tnod_s: float = 5696.70
    tnod_min: float = 94.945
    tnod_minus_tkep_s: float = 7.4
    tnod_minus_tkep_pct: float = 0.13
    v_circ_m_s: float = 7607.0
    raan_dot_deg_day: float = 0.985647
    argp_dot_j2_deg_day: float = -3.487  # só J2; líquido 0 com J3 (congelada)
    argp_dot_net_deg_day: float = 0.0
    dlambda_per_rev_deg: float = -23.736
    grid_deg: float = 3.956
    grid_km: float = 440.4  # no equador
    rev_per_day: float = 15.167
    h_inst_min_km: float = 502.5  # Sem4 §1: altitude instantânea
    h_inst_max_km: float = 517.4
    h_excursion_km: float = 14.8
    lat_max_deg: float = 82.56  # 180° − i
    beta_star_deg: float = 67.8  # Sem4 §5.3: |β| > β* ⇒ sem eclipse (510 km)
    eclipse_min_per_orbit: float = 35.0  # Sem4 §5.3 (LTDN 10h30)
    eclipse_fraction: float = 0.37
    beta_min_deg: float = -24.0  # Sem4 §5.4: β nessa órbita
    beta_max_deg: float = -17.0


DESIGN_ORBIT = DesignOrbit()

# ===============================================================================================
# A.3 Tabela 1 (Sem4 §4.3): SSO circulares
# ===============================================================================================


@dataclass(frozen=True)
class SsoCircularRow:
    h_km: float
    i_deg: float
    tnod_s: float
    rev_per_day: float
    mission: str = ""


SSO_CIRCULAR_TABLE1: tuple[SsoCircularRow, ...] = (
    SsoCircularRow(300, 96.672, 5438.8, 15.886),
    SsoCircularRow(400, 97.030, 5561.1, 15.536),
    SsoCircularRow(500, 97.402, 5684.4, 15.200, "SkySat, PlanetScope"),
    SsoCircularRow(510, 97.440, 5696.7, 15.167, "órbita de projeto"),
    SsoCircularRow(600, 97.788, 5808.5, 14.875),
    SsoCircularRow(705, 98.208, 5939.8, 14.546, "Landsat 8/9"),
    SsoCircularRow(786, 98.544, 6041.8, 14.300, "Sentinel-2"),
    SsoCircularRow(800, 98.603, 6059.5, 14.259),
    SsoCircularRow(1000, 99.479, 6313.9, 13.684),
    SsoCircularRow(1200, 100.420, 6571.9, 13.147),
)

# ===============================================================================================
# A.4 Tabela 3 (Sem4 §7.4): SSO de repetição. Linhas com gcd > 1 são degeneradas (sem grade).
# ===============================================================================================


@dataclass(frozen=True)
class RepeatOrbitRow:
    n_revs: int
    d_days: int
    gcd: int
    h_km: float
    i_deg: float
    grid_deg: float | None  # None nas degeneradas
    grid_km: float | None  # no equador; None nas degeneradas
    note: str = ""


SSO_REPEAT_TABLE3: tuple[RepeatOrbitRow, ...] = (
    RepeatOrbitRow(95, 6, 1, 314.80, 96.724, 3.789, 421.8),
    RepeatOrbitRow(93, 6, 3, 410.63, 97.069, None, None, "degenerada (ciclo real 2 dias)"),
    RepeatOrbitRow(91, 6, 1, 509.95, 97.440, 3.956, 440.4, "órbita de projeto"),
    RepeatOrbitRow(106, 7, 1, 517.19, 97.467, 3.396, 378.1),
    RepeatOrbitRow(121, 8, 1, 522.62, 97.488, 2.975, 331.2),
    RepeatOrbitRow(90, 6, 6, 560.99, 97.636, None, None, "degenerada (ciclo real 1 dia)"),
    RepeatOrbitRow(89, 6, 1, 612.99, 97.839, 4.045, 450.3),
    RepeatOrbitRow(175, 12, 1, 692.83, 98.159, 2.057, 229.0),
    RepeatOrbitRow(233, 16, 1, 699.59, 98.186, 1.545, 172.0, "Landsat"),
    RepeatOrbitRow(143, 10, 1, 786.12, 98.545, 2.517, 280.2, "Sentinel-2"),
)

# ===============================================================================================
# A.5 Taxas seculares de J2 (Sem4 §3.5) e números avulsos de Sem4 §3.4, §4.3, §6.1, §7.6
# ===============================================================================================


@dataclass(frozen=True)
class SecularRateRow:
    orbit: str
    i_deg: float
    raan_dot_deg_day: float
    argp_dot_deg_day: float


J2_SECULAR_RATES: tuple[SecularRateRow, ...] = (
    SecularRateRow("Molniya (inclinação crítica)", 63.435, -0.148, 0.0),
    SecularRateRow("Polar exata", 90.0, 0.0, -3.806),
    SecularRateRow("SSO de projeto", 97.44, 0.986, -3.487),  # ω-ponto só J2
    SecularRateRow("ISS (h = 420 km)", 51.64, -4.947, 3.690),
)


@dataclass(frozen=True)
class MiscSem4:
    # Sem4 §3.4: precessão do plano da ISS
    iss_raan_dot_deg_day_approx: float = -5.0
    iss_raan_full_cycle_days: float = 72.0
    # Sem4 §4.3: teto de altitude da SSO (i = 180°)
    sso_a_max_km: float = 12352.0
    sso_h_max_km: float = 5974.0
    # Sem4 §6.1, eq. 32: sensibilidade da deriva de LTAN
    tan_i_design: float = -7.66
    ltan_drift_min_per_year_per_dh_1km: float = -0.73
    ltan_drift_min_per_year_per_di_001deg: float = 1.9
    # Sem4 §7.6, eq. 40: deriva da grade para δa = −100 m (D = 6)
    grid_drift_da_m: float = -100.0
    grid_drift_deg_per_cycle: float = 0.047
    grid_drift_km_equator: float = 5.2
    # Sem4 §7.6, eq. 41: cadência de manobras com ȧ = −35 m/dia; banda [±km] → dias
    maneuver_adot_m_day: float = -35.0
    maneuver_cadence_days: tuple[tuple[float, float], ...] = ((1.0, 2.6), (5.0, 5.7), (20.0, 11.4))


MISC_SEM4 = MiscSem4()

# ===============================================================================================
# Sem3 (arrasto) — Apêndice A.5 do prompt; Sem3 §2.4 (a_D a 400 km), §3.4 (184 m/dia) e
#     §4.4 (tabela de valores de referência: ρ, decaimento e tempo de vida)
# ===============================================================================================


@dataclass(frozen=True)
class DragCase:
    h_km: float
    rho_kg_m3: float
    cd: float = 2.2
    area_over_mass_m2_kg: float = 0.005
    drag_accel_m_s2: float | None = None
    decay_m_day: float | None = None  # |ȧ|
    lifetime_years_approx: float | None = None


DRAG_SEM3: tuple[DragCase, ...] = (
    DragCase(400.0, 3.7e-12, drag_accel_m_s2=1.2e-6, decay_m_day=184.0, lifetime_years_approx=1.0),
    DragCase(500.0, 7.0e-13, decay_m_day=35.0, lifetime_years_approx=5.0),
)

# ===============================================================================================
# A.6 Constelação de 6 satélites — Sem4 §8.5, Tabela 4; verificação da grade Sem4 §9.2
# Todos compartilham a, e, i, Ω e ω da órbita de projeto; diferem apenas em ψ0.
# ===============================================================================================


@dataclass(frozen=True)
class ConstellationSat:
    name: str
    psi0_deg: float
    node_delay_s: float


CONSTELLATION_TABLE4: tuple[ConstellationSat, ...] = (
    ConstellationSat("S1", 0.0, 0.0),
    ConstellationSat("S2", 300.0, 949.5),
    ConstellationSat("S3", 240.0, 1898.9),
    ConstellationSat("S4", 180.0, 2848.4),
    ConstellationSat("S5", 120.0, 3797.8),
    ConstellationSat("S6", 60.0, 4747.3),
)


@dataclass(frozen=True)
class ConstellationOrbit:
    a_km: float = 6888.089  # Sem4 Tabela 4
    e: float = 0.001074
    i_deg: float = 97.4396
    walker: str = "97,44°: 6/1/0"


CONSTELLATION_ORBIT = ConstellationOrbit()


@dataclass(frozen=True)
class GridCheckSem4:
    # Sem4 §9.2: P = 1 em 6 dias e P = 6 em 1 dia → 91 nós, passo uniforme
    nodes_p6_1day: int = 91
    step_p6_deg: float = 3.9560
    # P = 7 em 1 dia → 107 nós, passo não uniforme
    nodes_p7_1day: int = 107
    step_p7_min_deg: float = 0.5651
    step_p7_max_deg: float = 3.3909


GRID_CHECK = GridCheckSem4()
