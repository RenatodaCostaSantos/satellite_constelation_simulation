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
    # Sem4 §1: altitude instantânea, como publicada. Diferença de arredondamento conhecida: pelos
    # elementos, a(1 ∓ e) − R⊕ = 502,554 e 517,350 km; o PDF traz 502,5 e 517,4 (diferença
    # 14,9 km vs. excursão 2ae = 14,8 km). Mantido como no PDF; ver o xfail
    # test_design_orbit_instantaneous_altitude_sem4_s1 em tests/test_golden_consistency.py.
    h_inst_min_km: float = 502.5
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

# ===============================================================================================
# Semana 2 — Apêndice A.1 do prompt da Semana 2: taxas seculares de J2 (Sem4 eqs. 23–25, 34) com
#     mais casas que a A.5 acima, recalculadas a partir dos elementos de cada caso.
#     ISS: a = R⊕ + 420 km. A inclinação crítica arccos(1/√5) vale para qualquer a e não entra
#     aqui (é fórmula, não número tabelado).
# ===============================================================================================


@dataclass(frozen=True)
class J2SecularCase:
    orbit: str
    a_km: float
    e: float
    i_deg: float
    raan_dot_deg_day: float
    argp_dot_deg_day: float


J2_SECULAR_CASES_S2: tuple[J2SecularCase, ...] = (
    J2SecularCase("ISS", 6378.137 + 420.0, 0.0, 51.64, -4.947, 3.690),
    J2SecularCase("SSO de projeto", 6888.089, 1.074e-3, 97.4396, 0.9856, -3.487),
    J2SecularCase("Polar exata", 6888.089, 0.0, 90.0, 0.0, -3.806),
    J2SecularCase("Molniya", 26561.75, 0.74, 63.4349, -0.148, 0.0),
)


@dataclass(frozen=True)
class J2SecularMiscS2:
    # SSO de projeto (A.1): ψ̇/n − 1 com ω̇ de J2 isolado. Sem4 §3.5 traz "+0,0065%" por erro de
    # digitação; o valor correto é −0,066% (≈ 3,8 s de Tnod − Tkep; ω̇ dá os outros ≈ 3,6 s).
    sso_psi_dot_over_n_minus_1: float = -6.6e-4
    # ISS (A.1): período de precessão do plano 360°/|Ω̇| e período nodal
    iss_plane_precession_days: float = 72.8
    iss_tnod_s: float = 5573.9
    iss_rev_per_day: float = 15.50


J2_SECULAR_MISC_S2 = J2SecularMiscS2()

# ===============================================================================================
# Semana 2 — Apêndice A.4 do prompt da Semana 2: referências de comportamento do Cowell, de um
#     experimento independente (Cowell DOP853, rtol = 1e-11, atol = 1e-3, J2 apenas, 10 dias a
#     60 s, raan0 = 30°, psi0 = 0). Diferenças de até ~30% são aceitáveis; os critérios de aceite
#     dos cards são os que valem.
#     argp0: SSO 90° e órbita de teste 30° (do prompt); ISS-like 0° *inferido* na S2-02, pois só
#     ele reproduz a tabela. Motivo: em órbita quase circular o semieixo osculador oscila como
#     a_osc − a_médio ≈ 1,5·J2·R⊕²/a·sin²i·cos 2u (u = ω + ν). Tomar os elementos nominais como
#     osculadores em u0 = 0° (nó) põe a_médio ≈ 6,0 km *abaixo* do nominal; em u0 = 90° (argp0 =
#     90°, latitude máxima), 6,0 km *acima*. Com argp0 = 90° a ISS-like dá a_médio − a_nominal =
#     +6,0 km, Ω̇ ajustado −4,933 °/dia (teoria@médio −4,929, +0,08%; @nominal −0,27%) e diferença
#     MeanJ2(médio) × Cowell de 3,4 km em 24 h: os critérios valem, mas os números não são os da
#     tabela. A escolha de argp0 para a ISS-like na S2-06 deve ser coerente com esta linha.
# ===============================================================================================


@dataclass(frozen=True)
class CowellReferenceCase:
    orbit: str
    a_km: float
    e: float
    i_deg: float
    argp0_deg: float
    a_mean_minus_nominal_km: float
    raan_dot_fit_deg_day: float
    raan_dot_theory_at_mean_deg_day: float
    diff_mean_j2_mean_elements_24h_km: float
    a_osc_peak_to_peak_km: float | None = None  # semieixo osculador numa revolução
    nominal_error_1_6_12_24h_km: tuple[float, ...] | None = None  # MeanJ2 nominal vs Cowell


COWELL_RAAN0_DEG: float = 30.0
COWELL_PSI0_DEG: float = 0.0

COWELL_REFERENCE_S2: tuple[CowellReferenceCase, ...] = (
    CowellReferenceCase(
        orbit="SSO",
        a_km=6888.089,
        e=1.074e-3,
        i_deg=97.4396,
        argp0_deg=90.0,
        a_mean_minus_nominal_km=9.5,
        raan_dot_fit_deg_day=0.9794,
        raan_dot_theory_at_mean_deg_day=0.9803,
        diff_mean_j2_mean_elements_24h_km=6.8,
        a_osc_peak_to_peak_km=19.0,
        nominal_error_1_6_12_24h_km=(79.0, 363.0, 688.0, 1333.0),
    ),
    CowellReferenceCase(
        orbit="ISS-like",
        a_km=6378.137 + 420.0,
        e=5e-4,
        i_deg=51.64,
        argp0_deg=0.0,
        a_mean_minus_nominal_km=-6.0,
        raan_dot_fit_deg_day=-4.968,
        raan_dot_theory_at_mean_deg_day=-4.9641,
        diff_mean_j2_mean_elements_24h_km=14.0,
        a_osc_peak_to_peak_km=12.0,
        nominal_error_1_6_12_24h_km=(47.0, 226.0, 449.0, 874.0),
    ),
    CowellReferenceCase(
        orbit="e = 0,02",
        a_km=7200.0,
        e=0.02,
        i_deg=45.0,
        argp0_deg=30.0,
        a_mean_minus_nominal_km=-2.5,
        raan_dot_fit_deg_day=-4.6251,
        raan_dot_theory_at_mean_deg_day=-4.6197,
        diff_mean_j2_mean_elements_24h_km=6.6,
    ),
)

# ω-ponto na órbita e = 0,02 (A.4): ajustado vs teoria no semieixo médio
COWELL_ARGP_DOT_FIT_E002_DEG_DAY: float = 4.9095
COWELL_ARGP_DOT_THEORY_E002_DEG_DAY: float = 4.9018

# ===============================================================================================
# Semana 2 — Apêndice A.3 do prompt da Semana 2: WGS84 e geometria. f com o valor completo do
#     WGS84 (NIMA TR8350.2); o texto do Sem4 usa 1/298,257. São Bento do Sapucaí a confirmar na
#     S3-04.
# ===============================================================================================


@dataclass(frozen=True)
class Wgs84Reference:
    a_m: float = 6378137.0
    inv_f: float = 298.257223563
    b_m: float = 6356752.314245  # a(1 − f)
    e2: float = 6.69437999e-3  # f(2 − f)
    equator_minus_pole_km: float = 21.385  # a − b
    geocentric_minus_geodetic_45_deg: float = -0.1924  # latitude a 45° geodésicos, h = 0
    sbs_lat_deg: float = -22.69  # São Bento do Sapucaí (≈)
    sbs_lon_deg: float = -45.73
    sbs_alt_m: float = 900.0


WGS84_REF = Wgs84Reference()

# ===============================================================================================
# Semana 2 — card S2-04 do prompt da Semana 2: traço no solo (MeanJ2Propagator, 24 h a 10 s).
#     Latitude máxima geocêntrica = 180° − i (SSO, DESIGN_ORBIT.lat_max_deg) ou i (ISS-like);
#     a geodésica é maior (≈ +0,05° na SSO). Deslocamento do nó por revolução e Tnod da SSO estão
#     em DESIGN_ORBIT (Sem4 §7.1, §7.5).
# ===============================================================================================


@dataclass(frozen=True)
class GroundTrackS2:
    sso_lat_max_geodetic_deg: float = 82.61  # ≈ (abs 0,1°)
    iss_like_i_deg: float = 51.64  # latitude máxima geocêntrica da ISS-like
    closure_91_6_ref_deg: float = 1.3e-4  # |λ(nó k+91) − λ(nó k)| após 6 dias (≈)
    # Cruzamentos do antimeridiano na SSO. O card S2-04 dizia {15, 16} (≈ 15,17 rev/dia), mas a
    # SSO é retrógrada: por revolução nodal a longitude do subponto recua 360° + 23,736°, logo
    # (360 + 23,736)/360 × 86400/5696,70 = 16,17 cruzamentos/dia → 16 ou 17 em 24 h, conforme Ω0;
    # no ciclo 91/6, (360 + 23,736)/360 × 91 = 97,0. Corrigido na revisão da S2-04 (dedução acima,
    # confirmada por simulação com vários Ω0).
    sso_antimeridian_crossings_24h: tuple[int, ...] = (16, 17)
    sso_antimeridian_crossings_per_day: float = 16.17
    sso_antimeridian_crossings_91_6_cycle: int = 97


GROUND_TRACK_S2 = GroundTrackS2()

# ===============================================================================================
# Semana 2 — Apêndice A.2 e card S2-05 do prompt da Semana 2: projeto de SSO (Sem4 §4.3, §6.1).
#     Tabela 1, teto, tan i e derivas de LTAN já estão em SSO_CIRCULAR_TABLE1 e MISC_SEM4.
#     A.5: convenção de ω̇ no ciclo 91/6 (pendência para a Semana 5).
# ===============================================================================================


@dataclass(frozen=True)
class SsoDesignS2:
    di_dh_510_deg_per_km: float = 0.0038  # di/dh da SSO circular em 510 km
    cos_i_rel_change_e001: float = -2.0e-4  # cos i ∝ (1 − e²)²: variação relativa para e = 0,01
    # A.5: solução de repetição 91/6 com ω̇ de J2 isolado (código do PDF) vs órbita congelada
    # (J3 cancela ω̇ → ω̇ líquido = 0); Tnod congelada para o mesmo a ≈ 5693,1 s
    repeat_91_6_j2_only_h_km: float = 509.952
    repeat_91_6_j2_only_i_deg: float = 97.4396
    repeat_91_6_frozen_h_km: float = 512.89
    repeat_91_6_frozen_i_deg: float = 97.4507
    tnod_frozen_same_a_s: float = 5693.1


SSO_DESIGN_S2 = SsoDesignS2()

# ===============================================================================================
# Semana 2 — card S2-07 do prompt da Semana 2: cena da GUI com subponto geodésico WGS84.
#     O card dizia "≈ 51,7° (abs 0,1°)" para a latitude geodésica máxima da ISS-like. Corrigido:
#     a diferença geodésica − geocêntrica na latitude φ e altitude h é ≈ f·sin 2φ·R⊕/(R⊕ + h)
#     [rad] → 0,1924°·sin(103,3°)·6378/6798 ≈ 0,176° a 51,64° e 420 km, logo 51,64 + 0,176 =
#     51,82° (a mesma conta dá os 82,61° da SSO). Confirmado pela conversão validada contra o
#     Skyfield (S2-03) e pelo traço no solo (S2-04: 51,815°).
# ===============================================================================================


@dataclass(frozen=True)
class SceneS2:
    iss_like_lat_max_geocentric_deg: float = 51.64  # = i
    iss_like_lat_max_geodetic_deg: float = 51.82
    sso_lat_max_geodetic_deg: float = 82.61
    max_abs_dlat_vs_spherical_deg: float = 0.2  # |lat geodésica − geocêntrica| ≤ 0,2°
    alt_vs_spherical_km: tuple[float, float] = (-25.0, 25.0)  # alt geodésica − (r − R⊕)
    iss_raan_regression_48h_deg: float = -9.9  # Ω̇·2 dias, ISS-like (−4,947 °/dia)


SCENE_S2 = SceneS2()

# ===============================================================================================
# Semana 3 — Apêndice A.1 do prompt da Semana 3: TLE de teste da ISS (congelado em
#     tests/data/iss_sample.tle). Época 2019-12-09 16:38:29 UTC (mínimo solar); checksums válidos.
#     B* = 0,38792e-4; movimento médio 15,50103472 rev/dia. Apêndice A.6: limiares de idade do TLE
#     (aviso > 3 dias, aviso forte > 14 dias) e a CelesTrak pede ≤ 1 consulta a cada 2 h.
# ===============================================================================================


@dataclass(frozen=True)
class IssTleA1:
    name: str = "ISS (ZARYA)"
    line1: str = "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991"
    line2: str = "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482"
    norad_id: int = 25544
    epoch_utc_iso: str = "2019-12-09T16:38:29+00:00"  # tolerância do card: abs = 1 s
    inclination_deg: float = 51.6439
    mean_motion_rev_day: float = 15.50103472
    bstar: float = 3.8792e-5


ISS_TLE_A1 = IssTleA1()


@dataclass(frozen=True)
class TleAgingA6:
    warning_days: float = 3.0  # aviso "precisão degradada"
    strong_warning_days: float = 14.0  # aviso forte
    celestrak_min_interval_h: float = 2.0  # padrão de max_cache_age


TLE_AGING = TleAgingA6()

# ===============================================================================================
# Semana 3 — Apêndice A.1 do prompt da Semana 3: vetores de regressão do SGP4 para o TLE de teste
#     (sgp4 2.27, WGS72, saída TEME em km e km/s; tolerâncias do card S3-02: 1 m e 1e-3 m/s).
#     Sanidade em t = 0: |r| ≈ 6787,2 km, |v| ≈ 7,665 km/s; período 2π/n ≈ 5573,9 s (coincide com
#     o Tnod da ISS na Semana 2). Skyfield 1.55 (ITRS, timescale builtin) × Rz(GMST)·TEME: diferença
#     máxima 86 m, média 70 m em 3 dias (UT1 ≈ UTC e sem movimento polar); critério < 150 m.
# ===============================================================================================


@dataclass(frozen=True)
class Sgp4RegressionRow:
    t_h: float
    r_km: tuple[float, float, float]
    v_km_s: tuple[float, float, float]


SGP4_REGRESSION_A1: tuple[Sgp4RegressionRow, ...] = (
    Sgp4RegressionRow(0, (3469.947984, -2690.388430, 5175.831925),
                      (5.810229142, 4.802261185, -1.388280333)),
    Sgp4RegressionRow(24, (-3591.826765, 2723.666458, -5090.448286),
                      (-5.927709568, -4.496384380, 1.785277102)),
    Sgp4RegressionRow(72, (-3873.542270, 2837.516435, -4813.891216),
                      (-6.055171616, -3.901758569, 2.582041156)),
    Sgp4RegressionRow(168, (-4346.138715, 3409.792942, -3962.984695),
                      (-5.869671636, -2.825677530, 4.019415520)),
)  # fmt: skip


@dataclass(frozen=True)
class Sgp4SanityA1:
    r0_km: float = 6787.2
    v0_km_s: float = 7.665
    osc_inclination_deg: float = 51.64  # inclinação osculadora em t = 0 (card: abs 0,05°)
    period_s: float = 5573.9  # 2π/n do movimento médio (card: abs 1 s)
    skyfield_ecef_max_diff_m: float = 86.0  # referência; critério do card: < 150 m em 3 dias
    skyfield_ecef_mean_diff_m: float = 70.0
    # O card S3-02 dizia "|lat| ≤ 51,7° (geodésica)". Corrigido: 51,7° limita a latitude
    # GEOCÊNTRICA (≤ i osculadora, 51,62–51,66° em 1 dia). A geodésica fica Δφ ≈ f·sin 2φ·R⊕/(R⊕+h)
    # ≈ 0,176° acima (mesma correção de SCENE_S2): 51,62 + 0,18 = 51,80° medido.
    lat_max_geocentric_deg: float = 51.7


SGP4_SANITY_A1 = Sgp4SanityA1()

# ===============================================================================================
# Semana 3 — Apêndice A.3 do prompt da Semana 3: atmosfera exponencial por partes (Vallado).
#     Cada faixa h0 ≤ h < próxima h0: ρ = ρ0·exp(−(h − h0)/H). Continuidade verificada no prompt:
#     salto relativo < 1e-4 em todas as fronteiras (critério do card: < 1e-3).
# ===============================================================================================

ATMOSPHERE_A3: tuple[tuple[float, float, float], ...] = (
    # (h0 [km], ρ0 [kg/m³], H [km])
    (100.0, 5.297e-7, 5.877), (110.0, 9.661e-8, 7.263), (120.0, 2.438e-8, 9.473),
    (130.0, 8.484e-9, 12.636), (140.0, 3.845e-9, 16.149), (150.0, 2.07e-9, 22.523),
    (180.0, 5.464e-10, 29.740), (200.0, 2.789e-10, 37.105), (250.0, 7.248e-11, 45.546),
    (300.0, 2.418e-11, 53.628), (350.0, 9.518e-12, 53.298), (400.0, 3.725e-12, 58.515),
    (450.0, 1.585e-12, 60.828), (500.0, 6.967e-13, 63.822), (600.0, 1.454e-13, 71.835),
    (700.0, 3.614e-14, 88.667), (800.0, 1.17e-14, 124.640), (900.0, 5.245e-15, 181.050),
    (1000.0, 3.019e-15, 268.000),
)  # fmt: skip


# ===============================================================================================
# Semana 3 — Sem3 Tabela 3 (C_D = 2,2, A/m = 0,005 → C_D·A/m = 0,011; atividade solar moderada),
#     valores arredondados no PDF, e Apêndice A.4 do prompt da Semana 3 (calculados com a tabela
#     A.3, altitude esférica, sem co-rotação). Tolerâncias do card S3-03: ρ e |dh/dt| ±5%.
# ===============================================================================================


@dataclass(frozen=True)
class DecayRowSem3:
    h_km: float
    rho_table3_kg_m3: float  # Sem3 Tabela 3 (arredondado)
    dhdt_table3_m_day: float  # |dh/dt| Sem3 Tabela 3 (arredondado)
    dhdt_a4_m_day: float  # |dh/dt| calculado no Apêndice A.4...
    dhdt_a4_decimals: int  # ...impresso com estas casas decimais


DECAY_TABLE3: tuple[DecayRowSem3, ...] = (
    DecayRowSem3(200.0, 2.8e-10, 1.4e4, 13573.0, 0),
    DecayRowSem3(300.0, 2.4e-11, 1.2e3, 1186.0, 0),
    DecayRowSem3(400.0, 3.7e-12, 184.0, 184.0, 1),
    DecayRowSem3(500.0, 7.0e-13, 35.0, 34.7, 1),
    DecayRowSem3(600.0, 1.5e-13, 7.0, 7.29, 2),
    DecayRowSem3(800.0, 1.2e-14, 0.6, 0.59, 2),
    DecayRowSem3(1000.0, 3.0e-15, 0.16, 0.16, 2),
)


@dataclass(frozen=True)
class DragDecayA4:
    cd: float = 2.2
    area_over_mass_m2_kg: float = 0.005
    cd_a_m_m2_kg: float = 0.011
    v_400_m_s: float = 7669.0  # Sem3 §2.4
    drag_accel_400_m_s2: float = 1.2e-6  # Sem3 §2.4 (card: rel 3%)
    g_400_m_s2: float = 8.68  # Sem3 §2.4
    drag_over_gravity_400: float = 1.4e-7  # Sem3 §2.4
    da_per_rev_400_m: float = -11.8  # Sem3 eq. 15 (A.4: −11,83; card: rel 3%)
    da_per_day_400_m: float = 184.0  # |Δa| em 1 dia (A.4: ≈ 184 m/dia; card: rel 5%)
    h_reentry_km: float = 120.0  # altitude que encerra a integração do tempo de vida (A.4)


DRAG_DECAY_A4 = DragDecayA4()


@dataclass(frozen=True)
class LifetimeCaseA4:
    h0_km: float
    cd_a_m_m2_kg: float
    lifetime_years: float  # integração da Sem3 eq. 22 com a tabela A.3, reentrada em 120 km
    pdf_years: float  # valor citado no PDF (Sem3 §5.1, Tabela 3), para comparação


LIFETIMES_A4: tuple[LifetimeCaseA4, ...] = (
    LifetimeCaseA4(500.0, 0.011, 4.64, 4.6),  # o PDF diz 4,6
    LifetimeCaseA4(500.0, 0.044, 1.16, 1.16),  # Sem3 Figura 2
    LifetimeCaseA4(500.0, 0.0022, 23.2, 23.2),  # Sem3 Figura 2
    LifetimeCaseA4(400.0, 0.011, 0.77, 1.0),  # tabela do PDF: "~1 ano"
    LifetimeCaseA4(600.0, 0.011, 23.7, 25.0),  # tabela do PDF: "~25 anos"
)


# ===============================================================================================
# Semana 3 — Apêndice A.5 do prompt da Semana 3: erro de posição 3D Cowell (J2 + arrasto) × SGP4
#     para o TLE de teste (DOP853, rtol = 1e-10, atol = 1e-3, partindo do estado do SGP4 em t = 0).
#     Critérios do card: ajustado ≤ 2 / 8 / 40 km; C_D·A/m ajustado em [5e-4, 5e-3]. Os números do
#     arrasto padrão NÃO são critério: são o resultado a explicar (mínimo solar em dez/2019).
# ===============================================================================================


@dataclass(frozen=True)
class CowellVsSgp4Row:
    label: str
    cd_a_m_m2_kg: float
    altitude_model: str
    err_24h_km: float
    err_72h_km: float
    err_7d_km: float


COWELL_VS_SGP4_A5: tuple[CowellVsSgp4Row, ...] = (
    CowellVsSgp4Row("Cowell J2 apenas", 0.0, "geodetic", 1.0, 8.7, 47.0),
    CowellVsSgp4Row("J2 + arrasto padrão (esférica)", 0.0088, "spherical", 6.6, 59.6, 325.0),
    CowellVsSgp4Row("J2 + arrasto padrão (geodésica)", 0.0088, "geodetic", 5.8, 52.4, 286.0),
    CowellVsSgp4Row("J2 + arrasto ajustado (72 h)", 0.00163, "geodetic", 0.3, 2.7, 14.5),
)


@dataclass(frozen=True)
class BallisticFitA5:
    fitted_cd_a_m_m2_kg: float = 0.00163  # A/m ≈ 7,4e-4 com C_D = 2,2
    fitted_range_m2_kg: tuple[float, float] = (5e-4, 5e-3)  # critério do card
    max_err_km: tuple[float, float, float] = (2.0, 8.0, 40.0)  # 24 h / 72 h / 7 dias (critério)


BALLISTIC_FIT_A5 = BallisticFitA5()

# ===============================================================================================
# Semana 3 — Apêndices A.2 e A.6 do prompt da Semana 3: estação de São Bento do Sapucaí
#     (coordenadas APROXIMADAS, a confirmar), convenções topocêntricas, footprint (Sem4 §9.4) e a
#     concordância do pipeline SGP4 → Rz(GMST) → geodésico/ENU com o Skyfield (3 dias, el > 0).
# ===============================================================================================


@dataclass(frozen=True)
class StationA6:
    lat_deg: float = -22.69
    lon_deg: float = -45.73
    alt_m: float = 900.0
    ecef_m: tuple[float, float, float] = (4110248.7, -4216342.7, -2445429.5)  # A.2
    # Skyfield 1.55 × satsim, TLE A.1, 3 dias a cada 30 s, só el > 0 (A.2):
    skyfield_max_daz_deg: float = 0.026  # critério do card: < 0,05°
    skyfield_max_del_deg: float = 0.009  # critério: < 0,05°
    skyfield_max_drange_m: float = 74.0  # critério: < 150 m


STATION_A6 = StationA6()


@dataclass(frozen=True)
class FootprintA6:
    # λ = arccos(R⊕/(R⊕ + h)·cos ε) − ε (Sem4 §9.4)
    # O prompt diz "h = 400 km, ε = 0° → 19,7°". Corrigido: arccos(6378,137/6778,137) = 19,782°
    # (o prompt truncou; arredondado seria 19,8°).
    h400_el0_deg: float = 19.782  # h = 400 km, ε = 0° (horizonte)
    h420_el10_deg: float = 12.5  # h = 420 km, ε = 10° (critério de sanidade, não de precisão)


FOOTPRINT_A6 = FootprintA6()

# ===============================================================================================
# Semana 3 — Apêndice A.2 do prompt da Semana 3: passagens da ISS (TLE de teste A.1) sobre São
#     Bento do Sapucaí, el_min = 10°, 3 dias após a época. Instantes em SEGUNDOS DESDE A ÉPOCA DO
#     TLE; UTC = época + t. Tolerâncias do card S3-05: AOS/TCA/LOS ±1 s; el. máx. ±0,05°;
#     azimutes ±0,1°. Skyfield find_events: TCA ≤ 0,1 s, AOS/LOS ≤ 0,19 s, el. máx. ≤ 0,01°.
#     Com el_min = 5°: 11 passagens (22 cruzamentos); as de 10° estão contidas nas de 5°.
#
#     Correção (passagem 9): o prompt traz AOS = 237528,4 s, incoerente com a própria linha —
#     a coluna UTC (10:37:28) e a duração (132 s) dão AOS = 237538,4 s; confirmado pelo satsim
#     (237538,7 s) e pelo Skyfield (237538,66 s). Erro de digitação (2 ↔ 3); ver
#     test_passes_a2_internal_consistency em tests/test_golden_consistency.py.
# ===============================================================================================


@dataclass(frozen=True)
class PassA2:
    aos_s: float
    tca_s: float
    los_s: float
    duration_s: float
    max_el_deg: float
    az_aos_deg: float
    az_tca_deg: float
    az_los_deg: float
    aos_utc: str  # coluna UTC (AOS → LOS) do A.2, "dd/mm HH:MM:SS"
    los_utc: str


PASSES_A2: tuple[PassA2, ...] = (
    PassA2(29258.7, 29457.1, 29657.6, 399, 54.72, 333.4, 51.9, 130.3,
           "10/12 00:46:08", "10/12 00:52:46"),
    PassA2(64638.7, 64833.5, 65025.6, 387, 46.88, 232.8, 307.7, 22.7,
           "10/12 10:35:48", "10/12 10:42:14"),
    PassA2(112786.2, 112954.4, 113124.0, 338, 24.27, 357.9, 54.3, 110.6,
           "10/12 23:58:15", "11/12 00:03:53"),
    PassA2(118618.4, 118755.1, 118892.7, 274, 16.69, 265.9, 223.6, 181.2,
           "11/12 01:35:27", "11/12 01:40:02"),
    PassA2(148126.5, 148328.2, 148526.8, 400, 71.54, 215.4, 130.3, 45.3,
           "11/12 09:47:15", "11/12 09:53:56"),
    PassA2(196377.3, 196454.1, 196531.0, 154, 11.70, 34.0, 56.4, 78.8,
           "11/12 23:11:26", "11/12 23:14:00"),
    PassA2(202056.8, 202241.1, 202427.1, 370, 31.53, 291.9, 226.7, 161.6,
           "12/12 00:46:06", "12/12 00:52:16"),
    PassA2(231634.9, 231818.7, 232000.0, 365, 31.40, 198.4, 133.3, 68.2,
           "12/12 08:59:04", "12/12 09:05:09"),
    PassA2(237538.4, 237604.9, 237670.9, 132, 11.26, 284.1, 303.5, 323.0,
           "12/12 10:37:28", "12/12 10:39:40"),  # AOS corrigido (prompt: 237528,4)
)  # fmt: skip


@dataclass(frozen=True)
class PassesA2Meta:
    el_min_deg: float = 10.0
    window_days: float = 3.0
    n_passes_10deg: int = 9
    n_passes_5deg: int = 11
    skyfield_max_dt_s: float = 0.19  # AOS/LOS (TCA ≤ 0,1 s); critério do card: < 1 s
    tol_time_s: float = 1.0
    tol_el_deg: float = 0.05
    tol_az_deg: float = 0.1


PASSES_A2_META = PassesA2Meta()

# ===============================================================================================
# Semana 4 — Apêndice A.2 do prompt da Semana 4: gravidade zonal J2–J4 (Sem4 eq. 6, card S4-01).
#     Experimento de referência independente (23/10/2026): forças escritas do zero e integradas no
#     CowellPropagator do projeto (DOP853). J2 1,08262668e-3, J3 −2,533e-6, J4 −1,620e-6; µ e R⊕
#     do projeto. Acelerações em m/s², em r = 7000 km; o ponto de 45° N é 7000 km·(cos 45°, 0,
#     sin 45°) (o prompt imprime (4949,7; 0; 4949,7) km).
# ===============================================================================================


@dataclass(frozen=True)
class ZonalPointA2S4:
    label: str
    r_km: tuple[float, float, float]
    a_j2: tuple[float, float, float]
    a_j3: tuple[float, float, float]
    a_j4: tuple[float, float, float]


_C45_KM = 7000.0 * 0.5**0.5

ZONAL_POINTS_A2_S4: tuple[ZonalPointA2S4, ...] = (
    ZonalPointA2S4("equador", (7000.0, 0.0, 0.0),
                   (-1.0967390000e-2, 0.0, 0.0),
                   (0.0, 0.0, -2.3380594778e-5),
                   (-1.7031041069e-5, 0.0, 0.0)),
    ZonalPointA2S4("45N", (_C45_KM, 0.0, _C45_KM),
                   (1.1632673762e-2, 0.0, -3.8775579205e-3),
                   (-9.7419144910e-6, 0.0, 2.5328977677e-5),
                   (9.0320734728e-6, 0.0, 1.7060583226e-5)),
    ZonalPointA2S4("polo N", (0.0, 0.0, 7000.0),
                   (0.0, 0.0, 2.1934780000e-2),
                   (0.0, 0.0, -6.2348252742e-5),
                   (0.0, 0.0, -4.5416109517e-5)),
)  # fmt: skip


@dataclass(frozen=True)
class ZonalA2S4:
    # Medidas do experimento (A.2): parte J2 × j2_acceleration (rel máx, 200 pontos de LEO),
    # gradiente numérico do potencial (h = 1 m) e energia no Cowell J2–J4 (10 dias, rtol 1e-11)
    j2_consistency_rel: float = 1.3e-15  # critério do card: ≤ 1e-12
    gradient_rel: float = 1.4e-9  # critério: ≤ 1e-7
    energy_rel_10d: float = 3.4e-11  # critério: ≤ 1e-9
    # Elementos de projeto tomados como osculadores (Ω0 = 30°, ψ0 = 0, 2026-10-23 00:00 UTC),
    # 30 dias a 60 s, fit_rate de Ω [°/dia]
    raan_dot_j2_deg_day: float = 0.979419
    raan_dot_j2_j3_deg_day: float = 0.979418
    raan_dot_j2_j4_deg_day: float = 0.977169
    raan_dot_ratio_j2j4_j2: float = 0.997703  # tolerância do card: abs 2e-5
    epoch_utc: str = "2026-10-23T00:00:00+00:00"
    raan0_deg: float = 30.0
    psi0_deg: float = 0.0


ZONAL_A2_S4 = ZonalA2S4()

# ===============================================================================================
# Semana 4 — Apêndice A.3 do prompt da Semana 4: Sol, Lua e terceiro corpo (card S4-02).
#     Quadro: equador verdadeiro e equinócio da data (= ECI/TEME do projeto). Posições do JPL
#     DE421 via Skyfield 1.55 (experimento de referência, 23/10/2026), em km, 00:00 UTC.
#     Erros em 2026 (a cada 6 h, 1460 instantes) × DE421: Sol (Astronomical Almanac) direção máx.
#     0,0140°, |Δr|/r ≤ 8,5e-5; Lua (Montenbruck & Gill §3.3.2, da data) direção máx. 0,085°,
#     |Δr| ≤ 445 km. Nota (S4-02): esses |Δr| são diferenças de DISTÂNCIA (|r| − |r_DE421|); a
#     diferença vetorial inclui o erro de direção (Lua em 2026-10-23: distância 439 km, vetor
#     494 km; Sol: o erro de direção sozinho já dá ~2,2e-4 relativo).
# ===============================================================================================


@dataclass(frozen=True)
class EphemerisPointA3:
    date_utc: str
    sun_km: tuple[float, float, float]
    moon_km: tuple[float, float, float]


EPHEMERIS_A3: tuple[EphemerisPointA3, ...] = (
    EphemerisPointA3("2026-03-20T00:00:00+00:00",
                     (148953369.6, -1445437.0, -626444.2), (362066.5, 61667.2, 46759.2)),
    EphemerisPointA3("2026-06-21T00:00:00+00:00",
                     (872345.9, 139472241.3, 60464721.9), (-375552.0, 72784.2, 20825.8)),
    EphemerisPointA3("2026-10-23T00:00:00+00:00",
                     (-129426274.1, -67475521.4, -29252187.5), (376788.1, -72295.8, -18052.9)),
)  # fmt: skip


@dataclass(frozen=True)
class ThirdBodyA3:
    # Tolerâncias do card S4-02 e erros máximos medidos em 2026 (A.3)
    sun_dir_tol_deg: float = 0.02
    sun_dist_rel_tol: float = 1e-4
    moon_dir_tol_deg: float = 0.1
    moon_dist_tol_km: float = 500.0
    sun_dir_max_2026_deg: float = 0.0140
    sun_dist_rel_max_2026: float = 8.5e-5
    moon_dir_max_2026_deg: float = 0.085
    moon_dist_max_2026_km: float = 445.0
    sun_mg_rejected_deg: float = 0.094  # Sol de M&G §3.3.2 (periélio fixo)
    no_precession_sun_deg: float = 0.47  # sem a correção para o equinócio da data
    no_precession_moon_deg: float = 0.46
    # 2026-10-23 00:00 UTC
    epoch_utc: str = "2026-10-23T00:00:00+00:00"
    sun_ra_aa_deg: float = 207.5235  # fórmula do AA (tolerância do card: abs 0,001°)
    sun_ra_de421_deg: float = 207.5350
    sun_distance_au: float = 0.995080
    ltan_raan_deg: float = 5.0235  # Ω que dá LTAN 22h30 nessa data
    ltan_h: float = 22.500  # tolerância: abs 0,005 h
    # Terceiro corpo em r = (R⊕ + 500 km, 0, 0), m/s² (tolerância do card: rel 1%)
    mu_sun_m3_s2: float = 1.32712440018e20  # IERS Conventions 2010
    mu_moon_m3_s2: float = 4.9028000661e12  # DE421
    au_m: float = 149597870700.0  # IAU 2012
    r_alt_km: float = 500.0
    a_sun_m_s2: float = 5.002320e-7
    a_moon_m_s2: float = 1.200458e-6
    # SSO calibrada (A.7), Cowell J2–J4, 30 dias: efeito de Sol/Lua em Ω̇. Discrepância com o
    # Sem4 §10, que estima a precessão nodal lunissolar em "~1e-3 grau/dia em LEO" (~40× maior).
    lunisolar_raan_dot_deg_day: float = -2.45e-5
    lunisolar_ltan_drift_min_year: float = -0.035
    sem4_s10_estimate_deg_day: float = 1e-3
    i_calibrated_with_sun_moon_deg: float = 97.46344  # +0,0002° sobre a calibrada só J2–J4
    # Medido no S4-02 (não substitui o A.3; a conversão médio → osculador só chega no S4-04):
    # elementos do A.7 tomados como osculadores (Ω0 = 5,0235°, LTAN 22h30), Cowell J2–J4 ± Sol/Lua,
    # 30 dias a 60 s. O ajuste linear de Ω dá −3,52e-4 °/dia (Sol sozinho −2,17e-4), porque o Sol
    # produz di/dt ≈ −9,8e-5 °/dia (deriva de inclinação da SSO) e, via J2, ΔΩ ∝ t²
    # (−0,0058° previsto, −0,0065° medido em 30 dias com o Sol). Nos primeiros 5 dias
    # ΔΩ ≈ −1e-4° (≈ −2e-5 °/dia), coerente com o A.3. Reavaliar com o estado convertido no S4-05.


THIRD_BODY_A3 = ThirdBodyA3()

# ===============================================================================================
# Semana 4 — Apêndice A.4 do prompt da Semana 4: pressão de radiação solar e eclipse (card S4-03).
#     P☉ = 4,56e-6 N/m² a 1 UA (Montenbruck & Gill §3.4). Sombra cilíndrica (Sem4 §5.3), fração de
#     eclipse de órbita circular (Sem4 eq. 30). SSO de demonstração do A.7: a 6888,089 km,
#     i 97,46325°, e 1,07721e-3, ω 90°, ψ0 = 0, Ω = 5,0235° (LTAN 22h30 em 2026-10-23 00:00 UTC).
# ===============================================================================================


@dataclass(frozen=True)
class SrpEclipseA4:
    p_sun_n_m2: float = 4.56e-6
    cr: float = 1.3
    area_over_mass_m2_kg: float = 0.01
    a_1au_m_s2: float = 5.928e-8  # P☉·C_R·A/m a 1 UA
    epoch_utc: str = "2026-10-23T00:00:00+00:00"
    sun_distance_au: float = 0.995080
    a_epoch_m_s2: float = 5.986762e-8  # tolerância do card: rel 1e-4
    # SSO de demonstração (LTAN 22h30), Cowell completo no A.4 (J2–J4 no teste; diferença
    # desprezível), 24 h a 20 s, sombra cilíndrica
    demo_a_km: float = 6888.089
    demo_i_deg: float = 97.46325
    demo_e: float = 1.07721e-3
    demo_argp_deg: float = 90.0
    demo_raan_deg: float = 5.0235
    demo_beta_min_deg: float = 23.41
    demo_beta_max_deg: float = 23.43
    demo_eclipse_fraction: float = 0.3603  # tolerância do card: abs 0,01
    demo_eclipse_fraction_eq30: float = 0.3650  # Sem4 eq. 30 com β médio
    demo_eclipse_min_per_rev: float = 34.7  # por volta de 5 696 s (Sem4 §5.3: "~35 min")
    demo_period_s: float = 5696.0
    # Caso meio-dia/meia-noite (Ω = 30°, LTAN ≈ 0h10), 24 h a 5 s
    midnight_raan_deg: float = 30.0
    midnight_beta_deg: float = -0.94
    midnight_eclipse_fraction: float = 0.3720
    midnight_eclipse_fraction_eq30: float = 0.3767


SRP_ECLIPSE_A4 = SrpEclipseA4()

# ===============================================================================================
# Semana 4 — Apêndice A.6 do prompt da Semana 4: conversão médio → osculador (card S4-04).
#     Órbita de projeto (a 6888,089 km, e 1,074e-3, i 97,4396°, ω 90°, Ω 30°, ψ 0; 2026-10-23
#     00:00 UTC), modelos J2–J4, médias por janela nodal com passo 10 s.
#
#     Discrepância (S4-04): o experimento do A.6 parece ter usado janelas do período KEPLERIANO
#     (5 689,30 s), não do nodal (5 696,70 s) que o card pede. Com a janela kepleriana o satsim
#     reproduz o A.6 (resíduos 9 495,97 / 74,48 / 0,2469 m; osculador a 6 878 667,3 m, e 7,2726e-4,
#     i 97,444713°). Com a janela nodal (implementada): resíduos 9 483,6 / 74,30 / 0,2460 /
#     2,26e-4 m; osculador a 6 878 679,4 m, e 7,2492e-4, i 97,444706°, ω 270,0°. A janela nodal é a
#     correta: os termos de curto período dependem de cos 2u, e u completa 2π num período nodal;
#     a kepleriana fica 7,4 s curta e enviesa ⟨a⟩ em ≈ A·δT/T ≈ 9,4 km × 7,4/5 697 ≈ 12 m. É
#     também o que explica o "a − alvo −0,0 … +16,6 m" do A.6: com a janela nodal o semieixo
#     médio do estado convertido fica constante (−0,0 m em 30 dias).
# ===============================================================================================


@dataclass(frozen=True)
class MeanOscA6:
    epoch_utc: str = "2026-10-23T00:00:00+00:00"
    raan0_deg: float = 30.0
    psi0_deg: float = 0.0
    # Resíduos por iteração (A.6, janela kepleriana — ver discrepância acima)
    residual_a_m: tuple[float, ...] = (9.495e3, 7.446e1, 2.468e-1, 2.252e-4, 1.985e-6)
    residual_ex: tuple[float, ...] = (3.67e-6, 3.68e-6, 1.52e-8, 2.20e-11, 8.84e-14)
    residual_ey: tuple[float, ...] = (1.81e-3, 5.25e-6, 1.57e-9, 5.72e-11, 2.77e-13)
    residual_i_deg: tuple[float, ...] = (5.12e-3, 5.39e-6, 4.65e-9, 8.37e-11, 4.45e-13)
    first_residual_a_tol_km: float = 0.1  # critério do card: Δa da 1ª iteração ≈ 9,5 km
    # Osculador resultante (A.6, janela kepleriana)
    osc_a_m: float = 6878668.405
    osc_e: float = 7.270362e-4
    osc_i_deg: float = 97.444712
    osc_argp_deg: float = 270.0
    # Medido no satsim com a janela nodal (S4-04)
    nodal_osc_a_m: float = 6878679.4
    nodal_first_residual_a_m: float = 9483.6
    # Congelamento: estado convertido, Cowell J2–J4, médias por janela nodal (A.6)
    frozen_10d_a_minus_target_m: tuple[float, float] = (-0.0, 16.6)
    frozen_10d_e: tuple[float, float] = (1.0740e-3, 1.0781e-3)
    frozen_10d_argp_deg: tuple[float, float] = (90.00, 90.16)
    frozen_30d_e: tuple[float, float] = (1.0740e-3, 1.0837e-3)
    frozen_30d_argp_deg: tuple[float, float] = (90.00, 90.29)
    # Critérios do card (10 dias)
    crit_argp_deg: tuple[float, float] = (89.5, 90.5)
    crit_e: tuple[float, float] = (1.06e-3, 1.10e-3)
    crit_a_minus_target_m: tuple[float, float] = (-50.0, 50.0)
    # Contraste: médio usado como osculador (A.6)
    naive_10d_a_minus_target_m: tuple[float, float] = (9417.0, 9511.0)
    naive_10d_e: tuple[float, float] = (2.755e-3, 2.883e-3)
    naive_10d_argp_deg: tuple[float, float] = (68.39, 89.93)
    naive_10d_argp_range_deg: float = 21.5  # critério do card: > 15°
    naive_30d_argp_deg: tuple[float, float] = (20.67, 89.93)
    # Ω̇ (30 dias a 60 s) e deriva do LTAN (A.6)
    naive_raan_dot_deg_day: float = 0.977169
    naive_ltan_drift_min_year: float = -12.39
    converted_raan_dot_deg_day: float = 0.982540
    converted_ltan_drift_min_year: float = -4.54  # i = 97,4396° vem da teoria só J2 (S4-05)
    theory_j2_j4_raan_dot_deg_day: float = 0.983318  # nos médios medidos (Cowell −0,079%)
    theory_j2_raan_dot_deg_day: float = 0.985603  # (−0,31%)


MEAN_OSC_A6 = MeanOscA6()

# ===============================================================================================
# Semana 4 — Apêndice A.5 do prompt da Semana 4: taxas médias de 1ª ordem com J2, J3 e J4
#     (card S4-05). Órbita de projeto: a 6888,089 km, i 97,4396°, ω 90°, e 1,074e-3. Fórmulas
#     conferidas pelo autor do prompt contra a dedução simbólica em 2000 pontos.
# ===============================================================================================


@dataclass(frozen=True)
class MeanZonalA5:
    a_km: float = 6888.089
    i_deg: float = 97.4396
    argp_deg: float = 90.0
    e: float = 1.074e-3
    # Ω̇ [°/d]
    raan_dot_j2_deg_day: float = 0.985653
    raan_dot_j2_j4_deg_day: float = 0.983375
    raan_dot_j2_j3_j4_deg_day: float = 0.983369
    # ω̇ [°/d]
    argp_dot_j2_deg_day: float = -3.48711
    argp_dot_j2_j4_deg_day: float = -3.47711
    argp_dot_j2_j3_j4_deg_day: float = 0.01039  # 4 algarismos no A.5: comparar a 5e-6 °/d
    # ψ̇ + ω̇ [°/d]
    u_dot_j2_deg_day: float = 5459.9997
    u_dot_j2_j4_deg_day: float = 5460.0070
    u_dot_j2_j3_j4_deg_day: float = 5460.0071
    # Excentricidade congelada (ω̇ = 0 em ω = 90°)
    frozen_e_j2_j3: float = 1.074119e-3
    frozen_e_eq33: float = 1.074113e-3  # Sem4 eq. 33: −(J3/2J2)(R⊕/a)·sin i
    frozen_e_j2_j4: float = 1.077210e-3
    # Medido no satsim (S4-05): Ω̇ 0,98565339 / 0,98337502 / 0,98336881 °/d; ω̇ J2–J4
    # +0,0103926 °/d; econg J2+J3 1,0741194e-3, J2–J4 1,0772101e-3.


MEAN_ZONAL_A5 = MeanZonalA5()

# ===============================================================================================
# Semana 4 — Apêndice A.7 do prompt da Semana 4: SSO calibrada, modelo médio × Cowell e o
#     comparador real (card S4-05). Época 2026-10-23 00:00 UTC.
#
#     Diferença esperada no período nodal medido (S4-05): o A.7 mede Tnod = 5 696,6659 s no Cowell
#     (+4,7 ppm da teoria); o satsim mede 5 696,6796 s (+2,3 ppm). É o mesmo viés da janela
#     kepleriana do A.6 (ver MEAN_OSC_A6): com ela o estado convertido fica com o semieixo médio
#     real ≈ 11 m abaixo do alvo, o que encurta o período em 1,5·11/6 888 089 ≈ 2,4 ppm e
#     reproduz os 5 696,6659 s. Com a janela nodal (implementada) o equivalente ao longo da órbita
#     é ≈ 45 km em 30 dias (A.7: ≈ 92 km). O Ω̇ (+0,085%) e a i calibrada (+0,00005°) quase não
#     mudam.
# ===============================================================================================


@dataclass(frozen=True)
class SsoCalibratedA7:
    epoch_utc: str = "2026-10-23T00:00:00+00:00"
    # SSO de projeto calibrada no Cowell J2–J4 (e = econg J2–J4, ω 90°)
    e_frozen: float = 1.07721e-3
    i_first_order_deg: float = 97.4396
    raan_dot_first_order_cowell_deg_day: float = 0.982540
    i_calibrated_deg: float = 97.46325
    i_calibrated_tol_deg: float = 0.0005  # critério do card
    ltan_drift_tol_min_year: float = 0.1  # critério do card
    # Medido no satsim (S4-05): i calibrada 97,463303°; Cowell 30 dias: Ω̇ 0,9856473 °/d, LTAN
    # −6e-5 min/ano.
    # SSO de demonstração: LTAN 22h30 → Ω
    ltan_h: float = 22.5
    ra_sun_deg: float = 207.5235
    raan_demo_deg: float = 5.0235
    # Modelo médio J2–J4 de 1ª ordem × Cowell (SSO calibrada, Ω 5°, 30 dias)
    raan0_mean_vs_cowell_deg: float = 5.0
    raan_dot_ratio_minus_1: float = 8.43e-4
    raan_dot_ratio_tol: float = 1e-3  # critério do card: 0,1%
    tnod_cowell_s: float = 5696.6659  # janela kepleriana: ver nota acima
    tnod_theory_s: float = 5696.6926
    tnod_theory_j2_s: float = 5696.7003
    tnod_ppm_tol: float = 10.0  # critério do card
    along_track_30d_km: float = 92.0
    # Medido no satsim (S4-05): Ω̇ +0,0849%; Tnod Cowell 5 696,6796 s (+2,28 ppm); teoria
    # 5 696,69257 s; ≈ 45 km ao longo da órbita em 30 dias.
    tnod_cowell_satsim_s: float = 5696.6796
    # Arrasto médio: C_D·A/m 0,011, sem co-rotação, 10 dias (ȧ [m/d])
    cd_a_m: float = 0.011
    a_dot_cowell_spherical_m_day: float = -27.76
    a_dot_cowell_geodetic_m_day: float = -23.63
    a_dot_eq14_spherical_m_day: float = -29.78
    a_dot_eq14_geodetic_m_day: float = -25.42
    a_dot_eq14_at_mean_altitude_m_day: float = -29.69
    a_dot_ratio_tol: float = 0.10  # critério do card: modelo médio a 10% do Cowell
    # Medido no satsim (S4-05): Cowell −27,85 (esférica) e −23,71 (geodésica) m/d; modelo médio
    # −29,85 e −25,47 m/d (razão Cowell/médio 0,933 e 0,931).
    # Sentinel-2A: médias do SGP4 do TLE do A.1
    sentinel_a_m: float = 7164249.6
    sentinel_e: float = 1.142e-3
    sentinel_argp_deg: float = 90.5
    sentinel_i_real_deg: float = 98.56361
    sentinel_i_calibrated_deg: float = 98.56939
    sentinel_i_first_order_j2_deg: float = 98.54454
    sentinel_i_tol_deg: float = 0.01  # critério do card
    # Medido no satsim (S4-05): i calibrada 98,569401° (+0,00579° do real).


SSO_CALIBRATED_A7 = SsoCalibratedA7()

# ===============================================================================================
# Semana 4 — Apêndice A.8 do prompt da Semana 4: Cowell × SGP4 (card S4-06). Erro de posição 3D
#     [km] em 24 h / 72 h / 7 dias. (A) partindo do estado do SGP4 em t = 0; (B) estado inicial
#     (e C_D·A/m) ajustado por mínimos quadrados em 72 h, amostras a cada 600 s. DOP853 rtol
#     1e-10, atol 1e-3. TLEs: tests/data/sentinel2a_sample.tle (A.1) e iss_sample.tle (Semana 3).
#     A tabela completa está em scripts/compare_high_fidelity.py (A8_KM) e no relatório
#     reports/s4_alta_fidelidade/relatorio.md, com os valores medidos no satsim.
# ===============================================================================================


@dataclass(frozen=True)
class HighFidelityA8:
    # (B) com forças completas (J2–J4 + Sol/Lua + SRP + arrasto): referência e critério do card
    sentinel_full_b_km: tuple[float, float, float] = (0.06, 0.09, 1.11)
    sentinel_full_b_crit_km: tuple[float, float, float] = (0.2, 0.3, 3.0)
    iss_full_b_km: tuple[float, float, float] = (0.14, 0.29, 3.39)
    iss_full_b_crit_km: tuple[float, float, float] = (0.5, 1.0, 8.0)
    # C_D·A/m efetivo ajustado em (B) [m²/kg]
    sentinel_cd_a_m_full: float = 0.00208
    iss_cd_a_m_full: float = 0.00134
    # Teste rápido do card: Sentinel-2A, J2–J4, sem arrasto, janela de 24 h a 600 s
    quick_window_s: float = 24 * 3600.0
    quick_rms_crit_km: float = 0.2
    quick_err_24h_crit_km: float = 0.3
    # Medido no satsim (S4-06), teste rápido: RMS 0,0071 km, erro em 24 h 0,0076 km (48
    # avaliações, ≈ 5 s); RMS partindo do estado do SGP4 em t = 0: 0,923 km.
    # Ajuste sintético (card): estado perturbado de 1 km e 1 m/s → recuperar ≤ 1 m e ≤ 1 mm/s
    synthetic_pos_tol_m: float = 1.0
    synthetic_vel_tol_m_s: float = 1e-3


HIGH_FIDELITY_A8 = HighFidelityA8()


@dataclass(frozen=True)
class DemoSceneA7:
    """SSO de demonstração do A.7 na cena (S4-06): LTAN 22h30 em 2026-10-23 00:00 UTC, 24 h a 20 s.

    Medido no satsim (build_demo, Cowell completo): LTAN 22,4997 … 22,5026 h; β 23,407 … 23,431°;
    fração de eclipse 0,3603 (idênticos ao A.7/A.4).
    """

    start_utc: str = "2026-10-23T00:00:00Z"
    ltan_h: float = 22.5
    ltan_tol_h: float = 0.01  # critério do card
    ltan_range_h: tuple[float, float] = (22.4997, 22.5026)
    eclipse_fraction: float = 0.36
    eclipse_fraction_tol: float = 0.02  # critério do card
    eclipse_fraction_a4: float = 0.3603
    beta_range_deg: tuple[float, float] = (23.41, 23.43)


DEMO_SCENE_A7 = DemoSceneA7()

# ===============================================================================================
# Semana 5 — Apêndice A.1 do prompt da Semana 5: constelação de referência Epic of Sun (Sem4
#     Tabela 4, com a calibração da Semana 4). Época 2026-10-23 00:00 UTC. Experimento de
#     referência (03/11/2026, satsim no commit 779a6d4): valores de regressão, não independentes.
# ===============================================================================================


@dataclass(frozen=True)
class ConstellationA1:
    epoch_utc: str = "2026-10-23T00:00:00+00:00"
    a_design_m: float = 6888089.0  # Sem4 §7.5
    a_repeat_calibrated_m: float = 6888108.4  # experimento (A.4)
    e_frozen: float = 1.07721e-3  # A.5 da Semana 4; ω = 90°
    argp_deg: float = 90.0
    i_calibrated_deg: float = 97.46325  # A.7 da Semana 4 (Cowell J2–J4)
    raan_deg: float = 5.0235  # A.7 da Semana 4 (LTAN 22h30)
    ltan_h: float = 22.5
    psi0_deg: tuple[float, ...] = (0.0, 300.0, 240.0, 180.0, 120.0, 60.0)  # Sem4 Tabela 4 / eq. 42
    node_delay_s: tuple[float, ...] = (0.0, 949.5, 1898.9, 2848.4, 3797.8, 4747.3)  # Tabela 4
    node_delay_tol_s: float = 2.0  # critério do card S5-01
    # Forças do Cowell completo (A.7 da Semana 4): J2–J4 + Sol/Lua + SRP + arrasto geodésico
    srp_cr: float = 1.3
    srp_area_over_mass: float = 0.01
    drag_cd_a_m: float = 0.011
    # Experimento: 6 conversões médio → osculador com i dada
    build_time_s: float = 0.33
    build_time_limit_s: float = 3.0  # critério do card S5-01
    # Primeiros nós ascendentes (Cowell completo, a de projeto) [s desde a época], S1 … S6
    first_node_s: tuple[float, ...] = (4274.46, 5223.99, 476.60, 1426.12, 2375.65, 3324.93)
    delay_s2_s1_s: float = 949.53
    # Medido no satsim (S5-01, nós por interpolação linear de z_ECI com passo 10 s): 4 274,463;
    # 5 223,988; 476,596; 1 426,122; 2 375,648; 3 324,934 s; construção dos 6 satélites 0,33 s.


CONSTELLATION_A1 = ConstellationA1()

# ===============================================================================================
# Semana 5 — Apêndice A.2 do prompt da Semana 5: órbitas de repetição só J2 (Sem4 Tabela 3,
#     e = 0, eqs. 27 + 37) e grade no modelo médio (Sem4 §9.2). Teoria (independente).
# ===============================================================================================


@dataclass(frozen=True)
class RepeatOrbitA2:
    n_revs: int
    d_days: int
    h_km: float
    i_deg: float
    tnod_s: float
    grid_deg: float
    grid_km: float  # no equador


REPEAT_ORBITS_A2: tuple[RepeatOrbitA2, ...] = (
    RepeatOrbitA2(91, 6, 509.95, 97.4396, 5696.70, 3.956, 440.4),
    RepeatOrbitA2(106, 7, 517.19, 97.4671, 5705.66, 3.396, 378.1),
    RepeatOrbitA2(121, 8, 522.62, 97.4879, 5712.40, 2.975, 331.2),
    RepeatOrbitA2(89, 6, 612.99, 97.8388, 5824.72, 4.045, 450.3),
    RepeatOrbitA2(143, 10, 786.12, 98.5446, 6041.96, 2.517, 280.2),
    RepeatOrbitA2(233, 16, 699.59, 98.1863, 5933.05, 1.545, 172.0),
)
"""Medido no satsim (S5-01, ``repeat_orbit``): h 509,9523 / 517,1863 / 522,6243 / 612,9856 /
786,1213 / 699,5852 km; i 97,43957 / 97,46711 / 97,48786 / 97,83883 / 98,54460 / 98,18629°."""


@dataclass(frozen=True)
class RepeatOrbitTolA2:
    h_km: float = 0.01  # critério do card S5-01
    i_deg: float = 0.001  # critério do card S5-01
    tnod_s: float = 0.01  # arredondamento do A.2 (2 casas)


REPEAT_ORBIT_TOL_A2 = RepeatOrbitTolA2()


@dataclass(frozen=True)
class MeanGridA2:
    """Grade no modelo médio (Sem4 §9.2): nós a cada Tnod, longitude −(ω⊕ − Ω̇)·t."""

    grid_ideal_deg: float = 3.9560
    nodes_p1_6days: int = 91
    nodes_p6_1day: int = 91
    nodes_p7_1day: int = 107  # depende de onde a janela corta (Cowell, [0, 24 h): 106)
    step_p7_min_deg: float = 0.5651
    step_p7_max_deg: float = 3.3909


MEAN_GRID_A2 = MeanGridA2()

# ===============================================================================================
# Semana 5 — Apêndice A.3 do prompt da Semana 5: nós ascendentes e grade no Cowell completo
#     (experimento de referência, constelação do A.1 com o a de projeto, janela [0, 24 h) salvo
#     indicação). Ideal: 360°/91 = 3,95604°. Afastamento ancorado no 1º nó de S1.
# ===============================================================================================


@dataclass(frozen=True)
class NodeGridA3:
    grid_ideal_deg: float = 3.95604
    # P = 6, 1 dia
    p6_nodes: int = 91
    p6_gap_min_deg: float = 3.95173
    p6_gap_max_deg: float = 3.96431
    p6_gap_tol_deg: float = 0.002  # critério do card S5-02
    p6_max_offset_deg: float = 0.0064
    p6_max_offset_limit_deg: float = 0.01  # critério do card S5-02
    # P = 1, 6 dias (91 primeiros nós)
    p1_6d_nodes: int = 91
    p1_6d_gap_min_deg: float = 3.9169
    p1_6d_gap_max_deg: float = 3.9674
    # P = 7 na mesma órbita (Δψ = 51,43°)
    p7_nodes: int = 106
    p7_gap_min_deg: float = 3.3834
    p7_gap_max_deg: float = 3.9585
    p7_max_offset_deg: float = 1.694
    p7_max_offset_limit_deg: float = 1.0  # critério do card: > 1°
    # Precisão do nó (S1, 15 nós em 24 h): interpolação linear × brentq em z_ECI = 0
    node_err_60s_deg: float = 7.3e-5
    node_err_60s_s: float = 1.1e-3
    node_err_30s_deg: float = 9.3e-6
    node_err_30s_s: float = 3e-4
    node_err_10s_deg: float = 4.8e-7
    node_err_30s_limit_deg: float = 2e-5  # critério do card S5-02
    node_err_30s_limit_s: float = 1e-3  # critério do card S5-02
    # Período nodal medido no Cowell completo em 6 dias (ajuste linear) e alvo da repetição
    tnod_cowell_full_6d_s: float = 5696.590
    tnod_repeat_target_s: float = 5696.704
    # Medido no satsim (S5-02, node_crossings a 30 s): P = 6: 91 nós, 3,951728 … 3,964321°,
    # afastamento 0,006447° (âncora S1; com a âncora no 1º nó no tempo, o de S3, 0,010672°);
    # P = 7: 106 nós, 3,383367 … 3,958479°, afastamento 1,694° (âncora S1; 1,702° no tempo).
    # Precisão: 60 s 7,31e-5° / 1,12e-3 s; 30 s 9,26e-6° / 2,7e-4 s; 10 s 5,17e-7° / 2,9e-5 s.
    # Os desvios à grade têm padrão de período 3 (s1/s4, s2/s5, s3/s6; ±0,0045°): o termo de
    # curto período de J2 em Ω e u, que a conversão médio → osculador (S4-04) não corrige.


NODE_GRID_A3 = NodeGridA3()

# ===============================================================================================
# Semana 5 — Apêndice A.4 do prompt da Semana 5: fechamento do ciclo de 6 dias e semieixo de
#     repetição calibrado (experimento; S1 da constelação do A.1, Cowell).
# ===============================================================================================


@dataclass(frozen=True)
class RepeatClosureA4:
    a_design_m: float = 6888089.0
    # J2–J4, a de projeto
    closure_design_j2j4_deg: float = 0.00895
    closure_design_j2j4_km: float = 1.00
    closure_design_tol_km: float = 0.2  # critério do card S5-02
    residual_design_s: float = -2.197  # 91·Tnod − 6·Tdia,nod
    tnod_design_s: float = 5696.6796
    # J2–J4, a calibrado
    a_calibrated_m: float = 6888108.37
    delta_a_m: float = 19.37
    delta_a_tol_m: float = 1.0  # critério do card S5-02
    closure_calibrated_j2j4_deg: float = -0.00020
    closure_calibrated_j2j4_km: float = -0.022
    closure_calibrated_limit_km: float = 0.1  # critério do card S5-02
    secant_residuals_s: tuple[float, float] = (1e-2, 3e-9)  # 2º ponto → secante
    # Cowell completo (com arrasto)
    closure_design_full_deg: float = 0.0422
    closure_design_full_km: float = 4.69
    closure_calibrated_full_deg: float = 0.0330
    closure_calibrated_full_km: float = 3.68
    # Ω medido (6 dias, 60 s) e dia nodal
    raan_dot_design_deg_day: float = 0.985632
    raan_dot_calibrated_deg_day: float = 0.985622
    t_day_nod_s: float = 86400.0065
    # Sem4 eq. 32: efeito de +19 m no LTAN
    ltan_drift_delta_a_min_year: float = -0.014
    # Medido no satsim (S5-02): ver o relatório reports/s5_constelacao/relatorio.md e o teste
    # tests/test_repeat_grid.py. J2–J4: a calibrado 6 888 108,374 m (Δa +19,374 m em relação a
    # 6 888 089 m), f: −2,1686 s (a0 = repeat_orbit, 6 888 089,253 m) → 1,0e-2 s → 1,4e-8 s;
    # fechamento +0,008952° (0,996 km) com o a de projeto e −0,000201° (−0,022 km) com o
    # calibrado. A maior parte do resíduo do calibrado vem de ω⊕ = 7,292115e-5 rad/s (usado em
    # Tdia,nod, como no card) ser 8,55e-12 rad/s menor que a taxa real do GMST que gira o ECEF
    # (7,2921158553e-5): em 6 dias, −2,54e-4° (−28 m) previstos; o resto (≈ +0,5e-4°) é o
    # ruído do ajuste de Tnod contra o nó 91 medido.


REPEAT_CLOSURE_A4 = RepeatClosureA4()

# ===============================================================================================
# Semana 5 — Apêndice A.5 do prompt da Semana 5: cobertura (Sem4 eq. 45 e experimento: Cowell
#     completo, a de projeto, η = 23,2°, alvos em |lat| ≤ 80°, grade 1°, passo 10 s; acessos a
#     menos de 30 min contam como uma passagem; o intervalo inclui as bordas da janela).
# ===============================================================================================


@dataclass(frozen=True)
class CoverageA5:
    a_m: float = 6888089.0
    # Sem4 eq. 45: rolagem η [°] → λ [°] e alcance transversal [km]
    roll_deg: tuple[float, ...] = (10.0, 20.0, 23.2, 30.0)
    swath_half_angle_deg: tuple[float, ...] = (0.8088, 1.6765, 1.9784, 2.6820)
    swath_km: tuple[float, ...] = (90.0, 186.6, 220.2, 298.6)
    swath_km_sem4: tuple[float, ...] = (90.0, 187.0, 220.0, 299.0)  # Sem4 §8.6
    swath_km_tol: float = 1.0  # critério do card S5-03
    roll_half_grid_deg: float = 23.1957
    roll_half_grid_tol_deg: float = 0.01  # critério do card S5-03
    # Experimento (grade 1°)
    roll_coverage_deg: float = 23.2
    fraction_p6_1d: float = 1.0
    fraction_p1_1d: float = 0.4848
    fraction_p1_1d_tol: float = 0.02  # critério do card S5-03
    fraction_p1_6d: float = 1.0
    revisit_p6_3d_lat60_h: float = 13.97
    revisit_p6_3d_lat60_limit_h: float = 15.0  # critério do card S5-03
    revisit_p6_3d_lat80_h: float = 19.04
    fraction_p6_3d_sun10: float = 0.9255
    revisit_p6_3d_sun10_lat60_h: float = 24.26
    revisit_p6_3d_sun10_lat60_limit_h: float = 25.0  # critério do card S5-03
    min_sun_elevation_deg: float = 10.0
    gap_merge_s: float = 1800.0
    # Medido no satsim (S5-03), grade 1°: λ 0,808783 / 1,676505 / 1,978441 / 2,682043°; meia
    # grade 23,195673°; P = 1 em 1 dia 48,82% (A.5: 48,48%); P = 6 3 dias: revisita 13,97 h
    # (|lat| ≤ 60°) e 19,04 h (≤ 80°); com Sol ≥ 10°: 92,55% e 24,26 h. O A.5 se reproduz
    # medindo o intervalo entre os INÍCIOS de passagens sucessivas (coverage_stats(...,
    # interval="revisit")); do fim de uma ao início da seguinte (interval="blackout") dá 13,69 h /
    # 16,68 h / 24,00 h: a diferença é a duração da passagem, que em latitude alta junta vários
    # satélites do trem (15,8 min entre eles).


COVERAGE_A5 = CoverageA5()


@dataclass(frozen=True)
class CoverageTestS5:
    """Referência própria dos testes de cobertura (S5-03): janelas de 3 e 6 dias com o Cowell
    só J2–J4 (≈ 7× mais rápido que o completo) e grade de 2°.

    Medido no satsim (S5-03): em 3 dias a cobertura no Cowell J2–J4 é idêntica à do completo a
    menos de 0,01 h (grade 1°: 13,97 / 19,04 / 24,26 h e 92,55% nos dois).
    """

    grid_deg: float = 2.0
    revisit_p6_3d_lat60_h: float = 13.967
    revisit_p6_3d_lat80_h: float = 19.039
    blackout_p6_3d_lat60_h: float = 13.694
    fraction_p6_3d_sun10: float = 0.9259
    revisit_p6_3d_sun10_lat60_h: float = 24.258
    regression_tol_h: float = 0.02  # ≈ 7 amostras de 10 s
    fraction_p1_1d_full_2deg: float = 0.4922


COVERAGE_TEST_S5 = CoverageTestS5()

# ===============================================================================================
# Semana 5 — Apêndice A.6 do prompt da Semana 5: passagens da constelação do A.1 (a de projeto,
#     Cowell completo) sobre São Bento do Sapucaí, 2026-10-23 00:00 UTC + 24 h, el ≥ 10°,
#     find_passes (experimento de referência).
# ===============================================================================================


@dataclass(frozen=True)
class ContactPassA6:
    sat: str
    aos_utc: str  # HH:MM:SS em 2026-10-23
    max_el_deg: float
    duration_s: float


CONTACT_PASSES_A6: tuple[ContactPassA6, ...] = (
    ContactPassA6("S5", "00:32:20", 12.83, 209.9),
    ContactPassA6("S6", "00:46:51", 21.42, 351.3),
    ContactPassA6("S1", "01:01:59", 36.04, 418.5),
    ContactPassA6("S2", "01:17:23", 65.20, 448.4),
    ContactPassA6("S3", "01:33:01", 69.55, 449.3),
    ContactPassA6("S4", "01:48:51", 38.48, 422.2),
    ContactPassA6("S5", "02:04:57", 22.98, 361.5),
    ContactPassA6("S6", "02:21:30", 14.08, 242.8),
    ContactPassA6("S3", "12:05:18", 12.27, 188.7),
    ContactPassA6("S4", "12:19:37", 20.13, 337.0),
    ContactPassA6("S5", "12:34:35", 33.21, 409.5),
    ContactPassA6("S6", "12:49:54", 59.02, 444.8),
    ContactPassA6("S1", "13:05:28", 76.87, 451.3),
    ContactPassA6("S2", "13:21:18", 42.03, 429.8),
    ContactPassA6("S3", "13:37:26", 24.53, 374.1),
    ContactPassA6("S4", "13:54:05", 14.80, 261.5),
)


@dataclass(frozen=True)
class ContactPlanA6:
    start_utc: str = "2026-10-23T00:00:00+00:00"
    duration_s: float = 86400.0
    el_min_deg: float = 10.0
    n_passes: int = 16
    n_passes_tol: int = 1  # critério do card S5-04
    passes_per_sat: tuple[int, ...] = (2, 2, 3, 3, 3, 3)  # S1 … S6
    total_contact_min: float = 96.68
    total_contact_tol_min: float = 3.0  # critério do card S5-04
    max_elevation_deg: float = 76.87
    n_overlaps: int = 0
    n_trains: int = 2
    # Trens (BRT = UTC−3): noturno 00:32–02:25 UTC = 21h32–23h25 BRT;
    # diurno 12:05–13:58 UTC = 9h05–10h58 BRT
    night_train_utc: tuple[str, str] = ("00:32", "02:25")
    day_train_utc: tuple[str, str] = ("12:05", "13:58")
    # Dentro de um trem: AOS consecutivos a 14–17 min (≈ Tnod/6 = 15,8 min + geometria)
    train_aos_spacing_min: tuple[float, float] = (14.0, 17.0)  # critério do card S5-04
    # Medido no satsim (S5-04, contact_plan): as 16 passagens idênticas ao A.6 (AOS no segundo,
    # el máx a 0,005°, duração a 0,05 s); contato total 96,68 min; trens 21:32–23:25 e
    # 09:05–10:58 BRT; Sol no TCA de −42,9° a −55,5° (noite) e de 50,2° a 73,5° (dia); maior
    # intervalo sem contato 10,03 h (do último LOS ao fim da janela).


CONTACT_PLAN_A6 = ContactPlanA6()

# ===============================================================================================
# Semana 5 — Apêndice A.7 do prompt da Semana 5: manutenção (Sem4 §6.1, §7.6 e experimento).
#     λ(t) = 3π|ȧ|t²/(2a) (Sem4 §7.6); v = 7 607,1 m/s; Δv = v·δa/(2a). Banda do traço ±5 km.
# ===============================================================================================


@dataclass(frozen=True)
class TrackMaintenanceRow:
    a_dot_m_day: float
    days_to_5km: float
    delta_a_m: float
    delta_v_per_maneuver_m_s: float
    maneuvers_per_year: float
    delta_v_per_year_m_s: float
    source: str


TRACK_MAINTENANCE_A7: tuple[TrackMaintenanceRow, ...] = (
    TrackMaintenanceRow(-23.71, 6.95, 164.8, 0.0910, 52.5, 4.78, "Cowell (A.7 da Semana 4)"),
    TrackMaintenanceRow(-25.47, 6.71, 170.8, 0.0943, 54.5, 5.14, "modelo médio, mesmo caso"),
    TrackMaintenanceRow(-35.0, 5.72, 200.3, 0.1106, 63.8, 7.06, "Sem4 §6.2"),
)


@dataclass(frozen=True)
class MaintenanceA7:
    a_m: float = 6888089.0
    i_deg: float = 97.46325
    v_m_s: float = 7607.1
    deadband_km: float = 5.0
    rel_tol: float = 0.02  # critério do card S5-05 ("A.7, rel 2%")
    # Verificação no Cowell: a calibrado + arrasto, fechamento de 6 dias (A.4) × fórmula
    closure_cowell_6d_km: float = 3.68
    closure_formula_6d_km: float = 3.72
    closure_formula_tol: float = 0.10  # critério do card S5-05: fórmula a ≤ 10% do Cowell
    # LTAN (Sem4 eq. 32)
    ltan_per_di_001deg_min_year: float = 1.92  # Sem4: +1,9
    ltan_per_dh_1km_min_year: float = -0.73
    di_dt_sun_deg_day: float = -9.8e-5  # S4-02
    di_dt_sun_deg_year: float = -0.0358
    ltan_accel_min_year2: float = -6.87  # → −3,43·t² min (t em anos)
    ltan_quadratic_coeff_min: float = -3.43
    ltan_band_min: float = 5.0
    years_to_ltan_band: float = 1.21
    inclination_delta_v_per_year_m_s: float = 4.75
    # Medido no satsim (S5-05): fórmulas idênticas ao A.7 (6,95 d, 164,8 m, 0,0910 m/s, 52,5/ano,
    # 4,78 m/s/ano; −6,867 min/ano², 1,207 ano, 4,752 m/s/ano; 1,919 e −0,732 min/ano).
    # ȧ da constelação de referência por janelas nodais (Cowell completo): −23,95 m/dia em 1 dia,
    # −24,01 em 6 dias (A.7: −23,71; 1,0–1,3%). Fechamento de 6 dias com o a calibrado e arrasto:
    # Cowell completo 3,673 km × fórmula 3,772 km (2,7%). Validação do LTAN no Cowell J2–J4 + Sol
    # (60 dias): di/dt −9,82e-5 °/d, aceleração do LTAN −7,13 min/ano² × fórmula −6,88 (4%).


MAINTENANCE_A7 = MaintenanceA7()


@dataclass(frozen=True)
class MaintenanceTestS5:
    """Referência própria do teste Cowell × fórmula (S5-05): S1 com o a calibrado, Cowell só
    J2–J4 + arrasto (≈ 2× mais rápido que o completo), 6 dias.

    Medido no satsim (S5-05): fechamento 3,838 km; ȧ em 6 dias −24,03 m/dia; fórmula 3,774 km
    (1,7%). No Cowell completo: 3,673 km e 3,772 km (2,7%).
    """

    a_calibrated_m: float = 6888108.374
    closure_6d_km: float = 3.838
    formula_6d_km: float = 3.774
    regression_tol_km: float = 0.01
    a_dot_6d_m_day: float = -24.03


MAINTENANCE_TEST_S5 = MaintenanceTestS5()


@dataclass(frozen=True)
class DemoSceneS5:
    """Constelação na cena da demo (S5-06): critérios do card e medidas de regressão.

    Critérios do card S5-06: 6 membros em ``meta.constellation``; 91 nós somados em 1 dia (A.3,
    ``NODE_GRID_A3.p6_nodes``); LTAN dos 6 satélites constante em 24 h, 22,500 ± 0,01 h, e iguais
    entre si (mesma tolerância); ``--days 6`` com arquivo < 15 MB.

    Medido no satsim (S5-06, build_demo --offline, 2026-10-23 00:00 UTC, Cowell completo, a
    calibrado): LTAN 22,4994 … 22,5029 h nos seis, diferença máxima entre eles 0,0006 h; 15 ou 16
    nós por satélite (91 no total); faixa λ = 1,9955° (rolagem de meia grade do A.5, 23,1957°, no
    raio médio da órbita); ȧ −23,94 … −24,01 m/dia. Arquivo: 1 dia (passo 20 s) 6,7 MB; 6 dias
    (passo 60 s) 11,2 MB, 8 641 amostras por satélite, 40 s para gerar.
    """

    start_utc: str = "2026-10-23T00:00:00Z"
    member_ids: tuple[str, ...] = ("s1", "s2", "s3", "s4", "s5", "s6")
    nodes_1d: int = 91  # A.3
    ltan_h: float = 22.5
    ltan_tol_h: float = 0.01  # critério do card
    ltan_spread_tol_h: float = 0.01  # critério do card ("iguais entre si")
    max_file_mb: float = 15.0  # critério do card
    step_6d_s: float = 60.0


DEMO_SCENE_S5 = DemoSceneS5()

# ===============================================================================================
# Semana 6 — Apêndice A do prompt da Semana 6: constelação de 18 satélites (P = 3·D) na órbita
#     91/6. Experimento de referência em 09/11/2026 com o satsim no commit e423447 (fim da S5-06):
#     sso_constellation(2026-10-23 00:00 UTC, Cowell completo do build_demo, n_sats=18,
#     a = 6 888 089 m de projeto, i = 97,46325°), LTAN 22h30. Valores de regressão e de ordem de
#     grandeza, não independentes do satsim; os da teoria (A.2: Sem4 eqs. 42, 43, 45) são.
#     Convenção: ψ0_j = −j·20°, ids s1 … s18 em ordem de fase; família f = j mod 3.
# ===============================================================================================


@dataclass(frozen=True)
class DesignRowA2S6:
    n_sats: int
    delta_u_deg: float
    tracks_per_day: int | None
    daily_grid_deg: float | None
    daily_grid_km: float | None  # no equador
    roll_daily_coverage_deg: float | None
    uniform: bool


DESIGN_A2_S6: tuple[DesignRowA2S6, ...] = (
    DesignRowA2S6(6, 60.0, 91, 3.95604, 440.0, 23.1957, True),  # P = D
    DesignRowA2S6(18, 20.0, 273, 1.31868, 147.0, 8.1835, True),  # P = 3D
    DesignRowA2S6(7, 51.43, None, None, None, None, False),  # A.2/A.3 da S5
)
"""A.2 do prompt da Semana 6 (teoria, Sem4 eqs. 42, 43 e 45; rolagem pela eq. 45 com
a = 6 888,089 km e λ = π/273). Critérios do card S6-01: rolagem 23,1957° (P = 6) e 8,1835° ±
0,01° (P = 18); P = 3 uniforme por P | D, com revisita da grade em 2 dias. Com η = 23,2°,
λ = 1,9784° = 1,5 grade diária de 18. Medido no satsim (S6-01, ``ConstellationDesign``):
23,195673° e 8,183479°."""

ROLL_DAILY_COVERAGE_TOL_DEG_S6 = 0.01  # critério do card S6-01


@dataclass(frozen=True)
class Constellation18A1S6:
    """A.1 do prompt da Semana 6: constelação de 18 e grade no Cowell completo (experimento).

    Medido no satsim (S6-01, ``DESIGN_91_6_P18.build``, nós por ``node_crossings`` a 30 s):
    construção 1,02 s; nós de 1 dia 21,3 s; 1ºs nós 4 274,46 / 4 591,00 / 4 907,52 / 5 223,99 s;
    atraso s2 − s1 316,542 s (vizinhos: 316,42 … 316,54 s); P = 18: 273 nós, 1,31556 …
    1,32196°, afastamento 0,00695° (âncora no 1º nó de s1); família 0: 91 nós, 3,95173 …
    3,96432°, afastamento 0,00645°; famílias 1 e 2: 91 nós cada, 3,9482 … 3,9625°.
    Modelo médio (MeanJ2, e = 0, a e i de ``repeat_orbit(91, 6)``, janela [−1 s, 1 dia)): 273 nós,
    1,318641 … 1,318698° (ideal 1,318681°).
    """

    n_sats: int = 18
    build_time_s: float = 0.98
    nodes_1d_time_s: float = 20.3
    node_delay_s2_s1_s: float = 316.54  # Tnod/18
    node_delay_tol_s: float = 1.0  # critério do card S6-01
    first_node_s: tuple[float, ...] = (4274.46, 4591.00, 4907.52, 5223.99)  # s1 … s4
    p18_nodes: int = 273
    p18_gap_min_deg: float = 1.31556
    p18_gap_max_deg: float = 1.32196
    gap_tol_deg: float = 0.002  # critério do card S6-01
    p18_max_offset_deg: float = 0.00695
    p18_max_offset_limit_deg: float = 0.01  # critério do card S6-01
    family0_nodes: int = 91  # = A.3 da Semana 5
    family0_gap_min_deg: float = 3.95173
    family0_gap_max_deg: float = 3.96432
    family0_max_offset_deg: float = 0.00645
    mean_grid_tol_deg: float = 1e-4  # arredondamento (5e-5°) + deriva do GMST em 1 dia (4,2e-5°)


CONSTELLATION18_A1_S6 = Constellation18A1S6()


@dataclass(frozen=True)
class Coverage18A3S6:
    """A.3 do prompt da Semana 6: cobertura e revisita (experimento; grade de 2° × 2°,
    |lat| ≤ 80°, passo 10 s; acessos a menos de 30 min contam como uma passagem; intervalo de
    início a início, com as bordas da janela).

    Medido no satsim (S6-01, grade 2°, passo 10 s), Cowell completo: P = 18 com 8,1835°: 100% em
    1 dia, 100% e 24,00 h em 3 dias; P = 18 com 23,2° em 3 dias: 100%, 13,794 h (|lat| ≤ 60°) e
    18,953 h (≤ 80°); P = 6 (família 0) com 23,2°: 13,967 h e 19,039 h; P = 6 com 8,1835° em 3 dias:
    74,348%, 72,0 h. No Cowell só J2–J4 (o dos testes de 3 dias): os mesmos valores, exceto P = 6
    com 8,1835°: 73,896%. 17 satélites (sem s1 ou sem s5) com 23,2°: 100% em 1 dia. Passagens por
    alvo em 1 dia (P = 18, 23,2°, |lat| ≤ 60°, fusão de 30 min): mín. 2, mediana 2, máx. 3;
    passagens por satélite (sem fundir satélites do trem): mín. 6, mediana 7, máx. 13.
    """

    grid_deg: float = 2.0
    step_s: float = 10.0
    gap_merge_s: float = 1800.0
    roll_p18_deg: float = 8.1835
    roll_p6_deg: float = 23.2
    fraction_p18_roll8_1d: float = 1.0
    fraction_p18_roll8_3d: float = 1.0
    revisit_p18_roll8_3d_lat60_h: float = 24.0
    revisit_p18_roll23_3d_lat60_h: float = 13.79
    revisit_p6_roll23_3d_lat60_h: float = 13.97
    revisit_tol_h: float = 0.5  # critério do card S6-01 (e diferença P = 18 × P = 6 < 0,5 h)
    fraction_p6_roll8_3d: float = 0.7435
    fraction_p6_roll8_3d_tol: float = 0.02  # critério do card S6-01
    fraction_p17_roll23_1d: float = 1.0
    passes_per_target_1d_min: int = 2
    passes_per_target_1d_median: int = 2
    passes_per_target_1d_max: int = 3
    train_opportunities_window_min: float = 10.0


COVERAGE18_A3_S6 = Coverage18A3S6()


@dataclass(frozen=True)
class Contacts18A4S6:
    """A.4 do prompt da Semana 6: contatos com São Bento do Sapucaí (experimento; 24 h desde
    2026-10-23 00:00 UTC, el ≥ 10°, ``contact_plan`` com step_s = 20 s).

    Medido no satsim (S6-01): 49 passagens (2 a 4 por satélite), 287,99 min, 2 trens (noite
    00:27:58–02:30:08 UTC, 24 passagens; dia 12:01:10–14:07:20 UTC, 25), AOS consecutivos a
    4,37–5,77 min (noite) e 4,14–6,28 min (dia), 36 sobreposições, no máximo 2 satélites
    simultâneos, maior intervalo sem contato 9,878 h, el. máx. 87,57°; 0,3 s de cálculo.
    """

    n_passes: int = 49
    n_passes_tol: int = 2  # critério do card S6-01
    total_contact_min: float = 287.99
    total_contact_tol_min: float = 6.0  # critério do card S6-01
    n_trains: int = 2
    night_train_utc: tuple[str, str] = ("00:28", "02:30")
    day_train_utc: tuple[str, str] = ("12:01", "14:07")
    night_train_passes: int = 24
    day_train_passes: int = 25
    train_aos_spacing_min: tuple[float, float] = (4.0, 6.5)  # critério do card S6-01
    a4_aos_spacing_min: tuple[float, float] = (4.1, 6.3)  # A.4 (≈ Tnod/18 = 5,3 min)
    n_overlaps: int = 36
    n_overlaps_tol: int = 3  # critério do card S6-01
    max_simultaneous: int = 2
    max_gap_h: float = 9.88
    passes_per_sat_range: tuple[int, int] = (2, 4)
    p6_n_passes: int = 16  # A.6 da Semana 5
    p6_total_contact_min: float = 96.68


CONTACTS18_A4_S6 = Contacts18A4S6()


@dataclass(frozen=True)
class Maintenance18S6:
    """Manutenção dos 18 (critério do card S6-01): mesmo a, mesmo arrasto → ``MaintenanceStatus``
    igual ao dos 6 a menos de 1% (ȧ medido em 1 dia).

    Medido no satsim (S6-01, ``measure_decay_rate`` em 1 dia, Cowell completo): ȧ dos 18 de
    −24,022 a −23,945 m/dia (dispersão 0,077 m/dia, 0,3%), contra −23,95 dos 6 (A.7 da S5).
    """

    rel_tol: float = 0.01  # critério do card S6-01
    a_dot_min_m_day: float = -24.022
    a_dot_max_m_day: float = -23.945


MAINTENANCE18_S6 = Maintenance18S6()


@dataclass(frozen=True)
class DemoSceneS6:
    """Demo com 6 ou 18 satélites e variantes (S6-02): critérios do card e medidas.

    Critérios do card S6-02: cada arquivo < 15 MB; flecha da interpolação do viewer ≤ 5 km;
    faixa da cena de 18 com a rolagem de cobertura diária, 8,1835° (A.2 da Semana 6).

    Medido no satsim (S6-02, ``build_demo --all --offline``, Cowell completo, 20 objetos com 18):
    index.html 6,7 MB (20 s); ciclo_6dias.html 11,2 MB (60 s); constelacao18.html 11,4 MB (20 s,
    precisão compacta); ciclo_6dias_18.html 11,6 MB (120 s, compacta); 2 min 45 s para os quatro.
    Com precisão completa, 18 satélites a 20 s em 1 dia: 13,4 MB; 6 dias a 60 s: 24,5 MB (compacta:
    20,6 MB; a 90 s, 14,5 MB). A estimativa do A.5 (17 MB em 1 dia a 20 s) é proporcional a 8
    objetos e ignora a parte fixa (2,3 MB de globe.gl e textura).
    Interpolação, S1 no Cowell completo, 1 dia, erro máximo entre amostras: linear em lat/lon (o
    método das Semanas 3–5) 3,28 / 7,37 / 13,06 / 29,16 / 51,3 / 79,1 / 112,1 km a 20 / 30 / 40 /
    60 / 80 / 100 / 120 s, concentrado perto da latitude máxima (a flecha r·θ²/8 do A.5 vale para a
    corda em ECI, não para lat/lon); Hermite cúbica no vetor unitário (S6-02): 0,031 / 0,070 /
    0,124 / 0,279 / 0,497 / 0,778 / 1,127 km (o máximo vem das pontas da janela, com tangente
    unilateral; no meio da janela, 0,27 km a 120 s).
    """

    max_file_mb: float = 15.0  # critério do card
    interp_limit_km: float = 5.0  # critério do card
    roll18_deg: float = 8.1835  # A.2 da Semana 6
    file_mb: tuple[tuple[str, float], ...] = (
        ("index.html", 6.7),
        ("ciclo_6dias.html", 11.2),
        ("constelacao18.html", 11.4),
        ("ciclo_6dias_18.html", 11.6),
    )
    interp_error_120s_km: float = 1.127


DEMO_SCENE_S6 = DemoSceneS6()
