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
