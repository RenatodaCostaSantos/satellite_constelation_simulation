"""Relatório da constelação 91/6 da Semana 5 (reports/s5_constelacao/).

Uso (na raiz do repositório, com o pacote instalado com o extra dev):

    python scripts/constellation_report.py                  # tudo (~2 min) e o relatório
    python scripts/constellation_report.py --from-json      # só refaz gráficos e relatório
    python scripts/constellation_report.py --only coverage  # recalcula uma seção, mantém as outras

Seções (uma por card; os próximos cards acrescentam as suas):

- **S5-02 · Grade e repetição.** Nós ascendentes da constelação de referência (Apêndice A.1 do
  prompt da Semana 5) no Cowell completo em 1 dia para P = 6 e P = 7 (equivalente à Figura 12 do
  Sem4), a mesma grade no modelo médio (Sem4 §9.2), a precisão do nó por interpolação linear, o
  fechamento do ciclo de 6 dias com o a de projeto e com o a calibrado
  (``repeat_semi_major_axis_calibrated``), com e sem arrasto, e o efeito de Δa no LTAN.
- **S5-03 · Cobertura e revisita.** Rolagem de ±23,2° (meia grade, Sem4 §8.6), alvos em
  |lat| ≤ 80° numa grade de 1°, passo de 10 s, Cowell completo: fração coberta por latitude em
  1 dia (P = 1 e P = 6, com e sem a restrição de luz), P = 1 em 6 dias e o mapa do intervalo
  máximo de revisita da constelação em 3 dias (sem e com luz).
- **S5-04 · Plano de contatos.** Passagens dos 6 satélites sobre São Bento do Sapucaí em 24 h
  (el ≥ 10°, ``contact_plan``): linha do tempo, trens, dia/noite, contato total e sobreposições.
- **S5-05 · Manutenção.** ȧ medido no Cowell, ``MaintenanceStatus`` dos 6 satélites, deriva do
  traço × tempo para 3 coeficientes balísticos (como a Figura 10 do Sem4) com a curva do Cowell
  completo, e a deriva do LTAN em 2 anos validada no Cowell J2–J4 + Sol (60 dias).

Toda constelação é criada por ``satsim.constellation.sso_constellation`` (também a de 1 satélite).
Gera grade_nos.png, cobertura_latitude.png, revisita_mapa.png, contatos_24h.png,
manutencao_traco.png, manutencao_ltan.png, relatorio.md e resultados.json
(números brutos, para --from-json).
"""

import argparse
import json
import logging
import math
import sys
import time
from datetime import UTC, datetime, timedelta
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # sem janela: só gera os arquivos

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from satsim.analysis import fit_rate, osculating_elements  # noqa: E402
from satsim.astrotime import datetime_to_jd, gmst_rate  # noqa: E402
from satsim.constants import OMEGA_EARTH, OMEGA_SUN, R_EARTH, SECONDS_PER_DAY  # noqa: E402
from satsim.constellation import (  # noqa: E402
    CONSTELLATION_PALETTE,
    Constellation,
    sso_constellation,
)
from satsim.coverage import (  # noqa: E402
    TargetGrid,
    access_times,
    coverage_stats,
    latitude_profile,
    roll_for_half_grid,
    swath_half_angle,
)
from satsim.elements import KeplerianElements  # noqa: E402
from satsim.ephemeris import sun_ra_dec  # noqa: E402
from satsim.forces import (  # noqa: E402
    DragAcceleration,
    ThirdBodyAcceleration,
    ZonalAcceleration,
)
from satsim.groundtrack import grid_metrics, node_crossings, repeat_closure  # noqa: E402
from satsim.maintenance import (  # noqa: E402
    MINUTES_PER_RAD_OF_RAAN,
    YEAR_S,
    ltan_offset,
    maintenance_status,
    measure_decay_rate,
    time_to_deadband,
    time_to_ltan_band,
    track_drift,
)
from satsim.mean_osc import window_mean_elements  # noqa: E402
from satsim.observer import SAO_BENTO_DO_SAPUCAI  # noqa: E402
from satsim.orbit_design import (  # noqa: E402
    constellation_phases,
    frozen_eccentricity,
    ltan_drift_min_per_year,
    raan_for_ltan,
    repeat_orbit,
    repeat_semi_major_axis_calibrated,
    sso_inclination,
)
from satsim.passes import contact_plan, sun_elevation_deg  # noqa: E402
from satsim.propagators import MeanJ2Propagator  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_demo import full_force_models  # noqa: E402

logger = logging.getLogger("constellation_report")

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "reports" / "s5_constelacao"

EPOCH = datetime(2026, 10, 23, tzinfo=UTC)  # época padrão do Apêndice A
DAY = SECONDS_PER_DAY
N_REVS, N_DAYS = 91, 6
A_DESIGN = 6888089.0  # Sem4 §7.5 (A.1)
I_CAL = math.radians(97.46325)  # A.7 da Semana 4 (A.1)
LTAN_H = 22.5
GRID_DEG = 360.0 / N_REVS
PALETTE_P7 = (*CONSTELLATION_PALETTE, "#1f3a93")  # 7ª cor só para o caso P = 7 do gráfico


def _deg(x: float) -> float:
    return math.degrees(x)


def _metrics_deg(lon: np.ndarray) -> dict:
    count, gap_min, gap_max, offset = grid_metrics(lon, N_REVS)
    return {"nodes": count, "gap_min_deg": _deg(gap_min), "gap_max_deg": _deg(gap_max),
            "max_offset_deg": _deg(offset)}  # fmt: skip


def _nodes_by_sat(sats, start: datetime, duration_s: float, step_s: float = 30.0) -> list:
    return [node_crossings(s, start, duration_s, step_s)[1] for s in sats]


# --- S5-02: grade e repetição ---------------------------------------------------------------------


def grid_section() -> dict:
    """Grade no modelo médio e no Cowell completo, P = 1 / 6 / 7."""
    out: dict = {}
    # Modelo médio (Sem4 §9.2): MeanJ2 na órbita de projeto, e = 0; janela de −1 s (nó em t = 0)
    a, i = repeat_orbit(N_REVS, N_DAYS)
    for n_sats, days in ((1, 6.0), (6, 1.0), (7, 1.0)):
        sats = [MeanJ2Propagator(KeplerianElements(a, 0.0, i, 0.3, 0.0, float(p)), EPOCH)
                for p in constellation_phases(n_sats)]  # fmt: skip
        lon = np.concatenate(_nodes_by_sat(sats, EPOCH - timedelta(seconds=1), days * DAY))
        out[f"mean_p{n_sats}_{days:g}d"] = _metrics_deg(lon)

    # Cowell completo (A.1): P = 6 e P = 7 em 1 dia; P = 1 (S1) em 6 dias, 91 primeiros nós
    models = full_force_models(EPOCH)
    for n_sats in (6, 7):
        c = sso_constellation(EPOCH, models, n_sats=n_sats, a=A_DESIGN, i=I_CAL,
                              color_palette=PALETTE_P7)  # fmt: skip
        lon_by_sat = _nodes_by_sat(c, EPOCH, DAY)
        out[f"cowell_p{n_sats}_1d"] = _metrics_deg(np.concatenate(lon_by_sat))
        out[f"cowell_p{n_sats}_1d"]["lon_deg_by_sat"] = [np.degrees(x).tolist() for x in lon_by_sat]
        out[f"cowell_p{n_sats}_1d"]["colors"] = [s.color for s in c]
        if n_sats == 6:
            s1 = c[0]
    lon_s1 = node_crossings(s1, EPOCH, 6 * DAY + 3000.0)[1][:N_REVS]
    out["cowell_p1_6d"] = _metrics_deg(lon_s1)
    out["gmst_minus_omega_earth_rad_s"] = float(gmst_rate(datetime_to_jd(EPOCH))) - OMEGA_EARTH
    return out


def _s1(a: float, models) -> object:
    """S1 da constelação (ψ0 = 0) com semieixo a, pela fábrica (n_sats = 1)."""
    return sso_constellation(EPOCH, models, n_sats=1, a=a, i=I_CAL)[0]


def _closure_and_rates(sat) -> dict:
    t_nodes, lon = node_crossings(sat, EPOCH, 6 * DAY + 2 * 5697.0)
    t = np.arange(0.0, 6 * DAY + 1.0, 60.0)
    raan_dot = fit_rate(t, osculating_elements(sat.propagator, t).raan)
    t_nod = float(np.polyfit(np.arange(N_REVS + 1), t_nodes[: N_REVS + 1], 1)[0])
    t_day = 2.0 * math.pi / (OMEGA_EARTH - raan_dot)
    closure = repeat_closure(lon, N_REVS)
    return {"closure_deg": _deg(closure), "closure_km": closure * R_EARTH / 1e3,
            "raan_dot_deg_day": _deg(raan_dot) * DAY, "tnod_s": t_nod, "t_day_nod_s": t_day,
            "residual_s": N_REVS * t_nod - N_DAYS * t_day}  # fmt: skip


def repeat_section() -> dict:
    """Semieixo calibrado (J2–J4) e o fechamento do ciclo com e sem arrasto."""
    zonal = (ZonalAcceleration(),)
    e = frozen_eccentricity(A_DESIGN, float(sso_inclination(A_DESIGN)))
    raan = float(raan_for_ltan(LTAN_H, datetime_to_jd(EPOCH)))
    t0 = time.perf_counter()
    a_cal = repeat_semi_major_axis_calibrated(N_REVS, N_DAYS, e, I_CAL, EPOCH, zonal, raan=raan)
    out: dict = {"a_design_m": A_DESIGN, "a_calibrated_m": a_cal, "delta_a_m": a_cal - A_DESIGN,
                 "calibration_s": time.perf_counter() - t0}  # fmt: skip
    full = full_force_models(EPOCH)
    for label, models in (("j2j4", zonal), ("full", full)):
        for which, a in (("design", A_DESIGN), ("calibrated", a_cal)):
            out[f"{label}_{which}"] = _closure_and_rates(_s1(a, models))
    # Sem4 eq. 32: Ω̇ ∝ a^(−7/2) → deriva de LTAN causada por Δa
    raan_dot = OMEGA_SUN * (1.0 - 3.5 * out["delta_a_m"] / A_DESIGN)
    out["ltan_drift_delta_a_min_year"] = float(ltan_drift_min_per_year(raan_dot))
    d_raan = out["j2j4_calibrated"]["raan_dot_deg_day"] - out["j2j4_design"]["raan_dot_deg_day"]
    out["ltan_drift_measured_min_year"] = d_raan * 365.2422 / 15.0 * 60.0
    return out


# --- S5-03: cobertura e revisita ----------------------------------------------------------------

ROLL_DEG = 23.2
SUN_MIN_DEG = 10.0


def coverage_section() -> dict:
    """Cobertura em 1 dia (P = 1, 6), P = 1 em 6 dias e revisita de P = 6 em 3 dias (grade 1°)."""
    roll = math.radians(ROLL_DEG)
    grid = TargetGrid(1.0, 1.0)
    c = sso_constellation(EPOCH, full_force_models(EPOCH), a=A_DESIGN, i=I_CAL)
    s1 = Constellation(EPOCH, [c[0]])
    lat60 = np.abs(grid.lat_deg) <= 60.0
    out: dict = {
        "roll_deg": ROLL_DEG,
        "swath_half_angle_deg": _deg(swath_half_angle(A_DESIGN, roll)),
        "roll_half_grid_deg": _deg(roll_for_half_grid(A_DESIGN, N_REVS)),
        "band_lat_deg": grid.band_lat_deg.tolist(),
        "n_lon": int(grid.band_lon_deg.size),
        "profile_1d": {},
    }
    for label, const in (("p1", s1), ("p6", c)):
        for light, sun in (("dark", None), ("sun", SUN_MIN_DEG)):
            acc = access_times(const, grid, EPOCH, DAY, roll, min_sun_elevation_deg=sun)
            frac, revisit = coverage_stats(acc, DAY)
            _, prof = latitude_profile(grid, revisit < DAY)
            out["profile_1d"][f"{label}_{light}"] = {"fraction": frac, "by_lat": prof.tolist()}
    frac6, _ = coverage_stats(access_times(s1, grid, EPOCH, 6 * DAY, roll), 6 * DAY)
    out["fraction_p1_6d"] = frac6
    for light, sun in (("dark", None), ("sun", SUN_MIN_DEG)):
        acc = access_times(c, grid, EPOCH, 3 * DAY, roll, min_sun_elevation_deg=sun)
        frac, revisit = coverage_stats(acc, 3 * DAY)
        _, blackout = coverage_stats(acc, 3 * DAY, interval="blackout")
        seen = revisit < 3 * DAY
        out[f"revisit_3d_{light}"] = {
            "fraction": frac,
            "lat60_h": float(revisit[lat60].max() / 3600),
            "lat80_h": float(revisit[seen].max() / 3600) if seen.any() else None,
            "blackout_lat60_h": float(blackout[lat60].max() / 3600),
            "blackout_lat80_h": float(blackout[seen].max() / 3600) if seen.any() else None,
            "unseen_min_lat_deg": float(grid.lat_deg[~seen].min()) if (~seen).any() else None,
            # mapa em minutos inteiros (−1 = nunca visto), para o JSON ficar pequeno
            "map_min": np.where(seen, np.round(revisit / 60), -1).astype(int).tolist(),
        }
    _, dec = sun_ra_dec(datetime_to_jd(EPOCH + timedelta(days=1.5)))
    out["sun_dec_deg"] = _deg(float(dec))
    return out


def plot_coverage(data: dict, out_lat: Path, out_map: Path) -> None:
    """Fração coberta por latitude (1 dia) e mapa da revisita de P = 6 em 3 dias."""
    cov = data["coverage"]
    lat = np.array(cov["band_lat_deg"])
    fig, ax = plt.subplots(figsize=(8, 4.2))
    styles = {
        "p6_dark": ("P = 6", "#0b6e8a", "-"),
        "p6_sun": ("P = 6, só com luz", "#0b6e8a", "--"),
        "p1_dark": ("P = 1", "#ff8a3d", "-"),
        "p1_sun": ("P = 1, só com luz", "#ff8a3d", "--"),
    }
    for key, (label, color, ls) in styles.items():
        p = cov["profile_1d"][key]
        ax.plot(lat, 100 * np.array(p["by_lat"]), color=color, ls=ls, lw=1.8,
                label=f"{label} ({_f(100 * p['fraction'], 1)}%)")  # fmt: skip
    ax.set_xlabel("latitude do alvo [°]")
    ax.set_ylabel("alvos vistos em 1 dia [%]")
    ax.set_xlim(-80, 80)
    ax.set_ylim(-3, 103)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="lower center")
    ax.set_title(f"Cobertura em 1 dia, rolagem ±{_f(ROLL_DEG, 1)}° (2026-10-23, Cowell completo)",
                 fontsize=10)  # fmt: skip
    fig.tight_layout()
    fig.savefig(out_lat, dpi=130)
    plt.close(fig)

    n_lon = cov["n_lon"]
    lon_edges = np.arange(-180.0, 180.0 + 1e-9, 360.0 / n_lon) - 0.5 * 360.0 / n_lon
    lat_edges = np.r_[lat - 0.5, lat[-1] + 0.5]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4), sharey=True)
    # Escala própria por painel: sem luz a revisita vai de ≈ 12 h (equador) a ≈ 19 h (±80°); com
    # luz fica em ≈ 24 h quase em toda parte (uma visita por dia), mais longa só na borda da
    # noite polar.
    panels = (
        (axes[0], "revisit_3d_dark", "sem restrição de luz", (12.0, 19.0)),
        (axes[1], "revisit_3d_sun", f"só com Sol ≥ {SUN_MIN_DEG:g}° no alvo", (23.5, 26.0)),
    )
    for ax, key, title, (vmin, vmax) in panels:
        m = np.array(cov[key]["map_min"], dtype=float).reshape(lat.size, n_lon) / 60.0
        masked = np.ma.masked_less(m, 0.0)
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad("0.75")
        mesh = ax.pcolormesh(lon_edges, lat_edges, masked, cmap=cmap, vmin=vmin, vmax=vmax,
                             shading="flat")  # fmt: skip
        r = cov[key]
        ax.set_title(f"P = 6, 3 dias, {title}\nmáx. {_f(r['lat60_h'], 2)} h em |lat| ≤ 60°",
                     fontsize=9.5)  # fmt: skip
        ax.set_xlabel("longitude [°]")
        ax.set_aspect("equal")
        fig.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.2, shrink=0.8, extend="max",
                     label="intervalo máximo entre visitas [h] (cinza: nunca visto)")  # fmt: skip
    axes[0].set_ylabel("latitude [°]")
    fig.savefig(out_map, dpi=120, bbox_inches="tight")
    plt.close(fig)


def coverage_report_lines(data: dict) -> list[str]:
    cov = data["coverage"]
    prof = cov["profile_1d"]
    dark, sun = cov["revisit_3d_dark"], cov["revisit_3d_sun"]
    lat_lim = 90.0 - SUN_MIN_DEG + cov["sun_dec_deg"]
    swath_km = math.radians(cov["swath_half_angle_deg"]) * R_EARTH / 1e3
    return [
        "## S5-03 · Cobertura diária e revisita",
        "",
        f"Rolagem η = ±{_f(ROLL_DEG, 1)}°: alcance transversal λ = "
        f"{_f(cov['swath_half_angle_deg'], 4)}° = {_f(swath_km, 1)} km "
        f"(Sem4 eq. 45; A.5: 1,9784° = 220,2 km). A rolagem de meia grade (λ = π/91) é "
        f"{_f(cov['roll_half_grid_deg'], 4)}° (A.5: 23,1957°). Alvos em |lat| ≤ 80° numa grade de "
        "1° × 1°, passo de 10 s, constelação de referência no Cowell completo (`access_times`, "
        "`coverage_stats`); acessos a menos de 30 min contam como uma passagem.",
        "",
        "![Cobertura por latitude](cobertura_latitude.png)",
        "",
        "| Caso | Alvos vistos |",
        "|---|---:|",
        f"| P = 6, 1 dia | {_f(100 * prof['p6_dark']['fraction'], 2)}% (100%) |",
        f"| P = 1, 1 dia | {_f(100 * prof['p1_dark']['fraction'], 2)}% (48,48%) |",
        f"| P = 1, 6 dias | {_f(100 * cov['fraction_p1_6d'], 2)}% (100%) |",
        f"| P = 6, 1 dia, só com luz | {_f(100 * prof['p6_sun']['fraction'], 2)}% |",
        f"| P = 1, 1 dia, só com luz | {_f(100 * prof['p1_sun']['fraction'], 2)}% |",
        f"| P = 6, 3 dias, só com luz | {_f(100 * sun['fraction'], 2)}% (92,55%) |",
        "",
        "- **A promessa da constelação:** com seis satélites e ±23,2° de rolagem, todo ponto entre "
        "±80° é alcançado no mesmo dia; um satélite sozinho vê metade da Terra por dia e precisa "
        "do ciclo inteiro de 6 dias. A curva de P = 1 sobe perto de ±80°, onde as trilhas "
        "convergem e as faixas se sobrepõem.",
        "",
        "### Revisita em 3 dias",
        "",
        "![Mapa do intervalo máximo de revisita](revisita_mapa.png)",
        "",
        "| P = 6, 3 dias | Revisita máx. \\|lat\\| ≤ 60° | Revisita máx. \\|lat\\| ≤ 80° "
        "| Tempo sem acesso máx. \\|lat\\| ≤ 60° |",
        "|---|---:|---:|---:|",
        f"| sem restrição de luz | {_f(dark['lat60_h'], 2)} h (13,97 h) | "
        f"{_f(dark['lat80_h'], 2)} h (19,04 h) | {_f(dark['blackout_lat60_h'], 2)} h |",
        f"| só com Sol ≥ {SUN_MIN_DEG:g}° | {_f(sun['lat60_h'], 2)} h (24,26 h) | "
        f"{_f(sun['lat80_h'], 2)} h (dos vistos) | {_f(sun['blackout_lat60_h'], 2)} h |",
        "",
        '- **Revisita × tempo sem acesso:** a revisita (`interval="revisit"`, a do A.5) vai do '
        'início de uma passagem ao início da seguinte; o tempo sem acesso (`"blackout"`) vai do '
        "fim de uma ao início da outra. A diferença é a duração da passagem: perto do equador, "
        "dois satélites vizinhos do trem (15,8 min) podem ver o mesmo alvo; em latitude alta, as "
        "faixas convergem e vários satélites, em voltas seguidas, o veem em sequência.",
        "- **Sem luz: ≈ 12 a 14 h.** Cada alvo é visto duas vezes por dia, numa passagem "
        "ascendente (nó às 22h30 local) e numa descendente (10h30), a ≈ 12 h uma da outra. No "
        "equador a revisita máxima é ≈ 12,1 h; ela cresce com a latitude (≈ 14 h a ±60°, ≈ 19 h a "
        "±80°) porque fora do equador as duas passagens deixam de ser simétricas: a hora local "
        "do sobrevoo se afasta de 10h30 e 22h30 em sentidos opostos, e o intervalo entre elas "
        "fica alternadamente mais curto e mais longo que 12 h.",
        "- **Com luz: ≈ 24 h.** Uma câmera óptica só usa a passagem descendente, de manhã (10h30): "
        "a ascendente é noturna. Resta uma visita por dia, e a revisita dobra. Em outubro a "
        f"declinação do Sol é {_f(cov['sun_dec_deg'], 1)}°: ao meio-dia solar, ele fica abaixo "
        f"de {SUN_MIN_DEG:g}° de elevação para latitudes acima de 90° − {SUN_MIN_DEG:g}° + δ ≈ "
        f"{_f(lat_lim, 1)}° N. Nenhum alvo acima de ≈ {_f(sun['unseen_min_lat_deg'], 0)}° N é "
        f"visto com luz: a fração cai para {_f(100 * sun['fraction'], 2)}% (A.5: 92,55%). No "
        "hemisfério sul é o contrário: a primavera austral deixa o Sol alto, e com P = 1 a "
        "fração com luz sobe para perto de −80° (gráfico por latitude).",
        "",
    ]


# --- S5-04: plano de contatos ---------------------------------------------------------------------


def contacts_section() -> dict:
    """Plano de contatos de 24 h sobre São Bento do Sapucaí (el ≥ 10°) e o Sol na estação."""
    c = sso_constellation(EPOCH, full_force_models(EPOCH), a=A_DESIGN, i=I_CAL)
    plan = contact_plan(c, SAO_BENTO_DO_SAPUCAI, EPOCH, DAY, el_min_deg=10.0)

    def h(dt: datetime) -> float:
        return (dt - EPOCH).total_seconds() / 3600.0

    contacts = [
        {
            "sat": ct.pass_.sat_id,
            "aos_h": h(ct.pass_.aos),
            "los_h": h(ct.pass_.los),
            "tca_h": h(ct.pass_.tca),
            "max_el_deg": ct.pass_.max_elevation_deg,
            "duration_s": ct.pass_.duration_s,
            "sun_el_deg": ct.sun_elevation_deg,
            "lighting": ct.lighting,
            "train": ct.train,
        }  # fmt: skip
        for ct in plan.contacts
    ]
    trains = []
    for tr in plan.trains:
        aos = [contacts[k]["aos_h"] for k in tr.contact_indices]
        trains.append({"start_h": h(tr.start), "end_h": h(tr.end), "lighting": tr.lighting,
                       "n": len(tr.contact_indices),
                       "aos_spacing_min": (np.diff(aos) * 60).tolist()})  # fmt: skip
    t_sun = [EPOCH + timedelta(minutes=10 * k) for k in range(145)]
    return {
        "contacts": contacts,
        "trains": trains,
        "colors": {s.id: s.color for s in c},
        "total_contact_min": plan.total_contact_s / 60,
        "max_gap_h": plan.max_gap_s / 3600,
        "n_overlaps": len(plan.overlaps),
        "per_sat": plan.passes_per_satellite,
        "sun_h": [h(t) for t in t_sun],
        "sun_el_deg": sun_elevation_deg(SAO_BENTO_DO_SAPUCAI, t_sun).tolist(),
    }


def plot_contacts(data: dict, out: Path) -> None:
    """Linha do tempo de 24 h dos contatos (um satélite por linha) e um zoom em cada trem."""
    ct = data["contacts"]
    sats = sorted(ct["colors"])
    row = {s: k for k, s in enumerate(sats)}
    sun_h, sun_el = np.array(ct["sun_h"]), np.array(ct["sun_el_deg"])
    fig = plt.figure(figsize=(12, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.15), hspace=0.45, wspace=0.12)
    top = fig.add_subplot(gs[0, :])
    zooms = [fig.add_subplot(gs[1, k]) for k in range(len(ct["trains"][:2]))]

    def draw(ax, xlim, label_el: bool) -> None:
        night = sun_el < -6.0
        ax.fill_between(sun_h, -0.6, len(sats) - 0.4, where=night, color="#1d2330", alpha=0.12,
                        step="mid", lw=0)  # fmt: skip
        for c in ct["contacts"]:
            y = row[c["sat"]]
            ax.barh(y, c["los_h"] - c["aos_h"], left=c["aos_h"], height=0.6,
                    color=ct["colors"][c["sat"]], edgecolor="0.25", lw=0.5)  # fmt: skip
            if label_el and xlim[0] <= c["tca_h"] <= xlim[1]:
                ax.text(c["tca_h"], y + 0.42, f"{c['max_el_deg']:.0f}°", ha="center", va="top",
                        fontsize=7, color="0.25", clip_on=True)  # fmt: skip
        ax.set_yticks(range(len(sats)), [s.upper() for s in sats])
        ax.set_ylim(-0.6, len(sats) - 0.4)
        ax.invert_yaxis()
        ax.set_xlim(*xlim)
        ax.grid(axis="x", alpha=0.3)

    draw(top, (0.0, 24.0), False)
    ticks = np.arange(0, 25, 3)
    top.set_xticks(ticks, [f"{t:02d}h\n{(t - 3) % 24:02d}h" for t in ticks])
    top.set_xlabel("UTC (em cima) / BRT (embaixo) — faixas escuras: noite na estação (Sol < −6°)")
    top.set_title(
        f"Contatos com São Bento do Sapucaí em 24 h (el ≥ 10°): {len(ct['contacts'])} passagens, "
        f"{_f(ct['total_contact_min'], 1)} min, {ct['n_overlaps']} sobreposições",
        fontsize=10,
    )
    for ax, tr in zip(zooms, ct["trains"], strict=False):
        draw(ax, (tr["start_h"] - 0.1, tr["end_h"] + 0.1), True)
        hh = np.arange(math.floor(tr["start_h"] * 4) / 4, tr["end_h"] + 0.25, 0.25)
        ax.set_xticks(hh, [f"{int(x) % 24:02d}:{int(round(x % 1 * 60)):02d}" for x in hh],
                      fontsize=8)  # fmt: skip
        spacing = tr["aos_spacing_min"]
        ax.set_title(
            f"trem {'noturno' if tr['lighting'] == 'noite' else 'diurno'}: {tr['n']} passagens, "
            f"AOS a {_f(min(spacing), 1)}–{_f(max(spacing), 1)} min",
            fontsize=9.5,
        )
        ax.set_xlabel("UTC (rótulos: elevação máxima)")
    fig.savefig(out, dpi=125, bbox_inches="tight")
    plt.close(fig)


def contacts_report_lines(data: dict) -> list[str]:
    ct = data["contacts"]
    night, day = ct["trains"][0], ct["trains"][1]

    def brt(hours: float) -> str:
        m = int((hours - 3.0) % 24 * 60)  # trunca no minuto, como o A.6
        return f"{m // 60:02d}h{m % 60:02d}"

    spacing = night["aos_spacing_min"] + day["aos_spacing_min"]
    per_sat = ", ".join(f"{k.upper()} {v}" for k, v in sorted(ct["per_sat"].items()))
    rows = []
    for k, c in enumerate(ct["contacts"], start=1):
        rows.append(
            f"| {k} | {c['train'] + 1} | {c['sat'].upper()} | {brt(c['aos_h'])} | "
            f"{_f(c['max_el_deg'], 1)}° | {_f(c['duration_s'] / 60, 1)} min | "
            f"{_f(c['sun_el_deg'], 1)}° | {c['lighting']} |"
        )
    return [
        "## S5-04 · Plano de contatos sobre São Bento do Sapucaí",
        "",
        "`contact_plan` (estação −22,69°, −45,73°, 900 m; el ≥ 10°; 2026-10-23 00:00 UTC + 24 h; "
        "constelação de referência no Cowell completo). Dia/noite pela elevação do Sol na estação "
        "no TCA (crepúsculo civil, −6°); trens = passagens com AOS a menos de 30 min da anterior.",
        "",
        "![Linha do tempo dos contatos](contatos_24h.png)",
        "",
        f"- **{len(ct['contacts'])} passagens** ({per_sat}), "
        f"**{_f(ct['total_contact_min'], 1)} min** de contato, **{ct['n_overlaps']} "
        "sobreposições** (uma antena basta) — A.6: 16, 96,68 min, nenhuma.",
        f"- **Dois trens:** noturno de {brt(night['start_h'])} a {brt(night['end_h'])} BRT e "
        f"diurno de {brt(day['start_h'])} a {brt(day['end_h'])} BRT (A.6: 21h32–23h25 e "
        f"9h05–10h58). O maior intervalo sem contato é de {_f(ct['max_gap_h'], 2)} h.",
        "",
        "| # | Trem | Satélite | AOS (BRT) | El. máx. | Duração | Sol no TCA | Luz |",
        "|---:|---:|---|---|---:|---:|---:|---|",
        *rows,
        "",
        "### Por que as passagens chegam em trens de 15,8 min",
        "",
        "Os seis satélites estão no **mesmo plano**, em fila, separados de 60° em argumento "
        "de latitude: cada um cruza qualquer ponto da órbita Tnod/6 = 5 696,7/6 s ≈ 15,8 min "
        "depois do anterior (Sem4 §8.5, Tabela 4). Nesse intervalo a Terra gira ≈ 3,96° (um "
        "passo da grade), então o satélite seguinte passa ≈ 440 km a oeste no equador (≈ 400 km "
        "na latitude de São Bento): ainda dentro do círculo de visibilidade da estação, que tem "
        "≈ 1 590 km de raio para el ≥ 10° a 510 km (`footprint_half_angle`: 14,2°). Por isso a "
        "estação vê o trem passar inteiro, um satélite a cada ≈ 15,8 min (AOS a "
        f"{_f(min(spacing), 1)}–{_f(max(spacing), 1)} min: o tempo entre AOS varia porque cada "
        "satélite corta o círculo numa corda diferente, as primeiras do trem perto de uma borda, "
        "as do meio perto do centro e as últimas perto da borda oposta). As passagens do meio do "
        "trem são as mais altas e longas; as das pontas, rasantes e curtas.",
        "",
        "O plano orbital acompanha o Sol (SSO): o trem passa pela estação **duas vezes por "
        "dia**, na metade ascendente da órbita (nó às 22h30 local, de noite) e na descendente "
        "(10h30, de manhã), **sempre nos mesmos horários solares**. Entre os dois trens não há "
        "nenhum contato: as outras ≈ 13 voltas de cada satélite passam longe demais de São Bento.",
        "",
    ]


# --- S5-05: manutenção ----------------------------------------------------------------------------

BALLISTIC = (0.0044, 0.011, 0.022)  # C_D·A/m [m²/kg]: 0,4×, 1× (A.1) e 2× o de referência
DEADBAND_KM = 5.0
LTAN_BAND_MIN = 5.0


def _lon_diff_km(lon_a: np.ndarray, lon_b: np.ndarray) -> np.ndarray:
    d = np.mod(lon_a - lon_b + np.pi, 2.0 * np.pi) - np.pi
    return d * R_EARTH / 1e3


def maintenance_section() -> dict:
    """ȧ medido, MaintenanceStatus, deriva do traço no Cowell × fórmula e validação do LTAN."""
    full = full_force_models(EPOCH)
    out: dict = {}
    # 1. Constelação de referência: ȧ por satélite (1 dia) e o indicador
    c = sso_constellation(EPOCH, full, a=A_DESIGN, i=I_CAL)
    out["status"] = [
        maintenance_status(s.id, measure_decay_rate(s, 1.0), A_DESIGN, I_CAL).as_dict() for s in c
    ]
    # 2. Deriva do traço no Cowell completo: nós com arrasto − nós sem arrasto (mesmo a calibrado)
    a_cal = data_a_calibrated()
    days = 14.0
    with_drag = sso_constellation(EPOCH, full, n_sats=1, a=a_cal, i=I_CAL)[0]
    no_drag = sso_constellation(EPOCH, full[:-1], n_sats=1, a=a_cal, i=I_CAL)[0]
    t_nodes, lon_drag = node_crossings(with_drag, EPOCH, days * DAY)
    _, lon_free = node_crossings(no_drag, EPOCH, days * DAY)
    n = min(lon_drag.size, lon_free.size)
    a_dot_cal = measure_decay_rate(with_drag, 6.0)
    _, closure_lon = node_crossings(with_drag, EPOCH, 6 * DAY + 2 * 5697.0)
    out["cowell_track"] = {
        "a_calibrated_m": a_cal,
        "t_days": (t_nodes[:n] / DAY).tolist(),
        "drift_km": _lon_diff_km(lon_drag[:n], lon_free[:n]).tolist(),
        "a_dot_m_day": a_dot_cal * DAY,
        "closure_6d_km": repeat_closure(closure_lon, N_REVS) * R_EARTH / 1e3,
        "formula_6d_km": float(track_drift(6 * DAY, a_dot_cal, a_cal)[1]),
    }
    # 3. Três coeficientes balísticos: ȧ no Cowell J2–J4 + arrasto (2 dias)
    out["ballistic"] = []
    for cda in BALLISTIC:
        models = (ZonalAcceleration(), DragAcceleration(cd=2.2, area_over_mass=cda / 2.2))
        sat = sso_constellation(EPOCH, models, n_sats=1, a=a_cal, i=I_CAL)[0]
        a_dot = measure_decay_rate(sat, 2.0)
        t_band = time_to_deadband(DEADBAND_KM, a_dot, a_cal) / DAY
        out["ballistic"].append({"cd_a_m": cda, "a_dot_m_day": a_dot * DAY, "days_to_band": t_band})
    # 4. LTAN no Cowell J2–J4 + Sol (e + Lua), 60 dias: di/dt, aceleração do LTAN × eq. 32
    out["ltan"] = {}
    for label, models in (
        ("sun", (ZonalAcceleration(), ThirdBodyAcceleration("sun", EPOCH))),
        ("sun_moon", (ZonalAcceleration(), ThirdBodyAcceleration("sun", EPOCH),
                      ThirdBodyAcceleration("moon", EPOCH))),
    ):  # fmt: skip
        out["ltan"][label] = _ltan_validation(models)
    return out


def data_a_calibrated() -> float:
    """O semieixo calibrado do S5-02 (do JSON, se já calculado; senão, recalcula)."""
    path = OUT_DIR / "resultados.json"
    if path.exists():
        saved = json.loads(path.read_text(encoding="utf-8"))
        if "repeat" in saved:
            return float(saved["repeat"]["a_calibrated_m"])
    return float(repeat_section()["a_calibrated_m"])


def _ltan_validation(models, days: float = 60.0) -> dict:
    """i(t) e Ω(t) médios por volta no Cowell: di/dt e a parte quadrática do LTAN × eq. 32."""
    sat = sso_constellation(EPOCH, models, n_sats=1, a=A_DESIGN, i=I_CAL)[0]
    t_nodes, _ = node_crossings(sat, EPOCH, days * DAY)
    period = float(np.mean(np.diff(t_nodes)))
    starts = t_nodes[:-1]
    inc = window_mean_elements(sat.propagator, starts, period, 30.0).i
    tt = (starts[:, None] + np.linspace(0.0, period, 60)[None, :]).ravel()
    raan = np.unwrap(osculating_elements(sat.propagator, tt).raan).reshape(starts.size, -1)
    tc = starts + period / 2.0
    ltan_min = (raan.mean(axis=1) - OMEGA_SUN * tc) * MINUTES_PER_RAD_OF_RAAN
    di_dt = float(np.polyfit(tc, inc, 1)[0])
    c2, c1, c0 = np.polyfit(tc, ltan_min, 2)
    return {
        "di_dt_deg_day": math.degrees(di_dt) * DAY,
        "ltan_accel_min_year2": 2.0 * c2 * YEAR_S**2,
        "ltan_accel_formula_min_year2": 2.0 * float(ltan_offset(YEAR_S, di_dt, I_CAL)),
        "t_days": (tc / DAY).tolist(),
        "ltan_min": ltan_min.tolist(),  # LTAN − nominal medido [min], com tendência linear
        "formula_min": np.asarray(ltan_offset(tc, di_dt, I_CAL)).tolist(),  # eq. 32, di/dt medido
    }


def plot_maintenance(data: dict, out_track: Path, out_ltan: Path) -> None:
    """Deriva do traço × tempo (3 coeficientes balísticos + Cowell) e deriva do LTAN em 2 anos."""
    m = data["maintenance"]
    a_cal = m["cowell_track"]["a_calibrated_m"]
    t = np.linspace(0.0, 20.0, 400)
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    colors = ("#7cff6b", "#0b6e8a", "#ff4d8d")
    for b, color in zip(m["ballistic"], colors, strict=True):
        _, km = track_drift(t * DAY, b["a_dot_m_day"] / DAY, a_cal)
        ax.plot(t, km, color=color, lw=2,
                label=f"C_D·A/m = {_f(b['cd_a_m'], 4)} m²/kg: ȧ = {_f(b['a_dot_m_day'], 1)} m/dia, "
                      f"±5 km em {_f(b['days_to_band'], 1)} dias")  # fmt: skip
    ct = m["cowell_track"]
    ax.plot(ct["t_days"], ct["drift_km"], ".", ms=2.5, color="0.15",
            label="Cowell completo (C_D·A/m = 0,011): nós com arrasto − sem arrasto")  # fmt: skip
    ax.axhline(DEADBAND_KM, color="0.4", ls="--", lw=1)
    ax.text(0.3, DEADBAND_KM + 0.4, "banda morta ±5 km", fontsize=8.5, color="0.3")
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 20)
    ax.set_xlabel("tempo desde a última manobra [dias]")
    ax.set_ylabel("desvio do traço no equador [km]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.8, loc="upper left")
    ax.set_title(
        "Deriva do traço por arrasto: λ(t) = 3π|ȧ|t²/(2a) (Sem4 §7.6) × Cowell", fontsize=10
    )
    fig.tight_layout()
    fig.savefig(out_track, dpi=130)
    plt.close(fig)

    st = m["status"][0]
    di_dt = math.radians(st["di_dt_deg_day"]) / DAY
    years = np.linspace(0.0, 2.0, 300)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3))
    ax = axes[0]
    label = f"ΔLTAN = ½·k·(di/dt)·t², di/dt = {_sci(st['di_dt_deg_day'], 1)} °/d (S4-02)"
    ax.plot(years, ltan_offset(years * YEAR_S, di_dt, I_CAL), color="#0b6e8a", lw=2, label=label)
    t_band = time_to_ltan_band(LTAN_BAND_MIN, di_dt, I_CAL) / YEAR_S
    ax.axhline(-LTAN_BAND_MIN, color="0.4", ls="--", lw=1)
    ax.axvline(t_band, color="0.6", ls=":", lw=1)
    ax.text(t_band + 0.03, -1.0, f"±5 min em {_f(t_band, 2)} ano", fontsize=8.5, color="0.3")
    ax.set_xlabel("tempo [anos]")
    ax.set_ylabel("afastamento do LTAN [min]")
    ax.set_xlim(0, 2)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    ax.set_title("Deriva do LTAN em 2 anos (Sem4 eq. 32)", fontsize=10)
    ax = axes[1]
    for key, color, label in (("sun", "#ff8a3d", "Cowell J2–J4 + Sol"),
                              ("sun_moon", "#8f7bff", "Cowell J2–J4 + Sol + Lua")):  # fmt: skip
        v = m["ltan"][key]
        t_d = np.array(v["t_days"])

        def detrend(y: list, t_d: np.ndarray = t_d) -> np.ndarray:
            """Tira a reta de melhor ajuste: sobra a parte quadrática (a mesma base nas duas)."""
            y = np.asarray(y)
            return y - np.polyval(np.polyfit(t_d, y, 1), t_d)

        ax.plot(t_d, detrend(v["ltan_min"]) * 60, ".", ms=2, color=color,
                label=f"{label}: {_f(v['ltan_accel_min_year2'], 2)} min/ano²")  # fmt: skip
        ax.plot(t_d, detrend(v["formula_min"]) * 60, color=color, lw=1.5,
                label=f"eq. 32 com o di/dt medido ({_sci(v['di_dt_deg_day'])} °/d): "
                      f"{_f(v['ltan_accel_formula_min_year2'], 2)}")  # fmt: skip
    ax.set_xlabel("tempo [dias]")
    ax.set_ylabel("parte quadrática do LTAN [s]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.2, loc="lower center")
    ax.set_title("Validação no Cowell (60 dias; tendência linear removida)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_ltan, dpi=125)
    plt.close(fig)


def maintenance_report_lines(data: dict) -> list[str]:
    m = data["maintenance"]
    ct = m["cowell_track"]
    st0 = m["status"][0]
    sun, sun_moon = m["ltan"]["sun"], m["ltan"]["sun_moon"]
    design_full = data.get("repeat", {}).get("full_design", {}).get("closure_km")
    rows = [
        f"| {s['sat_id'].upper()} | {_f(s['a_dot_m_day'], 2)} | "
        f"{_f(s['days_to_track_correction'], 2)} | {_f(s['track_delta_v_per_maneuver_m_s'], 4)} | "
        f"{_f(s['track_delta_v_per_year_m_s'], 2)} | {_f(s['ltan_offset_1y_min'], 2)} | "
        f"{_f(s['days_to_ltan_correction'] / 365.2422, 2)} | "
        f"{_f(s['inclination_delta_v_per_year_m_s'], 2)} |"
        for s in m["status"]
    ]
    err = ct["formula_6d_km"] / ct["closure_6d_km"] - 1.0
    lines = [
        "## S5-05 · Manutenção: deriva do traço e do LTAN",
        "",
        "Indicador de manutenção (`maintenance_status`), com ȧ medido no Cowell completo por "
        "médias em janelas nodais (`measure_decay_rate`, 1 dia) e di/dt = "
        f"{_sci(st0['di_dt_deg_day'], 1)} °/dia (Sol, S4-02). Bandas: traço ±{DEADBAND_KM:g} km no "
        f"equador, LTAN ±{LTAN_BAND_MIN:g} min. A.7: ȧ −23,71 m/dia → 6,95 dias, 0,0910 m/s por "
        "manobra, 4,78 m/s/ano; LTAN −3,43 min em 1 ano, 1,21 ano, 4,75 m/s/ano.",
        "",
        "| Satélite | ȧ [m/dia] | Próxima correção do traço [dias] | Δv por manobra [m/s] "
        "| Δv/ano traço [m/s] | LTAN em 1 ano [min] | LTAN fora da banda [anos] "
        "| Δv/ano inclinação [m/s] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
        *rows,
        "",
        "- **ȧ medido:** −24,0 m/dia nos seis satélites (A.7: −23,71, 1,0% de diferença). As "
        "janelas de média começam nos nós ascendentes medidos e duram o período nodal medido: com "
        "o período calculado do semieixo osculador (≈ 12 s curto) o resto do termo de 2u escorrega "
        "de janela em janela e o ȧ de 6 dias sairia −30 m/dia.",
        "",
        "### Traço: a fórmula do Sem4 contra o Cowell",
        "",
        "![Deriva do traço](manutencao_traco.png)",
        "",
        "- **Curva do Cowell:** diferença de longitude, nó a nó, entre o S1 com e sem arrasto "
        f"(mesmo a calibrado, {_f(ct['a_calibrated_m'], 1)} m, Cowell completo nos dois). Ela "
        "segue a parábola da fórmula com o ȧ medido: a deriva do traço é quadrática, como no "
        "Sem4 §7.6.",
        f"- **Fechamento de 6 dias com arrasto:** Cowell {_f(ct['closure_6d_km'], 3)} km × fórmula "
        f"{_f(ct['formula_6d_km'], 3)} km com ȧ = {_f(ct['a_dot_m_day'], 2)} m/dia "
        f"({_f(100 * err, 1)}%; critério ≤ 10%; A.4/A.7: 3,68 × 3,72 km). Com o a de projeto o "
        + (
            f"Cowell dá {_f(design_full, 2)} km: "
            if design_full is not None
            else "Cowell dá mais: "
        )
        + "≈ 1 km de viés de semieixo (S5-02) somado à deriva do arrasto.",
        "- **Três coeficientes balísticos:** o tempo até a banda cai com √|ȧ| (Sem4 eq. 41); "
        "dobrar o C_D·A/m não reduz o intervalo à metade, mas a 1/√2.",
        "",
        "### LTAN: a inclinação que o Sol puxa",
        "",
        "![Deriva do LTAN](manutencao_ltan.png)",
        "",
        "- **Fórmula:** δΩ̇ = −Ω☉·tan i·δi (Sem4 eq. 32) com δi(t) = (di/dt)·t: o LTAN se afasta "
        f"quadraticamente, {_f(st0['ltan_offset_1y_min'], 2)} min em 1 ano, e sai da banda de "
        f"±5 min em {_f(st0['days_to_ltan_correction'] / 365.2422, 2)} ano.",
        "- **Validação no Cowell J2–J4 + Sol (60 dias):** di/dt = "
        f"{_sci(sun['di_dt_deg_day'])} °/dia (S4-02: −9,8e-5) e aceleração do LTAN "
        f"{_f(sun['ltan_accel_min_year2'], 2)} min/ano² × "
        f"{_f(sun['ltan_accel_formula_min_year2'], 2)} "
        "pela eq. 32 com o di/dt medido. Com a Lua, os termos de ≈ 14 dias contaminam o ajuste em "
        f"60 dias: di/dt = {_sci(sun_moon['di_dt_deg_day'])} °/dia, "
        f"{_f(sun_moon['ltan_accel_min_year2'], 2)} × "
        f"{_f(sun_moon['ltan_accel_formula_min_year2'], 2)} "
        "min/ano². O di/dt do S4-02 é a parte do Sol; ele também varia com a estação, então o "
        "indicador é uma previsão de ordem de grandeza para um ano.",
        "",
        "### Dois orçamentos parecidos, por motivos diferentes",
        "",
        f"O traço consome ≈ {_f(st0['track_delta_v_per_year_m_s'], 1)} m/s por ano em ≈ "
        f"{_f(365.2422 / st0['days_to_track_correction'], 0)} manobras pequenas (uma por semana, "
        "no plano, repondo o semieixo que o arrasto tira); o LTAN, ≈ "
        f"{_f(st0['inclination_delta_v_per_year_m_s'], 1)} m/s por ano numa correção de "
        "inclinação (fora do plano). Os números são parecidos por coincidência da órbita de "
        "510 km; o que dita a operação é o traço, pela cadência (Sem4 §6.2). As manobras não são "
        "simuladas: o indicador prevê quando e quanto.",
        "",
    ]
    return lines


# --- Gráficos e relatório ------------------------------------------------------------------------


def plot_grid(data: dict, out: Path) -> None:
    """Figura 12 do Sem4 no Cowell completo: nós de 1 dia para P = 6 e P = 7 (0° a 120°)."""
    fig, axes = plt.subplots(2, 1, figsize=(11, 4.6), sharex=True)
    for ax, key, title in (
        (axes[0], "cowell_p6_1d", "P = 6: uma passagem por linha de grade"),
        (axes[1], "cowell_p7_1d", "P = 7 na mesma órbita: passagens fora da grade"),
    ):
        case = data["grid"][key]
        anchor = case["lon_deg_by_sat"][0][0]
        for k in range(-N_REVS, N_REVS):
            x = np.mod(anchor + k * GRID_DEG, 360.0)
            ax.axvline(x, color="0.82", lw=0.8, zorder=0)
        for j, (lon, color) in enumerate(zip(case["lon_deg_by_sat"], case["colors"], strict=True)):
            ax.scatter(np.mod(lon, 360.0), np.zeros(len(lon)), s=60, color=color,
                       edgecolor="0.2", lw=0.5, label=f"S{j + 1}", zorder=3)  # fmt: skip
        ax.set_title(
            f"{title} — {case['nodes']} nós/dia, espaçamento {_f(case['gap_min_deg'], 3)}° a "
            f"{_f(case['gap_max_deg'], 3)}°, afastamento máx. {_f(case['max_offset_deg'], 4)}°",
            fontsize=9.5,
        )
        ax.set_yticks([])
        ax.set_xlim(0.0, 120.0)
    axes[1].legend(ncol=7, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.32))
    axes[1].set_xlabel("longitude do nó ascendente [°] (linhas cinza: grade de 91 trilhas)")
    fig.suptitle("Nós ascendentes em 1 dia, Cowell completo (2026-10-23)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def _sci(x: float, nd: int = 2) -> str:
    """Notação científica no formato do relatório: −9,82e-5 (vírgula, sinal "−", expoente curto)."""
    mant, exp = f"{x:.{nd}e}".split("e")
    return f"{mant.replace('.', ',').replace('-', '−')}e{int(exp)}"


def _f(x: float, nd: int) -> str:
    """Número no formato do relatório: espaço nos milhares, vírgula decimal, sinal "−"."""
    return f"{x:,.{nd}f}".replace(",", " ").replace(".", ",").replace("-", "−")


def write_report(data: dict, out: Path, elapsed_s: float | None) -> None:
    g, r = data["grid"], data["repeat"]
    gen = f" em {_f(elapsed_s / 60, 1)} min" if elapsed_s is not None else ""
    lines = [
        "# Relatório da Semana 5 — Constelação 91/6",
        "",
        f"Gerado por `python scripts/constellation_report.py`{gen}. Época 2026-10-23 00:00 UTC. "
        "Constelação de referência do Apêndice A.1 do prompt da Semana 5 (`sso_constellation`: a "
        "6 888,089 km, e congelada, i 97,46325°, ω 90°, LTAN 22h30, ψ0 = 0°, 300°, …, 60°). "
        "Cowell completo = J2–J4 + Sol/Lua + SRP (C_R 1,3, A/m 0,01 m²/kg) + arrasto "
        "(C_D·A/m 0,011 m²/kg, geodésico). Valores do Apêndice entre parênteses.",
        "",
        "## S5-02 · Grade de 91 trilhas e repetição",
        "",
        "![Nós ascendentes em 1 dia](grade_nos.png)",
        "",
        "### Grade no equador",
        "",
        "Nós por `node_crossings` (interpolação linear da latitude, passo 30 s); métricas por "
        "`grid_metrics` (espaçamentos com a volta de 360°; afastamento à grade ideal de "
        f"360°/91 = {_f(GRID_DEG, 5)}°, ancorada no 1º nó de S1).",
        "",
        "| Caso | Nós | Espaçamento mín. … máx. [°] | Afastamento máx. [°] |",
        "|---|---:|---|---:|",
    ]
    rows = (
        ("mean_p1_6d", "Modelo médio, P = 1, 6 dias", "(91; 3,9560 uniforme)", ""),
        ("mean_p6_1d", "Modelo médio, P = 6, 1 dia", "(91; 3,9560 uniforme)", ""),
        ("mean_p7_1d", "Modelo médio, P = 7, 1 dia", "(107; 0,5651 … 3,3909)", ""),
        ("cowell_p6_1d", "Cowell completo, P = 6, 1 dia", "(91; 3,95173 … 3,96431)", "(0,0064)"),
        ("cowell_p1_6d", "Cowell completo, P = 1, 6 dias (91 primeiros)", "(3,9169 … 3,9674)", ""),
        ("cowell_p7_1d", "Cowell completo, P = 7, 1 dia", "(106; 3,3834 … 3,9585)", "(1,694)"),
    )
    for key, label, ref, ref_off in rows:
        c = g[key]
        lines.append(
            f"| {label} | {c['nodes']} | {_f(c['gap_min_deg'], 5)} … {_f(c['gap_max_deg'], 5)} "
            f"{ref} | {_f(c['max_offset_deg'], 4)} {ref_off} |"
        )
    drift6 = _deg(g["gmst_minus_omega_earth_rad_s"]) * 6 * DAY
    gmst_txt = f"{g['gmst_minus_omega_earth_rad_s']:.3e}".replace(".", ",")
    lines += [
        "",
        "- **Teorema P = D (Sem4 §8.3) no Cowell completo:** seis satélites em um dia depositam as "
        "91 trilhas, com espaçamentos a ±0,006° do ideal; sete satélites na mesma órbita deixam a "
        "grade (afastamento de 1,7°, 6/7 de passo), com quase-duplicatas e lacunas mantidas.",
        "- **Modelo médio:** com 1 satélite em 6 dias o passo varia "
        f"{_f(g['mean_p1_6d']['gap_max_deg'] - g['mean_p1_6d']['gap_min_deg'], 6)}°. Não é a "
        "teoria: o ECEF gira com o GMST real, "
        f"{gmst_txt} rad/s mais rápido que o ω⊕ = 7,292115e-5 rad/s "
        f"do Sem4, o que desloca a grade em {_f(drift6, 6)}° (≈ 28 m) em 6 dias.",
        "- **Afastamento no Cowell (P = 6):** os desvios à grade têm padrão de período 3 (S1/S4, "
        "S2/S5, S3/S6, ±0,0045°): é o termo de curto período de J2 em Ω e u, que a conversão "
        "médio → osculador (S4-04) não corrige (ela ajusta a, o vetor e e i). Ancorado no 1º nó no "
        "tempo (o de S3), o afastamento máximo seria 0,0107°.",
        "- **Precisão do nó:** contra nós refinados por `brentq` em z_ECI = 0 (S1, 15 nós em "
        "24 h), a interpolação linear erra 7,3e-5° / 1,1e-3 s com passo de 60 s, "
        "9,3e-6° / 2,7e-4 s com 30 s e 5,2e-7° com 10 s (`tests/test_repeat_grid.py`; A.3: "
        "7,3e-5°, 9,3e-6°, 4,8e-7°). O refinamento por `brentq` não é necessário.",
        "",
        "### Repetição de 6 dias e semieixo calibrado",
        "",
        "`repeat_semi_major_axis_calibrated` (Cowell J2–J4, conversão só zonal): secante em a "
        "até N·Tnod = D·Tdia,nod, Tnod pelo ajuste linear dos nós e Tdia,nod = 2π/(ω⊕ − Ω̇) com Ω̇ "
        f"por `fit_rate`. Resultado: **a = {_f(r['a_calibrated_m'], 2)} m, Δa = "
        f"+{_f(r['delta_a_m'], 2)} m** (A.4: 6 888 108,37 m, +19,37 m), em "
        f"{_f(r['calibration_s'], 1)} s.",
        "",
        "| Forças | a | Fechamento após 91 voltas | N·Tnod − D·Tdia,nod [s] | Tnod [s] | Ω̇ [°/d] |",
        "|---|---|---|---:|---:|---:|",
    ]
    refs = {
        "j2j4_design": "(+0,00895° = 1,00 km)",
        "j2j4_calibrated": "(−0,00020° = −0,022 km)",
        "full_design": "(+0,0422° = 4,69 km)",
        "full_calibrated": "(+0,0330° = 3,68 km)",
    }
    for key, ref in refs.items():
        c = r[key]
        forces = "J2–J4" if key.startswith("j2j4") else "completo (com arrasto)"
        a_label = "projeto" if key.endswith("design") else "calibrado"
        lines.append(
            f"| {forces} | {a_label} | {_f(c['closure_deg'], 5)}° = {_f(c['closure_km'], 3)} km "
            f"{ref} | {_f(c['residual_s'], 4)} | {_f(c['tnod_s'], 4)} | "
            f"{_f(c['raan_dot_deg_day'], 6)} |"
        )
    lines += [
        "",
        "- **Fechamento com o a calibrado (J2–J4):** o resíduo de "
        f"{_f(r['j2j4_calibrated']['closure_km'] * 1e3, 0)} m é, na maior parte, a mesma "
        "diferença ω⊕ × GMST real do modelo médio (−28 m previstos em 6 dias): a calibração usa o "
        "ω⊕ do Sem4 em Tdia,nod, como pede o card.",
        "- **Com arrasto o ciclo não fecha** mesmo com o a calibrado (≈ 3,7 km por ciclo): é a "
        "deriva quadrática do Sem4 §7.6, assunto do S5-05. A calibração é feita sem arrasto, "
        "como a da inclinação.",
        "- **Efeito no LTAN:** Δa = +19 m muda Ω̇ por Ω̇ ∝ a^(−7/2) (Sem4 eq. 32): "
        f"{_f(r['ltan_drift_delta_a_min_year'], 4)} min/ano (Sem4: −0,014). Medido no Cowell "
        "J2–J4 (diferença dos Ω̇ ajustados em 6 dias): "
        f"{_f(r['ltan_drift_measured_min_year'], 4)} min/ano. Desprezível diante da banda de "
        "±5 min.",
        "",
    ]
    if "coverage" in data:
        lines += coverage_report_lines(data)
    if "contacts" in data:
        lines += contacts_report_lines(data)
    if "maintenance" in data:
        lines += maintenance_report_lines(data)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


SECTIONS = {
    "grid": ("S5-02: grade", grid_section),
    "repeat": ("S5-02: repetição", repeat_section),
    "coverage": ("S5-03: cobertura", coverage_section),
    "contacts": ("S5-04: contatos", contacts_section),
    "maintenance": ("S5-05: manutenção", maintenance_section),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--from-json", action="store_true", help="refaz só gráficos e relatório")
    parser.add_argument("--only", action="append", choices=sorted(SECTIONS),
                        help="recalcula só esta seção; as outras vêm do JSON")  # fmt: skip
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUT_DIR / "resultados.json"
    elapsed = None
    data = json.loads(json_path.read_text(encoding="utf-8")) if json_path.exists() else {}
    if not args.from_json:
        t0 = time.perf_counter()
        for key in args.only or SECTIONS:
            label, compute = SECTIONS[key]
            logger.info("%s …", label)
            data[key] = compute()
        elapsed = time.perf_counter() - t0 if not args.only else None
        json_path.write_text(json.dumps(data), encoding="utf-8")
    plot_grid(data, OUT_DIR / "grade_nos.png")
    if "coverage" in data:
        plot_coverage(data, OUT_DIR / "cobertura_latitude.png", OUT_DIR / "revisita_mapa.png")
    if "contacts" in data:
        plot_contacts(data, OUT_DIR / "contatos_24h.png")
    if "maintenance" in data:
        plot_maintenance(data, OUT_DIR / "manutencao_traco.png", OUT_DIR / "manutencao_ltan.png")
    write_report(data, OUT_DIR / "relatorio.md", elapsed)
    logger.info("relatório em %s", OUT_DIR)


if __name__ == "__main__":
    main()
