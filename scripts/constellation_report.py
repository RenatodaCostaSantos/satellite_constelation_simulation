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

Toda constelação é criada por ``satsim.constellation.sso_constellation`` (também a de 1 satélite).
Gera grade_nos.png, cobertura_latitude.png, revisita_mapa.png, relatorio.md e resultados.json
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
from satsim.forces import ZonalAcceleration  # noqa: E402
from satsim.groundtrack import grid_metrics, node_crossings, repeat_closure  # noqa: E402
from satsim.orbit_design import (  # noqa: E402
    constellation_phases,
    frozen_eccentricity,
    ltan_drift_min_per_year,
    raan_for_ltan,
    repeat_orbit,
    repeat_semi_major_axis_calibrated,
    sso_inclination,
)
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
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


SECTIONS = {
    "grid": ("S5-02: grade", grid_section),
    "repeat": ("S5-02: repetição", repeat_section),
    "coverage": ("S5-03: cobertura", coverage_section),
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
    write_report(data, OUT_DIR / "relatorio.md", elapsed)
    logger.info("relatório em %s", OUT_DIR)


if __name__ == "__main__":
    main()
