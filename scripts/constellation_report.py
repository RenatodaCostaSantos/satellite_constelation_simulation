"""Relatório da constelação 91/6 da Semana 5 (reports/s5_constelacao/).

Uso (na raiz do repositório, com o pacote instalado com o extra dev):

    python scripts/constellation_report.py               # calcula tudo (~1 min) e gera o relatório
    python scripts/constellation_report.py --from-json   # só refaz gráficos e relatório

Seções (uma por card; os próximos cards acrescentam as suas):

- **S5-02 · Grade e repetição.** Nós ascendentes da constelação de referência (Apêndice A.1 do
  prompt da Semana 5) no Cowell completo em 1 dia para P = 6 e P = 7 (equivalente à Figura 12 do
  Sem4), a mesma grade no modelo médio (Sem4 §9.2), a precisão do nó por interpolação linear, o
  fechamento do ciclo de 6 dias com o a de projeto e com o a calibrado
  (``repeat_semi_major_axis_calibrated``), com e sem arrasto, e o efeito de Δa no LTAN.

Toda constelação é criada por ``satsim.constellation.sso_constellation`` (também a de 1 satélite).
Gera grade_nos.png, relatorio.md e resultados.json (números brutos, para --from-json).
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
from satsim.constellation import CONSTELLATION_PALETTE, sso_constellation  # noqa: E402
from satsim.elements import KeplerianElements  # noqa: E402
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
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--from-json", action="store_true", help="refaz só gráficos e relatório")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUT_DIR / "resultados.json"
    elapsed = None
    if args.from_json:
        data = json.loads(json_path.read_text(encoding="utf-8"))
    else:
        t0 = time.perf_counter()
        logger.info("S5-02: grade …")
        data = {"grid": grid_section()}
        logger.info("S5-02: repetição …")
        data["repeat"] = repeat_section()
        elapsed = time.perf_counter() - t0
        json_path.write_text(json.dumps(data), encoding="utf-8")
    plot_grid(data, OUT_DIR / "grade_nos.png")
    write_report(data, OUT_DIR / "relatorio.md", elapsed)
    logger.info("relatório em %s", OUT_DIR)


if __name__ == "__main__":
    main()
