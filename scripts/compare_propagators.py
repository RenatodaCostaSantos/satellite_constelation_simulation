"""Compara os propagadores Kepler, MeanJ2 e Cowell e gera o relatório da Semana 2.

Uso (na raiz do repositório, com o pacote instalado com o extra dev):

    python scripts/compare_propagators.py

Gera em reports/s2_comparacao/ cinco figuras (PNG) e relatorio.md com a tabela de números e as
conclusões. Órbitas, por 10 dias a 60 s:

- SSO de projeto (a = 6888,089 km, e = 1,074e-3, i = 97,4396°, Ω0 = 30°, ω0 = 90°);
- ISS-like (a = R⊕ + 420 km, e = 5e-4, i = 51,64°, Ω0 = 30°, ω0 = 0°);
- órbita de teste e = 0,02 (a = 7200 km, i = 45°, Ω0 = 30°, ω0 = 30°), para ω(t).

A lição central é a diferença entre elementos médios e osculadores: os elementos nominais,
usados como condição inicial (osculadora) do Cowell, correspondem a uma órbita média diferente.
"""

import argparse
import logging
import math
import time
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # sem janela: só gera os arquivos

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from satsim.analysis import fit_rate, mean_elements, osculating_elements  # noqa: E402
from satsim.constants import OMEGA_SUN, R_EARTH  # noqa: E402
from satsim.elements import KeplerianElements  # noqa: E402
from satsim.propagators import CowellPropagator, KeplerPropagator, MeanJ2Propagator  # noqa: E402
from satsim.secular import nodal_period, secular_rates  # noqa: E402

logger = logging.getLogger("compare_propagators")

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "reports" / "s2_comparacao"

EPOCH = datetime(2026, 9, 29, 12, 0, 0, tzinfo=UTC)
DAY = 86400.0
DURATION_S = 10 * DAY
STEP_S = 60.0
DEG_DAY = math.degrees(1.0) * DAY  # rad/s → °/dia

# Faixas do Apêndice A.4 / critérios do card S2-06
RAAN_DOT_MAX_REL_ERR = 5e-3  # Ω̇ ajustado vs teoria no semieixo médio
ERR_B_MAX_KM = 20.0  # MeanJ2 com elementos médios medidos vs Cowell, 24 h
ERR_A_MIN_KM = 500.0  # MeanJ2 nominal vs Cowell, 24 h

COLORS = {"kepler": "#8d99ae", "mean_j2": "#ef476f", "cowell": "#118ab2", "mean": "#ffb703"}


def _elements(
    a_km: float, e: float, i_deg: float, raan_deg: float, argp_deg: float
) -> KeplerianElements:
    return KeplerianElements(
        a=a_km * 1e3,
        e=e,
        i=math.radians(i_deg),
        raan=math.radians(raan_deg),
        argp=math.radians(argp_deg),
        psi=0.0,
    )


ORBITS = {
    "SSO": ("SSO de projeto", _elements(6888.089, 1.074e-3, 97.4396, 30.0, 90.0)),
    "ISS-like": ("ISS-like", _elements(R_EARTH / 1e3 + 420.0, 5e-4, 51.64, 30.0, 0.0)),
    "e = 0,02": ("órbita de teste e = 0,02", _elements(7200.0, 0.02, 45.0, 30.0, 30.0)),
}


@dataclass
class Run:
    """Propagações e medidas de uma órbita."""

    key: str
    label: str
    nominal: KeplerianElements
    t: np.ndarray
    osc: KeplerianElements  # elementos osculadores do Cowell
    a_mean: float
    e_mean: float
    i_mean: float
    raan_dot_fit: float
    argp_dot_fit: float
    raan_dot_theory_mean: float
    raan_dot_theory_nominal: float
    argp_dot_theory_mean: float
    raan_mean_j2: np.ndarray  # Ω(t) do MeanJ2 nominal
    err_a: np.ndarray  # |MeanJ2(nominal) − Cowell| [m]
    err_b: np.ndarray  # |MeanJ2(médios medidos) − Cowell| [m]
    err_kepler: np.ndarray  # |Kepler − Cowell| [m]


def run_orbit(key: str) -> Run:
    """Propaga a órbita com os três propagadores e mede elementos médios e taxas."""
    label, nominal = ORBITS[key]
    t = np.arange(0.0, DURATION_S + STEP_S / 2, STEP_S)
    cowell = CowellPropagator(nominal, EPOCH)
    osc = osculating_elements(cowell, t)
    period = nodal_period(nominal.a, nominal.e, nominal.i)
    a_mean, e_mean, i_mean = mean_elements(t, osc, period)
    th_mean = secular_rates(a_mean, e_mean, i_mean)
    th_nom = secular_rates(nominal.a, nominal.e, nominal.i)

    mean_j2_nominal = MeanJ2Propagator(nominal, EPOCH)
    mean_j2_measured = MeanJ2Propagator(replace(nominal, a=a_mean, e=e_mean, i=i_mean), EPOCH)
    r_cowell, _ = cowell.propagate(t)
    r_nom, _ = mean_j2_nominal.propagate(t)
    r_meas, _ = mean_j2_measured.propagate(t)
    r_kep, _ = KeplerPropagator(nominal, EPOCH).propagate(t)

    return Run(
        key=key,
        label=label,
        nominal=nominal,
        t=t,
        osc=osc,
        a_mean=a_mean,
        e_mean=e_mean,
        i_mean=i_mean,
        raan_dot_fit=fit_rate(t, osc.raan),
        argp_dot_fit=fit_rate(t, osc.argp),
        raan_dot_theory_mean=th_mean.raan_dot,
        raan_dot_theory_nominal=th_nom.raan_dot,
        argp_dot_theory_mean=th_mean.argp_dot,
        raan_mean_j2=mean_j2_nominal.elements_at(t).raan,
        err_a=np.linalg.norm(r_nom - r_cowell, axis=1),
        err_b=np.linalg.norm(r_meas - r_cowell, axis=1),
        err_kepler=np.linalg.norm(r_kep - r_cowell, axis=1),
    )


def _rel_err(fit: float, theory: float) -> float:
    return fit / theory - 1.0


def _max_24h_km(run: Run, err: np.ndarray) -> float:
    return float(np.max(err[run.t <= DAY])) / 1e3


# --- Figuras -----------------------------------------------------------------------------------


def fig_raan(runs: dict[str, Run], out: Path) -> None:
    """Figura 1: Ω(t) − Ω0 para Kepler (constante), MeanJ2 (reta) e Cowell (osculador)."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=110)
    for ax, key in zip(axes, ("SSO", "ISS-like"), strict=True):
        r = runs[key]
        days = r.t / DAY
        raan0 = r.nominal.raan
        ax.plot(
            days, np.zeros_like(days), color=COLORS["kepler"], lw=2, label="Kepler (Ω constante)"
        )
        ax.plot(days, np.degrees(np.unwrap(r.osc.raan) - raan0), color=COLORS["cowell"], lw=1.0,
                label="Cowell (osculador)")  # fmt: skip
        ax.plot(days, np.degrees(r.raan_mean_j2 - raan0), "--", color=COLORS["mean_j2"], lw=1.5,
                label="MeanJ2 (reta, Sem4 eq. 23)")  # fmt: skip
        ax.set_title(
            f"{r.label}: dΩ/dt ajustado = {r.raan_dot_fit * DEG_DAY:+.4f} °/dia", fontsize=10
        )
        ax.set_xlabel("tempo [dias]")
        ax.set_ylabel("Ω − Ω0 [°]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Figura 1 — regressão do nó: J2 gira o plano da órbita", fontsize=11)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_argp(run: Run, sso: Run, out: Path) -> None:
    """Figura 2: ω(t) na órbita e = 0,02, com a nota sobre e ≈ 1e-3."""
    fig, ax = plt.subplots(figsize=(8, 4.4), dpi=110)
    days = run.t / DAY
    argp0 = run.nominal.argp
    theory = run.argp_dot_theory_mean * run.t
    ax.plot(days, np.degrees(np.unwrap(run.osc.argp) - argp0), color=COLORS["cowell"], lw=1.0,
            label="Cowell (ω osculador)")  # fmt: skip
    ax.plot(days, np.degrees(theory), "--", color=COLORS["mean_j2"], lw=1.5,
            label="teoria (Sem4 eq. 24) no semieixo médio")  # fmt: skip
    ax.set_xlabel("tempo [dias]")
    ax.set_ylabel("ω − ω0 [°]")
    fit, th = run.argp_dot_fit * DEG_DAY, run.argp_dot_theory_mean * DEG_DAY
    ax.set_title(
        f"Figura 2 — rotação da linha dos apsides, {run.label}:\n"
        f"dω/dt ajustado {fit:+.4f} °/dia vs teoria {th:+.4f} °/dia",
        fontsize=10,
    )
    note = (
        "Por que não usar e ≈ 1e-3: o perigeu de uma órbita quase circular é mal definido.\n"
        "As oscilações de curto período de J2 no vetor excentricidade são da mesma ordem\n"
        f"que o próprio e; na SSO o ajuste de ω daria {sso.argp_dot_fit * DEG_DAY:+.2f} °/dia "
        f"(teoria {sso.argp_dot_theory_mean * DEG_DAY:+.2f})."
    )
    ax.text(0.02, 0.97, note, transform=ax.transAxes, va="top", fontsize=8,
            bbox={"boxstyle": "round", "fc": "white", "alpha": 0.9})  # fmt: skip
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_semi_major_axis(runs: dict[str, Run], out: Path) -> None:
    """Figura 3: semieixo osculador do Cowell, média sobre revoluções inteiras e nominal (12 h)."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=110)
    for ax, key in zip(axes, ("SSO", "ISS-like"), strict=True):
        r = runs[key]
        w = r.t <= 0.5 * DAY
        hours = r.t[w] / 3600.0
        ax.plot(
            hours, r.osc.a[w] / 1e3, color=COLORS["cowell"], lw=1.0, label="a osculador (Cowell)"
        )
        ax.axhline(r.a_mean / 1e3, color=COLORS["mean"], lw=2,
                   label=f"a médio (10 dias) = {r.a_mean / 1e3:.2f} km")  # fmt: skip
        ax.axhline(r.nominal.a / 1e3, color=COLORS["mean_j2"], ls="--", lw=1.5,
                   label=f"a nominal = {r.nominal.a / 1e3:.2f} km")  # fmt: skip
        p2p = (np.max(r.osc.a[w]) - np.min(r.osc.a[w])) / 1e3
        ax.set_title(f"{r.label}: oscilação pico a pico ≈ {p2p:.1f} km", fontsize=10)
        ax.set_xlabel("tempo [h]")
        ax.set_ylabel("semieixo maior [km]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="lower right")
    fig.suptitle(
        "Figura 3 — elementos osculadores oscilam a cada meia revolução em torno do valor médio",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_error_mean_j2(runs: dict[str, Run], out: Path) -> None:
    """Figura 4: |MeanJ2 − Cowell| em 24 h, (a) nominais como osculadores e (b) médios medidos."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=110, sharey=True)
    for ax, key in zip(axes, ("SSO", "ISS-like"), strict=True):
        r = runs[key]
        w = (r.t <= DAY) & (r.t > 0.0)
        hours = r.t[w] / 3600.0
        ax.semilogy(hours, r.err_a[w] / 1e3, color=COLORS["mean_j2"],
                    label="(a) elementos nominais tratados como osculadores")  # fmt: skip
        ax.semilogy(hours, r.err_b[w] / 1e3, color=COLORS["cowell"],
                    label="(b) elementos médios medidos na simulação")  # fmt: skip
        ax.set_title(
            f"{r.label}: máx 24 h (a) {_max_24h_km(r, r.err_a):.0f} km, "
            f"(b) {_max_24h_km(r, r.err_b):.1f} km",
            fontsize=10,
        )
        ax.set_xlabel("tempo [h]")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8, loc="lower right")
    axes[0].set_ylabel("|r_MeanJ2 − r_Cowell| [km]")
    fig.suptitle("Figura 4 — erro de posição MeanJ2 × Cowell: médio ≠ osculador", fontsize=11)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_error_kepler(runs: dict[str, Run], out: Path) -> None:
    """Figura 5: |Kepler − Cowell| em 10 dias."""
    fig, ax = plt.subplots(figsize=(8, 4.4), dpi=110)
    for key, color in (("SSO", COLORS["cowell"]), ("ISS-like", COLORS["mean_j2"])):
        r = runs[key]
        ax.plot(r.t / DAY, r.err_kepler / 1e3, color=color, lw=1.0, label=r.label)
    ax.set_xlabel("tempo [dias]")
    ax.set_ylabel("|r_Kepler − r_Cowell| [km]")
    ax.set_title("Figura 5 — ignorar J2: regressão do nó + deriva ao longo da órbita", fontsize=10)
    note = (
        "O Kepler mantém o plano fixo (o Cowell precessa) e anda com outro movimento\n"
        "médio (n nominal vs dψ/dt + dω/dt no semieixo médio). Na SSO domina a deriva ao\n"
        "longo da órbita, e o erro satura perto do diâmetro da órbita (≈ 13 800 km),\n"
        "quando os dois ficam em lados opostos."
    )
    ax.text(0.98, 0.04, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox={"boxstyle": "round", "fc": "white", "alpha": 0.9})  # fmt: skip
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


# --- Relatório ---------------------------------------------------------------------------------


def _ok(cond: bool) -> str:
    return "OK" if cond else "**FORA**"


def _n(value: float, fmt: str) -> str:
    """Número formatado com vírgula decimal (o relatório é em português)."""
    return format(value, fmt).replace(".", ",")


def _table_row(r: Run) -> str:
    err_mean = _rel_err(r.raan_dot_fit, r.raan_dot_theory_mean)
    err_nom = _rel_err(r.raan_dot_fit, r.raan_dot_theory_nominal)
    err_a, err_b = _max_24h_km(r, r.err_a), _max_24h_km(r, r.err_b)
    ok = abs(err_mean) < RAAN_DOT_MAX_REL_ERR and err_b < ERR_B_MAX_KM
    if r.key != "e = 0,02":  # o card pede (a) > 500 km para SSO e ISS-like
        ok = ok and err_a > ERR_A_MIN_KM
    cells = [
        r.label,
        _n((r.a_mean - r.nominal.a) / 1e3, "+.2f"),
        _n(r.raan_dot_fit * DEG_DAY, "+.4f"),
        f"{_n(r.raan_dot_theory_mean * DEG_DAY, '+.4f')} ({_n(100 * err_mean, '+.3f')}%)",
        f"{_n(r.raan_dot_theory_nominal * DEG_DAY, '+.4f')} ({_n(100 * err_nom, '+.3f')}%)",
        _n(err_a, ".0f"),
        _n(err_b, ".1f"),
        _n(float(np.max(r.err_kepler)) / 1e3, ".0f"),
        _ok(ok),
    ]
    return "| " + " | ".join(cells) + " |"


def _conclusions(runs: dict[str, Run]) -> str:
    """Parágrafo de conclusões (até 8 linhas), montado com os números calculados."""
    sso, iss, e002 = runs["SSO"], runs["ISS-like"], runs["e = 0,02"]
    worst = max(abs(_rel_err(r.raan_dot_fit, r.raan_dot_theory_mean)) for r in runs.values())
    return (
        f"J2 faz o plano da órbita precessar: Ω̇ = {_n(sso.raan_dot_fit * DEG_DAY, '+.3f')} °/dia "
        f"na SSO (a teoria nominal dá a taxa do Sol, {_n(OMEGA_SUN * DEG_DAY, '.4f')} °/dia; a "
        f"diferença de {_n(100 * _rel_err(sso.raan_dot_fit, sso.raan_dot_theory_nominal), '+.1f')}%"
        " vem do semieixo médio, ver abaixo) e "
        f"{_n(iss.raan_dot_fit * DEG_DAY, '+.3f')} °/dia na ISS-like; na órbita e = 0,02 ele "
        f"também gira a linha dos apsides (ω̇ = {_n(e002.argp_dot_fit * DEG_DAY, '+.3f')} °/dia). "
        f"Ignorar J2 (Kepler) leva a erros de {_n(float(np.max(sso.err_kepler)) / 1e3, '.0f')} km "
        f"(SSO) e {_n(float(np.max(iss.err_kepler)) / 1e3, '.0f')} km (ISS-like) em 10 dias. "
        f"O MeanJ2 reproduz a taxa medida no Cowell dentro de {_n(100 * worst, '.2f')}%, *desde "
        "que* a teoria use o semieixo médio medido: elementos médios e osculadores diferem "
        f"(a_médio − a_nominal = {_n((sso.a_mean - sso.nominal.a) / 1e3, '+.1f')} km na SSO, "
        f"{_n((iss.a_mean - iss.nominal.a) / 1e3, '+.1f')} km na ISS-like). Tratar os elementos "
        f"nominais como osculadores dá {_n(_max_24h_km(sso, sso.err_a), '.0f')} km de erro em 24 h "
        f"na SSO, contra {_n(_max_24h_km(sso, sso.err_b), '.1f')} km com os médios medidos; o "
        "resíduo são as oscilações de curto período (Figura 3), que o MeanJ2 não representa por "
        "construção."
    )


def write_report(runs: dict[str, Run], elapsed_s: float, out: Path) -> None:
    """Escreve relatorio.md com a tabela de números, as figuras e as conclusões."""
    e002 = runs["e = 0,02"]
    argp_err = _rel_err(e002.argp_dot_fit, e002.argp_dot_theory_mean)
    lines = [
        "# Relatório da Semana 2 — Kepler × MeanJ2 × Cowell",
        "",
        f"Gerado por `python scripts/compare_propagators.py` em {elapsed_s:.0f} s. Época "
        f"{EPOCH:%Y-%m-%d %H:%M} UTC; 10 dias a 60 s; Cowell DOP853 (rtol = 1e-11, atol = 1e-6), "
        "só J2. Ω0 = 30° e ψ0 = 0 em todas as órbitas; ω0 = 90° (SSO), 0° (ISS-like, ver nota) e "
        "30° (e = 0,02). Faixas de referência: Apêndice A.4 do prompt da Semana 2.",
        "",
        "## Números",
        "",
        "| Órbita | a_médio − a_nominal [km] | Ω̇ ajustado [°/dia] | Ω̇ teoria@médio (erro) "
        "| Ω̇ teoria@nominal (erro) | erro (a) 24 h [km] | erro (b) 24 h [km] "
        "| \\|Kepler − Cowell\\| 10 d [km] | faixas A.4 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|:---:|",
        *(_table_row(r) for r in runs.values()),
        "",
        f"**ω̇ na órbita e = 0,02:** ajustado {_n(e002.argp_dot_fit * DEG_DAY, '+.4f')} °/dia vs "
        f"teoria no semieixo médio {_n(e002.argp_dot_theory_mean * DEG_DAY, '+.4f')} °/dia "
        f"({_n(100 * argp_err, '+.3f')}%; critério < 1%: {_ok(abs(argp_err) < 0.01)}).",
        "",
        "Erro (a): MeanJ2 com os elementos nominais tratados como osculadores. Erro (b): MeanJ2 "
        "com a, e, i médios medidos no Cowell (ângulos nominais). Critérios: Ω̇ teoria@médio "
        "< 0,5%; (b) < 20 km; (a) > 500 km (SSO e ISS-like).",
        "",
        "**Nota (ISS-like, ω0 = 0°):** o semieixo osculador oscila como "
        "1,5·J2·R⊕²/a·sin²i·cos 2u (u = ω + ν), então o semieixo médio depende de onde a órbita "
        "começa. Com ω0 = 0° (início no nó) reproduz-se a linha ISS-like do Apêndice A.4; com "
        "ω0 = 90° o sinal de a_médio − a_nominal se inverte (ver o golden file, S2-02).",
        "",
        "## Figuras",
        "",
        "1. [Ω(t): Kepler, MeanJ2 e Cowell](fig1_raan.png)",
        "2. [ω(t) na órbita e = 0,02](fig2_argp_e002.png)",
        "3. [Semieixo osculador e médio](fig3_semieixo.png)",
        "4. [\\|MeanJ2 − Cowell\\| em 24 h, (a) e (b)](fig4_erro_meanj2.png)",
        "5. [\\|Kepler − Cowell\\| em 10 dias](fig5_erro_kepler.png)",
        "",
        "## Conclusões",
        "",
        _conclusions(runs),
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR, help="pasta de saída")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    start = time.perf_counter()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    runs = {key: run_orbit(key) for key in ORBITS}
    fig_raan(runs, args.out_dir / "fig1_raan.png")
    fig_argp(runs["e = 0,02"], runs["SSO"], args.out_dir / "fig2_argp_e002.png")
    fig_semi_major_axis(runs, args.out_dir / "fig3_semieixo.png")
    fig_error_mean_j2(runs, args.out_dir / "fig4_erro_meanj2.png")
    fig_error_kepler(runs, args.out_dir / "fig5_erro_kepler.png")
    elapsed = time.perf_counter() - start
    write_report(runs, elapsed, args.out_dir / "relatorio.md")
    logger.info("Relatório gerado em %s (%.1f s)", args.out_dir, time.perf_counter() - start)


if __name__ == "__main__":
    main()
