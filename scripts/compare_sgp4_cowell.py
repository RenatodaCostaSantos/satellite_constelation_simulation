"""Compara o Cowell (J2 e J2 + arrasto) com o SGP4 da ISS e gera o relatório do S3-03.

Uso (na raiz do repositório, com o pacote instalado com o extra dev):

    python scripts/compare_sgp4_cowell.py                 # TLE de teste (Apêndice A.1)
    python scripts/compare_sgp4_cowell.py --tle-file meu.tle

Gera em reports/s3_arrasto/ o gráfico erro_sgp4_cowell.png (erro de posição 3D vs tempo, 7 dias,
escala log), o gráfico erro_vs_coeficiente.png (erro em 24 h, 72 h e 7 dias em função do C_D·A/m
usado: um "V" com mínimo no valor ajustado) e relatorio.md com as tabelas, o C_D·A/m ajustado e
a explicação de por que o arrasto padrão piora o resultado.

Todas as configurações partem do estado do SGP4 em t = 0 (CowellPropagator.from_state):

(a) Cowell só com J2;
(b) J2 + arrasto com os parâmetros padrão (C_D = 2,2, A/m = 0,004 m²/kg), altitude esférica e
    geodésica;
(c) J2 + arrasto com C_D·A/m ajustado por mínimos quadrados ao SGP4 numa janela de 72 h.

O TLE padrão é o congelado de tests/data/iss_sample.tle (época 2019-12-09), para que os números
reproduzam o Apêndice A.5 do prompt da Semana 3.
"""

import argparse
import logging
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # sem janela: só gera os arquivos

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from satsim.analysis import (  # noqa: E402
    cowell_j2_drag,
    fit_ballistic_coefficient,
    position_error,
)
from satsim.propagators import SGP4Propagator  # noqa: E402
from satsim.tle import Tle  # noqa: E402

logger = logging.getLogger("compare_sgp4_cowell")

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "reports" / "s3_arrasto"
DEFAULT_TLE = ROOT / "tests" / "data" / "iss_sample.tle"

HOUR = 3600.0
DURATION_S = 7 * 24 * HOUR
STEP_S = 600.0
CHECK_H = (24.0, 72.0, 168.0)
FIT_WINDOW_S = 72 * HOUR

CD = 2.2
DEFAULT_AREA_OVER_MASS = 0.004  # ISS, Sem3 §2.4
MAX_FITTED_ERR_KM = (2.0, 8.0, 40.0)  # critério do card S3-03 (24 h / 72 h / 7 dias)
# Varredura de C_D·A/m [m²/kg]: curva do gráfico e linhas da tabela do "V" (mais o ajustado e o
# padrão, acrescentados em run)
SWEEP_CURVE = tuple(np.linspace(0.0, 0.009, 10))
SWEEP_TABLE = (0.0, 0.0008, 0.0010, 0.0012, 0.0014, 0.0018, 0.0020, 0.0025, 0.0035, 0.005)
A5_KM = {  # Apêndice A.5 do prompt da Semana 3 (referência)
    "j2": (1.0, 8.7, 47.0),
    "default_spherical": (6.6, 59.6, 325.0),
    "default_geodetic": (5.8, 52.4, 286.0),
    "fitted": (0.3, 2.7, 14.5),
}


@dataclass(frozen=True)
class Case:
    key: str
    label: str
    cd_a_m: float
    altitude_model: str
    color: str
    style: str


@dataclass(frozen=True)
class Result:
    case: Case
    err_m: np.ndarray  # erro em cada t da grade

    def at_hours(self, t: np.ndarray, hours: float) -> float:
        return float(self.err_m[np.argmin(np.abs(t - hours * HOUR))])


@dataclass(frozen=True)
class Sweep:
    cd_a_m: np.ndarray  # (M,) valores varridos, crescentes
    rms72_km: np.ndarray  # (M,) RMS do erro em [0, 72 h]: o que o ajuste minimiza
    rms7d_km: np.ndarray  # (M,) RMS do erro em [0, 7 dias]
    point_km: np.ndarray  # (M, 3) erro pontual em 24 h, 72 h e 7 dias

    def index(self, cd_a_m: float) -> int:
        return int(np.argmin(np.abs(self.cd_a_m - cd_a_m)))


def sweep(sgp4: SGP4Propagator, t: np.ndarray, values: list[float]) -> Sweep:
    """Erro do Cowell (J2 + arrasto geodésico) × SGP4 para cada C_D·A/m da lista."""
    vals = np.array(sorted(set(values)))
    in_fit = t <= FIT_WINDOW_S + 1e-9
    checks = [int(np.argmin(np.abs(t - h * HOUR))) for h in CHECK_H]
    rms72, rms7d, point = [], [], []
    for c in vals:
        err = position_error(cowell_j2_drag(sgp4, c, cd=CD), sgp4, t) / 1e3
        rms72.append(np.sqrt(np.mean(err[in_fit] ** 2)))
        rms7d.append(np.sqrt(np.mean(err**2)))
        point.append(err[checks])
    return Sweep(vals, np.array(rms72), np.array(rms7d), np.array(point))


def run(tle: Tle) -> tuple[np.ndarray, list[Result], float, float, Sweep]:
    sgp4 = SGP4Propagator(tle)
    t = np.arange(0.0, DURATION_S + 0.5 * STEP_S, STEP_S)

    start = time.perf_counter()
    fitted = fit_ballistic_coefficient(sgp4, window_s=FIT_WINDOW_S, step_s=STEP_S, cd=CD)
    fit_s = time.perf_counter() - start
    logger.info("C_D·A/m ajustado = %.5f m²/kg (%.1f s)", fitted, fit_s)

    default = CD * DEFAULT_AREA_OVER_MASS
    cases = [
        Case("j2", "(a) J2 apenas", 0.0, "geodetic", "#118ab2", "-"),
        Case("default_spherical", "(b) J2 + arrasto padrão, alt. esférica", default,
             "spherical", "#ef476f", ":"),
        Case("default_geodetic", "(b) J2 + arrasto padrão, alt. geodésica", default,
             "geodetic", "#ef476f", "-"),
        Case("fitted", f"(c) J2 + arrasto ajustado (C_D·A/m = {_dec(fitted, '.5f')})", fitted,
             "geodetic", "#06d6a0", "-"),
    ]  # fmt: skip
    results = []
    for case in cases:
        cowell = cowell_j2_drag(sgp4, case.cd_a_m, cd=CD, altitude_model=case.altitude_model)
        results.append(Result(case, position_error(cowell, sgp4, t)))
    sweep_result = sweep(sgp4, t, [*SWEEP_CURVE, *SWEEP_TABLE, fitted, default])
    return t, results, fitted, fit_s, sweep_result


def plot(t: np.ndarray, results: list[Result], tle: Tle, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.2), dpi=120)
    days = t[1:] / (24 * HOUR)  # t = 0 tem erro nulo: fora da escala log
    for res in results:
        ax.plot(days, res.err_m[1:] / 1e3, res.case.style, color=res.case.color, lw=1.6,
                label=res.case.label)  # fmt: skip
    for hours in CHECK_H:
        ax.axvline(hours / 24, color="#adb5bd", lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_xlim(0, DURATION_S / (24 * HOUR))
    ax.set_xlabel("tempo desde a época do TLE [dias]")
    ax.set_ylabel("erro de posição 3D |r_Cowell − r_SGP4| [km]")
    name = tle.name or tle.norad_id
    ax.set_title(f"Cowell × SGP4 — {name}, época {tle.epoch:%Y-%m-%d %H:%M} UTC")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_sweep(sw: Sweep, fitted: float, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.2), dpi=120)
    x = sw.cd_a_m * 1e3
    ax.plot(x, sw.rms72_km, "o-", color="#118ab2", ms=3.5, lw=1.8,
            label="RMS em 0–72 h (o que o ajuste minimiza)")  # fmt: skip
    ax.plot(x, sw.rms7d_km, "o-", color="#073b4c", ms=3.5, lw=1.8, label="RMS em 0–7 dias")
    ax.plot(x, sw.point_km[:, 2], "--", color="#adb5bd", lw=1.2,
            label="erro pontual em t = 7 dias (oscila com a órbita)")  # fmt: skip
    default = CD * DEFAULT_AREA_OVER_MASS
    ax.axvline(fitted * 1e3, color="#06d6a0", lw=1.2, ls="--")
    ax.axvline(default * 1e3, color="#adb5bd", lw=1.2, ls="--")
    ax.text(fitted * 1e3, 0.6, f"  ajustado (72 h)\n  {_dec(fitted * 1e3, '.2f')}",
            color="#06a77d", fontsize=9, va="bottom")  # fmt: skip
    ax.text(default * 1e3, 2.5, "padrão  \n8,80  ", color="#6c757d", fontsize=9, va="bottom",
            ha="right")  # fmt: skip
    ax.set_yscale("log")
    ax.set_ylim(0.5, None)
    ax.set_xlabel("C_D·A/m usado no Cowell [10⁻³ m²/kg]   (0 = só J2)")
    ax.set_ylabel("erro de posição 3D contra o SGP4 [km]")
    ax.set_title("O erro depende de quão errado está o arrasto (forma de V)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def _dec(x: float, fmt: str) -> str:
    """Número com vírgula decimal (relatório em português)."""
    return format(x, fmt).replace(".", ",")


def _km(x_m: float) -> str:
    return _dec(x_m / 1e3, ".1f")


def write_report(
    t: np.ndarray, results: list[Result], tle: Tle, fitted: float, fit_s: float, sw: Sweep,
    elapsed_s: float, out: Path,
) -> None:  # fmt: skip
    by_key = {r.case.key: r for r in results}
    fitted_err = [by_key["fitted"].at_hours(t, h) / 1e3 for h in CHECK_H]
    ok = all(e <= lim for e, lim in zip(fitted_err, MAX_FITTED_ERR_KM, strict=True))
    j2 = [by_key["j2"].at_hours(t, h) for h in CHECK_H]
    dflt = [by_key["default_geodetic"].at_hours(t, h) for h in CHECK_H]
    lines = [
        "# Relatório do S3-03 — Cowell (J2 + arrasto) × SGP4",
        "",
        f"Gerado por `python scripts/compare_sgp4_cowell.py` em {elapsed_s:.0f} s (ajuste: "
        f"{fit_s:.0f} s). TLE: {tle.name or '—'} (catálogo {tle.norad_id}), época "
        f"{tle.epoch:%Y-%m-%d %H:%M:%S} UTC. Todas as configurações partem do estado do SGP4 em "
        "t = 0; Cowell DOP853 (rtol = 1e-10, atol = 1e-3); atmosfera exponencial por partes de "
        "Vallado (`satsim.forces.atmosphere`), co-rotante. Referência: Apêndice A.5 do prompt da "
        "Semana 3.",
        "",
        "## Erro de posição 3D contra o SGP4",
        "",
        "| Configuração | C_D·A/m [m²/kg] | 24 h [km] | 72 h [km] | 7 dias [km] | A.5 (24 h / 72 h / 7 d) |",  # noqa: E501
        "|---|---:|---:|---:|---:|---|",
    ]
    for res in results:
        errs = " | ".join(_km(res.at_hours(t, h)) for h in CHECK_H)
        ref = " / ".join(_dec(x, "g") for x in A5_KM[res.case.key])
        lines.append(f"| {res.case.label} | {_dec(res.case.cd_a_m, '.5f')} | {errs} | {ref} |")
    lines += [
        "",
        f"**C_D·A/m ajustado:** {_dec(fitted, '.5f')} m²/kg (janela de 72 h, amostras a cada "
        "600 s, mínimos quadrados no RMS do erro de posição). Com C_D = 2,2 equivale a "
        f"A/m = {_dec(fitted / CD, '.2g')} m²/kg. Referência A.5: 0,00163; faixa aceita: "
        "0,0005 a 0,005.",
        "",
        "**Critério do card (ajustado ≤ 2 / 8 / 40 km):** "
        f"{' / '.join(_dec(e, '.1f') for e in fitted_err)} km → {'OK' if ok else 'FALHA'}.",
        "",
        "![Erro de posição Cowell × SGP4](erro_sgp4_cowell.png)",
        "",
        "## O erro em função do C_D·A/m: um “V”",
        "",
        "Seria natural esperar que acrescentar arrasto sempre aproximasse o Cowell do SGP4. A "
        "varredura abaixo mostra que não: o Cowell (J2 + arrasto, altitude geodésica) foi rodado "
        "com vários valores de C_D·A/m, do zero (só J2) ao padrão. As duas primeiras colunas são "
        "o RMS do erro numa janela (variam suavemente); as três últimas, o erro num instante "
        "exato (oscilam com a órbita, como no gráfico acima).",
        "",
        "| C_D·A/m usado [m²/kg] | RMS 0–72 h [km] | RMS 0–7 d [km] | 24 h [km] | 72 h [km] | 7 dias [km] |",  # noqa: E501
        "|---|---:|---:|---:|---:|---:|",
    ]
    default = CD * DEFAULT_AREA_OVER_MASS
    for c in sorted({*SWEEP_TABLE, fitted, default}):
        k = sw.index(c)
        tag = " (só J2)" if c == 0.0 else " (**ajustado**)" if c == fitted else (
            " (padrão)" if c == default else "")  # fmt: skip
        rms72 = _dec(sw.rms72_km[k], ".2f")
        if c == fitted:
            rms72 = f"**{rms72}**"
        cells = [rms72, _dec(sw.rms7d_km[k], ".2f"), *(_dec(x, ".1f") for x in sw.point_km[k])]
        lines.append(f"| {_dec(c, '.5f')}{tag} | {' | '.join(cells)} |")
    k0, kf, kd = sw.index(0.0), sw.index(fitted), sw.index(default)
    k7 = int(np.argmin(sw.rms7d_km))
    lines += [
        "",
        "![Erro em função do C_D·A/m](erro_vs_coeficiente.png)",
        "",
        "Em todas as colunas o erro tem forma de **V**: acrescentar arrasto ajuda enquanto a "
        "quantidade se aproxima da certa; passando dela, o erro volta a crescer, e rápido. O que "
        "importa é **quão errado está o arrasto**, não se ele está presente:",
        "",
        f"- **Só J2** erra para *menos*: falta arrasto (0 em vez de ≈ {_dec(fitted, '.4f')}); "
        f"RMS de {_dec(sw.rms7d_km[k0], '.1f')} km em 7 dias.",
        f"- **O arrasto padrão** erra para *mais*: {_dec(default, '.4f')}, cerca de "
        f"{_dec(default / fitted, '.1f')} vezes o arrasto que o SGP4 “enxerga”; RMS de "
        f"{_dec(sw.rms7d_km[kd], '.1f')} km em 7 dias.",
        "",
        "Errar para mais por esse fator é pior do que não ter arrasto algum: um valor só ≈ 50% "
        "acima do ajustado (0,0025) já empata com o caso sem arrasto no erro em 7 dias.",
        "",
        "**Onde fica o mínimo depende da janela.** O RMS em 0–72 h, que é o que "
        f"`fit_ballistic_coefficient` minimiza, tem o mínimo no valor ajustado "
        f"({_dec(fitted, '.5f')}, {_dec(sw.rms72_km[kf], '.2f')} km). Em 0–7 dias o melhor valor é "
        f"outro (≈ {_dec(sw.cd_a_m[k7], '.4f')}, RMS de {_dec(sw.rms7d_km[k7], '.2f')} km contra "
        f"{_dec(sw.rms7d_km[kf], '.2f')} km do ajustado). E os erros pontuais (24 h, 72 h, 7 dias) "
        "oscilam com a fase da órbita, então o seu mínimo pode cair ao lado do ajustado. Se o "
        "arrasto fosse o único erro do modelo, todas as janelas apontariam para o mesmo valor; "
        "como não apontam, o C_D·A/m ajustado é um parâmetro *efetivo*, que também compensa o que "
        "o Cowell não tem (ver a seção seguinte).",
        "",
        "**Por que o erro de arrasto cresce tão depressa.** O arrasto tira energia, a órbita baixa "
        "e uma órbita mais baixa tem período menor: o satélite *adianta* ao longo da trajetória. "
        "O semieixo cai de forma aproximadamente linear no tempo, mas o adiantamento acumula essa "
        "diferença de velocidade angular, de modo que o erro ao longo da órbita cresce com o "
        "**quadrado do tempo**. Um erro pequeno em C_D·A/m vira centenas de km em uma semana.",
        "",
        "**Por que o SGP4 “vê” tão pouco arrasto.** O arrasto do SGP4 vem do **B\\*** do TLE, que "
        "o NORAD ajustou às medições de radar daquela semana. Em dezembro de 2019 a atividade "
        "solar estava no mínimo, a termosfera estava contraída e o arrasto real era fraco; o B\\* "
        "registra isso. A tabela de Vallado representa atividade *moderada*: não sabe que era "
        "mínimo solar e põe muito mais ar a 400 km do que havia (ver a seção seguinte).",
        "",
        "**Com o valor certo, o arrasto ajuda.** Com o C_D·A/m ajustado o RMS em 7 dias cai de "
        f"{_dec(sw.rms7d_km[k0], '.1f')} km (só J2) para {_dec(sw.rms7d_km[kf], '.1f')} km, "
        f"{_dec(sw.rms7d_km[k0] / sw.rms7d_km[kf], '.1f')} vezes menos. Mais física deixa o "
        "modelo mais fiel *se os parâmetros forem certos*; com parâmetros genéricos, um modelo "
        "mais completo pode ser pior do que um mais simples.",
        "",
        "## Por que o arrasto padrão piora o resultado",
        "",
        f"Com os parâmetros padrão (C_D = 2,2, A/m = 0,004) o Cowell fica **mais longe** do SGP4 "
        f"do que só com J2: {_km(dflt[0])} / {_km(dflt[1])} / {_km(dflt[2])} km contra "
        f"{_km(j2[0])} / {_km(j2[1])} / {_km(j2[2])} km em 24 h / 72 h / 7 dias. Não é erro de "
        "código: o arrasto está implementado e validado contra o Sem3 (1,2e-6 m/s² a 400 km, "
        "−11,8 m por revolução, Tabela 3 e tempos de vida). O problema é o *valor* do arrasto:",
        "",
        "- **Mínimo solar.** Dezembro de 2019 foi o fundo do mínimo do ciclo solar 24/25. Com "
        "pouca radiação EUV a termosfera esfria e se contrai, e a densidade a 400 km fica várias "
        "vezes abaixo da média. A tabela de Vallado representa atividade solar *moderada* "
        "(Sem3 Tabela 3: ρ varia por fatores de 5 a 10 ao longo do ciclo), então superestima o "
        "arrasto real daquela semana.",
        "- **A/m efetivo da ISS.** A área frontal depende da atitude e da orientação dos painéis "
        "solares, que giram para seguir o Sol; 0,004 m²/kg é um valor nominal.",
        "- **O que o ajuste absorve.** O C_D·A/m ajustado "
        f"({_dec(fitted, '.5f')}, cerca de {fitted / (CD * DEFAULT_AREA_OVER_MASS):.0%} do "
        "padrão) é um "
        "parâmetro *efetivo*: compensa a atmosfera do dia, a atitude e também parte do que o "
        "Cowell não modela e o SGP4 sim (J3, J4 e, em menor grau, Sol e Lua). Por isso o erro "
        "ajustado não vai a zero e cresce com o tempo.",
        "",
        "A diferença entre as linhas esférica e geodésica do arrasto padrão vem da altitude usada "
        "na densidade: medida sobre o elipsoide WGS84, a altitude da ISS é maior que pela esfera "
        "de raio R⊕ (de 0 no equador a ≈ 13 km nas latitudes máximas, 51,6°), e a densidade ali é "
        "menor, então o arrasto geodésico é um pouco mais fraco.",
        "",
        "**Lição:** um modelo de arrasto correto com parâmetros genéricos pode ser pior do que "
        "nenhum arrasto. Na prática, operadores ajustam um coeficiente balístico efetivo contra "
        "dados de rastreio (Sem3 §7, item 4), exatamente o que `fit_ballistic_coefficient` faz "
        "aqui com o SGP4 no lugar dos radares. A observação sobre atividade solar volta no "
        "documento da Semana 3.",
        "",
        "## Fora do escopo",
        "",
        "- Atmosferas empíricas (NRLMSISE-00, JB2008) e variação com o índice F10.7.",
        "- Tempestades geomagnéticas.",
        "- J3/J4 e atração lunissolar no Cowell.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tle-file", type=Path, default=DEFAULT_TLE, help="TLE (2 ou 3 linhas)")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    start = time.perf_counter()
    tle = Tle.from_text(args.tle_file.read_text(encoding="utf-8"))
    t, results, fitted, fit_s, sw = run(tle)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    plot(t, results, tle, args.out_dir / "erro_sgp4_cowell.png")
    plot_sweep(sw, fitted, args.out_dir / "erro_vs_coeficiente.png")
    write_report(t, results, tle, fitted, fit_s, sw, time.perf_counter() - start,
                 args.out_dir / "relatorio.md")  # fmt: skip
    for res in results:
        errs = ", ".join(f"{res.at_hours(t, h) / 1e3:.2f}" for h in CHECK_H)
        logger.info("%-50s %s km", res.case.label, errs)
    logger.info("gerado em %s (%.1f s)", args.out_dir, time.perf_counter() - start)


if __name__ == "__main__":
    main()
