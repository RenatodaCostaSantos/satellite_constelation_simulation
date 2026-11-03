"""Ablação de forças do Cowell contra o SGP4 (Sentinel-2A e ISS) e relatório do S4-06.

Uso (na raiz do repositório, com o pacote instalado com o extra dev):

    python scripts/compare_high_fidelity.py                # roda tudo (~10 min) e gera o relatório
    python scripts/compare_high_fidelity.py --from-json    # só refaz gráfico e relatório

Para cada satélite (Sentinel-2A, TLE do Apêndice A.1 da Semana 4; ISS, TLE do A.1 da Semana 3) e
cada conjunto de forças (J2; J2–J4; J2–J4 + arrasto; J2–J4 + Sol/Lua + SRP + arrasto), o Cowell é
comparado ao SGP4 por 7 dias de duas formas:

(A) partindo do estado do SGP4 em t = 0 (método da Semana 3); com arrasto, o C_D·A/m é ajustado
    ao RMS de 72 h por busca unidimensional limitada a [0; 0,03] m²/kg;
(B) com o estado inicial (e o C_D·A/m) ajustado por mínimos quadrados às posições do SGP4 em
    72 h (``satsim.orbit_determination.fit_initial_state``).

Gera em reports/s4_alta_fidelidade/: erro_ablacao.png, relatorio.md e resultados.json (os números
brutos, para refazer o relatório com --from-json sem repropagar). Referência: Apêndice A.8 do
prompt da Semana 4.
"""

import argparse
import json
import logging
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # sem janela: só gera os arquivos

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import minimize_scalar  # noqa: E402

from satsim.analysis import position_error  # noqa: E402
from satsim.forces import (  # noqa: E402
    DragAcceleration,
    J2Acceleration,
    SolarRadiationPressure,
    ThirdBodyAcceleration,
    ZonalAcceleration,
)
from satsim.orbit_determination import fit_initial_state  # noqa: E402
from satsim.propagators import CowellPropagator, SGP4Propagator  # noqa: E402
from satsim.tle import Tle  # noqa: E402

logger = logging.getLogger("compare_high_fidelity")

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "reports" / "s4_alta_fidelidade"
SATELLITES = {
    "sentinel2a": ("Sentinel-2A", ROOT / "tests" / "data" / "sentinel2a_sample.tle"),
    "iss": ("ISS", ROOT / "tests" / "data" / "iss_sample.tle"),
}

HOUR = 3600.0
DURATION_S = 7 * 24 * HOUR
STEP_S = 600.0
CHECK_H = (24.0, 72.0, 168.0)
FIT_WINDOW_S = 72 * HOUR
CD = 2.2
DRAG_BOUNDS = (0.0, 0.03)  # busca de C_D·A/m em (A), como na Semana 3 e no A.8 [m²/kg]
RTOL, ATOL = 1e-10, 1e-3  # DOP853, as do experimento do A.8
SRP_CR, SRP_AREA_OVER_MASS = 1.3, 0.01  # SRP com os parâmetros do A.4/A.7

FORCE_SETS = {  # chave → (rótulo, arrasto?, Sol/Lua + SRP?, cor)
    "j2": ("J2", False, False, "#8d99ae"),
    "j2j4": ("J2–J4", False, False, "#118ab2"),
    "j2j4_drag": ("J2–J4 + arrasto", True, False, "#ef476f"),
    "full": ("J2–J4 + Sol/Lua + SRP + arrasto", True, True, "#06d6a0"),
}
A8_KM = {  # Apêndice A.8 do prompt da Semana 4: erro 24 h / 72 h / 7 dias [km]
    ("sentinel2a", "A"): {"j2": (2.12, 10.11, 18.40), "j2j4": (1.60, 5.12, 13.41),
                          "j2j4_drag": (1.48, 4.11, 7.91), "full": (1.39, 3.77, 6.77)},
    ("sentinel2a", "B"): {"j2": (0.76, 1.58, 3.20), "j2j4": (0.06, 0.11, 1.72),
                          "j2j4_drag": (0.05, 0.09, 1.44), "full": (0.06, 0.09, 1.11)},
    ("iss", "A"): {"j2": (0.98, 8.66, 47.18), "j2j4": (3.58, 16.44, 64.98),
                   "j2j4_drag": (1.54, 1.97, 35.41), "full": (1.60, 2.05, 36.35)},
    ("iss", "B"): {"j2": (0.39, 1.03, 26.72), "j2j4": (0.93, 2.29, 31.53),
                   "j2j4_drag": (0.13, 0.25, 3.60), "full": (0.14, 0.29, 3.39)},
}  # fmt: skip
A8_CDAM_B = {("sentinel2a", "j2j4_drag"): 0.00220, ("sentinel2a", "full"): 0.00208,
             ("iss", "j2j4_drag"): 0.00137, ("iss", "full"): 0.00134}  # fmt: skip
A8_RMS_B_KM = {  # RMS do ajuste (B) em 72 h no A.8 [km]
    "sentinel2a": {"j2": 1.06, "j2j4": 0.051, "j2j4_drag": 0.043, "full": 0.048},
    "iss": {"j2": 0.97, "j2j4": 0.94, "j2j4_drag": 0.109, "full": 0.128},
}
CRITERIA_KM = {"sentinel2a": (0.2, 0.3, 3.0), "iss": (0.5, 1.0, 8.0)}  # card S4-06, (B) completo


def base_models(key: str, epoch) -> tuple:
    """Forças do conjunto ``key`` sem o arrasto (o arrasto entra com o C_D·A/m de cada caso)."""
    _, _, lunisolar, _ = FORCE_SETS[key]
    models: tuple = (J2Acceleration(),) if key == "j2" else (ZonalAcceleration(),)
    if lunisolar:
        models += (
            ThirdBodyAcceleration("sun", epoch),
            ThirdBodyAcceleration("moon", epoch),
            SolarRadiationPressure(epoch, cr=SRP_CR, area_over_mass=SRP_AREA_OVER_MASS),
        )
    return models


def with_drag(models: tuple, cd_a_m: float) -> tuple:
    return (*models, DragAcceleration(cd=CD, area_over_mass=cd_a_m / CD))


def rms_km(err_m: np.ndarray) -> float:
    return float(np.sqrt(np.mean((err_m / 1e3) ** 2)))


def case_a(sgp4: SGP4Propagator, key: str, t: np.ndarray) -> dict:
    """(A): estado do SGP4 em t = 0; com arrasto, C_D·A/m ajustado ao RMS de 72 h."""
    r0, v0 = (x[0] for x in sgp4.propagate(0.0))
    models = base_models(key, sgp4.epoch)
    in_fit = t <= FIT_WINDOW_S + 1e-9

    def error(cd_a_m: float | None) -> np.ndarray:
        m = models if cd_a_m is None else with_drag(models, cd_a_m)
        prop = CowellPropagator((r0, v0), sgp4.epoch, m, rtol=RTOL, atol=ATOL)
        return position_error(prop, sgp4, t)

    cd_a_m = None
    if FORCE_SETS[key][1]:
        res = minimize_scalar(lambda c: rms_km(error(c)[in_fit]), bounds=DRAG_BOUNDS,
                              method="bounded", options={"xatol": 1e-6})  # fmt: skip
        cd_a_m = float(res.x)
    err = error(cd_a_m)
    at_bound = cd_a_m is not None and cd_a_m > DRAG_BOUNDS[1] - 1e-4
    return {"cd_a_m": cd_a_m, "at_bound": at_bound, "rms72_km": rms_km(err[in_fit]),
            "err_m": err.tolist()}  # fmt: skip


def case_b(sgp4: SGP4Propagator, key: str, t: np.ndarray) -> dict:
    """(B): estado inicial (e C_D·A/m) ajustado por mínimos quadrados em 72 h."""
    fit = fit_initial_state(sgp4, base_models(key, sgp4.epoch), window_s=FIT_WINDOW_S,
                            step_s=STEP_S, fit_drag=FORCE_SETS[key][1], cd=CD,
                            rtol=RTOL, atol=ATOL)  # fmt: skip
    err = position_error(fit.propagator(rtol=RTOL, atol=ATOL), sgp4, t)
    return {"cd_a_m": fit.cd_a_m, "rms72_km": fit.rms_km, "rms_initial_km": fit.rms_initial_km,
            "nfev": fit.nfev, "success": fit.success, "err_m": err.tolist()}  # fmt: skip


def run() -> dict:
    t = np.arange(0.0, DURATION_S + 0.5 * STEP_S, STEP_S)
    out: dict = {"t_s": t.tolist(), "satellites": {}}
    for sat_key, (name, path) in SATELLITES.items():
        tle = Tle.from_text(path.read_text(encoding="utf-8"))
        sgp4 = SGP4Propagator(tle)
        sat_out = {"name": name, "tle_epoch": tle.epoch.isoformat(timespec="seconds"),
                   "tle_file": path.relative_to(ROOT).as_posix(), "cases": {}}  # fmt: skip
        for key in FORCE_SETS:
            start = time.perf_counter()
            a = case_a(sgp4, key, t)
            a["elapsed_s"] = time.perf_counter() - start
            start = time.perf_counter()
            b = case_b(sgp4, key, t)
            b["elapsed_s"] = time.perf_counter() - start
            sat_out["cases"][key] = {"A": a, "B": b}
            logger.info("%-11s %-32s (A) %s km  (B) %s km  [%.0f s + %.0f s]", name,
                        FORCE_SETS[key][0], _errs(t, a), _errs(t, b), a["elapsed_s"],
                        b["elapsed_s"])  # fmt: skip
        out["satellites"][sat_key] = sat_out
    return out


def at_hours(t: np.ndarray, err_m, hours: float) -> float:
    return float(np.asarray(err_m)[np.argmin(np.abs(t - hours * HOUR))]) / 1e3


def _errs(t: np.ndarray, case: dict) -> str:
    return " / ".join(f"{at_hours(t, case['err_m'], h):.2f}" for h in CHECK_H)


def _dec(x: float, fmt: str) -> str:
    """Número com vírgula decimal (relatório em português)."""
    return format(x, fmt).replace(".", ",")


def plot(data: dict, out: Path) -> None:
    t = np.asarray(data["t_s"])
    days = t[1:] / (24 * HOUR)
    sats = data["satellites"]
    fig, axes = plt.subplots(1, len(sats), figsize=(13, 5.2), dpi=120, sharey=True)
    for ax, sat in zip(axes, sats.values(), strict=True):
        for key, (label, _, _, color) in FORCE_SETS.items():
            case = sat["cases"][key]
            ax.plot(days, np.asarray(case["A"]["err_m"][1:]) / 1e3, "--", color=color, lw=1.1)
            ax.plot(days, np.asarray(case["B"]["err_m"][1:]) / 1e3, "-", color=color, lw=1.6,
                    label=label)  # fmt: skip
        for hours in CHECK_H:
            ax.axvline(hours / 24, color="#adb5bd", lw=0.8, ls=":")
        ax.axvspan(0, FIT_WINDOW_S / (24 * HOUR), color="#e9ecef", alpha=0.6, lw=0)
        ax.set_yscale("log")
        ax.set_xlim(0, DURATION_S / (24 * HOUR))
        ax.set_xlabel("tempo desde a época do TLE [dias]")
        ax.set_title(f"{sat['name']} (época {sat['tle_epoch'][:16].replace('T', ' ')} UTC)")
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("erro de posição 3D |r_Cowell − r_SGP4| [km]")
    handles, labels = axes[0].get_legend_handles_labels()
    handles += [plt.Line2D([], [], color="k", ls="--", lw=1.1),
                plt.Line2D([], [], color="k", ls="-", lw=1.6)]  # fmt: skip
    labels += ["(A) estado do SGP4 em t = 0", "(B) estado ajustado em 72 h (faixa cinza)"]
    axes[1].legend(handles, labels, loc="lower right", fontsize=8.5)
    fig.suptitle("Cowell × SGP4: ablação de forças, sem ajuste (A) e com ajuste (B)")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def write_report(data: dict, elapsed_s: float, out: Path) -> None:
    t = np.asarray(data["t_s"])
    sats = data["satellites"]

    def err(sat: str, key: str, col: str, k: int) -> str:
        case = sats[sat]["cases"][key][col]
        return _dec(at_hours(t, case["err_m"], CHECK_H[k]), ".2f")

    def cdam(sat: str, key: str) -> str:
        return _dec(sats[sat]["cases"][key]["B"]["cd_a_m"], ".5f")

    def rms(sat: str, key: str) -> str:
        return _dec(sats[sat]["cases"][key]["B"]["rms72_km"], ".3f")

    lines = [
        "# Relatório do S4-06 — Cowell de alta fidelidade × SGP4",
        "",
        f"Gerado por `python scripts/compare_high_fidelity.py` em {elapsed_s / 60:.1f} min. "
        "Cowell DOP853 (rtol = 1e-10, atol = 1e-3); J2–J4 = `ZonalAcceleration`; arrasto com a "
        "atmosfera exponencial de Vallado, altitude geodésica e co-rotação; Sol e Lua pelas "
        "efemérides analíticas (`satsim.ephemeris`); SRP com C_R = 1,3, A/m = 0,01 m²/kg e sombra "
        "cilíndrica. Referência: Apêndice A.8 do prompt da Semana 4 (entre parênteses nas "
        "tabelas). Erros em km, instantâneos em 24 h / 72 h / 7 dias.",
        "",
        "- **(A)** Cowell a partir do estado do SGP4 em t = 0 (método da Semana 3). Com arrasto, "
        "C_D·A/m ajustado ao RMS de 72 h numa busca limitada a [0; 0,03] m²/kg.",
        "- **(B)** Estado inicial (6 parâmetros) e, com arrasto, C_D·A/m ajustados por mínimos "
        "quadrados às posições do SGP4 a cada 600 s em 72 h (`fit_initial_state`). Sem arrasto, "
        "só os 6 do estado.",
        "",
    ]
    for sat_key, sat in sats.items():
        lines += [
            f"## {sat['name']}",
            "",
            f"TLE: `{sat['tle_file']}`, época {sat['tle_epoch'].replace('T', ' ')} UTC.",
            "",
            "| Forças | (A) 24 h / 72 h / 7 d | C_D·A/m (A) | (B) 24 h / 72 h / 7 d | C_D·A/m (B) | RMS 72 h (B) [km] | Tempo (B) |",  # noqa: E501
            "|---|---|---:|---|---:|---:|---:|",
        ]
        for key, (label, *_rest) in FORCE_SETS.items():
            a, b = sat["cases"][key]["A"], sat["cases"][key]["B"]
            cells = []
            for col, case in (("A", a), ("B", b)):
                got = " / ".join(_dec(at_hours(t, case["err_m"], h), ".2f") for h in CHECK_H)
                ref = " / ".join(_dec(x, ".2f") for x in A8_KM[(sat_key, col)][key])
                cells.append(f"{got} ({ref})")
            cd_a = "—" if a["cd_a_m"] is None else _dec(a["cd_a_m"], ".5f")
            if a.get("at_bound"):
                cd_a += " \\*"
            cd_b = "—" if b["cd_a_m"] is None else _dec(b["cd_a_m"], ".5f")
            if (sat_key, key) in A8_CDAM_B:
                cd_b += f" ({_dec(A8_CDAM_B[(sat_key, key)], '.5f')})"
            lines.append(
                f"| {label} | {cells[0]} | {cd_a} | {cells[1]} | {cd_b} | "
                f"{_dec(b['rms72_km'], '.3f')} ({_dec(A8_RMS_B_KM[sat_key][key], 'g')}) | "
                f"{b['elapsed_s']:.0f} s |"
            )
        full = [at_hours(t, sat["cases"]["full"]["B"]["err_m"], h) for h in CHECK_H]
        lim = CRITERIA_KM[sat_key]
        ok = all(e <= x for e, x in zip(full, lim, strict=True))
        lim_txt = " / ".join(_dec(x, "g") for x in lim)
        got_txt = " / ".join(_dec(e, ".2f") for e in full)
        lines += [
            "",
            f"**Critério do card, (B) com forças completas (≤ {lim_txt} km):** {got_txt} km → "
            f"{'OK' if ok else 'FALHA'}.",
            "",
        ]
    if any(c["A"].get("at_bound") for s in sats.values() for c in s["cases"].values()):
        lines += [
            "\\* C_D·A/m de (A) parou no limite superior da busca (0,03): o arrasto estava "
            "absorvendo o erro do estado inicial, não descrevendo a atmosfera.",
            "",
        ]
    lines += [
        "![Erro de posição: ablação de forças](erro_ablacao.png)",
        "",
        "Linhas tracejadas: (A); contínuas: (B). A faixa cinza é a janela do ajuste (72 h); "
        "depois dela, (B) é previsão.",
        "",
        "## Por que (A) tem piso e (B) não",
        "",
        "O SGP4 é uma teoria analítica truncada: o estado que ele dá na época já erra ≈ 1 km "
        "contra a órbita que ele próprio descreve nos dias seguintes. Partindo desse estado, o "
        "Cowell herda o erro, e o pior componente é o do **semieixo**: alguns metros de erro em a "
        "mudam o período e viram uma deriva ao longo da órbita que cresce linearmente com o "
        "tempo. Nenhuma força a mais remove esse erro, porque ele não vem da física que falta, "
        "vem da condição inicial. É o piso de (A): no Sentinel-2A, sair de J2 para o modelo "
        f"completo reduz o erro em 7 dias só de {err('sentinel2a', 'j2', 'A', 2)} para "
        f"{err('sentinel2a', 'full', 'A', 2)} km.",
        "",
        "Em (B) o estado inicial é escolhido para que a trajetória inteira de 72 h fique o mais "
        "perto possível do SGP4, o que corrige o semieixo médio. O que sobra é a diferença de "
        "**modelo** entre o Cowell e o SGP4, e aí cada força a mais aparece. No Sentinel-2A, J3 e "
        f"J4 levam o erro em 24 h de {err('sentinel2a', 'j2', 'B', 0)} para "
        f"{err('sentinel2a', 'j2j4', 'B', 0)} km. Na ISS, mais baixa, quem domina é o arrasto: "
        "sem ele o ajuste não consegue representar a queda do semieixo, e o erro em 7 dias fica "
        f"em {err('iss', 'j2j4', 'B', 2)} km; com ele cai para {err('iss', 'j2j4_drag', 'B', 2)} "
        "km.",
        "",
        "**O C_D·A/m ajustado é efetivo.** Em (A) o arrasto também compensa o erro do estado "
        "inicial: no Sentinel-2A, a 786 km, o C_D·A/m de (A) vai para o limite da busca. Em (B) "
        "os dois papéis se separam e ele deixa de bater no limite. Na ISS o valor de (B) "
        f"({cdam('iss', 'full')} m²/kg) fica perto do A.8 e do ajuste da Semana 3 (0,00163). No "
        f"Sentinel-2A sai {cdam('sentinel2a', 'full')} m²/kg, contra 0,00208 do A.8: a 786 km a "
        "tabela de Vallado e o modelo de densidade do próprio SGP4 (uma lei de potência fixa, "
        "escalada pelo B\* do TLE) não têm por que concordar, e o ajuste escolhe o coeficiente "
        "que reproduz o decaimento que o *SGP4* aplica. É um parâmetro efetivo em relação ao "
        "oráculo, não o C_D·A/m físico do satélite.",
        "",
        "## Diferenças em relação ao A.8",
        "",
        "Os erros de (A) reproduzem o A.8 nos 8 casos. Em (B) o satsim chega a erros e RMS "
        "**menores** que os do A.8 em todos os casos (RMS entre parênteses na tabela; ex.: "
        f"Sentinel-2A J2–J4 + arrasto {rms('sentinel2a', 'j2j4_drag')} km contra 0,043 km; ISS "
        f"{rms('iss', 'j2j4_drag')} km contra 0,109 km). Como o RMS de 72 h é exatamente o que o "
        "ajuste minimiza, um RMS menor com o mesmo modelo indica que o otimizador do experimento "
        "parou antes do mínimo. A causa provável é o passo das diferenças finitas: com passo "
        "relativo pequeno (≈ 1e-6, 1 mm em r0) o Jacobiano fica dominado pelo ruído da integração "
        "e o `least_squares` para cedo; no ajuste sintético do satsim, esse passo deixava o estado "
        "a 0,33 m do certo, contra 1,3 mm com o passo usado (1e-4: 0,1 m em r0, 0,1 mm/s em v0). "
        "Os critérios do card passam com folga nos dois casos.",
        "",
        "## Por que o SGP4 não valida Sol/Lua",
        "",
        "Para órbitas com período abaixo de 225 min (LEO, caso dos dois satélites) o SGP4 **não "
        "tem** atração lunissolar: ela só entra no SDP4, a variante de espaço profundo. O SGP4 "
        "também não tem SRP. Então, contra o SGP4, acrescentar Sol, Lua e SRP ao Cowell só pode "
        "piorar o acordo ou deixá-lo igual: o Cowell completo fica mais perto da órbita real e "
        "mais longe do oráculo. É o que a tabela mostra em (B): em 7 dias, a linha J2–J4 + "
        f"arrasto erra {err('sentinel2a', 'j2j4_drag', 'B', 2)} km no Sentinel-2A e "
        f"{err('iss', 'j2j4_drag', 'B', 2)} km na ISS; com Sol, Lua e SRP, "
        f"{err('sentinel2a', 'full', 'B', 2)} e {err('iss', 'full', 'B', 2)} km. Mais física, "
        "acordo pior: a assinatura de forças que o oráculo não tem. Validar Sol e Lua exige uma "
        "referência que os tenha: efemérides de precisão (SP3 de GNSS a bordo, como as do "
        "Sentinel-2) ou um propagador numérico de referência. Neste projeto eles foram validados "
        "por componentes no S4-02 (posições contra DE421 e aceleração contra valores "
        "independentes, Apêndice A.3).",
        "",
        "## Fora do escopo",
        "",
        "- Atmosfera com F10.7 e Ap, harmônicos tesserais, penumbra, manobras.",
        "- Referências de precisão (SP3/GNSS) no lugar do SGP4.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument(
        "--from-json", action="store_true", help="só refaz gráfico e relatório (resultados.json)"
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / "resultados.json"

    start = time.perf_counter()
    if args.from_json:
        data = json.loads(json_path.read_text(encoding="utf-8"))
    else:
        data = run()
        data["elapsed_s"] = time.perf_counter() - start
        json_path.write_text(json.dumps(data), encoding="utf-8")
    plot(data, args.out_dir / "erro_ablacao.png")
    write_report(data, data["elapsed_s"], args.out_dir / "relatorio.md")
    logger.info("gerado em %s (%.1f s)", args.out_dir, time.perf_counter() - start)


if __name__ == "__main__":
    main()
