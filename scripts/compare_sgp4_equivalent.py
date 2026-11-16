"""Ablação a partir do "Cowell equivalente ao SGP4": quanto da diferença Cowell × SGP4 é física.

Uso (na raiz do repositório, com o pacote instalado com o extra dev):

    python scripts/compare_sgp4_equivalent.py               # roda tudo (~10 min) e gera o relatório
    python scripts/compare_sgp4_equivalent.py --from-json   # só refaz gráfico e relatório

A ablação do S4-06 (``compare_high_fidelity.py``) parte de um Cowell que difere do SGP4 em várias
coisas ao mesmo tempo (constantes, atmosfera, co-rotação, Sol/Lua/SRP), e o erro que sobra mistura
todas. Aqui o ponto de partida é o Cowell com **as mesmas premissas físicas do SGP4**
(``satsim.sgp4_equivalent``): o erro do passo 0 mede só a truncagem da teoria analítica do SGP4
e o acordo do integrador numérico com ela. Os passos seguintes acrescentam uma diferença por vez:

0. premissas do SGP4: zonais J2–J4 WGS-72, µ WGS-72, arrasto do B* (lei de potência, sem
   co-rotação), sem Sol/Lua/SRP;
1. constantes do projeto (µ, R⊕, J2–J4 de ``satsim.constants``), mesmo arrasto do B*;
2. atmosfera do projeto (Vallado, altitude geodésica, co-rotação) com C_D·A/m ajustado: é o
   "J2–J4 + arrasto" do S4-06;
3. + Sol, Lua e SRP: é o "completo" do S4-06.

Satélites e métodos como no S4-06 (Sentinel-2A e ISS, 7 dias a cada 600 s): (A) estado do SGP4 em
t = 0 (C_D·A/m por busca 1D nos passos com atmosfera do projeto); (B) estado inicial ajustado por
mínimos quadrados às posições do SGP4 em 72 h (``fit_initial_state``), com o C_D·A/m ajustado
junto nos passos 2 e 3 e o arrasto fixo pelo B* nos passos 0 e 1.

Gera em reports/sgp4_equivalente/: erro_sgp4_equivalente.png, relatorio.md e resultados.json.
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
from satsim.constants import MU  # noqa: E402
from satsim.forces import (  # noqa: E402
    DragAcceleration,
    SolarRadiationPressure,
    ThirdBodyAcceleration,
    ZonalAcceleration,
)
from satsim.orbit_determination import fit_initial_state  # noqa: E402
from satsim.propagators import CowellPropagator, SGP4Propagator  # noqa: E402
from satsim.sgp4_equivalent import (  # noqa: E402
    WGS72_MU,
    perigee_altitude_from_tle,
    sgp4_drag,
    sgp4_gravity,
)
from satsim.tle import Tle  # noqa: E402

logger = logging.getLogger("compare_sgp4_equivalent")

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "reports" / "sgp4_equivalente"
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
DRAG_BOUNDS = (0.0, 0.03)  # busca de C_D·A/m em (A), como no S4-06 [m²/kg]
RTOL, ATOL = 1e-10, 1e-3  # DOP853, as do S4-06
SRP_CR, SRP_AREA_OVER_MASS = 1.3, 0.01

STEPS = {  # chave → (rótulo, cor)
    "sgp4eq": ("0 · premissas do SGP4", "#2a78d6"),
    "consts": ("1 · + constantes do projeto", "#eb6834"),
    "atmos": ("2 · + atmosfera do projeto", "#1baf7a"),
    "full": ("3 · + Sol, Lua e SRP", "#eda100"),
}
FITTED_DRAG = {"atmos", "full"}  # passos com C_D·A/m ajustado (atmosfera do projeto)


def step_models(key: str, tle: Tle) -> tuple[tuple, float]:
    """Forças (sem o arrasto ajustado dos passos 2 e 3) e µ do termo central do passo ``key``."""
    epoch = tle.epoch
    if key == "sgp4eq":
        perigee = perigee_altitude_from_tle(tle.mean_motion_rev_day, tle.eccentricity)
        return (sgp4_gravity(), sgp4_drag(tle.bstar, perigee)), WGS72_MU
    if key == "consts":
        perigee = perigee_altitude_from_tle(tle.mean_motion_rev_day, tle.eccentricity)
        return (ZonalAcceleration(), sgp4_drag(tle.bstar, perigee)), MU
    models: tuple = (ZonalAcceleration(),)
    if key == "full":
        models += (
            ThirdBodyAcceleration("sun", epoch),
            ThirdBodyAcceleration("moon", epoch),
            SolarRadiationPressure(epoch, cr=SRP_CR, area_over_mass=SRP_AREA_OVER_MASS),
        )
    return models, MU


def with_drag(models: tuple, cd_a_m: float) -> tuple:
    return (*models, DragAcceleration(cd=CD, area_over_mass=cd_a_m / CD))


def rms_km(err_m: np.ndarray) -> float:
    return float(np.sqrt(np.mean((err_m / 1e3) ** 2)))


def case_a(sgp4: SGP4Propagator, tle: Tle, key: str, t: np.ndarray) -> dict:
    """(A): estado do SGP4 em t = 0; nos passos 2 e 3, C_D·A/m ajustado ao RMS de 72 h."""
    r0, v0 = (x[0] for x in sgp4.propagate(0.0))
    models, mu = step_models(key, tle)
    in_fit = t <= FIT_WINDOW_S + 1e-9

    def error(cd_a_m: float | None) -> np.ndarray:
        m = models if cd_a_m is None else with_drag(models, cd_a_m)
        prop = CowellPropagator((r0, v0), sgp4.epoch, m, mu=mu, rtol=RTOL, atol=ATOL)
        return position_error(prop, sgp4, t)

    cd_a_m = None
    if key in FITTED_DRAG:
        res = minimize_scalar(lambda c: rms_km(error(c)[in_fit]), bounds=DRAG_BOUNDS,
                              method="bounded", options={"xatol": 1e-6})  # fmt: skip
        cd_a_m = float(res.x)
    err = error(cd_a_m)
    return {"cd_a_m": cd_a_m, "rms72_km": rms_km(err[in_fit]), "err_m": err.tolist()}


def case_b(sgp4: SGP4Propagator, tle: Tle, key: str, t: np.ndarray) -> dict:
    """(B): estado inicial (e, nos passos 2 e 3, C_D·A/m) ajustado por mínimos quadrados."""
    models, mu = step_models(key, tle)
    fit = fit_initial_state(sgp4, models, window_s=FIT_WINDOW_S, step_s=STEP_S,
                            fit_drag=key in FITTED_DRAG, cd=CD, rtol=RTOL, atol=ATOL,
                            mu=mu)  # fmt: skip
    err = position_error(fit.propagator(rtol=RTOL, atol=ATOL), sgp4, t)
    return {"cd_a_m": fit.cd_a_m, "rms72_km": fit.rms_km, "rms_initial_km": fit.rms_initial_km,
            "nfev": fit.nfev, "success": fit.success, "err_m": err.tolist()}  # fmt: skip


def run() -> dict:
    t = np.arange(0.0, DURATION_S + 0.5 * STEP_S, STEP_S)
    out: dict = {"t_s": t.tolist(), "satellites": {}}
    for sat_key, (name, path) in SATELLITES.items():
        tle = Tle.from_text(path.read_text(encoding="utf-8"))
        sgp4 = SGP4Propagator(tle)
        perigee = perigee_altitude_from_tle(tle.mean_motion_rev_day, tle.eccentricity)
        sat_out = {"name": name, "tle_epoch": tle.epoch.isoformat(timespec="seconds"),
                   "tle_file": path.relative_to(ROOT).as_posix(), "bstar": tle.bstar,
                   "perigee_km": perigee / 1e3, "cases": {}}  # fmt: skip
        for key in STEPS:
            start = time.perf_counter()
            a = case_a(sgp4, tle, key, t)
            a["elapsed_s"] = time.perf_counter() - start
            start = time.perf_counter()
            b = case_b(sgp4, tle, key, t)
            b["elapsed_s"] = time.perf_counter() - start
            sat_out["cases"][key] = {"A": a, "B": b}
            logger.info("%-11s %-28s (A) %s km  (B) %s km  [%.0f s + %.0f s]", name,
                        STEPS[key][0], _errs(t, a), _errs(t, b), a["elapsed_s"],
                        b["elapsed_s"])  # fmt: skip
        out["satellites"][sat_key] = sat_out
    return out


def at_hours(t: np.ndarray, err_m, hours: float) -> float:
    """Erro [km] no instante mais próximo de ``hours``."""
    return float(np.asarray(err_m)[np.argmin(np.abs(t - hours * HOUR))]) / 1e3


def _errs(t: np.ndarray, case: dict) -> str:
    return " / ".join(f"{at_hours(t, case['err_m'], h):.3f}" for h in CHECK_H)


def _dec(x: float, fmt: str) -> str:
    """Número com vírgula decimal (relatório em português)."""
    return format(x, fmt).replace(".", ",")


def _km(x_km: float) -> str:
    """Distância legível: metros abaixo de 1 km, km com 2 casas acima."""
    return f"{_dec(x_km * 1e3, '.0f')} m" if x_km < 1.0 else f"{_dec(x_km, '.2f')} km"


def plot(data: dict, out: Path) -> None:
    t = np.asarray(data["t_s"])
    days = t[1:] / (24 * HOUR)
    sats = data["satellites"]
    fig, axes = plt.subplots(1, len(sats), figsize=(13, 5.8), dpi=120, sharey=True)
    for ax, sat in zip(axes, sats.values(), strict=True):
        for key, (label, color) in STEPS.items():
            case = sat["cases"][key]
            ax.plot(days, np.asarray(case["A"]["err_m"][1:]) / 1e3, "--", color=color, lw=1.1)
            ax.plot(days, np.asarray(case["B"]["err_m"][1:]) / 1e3, "-", color=color, lw=1.8,
                    label=label)  # fmt: skip
        ax.axvspan(0, FIT_WINDOW_S / (24 * HOUR), color="#e9ecef", alpha=0.6, lw=0)
        ax.set_yscale("log")
        ax.set_xlim(0, DURATION_S / (24 * HOUR))
        ax.set_xlabel("tempo desde a época do TLE [dias]")
        ax.set_title(f"{sat['name']} (época {sat['tle_epoch'][:16].replace('T', ' ')} UTC)")
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("erro de posição 3D |r_Cowell − r_SGP4| [km]")
    handles, labels = axes[0].get_legend_handles_labels()
    handles += [plt.Line2D([], [], color="k", ls="--", lw=1.1),
                plt.Line2D([], [], color="k", ls="-", lw=1.8)]  # fmt: skip
    labels += ["(A) estado do SGP4 em t = 0", "(B) estado ajustado em 72 h (faixa cinza)"]
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9, frameon=False)
    fig.suptitle("Cowell × SGP4: ablação a partir das premissas do SGP4")
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(out)
    plt.close(fig)


def write_report(data: dict, elapsed_s: float, out: Path) -> None:
    t = np.asarray(data["t_s"])
    sats = data["satellites"]

    def e(sat: str, key: str, col: str, k: int) -> str:
        return _km(at_hours(t, sats[sat]["cases"][key][col]["err_m"], CHECK_H[k]))

    def rms(sat: str, key: str) -> str:
        return _km(sats[sat]["cases"][key]["B"]["rms72_km"])

    lines = [
        "# Cowell equivalente ao SGP4: ablação a partir das premissas do SGP4",
        "",
        f"Gerado por `python scripts/compare_sgp4_equivalent.py` em {elapsed_s / 60:.1f} min. "
        "Cowell DOP853 (rtol = 1e-10, atol = 1e-3) contra o SGP4 dos TLEs de teste, 7 dias a cada "
        "600 s. Erros de posição 3D instantâneos em 24 h / 72 h / 7 dias.",
        "",
        "**Pergunta.** Quando o Cowell assume exatamente a mesma física que o SGP4, ele recupera o "
        "SGP4? Se sim, a diferença que sobra na ablação do S4-06 é física (o que um modelo tem e o "
        "outro não), e não erro do integrador.",
        "",
        "Passos (cada um acrescenta uma diferença ao anterior):",
        "",
        "0. **Premissas do SGP4** (`satsim.sgp4_equivalent`): zonais J2–J4 e µ do WGS-72; arrasto "
        "com a atmosfera de lei de potência do SGP4, ρ ∝ [(q0 − s)/(r − s)]⁴ (q0 = R⊕ + 120 km, "
        "s = R⊕ + 78 km), esférica, sem co-rotação, intensidade dada pelo B\\* do TLE; sem Sol, "
        "Lua e SRP.",
        "1. **Constantes do projeto**: µ, R⊕ e J2–J4 de `satsim.constants` no lugar dos do WGS-72; "
        "o mesmo arrasto do B\\*.",
        "2. **Atmosfera do projeto**: exponencial de Vallado, altitude geodésica e co-rotação, com "
        'C_D·A/m ajustado (é o "J2–J4 + arrasto" do S4-06).',
        '3. **+ Sol, Lua e SRP** (C_R = 1,3, A/m = 0,01 m²/kg; é o "completo" do S4-06).',
        "",
        "Métodos: **(A)** estado do SGP4 em t = 0; **(B)** estado inicial ajustado por mínimos "
        "quadrados às posições do SGP4 em 72 h (`fit_initial_state`). O C_D·A/m é ajustado só "
        "nos passos 2 e 3; nos passos 0 e 1 o arrasto é o do B\\*, sem parâmetro livre.",
        "",
    ]
    for sat in sats.values():
        lines += [
            f"## {sat['name']}",
            "",
            f"TLE: `{sat['tle_file']}`, época {sat['tle_epoch'].replace('T', ' ')} UTC, "
            f"B\\* = {_dec(sat['bstar'], '.4e')} /ER, perigeu ≈ {sat['perigee_km']:.0f} km.",
            "",
            "| Passo | (A) 24 h / 72 h / 7 d | (B) 24 h / 72 h / 7 d | RMS 72 h (B) | C_D·A/m (B) [m²/kg] |",  # noqa: E501
            "|---|---|---|---:|---:|",
        ]
        for key, (label, _) in STEPS.items():
            a, b = sat["cases"][key]["A"], sat["cases"][key]["B"]
            ca = " / ".join(_km(at_hours(t, a["err_m"], h)) for h in CHECK_H)
            cb = " / ".join(_km(at_hours(t, b["err_m"], h)) for h in CHECK_H)
            cd = "B\\* do TLE" if b["cd_a_m"] is None else _dec(b["cd_a_m"], ".5f")
            lines.append(f"| {label} | {ca} | {cb} | {_km(b['rms72_km'])} | {cd} |")
        lines.append("")
    lines += [
        "![Erro de posição: ablação a partir das premissas do SGP4](erro_sgp4_equivalente.png)",
        "",
        "Linhas tracejadas: (A); contínuas: (B). A faixa cinza é a janela do ajuste (72 h); "
        "depois dela, (B) é previsão.",
        "",
        "## O integrador recupera o SGP4",
        "",
        "Com as mesmas premissas (passo 0) e o estado ajustado, o Cowell fica a "
        f"{rms('sentinel2a', 'sgp4eq')} (RMS em 72 h) do SGP4 no Sentinel-2A e a "
        f"{rms('iss', 'sgp4eq')} na ISS, e em 7 dias, já fora da janela do ajuste, a "
        f"{e('sentinel2a', 'sgp4eq', 'B', 2)} e {e('iss', 'sgp4eq', 'B', 2)}. É o resultado "
        "esperado de duas soluções das mesmas equações: o que sobra é a truncagem da teoria "
        "analítica do SGP4 (J2 de 1ª ordem nos termos de curto período, J3 só no longo período, "
        "J4 só no secular, arrasto só como deriva secular), e não erro do integrador. Na ISS, a "
        "≈ 410 km, o resíduo de 7 dias é maior porque o arrasto, que o SGP4 trata de forma mais "
        "aproximada que a gravidade, pesa mais.",
        "",
        "No método (A) o passo 0 também é o melhor (Sentinel-2A em 7 dias: "
        f"{e('sentinel2a', 'sgp4eq', 'A', 2)}, contra {e('sentinel2a', 'consts', 'A', 2)} no passo "
        "1), mas continua com o piso do S4-06: o estado osculador que o SGP4 entrega em t = 0 é "
        "uma série truncada, e o erro de semieixo que ele carrega vira deriva ao longo da órbita. "
        "É a razão de o ajuste (B) ser o método de referência.",
        "",
        "## O que cada passo acrescenta",
        "",
        "- **Constantes (passo 1).** É o maior efeito em (A): Sentinel-2A "
        f"{e('sentinel2a', 'sgp4eq', 'A', 2)} → {e('sentinel2a', 'consts', 'A', 2)}, ISS "
        f"{e('iss', 'sgp4eq', 'A', 2)} → {e('iss', 'consts', 'A', 2)} em 7 dias. O estado do SGP4 "
        "em t = 0 é coerente com o µ do WGS-72; com o µ do projeto (≈ 0,9 ppm menor) o mesmo "
        "(r, v) tem outro semieixo e outro período, e a diferença vira deriva ao longo da órbita "
        "(≈ 570 m em 12 h no Sentinel-2A). Em (B) o ajuste absorve o µ no estado inicial e sobra "
        "o efeito dos zonais, sobretudo J4 (−1,656e-6 → −1,620e-6, 2%; J3 muda 0,2%): "
        f"Sentinel-2A {e('sentinel2a', 'sgp4eq', 'B', 2)} "
        f"→ {e('sentinel2a', 'consts', 'B', 2)}, ISS {e('iss', 'sgp4eq', 'B', 2)} → "
        f"{e('iss', 'consts', 'B', 2)}.",
        "- **Atmosfera (passo 2).** Troca o arrasto fixo do B\\* pela atmosfera de Vallado com "
        "co-rotação **e um parâmetro livre**, o C_D·A/m ajustado. No Sentinel-2A, a 786 km, nada "
        f"muda ({e('sentinel2a', 'consts', 'B', 2)} → {e('sentinel2a', 'atmos', 'B', 2)}): o "
        "arrasto é desprezível e o C_D·A/m ajustado não tem o que representar (sai "
        f"{_dec(sats['sentinel2a']['cases']['atmos']['B']['cd_a_m'], '.4f')} m²/kg, sem sentido "
        "físico). Na ISS o erro **cai** "
        f"({e('iss', 'consts', 'B', 2)} → {e('iss', 'atmos', 'B', 2)}). Não é física melhor: o "
        "parâmetro livre absorve parte da truncagem do arrasto do SGP4, que no passo 0 fica "
        "visível porque o arrasto está preso ao B\\*.",
        "- **Sol, Lua e SRP (passo 3).** Forças que o SGP4 de órbitas próximas não tem. Contra "
        "ele, a mudança mede a discordância com a referência, não o erro do Cowell: Sentinel-2A "
        f"{e('sentinel2a', 'atmos', 'B', 2)} → {e('sentinel2a', 'full', 'B', 2)}, ISS "
        f"{e('iss', 'atmos', 'B', 2)} → {e('iss', 'full', 'B', 2)} em 7 dias. É o maior efeito "
        "em (B), e o único que vem de física que o SGP4 não tem.",
        "",
        "## Como ler o resultado",
        "",
        "O SGP4 não é a verdade: é outro modelo, mais simples, com cerca de 1 km de erro contra a "
        "órbita real já na época. O passo 0 **verifica** o Cowell (com a mesma física, ele "
        "concorda com uma solução independente das mesmas equações); os passos 1 a 3 mostram "
        "quanto cada escolha de modelagem do projeto afasta o resultado do SGP4. **Validar** essas "
        "escolhas (dizer se elas aproximam da órbita real) exige uma referência medida, como as "
        "órbitas precisas (SP3) do Sentinel-2.",
        "",
        "## Limites da equivalência",
        "",
        "- O s da atmosfera do SGP4 depende do perigeu (muda abaixo de 156 km); aqui ele é "
        "decidido pelo perigeu do TLE, calculado com o movimento médio de Kozai. Os dois "
        "satélites ficam bem acima do limite.",
        "- O SGP4 usa elementos e semieixo de Brouwer e ignora as derivadas ṅ e n̈ do TLE; o "
        "Cowell não usa elementos médios em nenhum ponto, só as posições do SGP4.",
        "- A simplificação de arrasto do SGP4 para perigeus abaixo de 220 km (`isimp`) não é "
        "reproduzida.",
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
    plot(data, args.out_dir / "erro_sgp4_equivalente.png")
    write_report(data, data["elapsed_s"], args.out_dir / "relatorio.md")
    logger.info("gerado em %s (%.1f s)", args.out_dir, time.perf_counter() - start)


if __name__ == "__main__":
    main()
