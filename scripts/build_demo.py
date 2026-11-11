"""Gera demo/index.html: o viewer do satsim autocontido, que abre com duplo clique (file://).

Uso (na raiz do repositório, com o pacote instalado):

    python scripts/build_demo.py                      # ISS real, 24 h a partir de agora
    python scripts/build_demo.py --offline            # sem rede: cache ou TLE embarcado
    python scripts/build_demo.py --tle-file tests/data/iss_sample.tle --start tle-epoch
    python scripts/build_demo.py --days 6              # o ciclo de repetição inteiro (S5-06)
    python scripts/build_demo.py --sats 18             # 18 satélites → demo/constelacao18.html
    python scripts/build_demo.py --all --offline       # as quatro cenas, ligadas pelo seletor

Desde a S3-06 a demo mostra a **ISS real**: TLE de satsim.tle.get_tle (cache com menos de 2 h
→ CelesTrak → cache antigo → TLE embarcado) ou de --tle-file, propagado com SGP4Propagator. A
janela é de 24 h a partir de agora (UTC, arredondado ao minuto), ou da época do TLE com
--start tle-epoch (determinístico, para testes e demonstração), com passo de 20 s (1 dia).

Desde a S4-06:

- o **Sentinel-2A real** (SGP4, TLE de ``get_tle(40697)``, com TLE embarcado para o modo
  offline) entra na cena, desde que o TLE esteja a menos de ``MAX_TLE_OFFSET_DAYS`` do início
  da cena (com ``--start tle-epoch`` e o TLE de 2019 da ISS, ele fica de fora);
- a cena traz LTAN, β e eclipse de cada satélite (``build_scene(..., sun_geometry=True)``).

Desde a S5-06, a SSO-1 dá lugar à **constelação 91/6** (``build_constellation``): os 6
satélites de ``sso_constellation`` num único plano, no **Cowell completo** (J2–J4 + Sol/Lua +
SRP + arrasto), com a i calibrada do A.7 da Semana 4 (97,46325°), o a calibrado do S5-02
(``A_CALIBRATED_M``, que fecha o ciclo de 6 dias no Cowell) e LTAN 22h30 no início da cena. A
cena ganha ``meta.constellation`` e, em cada membro, os nós ascendentes, a faixa alcançável com
a rolagem de meia grade (≈ 23,2°) e o indicador de manutenção (S5-05). ISS e Sentinel-2A
continuam na cena; o viewer os deixa desligados na legenda quando há constelação.

``--days`` (1 a 6) estende a janela até o ciclo inteiro; o passo cresce com a janela
(``STEP_BY_DAYS``) para o arquivo ficar abaixo de 15 MB.

A cena traz também a estação de São Bento do Sapucaí, o que a ISS faz vista de lá (azimute,
elevação, distância), o círculo de visibilidade (footprint) a el_min = 10° e as passagens. Tudo
é calculado em Python e embutido no HTML junto com o CSS, o JS do viewer, o bundle UMD do
globe.gl e a textura da Terra (data URI). Nada é carregado da rede ao abrir.

Desde a S6-02, ``--sats {6,18}`` escolhe o desenho em ``CONSTELLATION_DESIGNS`` (padrão 6). Os
nomes de saída padrão são demo/index.html (6) e demo/constelacao18.html (18); com ``--days 6``,
demo/ciclo_6dias.html e demo/ciclo_6dias_18.html. ``--all`` gera os quatro arquivos, com a mesma
época, e os liga como variantes (``meta.variants``): o seletor "6 satélites | 18 satélites" do
viewer navega entre os irmãos de mesma duração. O passo vem de ``STEP_BY_SATS_DAYS`` e a precisão
de ``PRECISION_BY_SATS``, escolhidos pelo tamanho medido (cada arquivo < 15 MB). A faixa usa a
rolagem de cobertura diária do desenho (23,1957° com 6, 8,1835° com 18).

Para adicionar um satélite à demo, acrescente-o à lista montada em ``build_demo``; para mudar a
constelação, ``build_constellation``.
"""

import argparse
import base64
import json
import logging
import math
from datetime import UTC, datetime
from pathlib import Path

from satsim.app.scene_export import (
    COMPACT_PRECISION,
    FULL_PRECISION,
    SceneConstellation,
    ScenePrecision,
    SceneSatellite,
    SceneStation,
    build_scene,
    round_ms,
)
from satsim.constellation import (
    CONSTELLATION_DESIGNS,
    DESIGN_SEMI_MAJOR_AXIS_M,
    Constellation,
)
from satsim.forces import (
    DragAcceleration,
    SolarRadiationPressure,
    ThirdBodyAcceleration,
    ZonalAcceleration,
)
from satsim.maintenance import maintenance_status, measure_decay_rate
from satsim.observer import SAO_BENTO_DO_SAPUCAI
from satsim.propagators import SGP4Propagator
from satsim.tle import ISS_NORAD_ID, SENTINEL_2A_NORAD_ID, Tle, TleError, get_tle

logger = logging.getLogger("build_demo")

ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "src" / "satsim" / "app" / "static"
DEMO_DIR = ROOT / "demo"
OUTPUT = DEMO_DIR / "index.html"

DAY_S = 24 * 3600
STEP_BY_DAYS = {1: 20, 2: 20, 3: 30, 4: 40, 5: 50, 6: 60}
"""Passo da amostragem [s] por duração da janela [dias]: o arquivo fica abaixo de 15 MB (com 8
satélites, ≈ 1,2 MB por 1 000 amostras de cada; de 2 a 6 dias, 8 641 amostras ≈ 12 MB)."""
STEP_BY_SATS_DAYS: dict[tuple[int, int], int] = {
    **{(6, days): step for days, step in STEP_BY_DAYS.items()},
    **{(18, days): 20 * days for days in STEP_BY_DAYS},
}
"""Passo [s] por (P, dias) (S6-02), escolhido pelo tamanho medido (cada arquivo < 15 MB). P = 6:
o ``STEP_BY_DAYS`` da S5-06 (cenas idênticas). P = 18 (20 objetos, ``COMPACT_PRECISION``): 20 s ×
dias, 4 321 amostras por objeto, ≈ 11,5 MB em qualquer duração (6 dias a 60 s: 20,6 MB; a 90 s:
14,5 MB). Com 120 s a interpolação do viewer (Hermite cúbica no vetor unitário) erra ≈ 0,3 km."""
PRECISION_BY_SATS: dict[int, ScenePrecision] = {6: FULL_PRECISION, 18: COMPACT_PRECISION}
"""Casas decimais da cena por P: a de 6 continua com as da S5-06; a de 18, sem as casas que não
aparecem na tela (≈ 18% menos bytes)."""
OUTPUT_NAMES: dict[tuple[int, int], str] = {
    (6, 1): "index.html",
    (18, 1): "constelacao18.html",
    (6, 6): "ciclo_6dias.html",
    (18, 6): "ciclo_6dias_18.html",
}
"""Nomes de saída padrão por (P, dias) (S6-02); as outras durações usam ``ciclo_{d}dias[_18]``."""
VARIANT_LABELS = {6: "6 satélites", 18: "18 satélites"}
ALL_DAYS = (1, 6)  # --all: 1 dia e o ciclo de 6 dias, cada um com 6 e 18 satélites
EL_MIN_DEG = 10.0  # elevação mínima das passagens e do footprint
ISS_COLOR = "#ffb020"
SENTINEL_COLOR = "#c084fc"
MAX_TLE_OFFSET_DAYS = 30.0  # TLE mais distante que isso do início da cena: satélite omitido

# Constelação 91/6 (S5-01 a S5-06): Sem4 §7.5 e §8, Apêndices A.1, A.4 e A.7
CONSTELLATION_NAME = "Constelação 91/6"
A_CALIBRATED_M = 6888108.374  # a que fecha o ciclo no Cowell completo (S5-02, Δa = +19,37 m)
I_CALIBRATED_DEG = 97.46325  # A.7 da Semana 4 (Cowell J2–J4)
LTAN_H = 22.5
MAINTENANCE_WINDOW_DAYS = 1.0  # janela da medida de ȧ (S5-05)
# Cowell completo do A.7: SRP C_R 1,3, A/m 0,01 m²/kg; arrasto C_D·A/m 0,011 m²/kg
SRP_CR, SRP_AREA_OVER_MASS = 1.3, 0.01
DRAG_CD, DRAG_AREA_OVER_MASS = 2.2, 0.005


def full_force_models(epoch: datetime) -> tuple:
    """Cowell completo do A.7: J2–J4 + Sol/Lua + SRP + arrasto."""
    return (
        ZonalAcceleration(),
        ThirdBodyAcceleration("sun", epoch),
        ThirdBodyAcceleration("moon", epoch),
        SolarRadiationPressure(epoch, cr=SRP_CR, area_over_mass=SRP_AREA_OVER_MASS),
        DragAcceleration(cd=DRAG_CD, area_over_mass=DRAG_AREA_OVER_MASS),
    )


def build_constellation(start: datetime, n_sats: int = 6) -> Constellation:
    """A constelação 91/6 de ``n_sats`` satélites (``CONSTELLATION_DESIGNS``) no Cowell
    completo, época = início da cena (LTAN 22h30 nesse instante)."""
    design = CONSTELLATION_DESIGNS[n_sats]
    return design.build(start, full_force_models(start), a=A_CALIBRATED_M, ltan_h=LTAN_H,
                        i=math.radians(I_CALIBRATED_DEG))  # fmt: skip


def constellation_info(constellation: Constellation) -> SceneConstellation:
    """``SceneConstellation`` da demo: nós, famílias, faixa da rolagem de cobertura diária do
    desenho e manutenção (ȧ medido em 1 dia)."""
    design = CONSTELLATION_DESIGNS[len(constellation)]
    i = math.radians(I_CALIBRATED_DEG)
    maintenance = {
        sat.id: maintenance_status(
            sat.id, measure_decay_rate(sat, MAINTENANCE_WINDOW_DAYS), A_CALIBRATED_M, i
        )
        for sat in constellation
    }
    # rolagem de cobertura diária no a de projeto: 23,1957° (P = 6, a do S5-03) e 8,1835° (P = 18)
    return SceneConstellation.from_design(
        design, CONSTELLATION_NAME, constellation.ids,
        roll_rad=design.roll_for_daily_coverage(DESIGN_SEMI_MAJOR_AXIS_M), maintenance=maintenance,
    )  # fmt: skip


def default_output(n_sats: int, days: int, out_dir: Path = DEMO_DIR) -> Path:
    """Arquivo de saída padrão de (P, dias): ``OUTPUT_NAMES`` ou ``ciclo_{d}dias[_18].html``."""
    name = OUTPUT_NAMES.get((n_sats, days))
    if name is None:
        name = f"ciclo_{days}dias" + ("" if n_sats == 6 else f"_{n_sats}") + ".html"
    return out_dir / name


def variants_for(days: int) -> list[dict[str, str]]:
    """``meta.variants`` das cenas de ``days`` dias: uma por P, com ``href`` relativo (irmãos na
    mesma pasta)."""
    return [
        {"id": f"p{n}", "label": VARIANT_LABELS[n], "href": default_output(n, days).name}
        for n in sorted(CONSTELLATION_DESIGNS)
    ]


STATIONS = [SceneStation.from_ground_station("sbs", SAO_BENTO_DO_SAPUCAI)]


def load_tle(tle_file: Path | None, offline: bool, cache_dir: str | None) -> tuple[Tle, str]:
    """TLE da ISS de um arquivo ou de ``get_tle`` (rede → cache → embarcado) e a sua fonte."""
    if tle_file is not None:
        return Tle.from_text(tle_file.read_text(encoding="utf-8")), f"arquivo {tle_file.name}"
    result = get_tle(ISS_NORAD_ID, offline=offline, cache_dir=cache_dir)
    return result.tle, result.source


def sentinel_2a(start: datetime, offline: bool, cache_dir: str | None) -> SceneSatellite | None:
    """Sentinel-2A real (SGP4, ``get_tle(40697)``), ou ``None`` se o TLE estiver longe da cena."""
    try:
        result = get_tle(SENTINEL_2A_NORAD_ID, offline=offline, cache_dir=cache_dir)
    except TleError as exc:
        logger.warning("Sentinel-2A omitido: %s", exc)
        return None
    offset_days = abs((start - result.tle.epoch).total_seconds()) / 86400.0
    if offset_days > MAX_TLE_OFFSET_DAYS:
        logger.warning("Sentinel-2A omitido: TLE a %.0f dias do início da cena", offset_days)
        return None
    logger.info("Sentinel-2A: TLE época %s, fonte %s",
                result.tle.epoch.isoformat(timespec="seconds"), result.source)  # fmt: skip
    return SceneSatellite("s2a", "Sentinel-2A", SENTINEL_COLOR, SGP4Propagator(result.tle))


def scene_start(text: str | None, tle: Tle) -> datetime:
    """Início da cena: agora (arredondado ao minuto), ``"tle-epoch"`` ou um instante ISO 8601."""
    if text is None:
        return datetime.now(UTC).replace(second=0, microsecond=0)
    if text == "tle-epoch":
        return tle.epoch
    start = datetime.fromisoformat(text.replace("Z", "+00:00"))
    if start.tzinfo is None:
        raise SystemExit("--start precisa de fuso horário (ex.: 2026-10-23T00:00:00Z)")
    return start


def _inline_script(code: str) -> str:
    # "</script" dentro do conteúdo fecharia a tag; o escape "<\/" é equivalente em JS.
    return "<script>\n" + code.replace("</script", "<\\/script") + "\n</script>"


def render_html(scene: dict) -> str:
    """Troca as referências externas de viewer.html pelo conteúdo embutido."""
    html = (STATIC / "viewer.html").read_text(encoding="utf-8")
    css = (STATIC / "viewer.css").read_text(encoding="utf-8")
    vendor = (STATIC / "vendor" / "globe.gl.min.js").read_text(encoding="utf-8")
    viewer = (STATIC / "viewer.js").read_text(encoding="utf-8")
    texture_b64 = base64.b64encode((STATIC / "textures" / "earth_2048.jpg").read_bytes())
    texture_uri = "data:image/jpeg;base64," + texture_b64.decode("ascii")

    scene_js = (
        "window.SATSIM_SCENE = "
        + json.dumps(scene, ensure_ascii=False, separators=(",", ":"))
        + ";\nwindow.SATSIM_TEXTURE = "
        + json.dumps(texture_uri)
        + ";"
    )
    replacements = {
        '<link rel="stylesheet" href="viewer.css">': "<style>\n" + css + "\n</style>",
        '<script src="vendor/globe.gl.min.js"></script>': _inline_script(vendor),
        '<script src="scene.js"></script>': _inline_script(scene_js),
        '<script src="viewer.js"></script>': _inline_script(viewer),
    }
    for tag, content in replacements.items():
        if tag not in html:
            raise RuntimeError(f"viewer.html não contém a tag esperada: {tag}")
        html = html.replace(tag, content)
    return html


def build_demo(
    output: Path,
    tle: Tle,
    source: str,
    start: datetime,
    constellation: Constellation,
    sentinel: SceneSatellite | None,
    days: int,
    variants: list[dict[str, str]] | None = None,
) -> dict:
    """Monta a cena de ``days`` dias com a constelação dada, grava o HTML e devolve a cena."""
    n_sats = len(constellation)
    iss = SceneSatellite("iss", "ISS", ISS_COLOR, SGP4Propagator(tle))
    satellites = [iss, *constellation]
    if sentinel is not None:
        satellites.append(sentinel)
    duration_s, step_s = days * DAY_S, STEP_BY_SATS_DAYS[(n_sats, days)]
    scene = build_scene(satellites, STATIONS, start, duration_s, step_s, el_min_deg=EL_MIN_DEG,
                        tle=tle, tle_source=source, sun_geometry=True,
                        constellation_info=constellation_info(constellation), variants=variants,
                        precision=PRECISION_BY_SATS[n_sats])  # fmt: skip
    html = render_html(scene)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    n_samples = len(scene["satellites"][0]["lla"]) if scene["satellites"] else 0
    logger.info(
        "Demo gerada: %s (%.1f MB, %d satélites, %d amostras cada a %d s, %d passagens)",
        output,
        len(html.encode("utf-8")) / 1e6,
        len(scene["satellites"]),
        n_samples,
        step_s,
        len(scene["passes"]),
    )
    return scene


def build_all(
    out_dir: Path, tle: Tle, source: str, start: datetime, sentinel: SceneSatellite | None
) -> dict[tuple[int, int], Path]:
    """As quatro cenas de ``--all`` (P em ``CONSTELLATION_DESIGNS`` × ``ALL_DAYS``) em
    ``out_dir``, ligadas como variantes; devolve o arquivo de cada (P, dias). Usada também por
    scripts/build_site.py (S6-03)."""
    outputs: dict[tuple[int, int], Path] = {}
    for n_sats in sorted(CONSTELLATION_DESIGNS):
        # uma constelação por P: a cena de 6 dias reaproveita a integração da de 1 dia
        constellation = build_constellation(start, n_sats=n_sats)
        for days in ALL_DAYS:
            output = default_output(n_sats, days, out_dir)
            build_demo(output, tle, source, start, constellation, sentinel, days,
                       variants_for(days))  # fmt: skip
            outputs[(n_sats, days)] = output
    return outputs


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="arquivo de saída (padrão: OUTPUT_NAMES em --out-dir)")  # fmt: skip
    parser.add_argument("--out-dir", type=Path, default=DEMO_DIR,
                        help="pasta dos nomes padrão e de --all (padrão: demo/)")  # fmt: skip
    parser.add_argument("--tle-file", type=Path, default=None, help="TLE da ISS de um arquivo")
    parser.add_argument("--offline", action="store_true", help="não acessa a rede")
    parser.add_argument("--start", default=None, help="agora (padrão), ISO 8601 ou tle-epoch")
    parser.add_argument(
        "--cache-dir", default=None, help="padrão: $SATSIM_CACHE_DIR ou data/cache/tle"
    )
    parser.add_argument("--days", type=int, default=1, choices=sorted(STEP_BY_DAYS),
                        help="duração da janela em dias (6 = o ciclo inteiro)")  # fmt: skip
    parser.add_argument("--sats", type=int, default=6, choices=sorted(CONSTELLATION_DESIGNS),
                        help="satélites da constelação (S6-02)")  # fmt: skip
    parser.add_argument("--all", action="store_true",
                        help="as quatro cenas (6 e 18 satélites × 1 e 6 dias), ligadas como "
                             "variantes; ignora --sats e --days")  # fmt: skip
    args = parser.parse_args(argv)
    if args.all and args.output is not None:
        parser.error("--all grava em --out-dir; não use -o")
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    tle, source = load_tle(args.tle_file, args.offline, args.cache_dir)
    # arredondado ao milissegundo, como build_scene faz com a época da cena: o Cowell da constelação
    # (época = início da cena) não pode ser amostrado em t < 0
    start = round_ms(scene_start(args.start, tle))
    sentinel = sentinel_2a(start, args.offline, args.cache_dir)
    if args.all:
        build_all(args.out_dir, tle, source, start, sentinel)
        output = default_output(6, 1, args.out_dir)
    else:
        output = args.output or default_output(args.sats, args.days, args.out_dir)
        constellation = build_constellation(start, n_sats=args.sats)
        build_demo(output, tle, source, start, constellation, sentinel, args.days)
    logger.info("TLE: época %s, fonte %s, idade %.2f dias no início",
                tle.epoch.isoformat(timespec="seconds"), source,
                (start - tle.epoch).total_seconds() / 86400.0)  # fmt: skip
    logger.info("Abra com duplo clique ou: xdg-open %s", output)


if __name__ == "__main__":
    main()
