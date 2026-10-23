"""Gera demo/index.html: o viewer do satsim autocontido, que abre com duplo clique (file://).

Uso (na raiz do repositório, com o pacote instalado):

    python scripts/build_demo.py                      # ISS real, 24 h a partir de agora
    python scripts/build_demo.py --offline            # sem rede: cache ou TLE embarcado
    python scripts/build_demo.py --tle-file tests/data/iss_sample.tle --start tle-epoch

Desde a S3-06 a demo mostra a **ISS real**: TLE de satsim.tle.get_tle (cache com menos de 2 h
→ CelesTrak → cache antigo → TLE embarcado) ou de --tle-file, propagado com SGP4Propagator. A
janela é de 24 h a partir de agora (UTC, arredondado ao minuto), ou da época do TLE com
--start tle-epoch (determinístico, para testes e demonstração), com passo de 20 s. A SSO de
projeto (MeanJ2) segue como segundo satélite.

A cena traz também a estação de São Bento do Sapucaí, o que a ISS faz vista de lá (azimute,
elevação, distância), o círculo de visibilidade (footprint) a el_min = 10° e as passagens. Tudo
é calculado em Python e embutido no HTML junto com o CSS, o JS do viewer, o bundle UMD do
globe.gl e a textura da Terra (data URI). Nada é carregado da rede ao abrir.

Para adicionar um satélite à demo, acrescente UMA linha à lista SATELLITES abaixo.
"""

import argparse
import base64
import json
import logging
import math
from datetime import UTC, datetime
from pathlib import Path

from satsim.app.scene_export import SceneSatellite, SceneStation, build_scene
from satsim.elements import KeplerianElements
from satsim.observer import SAO_BENTO_DO_SAPUCAI
from satsim.propagators import MeanJ2Propagator, Propagator, SGP4Propagator
from satsim.tle import ISS_NORAD_ID, Tle, get_tle

logger = logging.getLogger("build_demo")

ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "src" / "satsim" / "app" / "static"
OUTPUT = ROOT / "demo" / "index.html"

DESIGN_EPOCH = datetime(2026, 9, 29, 12, 0, 0, tzinfo=UTC)  # época dos satélites de projeto
DURATION_S = 24 * 3600  # janela de 24 h
STEP_S = 20  # passo de 20 s
EL_MIN_DEG = 10.0  # elevação mínima das passagens e do footprint
ISS_COLOR = "#ffb020"


def sat(
    sat_id: str,
    name: str,
    color: str,
    propagator: type[Propagator],
    *,
    a_km: float,
    e: float,
    i_deg: float,
    raan_deg: float,
    argp_deg: float,
    psi_deg: float,
) -> SceneSatellite:
    """Satélite de projeto com elementos em km e graus, época ``DESIGN_EPOCH``.

    ``propagator`` é a classe do propagador (ex.: KeplerPropagator ou MeanJ2Propagator).
    """
    elements = KeplerianElements(
        a=a_km * 1e3,
        e=e,
        i=math.radians(i_deg),
        raan=math.radians(raan_deg),
        argp=math.radians(argp_deg),
        psi=math.radians(psi_deg),
    )
    return SceneSatellite(sat_id, name, color, propagator(elements, DESIGN_EPOCH))


# A ISS real (TLE/SGP4) entra primeiro, em main(); aqui ficam os satélites de projeto.
# fmt: off
SATELLITES = [
    # SSO de projeto 91/6 (Sem4 §7.5): a = 6888,089 km, e = 1,074e-3, i = 97,4396°, ω = 90°
    sat("sso1", "SSO-1 (J2)", "#00e5ff", MeanJ2Propagator, a_km=6888.089, e=1.074e-3, i_deg=97.4396, raan_deg=30.0, argp_deg=90.0, psi_deg=0.0),  # noqa: E501
]
# fmt: on

STATIONS = [SceneStation.from_ground_station("sbs", SAO_BENTO_DO_SAPUCAI)]


def load_tle(tle_file: Path | None, offline: bool, cache_dir: str | None) -> tuple[Tle, str]:
    """TLE da ISS de um arquivo ou de ``get_tle`` (rede → cache → embarcado) e a sua fonte."""
    if tle_file is not None:
        return Tle.from_text(tle_file.read_text(encoding="utf-8")), f"arquivo {tle_file.name}"
    result = get_tle(ISS_NORAD_ID, offline=offline, cache_dir=cache_dir)
    return result.tle, result.source


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


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-o", "--output", type=Path, default=OUTPUT, help="arquivo de saída")
    parser.add_argument("--tle-file", type=Path, default=None, help="TLE da ISS de um arquivo")
    parser.add_argument("--offline", action="store_true", help="não acessa a rede")
    parser.add_argument("--start", default=None, help="agora (padrão), ISO 8601 ou tle-epoch")
    parser.add_argument(
        "--cache-dir", default=None, help="padrão: $SATSIM_CACHE_DIR ou data/cache/tle"
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    tle, source = load_tle(args.tle_file, args.offline, args.cache_dir)
    start = scene_start(args.start, tle)
    iss = SceneSatellite("iss", "ISS", ISS_COLOR, SGP4Propagator(tle))
    scene = build_scene([iss, *SATELLITES], STATIONS, start, DURATION_S, STEP_S,
                        el_min_deg=EL_MIN_DEG, tle=tle, tle_source=source)  # fmt: skip
    html = render_html(scene)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html, encoding="utf-8")
    n_samples = len(scene["satellites"][0]["lla"]) if scene["satellites"] else 0
    logger.info(
        "Demo gerada: %s (%.1f MB, %d satélites, %d amostras cada, %d passagens)",
        args.output,
        len(html.encode("utf-8")) / 1e6,
        len(scene["satellites"]),
        n_samples,
        len(scene["passes"]),
    )
    logger.info("TLE: época %s, fonte %s, idade %.2f dias no início",
                tle.epoch.isoformat(timespec="seconds"), source,
                scene["meta"]["tle"]["age_days_at_start"])  # fmt: skip
    logger.info("Abra com duplo clique ou: xdg-open %s", args.output)


if __name__ == "__main__":
    main()
