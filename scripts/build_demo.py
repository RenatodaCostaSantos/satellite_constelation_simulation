"""Gera demo/index.html: o viewer do satsim autocontido, que abre com duplo clique (file://).

Uso (na raiz do repositório, com o pacote instalado):

    python scripts/build_demo.py

A cena é calculada em Python (KeplerPropagator) e embutida no HTML junto com o CSS, o JS do
viewer, o bundle UMD do globe.gl e a textura da Terra (data URI). Nada é carregado da rede.

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
from satsim.propagators import KeplerPropagator

logger = logging.getLogger("build_demo")

ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "src" / "satsim" / "app" / "static"
OUTPUT = ROOT / "demo" / "index.html"

EPOCH = datetime(2026, 9, 29, 12, 0, 0, tzinfo=UTC)
DURATION_S = 6 * 3600  # janela de 6 h
STEP_S = 10  # passo de 10 s


def kepler_sat(
    sat_id: str,
    name: str,
    color: str,
    *,
    a_km: float,
    e: float,
    i_deg: float,
    raan_deg: float,
    argp_deg: float,
    psi_deg: float,
) -> SceneSatellite:
    """Satélite kepleriano com elementos em km e graus, época = início da cena."""
    elements = KeplerianElements(
        a=a_km * 1e3,
        e=e,
        i=math.radians(i_deg),
        raan=math.radians(raan_deg),
        argp=math.radians(argp_deg),
        psi=math.radians(psi_deg),
    )
    return SceneSatellite(sat_id, name, color, KeplerPropagator(elements, EPOCH))


# fmt: off
SATELLITES = [
    # SSO de projeto 91/6 (Sem4 §7.5): a = 6888,089 km, e = 1,074e-3, i = 97,4396°, ω = 90°
    kepler_sat("sso1", "SSO-1", "#00e5ff", a_km=6888.089, e=1.074e-3, i_deg=97.4396, raan_deg=30.0, argp_deg=90.0, psi_deg=0.0),  # noqa: E501
    # Só para provar o multi-satélite: NÃO é a ISS real (a ISS via TLE/SGP4 entra na Semana 3)
    kepler_sat("kep1", "Exemplo kepleriano (i = 51,64°)", "#ffb020", a_km=6798.0, e=0.0, i_deg=51.64, raan_deg=120.0, argp_deg=0.0, psi_deg=0.0),  # noqa: E501
]
# fmt: on

STATIONS = [SceneStation("sbs", "São Bento do Sapucaí", -22.69, -45.73, 900.0)]


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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-o", "--output", type=Path, default=OUTPUT, help="arquivo de saída")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    scene = build_scene(SATELLITES, STATIONS, EPOCH, DURATION_S, STEP_S)
    html = render_html(scene)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html, encoding="utf-8")
    n_samples = len(scene["satellites"][0]["lla"]) if scene["satellites"] else 0
    logger.info(
        "Demo gerada: %s (%.1f MB, %d satélites, %d amostras cada)",
        args.output,
        len(html.encode("utf-8")) / 1e6,
        len(scene["satellites"]),
        n_samples,
    )
    logger.info("Abra com duplo clique ou: xdg-open %s", args.output)


if __name__ == "__main__":
    main()
