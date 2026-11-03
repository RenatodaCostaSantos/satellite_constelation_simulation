"""Gera demo/index.html: o viewer do satsim autocontido, que abre com duplo clique (file://).

Uso (na raiz do repositório, com o pacote instalado):

    python scripts/build_demo.py                      # ISS real, 24 h a partir de agora
    python scripts/build_demo.py --offline            # sem rede: cache ou TLE embarcado
    python scripts/build_demo.py --tle-file tests/data/iss_sample.tle --start tle-epoch

Desde a S3-06 a demo mostra a **ISS real**: TLE de satsim.tle.get_tle (cache com menos de 2 h
→ CelesTrak → cache antigo → TLE embarcado) ou de --tle-file, propagado com SGP4Propagator. A
janela é de 24 h a partir de agora (UTC, arredondado ao minuto), ou da época do TLE com
--start tle-epoch (determinístico, para testes e demonstração), com passo de 20 s.

Desde a S4-06:

- a **SSO-1** é a SSO de projeto no **Cowell completo** (J2–J4 + Sol/Lua + SRP + arrasto), a
  partir de ``mean_to_osculating`` dos elementos médios calibrados: i = 97,46325° (Apêndice A.7),
  e congelada (``frozen_eccentricity``), ω = 90° e Ω de ``raan_for_ltan(22,5 h)`` no início da
  cena, isto é, LTAN 22h30. A conversão usa só a gravidade zonal (``conversion_models``);
- o **Sentinel-2A real** (SGP4, TLE de ``get_tle(40697)``, com TLE embarcado para o modo
  offline) entra como terceiro satélite, desde que o TLE esteja a menos de
  ``MAX_TLE_OFFSET_DAYS`` do início da cena (com ``--start tle-epoch`` e o TLE de 2019 da ISS,
  ele fica de fora);
- a cena traz LTAN, β e eclipse de cada satélite (``build_scene(..., sun_geometry=True)``).

A cena traz também a estação de São Bento do Sapucaí, o que a ISS faz vista de lá (azimute,
elevação, distância), o círculo de visibilidade (footprint) a el_min = 10° e as passagens. Tudo
é calculado em Python e embutido no HTML junto com o CSS, o JS do viewer, o bundle UMD do
globe.gl e a textura da Terra (data URI). Nada é carregado da rede ao abrir.

Para adicionar um satélite de projeto à demo, acrescente-o em ``design_satellites``.
"""

import argparse
import base64
import json
import logging
import math
from datetime import UTC, datetime
from pathlib import Path

from satsim.app.scene_export import SceneSatellite, SceneStation, build_scene, round_ms
from satsim.astrotime import datetime_to_jd
from satsim.elements import KeplerianElements
from satsim.forces import (
    DragAcceleration,
    SolarRadiationPressure,
    ThirdBodyAcceleration,
    ZonalAcceleration,
)
from satsim.mean_osc import conversion_models, mean_to_osculating
from satsim.observer import SAO_BENTO_DO_SAPUCAI
from satsim.orbit_design import frozen_eccentricity, raan_for_ltan
from satsim.propagators import CowellPropagator, SGP4Propagator
from satsim.tle import ISS_NORAD_ID, SENTINEL_2A_NORAD_ID, Tle, TleError, get_tle

logger = logging.getLogger("build_demo")

ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "src" / "satsim" / "app" / "static"
OUTPUT = ROOT / "demo" / "index.html"

DURATION_S = 24 * 3600  # janela de 24 h
STEP_S = 20  # passo de 20 s
EL_MIN_DEG = 10.0  # elevação mínima das passagens e do footprint
ISS_COLOR = "#ffb020"
SENTINEL_COLOR = "#c084fc"
MAX_TLE_OFFSET_DAYS = 30.0  # TLE mais distante que isso do início da cena: satélite omitido

# SSO de projeto 91/6 (Sem4 §7.5) calibrada no Cowell J2–J4 (S4-05, Apêndice A.7)
SSO_A_M = 6888.089e3
SSO_I_DEG = 97.46325
SSO_LTAN_H = 22.5
SSO_COLOR = "#00e5ff"
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


def sso_design(start: datetime) -> SceneSatellite:
    """SSO-1: elementos médios calibrados → osculador → Cowell completo, época = início da cena."""
    i = math.radians(SSO_I_DEG)
    mean = KeplerianElements(
        a=SSO_A_M,
        e=frozen_eccentricity(SSO_A_M, i),
        i=i,
        raan=float(raan_for_ltan(SSO_LTAN_H, datetime_to_jd(start))),
        argp=math.pi / 2.0,
        psi=0.0,
    )
    models = full_force_models(start)
    osc = mean_to_osculating(mean, start, conversion_models(models)).elements
    return SceneSatellite("sso1", "SSO-1 (Cowell)", SSO_COLOR, CowellPropagator(osc, start, models))


def design_satellites(start: datetime) -> list[SceneSatellite]:
    """Satélites de projeto da demo (todos com época no início da cena)."""
    return [sso_design(start)]


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
    # arredondado ao milissegundo, como build_scene faz com a época da cena: o Cowell da SSO-1
    # (época = início da cena) não pode ser amostrado em t < 0
    start = round_ms(scene_start(args.start, tle))
    iss = SceneSatellite("iss", "ISS", ISS_COLOR, SGP4Propagator(tle))
    satellites = [iss, *design_satellites(start)]
    sentinel = sentinel_2a(start, args.offline, args.cache_dir)
    if sentinel is not None:
        satellites.append(sentinel)
    scene = build_scene(satellites, STATIONS, start, DURATION_S, STEP_S, el_min_deg=EL_MIN_DEG,
                        tle=tle, tle_source=source, sun_geometry=True)  # fmt: skip
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
