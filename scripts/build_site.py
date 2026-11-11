"""Monta site/: o simulador como site estático (página inicial + as quatro cenas), para o GitHub
Pages ou qualquer servidor de arquivos (S6-03).

Uso (na raiz do repositório, com o pacote instalado):

    python scripts/build_site.py --offline          # sem rede: cache ou TLEs embutidos
    python scripts/build_site.py                    # TLEs da CelesTrak (fallback: embutidos)
    python -m http.server -d site                   # e abra http://localhost:8000/

Estrutura gerada (``site/`` é ignorado pelo git; o workflow .github/workflows/pages.yml publica):

    site/index.html                 página inicial leve (sem cena), em linguagem simples
    site/img/demo_constelacao18.png imagem da demo (docs/img)
    site/cenas/index.html           6 satélites · 1 dia       ┐ geradas por build_demo.build_all
    site/cenas/constelacao18.html   18 satélites · 1 dia      │ (o mesmo que build_demo.py --all),
    site/cenas/ciclo_6dias.html     6 satélites · 6 dias      │ ligadas pelo seletor 6/18 do
    site/cenas/ciclo_6dias_18.html  18 satélites · 6 dias     ┘ viewer (hrefs relativos)

Todos os links são relativos (o site funciona em https://<usuário>.github.io/<repositório>/, num
http.server local ou em file://) e nada é carregado da rede: as cenas embutem CSS, JS, globe.gl e
a textura, como a demo; a página inicial só usa CSS embutido e fontes do sistema.

Os TLEs da ISS e do Sentinel-2A vêm de ``satsim.tle.get_tle`` (CelesTrak → cache → embutidos no
pacote). Se a CelesTrak falhar (por exemplo, no CI), o site é montado com o fallback e a página
inicial avisa de onde vieram os TLEs e a idade deles.
"""

import argparse
import html
import logging
import shutil
import sys
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_demo as demo  # noqa: E402

from satsim.app.scene_export import SceneSatellite, round_ms  # noqa: E402
from satsim.tle import Tle  # noqa: E402

logger = logging.getLogger("build_site")

ROOT = Path(__file__).resolve().parents[1]
SITE_DIR = ROOT / "site"
SCENES_SUBDIR = "cenas"
IMAGE = ROOT / "docs" / "img" / "demo_constelacao18.png"
IMAGE_SUBPATH = "img/demo_constelacao18.png"
MAX_HOME_BYTES = 1_000_000  # página inicial < 1 MB (critério do S6-03)
MAX_SCENE_BYTES = 15_000_000  # cada cena < 15 MB (S6-02)

BUTTONS: tuple[tuple[int, int, str], ...] = (
    (6, 1, "6 satélites · 1 dia"),
    (18, 1, "18 satélites · 1 dia"),
    (6, 6, "6 satélites · ciclo de 6 dias"),
    (18, 6, "18 satélites · ciclo de 6 dias"),
)
"""Os quatro botões da página inicial: (P, dias, rótulo), na ordem da grade 2 × 2."""

SOURCE_LABELS = {
    "network": "CelesTrak, baixados na montagem do site",
    "cache": "cache local",
    "bundled": "TLEs embutidos no pacote",
}

_CSS = """
:root { --bg: #0b1020; --panel: #121a30; --ink: #e6ecf5; --muted: #9aa7bd;
        --accent: #4cc3ff; --warn: #ffb020; --line: #24304d; }
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink); line-height: 1.55;
       font: 16px/1.55 system-ui, -apple-system, "Segoe UI", Roboto, "DejaVu Sans", sans-serif; }
main { max-width: 860px; margin: 0 auto; padding: 32px 16px 48px; }
h1 { font-size: 2rem; margin: 0 0 4px; letter-spacing: .02em; }
h1 span { color: var(--muted); font-weight: 400; font-size: 1.1rem; }
h2 { font-size: 1.15rem; margin: 32px 0 8px; color: var(--accent); }
p { margin: 8px 0; }
.lead { font-size: 1.1rem; color: var(--ink); }
.warn { border-left: 4px solid var(--warn); background: #2a2210; padding: 10px 14px;
        margin: 20px 0; border-radius: 4px; }
.launch { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 12px;
          margin: 24px 0; }
.launch a { display: block; padding: 16px; border: 1px solid var(--line); border-radius: 8px;
            background: var(--panel); color: var(--ink); text-decoration: none; }
.launch a:hover, .launch a:focus { border-color: var(--accent); outline: none; }
.launch b { display: block; font-size: 1.05rem; }
.launch small { color: var(--muted); }
figure { margin: 24px 0; }
figure img { width: 100%; height: auto; border: 1px solid var(--line); border-radius: 6px; }
figcaption { color: var(--muted); font-size: .9rem; margin-top: 6px; }
ul { padding-left: 20px; }
footer { color: var(--muted); font-size: .9rem; border-top: 1px solid var(--line);
         margin-top: 40px; padding-top: 12px; }
@media (max-width: 560px) { .launch { grid-template-columns: 1fr; } }
"""


def _fmt_mb(n_bytes: int) -> str:
    if n_bytes < 100_000:
        return f"{n_bytes / 1e3:.0f} kB"
    return f"{n_bytes / 1e6:.1f} MB".replace(".", ",")


def _fmt_utc(t: datetime) -> str:
    return f"{t:%d/%m/%Y %H:%M} UTC"


def render_home(
    tle: Tle,
    source: str,
    start: datetime,
    sentinel_included: bool,
    scene_bytes: dict[tuple[int, int], int],
    scene_hrefs: dict[tuple[int, int], str],
    image_href: str | None,
) -> str:
    """HTML da página inicial: texto para leigos, os quatro botões, a imagem e a origem dos TLEs.

    ``scene_hrefs`` são caminhos relativos à página (``cenas/…``); ``scene_bytes`` vai no botão
    para quem tem conexão lenta. Se ``source`` não for ``"network"``, um aviso diz de onde os
    TLEs vieram."""
    e = html.escape
    age_days = (start - tle.epoch).total_seconds() / 86400.0
    buttons = "\n".join(
        f'  <a href="{e(scene_hrefs[(n, d)])}"><b>{e(label)}</b>'
        f"<small>{_fmt_mb(scene_bytes[(n, d)])} · abre no navegador, sem instalar nada</small></a>"
        for n, d, label in BUTTONS
    )
    warn = ""
    if source != "network":
        warn = (
            '<div class="warn" role="note"><b>Aviso:</b> os TLEs desta versão não foram baixados '
            "da CelesTrak na montagem do site (fonte: "
            f"{e(SOURCE_LABELS.get(source, source))}). O TLE da ISS tem época "
            f"{e(_fmt_utc(tle.epoch))}, {age_days:.0f} dias antes do início da cena: a posição "
            "da ISS (e do Sentinel-2A) fica menos precisa quanto mais velho o TLE. A constelação "
            "91/6 não depende de TLE.</div>"
        )
    figure = ""
    if image_href is not None:
        figure = (
            f'<figure><img src="{e(image_href)}" width="1442" height="454" loading="lazy" '
            'alt="Mapas 2D da demo no mesmo instante: com 6 satélites, 61 de 91 trilhas da '
            'grade percorridas; com 18, 182 de 273."><figcaption>O mesmo instante com 6 e 18 '
            "satélites: os traços no equador são as trilhas da grade do dia.</figcaption></figure>"
        )
    sentinel_text = (
        "Os satélites reais ISS e Sentinel-2A também estão na cena (desligados na legenda)."
        if sentinel_included
        else "A ISS real também está na cena (desligada na legenda); o Sentinel-2A ficou de fora "
        "porque o TLE disponível estava longe da data da cena."
    )
    end_1d = start + timedelta(days=1)
    return f"""<!doctype html>
<html lang="pt-BR">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>satsim — simulador de satélites</title>
<link rel="icon" href="data:,">
<meta name="description" content="Constelação de satélites em órbita 91/6 no navegador: 6 ou 18
satélites, globo 3D e mapa 2D.">
<style>{_CSS}</style>
</head>
<body>
<main>
<h1>satsim <span>simulador de satélites</span></h1>
<p class="lead">Veja no navegador uma constelação de satélites de observação da Terra girando
sobre o globo, com mapa 2D, relógio e controles de tempo. Escolha uma das quatro versões:</p>
{warn}
<nav class="launch" aria-label="Abrir o simulador">
{buttons}
</nav>
<p>Dentro do simulador, o seletor <b>6 satélites | 18 satélites</b> na barra superior troca de
versão no mesmo instante. Arraste o globo para girar, use a linha do tempo embaixo e clique num
satélite para ver os dados dele.</p>
{figure}
<h2>O que é</h2>
<p>O satsim é um simulador escrito em Python. Ele calcula onde cada satélite está, minuto a
minuto, com a gravidade da Terra achatada, a atração do Sol e da Lua, a pressão da luz do Sol e o
atrito com o resto da atmosfera. O resultado é guardado nestas páginas: o navegador só desenha,
sem servidor e sem calcular nada. {sentinel_text}</p>
<h2>A órbita 91/6</h2>
<p>Os satélites voam a cerca de 510 km de altitude, passando perto dos polos, sempre na mesma
hora solar (órbita heliossíncrona). A altitude foi escolhida para que cada satélite dê
exatamente 91 voltas na Terra em 6 dias e então repita o mesmo caminho sobre o solo. Todos ficam
no mesmo plano, igualmente espaçados, como contas num colar.</p>
<h2>6 ou 18 satélites</h2>
<ul>
<li><b>6 satélites</b>, a 60° um do outro: juntos percorrem as 91 trilhas do ciclo em um único
dia, a cerca de 440 km uma da outra no equador. Para ver qualquer ponto todo dia, a câmera
precisa inclinar até 23,2° para os lados.</li>
<li><b>18 satélites</b>, a 20° um do outro: são 273 trilhas por dia, a cerca de 147 km. Basta
inclinar a câmera 8,2°, então as imagens saem mais de frente.</li>
<li>O intervalo entre visitas a um mesmo lugar quase não muda (cerca de 14 h nos dois casos):
os satélites vizinhos passam com poucos minutos de diferença, em "trens".</li>
</ul>
<p>As versões de 1 dia são mais leves e detalhadas; as de 6 dias mostram o ciclo inteiro até o
caminho se repetir.</p>
<footer>
<p>Cenas de {e(_fmt_utc(start))} a {e(_fmt_utc(end_1d))} (versões de 1 dia) ou 6 dias a partir
do mesmo instante. TLE da ISS: época {e(_fmt_utc(tle.epoch))}, fonte:
{e(SOURCE_LABELS.get(source, source))}. TLEs de CelesTrak (celestrak.org), textura da Terra e
créditos em docs/creditos.md do repositório.</p>
<p>Gerado por scripts/build_site.py.</p>
</footer>
</main>
</body>
</html>
"""


def build_site(
    out_dir: Path,
    tle: Tle,
    source: str,
    start: datetime,
    sentinel: SceneSatellite | None,
    image: Path | None = IMAGE,
) -> dict[str, int]:
    """Monta ``out_dir`` (cenas em ``out_dir/cenas``, imagem e página inicial) e devolve o
    tamanho em bytes de cada arquivo, por caminho relativo."""
    scenes = demo.build_all(out_dir / SCENES_SUBDIR, tle, source, start, sentinel)
    hrefs = {key: f"{SCENES_SUBDIR}/{path.name}" for key, path in scenes.items()}
    scene_bytes = {key: path.stat().st_size for key, path in scenes.items()}
    image_href = None
    if image is not None and image.exists():
        (out_dir / IMAGE_SUBPATH).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(image, out_dir / IMAGE_SUBPATH)
        image_href = IMAGE_SUBPATH
    home = render_home(tle, source, start, sentinel is not None, scene_bytes, hrefs, image_href)
    (out_dir / "index.html").write_text(home, encoding="utf-8")
    files = ["index.html", *([image_href] if image_href else []), *hrefs.values()]
    return {name: (out_dir / name).stat().st_size for name in files}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-o", "--output", type=Path, default=SITE_DIR,
                        help="pasta do site (padrão: site/)")  # fmt: skip
    parser.add_argument("--offline", action="store_true", help="não acessa a rede")
    parser.add_argument("--start", default=None, help="agora (padrão), ISO 8601 ou tle-epoch")
    parser.add_argument(
        "--cache-dir", default=None, help="padrão: $SATSIM_CACHE_DIR ou data/cache/tle"
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    tle, source = demo.load_tle(None, args.offline, args.cache_dir)
    start = round_ms(demo.scene_start(args.start, tle))
    sentinel = demo.sentinel_2a(start, args.offline, args.cache_dir)
    sizes = build_site(args.output, tle, source, start, sentinel)

    print(f"Site em {args.output}:")
    for name, n_bytes in sizes.items():
        print(f"  {name:<34} {_fmt_mb(n_bytes):>9}")
    print(f"  {'total':<34} {_fmt_mb(sum(sizes.values())):>9}")
    if sizes["index.html"] >= MAX_HOME_BYTES:
        raise SystemExit("página inicial com 1 MB ou mais")
    too_big = [n for n in sizes if n.startswith(SCENES_SUBDIR) and sizes[n] >= MAX_SCENE_BYTES]
    if too_big:
        raise SystemExit(f"cenas com 15 MB ou mais: {', '.join(too_big)}")
    if source != "network":
        logger.warning("TLE da ISS de fonte %r: a página inicial mostra o aviso", source)


if __name__ == "__main__":
    main()
