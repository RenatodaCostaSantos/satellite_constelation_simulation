"""S6-03: site estático (scripts/build_site.py) e workflow do GitHub Pages.

Tudo offline (rede bloqueada por monkeypatch), em tmp_path, com a constelação barata de
test_scene_export (o tamanho e os links não dependem do propagador). O workflow é verificado por
texto: o ambiente não tem PyYAML e o S6-03 não adiciona dependência.
"""

import contextlib
import importlib.util
import io
import re
import urllib.error
import urllib.request
from datetime import UTC, datetime
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest
from golden.reference_values import DEMO_SCENE_S5
from test_scene_export import _fast_constellation, _scene_from_html

from satsim.tle import ISS_NORAD_ID, load_bundled_tle

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "pages.yml"
SCENES = ("cenas/index.html", "cenas/constelacao18.html", "cenas/ciclo_6dias.html",
          "cenas/ciclo_6dias_18.html")  # fmt: skip
MAX_SCENE_MB = 15.0
MAX_HOME_MB = 1.0


def _load_build_site():
    spec = importlib.util.spec_from_file_location("build_site", ROOT / "scripts" / "build_site.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class _Refs(HTMLParser):
    """src/href de todas as tags; o conteúdo de <script> e <style> não é lido como HTML."""

    def __init__(self) -> None:
        super().__init__()
        self.refs: list[tuple[str, str, str]] = []

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ("src", "href", "srcset", "poster", "action") and value is not None:
                self.refs.append((tag, key, value))


def _refs(path: Path) -> list[tuple[str, str, str]]:
    parser = _Refs()
    parser.feed(path.read_text(encoding="utf-8"))
    return parser.refs


@pytest.fixture(scope="module")
def site(tmp_path_factory):
    """site/ montado uma vez por ``main --offline`` (rede bloqueada), com a saída do terminal."""
    tmp = tmp_path_factory.mktemp("site_s6")
    out = tmp / "site"

    def no_network(*args, **kwargs):
        raise urllib.error.URLError("sem rede (teste)")

    module = _load_build_site()
    stdout = io.StringIO()
    with pytest.MonkeyPatch.context() as mp, contextlib.redirect_stdout(stdout):
        mp.setattr(urllib.request, "urlopen", no_network)
        mp.setattr(module.demo, "build_constellation", _fast_constellation)
        module.main(["--offline", "--cache-dir", str(tmp / "cache"), "-o", str(out),
                     "--start", DEMO_SCENE_S5.start_utc])  # fmt: skip
    return out, stdout.getvalue(), module


def test_s6_site_has_home_and_four_scenes_under_size_limits(site) -> None:
    out, printed, _ = site
    assert (out / "index.html").is_file()
    for name in SCENES:
        assert (out / name).is_file(), name
        assert (out / name).stat().st_size / 1e6 < MAX_SCENE_MB, name
    assert (out / "index.html").stat().st_size / 1e6 < MAX_HOME_MB
    # a página inicial é leve: não embute cena, viewer nem textura
    home = (out / "index.html").read_text(encoding="utf-8")
    assert "SATSIM_SCENE" not in home and "<script" not in home and "base64" not in home
    # imprime o tamanho de cada arquivo e o total
    for name in ("index.html", *SCENES):
        assert re.search(rf"^\s+{re.escape(name)}\s+\d+(,\d)? [kM]B$", printed, re.M), name
    assert re.search(r"^\s+total\s+\d+,\d MB$", printed, re.M)


def test_s6_site_links_are_relative_and_resolve(site) -> None:
    out, _, _ = site
    pages = sorted(out.rglob("*.html"))
    assert len(pages) == 5
    for page in pages:
        refs = _refs(page)
        for tag, key, value in refs:
            parts = urlsplit(value)
            if value.startswith("#") or parts.scheme == "data":
                continue  # âncora na própria página ou recurso embutido
            # nada absoluto: http(s)://, //host ou /raiz (quebraria em <usuário>.github.io/<repo>/)
            assert not parts.scheme and not parts.netloc, (page.name, tag, key, value)
            assert not value.startswith("/"), (page.name, tag, key, value)
            target = (page.parent / unquote(parts.path)).resolve()
            assert target.is_file(), (page.name, tag, key, value)
            assert out.resolve() in target.parents, (page.name, value)  # não sai do site
    # as cenas não têm tags com src/href externos (tudo embutido); as variantes do seletor
    # (meta.variants) são irmãos relativos na mesma pasta
    for name in SCENES:
        scene_page = out / name
        assert _refs(scene_page) == [], name
        meta = _scene_from_html(scene_page.read_text(encoding="utf-8"))["meta"]
        hrefs = [v["href"] for v in meta["variants"]]
        assert len(hrefs) == 2
        for href in hrefs:
            assert "/" not in href and ":" not in href
            assert (scene_page.parent / href).is_file()


def test_s6_home_page_content(site) -> None:
    out, _, module = site
    home = (out / "index.html").read_text(encoding="utf-8")
    links = [v for tag, key, v in _refs(out / "index.html") if tag == "a" and key == "href"]
    assert links[:4] == list(SCENES)  # quatro botões: 6/18 × 1 dia/ciclo de 6 dias
    for label in ("6 satélites · 1 dia", "18 satélites · 1 dia", "6 satélites · ciclo de 6 dias",
                  "18 satélites · ciclo de 6 dias"):  # fmt: skip
        assert label in home
    assert 'src="img/demo_constelacao18.png"' in home and (out / module.IMAGE_SUBPATH).is_file()
    for text in ("91 voltas", "6 dias", "273 trilhas", "440 km", "147 km", "23,2°", "8,2°"):
        assert text in home, text
    # data e fonte dos TLEs; offline sem cache = TLE embutido, com aviso
    tle = load_bundled_tle(ISS_NORAD_ID)
    assert f"{tle.epoch:%d/%m/%Y %H:%M} UTC" in home
    assert 'class="warn"' in home and module.SOURCE_LABELS["bundled"] in home
    assert "23/10/2026 00:00 UTC" in home  # início das cenas (--start)


def test_s6_home_page_without_warning_when_tles_come_from_network(site) -> None:
    _, _, module = site
    tle = load_bundled_tle(ISS_NORAD_ID)
    keys = [(n, d) for n, d, _ in module.BUTTONS]
    start = datetime(2026, 10, 23, tzinfo=UTC)
    args = (tle, "network", start, True, dict.fromkeys(keys, 10**7),
            {k: f"cenas/{k[0]}_{k[1]}.html" for k in keys}, None)  # fmt: skip
    home = module.render_home(*args)
    assert 'class="warn"' not in home and module.SOURCE_LABELS["network"] in home
    assert "<figure>" not in home  # sem imagem, sem figura (nenhum src quebrado)
    assert "Sentinel-2A também" in home
    bundled = module.render_home(tle, "bundled", *args[2:])
    assert 'class="warn"' in bundled


def test_s6_pages_workflow_by_text() -> None:
    text = WORKFLOW.read_text(encoding="utf-8")
    lines = [ln.split("#", 1)[0].rstrip() for ln in text.splitlines()]
    body = "\n".join(lines)
    # gatilhos: push na main e manual
    assert re.search(r"^on:\n  push:\n    branches: \[main\]\n  workflow_dispatch:", body, re.M)
    # dois jobs; o deploy depende dos testes
    assert re.search(r"^  test:$", body, re.M) and re.search(r"^  deploy:$", body, re.M)
    test_job, deploy_job = body.split("\n  deploy:\n")
    assert re.search(r"^    needs: test$", deploy_job, re.M)
    # job test: Python 3.11, instala com [dev], ruff e pytest
    for step in ('python-version: "3.11"', 'pip install -e ".[dev]"', "run: ruff check .",
                 "run: ruff format --check .", "run: pytest"):  # fmt: skip
        assert step in test_job, step
    # job deploy: build do site com rede (sem --offline: fallback no get_tle) e Pages
    assert "python scripts/build_site.py -o site" in deploy_job
    assert "--offline" not in deploy_job
    for action in ("actions/upload-pages-artifact@", "actions/deploy-pages@"):
        assert action in deploy_job
    assert re.search(r"^          path: site$", deploy_job, re.M)
    assert re.search(r"^      pages: write$", deploy_job, re.M)
    assert re.search(r"^      id-token: write$", deploy_job, re.M)
    # só actions oficiais e nenhum segredo
    uses = re.findall(r"^\s*(?:- )?uses: (\S+)$", body, re.M)
    assert uses and all(u.startswith("actions/") for u in uses), uses
    assert "secrets." not in text
