"""Desenha o traço no solo (ground track) de uma órbita sobre o mapa da Terra e salva em PNG.

Uso (na raiz do repositório, com o pacote instalado com o extra dev):

    python scripts/plot_ground_track.py --preset sso
    python scripts/plot_ground_track.py --preset iss-like --duration-h 12
    python scripts/plot_ground_track.py --elements 7200 0.01 45 30 0 0 --out /tmp/traco.png

O subponto é calculado pela biblioteca (``satsim.groundtrack``: ECI → ECEF pelo GMST →
geodésico WGS84) e quebrado no antimeridiano, para que nenhuma linha atravesse o mapa. O fundo é
a textura do viewer (src/satsim/app/static/textures/), em projeção equiretangular.
"""

import argparse
import logging
import math
from datetime import UTC, datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # sem janela: só gera o arquivo

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from satsim.constants import R_EARTH  # noqa: E402
from satsim.elements import KeplerianElements  # noqa: E402
from satsim.groundtrack import (  # noqa: E402
    ascending_node_crossings,
    split_at_antimeridian,
    subpoint,
)
from satsim.propagators import (  # noqa: E402
    CowellPropagator,
    KeplerPropagator,
    MeanJ2Propagator,
    Propagator,
)

logger = logging.getLogger("plot_ground_track")

ROOT = Path(__file__).resolve().parents[1]
TEXTURE = ROOT / "src" / "satsim" / "app" / "static" / "textures" / "earth_2048.jpg"
REPORTS = ROOT / "reports" / "s2_comparacao"

EPOCH = datetime(2026, 9, 29, 12, 0, 0, tzinfo=UTC)
STATION = ("São Bento do Sapucaí", -22.69, -45.73)  # lat°, lon° (Apêndice A.3)

# Elementos (a [km], e, i [°], Ω0 [°], ω0 [°], ψ0 [°]); Ω0 = 30° e ψ0 = 0 como no Apêndice A.4.
PRESETS: dict[str, tuple[str, tuple[float, float, float, float, float, float]]] = {
    "sso": ("SSO de projeto 91/6", (6888.089, 1.074e-3, 97.4396, 30.0, 90.0, 0.0)),
    "iss-like": ("ISS-like (a = R⊕ + 420 km)", (R_EARTH / 1e3 + 420.0, 0.0, 51.64, 30.0, 0.0, 0.0)),
}

PROPAGATORS = {"mean-j2": MeanJ2Propagator, "kepler": KeplerPropagator, "cowell": CowellPropagator}


def make_propagator(
    kind: str, a_km: float, e: float, i_deg: float, raan_deg: float, argp_deg: float, psi_deg: float
) -> Propagator:
    """Propagador do tipo pedido, com elementos em km e graus e época EPOCH."""
    elements = KeplerianElements(
        a=a_km * 1e3,
        e=e,
        i=math.radians(i_deg),
        raan=math.radians(raan_deg),
        argp=math.radians(argp_deg),
        psi=math.radians(psi_deg),
    )
    return PROPAGATORS[kind](elements, EPOCH)


def plot_ground_track(prop: Propagator, duration_s: float, step_s: float, title: str, out: Path):
    """Calcula o subponto, desenha o traço sobre o mapa e salva em ``out``."""
    t = np.arange(0.0, duration_s + step_s / 2, step_s)
    lat, lon, _ = subpoint(prop, t)
    _, lon_nodes = ascending_node_crossings(t, lat, lon)
    segments = split_at_antimeridian(lon, lat)

    fig, ax = plt.subplots(figsize=(12, 6.4), dpi=110)
    if TEXTURE.exists():
        ax.imshow(plt.imread(TEXTURE), extent=(-180, 180, -90, 90), aspect="auto", alpha=0.9)
    for k, (seg_lon, seg_lat) in enumerate(segments):
        ax.plot(np.degrees(seg_lon), np.degrees(seg_lat), color="#ffd166", lw=1.2,
                label="traço no solo" if k == 0 else None)  # fmt: skip
    ax.scatter(np.degrees(lon_nodes), np.zeros_like(lon_nodes), s=18, color="#ef476f", zorder=3,
               label="nós ascendentes")  # fmt: skip
    ax.plot(np.degrees(lon[0]), np.degrees(lat[0]), "o", color="#06d6a0", ms=7, label="início")
    ax.plot(STATION[2], STATION[1], "*", color="#ff70d0", ms=12, mec="black", label=STATION[0])

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(np.arange(-180, 181, 30))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.grid(color="white", alpha=0.35, lw=0.6)
    ax.set_xlabel("longitude [°]")
    ax.set_ylabel("latitude geodésica [°]")
    ax.set_title(title, fontsize=11)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.85)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    logger.info(
        "Traço salvo: %s (%d pontos, %d segmentos, %d nós ascendentes, lat máx %.3f°)",
        out,
        t.size,
        len(segments),
        lon_nodes.size,
        np.degrees(np.max(np.abs(lat))),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--preset", choices=sorted(PRESETS), default="sso", help="órbita pronta")
    group.add_argument(
        "--elements",
        nargs=6,
        type=float,
        metavar=("A_KM", "E", "I_DEG", "RAAN_DEG", "ARGP_DEG", "PSI_DEG"),
        help="elementos da órbita (km e graus), em vez de --preset",
    )
    parser.add_argument("--propagator", choices=sorted(PROPAGATORS), default="mean-j2")
    parser.add_argument("--duration-h", type=float, default=24.0, help="duração [h] (padrão 24)")
    parser.add_argument("--step-s", type=float, default=10.0, help="passo [s] (padrão 10)")
    parser.add_argument("--out", type=Path, help="arquivo PNG de saída")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.elements is not None:
        name, elements, stem = "órbita personalizada", tuple(args.elements), "custom"
    else:
        name, elements = PRESETS[args.preset]
        stem = args.preset.replace("-", "_")
    a_km, e, i_deg = elements[:3]
    prop = make_propagator(args.propagator, *elements)
    title = (
        f"Traço no solo — {name} · {args.propagator} · {args.duration_h:g} h a {args.step_s:g} s\n"
        f"a = {a_km:.3f} km, e = {e:g}, i = {i_deg:g}° · época {EPOCH:%Y-%m-%d %H:%M} UTC"
    )
    out = args.out or REPORTS / f"ground_track_{stem}.png"
    plot_ground_track(prop, args.duration_h * 3600.0, args.step_s, title, out)


if __name__ == "__main__":
    main()
