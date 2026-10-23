"""Lista as passagens de um satélite (padrão: ISS) sobre São Bento do Sapucaí.

Uso (na raiz do repositório, com o pacote instalado):

    python scripts/list_passes.py                          # ISS, próximos 3 dias, el ≥ 10°
    python scripts/list_passes.py --days 7 --el-min 5 --offline
    python scripts/list_passes.py --csv passagens.csv
    python scripts/list_passes.py --tle-file tests/data/iss_sample.tle --start tle-epoch

O TLE vem de satsim.tle.get_tle (cache com menos de 2 h → CelesTrak → cache antigo → TLE
embarcado), ou de --tle-file. A janela começa agora (UTC, arredondado ao minuto), num instante
ISO 8601 (--start 2026-10-23T00:00:00Z) ou na época do TLE (--start tle-epoch, reprodutível).

O script imprime a fonte e a idade do TLE no início da janela, com aviso acima de 3 dias, e a
tabela de passagens em UTC e no horário de Brasília (BRT = UTC−03:00). "Passagem" aqui é o
satélite acima da elevação mínima (visibilidade por rádio), não visível a olho nu.
"""

import argparse
import csv
import logging
from datetime import UTC, datetime, timedelta
from pathlib import Path

from satsim.observer import SAO_BENTO_DO_SAPUCAI
from satsim.passes import Pass, find_passes, format_passes
from satsim.propagators import SGP4Propagator
from satsim.tle import AGE_WARNING, ISS_NORAD_ID, Tle, get_tle

logger = logging.getLogger("list_passes")


def parse_start(text: str | None, tle: Tle) -> datetime:
    """Início da janela: ``None`` → agora (arredondado ao minuto); ``"tle-epoch"``; ou ISO 8601."""
    if text is None:
        return datetime.now(UTC).replace(second=0, microsecond=0)
    if text == "tle-epoch":
        return tle.epoch
    start = datetime.fromisoformat(text.replace("Z", "+00:00"))
    if start.tzinfo is None:
        raise SystemExit("--start precisa de fuso horário (ex.: 2026-10-23T00:00:00Z)")
    return start.astimezone(UTC)


def write_csv(passes: list[Pass], path: Path) -> None:
    fields = ["sat_id", "aos_utc", "tca_utc", "los_utc", "max_elevation_deg", "az_aos_deg",
              "az_tca_deg", "az_los_deg", "duration_s", "starts_inside", "ends_inside"]  # fmt: skip
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        for p in passes:
            writer.writerow([
                p.sat_id, p.aos.isoformat(), p.tca.isoformat(), p.los.isoformat(),
                f"{p.max_elevation_deg:.3f}", f"{p.az_aos_deg:.2f}", f"{p.az_tca_deg:.2f}",
                f"{p.az_los_deg:.2f}", f"{p.duration_s:.1f}", p.starts_inside, p.ends_inside,
            ])  # fmt: skip


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--norad", type=int, default=ISS_NORAD_ID, help="nº de catálogo NORAD")
    parser.add_argument("--days", type=float, default=3.0, help="duração da janela [dias]")
    parser.add_argument("--el-min", type=float, default=10.0, help="elevação mínima [°]")
    parser.add_argument("--step", type=float, default=20.0, help="passo da varredura [s]")
    parser.add_argument("--start", default=None, help="agora (padrão), ISO 8601 ou tle-epoch")
    parser.add_argument("--tle-file", type=Path, default=None, help="TLE de um arquivo")
    parser.add_argument("--offline", action="store_true", help="não acessa a rede")
    parser.add_argument(
        "--cache-dir", default=None, help="padrão: $SATSIM_CACHE_DIR ou data/cache/tle"
    )
    parser.add_argument("--csv", type=Path, default=None, help="grava as passagens em CSV")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.WARNING, format="aviso: %(message)s")

    if args.tle_file is not None:
        tle = Tle.from_text(args.tle_file.read_text(encoding="utf-8"))
        source = f"arquivo {args.tle_file}"
    else:
        result = get_tle(args.norad, offline=args.offline, cache_dir=args.cache_dir)
        tle, source = result.tle, result.source

    start = parse_start(args.start, tle)
    age = start - tle.epoch
    station = SAO_BENTO_DO_SAPUCAI
    print(f"TLE:      {tle.name or '—'} (catálogo {tle.norad_id}), fonte: {source}")
    print(f"época:    {tle.epoch:%Y-%m-%d %H:%M:%S} UTC")
    print(f"idade:    {age / timedelta(days=1):.2f} dias no início da janela")
    if age > AGE_WARNING:
        print("AVISO:    TLE com mais de 3 dias: a previsão perde precisão (reboosts e arrasto)")
    print(f"estação:  {station.name} ({station.lat_deg:.2f}°, {station.lon_deg:.2f}°, "
          f"{station.alt:.0f} m)")  # fmt: skip
    print(f"janela:   {start:%Y-%m-%d %H:%M} UTC + {args.days:g} dia(s), el ≥ {args.el_min:g}°")
    print()

    sat_id = tle.name or str(tle.norad_id)
    passes = find_passes(SGP4Propagator(tle), station, start, timedelta(days=args.days),
                         el_min_deg=args.el_min, step_s=args.step, sat_id=sat_id)  # fmt: skip
    print(format_passes(passes))
    if args.csv is not None:
        write_csv(passes, args.csv)
        print(f"\n{len(passes)} passagens gravadas em {args.csv}")


if __name__ == "__main__":
    main()
