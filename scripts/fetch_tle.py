"""Obtém o TLE de um objeto (padrão: ISS) com cache e fallback offline, e mostra fonte e idade.

Uso (na raiz do repositório, com o pacote instalado):

    python scripts/fetch_tle.py                  # rede se o cache tiver > 2 h
    python scripts/fetch_tle.py --offline        # nunca acessa a rede
    python scripts/fetch_tle.py --norad 25544 --cache-dir /tmp/tle

Ordem (satsim.tle.get_tle): cache com menos de 2 h → CelesTrak (grava o cache) → cache de
qualquer idade → TLE embarcado no pacote. Os avisos (rede indisponível, TLE antigo) saem no
stderr.
"""

import argparse
import logging
from datetime import UTC, datetime

from satsim.tle import ISS_NORAD_ID, get_tle


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--norad", type=int, default=ISS_NORAD_ID, help="nº de catálogo NORAD")
    parser.add_argument("--offline", action="store_true", help="não acessa a rede")
    parser.add_argument(
        "--cache-dir", default=None, help="padrão: $SATSIM_CACHE_DIR ou data/cache/tle"
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.WARNING, format="aviso: %(message)s")

    now = datetime.now(UTC)
    result = get_tle(args.norad, offline=args.offline, cache_dir=args.cache_dir, now=now)
    tle = result.tle
    fetched = result.fetched_at.isoformat(timespec="seconds") if result.fetched_at else "—"
    print(tle.to_text(), end="")
    print(f"fonte:      {result.source}")
    print(f"baixado em: {fetched}")
    print(f"época:      {tle.epoch.isoformat(timespec='seconds')}")
    print(f"idade:      {result.age.total_seconds() / 86400:.2f} dias")


if __name__ == "__main__":
    main()
