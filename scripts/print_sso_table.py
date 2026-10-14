"""Imprime a Tabela 1 do Sem4 (§4.3), órbitas heliossíncronas circulares, em Markdown.

Uso (na raiz do repositório, com o pacote instalado):

    python scripts/print_sso_table.py
    python scripts/print_sso_table.py --altitudes 300 510 800

Cada linha vem de ``satsim.orbit_design.sso_row``: inclinação heliossíncrona (Sem4 eq. 27),
período nodal (eq. 34) e revoluções por dia nodal (eq. 37), com teoria secular de J2.
"""

import argparse

from satsim.constants import R_EARTH
from satsim.orbit_design import sso_max_semi_major_axis, sso_row

TABLE1_ALTITUDES_KM = (300, 400, 500, 510, 600, 705, 786, 800, 1000, 1200)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--altitudes",
        nargs="+",
        type=float,
        default=TABLE1_ALTITUDES_KM,
        metavar="H_KM",
        help="altitudes [km] (padrão: as da Tabela 1)",
    )
    args = parser.parse_args()

    print("| h [km] | i [°] | Tnod [s] | rev/dia |")
    print("|---:|---:|---:|---:|")
    for h_km in args.altitudes:
        row = sso_row(h_km * 1e3)
        print(
            f"| {row.h_km:g} | {row.i_deg:.3f} | {row.nodal_period_s:.1f} | {row.rev_per_day:.3f} |"
        )
    a_max = sso_max_semi_major_axis()
    print(
        f"\nTeto (i = 180°): a_max = {a_max / 1e3:.1f} km, h_max = {(a_max - R_EARTH) / 1e3:.1f} km"
    )


if __name__ == "__main__":
    main()
