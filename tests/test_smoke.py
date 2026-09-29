"""Teste de fumaça: o pacote importa e expõe a versão."""

import satsim


def test_import_and_version() -> None:
    assert satsim.__version__ == "0.1.0"
