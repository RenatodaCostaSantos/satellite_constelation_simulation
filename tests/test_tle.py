"""Testes da ingestão de TLE (S3-01): parse, validação, cache, fallback e idade.

Referências: tests/golden/reference_values.py (ISS_TLE_A1, Apêndice A.1 do prompt da Semana 3;
TLE_AGING, Apêndice A.6). Nenhum teste acessa a rede: ``urllib.request.urlopen`` é substituído
por uma função que falha, e cada teste que simula a CelesTrak instala a sua.
"""

import io
import json
import logging
import os
import urllib.error
import urllib.request
from datetime import datetime, timedelta
from pathlib import Path

import pytest
from golden.reference_values import ISS_TLE_A1, TLE_AGING

from satsim import tle as tle_mod
from satsim.tle import Tle, TleError, fetch_celestrak, get_tle, load_bundled_tle, tle_checksum

SAMPLE_PATH = Path(__file__).parent / "data" / "iss_sample.tle"
A1_TEXT = f"{ISS_TLE_A1.name}\n{ISS_TLE_A1.line1}\n{ISS_TLE_A1.line2}\n"
A1_EPOCH = datetime.fromisoformat(ISS_TLE_A1.epoch_utc_iso)
NOW = A1_EPOCH + timedelta(days=1)  # TLE com 1 dia: nenhum aviso de idade


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    """Garante que nenhum teste acessa a rede, e isola o diretório de cache do ambiente."""
    calls = []

    def forbidden(*args, **kwargs):
        calls.append(args)
        raise AssertionError("teste tentou acessar a rede")

    monkeypatch.setattr(urllib.request, "urlopen", forbidden)
    monkeypatch.delenv("SATSIM_CACHE_DIR", raising=False)
    return calls


class FakeCelestrak:
    """Substituto de ``urlopen``: devolve ``body`` ou levanta ``error``; conta as chamadas."""

    def __init__(self, body: bytes = b"", error: Exception | None = None):
        self.body, self.error, self.requests = body, error, []

    def __call__(self, request, timeout=None):
        self.requests.append(request)
        if self.error is not None:
            raise self.error
        return io.BytesIO(self.body)


def with_checksum(line: str) -> str:
    """Recalcula o checksum (último caractere) de uma linha de 69 caracteres."""
    return line[:68] + str(tle_checksum(line))


def write_cache(cache_dir: Path, fetched_at: datetime, text: str = A1_TEXT) -> Path:
    path = cache_dir / "25544.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"tle": text, "fetched_at": fetched_at.isoformat()}))
    return path


# --- Parse do TLE de exemplo (A.1) ------------------------------------------------------------


def test_sample_tle_fields_a1():
    tle = Tle.from_text(SAMPLE_PATH.read_text())
    assert tle.name == ISS_TLE_A1.name
    assert tle.norad_id == ISS_TLE_A1.norad_id
    assert abs((tle.epoch - A1_EPOCH).total_seconds()) <= 1.0
    assert tle.epoch.tzinfo is not None
    assert tle.inclination_deg == pytest.approx(ISS_TLE_A1.inclination_deg, abs=1e-12)
    assert tle.mean_motion_rev_day == pytest.approx(ISS_TLE_A1.mean_motion_rev_day, abs=1e-12)
    assert tle.bstar == pytest.approx(ISS_TLE_A1.bstar, rel=1e-12)


def test_three_and_two_lines_give_same_tle():
    three = Tle.from_text(A1_TEXT)
    two = Tle.from_text(f"{ISS_TLE_A1.line1}\n{ISS_TLE_A1.line2}\n")
    assert two.name == ""
    assert (two.line1, two.line2) == (three.line1, three.line2)
    assert two.epoch == three.epoch and two.bstar == three.bstar


def test_from_text_tolerates_crlf_and_padded_name():
    # formato real da CelesTrak: nome com espaços até 24 caracteres e fim de linha CRLF
    text = f"{ISS_TLE_A1.name:<24}\r\n{ISS_TLE_A1.line1}\r\n{ISS_TLE_A1.line2}\r\n"
    assert Tle.from_text(text) == Tle.from_text(A1_TEXT)


def test_to_text_round_trip():
    tle = Tle.from_text(A1_TEXT)
    assert Tle.from_text(tle.to_text()) == tle


# --- Rejeições --------------------------------------------------------------------------------


def _bump_digit(line: str, col: int) -> str:
    return line[:col] + str((int(line[col]) + 1) % 10) + line[col + 1 :]


@pytest.mark.parametrize(
    ("line1", "line2"),
    [
        pytest.param(_bump_digit(ISS_TLE_A1.line1, 40), ISS_TLE_A1.line2, id="checksum-l1"),
        pytest.param(ISS_TLE_A1.line1, _bump_digit(ISS_TLE_A1.line2, 19), id="checksum-l2"),
        pytest.param(ISS_TLE_A1.line1[:68], ISS_TLE_A1.line2, id="68-chars"),
        pytest.param(ISS_TLE_A1.line1 + "0", ISS_TLE_A1.line2, id="70-chars"),
        pytest.param(ISS_TLE_A1.line2, ISS_TLE_A1.line1, id="linhas-trocadas"),
        pytest.param(
            ISS_TLE_A1.line1,
            with_checksum(ISS_TLE_A1.line2[:2] + "25545" + ISS_TLE_A1.line2[7:]),
            id="catalogo-diferente",
        ),
        pytest.param(
            with_checksum(ISS_TLE_A1.line1[:20] + "x" + ISS_TLE_A1.line1[21:]),
            ISS_TLE_A1.line2,
            id="campo-nao-numerico",
        ),
    ],
)
def test_invalid_tle_raises(line1, line2):
    with pytest.raises(TleError):
        Tle("ISS", line1, line2)


@pytest.mark.parametrize("text", ["", ISS_TLE_A1.line1, A1_TEXT + A1_TEXT])
def test_from_text_wrong_line_count_raises(text):
    with pytest.raises(TleError):
        Tle.from_text(text)


def test_checksum_counts_minus_as_one():
    assert tle_checksum("-" * 68 + "0") == 68 % 10
    assert tle_checksum("1 2+A." + " " * 63) == 3


# --- Regra do ano de 2 dígitos ----------------------------------------------------------------


@pytest.mark.parametrize(("yy", "year"), [("56", 2056), ("57", 1957), ("00", 2000), ("99", 1999)])
def test_epoch_two_digit_year_rule(yy, year):
    line1 = with_checksum(ISS_TLE_A1.line1[:18] + yy + ISS_TLE_A1.line1[20:])
    tle = Tle("", line1, ISS_TLE_A1.line2)
    assert tle.epoch.year == year
    assert (tle.epoch.month, tle.epoch.day) == (12, 9 if year % 4 else 8)  # dia 343


# --- Idade ------------------------------------------------------------------------------------


def test_age_uses_injected_now():
    tle = Tle.from_text(A1_TEXT)
    now = tle.epoch + timedelta(days=2, hours=3)
    assert tle.age(now) == timedelta(days=2, hours=3)
    with pytest.raises(ValueError, match="tz-aware"):
        tle.age(now.replace(tzinfo=None))


@pytest.mark.parametrize(
    ("age_days", "expected"),
    [(1.0, None), (TLE_AGING.warning_days + 0.5, "precisão degradada"),
     (TLE_AGING.strong_warning_days + 1.0, "MUITO DESATUALIZADO")],
)  # fmt: skip
def test_age_warnings(tmp_path, caplog, age_days, expected):
    write_cache(tmp_path, A1_EPOCH)  # cache válido: só a idade do TLE gera aviso
    now = Tle.from_text(A1_TEXT).epoch + timedelta(days=age_days)
    with caplog.at_level(logging.WARNING, logger="satsim.tle"):
        result = get_tle(cache_dir=tmp_path, offline=True, now=now)
    assert result.age == timedelta(days=age_days)
    age_msgs = [r.message for r in caplog.records if "dias" in r.message]
    if expected is None:
        assert age_msgs == []
    else:
        assert len(age_msgs) == 1 and expected in age_msgs[0]


# --- get_tle com urlopen simulado -------------------------------------------------------------


def test_get_tle_network_success_writes_cache(tmp_path, monkeypatch):
    fake = FakeCelestrak(A1_TEXT.encode())
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    result = get_tle(cache_dir=tmp_path, now=NOW)
    assert result.source == "network"
    assert result.fetched_at == NOW
    assert result.age == NOW - result.tle.epoch
    assert result.tle == Tle.from_text(A1_TEXT)
    (request,) = fake.requests
    assert request.full_url == "https://celestrak.org/NORAD/elements/gp.php?CATNR=25544&FORMAT=TLE"
    assert request.get_header("User-agent").startswith("satsim/")
    cached = json.loads((tmp_path / "25544.json").read_text())
    assert Tle.from_text(cached["tle"]) == result.tle
    assert datetime.fromisoformat(cached["fetched_at"]) == NOW


def test_get_tle_second_call_within_max_age_uses_cache(tmp_path, monkeypatch):
    fake = FakeCelestrak(A1_TEXT.encode())
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    get_tle(cache_dir=tmp_path, now=NOW)
    result = get_tle(cache_dir=tmp_path, now=NOW + timedelta(hours=1, minutes=59))
    assert result.source == "cache"
    assert result.fetched_at == NOW
    assert len(fake.requests) == 1  # a 2ª chamada não consultou a rede


def test_get_tle_expired_cache_refetches(tmp_path, monkeypatch):
    write_cache(tmp_path, NOW - timedelta(hours=3))
    fake = FakeCelestrak(A1_TEXT.encode())
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    result = get_tle(cache_dir=tmp_path, now=NOW)
    assert result.source == "network" and len(fake.requests) == 1


@pytest.mark.parametrize(
    "error",
    [urllib.error.URLError("sem rede"), TimeoutError("timeout"), ConnectionResetError()],
)
def test_get_tle_network_failure_uses_old_cache(tmp_path, monkeypatch, caplog, error):
    fetched_at = NOW - timedelta(days=30)  # cache de qualquer idade serve como fallback
    write_cache(tmp_path, fetched_at)
    monkeypatch.setattr(urllib.request, "urlopen", FakeCelestrak(error=error))
    with caplog.at_level(logging.WARNING, logger="satsim.tle"):
        result = get_tle(cache_dir=tmp_path, now=NOW)
    assert result.source == "cache"
    assert result.fetched_at == fetched_at
    assert any("falha ao baixar" in r.message for r in caplog.records)
    assert any("usando TLE do cache" in r.message for r in caplog.records)


def test_get_tle_network_failure_without_cache_uses_bundled(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(urllib.request, "urlopen", FakeCelestrak(error=OSError("sem rede")))
    with caplog.at_level(logging.WARNING, logger="satsim.tle"):
        result = get_tle(cache_dir=tmp_path, now=NOW)
    assert result.source == "bundled"
    assert result.fetched_at is None
    assert result.tle == load_bundled_tle()
    assert any("embarcado" in r.message for r in caplog.records)
    assert not (tmp_path / "25544.json").exists()


def test_get_tle_offline_never_calls_network(tmp_path, no_network):
    assert get_tle(cache_dir=tmp_path, offline=True, now=NOW).source == "bundled"
    write_cache(tmp_path, NOW - timedelta(days=10))
    assert get_tle(cache_dir=tmp_path, offline=True, now=NOW).source == "cache"
    assert no_network == []


@pytest.mark.parametrize(
    "body",
    [
        b"<html><body><h1>500 Internal Server Error</h1></body></html>",
        b"No GP data found",
        b"\xff\xfe\x00",
        # TLE válido, mas de outro catálogo
        (
            with_checksum(ISS_TLE_A1.line1.replace("25544", "25545"))
            + "\n"
            + with_checksum(ISS_TLE_A1.line2.replace("25544", "25545"))
        ).encode(),
    ],
    ids=["html", "no-data", "binario", "outro-catalogo"],
)
def test_get_tle_invalid_response_is_a_failure(tmp_path, monkeypatch, body):
    monkeypatch.setattr(urllib.request, "urlopen", FakeCelestrak(body))
    result = get_tle(cache_dir=tmp_path, now=NOW)
    assert result.source == "bundled"
    assert not (tmp_path / "25544.json").exists()


def test_fetch_celestrak_rejects_wrong_catalog(monkeypatch):
    other = Tle(
        "", with_checksum(ISS_TLE_A1.line1.replace("25544", "25545")),
        with_checksum(ISS_TLE_A1.line2.replace("25544", "25545")),
    )  # fmt: skip
    monkeypatch.setattr(urllib.request, "urlopen", FakeCelestrak(other.to_text().encode()))
    with pytest.raises(TleError, match="catálogo"):
        fetch_celestrak(25544)


def test_get_tle_unknown_object_without_any_source_raises(tmp_path):
    with pytest.raises(TleError, match="sem TLE embarcado"):
        get_tle(norad_id=99999, cache_dir=tmp_path, offline=True, now=NOW)


def test_get_tle_cache_dir_from_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("SATSIM_CACHE_DIR", str(tmp_path / "env"))
    write_cache(tmp_path / "env", NOW)
    assert get_tle(offline=True, now=NOW).source == "cache"


def test_get_tle_rejects_naive_now(tmp_path):
    with pytest.raises(ValueError, match="tz-aware"):
        get_tle(cache_dir=tmp_path, offline=True, now=NOW.replace(tzinfo=None))


# --- Cache corrompido e escrita atômica -------------------------------------------------------


@pytest.mark.parametrize(
    "content",
    [
        "{ isto não é json",
        json.dumps({"tle": A1_TEXT}),  # sem fetched_at
        json.dumps({"tle": "lixo", "fetched_at": NOW.isoformat()}),
        json.dumps({"tle": A1_TEXT, "fetched_at": "2019-12-10T00:00:00"}),  # naive
        json.dumps([1, 2, 3]),
    ],
    ids=["json-invalido", "sem-fetched_at", "tle-invalido", "fetched_at-naive", "lista"],
)
def test_corrupted_cache_is_ignored(tmp_path, caplog, content):
    (tmp_path / "25544.json").write_text(content)
    with caplog.at_level(logging.WARNING, logger="satsim.tle"):
        result = get_tle(cache_dir=tmp_path, offline=True, now=NOW)
    assert result.source == "bundled"
    assert any("corrompido" in r.message for r in caplog.records)


def test_atomic_write_keeps_previous_cache_on_failure(tmp_path, monkeypatch, caplog):
    old_fetched = NOW - timedelta(hours=5)
    path = write_cache(tmp_path, old_fetched)
    before = path.read_bytes()

    def fail_midway(obj, fp, **kwargs):
        fp.write('{"tle": "parcial')  # escreve metade e falha
        raise OSError("disco cheio")

    monkeypatch.setattr(tle_mod.json, "dump", fail_midway)
    monkeypatch.setattr(urllib.request, "urlopen", FakeCelestrak(A1_TEXT.encode()))
    with caplog.at_level(logging.WARNING, logger="satsim.tle"):
        result = get_tle(cache_dir=tmp_path, now=NOW)
    assert result.source == "network"  # falha ao gravar o cache não derruba o get_tle
    assert path.read_bytes() == before  # cache anterior íntegro
    assert os.listdir(tmp_path) == ["25544.json"]  # temporário removido
    assert any("não foi possível gravar" in r.message for r in caplog.records)


# --- TLE embarcado ----------------------------------------------------------------------------


def test_bundled_tle_is_valid_iss():
    tle = load_bundled_tle()
    assert tle.norad_id == ISS_TLE_A1.norad_id
    assert 51.0 < tle.inclination_deg < 52.0
    assert 15.0 < tle.mean_motion_rev_day < 16.0
    assert tle.epoch.year >= 2019
