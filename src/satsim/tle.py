"""TLE (Two-Line Element set): parse, validação, download da CelesTrak, cache em disco e fallback.

O TLE é o formato em que o NORAD publica os elementos médios de um objeto em órbita; é a entrada
do SGP4 (S3-02). Formato de colunas (1-indexadas) usado aqui:

    linha 1: 1 | 3–7 nº de catálogo | 8 classificação | 19–20 ano | 21–32 dia do ano fracionário
             | 54–61 B* (ex.: " 38792-4" = 0,38792e-4) | 69 checksum
    linha 2: 2 | 3–7 nº de catálogo | 9–16 inclinação [°] | 53–63 movimento médio [rev/dia]
             | 69 checksum

``get_tle`` busca o TLE em camadas, para o simulador funcionar sem internet:

    1. cache em disco mais novo que ``max_cache_age`` (padrão 2 h) → sem rede;
    2. rede (CelesTrak) → grava o cache;
    3. falha de rede ou TLE inválido → cache de qualquer idade, com aviso;
    4. sem cache → TLE embarcado no pacote (``satsim/data/iss_fallback.tle``), com aviso.

A CelesTrak pede no máximo 1 consulta a cada 2 h por conjunto de dados, daí o padrão de
``max_cache_age``. A precisão do SGP4 cai com a idade do TLE (reboosts, arrasto): acima de
3 dias é emitido um aviso de precisão degradada, acima de 14 dias um aviso forte.

Simplificação: o número de catálogo é lido como inteiro de 5 dígitos; o formato Alpha-5
(letra no 1º dígito, catálogos ≥ 100 000) é rejeitado com ``TleError``.

Papel no simulador: este módulo fornece o TLE de onde saem as condições iniciais da simulação.
Ele não calcula posições: o SGP4 (``satsim.propagators.sgp4_ref``) transforma os elementos médios
do TLE em coordenadas (r, v), e o estado em t = 0 é a condição inicial do Cowell.

A ISS é o satélite de teste porque o TLE dela é público e atualizado várias vezes por dia (a
CelesTrak o redistribui do NORAD), o que permite comparar a simulação com a órbita real.
"""

import json
import logging
import os
import tempfile
import urllib.request
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from http.client import HTTPException
from importlib import resources
from pathlib import Path
from typing import Literal

from satsim import __version__

logger = logging.getLogger(__name__)

TLE_LINE_LENGTH = 69
ISS_NORAD_ID = 25544
CELESTRAK_URL = "https://celestrak.org/NORAD/elements/gp.php?CATNR={id}&FORMAT=TLE"
DEFAULT_CACHE_DIR = Path("data") / "cache" / "tle"
CACHE_DIR_ENV = "SATSIM_CACHE_DIR"
BUNDLED_TLE = "data/iss_fallback.tle"
AGE_WARNING = timedelta(days=3)
AGE_STRONG_WARNING = timedelta(days=14)

TleSource = Literal["network", "cache", "bundled"]


class TleError(ValueError):
    """TLE malformado (comprimento, checksum, catálogo, campos) ou indisponível."""


def tle_checksum(line: str) -> int:
    """Checksum módulo 10 de uma linha de TLE (sem o último caractere).

    Soma dos dígitos; cada ``-`` conta 1; os demais caracteres contam 0.
    """
    total = 0
    for ch in line[:68]:
        if ch.isdigit():
            total += int(ch)
        elif ch == "-":
            total += 1
    return total % 10


def _parse_epoch(field: str) -> datetime:
    """Campo de época "YYDDD.DDDDDDDD" → ``datetime`` UTC (ano < 57 → 20xx, senão 19xx)."""
    yy = int(field[:2])
    year = 2000 + yy if yy < 57 else 1900 + yy
    day_of_year = float(field[2:])
    return datetime(year, 1, 1, tzinfo=UTC) + timedelta(days=day_of_year - 1.0)


def _parse_implied_decimal(field: str) -> float:
    """Campo com ponto decimal implícito e expoente, ex.: " 38792-4" → 0,38792e-4."""
    s = field.strip()
    sign = "-" if s.startswith("-") else ""
    s = s.lstrip("+-")
    mantissa, exp_sign, exponent = s[:-2], s[-2], s[-1]
    return float(f"{sign}0.{mantissa}e{exp_sign}{exponent}")


@dataclass(frozen=True)
class Tle:
    """Um TLE validado. As linhas são guardadas como texto; as propriedades as interpretam.

    Raises:
        TleError: se as linhas não forem um TLE válido (ver ``_validate``).
    """

    name: str
    line1: str
    line2: str

    def __post_init__(self) -> None:
        self._validate()

    def _validate(self) -> None:
        for number, line in (("1", self.line1), ("2", self.line2)):
            if len(line) != TLE_LINE_LENGTH:
                raise TleError(
                    f"linha {number} com {len(line)} caracteres (esperado {TLE_LINE_LENGTH})"
                )
            if line[0] != number or line[1] != " ":
                raise TleError(f"linha {number} deve começar com '{number} ': {line[:2]!r}")
            if not line[68].isdigit() or tle_checksum(line) != int(line[68]):
                raise TleError(
                    f"checksum da linha {number} inválido: {line[68]!r}, "
                    f"calculado {tle_checksum(line)}"
                )
        if self.line1[2:7] != self.line2[2:7]:
            raise TleError(
                f"número de catálogo difere entre as linhas: "
                f"{self.line1[2:7]!r} × {self.line2[2:7]!r}"
            )
        try:
            self.norad_id  # noqa: B018 (força o parse de todos os campos)
            self.epoch  # noqa: B018
            self.inclination_deg  # noqa: B018
            self.mean_motion_rev_day  # noqa: B018
            self.bstar  # noqa: B018
        except (ValueError, IndexError) as exc:
            raise TleError(f"campo numérico inválido no TLE: {exc}") from exc

    @classmethod
    def from_text(cls, text: str) -> "Tle":
        """TLE a partir de 2 linhas (sem nome) ou 3 linhas (a 1ª é o nome).

        Linhas em branco e espaços/``\\r`` no fim de cada linha são ignorados.
        """
        lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
        if len(lines) == 3:
            name, line1, line2 = lines
        elif len(lines) == 2:
            name, (line1, line2) = "", lines
        else:
            raise TleError(f"esperadas 2 ou 3 linhas de TLE, encontradas {len(lines)}")
        return cls(name.strip(), line1, line2)

    def to_text(self) -> str:
        """Texto em 3 linhas (ou 2, sem nome), terminado em nova linha."""
        lines = [self.name] if self.name else []
        return "\n".join([*lines, self.line1, self.line2]) + "\n"

    @property
    def norad_id(self) -> int:
        """Número de catálogo NORAD (colunas 3–7)."""
        return int(self.line1[2:7])

    @property
    def epoch(self) -> datetime:
        """Época do TLE, UTC tz-aware (linha 1, colunas 19–32)."""
        return _parse_epoch(self.line1[18:32])

    @property
    def inclination_deg(self) -> float:
        """Inclinação média [°] (linha 2, colunas 9–16)."""
        return float(self.line2[8:16])

    @property
    def mean_motion_rev_day(self) -> float:
        """Movimento médio [rev/dia] (linha 2, colunas 53–63)."""
        return float(self.line2[52:63])

    @property
    def bstar(self) -> float:
        """Termo de arrasto B* [1/raio terrestre] (linha 1, colunas 54–61)."""
        return _parse_implied_decimal(self.line1[53:61])

    def age(self, now: datetime | None = None) -> timedelta:
        """Idade do TLE: ``now − epoch`` (``now`` padrão = agora, UTC).

        Raises:
            ValueError: se ``now`` for naive.
        """
        if now is None:
            now = datetime.now(UTC)
        elif now.tzinfo is None or now.utcoffset() is None:
            raise ValueError("now deve ser tz-aware (ex.: tzinfo=timezone.utc)")
        return now - self.epoch


@dataclass(frozen=True)
class TleResult:
    """Resultado de ``get_tle``: o TLE, de onde veio, quando foi baixado e a idade.

    ``fetched_at`` é ``None`` para o TLE embarcado (a data de download está em
    ``docs/creditos.md``); ``age`` é a idade do TLE (agora − época), não a do download.
    """

    tle: Tle
    source: TleSource
    fetched_at: datetime | None
    age: timedelta


def fetch_celestrak(norad_id: int = ISS_NORAD_ID, timeout: float = 10) -> Tle:
    """Baixa o TLE de ``norad_id`` da CelesTrak (GP, formato TLE) e o valida.

    Raises:
        OSError: falha de rede (``urllib.error.URLError``, ``TimeoutError``...).
        TleError: resposta que não é um TLE válido do objeto pedido (ex.: página de erro).
    """
    url = CELESTRAK_URL.format(id=norad_id)
    request = urllib.request.Request(url, headers={"User-Agent": f"satsim/{__version__}"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        raw = response.read()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise TleError(f"resposta da CelesTrak não é texto UTF-8: {exc}") from exc
    tle = Tle.from_text(text)
    if tle.norad_id != norad_id:
        raise TleError(f"CelesTrak devolveu o catálogo {tle.norad_id}, pedido {norad_id}")
    return tle


def _resolve_cache_dir(cache_dir: str | Path | None) -> Path:
    if cache_dir is not None:
        return Path(cache_dir)
    return Path(os.environ.get(CACHE_DIR_ENV) or DEFAULT_CACHE_DIR)


def _read_cache(path: Path, norad_id: int) -> tuple[Tle, datetime] | None:
    """Lê o cache; arquivo ausente → ``None``; corrompido → aviso e ``None`` (nunca levanta)."""
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        tle = Tle.from_text(data["tle"])
        fetched_at = datetime.fromisoformat(data["fetched_at"])
        if fetched_at.tzinfo is None:
            raise ValueError("fetched_at sem fuso horário")
        if tle.norad_id != norad_id:
            raise ValueError(f"catálogo {tle.norad_id} no cache de {norad_id}")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        logger.warning("cache de TLE corrompido em %s, ignorado: %s", path, exc)
        return None
    return tle, fetched_at


def _write_cache(path: Path, tle: Tle, fetched_at: datetime) -> None:
    """Grava o cache de forma atômica: arquivo temporário no mesmo diretório + ``os.replace``.

    Uma falha no meio deixa o cache anterior intacto. Falhas de escrita só geram aviso.
    """
    payload = {
        "norad_id": tle.norad_id,
        "tle": tle.to_text(),
        "fetched_at": fetched_at.isoformat(),
        "source_url": CELESTRAK_URL.format(id=tle.norad_id),
    }
    tmp_name = None
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as tmp:
            tmp_name = tmp.name
            json.dump(payload, tmp, indent=2)
        os.replace(tmp_name, path)
        tmp_name = None
    except OSError as exc:
        logger.warning("não foi possível gravar o cache de TLE em %s: %s", path, exc)
    finally:
        if tmp_name is not None:
            Path(tmp_name).unlink(missing_ok=True)


def load_bundled_tle() -> Tle:
    """TLE da ISS embarcado no pacote (``satsim/data/iss_fallback.tle``)."""
    text = resources.files("satsim").joinpath(BUNDLED_TLE).read_text(encoding="utf-8")
    return Tle.from_text(text)


def _warn_if_old(tle: Tle, age: timedelta, source: str) -> None:
    days = age / timedelta(days=1)
    if age > AGE_STRONG_WARNING:
        logger.warning(
            "TLE MUITO DESATUALIZADO (fonte %s): %.1f dias desde a época %s; "
            "posições podem errar dezenas a centenas de km",
            source,
            days,
            tle.epoch.isoformat(),
        )
    elif age > AGE_WARNING:
        logger.warning(
            "TLE com %.1f dias (fonte %s, época %s): precisão degradada",
            days,
            source,
            tle.epoch.isoformat(),
        )


def get_tle(
    norad_id: int = ISS_NORAD_ID,
    max_cache_age: timedelta = timedelta(hours=2),
    offline: bool = False,
    cache_dir: str | Path | None = None,
    now: datetime | None = None,
) -> TleResult:
    """TLE de ``norad_id`` com cache e fallback em camadas (ver docstring do módulo).

    Args:
        norad_id: número de catálogo NORAD (padrão: ISS, 25544).
        max_cache_age: cache mais novo que isto é usado sem consultar a rede.
        offline: ``True`` nunca acessa a rede (cache de qualquer idade, depois o embarcado).
        cache_dir: diretório do cache; padrão ``$SATSIM_CACHE_DIR`` ou ``data/cache/tle``
            (relativo ao diretório de trabalho). Arquivo: ``<cache_dir>/<norad_id>.json``.
        now: instante de referência (UTC tz-aware) para idade do cache e do TLE; padrão agora.

    Raises:
        TleError: sem rede, sem cache e sem TLE embarcado para ``norad_id`` (só há o da ISS).
    """
    if now is None:
        now = datetime.now(UTC)
    elif now.tzinfo is None or now.utcoffset() is None:
        raise ValueError("now deve ser tz-aware (ex.: tzinfo=timezone.utc)")
    cache_path = _resolve_cache_dir(cache_dir) / f"{norad_id}.json"
    cached = _read_cache(cache_path, norad_id)

    result: TleResult | None = None
    if cached is not None and now - cached[1] <= max_cache_age:
        result = TleResult(cached[0], "cache", cached[1], cached[0].age(now))
    elif not offline:
        try:
            tle = fetch_celestrak(norad_id)
        except (OSError, HTTPException, TleError) as exc:
            logger.warning("falha ao baixar o TLE %d da CelesTrak: %s", norad_id, exc)
        else:
            _write_cache(cache_path, tle, now)
            result = TleResult(tle, "network", now, tle.age(now))

    if result is None and cached is not None:
        tle, fetched_at = cached
        logger.warning(
            "usando TLE do cache baixado em %s (%s)",
            fetched_at.isoformat(),
            "modo offline" if offline else "rede indisponível",
        )
        result = TleResult(tle, "cache", fetched_at, tle.age(now))

    if result is None:
        if norad_id != ISS_NORAD_ID:
            raise TleError(f"sem rede, sem cache e sem TLE embarcado para o catálogo {norad_id}")
        tle = load_bundled_tle()
        logger.warning("sem cache de TLE: usando o TLE embarcado no pacote (época %s)", tle.epoch)
        result = TleResult(tle, "bundled", None, tle.age(now))

    _warn_if_old(result.tle, result.age, result.source)
    return result
