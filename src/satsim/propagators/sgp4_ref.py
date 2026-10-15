"""Propagador SGP4 (oráculo): a biblioteca ``sgp4`` sob a interface ``Propagator``.

O SGP4 é o modelo analítico para o qual o TLE foi gerado: ele lê os elementos *médios* do TLE e
devolve o estado em qualquer instante, com J2–J4, arrasto (via B*) e, para períodos longos,
Sol e Lua. É a referência ("oráculo") da posição da ISS real, contra a qual o Cowell é medido
(S3-03).

Três pontos a lembrar:

- **Constantes WGS72.** O SGP4 usa internamente as constantes WGS72 (R⊕ = 6378,135 km,
  µ = 398 600,8 km³/s², J2 = 1,082616e-3), porque é assim que o NORAD ajusta os TLEs. Não as
  troque por WGS84: os elementos médios seriam interpretados de forma ligeiramente errada.
  Fora daqui o projeto usa WGS84 / Sem4; o ponto de encontro é a saída em SI.
- **TEME = ECI do projeto, sem conversão.** O SGP4 devolve o estado no quadro TEME (True
  Equator, Mean Equinox: equador verdadeiro da data, eixo x no equinócio médio da data). É
  exatamente o "ECI" deste projeto (ADR 0001, Quadros), e a passagem para ECEF é Rz(GMST), como
  para os outros propagadores.
- **Precisão.** Para um TLE da ISS, ≈ 1 km na época, crescendo alguns km por dia (reboosts e
  arrasto mal modelado). Contra o Skyfield (mesmo SGP4, com UT1 e movimento polar), a posição
  ECEF difere ≈ 86 m no máximo em 3 dias (Apêndice A.1 do prompt da Semana 3).

A biblioteca trabalha em km e km/s; a conversão para SI é feita aqui, na fronteira.
"""

import logging
import math
from datetime import datetime

import numpy as np
from numpy.typing import ArrayLike, NDArray
from sgp4.api import SGP4_ERRORS, WGS72, Satrec

from satsim.constants import SECONDS_PER_DAY
from satsim.propagators.base import Propagator, as_time_array
from satsim.tle import Tle

logger = logging.getLogger(__name__)


class SGP4Propagator(Propagator):
    """Propaga um TLE com SGP4 (biblioteca ``sgp4``, constantes WGS72); saída em ECI (= TEME), SI.

    ``t`` pode ser negativo (antes da época do TLE): o SGP4 o permite. Instantes em que o SGP4
    falha (ex.: código 6, satélite decaído) saem como linhas ``NaN``, com um único
    ``logging.warning`` por chamada e os códigos em ``last_error_codes``; nenhuma exceção é
    levantada no meio de uma série.

    Args:
        tle: TLE validado (``satsim.tle.Tle``).
    """

    def __init__(self, tle: Tle) -> None:
        self._tle = tle
        self._satrec = Satrec.twoline2rv(tle.line1, tle.line2, WGS72)
        self.last_error_codes: NDArray[np.int_] = np.zeros(0, dtype=int)

    @classmethod
    def from_tle_text(cls, text: str) -> "SGP4Propagator":
        """Atalho: ``SGP4Propagator(Tle.from_text(text))`` (2 ou 3 linhas)."""
        return cls(Tle.from_text(text))

    @property
    def tle(self) -> Tle:
        """O TLE propagado."""
        return self._tle

    @property
    def epoch(self) -> datetime:
        """Época do TLE (t = 0), UTC tz-aware.

        Lida do TLE; difere do JD interno da biblioteca em ~10 µs (arredondamento do float).
        """
        return self._tle.epoch

    @property
    def mean_motion(self) -> float:
        """Movimento médio do TLE n [rad/s] (``no_kozai`` da biblioteca, em rad/min / 60)."""
        return self._satrec.no_kozai / 60.0

    @property
    def period(self) -> float:
        """Período 2π/n [s] do movimento médio do TLE."""
        return 2.0 * math.pi / self.mean_motion

    def propagate(self, t: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estado ECI (= TEME) em t segundos desde a época do TLE, vetorizado.

        Chama ``Satrec.sgp4_array(jd, fr)`` com ``jd = jdsatepoch`` e
        ``fr = jdsatepochF + t/86400``: a parte inteira do JD fica constante e só a fração
        carrega t, o que preserva a resolução do tempo.

        Args:
            t: escalar ou sequência (N,) de tempos [s], em qualquer ordem; pode ser negativo.

        Returns:
            (r, v) em ECI [m, m/s], cada um com shape (N, 3); linhas ``NaN`` onde o SGP4 falhou.
        """
        t = as_time_array(t)
        sat = self._satrec
        jd = np.full(t.shape, sat.jdsatepoch)
        fr = sat.jdsatepochF + t / SECONDS_PER_DAY
        errors, r_km, v_km_s = sat.sgp4_array(jd, fr)
        errors = np.asarray(errors, dtype=int)
        r = np.asarray(r_km, dtype=float).reshape(t.size, 3) * 1e3
        v = np.asarray(v_km_s, dtype=float).reshape(t.size, 3) * 1e3
        self.last_error_codes = errors
        failed = errors != 0
        if np.any(failed):
            r[failed] = np.nan
            v[failed] = np.nan
            codes = sorted(set(errors[failed].tolist()))
            logger.warning(
                "SGP4 falhou em %d de %d instantes (catálogo %d): %s",
                int(failed.sum()),
                t.size,
                self._tle.norad_id,
                "; ".join(f"código {c}: {SGP4_ERRORS.get(c, 'desconhecido')}" for c in codes),
            )
        return r, v
