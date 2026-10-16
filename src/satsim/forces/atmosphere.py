"""Atmosfera exponencial por partes (Sem3 §4, eq. 19; tabela de Vallado).

Em cada faixa de altitude [h0, próxima h0) a densidade é uma exponencial local (Sem3 eq. 19)

    ρ(h) = ρ0 · exp(−(h − h0)/H),

com ρ0 e a escala de altura H ajustados por faixa. Um único H global não serve: entre 200 e
1000 km H vai de 37 a 268 km (Sem3 Tabela 1 e §4.3), e extrapolar uma só exponencial erra a
densidade por ordens de grandeza (Sem3 Figura 1). A tabela por partes é a discretização da
integral exata da equação hidrostática (Sem3 eq. 18).

Limites do modelo
-----------------
- **Atividade solar moderada, estática.** ρ real a 400 km varia por um fator ≳ 10 ao longo do
  ciclo solar e sobe em horas numa tempestade geomagnética (Sem3 §7); nada disso está aqui.
- **Abaixo de 100 km não há dados**: devolve a densidade da primeira faixa (100 km), com aviso.
  O modelo é para LEO; abaixo de ~120 km a reentrada é questão de minutos.
- **Acima de 1000 km** usa a última faixa (H = 268 km), extrapolada.
"""

import logging
from collections.abc import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

logger = logging.getLogger(__name__)

FloatOrArray = float | NDArray[np.float64]

VALLADO: tuple[tuple[float, float, float], ...] = (
    # (h0 [km], ρ0 [kg/m³], H [km]) — Vallado, exponencial por partes (Apêndice A.3 do prompt
    # da Semana 3; reproduz a Tabela 3 do Sem3, arredondada no PDF)
    (100.0, 5.297e-7, 5.877),
    (110.0, 9.661e-8, 7.263),
    (120.0, 2.438e-8, 9.473),
    (130.0, 8.484e-9, 12.636),
    (140.0, 3.845e-9, 16.149),
    (150.0, 2.07e-9, 22.523),
    (180.0, 5.464e-10, 29.740),
    (200.0, 2.789e-10, 37.105),
    (250.0, 7.248e-11, 45.546),
    (300.0, 2.418e-11, 53.628),
    (350.0, 9.518e-12, 53.298),
    (400.0, 3.725e-12, 58.515),
    (450.0, 1.585e-12, 60.828),
    (500.0, 6.967e-13, 63.822),
    (600.0, 1.454e-13, 71.835),
    (700.0, 3.614e-14, 88.667),
    (800.0, 1.17e-14, 124.640),
    (900.0, 5.245e-15, 181.050),
    (1000.0, 3.019e-15, 268.000),
)


class ExponentialAtmosphere:
    """Densidade ρ(h) por uma tabela exponencial por partes (Sem3 eq. 19).

    Args:
        table: faixas ``(h0 [km], ρ0 [kg/m³], H [km])`` com h0 crescente (padrão: Vallado).

    Raises:
        ValueError: tabela vazia ou com h0 não estritamente crescente.
    """

    def __init__(self, table: Sequence[tuple[float, float, float]] = VALLADO) -> None:
        arr = np.asarray(table, dtype=float)
        if arr.ndim != 2 or arr.shape[0] == 0 or arr.shape[1] != 3:
            raise ValueError("a tabela deve ter linhas (h0 [km], ρ0 [kg/m³], H [km])")
        if np.any(np.diff(arr[:, 0]) <= 0.0):
            raise ValueError("h0 da tabela deve ser estritamente crescente")
        self._h0 = arr[:, 0] * 1e3  # m
        self._rho0 = arr[:, 1]
        self._scale_height = arr[:, 2] * 1e3  # m
        self._warned_below = False

    @property
    def h_min(self) -> float:
        """Base da primeira faixa [m] (abaixo dela não há dados)."""
        return float(self._h0[0])

    def density(self, h_m: ArrayLike) -> FloatOrArray:
        """Densidade ρ [kg/m³] na altitude ``h_m`` [m], vetorizada.

        Abaixo de ``h_min`` devolve ρ0 da primeira faixa (com um aviso por instância); acima da
        última faixa, extrapola a última exponencial.

        Returns:
            ρ com o shape de ``h_m`` (float para escalar).
        """
        h = np.asarray(h_m, dtype=float)
        idx = np.clip(np.searchsorted(self._h0, h, side="right") - 1, 0, self._h0.size - 1)
        below = h < self._h0[0]
        if np.any(below) and not self._warned_below:
            logger.warning(
                "altitude abaixo de %.0f km (mín. %.1f km): a tabela não tem dados; usando a "
                "densidade da primeira faixa",
                self._h0[0] / 1e3,
                float(np.min(h)) / 1e3,
            )
            self._warned_below = True
        dh = np.where(below, 0.0, h - self._h0[idx])
        rho = self._rho0[idx] * np.exp(-dh / self._scale_height[idx])
        return float(rho) if rho.ndim == 0 else rho
