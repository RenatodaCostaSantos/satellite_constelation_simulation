"""Constantes físicas e astronômicas do satsim (unidades SI).

Os valores são idênticos aos dos PDFs teóricos (Sem4: "Como se comporta um satélite em órbita
heliossíncrona?", tabela de símbolos, §2.2 e código §9.1).

Armadilha da normalização de J2 (Sem4 §2.1)
-------------------------------------------
Modelos de geopotencial como o EGM2008 tabulam coeficientes *totalmente normalizados*:
C̄20 = −4,8417e-4. O fator de normalização do grau 2, ordem 0, é N20 = √5 = 2,2361, logo

    C20 = N20 · C̄20 = −1,0826e-3   e   J2 = −C20 = 1,0826e-3.

Usar C̄20 diretamente no lugar de C20 erra J2 por um fator 2,24 e a inclinação heliossíncrona
sai completamente fora. Aqui J2 já é o coeficiente *não normalizado*.
"""

import math
from typing import Final

# --- Terra (Sem4, tabela de símbolos) ---------------------------------------------------------

MU: Final[float] = 3.986004418e14
"""Parâmetro gravitacional da Terra µ [m³/s²] (Sem4, tabela de símbolos)."""

R_EARTH: Final[float] = 6378.137e3
"""Raio equatorial da Terra R⊕ [m] (Sem4, tabela de símbolos)."""

OMEGA_EARTH: Final[float] = 7.292115e-5
"""Velocidade de rotação da Terra ω⊕ [rad/s] (Sem4, tabela de símbolos)."""

# --- Elipsoide WGS84 (NIMA TR8350.2) ----------------------------------------------------------

WGS84_A: Final[float] = R_EARTH
"""Semieixo maior do elipsoide WGS84 a [m] (= R⊕, raio equatorial)."""

WGS84_F: Final[float] = 1.0 / 298.257223563
"""Achatamento do WGS84 f [adimensional] (valor completo; o texto do Sem4 usa 1/298,257)."""

WGS84_B: Final[float] = WGS84_A * (1.0 - WGS84_F)
"""Semieixo menor (polar) do WGS84 b = a(1 − f) [m] ≈ 6 356 752,314245 m."""

WGS84_E2: Final[float] = WGS84_F * (2.0 - WGS84_F)
"""Excentricidade ao quadrado do WGS84 e² = f(2 − f) [adimensional] ≈ 6,69437999e-3."""

# --- Harmônicos zonais, não normalizados (Sem4 §2.2 e código §9.1) ----------------------------

J2: Final[float] = 1.08262668e-3
"""Harmônico zonal J2 [adimensional] (Sem4 código §9.1; o texto do §2.2 usa 1,082627e-3)."""

J3: Final[float] = -2.533e-6
"""Harmônico zonal J3 [adimensional] (Sem4 §2.2)."""

J4: Final[float] = -1.620e-6
"""Harmônico zonal J4 [adimensional] (Sem4 §2.2)."""

# --- Tempo e Sol --------------------------------------------------------------------------------

T_TROP_DAYS: Final[float] = 365.2422
"""Ano trópico [dias] (Sem4, tabela de símbolos; o código do PDF usa 365,24219, diferença
< 3e-8 relativa)."""

SECONDS_PER_DAY: Final[float] = 86400.0
"""Segundos por dia [s/dia]."""

OMEGA_SUN: Final[float] = 2.0 * math.pi / (T_TROP_DAYS * SECONDS_PER_DAY)
"""Taxa heliossíncrona alvo: movimento médio aparente do Sol [rad/s] (Sem4 eq. 26).
Vale ≈ 1,991064e-7 rad/s = 0,985647 °/dia."""

# --- Conversão de ângulos -----------------------------------------------------------------------

DEG2RAD: Final[float] = math.pi / 180.0
"""Fator graus → radianos."""

RAD2DEG: Final[float] = 180.0 / math.pi
"""Fator radianos → graus."""
