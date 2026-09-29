# ADR 0001 — Arquitetura do satsim

## Status

Proposto — 2026-09-29.

## Contexto

O satsim simula a posição de satélites e seus horários de visibilidade sobre uma estação
terrestre. Começa com um satélite (a ISS sobre São Bento do Sapucaí, semanas 1–4), passa a uma
constelação heliossíncrona de 6 satélites (91/6, um plano, semana 5) e, se houver tempo, a uma
aplicação web com 18 satélites (semanas 6–7). O código precisa ser numericamente verificável
contra os documentos teóricos (Sem3: decaimento em LEO; Sem4: órbita heliossíncrona) e contra
o Skyfield, e crescer de 1 para N satélites sem reescrita.

## Decisão

1. **Toda propagação passa por uma única interface**: `Propagator.propagate(t) → (r, v)`,
   vetorizada em `t`, devolvendo estado ECI em SI com shape `(N, 3)`. Kepler, J2 médio, Cowell e
   SGP4 são implementações intercambiáveis dessa ABC.
2. **Satélite = dados + propagador.** Adicionar satélites é apenas instanciar mais `Satellite`;
   `Constellation` empilha os resultados em `(P, N, 3)`. Nenhum código de análise (visibilidade,
   cobertura) pode assumir um único satélite.
3. **Arrays NumPy vetorizados** em todo o núcleo numérico; nada de laços Python por instante.
4. **Runtime mínimo**: apenas `numpy` e `scipy`. Skyfield entra só como dependência de
   desenvolvimento, para validação nos testes.
5. **Notação dos documentos** (ψ, ℓ, h) nos nomes do código, para que cada fórmula seja
   rastreável até a equação de origem, citada na docstring (ex.: "Sem4 eq. 23").
6. O módulo de tempo chama-se `astrotime.py` (o plano original dizia `time.py`), para não
   sombrear o módulo `time` da biblioteca padrão.

## Estrutura de pacotes

```
satsim/                          (raiz do repositório)
├─ pyproject.toml · README.md · .gitignore · .github/workflows/ci.yml
├─ docs/
│  ├─ adr/0001-arquitetura.md
│  └─ referencias/               (PDFs teóricos Sem3 e Sem4)
├─ src/satsim/
│  ├─ __init__.py · constants.py · astrotime.py · elements.py        (semana 1)
│  ├─ propagators/
│  │  ├─ __init__.py · base.py (ABC Propagator) · kepler.py          (semana 1)
│  │  └─ mean_j2.py · cowell.py · sgp4_ref.py                        (semanas 2–3)
│  ├─ frames.py · geodesy.py · sun.py · observer.py · visibility.py  (semanas 2–4)
│  ├─ orbit_design.py · satellite.py · constellation.py · coverage.py (semanas 2–5)
│  └─ app/                                                           (semana 6)
└─ tests/
   ├─ test_smoke.py · test_constants.py · test_elements.py
   ├─ test_kepler_propagator.py · test_astrotime.py
   └─ golden/reference_values.py · test_golden_consistency.py
```

## Convenções

- **Unidades**: SI — m, s, kg, m/s. Graus só em entrada/saída e em nomes com sufixo `_deg`
  (analogamente `_km`, `_s`).
- **Ângulos**: radianos internamente.
- **Tempo**: `datetime` UTC *tz-aware* (naive é rejeitado com `ValueError`). Na propagação,
  `t` = segundos (float) desde a `epoch` do propagador. Datas Julianas via `astrotime`.
- **Quadros**: ECI = equatorial inercial (x → equinócio vernal, z → eixo polar);
  ECEF = Rz(θG)·ECI, com θG = GMST. ENU local no observador.
- **Shapes**: tempos `(N,)`; posições e velocidades `(N, 3)`; constelação `(P, N, 3)`.
  Entradas escalares são promovidas a `(1,)`.

| Símbolo | Código | Significado | Literatura usual |
|---|---|---|---|
| ψ | `psi` | anomalia média | M |
| ℓ | `ell` | momento angular específico r × v | h |
| h | `h` | altitude | — |
| ν | `nu` | anomalia verdadeira | ν, f |
| ω | `argp` | argumento do perigeu | ω |
| Ω | `raan` | ascensão reta do nó ascendente | Ω |
| u | `u` | argumento de latitude, ω + ν | u |
| ω⊕ | `OMEGA_EARTH` | rotação da Terra | ω⊕ |

## Interfaces

Assinaturas ilustrativas. Nesta semana só `Propagator` é implementado (S1-05); as demais serão
implementadas nas semanas 3–5.

```python
class Propagator(ABC):
    @property
    @abstractmethod
    def epoch(self) -> datetime: ...  # UTC, tz-aware

    @abstractmethod
    def propagate(self, t: ArrayLike) -> tuple[NDArray, NDArray]:
        """t: segundos desde `epoch`, shape () ou (N,).
        Retorna (r, v) em ECI, SI, sempre com shape (N, 3)."""


@dataclass(frozen=True)
class Satellite:
    id: str
    name: str
    propagator: Propagator
    cd: float | None = None  # coeficiente de arrasto
    area_over_mass: float | None = None  # A/m [m²/kg]

    def propagate(self, t):
        return self.propagator.propagate(t)


class Constellation:
    epoch: datetime  # época comum; converte t para cada satélite
    satellites: Sequence[Satellite]

    def propagate(self, t) -> tuple[NDArray, NDArray]: ...  # cada um com shape (P, N, 3)


@dataclass(frozen=True)
class GroundStation:
    name: str
    lat: float
    lon: float
    alt: float  # rad, rad, m (WGS84)

    @classmethod
    def from_degrees(cls, name, lat_deg, lon_deg, alt_m): ...
    def ecef(self) -> NDArray: ...  # (3,)
    def look_angles(self, r_ecef) -> tuple[NDArray, NDArray, NDArray]: ...  # az, el, alcance
```

`Constellation.propagate` converte o `t` da época comum para cada satélite:
`t_i = t + (constellation.epoch − sat.propagator.epoch).total_seconds()`.

## Fluxo de dados

```
Propagator ─(r, v) ECI─► rotação por GMST θG ─► ECEF ─► geodésico / observador (ENU)
          ─► elevação (az, el, alcance) ─► visibilidade (janelas de passagem)
```

Cada etapa é uma função pura sobre arrays; o eixo de satélites `P`, quando existe, atravessa
todas as etapas sem laços explícitos.

## Alternativas rejeitadas

- **Lista de objetos Python por instante** (um objeto "estado" por t): legível, porém ordens de
  grandeza mais lento (laço Python por instante) e difícil de vetorizar em N instantes × P satélites. Arrays `(P, N, 3)`
  resolvem os dois.
- **Skyfield/astropy no runtime**: dependências pesadas que esconderiam justamente o que
  queremos implementar e entender. Ficam como oráculo de validação nos testes.
- **Notação da literatura (M, h)**: colidiria com `h` = altitude dos documentos e quebraria a
  rastreabilidade das equações. Usamos ψ, ℓ e h como nos PDFs, com a tabela acima como ponte.
- **Walker Delta (vários planos)**: a heliossíncrona exige o mesmo LTAN para todos os
  satélites, logo um único plano orbital; os satélites se distinguem apenas por ψ0.

## Simplificações e limitações conhecidas

- **ECI simplificado**: sem precessão nem nutação; o "ECI" é o equatorial da data, aproximado.
- **UT1 ≈ UTC**: erro de até 0,9 s, ≈ 0,4 km no equador.
- **Segundos intercalares** são ignorados pelo `datetime`.
- **TEME ≠ ECI**: o SGP4 devolve estados em TEME; a conversão será tratada na Semana 3 (S3-02).
- **Elementos médios de primeira ordem em J2**: bom para projeto de órbita, mas não serve para a
  ISS por mais de poucos dias sem TLE fresco.
- **Taxa de GMST vs. `OMEGA_EARTH`**: diferem levemente (≈ 5 m/dia no equador); a conversão
  ECI → ECEF usa a fórmula completa do GMST.

## Consequências

- Positivas: trocar de modelo de propagação ou passar de 1 para 18 satélites não muda o código
  de análise; cada número é verificável contra os PDFs e o Skyfield; runtime leve, fácil de
  empacotar para a web.
- Negativas: toda função nova precisa respeitar os contratos de shape `(N, 3)` / `(P, N, 3)`,
  o que exige testes de vetorização em cada card; as simplificações de quadro e tempo limitam a
  precisão a centenas de metros, aceitável para visibilidade, mas não para determinação de
  órbita precisa.
