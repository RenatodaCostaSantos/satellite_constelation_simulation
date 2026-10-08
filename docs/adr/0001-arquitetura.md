# ADR 0001 — Arquitetura do satsim

## Status

Proposto — 2026-09-29.

## Contexto

O satsim é um **simulador com interface gráfica**: qualquer pessoa abre e vê os satélites sobre
um globo 3D e, sincronizado na mesma tela, o traço no solo num mapa 2D. Evolui de 1 satélite (ISS
sobre São Bento do Sapucaí, semanas 1–4) para 6 (SSO 91/6, semana 5) e 18 com deploy web
(semanas 6–7). A física deve ser verificável contra os PDFs (Sem3, Sem4) e o Skyfield; o núcleo
numérico em Python existe para alimentar a GUI.

## Decisão

1. **Interface única de propagação**: `Propagator.propagate(t) → (r, v)`, vetorizada, estado ECI
   em SI com shape `(N, 3)`. Kepler, J2 médio, Cowell e SGP4 são implementações intercambiáveis.
2. **Satélite = dados + propagador.** Adicionar satélites é instanciar mais `Satellite`;
   `Constellation` empilha em `(P, N, 3)`. Nem a análise nem o viewer assumem um único satélite.
3. **Arrays NumPy vetorizados** no núcleo numérico; runtime só `numpy` e `scipy` (Skyfield só
   nos testes, como oráculo).
4. **A física fica em Python; o front-end só desenha**: Python exporta a "cena" (JSON amostrado
   no tempo); o viewer só interpola e projeta. Nenhuma equação orbital em JavaScript.
5. **Notação dos documentos** (ψ, ℓ, h) no código; docstrings citam a equação (ex.: "Sem4 eq. 23").
6. O módulo de tempo chama-se `astrotime.py` (não `time.py`), para não sombrear o `time` da stdlib.

## Estrutura de pacotes

```
satsim/                          (raiz do repositório)
├─ pyproject.toml · README.md · .gitignore (inclui demo/) · .github/workflows/ci.yml
├─ scripts/build_demo.py         (gera demo/index.html)
├─ demo/index.html               (GERADO, não versionado; abre com duplo clique)
├─ docs/
│  ├─ adr/0001-arquitetura.md
│  ├─ creditos.md                (texturas e bibliotecas JS: fonte, versão e licença)
│  └─ referencias/               (PDFs teóricos Sem3 e Sem4)
├─ src/satsim/
│  ├─ __init__.py · constants.py · astrotime.py · elements.py        (semana 1)
│  ├─ secular.py                 (taxas seculares de J2, período nodal) (semana 2, S2-01)
│  ├─ forces/ __init__.py · j2.py (acelerações para o Cowell)         (semana 2, S2-02)
│  ├─ propagators/
│  │  ├─ __init__.py · base.py (ABC Propagator) · kepler.py          (semana 1)
│  │  └─ mean_j2.py · cowell.py · sgp4_ref.py                        (semanas 2–3)
│  ├─ frames.py · geodesy.py · sun.py · observer.py · visibility.py  (semanas 2–4)
│  ├─ groundtrack.py             (subponto, segmentos, nós ascendentes) (semana 2, S2-04)
│  ├─ analysis.py                (elementos osculadores e médios, ajustes) (semana 2, S2-02)
│  ├─ orbit_design.py · satellite.py · constellation.py · coverage.py (semanas 2–5)
│  └─ app/                       (GUI)
│     ├─ __init__.py · scene_export.py                               (semana 1, S1-08)
│     ├─ static/ viewer.html · viewer.js · viewer.css
│     │         vendor/ (globe.gl UMD) · textures/ (Terra)          (semana 1, S1-08)
│     └─ serviço web Streamlit/FastAPI                               (semana 6)
└─ tests/
   ├─ test_smoke.py · test_constants.py · test_elements.py
   ├─ test_kepler_propagator.py · test_astrotime.py · test_scene_export.py
   ├─ test_secular.py · test_mean_j2.py · test_j2_force.py · test_cowell.py   (semana 2)
   ├─ test_frames_geodesy.py · test_groundtrack.py · test_orbit_design_sso.py (semana 2)
   └─ golden/reference_values.py · test_golden_consistency.py
```

## Convenções

- **Unidades**: SI (m, s, kg, m/s). Graus só em entrada/saída e em nomes com sufixo `_deg`
  (analogamente `_km`, `_s`). A cena JSON usa km e graus, com o sufixo no nome da chave.
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
          └─► scene_export: amostras no tempo (eci_km, lla) + estações + passagens ─► JSON
              ─► viewer (globo 3D + mapa 2D + controles de tempo)
```

Cada etapa é uma função pura sobre arrays; o eixo de satélites `P`, quando existe, atravessa
todas as etapas sem laços explícitos.

## Interface gráfica (GUI)

### Requisitos

| # | Requisito |
|---|---|
| R1 | Globo 3D com Terra texturizada (textura pública), fundo espacial com estrelas, brilho atmosférico e câmera controlável (arrastar para girar, roda para zoom). |
| R2 | Satélites em movimento contínuo (interpolação fluida, ≥ 30 fps), como marcadores luminosos com nome, cor própria, trilha recente e órbita prevista. |
| R3 | Mapa 2D (equiretangular, mesma textura) na mesma tela: subponto atual, traço passado (cheio) e futuro (tracejado), quebra correta no antimeridiano (±180°, sem linha atravessando o mapa) e grade de lat/lon. |
| R4 | Sincronia: uma única fonte de tempo (`simTime`) comanda 3D e 2D; arrastar a linha do tempo atualiza os dois painéis. |
| R5 | Controles: play/pause, velocidade (1x, 10x, 60x, 300x, 1000x), linha do tempo, relógio UTC visível, botão "agora/início", mostrar/ocultar cada satélite. |
| R6 | Estação de solo (São Bento do Sapucaí, ≈ −22,69°, −45,73°) marcada nos dois painéis; a partir da Semana 4, destaque quando o satélite estiver visível. |
| R7 | Escala 1 → 6 → 18 satélites sem alterar o código do viewer: satélites entram por dados/configuração. |
| R8 | Impacto visual: tema escuro, cores vivas, animação fluida, tipografia limpa; layout responsivo (3D e 2D lado a lado em tela larga, empilhados em tela estreita). |
| R9 | Uso simples: abre com um comando ou duplo clique, sem configuração; legenda e dicas na tela. |
| R10 | Offline após clonar: bibliotecas JS e texturas versionadas (`static/vendor`, `static/textures`), com créditos e licenças em `docs/creditos.md`. |

### Layout-alvo

```
+--------------------------------------------------------------------------+
| satsim   UTC 2026-09-29 12:00:00  [>][||] [1x 10x 60x 300x 1000x] [agora] |
+-------------------------------------+------------------------------------+
|                                     | MAPA 2D — traço no solo            |
|            GLOBO 3D                 | subponto + passado/futuro          |
|     satélites em movimento          | estação de solo marcada            |
|     órbitas, trilhas, estação       +------------------------------------+
|                                     | Legenda · satélites · informações  |
+-------------------------------------+------------------------------------+
| |------------------o-----------------------------|  linha do tempo       |
+--------------------------------------------------------------------------+
```

### Contrato de dados da cena (JSON)

Python calcula; o viewer só desenha. Gerado por `satsim.app.scene_export.build_scene` (S1-08).

```jsonc
{
  "meta": {"schema_version": 1,
           "epoch_utc": "2026-09-29T12:00:00Z", "duration_s": 21600, "step_s": 10,
           "frame": "ECI simplificado (x -> equinócio vernal)"},
  "earth": {"gmst0_rad": 1.234, "gmst_rate_rad_s": 7.2921158553e-5},
  "satellites": [
    {"id": "sso1", "name": "SSO-1", "color": "#00e5ff",
     "eci_km": [[x, y, z], ...],                  // N pontos, um por passo
     "lla":    [[lat_deg, lon_deg, alt_km], ...]} // subponto (aprox. esférica na Semana 1)
  ],
  "stations": [{"id": "sbs", "name": "São Bento do Sapucaí",
                "lat_deg": -22.69, "lon_deg": -45.73, "alt_m": 900}],
  "passes": []                                    // preenchido a partir da Semana 3
}
```

- Amostras uniformes: `N = duration_s / step_s + 1`, instante `k` = `epoch_utc + k·step_s`;
  `lon_deg ∈ [−180°, 180°)`. O viewer interpola tratando o salto de ±180° e quebra o traço 2D no
  antimeridiano. `earth` permite girar a Terra (modo inercial) com θG ≈ `gmst0 + rate·t`.
- **Ajuste ao prompt:** `meta.schema_version`, para o viewer detectar cenas de formato antigo
  quando o contrato evoluir (passagens na Semana 3).

### Decisão de stack

| Alternativa | Impacto visual | Offline (`file://`) | Esforço | Deploy web | Licença / tokens |
|---|---|---|---|---|---|
| **globe.gl (three.js) + canvas 2D** | Alto: globo texturizado, atmosfera, estrelas, arcos e pontos prontos | Sim, com bundle UMD e textura em data URI | Baixo–médio | Estático; qualquer host | MIT; sem tokens |
| CesiumJS | Muito alto (globo geoespacial completo) | Difícil: web workers e assets não carregam de `file://`; bundle grande | Médio–alto | Estático, pesado | Apache-2.0; imagens/terreno padrão exigem token Cesium ion |
| Plotly (JS/Python) | Médio: globo `scattergeo` sem textura realista | Sim (`plotly.min.js` embutido) | Baixo | Fácil | MIT; sem tokens |
| PyVista/VTK | Alto em 3D | Sim, mas só desktop | Médio | Ruim: exige trame/servidor | MIT/BSD; sem tokens |
| Streamlit puro | Baixo–médio: widgets, sem animação 3D a ≥ 30 fps | Não: exige servidor Python | Baixo | Fácil (Streamlit Cloud) | Apache-2.0; sem tokens |

**Decisão:** viewer HTML/JS estático e autocontido: **globe.gl (UMD, versão fixada)** no 3D e
**canvas 2D** no mapa equiretangular com a mesma textura, alimentados pelo JSON do Python. É a
única opção que atende R1, R2 e R9–R10 juntos, sem servidor nem token. Na **Semana 6**, um
serviço Python (Streamlit ou FastAPI) hospeda o mesmo viewer, sem mudá-lo.
**Restrição de `file://`:** navegadores bloqueiam módulos ES e recursos carregados de `file://`;
usamos scripts clássicos (UMD) e embutimos textura (data URI) e cena (`<script>`) no
`demo/index.html`. **Vista 3D:** Terra fixa (ECEF) a partir de `lla` nesta fase; modo inercial
(ECI) no backlog, já suportado pelo contrato (`eci_km`, `earth`).

## Alternativas rejeitadas

- **Lista de objetos Python por instante** (um objeto "estado" por t): legível, porém ordens de
  grandeza mais lento (laço Python por instante) e difícil de vetorizar em N instantes × P
  satélites. Arrays `(P, N, 3)` resolvem os dois.
- **Skyfield/astropy no runtime**: dependências pesadas que esconderiam justamente o que
  queremos implementar e entender. Ficam como oráculo de validação nos testes.
- **Notação da literatura (M, h)**: colidiria com `h` = altitude dos documentos e quebraria a
  rastreabilidade das equações. Usamos ψ, ℓ e h como nos PDFs, com a tabela acima como ponte.
- **Walker Delta (vários planos)**: a heliossíncrona exige o mesmo LTAN para todos os
  satélites, logo um único plano orbital; os satélites se distinguem apenas por ψ0.
- **Física em JavaScript** (propagar no navegador): duplicaria o código e quebraria a
  verificação contra PDFs e Skyfield, que só existe em Python.

## Simplificações e limitações conhecidas

- **ECI simplificado**: sem precessão nem nutação; o "ECI" é o equatorial da data, aproximado.
- **UT1 ≈ UTC**: erro de até 0,9 s, ≈ 0,4 km no equador.
- **Segundos intercalares** são ignorados pelo `datetime`.
- **TEME ≠ ECI**: o SGP4 devolve estados em TEME; a conversão será tratada na Semana 3 (S3-02).
- **Elementos médios de primeira ordem em J2**: bom para projeto de órbita, mas não serve para a
  ISS por mais de poucos dias sem TLE fresco.
- **Taxa de GMST vs. `OMEGA_EARTH`**: diferem levemente (≈ 5 m/dia no equador); a conversão
  ECI → ECEF usa a fórmula completa do GMST.
- **Subponto esférico na Semana 1**: `lla` usa lat = arcsin(z/r) e alt = r − R⊕; o geodésico
  WGS84 entra na S2-03.
- **Cena pré-amostrada**: o viewer só mostra a janela exportada (ex.: 6 h a cada 10 s); tempos
  fora dela exigem gerar nova cena.

## Consequências

- Positivas: trocar de propagador ou ir de 1 a 18 satélites não muda a análise nem o viewer;
  números verificáveis contra PDFs e Skyfield; o mesmo viewer serve à demo offline e ao deploy.
- Negativas: contratos de shape e da cena exigem testes de vetorização e de export em cada card;
  o HTML da demo cresce com a janela e o número de satélites; as simplificações de quadro e tempo
  limitam a precisão a centenas de metros (aceitável para visibilidade).
