# ADR 0001 — Arquitetura do satsim

## Status

Proposto — 2026-09-29. Atualizado na S5-06 (2026-11-06): `Satellite` e `Constellation`
implementados (S5-01); constelação na cena (S5-06). Atualizado na S6-02 (2026-11-10):
`ConstellationDesign` (S6-01), R7 concluído (6 ou 18 satélites pelo seletor do viewer), variantes
e famílias na cena, interpolação Hermite no viewer. Atualizado na S6-03 (2026-11-11): deploy
como site estático no GitHub Pages (seção Deploy).

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
- **TEME ≡ ECI do projeto (S3-02)**: o "ECI" do projeto é o equador verdadeiro e o equinócio
  médio *da data*, sem precessão/nutação, que é exatamente o quadro TEME em que o SGP4 devolve o
  estado. Não há conversão TEME → ECI: o `SGP4Propagator` entrega o estado já no quadro do
  projeto, e a passagem ao ECEF é Rz(GMST), como para os outros propagadores (diferença para o
  ITRS do Skyfield ≤ 86 m em 3 dias, por UT1 ≈ UTC e ausência de movimento polar).
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

Assinaturas resumidas. Implementados: `Propagator` (S1-05; `propagate_utc` na S3-02),
`GroundStation` (S3-03, `satsim.observer`), `Satellite` (S5-01, `satsim.satellite`) e
`Constellation` (S5-01, `satsim.constellation`, com a fábrica `sso_constellation` da receita da
Semana 4).

```python
class Propagator(ABC):
    @property
    @abstractmethod
    def epoch(self) -> datetime: ...  # UTC, tz-aware

    @abstractmethod
    def propagate(self, t: ArrayLike) -> tuple[NDArray, NDArray]:
        """t: segundos desde `epoch`, shape () ou (N,).
        Retorna (r, v) em ECI, SI, sempre com shape (N, 3)."""

    def propagate_utc(self, epoch: datetime, t: ArrayLike) -> tuple[NDArray, NDArray]:
        """S3-02 (não abstrato): estado em `epoch + t`, para qualquer `epoch`. Converte para
        segundos desde `self.epoch` e chama `propagate`. O Cowell exige t ≥ 0 (ValueError antes
        da sua época); Kepler, MeanJ2 e SGP4 aceitam t negativo."""


@dataclass(frozen=True)
class Satellite:  # S5-01
    id: str
    name: str
    propagator: Propagator
    color: str = "#00e5ff"  # cor CSS na cena
    cd: float | None = None  # coeficiente de arrasto (informativo; a física está no propagador)
    area_over_mass: float | None = None  # A/m [m²/kg]

    @property
    def epoch(self) -> datetime: ...  # a do propagador
    def propagate(self, t):
        return self.propagator.propagate(t)

    def propagate_utc(self, epoch, t):
        return self.propagator.propagate_utc(epoch, t)


class Constellation:  # S5-01
    def __init__(self, epoch: datetime, satellites: Sequence[Satellite]): ...  # ids únicos

    epoch: datetime  # época comum; converte t para cada satélite
    satellites: tuple[Satellite, ...]
    ids: list[str]  # também len(), iteração e índice

    def propagate(self, t) -> tuple[NDArray, NDArray]: ...  # cada um com shape (P, N, 3)
    def propagate_utc(self, epoch, t) -> tuple[NDArray, NDArray]: ...  # idem, (P, N, 3)


def sso_constellation(
    epoch, models, n_sats=6, a=..., ltan_h=22.5, i=None, e=None, color_palette=...
) -> Constellation: ...  # 91/6 num plano (Sem4 §8)


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
| R7 | Escala 1 → 6 → 18 satélites sem alterar o código do viewer: satélites entram por dados/configuração. **Concluído (S6-02):** `CONSTELLATION_DESIGNS` (6 e 18) é a única fonte de P; a demo gera as cenas de 6 e 18 e o viewer troca entre elas pelo seletor, desenhando qualquer número de membros e famílias a partir do JSON. |
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

Python calcula; o viewer só desenha. Gerado por `satsim.app.scene_export.build_scene` (S1-08;
subponto geodésico WGS84 desde a S2-07; estação, footprint, passagens e TLE desde a S3-06; LTAN,
β e eclipse desde a S4-06; constelação, nós, faixa e manutenção desde a S5-06; famílias,
grade diária, variantes e precisão compacta desde a S6-02). Aceita
`Satellite`, `SceneSatellite` ou uma `Constellation` inteira (S5-01).
Aceita qualquer `Propagator` (Kepler, MeanJ2, Cowell, SGP4), cada um com a sua época.

```jsonc
{
  "meta": {"schema_version": 1,
           "epoch_utc": "2026-09-29T12:00:00Z", "duration_s": 21600, "step_s": 10,
           "frame": "ECI simplificado (x -> equinócio vernal)",
           "geodesy": "WGS84",                    // S2-07: lla geodésico (campo novo)
           "el_min_deg": 10,                      // S3-06 (opcional): elevação mínima
           "tle": {"name": "ISS (ZARYA)", "norad_id": 25544,          // S3-06 (opcional)
                   "epoch_utc": "2026-10-22T23:36:04.211Z",
                   "age_days_at_start": 0.69, "source": "network"},   // network|cache|bundled
           "constellation": {"name": "Constelação 91/6",                 // S5-06 (opcional)
                             "n_revs": 91, "n_days": 6, "n_sats": 6, "grid_deg": 3.956044,
                             "member_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
                             "roll_deg": 23.1957,                       // null sem faixa
                             "tracks_per_day": 91, "daily_grid_deg": 3.956044,  // S6-02
                             "families": 1},                            //   (opcionais)
           "variants": [{"id": "p6", "label": "6 satélites", "href": "index.html"},  // S6-02
                        {"id": "p18", "label": "18 satélites",                       // (opcional)
                         "href": "constelacao18.html"}]},
  "earth": {"gmst0_rad": 1.234, "gmst_rate_rad_s": 7.2921158553e-5},
  "satellites": [
    {"id": "sso1", "name": "SSO-1", "color": "#00e5ff",
     "eci_km": [[x, y, z], ...],                  // N pontos, um por passo
     "lla":    [[lat_deg, lon_deg, alt_km], ...], // subponto geodésico WGS84 (S2-07)
     "footprint_deg": [12.5, ...],                // S3-06 (opcional): raio do footprint
     "ltan_h": [22.5, ...],                       // S4-06 (opcional): LTAN [h] em [0, 24)
     "beta_deg": [23.42, ...],                    // S4-06 (opcional): ângulo β [°]
     "eclipse": [false, ...],                     // S4-06 (opcional): na sombra da Terra
     "nodes": {"t_s": [1234.567, ...],            // S5-06 (opcional, só membros): nós
               "lon_deg": [-12.345678, ...]},     //   ascendentes na janela
     "swath_deg": 1.9955,                         // S5-06 (opcional): λ da faixa (eq. 45)
     "family": 0,                                 // S6-02 (opcional, só membros): j mod m
     "maintenance": {"a_dot_m_day": -23.95,       // S5-06 (opcional): resumo do
                     "days_to_track_correction": 6.92, ...}}  //   MaintenanceStatus (S5-05)
  ],
  "stations": [{"id": "sbs", "name": "São Bento do Sapucaí",
                "lat_deg": -22.69, "lon_deg": -45.73, "alt_m": 900,
                "looks": {"sso1": {"az_deg": [...], "el_deg": [...],   // S3-06 (opcional):
                                   "range_km": [...]}}}],            // por satélite e amostra
  "passes": [{"sat_id": "iss", "station_id": "sbs",                    // S3-06
              "aos_utc": "2026-10-24T00:10:18.123Z", "tca_utc": "...", "los_utc": "...",
              "max_el_deg": 18.2, "az_aos_deg": 236.1, "az_tca_deg": 180.3,
              "az_los_deg": 124.6, "duration_s": 287.4}]
}
```

- Amostras uniformes: `N = duration_s / step_s + 1`, instante `k` = `epoch_utc + k·step_s`;
  `lon_deg ∈ [−180°, 180°)`. O viewer interpola tratando o salto de ±180° e quebra o traço 2D no
  antimeridiano. `earth` permite girar a Terra (modo inercial) com θG ≈ `gmst0 + rate·t`.
- **Ajuste ao prompt:** `meta.schema_version`, para o viewer detectar cenas de formato antigo
  quando o contrato evoluir (passagens na Semana 3).
- **S2-07 (compatível, `schema_version` continua 1):** `lla` passa a ser geodésico WGS84
  (`frames.eci_to_ecef` + `geodesy.ecef_to_geodetic`): latitude da normal ao elipsoide e
  altitude elipsoidal. Novo campo `meta.geodesy = "WGS84"`; sem ele (cenas da Semana 1) o
  subponto é esférico (lat = arcsin(z/r), alt = r − R⊕). A diferença chega a ~0,2° em latitude e
  ~21 km em altitude (polos). O viewer só muda o texto que exibe.
- **S3-06 (compatível, `schema_version` continua 1; todos os campos novos são opcionais e só
  aparecem quando `build_scene` recebe `el_min_deg`):**
  - `satellites[i].footprint_deg[]`: raio angular geocêntrico do círculo de visibilidade a
    `el_min`, por amostra, λ = arccos(R⊕/(R⊕ + h)·cos ε) − ε (`observer.footprint_half_angle`,
    Sem4 §9.4), com h = a altitude publicada em `lla`. O viewer desenha o círculo com geometria
    de esfera (destino a distância angular λ), não com equações orbitais;
  - `stations[j].looks[sat_id] = {az_deg[], el_deg[], range_km[]}`: o satélite visto da estação
    (`GroundStation.look_angles` sobre o mesmo r_ECEF do subponto), por amostra;
  - `passes[]`: passagens acima de `el_min` (`passes.find_passes`, as mesmas de
    `scripts/list_passes.py`), instantes ISO 8601 UTC arredondados ao milissegundo;
  - `meta.el_min_deg` e `meta.tle` (época, idade no início da cena e fonte do TLE);
  - `meta.epoch_utc` ganha milissegundos quando a época não é um segundo inteiro (ex.: a época
    de um TLE); a época da cena é arredondada ao milissegundo antes da amostragem, para que
    amostras, passagens e o relógio do viewer partam do mesmo instante.
- **S4-06 (compatível, `schema_version` continua 1; campos opcionais que só aparecem quando
  `build_scene` recebe `sun_geometry=True`), por satélite e amostra:**
  - `ltan_h[]`: hora local do nó ascendente [h] em [0, 24), do Ω osculador e da ascensão reta do
    Sol (`ephemeris.ltan_hours`, Sem4 eq. 28). O viewer interpola pelo caminho curto (23,9 h →
    0,1 h passa por 24 h);
  - `beta_deg[]`: ângulo entre o plano da órbita e a direção do Sol (`ephemeris.beta_angle`,
    Sem4 eq. 29);
  - `eclipse[]`: booleano, na sombra cilíndrica da Terra (`ephemeris.in_shadow_cylindrical`).
    O viewer usa a amostra mais próxima (não interpola) e esmaece o marcador do satélite.
  - Sem os campos (cenas antigas) o painel Informações mostra "—" nas três linhas.
- **S5-06 (compatível, `schema_version` continua 1; campos opcionais que só aparecem quando
  `build_scene` recebe `constellation_info`, um `SceneConstellation`):**
  - `meta.constellation`: nome, repetição N/D, número de membros, passo da grade no equador
    (360°/N), ids dos membros e a rolagem da faixa (`roll_deg`, `null` sem faixa);
  - em cada membro, `nodes = {t_s[], lon_deg[]}`: os nós ascendentes na janela
    (`groundtrack.node_crossings`, S5-02), em segundos desde o início da cena;
  - em cada membro, `swath_deg`: o alcance transversal λ = arcsin((r/R⊕)·sin η) − η da eq. 45
    do Sem4 (`coverage.swath_half_angle`) para a rolagem η, no raio médio r da órbita na janela;
  - em cada membro com indicador, `maintenance`: os campos de `MAINTENANCE_FIELDS` do
    `MaintenanceStatus` (S5-05); valores infinitos (sem decaimento) viram `null`;
  - satélites fora de `member_ids` não ganham esses campos. Sem `meta.constellation` (cenas
    antigas) o viewer se comporta como antes: legenda simples, todos os satélites ligados,
    sem as camadas do mapa e com "—" nas linhas da constelação. Com ela, o viewer agrupa a
    legenda, desliga os não membros no início, desenha a grade do dia (nós com `t_s` ≤ t) e a
    faixa (polígono a ±λ do traço recente) e calcula a fase Δu a partir de `eci_km`.
- **S6-02 (compatível, `schema_version` continua 1; campos opcionais):**
  - com `SceneConstellation.families` (via `SceneConstellation.from_design`, a partir de um
    `ConstellationDesign`), `meta.constellation` ganha `tracks_per_day` (N·P/D; `null` se D não
    divide P), `daily_grid_deg` (360°·D/(N·P)) e `families` (m = P/D), e cada membro ganha
    `family` (j mod m). A rolagem da faixa da demo passa a ser a de cobertura diária do desenho
    (`ConstellationDesign.roll_for_daily_coverage`): 23,1957° com 6 (a mesma da S5-06) e
    8,1835° com 18;
  - `meta.variants = [{id, label, href}]`: as cenas irmãs, com `href` relativo (na mesma pasta;
    `build_scene` recusa caminho absoluto ou URL). Com elas, o viewer mostra o seletor
    "6 satélites | 18 satélites"; a troca navega para o irmão levando no hash da URL o instante
    (UTC absoluto), a velocidade, as camadas, os satélites fora da constelação ligados e o
    selecionado (`#t=…&v=…&grid=…&swath=…&play=…&on=…&sel=…&p=…`; um membro vai para o de fase
    mais próxima na outra cena: s4 de 18 ↔ S2 de 6). Funciona em `file://` e por HTTP;
  - `build_scene(..., precision=ScenePrecision)`: casas decimais das séries por amostra. O
    padrão, `FULL_PRECISION`, é o das semanas anteriores (a cena de 6 da demo não muda);
    `COMPACT_PRECISION` (1 m em posição, 0,01° em az/el; ≈ 18% menos bytes) é usado nas cenas de
    18 para cada arquivo ficar abaixo de 15 MB;
  - sem os campos, o viewer se comporta como antes: legenda da constelação sem subdivisões,
    contador "X de N" com a grade de 360°/N e sem seletor.
- **Interpolação no viewer (S6-02):** a posição entre amostras é uma Hermite cúbica
  (Catmull-Rom, tangentes por diferença central das amostras vizinhas) no vetor unitário de
  (lat, lon) e na altitude, projetada de volta para lat/lon: geometria de desenho, sem equação
  orbital. A interpolação linear em lat/lon das semanas anteriores erra muito mais que a flecha
  da corda (r·θ²/8), porque a longitude varia rápido perto da latitude máxima: 3,3 km a 20 s e
  29 km a 60 s no Cowell completo. Com Hermite no vetor unitário: 0,03 km a 20 s, 0,28 km a 60 s
  e 1,1 km a 120 s (o passo da cena de 6 dias com 18). Com passo > 30 s, os traços ganham pontos
  interpolados a cada 30 s. `tests/test_scene_export.py` reimplementa o método e mede o erro.
- **Desempenho (S6-02, Chrome headless, rasterização por software, `#bench=300`):** quadro 2D
  com tudo ligado de 3,0 ms (8 objetos) e 6,3 ms (20 objetos, 1 dia) a 9,1 ms (20 objetos,
  6 dias), máximo 14 ms; só o JavaScript do quadro 2D, 0,47 ms (8) e 3,9 ms (20); montagem dos
  dados do globo ≤ 0,13 ms. O globo WebGL não roda no headless deste ambiente e não foi medido.

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
única opção que atende R1, R2 e R9–R10 juntos, sem servidor nem token. O plano original era
hospedar o viewer num serviço Python (Streamlit ou FastAPI) na Semana 6; a S6-03 o publica como
**site estático** (seção Deploy), sem servidor de aplicação. **Atualização
S6-02:** as cenas de 6 e 18 satélites são **arquivos irmãos**, não um arquivo com as duas (juntas
passariam de 15 MB), ligados por `meta.variants` e links relativos; navegar entre arquivos
locais funciona offline, sem servidor.
**Restrição de `file://`:** navegadores bloqueiam módulos ES e recursos carregados de `file://`;
usamos scripts clássicos (UMD) e embutimos textura (data URI) e cena (`<script>`) no
`demo/index.html`. **Vista 3D:** Terra fixa (ECEF) a partir de `lla` nesta fase; modo inercial
(ECI) no backlog, já suportado pelo contrato (`eci_km`, `earth`).

### Deploy (S6-03)

**Decisão:** o simulador é publicado como **site estático no GitHub Pages**. Não há servidor de
aplicação: a física já sai de Python como cena, e o viewer só desenha.

- **Build:** `scripts/build_site.py [--offline] [-o site/]` monta `site/`, que não é versionado:
  - `index.html`: uma página inicial leve (≈ 6 kB, sem cena nem script), com o texto para
    leigos, quatro botões, a imagem da demo e a data e a fonte dos TLEs;
  - `img/demo_constelacao18.png`: a imagem da demo;
  - `cenas/`: as quatro cenas de `build_demo.build_all`. São as mesmas da demo `--all`, ligadas
    por `meta.variants`.

  O script imprime o tamanho de cada arquivo e o total, e falha se a página inicial tiver
  ≥ 1 MB ou se alguma cena tiver ≥ 15 MB.
- **Links relativos, nada externo:** os `href`/`src` são relativos e as cenas embutem CSS, JS,
  globe.gl e textura. O mesmo `site/` funciona em `https://<usuário>.github.io/<repositório>/`
  (um subcaminho, por isso nada começa com `/`), num `python -m http.server` local e em
  `file://`. `tests/test_build_site.py` verifica que todo link interno aponta para um arquivo
  existente dentro do site.
- **TLEs:** sem `--offline`, o build usa `get_tle`, na ordem CelesTrak → cache → embutidos no
  pacote. Se a CelesTrak falhar no CI, o site é montado com o fallback e a página inicial avisa,
  com a fonte e a época do TLE. Para a constelação 91/6 isso não importa, porque ela não depende
  de TLE.
- **Workflow** `.github/workflows/pages.yml`, que roda em push na `main` e por `workflow_dispatch`:
  - job `test`: Python 3.11, `pip install -e .[dev]`, ruff (check e format) e pytest;
  - job `deploy`, com `needs: test`: `build_site.py`, depois `actions/configure-pages`,
    `actions/upload-pages-artifact` e `actions/deploy-pages`.

  Usa só actions oficiais (`actions/*`), sem segredos. As permissões `pages: write` e
  `id-token: write` ficam apenas no job de deploy, e `concurrency: pages` evita dois deploys ao
  mesmo tempo. O `ci.yml` continua rodando lint e testes em 3.11 e 3.12 em todo push.
- **Tamanho:** cerca de 41 MB por deploy (cenas de 6,7 a 11,6 MB), dentro dos limites do Pages
  (100 MB por arquivo, 1 GB por site). O Pages serve com compressão (gzip), então a transferência
  é menor que o arquivo.
- **Fora do escopo:** servidor de aplicação, propagação no navegador, contas de usuário e
  domínio próprio. A ativação do Pages (Source: GitHub Actions) e o push ficam com o mantenedor.

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
- **TEME ≡ ECI** (S3-02): o SGP4 devolve estados em TEME, que é o próprio ECI do projeto
  (equinócio médio da data, sem precessão/nutação); ver Convenções. O SGP4 usa internamente as
  constantes WGS72 (o TLE é ajustado com elas); o resto do projeto usa WGS84 / Sem4.
- **Elementos médios de primeira ordem em J2**: bom para projeto de órbita, mas não serve para a
  ISS por mais de poucos dias sem TLE fresco.
- **Taxa de GMST vs. `OMEGA_EARTH`**: diferem levemente (≈ 5 m/dia no equador); a conversão
  ECI → ECEF usa a fórmula completa do GMST.
- **Subponto**: geodésico WGS84 desde a S2-07 (na Semana 1 era esférico). O globo 3D continua
  sendo uma esfera: desenha o `lla` geodésico sobre ela, o que é só uma escolha de exibição.
- **Cena pré-amostrada**: o viewer só mostra a janela exportada (ex.: 48 h a cada 20 s); tempos
  fora dela exigem gerar nova cena. Com 8 satélites o HTML cresce ≈ 1,2 MB por 1 000 amostras;
  a demo de 6 dias usa passo de 60 s (11,2 MB).
- **Faixa no mapa**: desenhada com geometria de esfera a ±λ da perpendicular ao traço; é a área
  alcançável com a rolagem, não o apontamento real da câmera (fora do escopo).

## Consequências

- Positivas: trocar de propagador ou ir de 1 a 18 satélites não muda a análise nem o viewer;
  números verificáveis contra PDFs e Skyfield; o mesmo viewer e as mesmas cenas servem à demo offline e ao site no GitHub Pages (S6-03).
- Negativas: contratos de shape e da cena exigem testes de vetorização e de export em cada card;
  o HTML da demo cresce com a janela e o número de satélites; as simplificações de quadro e tempo
  limitam a precisão a centenas de metros (aceitável para visibilidade).
