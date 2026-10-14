# satellite_constelation_simulation

**satsim** — simulador de posicionamento de satélites em Python, modular e fácil de estender
para mais satélites.

- **Fase 1 (semanas 1–4):** um satélite (a ISS), com posição e horários de visibilidade sobre
  São Bento do Sapucaí (SP).
- **Fase 2 (semana 5):** generalização para 6 satélites (constelação heliossíncrona 91/6, um
  plano).
- **Fase 3 (semanas 6–7, se houver tempo):** aplicação web com deploy e 18 satélites.

Todo propagador implementa `propagate(t) -> (r, v)` vetorizado; adicionar satélites é apenas
instanciar mais objetos. Unidades SI, ângulos em radianos internamente, tempo em UTC
(`datetime` com fuso).

## Instalação

Requer Python ≥ 3.11.

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

## Comandos

```bash
pytest               # testes
ruff check .         # lint
ruff format .        # formatação (use --check para só verificar)
```

## Estrutura

```
pyproject.toml               empacotamento, pytest e ruff
.github/workflows/ci.yml     CI: lint + testes em Python 3.11 e 3.12
docs/adr/                    decisões de arquitetura (ADR)
docs/creditos.md             créditos e licenças dos assets do viewer
src/satsim/                  código da biblioteca
src/satsim/app/              GUI: exportação da cena (Python) e viewer (static/)
scripts/build_demo.py        gera a demo visual demo/index.html
scripts/                     também: relatório da Semana 2, traço no solo, Tabela 1 (SSO),
                             TLE da ISS (fetch_tle.py)
reports/s2_comparacao/       figuras e relatório da Semana 2 (versionados)
tests/                       testes (pytest)
```

## Demo visual

![Demo do satsim: globo 3D e traço no solo 2D, fim da janela de 48 h, com as órbitas tipo ISS com e sem J2 deslocadas](docs/img/demo.png)

Globo 3D e mapa 2D com o traço no solo, na mesma tela e sincronizados. As órbitas são calculadas
em Python e o viewer só desenha. A demo cobre **48 h** (passo de 20 s) com três satélites:

| Satélite | Propagador | O que mostra |
|---|---|---|
| SSO-1 (J2) | `MeanJ2Propagator` | a órbita heliossíncrona de projeto 91/6 |
| ISS-like (exemplo J2) | `MeanJ2Propagator` | órbita tipo ISS (a = R⊕ + 420 km, i = 51,64°) com J2 |
| ISS-like sem J2 (Kepler) | `KeplerPropagator` | a mesma órbita sem J2, em cinza |

**O efeito de J2 fica visível:** as duas ISS-like começam juntas e se separam no mapa, porque J2
faz o plano da órbita girar (regressão do nó, ≈ −9,9° em 2 dias). Arraste a linha do tempo até o
fim da janela para ver os dois traços deslocados. A ISS-like é só um exemplo: a ISS real (TLE/SGP4)
entra na Semana 3.

**Gerar** (na raiz do repositório, com o venv ativo):

```bash
python scripts/build_demo.py      # cria demo/index.html (~4 MB, autocontido)
```

**Abrir:** dê dois cliques em `demo/index.html`, ou rode `xdg-open demo/index.html` no Linux
(`open` no macOS). Não precisa de servidor nem de internet: biblioteca, textura e cena vão
embutidas no arquivo. A pasta `demo/` não é versionada; gere de novo quando o código mudar.

**Adicionar um satélite:** acrescente uma linha à lista `SATELLITES` em
[`scripts/build_demo.py`](scripts/build_demo.py), por exemplo

```text
sat("kep2", "Quarto (i = 30°)", "#b26bff", MeanJ2Propagator, a_km=7200.0, e=0.01, i_deg=30.0, raan_deg=200.0, argp_deg=0.0, psi_deg=90.0),
```

(o quarto argumento é o propagador: `MeanJ2Propagator` com J2, `KeplerPropagator` sem J2)

e gere a demo de novo. O viewer lê quantos satélites houver na cena (testado com 18).

**Controles:**

| Controle | Ação |
|---|---|
| ▶ / ❚❚ (ou barra de espaço) | reproduzir / pausar |
| 1x · 10x · 60x · 300x · 1000x | velocidade da simulação |
| linha do tempo (embaixo) | arrastar para qualquer instante da janela de 48 h |
| ⏮ início | voltar ao início da janela |
| caixas na lista de satélites | mostrar / ocultar cada satélite |
| nome na lista de satélites | selecionar: o painel **Informações** mostra latitude, longitude e altitude geodésicas atuais |
| arrastar / roda do mouse no globo | girar / zoom |

No mapa, a linha cheia é o traço já percorrido e a tracejada é o futuro. São Bento do Sapucaí
aparece em rosa nos dois painéis. A vista mantém a Terra fixa (ECEF); o subponto é geodésico
WGS84 desde a Semana 2 (latitude da vertical local, altitude acima do elipsoide). Os créditos e
licenças dos assets estão em [`docs/creditos.md`](docs/creditos.md).

## Relatório da Semana 2

Comparação Kepler × MeanJ2 × Cowell (J2): o que J2 muda na órbita e o que cada propagador
entrega. Um único comando gera as figuras e o relatório em
[`reports/s2_comparacao/`](reports/s2_comparacao/):

```bash
python scripts/compare_propagators.py      # ~5 s; 5 figuras PNG + relatorio.md
```

O [`relatorio.md`](reports/s2_comparacao/relatorio.md) traz, para a SSO de projeto, uma órbita
tipo ISS e uma órbita de teste com e = 0,02 (10 dias a 60 s): a regressão do nó medida no Cowell
contra a teoria secular, a rotação do perigeu, a oscilação de curto período do semieixo, o erro
de posição do MeanJ2 e do Kepler contra o Cowell, e um parágrafo de conclusões. A lição central é
a diferença entre **elementos médios e osculadores**: os mesmos números nominais usados como
condição inicial do Cowell descrevem uma órbita média alguns km diferente, o que separa os dois
propagadores em ~1 000 km em 24 h; com os elementos médios medidos, a diferença cai para < 15 km.

Também na pasta: o traço no solo da SSO e da órbita tipo ISS
(`python scripts/plot_ground_track.py --preset sso` ou `--preset iss-like`), e a Tabela 1 do
Sem4 (órbitas heliossíncronas) sai de `python scripts/print_sso_table.py`.

## TLE da ISS (Semana 3)

A ISS real entra pelo seu TLE (Two-Line Element set), baixado da
[CelesTrak](https://celestrak.org/). Para o simulador funcionar sem internet,
`satsim.tle.get_tle()` busca o TLE em camadas:

1. **cache** em disco com menos de 2 h (a CelesTrak pede no máximo 1 consulta a cada 2 h);
2. **rede** (CelesTrak), gravando o cache;
3. falha de rede ou resposta inválida → **cache de qualquer idade**, com aviso;
4. sem cache → **TLE embarcado** no pacote (`src/satsim/data/iss_fallback.tle`), com aviso.

O cache fica em `data/cache/tle/<norad>.json` (relativo ao diretório de trabalho, não
versionado); mude com `--cache-dir` ou com a variável de ambiente `SATSIM_CACHE_DIR`. A precisão
cai com a idade do TLE: acima de 3 dias sai um aviso de precisão degradada, acima de 14 dias um
aviso forte.

```bash
python scripts/fetch_tle.py              # usa o cache (< 2 h) ou baixa da CelesTrak
python scripts/fetch_tle.py --offline    # nunca acessa a rede (cache ou TLE embarcado)
python scripts/fetch_tle.py --norad 25544 --cache-dir /tmp/tle
```

Saída (o TLE, depois fonte, data do download, época e idade; avisos no stderr):

```text
ISS (ZARYA)
1 25544U 98067A   26286.82178372  .00006592  00000+0  12831-3 0  9991
2 25544  51.6325  73.0089 0006710 259.3258 100.6975 15.48855783590140
fonte:      network
baixado em: 2026-10-14T11:47:21+00:00
época:      2026-10-13T19:43:22+00:00
idade:      0.67 dias
```

Os testes não acessam a rede: usam o TLE congelado de `tests/data/iss_sample.tle` (época
2019-12-09) e simulam a CelesTrak.

## Golden tests

[`tests/golden/reference_values.py`](tests/golden/reference_values.py) reúne, num único arquivo
versionado, os números de referência dos PDFs teóricos (`docs/referencias/`: Sem3, decaimento
em LEO; Sem4, órbita heliossíncrona): constantes, órbita de projeto 91/6, Tabelas 1, 3 e 4 do
Sem4, taxas seculares de J2 e valores de arrasto do Sem3. Os testes de cada semana importam
esses valores em vez de repetir números soltos, por exemplo
`from golden.reference_values import DESIGN_ORBIT`.

[`tests/test_golden_consistency.py`](tests/test_golden_consistency.py) confere a coerência
interna dos números (por exemplo, `86400 / Tnod ≈ rev/dia`, `360° / N ≈ grade`, e o mesmo valor
igual em tabelas diferentes), para pegar erros de transcrição. Uma discrepância encontrada no
próprio PDF fica registrada como `xfail(strict=True)`, com a explicação no `reason`.

Para adicionar um valor:

1. coloque-o no grupo certo (ou crie um grupo novo, como `@dataclass(frozen=True)`), exatamente
   como aparece no documento, com o sufixo da unidade no nome (`h_km`, `i_deg`, `tnod_s`);
2. cite a fonte num comentário (documento, seção, tabela ou equação);
3. se o número puder ser deduzido de outros, acrescente um teste de consistência.

## Roadmap

| Semana | Entrega |
|---|---|
| 1 | Fundação e núcleo kepleriano: CI, ADR, constantes, elementos ↔ estado, propagador kepleriano, tempo (JD/GMST), golden tests |
| 2 | Perturbação J2 (elementos médios), quadros ECI/ECEF, geodésia e projeto de órbita |
| 3 | Propagação numérica (Cowell, arrasto) e validação contra SGP4/TLE da ISS |
| 4 | Sol, observador e visibilidade da ISS sobre São Bento do Sapucaí |
| 5 | Constelação heliossíncrona de 6 satélites (91/6, um plano) e cobertura |
| 6 | Aplicação web e deploy |
| 7 | Extensão para 18 satélites e acabamento (se houver tempo) |
