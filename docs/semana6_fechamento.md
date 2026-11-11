# Fechamento da Semana 6 — Constelação de 18 satélites, escolha 6/18 e site na web

Branch `semana-6`, criada a partir de `semana-5` (`1dbd2e8`). Os três cards (S6-01 a S6-03) foram
concluídos e aprovados. Este documento segue a seção 7 do prompt da Semana 6 e fecha também o
projeto (Semanas 1–6).

## 1. O que foi entregue

O requisito R7 do ADR está concluído: o mesmo viewer mostra 6 ou 18 satélites, que entram só por
dados. O simulador está pronto para ir à web como site estático, publicado por um workflow do
GitHub Pages. A afirmação da semana foi medida no Cowell completo, e não só tirada da teoria:
com 18 satélites num plano, o traço forma **273 trilhas por dia**, e não 3 passagens sobre as
mesmas 91, e a cobertura diária pede **±8,18°** de rolagem, no lugar de 23,2°. A revisita quase
não muda.

| Card | Entrega | Commit |
|---|---|---|
| S6-01 · Constelação de 18 satélites: P = 3·D no Cowell completo | `ConstellationDesign` (N, D, P), com as grandezas do Sem4 §8.2–8.3 (trilhas por dia, grade diária, famílias, rolagem para cobertura diária, `build`); presets `DESIGN_91_6_P6` e `DESIGN_91_6_P18`; paleta por família (com P = 6, igual à da S5); `constellation_report.py --sats 18` e `reports/s6_constelacao18/` | `864d003` |
| S6-02 · Visualizador com escolha 6 ou 18 satélites | `build_demo.py --sats {6,18}` e `--all`, com quatro arquivos irmãos ligados por `meta.variants`; seletor no viewer que preserva o instante, a velocidade, as camadas e a seleção; legenda por família; contador "X de 273"; `ScenePrecision` (precisão compacta só com 18); interpolação Hermite no viewer; modo `#bench` | `fae7bcb` |
| S6-03 · Site estático e deploy web | `scripts/build_site.py`: página inicial leve com quatro botões e as cenas em `site/cenas/`; `.github/workflows/pages.yml` (test → deploy, só `actions/*`, sem segredos); `tests/test_build_site.py`; README (publicação), ADR (Deploy), resumo dos cards | `3bd4219` |

Nenhum card precisou de "ajustes de revisão".

**Validação principal (contra o Apêndice A do prompt da Semana 6):**

| O quê | Resultado | Referência |
|---|---|---|
| Construir os 18 satélites (i dada) | ≈ 1 s | A.1: 0,98 s |
| Atraso no nó entre S1 e S2 | 316,54 s (vizinhos: 316,42 a 316,54 s) | A.1: 316,54 s (Tnod/18) |
| Primeiros nós s1…s4 | 4 274,46; 4 591,00; 4 907,52; 5 223,99 s | A.1: os mesmos |
| Grade P = 18, 1 dia | 273 nós, 1,31556…1,32196°, afastamento 0,00695° | A.1: os mesmos |
| Família 0 (s1, s4, …, s16) | 91 nós, 3,95173…3,96432°, afastamento 0,00645° (= constelação de 6) | A.1: os mesmos |
| Rolagem para cobertura diária (eq. 45) | 23,1957° (P = 6) e 8,1835° (P = 18) | A.2: os mesmos |
| Cobertura P = 18, η = 8,18°, 1 dia | 100% em \|lat\| ≤ 80° | A.3: 100% |
| Cobertura P = 6, η = 8,18°, 3 dias | 74,35%, intervalo de 72,0 h | A.3: 74,35%; 72,0 h |
| Revisita, 3 dias, η = 23,2°, \|lat\| ≤ 60° | 13,79 h (P = 18) e 13,97 h (P = 6) | A.3: os mesmos |
| Passagens por alvo em 1 dia (P = 18, η = 23,2°) | mín. 2, mediana 2, máx. 3 | A.3: 2 / 2 / 3 |
| Contatos com São Bento do Sapucaí, 24 h | 49 passagens, 287,99 min, 2 trens (24 e 25), 36 sobreposições, até 2 simultâneos, maior intervalo 9,88 h | A.4: os mesmos |
| ȧ dos 18 satélites | −24,02 a −23,95 m/dia | o mesmo da S5 (−23,95) |
| Cenas da demo (`--all`) | 6,7 / 11,2 / 11,4 / 11,6 MB | cada uma < 15 MB |
| Interpolação do viewer, erro máximo de posição | 0,07 km (30 s), 0,28 km (60 s), 1,13 km (120 s) | ≤ 5 km |
| Site (`build_site.py --offline`) | página inicial 6 kB; total 41,3 MB | página < 1 MB; cenas < 15 MB |

## 2. Testes e lint

`ruff check .`:

```text
All checks passed!
```

`ruff format --check .`: `93 files already formatted`.

`pytest --cov=satsim`: **608 passed, 1 xfailed** em 269,8 s com a medição de cobertura (151,7 s sem ela; o limite de ~3 min vale para a suíte sem cobertura), contra 568 testes na Semana 5. Cobertura total de **99%**. O xfail é a discrepância conhecida do Sem4 §1, registrada na Semana 1. Os 2 avisos são `DeprecationWarning` internos do Skyfield com o NumPy 2.5, não do nosso código. Na Semana 6 entraram 40 testes: 28 em `test_constellation18.py`, 7 em `test_scene_export.py` e 5 em `test_build_site.py`. Os módulos alterados na semana (`constellation.py` e `app/scene_export.py`) têm 100% de cobertura.

```text
============================= test session starts ==============================
platform linux -- Python 3.12.3, pytest-9.1.1, pluggy-1.6.0
rootdir: /home/kunumi/projetos/constelacao_satelites/satellite_constelation_simulation
configfile: pyproject.toml
testpaths: tests
plugins: cov-7.1.0
collected 609 items

tests/test_astrotime.py ..............                                   [  2%]
tests/test_atmosphere_drag.py .......................................... [  9%]
                                                                         [  9%]
tests/test_build_site.py .....                                           [ 10%]
tests/test_constants.py ..............                                   [ 12%]
tests/test_constellation.py .......................................      [ 18%]
tests/test_constellation18.py ............................               [ 23%]
tests/test_contact_plan.py .............                                 [ 25%]
tests/test_coverage.py .................                                 [ 28%]
tests/test_cowell.py ...........................                         [ 32%]
tests/test_elements.py ..........................                        [ 36%]
tests/test_ephemeris_third_body.py ..................                    [ 39%]
tests/test_frames_geodesy.py .............                               [ 42%]
tests/test_golden_consistency.py ..................................x.... [ 48%]
..                                                                       [ 48%]
tests/test_groundtrack.py ..................                             [ 51%]
tests/test_j2_force.py .....                                             [ 52%]
tests/test_kepler_propagator.py ................                         [ 55%]
tests/test_maintenance.py ............                                   [ 57%]
tests/test_mean_j2.py ...........                                        [ 58%]
tests/test_mean_osc.py ........                                          [ 60%]
tests/test_mean_zonal.py .................                               [ 63%]
tests/test_observer.py ...............                                   [ 65%]
tests/test_orbit_design_sso.py ......................                    [ 69%]
tests/test_orbit_determination.py ...                                    [ 69%]
tests/test_passes.py .................                                   [ 72%]
tests/test_repeat_grid.py ..............                                 [ 74%]
tests/test_scene_export.py ............................................  [ 81%]
tests/test_secular.py .......................                            [ 85%]
tests/test_sgp4_propagator.py .............                              [ 87%]
tests/test_smoke.py .                                                    [ 88%]
tests/test_srp_eclipse.py .......                                        [ 89%]
tests/test_tle.py ...............................................        [ 96%]
tests/test_zonal.py ...................                                  [100%]

=============================== warnings summary ===============================
tests/test_passes.py::test_matches_skyfield_find_events[10.0]
tests/test_passes.py::test_matches_skyfield_find_events[5.0]
  .venv/lib/python3.12/site-packages/skyfield/searchlib.py:224: DeprecationWarning: Setting the dtype on a NumPy array has been deprecated in NumPy 2.5.
  Instead of changing the dtype on an array x, create a new array with x.view(new_dtype)
    y.dtype = int8

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
================================ tests coverage ================================
_______________ coverage: platform linux, python 3.12.3-final-0 ________________

Name                                   Stmts   Miss  Cover
----------------------------------------------------------
src/satsim/__init__.py                     1      0   100%
src/satsim/analysis.py                    43      0   100%
src/satsim/app/__init__.py                 0      0   100%
src/satsim/app/scene_export.py           224      0   100%
src/satsim/astrotime.py                   32      0   100%
src/satsim/constants.py                   40      0   100%
src/satsim/constellation.py              154      0   100%
src/satsim/coverage.py                   146      0   100%
src/satsim/decay.py                       39      1    97%
src/satsim/elements.py                   114      1    99%
src/satsim/ephemeris.py                   61      0   100%
src/satsim/forces/__init__.py              7      0   100%
src/satsim/forces/atmosphere.py           31      2    94%
src/satsim/forces/drag.py                 40      1    98%
src/satsim/forces/j2.py                   26      0   100%
src/satsim/forces/srp.py                  34      1    97%
src/satsim/forces/third_body.py           28      0   100%
src/satsim/forces/zonal.py                56      0   100%
src/satsim/frames.py                      35      0   100%
src/satsim/geodesy.py                     36      1    97%
src/satsim/groundtrack.py                 53      0   100%
src/satsim/maintenance.py                 78      0   100%
src/satsim/mean_osc.py                   101      0   100%
src/satsim/observer.py                    65      0   100%
src/satsim/orbit_design.py               127      7    94%
src/satsim/orbit_determination.py         65      4    94%
src/satsim/passes.py                     203      3    99%
src/satsim/propagators/__init__.py         7      0   100%
src/satsim/propagators/base.py            20      0   100%
src/satsim/propagators/cowell.py         100      1    99%
src/satsim/propagators/kepler.py          34      2    94%
src/satsim/propagators/mean_j2.py         34      0   100%
src/satsim/propagators/mean_zonal.py      78      4    95%
src/satsim/propagators/sgp4_ref.py        47      0   100%
src/satsim/satellite.py                   22      0   100%
src/satsim/secular.py                     76      2    97%
src/satsim/tle.py                        195      2    99%
----------------------------------------------------------
TOTAL                                   2452     32    99%
============ 608 passed, 1 xfailed, 2 warnings in 269.80s (0:04:29) ============
```

## 3. Discrepâncias encontradas (PDFs/prompt × código)

Os números novos estão no golden (`tests/golden/reference_values.py`), com a fonte.

1. **Sem4 §8.3: "m passagens por dia em cada trilha" (S6-01).** O Sem4 diz que, com P = m·D,
   cada trilha recebe m passagens por dia. Para 91/6 com P = 18 (m = 3), seriam 3 passagens
   sobre cada uma das 91 trilhas. No Cowell completo, o vizinho cruza o equador 316,5 s depois
   (Tnod/18) e, nesse tempo, a Terra gira um terço de passo de grade (eq. 43: D/P = 1/3). Por
   isso os 18 depositam **273 trilhas diferentes**, com 1,31556 a 1,32196° de espaçamento
   (uniforme: 1,31868°), uma passagem em cada. O "m" do teorema multiplica a densidade da grade,
   não as passagens por trilha. As três famílias f = j mod 3 são três grades de 91 intercaladas
   (afastamentos 0,00645°, 0,00785° e 0,00784°).
2. **Passagens por dia com P = m·D (S6-01).** Num alvo, com η = 23,2° e \|lat\| ≤ 60° em 1 dia:
   - fundindo acessos a menos de 30 min, são **2 / 2 / 3** passagens (mín. / mediana / máx.)
     tanto com 6 quanto com 18 satélites, como no A.3;
   - contadas satélite a satélite, são 2 / 2 / 5 com 6 e **6 / 7 / 13** com 18: três vezes mais
     oportunidades, mas no mesmo trem.

   Por isso a revisita não cai: 13,79 h com 18, contra 13,97 h com 6. A janela de um trem sobre
   o alvo, do 1º ao último acesso, dá mediana de **16,2 min**, contra "≈ 10 min" no A.3. A
   diferença é de definição: em latitude alta a passagem fundida junta mais satélites (máx.
   151,7 min).
3. **A.5, tamanho da cena (S6-02).** A estimativa por proporção do A.5 ignora a parte fixa de
   cada arquivo, ≈ 2,3 MB de globe.gl e textura, que não cresce com o número de objetos.

   | Cena | Estimativa do A.5 | Precisão completa | Precisão compacta |
   |---|---:|---:|---:|
   | 18 satélites, 1 dia, passo de 20 s | 17 MB | 13,4 MB | 11,4 MB (adotada) |
   | 18 satélites, 6 dias, passo de 60 s | 28 MB | 24,5 MB | 20,6 MB |
   | 18 satélites, 6 dias, passo de 120 s | 14 MB | — | 11,6 MB (adotada) |

   A cena de 6 satélites manteve a precisão completa, para continuar idêntica à da S5.
4. **A.5, flecha da interpolação (S6-02).** A flecha da corda r·θ²/8 (0,9 / 3,8 / 15 km a 30 /
   60 / 120 s) vale para uma corda em ECI. O viewer das Semanas 3–5 interpolava linearmente em
   latitude/longitude, e perto da latitude máxima o erro real era muito maior: **7,37 / 29,16 /
   112,1 km**. Com a Hermite cúbica no vetor unitário, adotada em todas as cenas, o erro cai
   para 0,070 / 0,279 / 1,127 km. O máximo de 1,127 km a 120 s vem das tangentes de um lado só
   nas bordas da janela.
5. **Página inicial e rede no CI (S6-03).** O card pede que "a página inicial avise" quando a
   CelesTrak falhar. `build_site.py` avisa sempre que o TLE da ISS não veio da rede nesta
   montagem, seja do cache, seja embutido. O `sentinel_2a()` da S4 não expõe a fonte do TLE dele,
   então a página só diz se o Sentinel-2A entrou ou não na cena.
6. **Verificação visual limitada (S6-02, S6-03).** O WebGL por software do Chrome headless deste
   ambiente não desenha a Terra no globo, nem com o viewer antigo. O seletor, o relógio, a
   legenda e o mapa 2D foram verificados no headless, por `file://` e por HTTP. O globo com 18
   marcadores ficou para a sua verificação visual, e a imagem do README é a composição dos dois
   mapas 2D.

## 4. Alterações em código de semanas anteriores

Todas compatíveis. Nenhum teste antigo foi removido nem teve a tolerância afrouxada, e com P = 6
os resultados da Semana 5 continuam idênticos (`test_s6_scene_6_identical_except_new_optional_fields`).

- `constellation.py` (S6-01): só acréscimos (`ConstellationDesign`, presets,
  `constellation_palette`). `sso_constellation` e `CONSTELLATION_PALETTE` não mudaram.
- `app/scene_export.py` (S6-02):
  - parâmetros novos `variants` e `precision`, com padrão `FULL_PRECISION`, igual às casas da S5;
  - `SceneConstellation` ganhou `families` e `from_design`;
  - campos opcionais novos: `meta.variants`, `tracks_per_day`, `daily_grid_deg`, `families` e
    `family`. O `schema_version` continua 1.
- `app/static/viewer.js`, `viewer.html`, `viewer.css` (S6-02):
  - interpolação Hermite no lugar da linear em lat/lon, o que melhora também as cenas antigas;
  - estado na URL (`#t=…&v=…`), seletor de variantes, legenda por família, contador de trilhas
    lido da cena (antes fixo em 91);
  - modo `#bench=N`.
  - Cenas sem os campos novos abrem como antes.
- `scripts/build_demo.py`:
  - `build_constellation(start, n_sats=6)` e `constellation_info` via `from_design` (S6-02);
  - `main` separado em `build_demo()` (S6-02) e `build_all()` (S6-03);
  - novos `--sats`, `--all` e `--out-dir`. Sem argumentos, gera a mesma `demo/index.html`.
- `scripts/constellation_report.py` (S6-01): novo `--sats {6,18}`. O padrão, 6, gera o
  relatório da S5 sem mudança.
- `tests/test_scene_export.py`: os `monkeypatch` antigos de `build_constellation` passaram a
  aceitar `n_sats`, e `_fast_constellation` foi generalizado. Só a assinatura mudou.
- `tests/conftest.py`: fixtures de sessão da constelação de 18 (construção e nós de 1 dia).
- Documentação: README (Demo, "6 ou 18 satélites", "Site na web"), ADR (status, R7, contrato,
  stack, Deploy), `.gitignore` (`site/`), `docs/cards_summary/cards_summary.html` (agora da
  Semana 6).

## 5. Fechamento do projeto (Semanas 1–6)

### O que o simulador faz (para leigos)

O **satsim** responde a uma pergunta prática: se colocarmos uma constelação de satélites de
observação da Terra numa certa órbita, onde cada um vai estar, o que vai enxergar e quando vai
falar com a estação de solo?

Para isso, ele calcula o movimento dos satélites com as forças que importam a 500 km de altura:
- a Terra não é uma esfera perfeita (é achatada nos polos e tem uma leve forma de pera), e isso
  gira o plano da órbita;
- o Sol e a Lua puxam o satélite;
- a luz do Sol empurra os painéis;
- o resto da atmosfera freia o satélite e o faz descer cerca de 24 m por dia.

A órbita escolhida é **heliossíncrona**: o satélite cruza o equador sempre à mesma hora solar
(22h30 na subida), e por isso as imagens têm sempre a mesma iluminação. Ela também é de
**repetição 91/6**: depois de 91 voltas, em exatos 6 dias, o satélite repassa sobre o mesmo
caminho no solo.

Com **6 satélites** em fila no mesmo plano, as 91 trilhas do ciclo são percorridas todos os
dias. Com **18**, são 273 trilhas por dia, e a câmera precisa inclinar bem menos para ver
qualquer lugar. O simulador também mostra o que **não** melhora: o tempo entre duas fotos do
mesmo lugar continua em cerca de 14 h, porque os satélites passam juntos, como vagões de um
trem.

Para conferir se as contas estão certas, o simulador compara seus resultados com satélites
reais: a **Estação Espacial Internacional** e o **Sentinel-2A** europeu, cujas posições oficiais
(TLEs) são públicas. Ele calcula as passagens sobre a estação de **São Bento do Sapucaí** (SP),
a cobertura da superfície e o combustível para manter a órbita (≈ 4,8 m/s por ano por satélite).

Tudo pode ser visto num **navegador**: um globo 3D e um mapa sincronizados, com relógio,
velocidades de 1x a 5000x, a grade de trilhas se formando e a troca entre 6 e 18 satélites no
mesmo instante. O navegador não calcula nada; ele só desenha o que o Python já calculou. Por
isso o mesmo arquivo abre com duplo clique, sem internet, ou num site estático.

### Validações principais

| O quê | Resultado | Critério / referência | Semana |
|---|---|---|---|
| SGP4 × vetores de regressão (ISS) | iguais a 1 m e 1e-3 m/s | 1 m / 1e-3 m/s | S3 |
| ISS: SGP4 → ECEF × Skyfield, 3 dias | máx. 85,6 m | < 150 m | S3 |
| ISS: azimute / elevação / distância × Skyfield | 0,026° / 0,009° / 74 m | < 0,05° / 0,05° / 150 m | S3 |
| ISS: passagens sobre São Bento × Skyfield | ≤ 0,19 s | < 1 s | S3 |
| ISS: Cowell (J2 + arrasto ajustado) × SGP4 | 0,3 / 2,7 / 14,5 km em 24 h / 72 h / 7 d | ≤ 2 / 8 / 40 km | S3 |
| ISS: Cowell completo ajustado × SGP4 | 0,03 / 0,06 / 0,49 km em 24 h / 72 h / 7 d | ≤ 0,5 / 1 / 8 km | S4 |
| Sentinel-2A: Cowell completo ajustado × SGP4 | 0,01 / 0,02 / 0,36 km em 24 h / 72 h / 7 d | ≤ 0,2 / 0,3 / 3 km | S4 |
| Sentinel-2A: inclinação heliossíncrona calibrada | 98,56940° (real 98,56361°) | < 0,01° do real | S4 |
| Sol / Lua × DE421 (direção) | 0,014° / 0,085° | ≤ 0,02° / 0,1° | S4 |
| Repetição 91/6: fechamento de 6 dias (J2–J4, a calibrado) | −0,022 km | A.4 S5: −0,022 km | S5 |
| Grade de 91 trilhas (P = 6, Cowell completo) | afastamento 0,0064° | A.3 S5: 0,0064° | S5 |
| Grade de 273 trilhas (P = 18, Cowell completo) | afastamento 0,00695° | A.1 S6: 0,00695° | S6 |
| Cobertura diária: P = 6 com 23,2° e P = 18 com 8,18° | 100% / 100% | 100% | S5, S6 |
| Revisita, \|lat\| ≤ 60°, 3 dias | 13,97 h (P = 6) e 13,79 h (P = 18) | ≤ 15 h; A.3 S6 | S5, S6 |
| Contatos com São Bento, 24 h | 16 passagens, 96,68 min (P = 6); 49, 287,99 min (P = 18) | A.6 S5; A.4 S6 | S5, S6 |
| Manutenção: ȧ medido × A.7 | −23,95 × −23,71 m/dia (1,0%) | rel 2% | S5 |
| Manutenção: deriva do traço, Cowell × fórmula, 6 dias | 3,673 × 3,772 km (2,7%) | ≤ 10% | S5 |
| LTAN da constelação na cena, 24 h | 22,4994 … 22,5029 h | 22,500 ± 0,01 h | S5 |

A referência "real" da ISS e do Sentinel-2A é o SGP4, que tem erro próprio de cerca de 1 km na
época. Erros do Cowell abaixo disso não dão para distinguir (README, Limitações).

### Limitações conhecidas

- **Esfera na cobertura:** `access_times` usa a Terra esférica e o ângulo geocêntrico da eq. 45.
  O achatamento (21 km entre os raios equatorial e polar) não entra no alcance, e esse erro não
  foi medido.
- **Manobras não simuladas:** o `MaintenanceStatus` prevê quando e quanto corrigir (traço ≈ 4,8
  m/s/ano, inclinação ≈ 4,75 m/s/ano), mas o Cowell nunca aplica a correção. Na cena de 6 dias,
  a correção aparece "vencida" depois do prazo.
- **Interpolação do viewer:** as cenas são pré-amostradas (20 a 120 s). A Hermite erra até
  1,13 km a 120 s. Fora da janela exportada, é preciso gerar uma nova cena.
- **Tamanho das cenas:** de 6,7 a 11,6 MB por arquivo, 41 MB no site. Com mais objetos ou
  janelas mais longas, seria preciso aumentar o passo ou trocar o formato (binário, ou
  propagação no navegador).
- **Um único plano:** todas as passagens acontecem no mesmo horário solar. Por isso 18 satélites
  não encurtam a revisita, e a estação recebe trens com até 2 satélites simultâneos. O
  escalonamento da antena foi quantificado, mas não implementado.
- **Referência é o SGP4:** "mais correto" quer dizer mais perto do SGP4, que erra cerca de 1 km.
  Sol, Lua e SRP foram validados por componentes contra o DE421, não pela órbita.
- **Modelo de tempo e quadros:** sem precessão e nutação, UT1 ≈ UTC e atmosfera estática sem
  F10.7/Ap. O erro fica em centenas de metros, aceitável para visibilidade e projeto.

### Próximos passos possíveis

1. **Três planos com LTAN diferentes** (por exemplo, 10h30, 13h30 e 22h30 em vez de um plano com
   18): revisita em horários solares distintos, ao custo do horário solar único. O
   `ConstellationDesign` já separa N, D e P; faltaria o RAAN por plano.
2. **Manobras simuladas:** impulsos no Cowell quando o traço sai da banda de ±5 km, com o Δv
   acumulado comparado ao previsto pelo `MaintenanceStatus`.
3. **Apontamento real da câmera:** campo de visada e rolagem por alvo, com a Terra elipsoidal no
   alcance, no lugar da faixa ±λ em esfera.
4. **Escalonamento da estação:** escolher um satélite por vez nos trens com sobreposição, ou
   dimensionar duas antenas.
5. **Referência mais precisa que o SGP4:** efemérides POD/SP3 do Sentinel-2A como `Propagator`
   no `fit_initial_state`, sem mudar a determinação de órbita.

## 6. Lembrete

O `git push` da branch `semana-6`, a abertura do PR para a `main` e a ativação do GitHub Pages
(Settings → Pages → Source: "GitHub Actions") ficam com você. Depois do merge na `main`, o
workflow `pages.yml` roda os testes e publica o site em
<https://renatodacostasantos.github.io/satellite_constelation_simulation/>.
