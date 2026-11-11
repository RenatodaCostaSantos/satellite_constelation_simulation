# Fechamento da Semana 5 — Constelação heliossíncrona 91/6 de 6 satélites

Branch `semana-5`. Os seis cards (S5-01 a S5-06) foram concluídos e aprovados. Este documento segue
a seção 7 do prompt da Semana 5.

## 1. O que foi entregue

A constelação "Epic of Sun" do Sem4 §8 existe no simulador e roda no Cowell completo. São seis
satélites num único plano, separados de 60°, numa órbita que repete o traço a cada 91 voltas em
6 dias. A afirmação central do projeto foi medida, e não só calculada pela teoria: em 1 dia, os
seis satélites cobrem as 91 trilhas da grade, sempre no mesmo horário solar. O custo de manter a
constelação também está quantificado.

| Card | Entrega | Commits |
|---|---|---|
| S5-01 · Satellite, Constellation e a constelação 91/6 | `satellite.py` (`Satellite`), `constellation.py` (`Constellation` com saída (P, N, 3), `sso_constellation` com a receita da Semana 4), `orbit_design.repeat_orbit` (Sem4 Tabela 3) e `constellation_phases` (eq. 42); `build_scene` aceita `Satellite` e `Constellation` | `e9e4988`, `d547004` |
| S5-02 · Nós ascendentes, grade de 91 trilhas e repetição calibrada | `groundtrack.node_crossings`, `grid_metrics`, `repeat_closure`; `orbit_design.repeat_semi_major_axis_calibrated` (Δa = +19,37 m); `scripts/constellation_report.py` e `reports/s5_constelacao/` | `3bcc86b`, `fafffec`, `ff0573d` |
| S5-03 · Cobertura diária e revisita | `coverage.py`: `swath_half_angle` (eq. 45), `roll_for_half_grid`, `TargetGrid`, `access_times` (vetorizado por faixa de latitude, opção de luz), `coverage_stats`, `latitude_profile`; mapa de revisita e cobertura por latitude no relatório | `e58a2b7`, `2bd04e2` |
| S5-04 · Plano de contatos com a estação | `passes.contact_plan` (`ContactPlan`, trens, dia/noite pelo crepúsculo civil, sobreposições), `sun_elevation_deg`, `format_contact_plan`; `scripts/list_passes.py --constellation` | `9058782`, `8eec46b` |
| S5-05 · Manutenção: deriva do traço e do LTAN | `maintenance.py`: `track_drift`, `time_to_deadband`, Δv por manobra e por ano, `ltan_drift_from_inclination`, `ltan_offset`, `time_to_ltan_band`, `measure_decay_rate`, `MaintenanceStatus` | `e8dca05`, `09d5d3e` |
| S5-06 · GUI da constelação | Demo com os 6 satélites no Cowell completo (`--days` 1–6); na cena, `meta.constellation`, `nodes`, `swath_deg` e `maintenance` (opcionais); no viewer, legenda agrupada, camadas "grade do dia" e "faixa", painel com fase, correção e Δv, e velocidade 5000x | `e423447` |

Também na branch: `4848693` e `099624e` (resumo dos cards em `docs/cards_summary/`) e as seções
da Semana 5 em `docs/codigo_explicado`.

**Validação principal (contra o Apêndice A do prompt):**

| O quê | Resultado | Referência / critério |
|---|---|---|
| `repeat_orbit`, Sem4 Tabela 3 (6 órbitas) | h a ≤ 0,005 km, i a ≤ 0,0001° | 0,01 km / 0,001° |
| Atraso no nó entre S1 e S2 (Cowell completo) | 949,53 s | 949,5 ± 2 s |
| Construir os 6 satélites (i dada) | 0,33 s | < 3 s (A.1: 0,33 s) |
| Grade P = 6 em 1 dia (Cowell completo) | 91 nós, 3,95173…3,96432°, afastamento 0,0064° | A.3: 3,95173…3,96431°, 0,0064° |
| Grade P = 7 na mesma órbita | 106 nós, afastamento 1,694° | > 1° (A.3: 1,694°) |
| Fechamento de 6 dias, J2–J4, a de projeto / calibrado | 0,996 km / −0,022 km | A.4: 1,00 km / −0,022 km |
| Semieixo de repetição calibrado | 6 888 108,374 m (Δa +19,37 m) | A.4: +19,4 ± 1 m |
| Rolagem para meia grade | 23,1957° | A.5: 23,1957° |
| Cobertura P = 6, 1 dia, η = 23,2° | 100% em \|lat\| ≤ 80° | 100% |
| Revisita P = 6, 3 dias, \|lat\| ≤ 60° (sem / com luz) | 13,97 h / 24,26 h | ≤ 15 h / ≤ 25 h |
| Contatos com São Bento do Sapucaí, 24 h | 16 passagens idênticas ao A.6, 96,68 min, 0 sobreposições, 2 trens | A.6 |
| ȧ da constelação (Cowell completo) | −23,95 m/dia (1 dia), −24,01 (6 dias) | A.7: −23,71 (rel 2%) |
| Deriva do traço, Cowell × fórmula, 6 dias | 3,673 × 3,772 km (2,7%) | ≤ 10% (A.7: 3,68 × 3,72) |
| LTAN dos 6 satélites na cena, 24 h | 22,4994…22,5029 h, dispersão 0,0006 h | 22,500 ± 0,01 h |

## 2. Testes e lint

`ruff check .`:

```text
All checks passed!
```

`ruff format --check .`: `89 files already formatted`.

`pytest --cov=satsim`: **568 passed, 1 xfailed** em 191,9 s com a medição de cobertura (117,1 s sem ela), contra 463 testes na Semana 4. O xfail é a discrepância conhecida do Sem4 §1, registrada na Semana 1. Os 2 avisos são `DeprecationWarning` internos do Skyfield com o NumPy 2.5, não do nosso código. Os módulos novos da Semana 5 têm 100% de cobertura, exceto `orbit_design.py` (94%: ramos de erro do solver) e `passes.py` (99%).

```text
Name                                   Stmts   Miss  Cover
----------------------------------------------------------
src/satsim/__init__.py                     1      0   100%
src/satsim/analysis.py                    43      0   100%
src/satsim/app/__init__.py                 0      0   100%
src/satsim/app/scene_export.py           175      0   100%
src/satsim/astrotime.py                   32      0   100%
src/satsim/constants.py                   40      0   100%
src/satsim/constellation.py               72      0   100%
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
TOTAL                                   2321     32    99%
============ 568 passed, 1 xfailed, 2 warnings in 191.93s (0:03:11) ============
```

## 3. Discrepâncias encontradas (prompt × código)

Todas estão registradas no golden (`tests/golden/reference_values.py`) com os valores medidos.

1. **ω⊕ do Sem4 × rotação real do ECEF (S5-02).** O simulador gira o ECEF com o GMST real,
   7,2921158553e-5 rad/s. A teoria do Sem4 e o card usam ω⊕ = 7,292115e-5 rad/s, 8,55e-12 rad/s
   menor. Em 6 dias isso desloca a grade em 2,54e-4° (≈ 28 m). Daí vêm o passo "não uniforme" de
   0,000285° no modelo médio com P = 1 e a maior parte do resíduo de −22 m do a calibrado (A.4).
   A calibração seguiu o card e usa o ω⊕ do Sem4.
2. **Âncora do afastamento da grade (S5-02).** O A.3 dá 0,0064° com a grade ancorada no 1º nó de
   S1. Ancorada no 1º nó no tempo, que é o de S3, dá 0,0107°. Os desvios têm padrão de período 3
   (S1/S4, S2/S5, S3/S6; ±0,0045°): é o termo de curto período de J2 em Ω e u, que a conversão
   médio → osculador da S4-04 não corrige. Fica a âncora em S1, como no A.3.
3. **A.5, P = 1 em 1 dia (S5-03).** O satsim mede 48,82% de cobertura contra 48,48% do A.5 (grade
   de 1°), dentro do ±2% do card.
4. **A.5, definição de "intervalo entre visitas" (S5-03).** O A.5 se reproduz medindo de um
   **início** de passagem ao início seguinte (`interval="revisit"`, o padrão): 13,97 h, 19,04 h e
   24,26 h. Do **fim** de uma passagem ao início da seguinte (`interval="blackout"`) dá 13,69 h,
   16,68 h e 24,00 h. A diferença é a duração da passagem, que em latitude alta junta vários
   satélites do trem (15,8 min entre eles). As duas definições estão implementadas.
5. **A.7, ȧ medido (S5-05).** No Cowell completo da constelação, ȧ = −23,95 m/dia em 1 dia e
   −24,01 em 6 dias, contra −23,71 do A.7: 1,0 a 1,3% acima, dentro do rel 2% do card. Com esse ȧ,
   a próxima correção do traço cai em 6,92 dias (A.7: 6,95).
6. **Deriva do traço, Cowell × fórmula (S5-05).** O fechamento de 6 dias com o a calibrado e
   arrasto deu 3,673 km no Cowell completo (A.4: 3,68) e 3,772 km na fórmula com o ȧ medido
   (A.7: 3,72). A diferença de 2,7% está dentro dos 10% do card. A aceleração do LTAN medida no
   Cowell J2–J4 + Sol (60 dias) foi −7,13 min/ano², contra −6,88 da fórmula (4%).
7. **Faixa na cena (S5-06).** A faixa usa a rolagem de meia grade calculada no a de projeto
   (23,1957°). Com o raio médio real da órbita, cerca de 4 km acima de a, λ sai 1,9955° em vez de
   π/91 = 1,9780°. Por isso as faixas vizinhas se sobrepõem um pouco no mapa, em vez de só se
   encostarem.

## 4. Dívida técnica e pontos frágeis

| Ponto | Situação | Impacto / o que fazer |
|---|---|---|
| Tempo dos testes com 6 satélites | suíte em 117,1 s (191,9 s com cobertura), contra 84,8 s na Semana 4. A constelação de referência no Cowell completo é fixture de sessão (`conftest.py`) e serve S5-02 a S5-06. Mais lentos: fixture de taxas do `test_zonal` (16,0 s, Semana 4), ajuste do coeficiente balístico (11,8 s, Semana 3), Cowell de 30 dias do `test_mean_zonal` (10,9 s), grade com P = 7 no Cowell (8,7 s, S5-02) e a construção da constelação de referência, que roda uma vez por sessão (5,9 s, no primeiro teste que a usa) | o limite passou de ~2 para **~3 min** (decisão de 06/11/2026). A Semana 6, com 18 satélites, deve construir a constelação de 18 uma vez por sessão (nós de 1 dia ≈ 20 s no Cowell) |
| Esfera na cobertura | `access_times` usa a Terra esférica de raio R⊕ e o ângulo geocêntrico da eq. 45, como no Sem4 | o achatamento (21 km entre os raios equatorial e polar) não entra no alcance; o erro não foi medido. Basta para a cobertura diária com folga de rolagem; usar altitude geodésica e o alcance real da câmera se a precisão pedir |
| Manobras não simuladas | o `MaintenanceStatus` prevê quando e quanto (≈ 4,8 m/s/ano de traço + 4,75 de inclinação), mas o Cowell nunca corrige | na demo de 6 dias, a contagem da próxima correção não é zerada: passado o prazo, o painel mostra "vencida há X dias". Simular manobras impulsivas no Cowell fica fora do escopo |
| Tamanho da cena de 6 dias | 11,2 MB com passo de 60 s e 8 objetos (1 dia, 20 s: 6,7 MB); viewer interpola linearmente | com 18 satélites, a cena estimada de 1 dia a 20 s passa de 15 MB; a flecha da interpolação linear é 3,8 km a 60 s e 15 km a 120 s. Reduzir casas decimais antes de aumentar o passo, ou interpolar por Hermite |
| Termo de curto período de J2 | desvios de ±0,0045° da grade com período 3 (seção 3, item 2) | sem efeito na cobertura; corrigir exigiria incluir Ω e u na conversão médio → osculador |
| Banner do TLE da ISS | o aviso de TLE embarcado aparece mesmo com a ISS desligada na legenda | cosmético; pode virar condicional à visibilidade |
| Avisos do Skyfield | 2 `DeprecationWarning` internos do Skyfield 1.55 com o NumPy 2.5 em `test_passes.py` (Semana 3) | não é código nosso; versões inalteradas desde 29/09 |

## 5. Alterações em código de semanas anteriores

Todas compatíveis: nenhum teste antigo foi removido nem teve a tolerância afrouxada.

- `app/scene_export.py`:
  - `build_scene` aceita `Satellite`, `SceneSatellite` ou uma `Constellation` (S5-01). O
    `SceneSatellite` provisório continua funcionando.
  - Novos campos opcionais da cena, com o novo parâmetro `constellation_info` (S5-06):
    `meta.constellation`, `nodes`, `swath_deg` e `maintenance`. O `schema_version` continua 1.
- `groundtrack.py` (S5-02): acréscimo de `node_crossings`, `grid_metrics` e `repeat_closure`.
  `subpoint` e `ascending_node_crossings` não mudaram.
- `orbit_design.py` (S5-01, S5-02): acréscimo de `repeat_orbit`, `constellation_phases` e
  `repeat_semi_major_axis_calibrated`.
- `passes.py` (S5-04): acréscimo de `contact_plan` e tipos associados, `sun_elevation_deg` e
  `format_contact_plan`. `find_passes` e `find_passes_many` não mudaram.
- `scripts/build_demo.py` (S5-06):
  - A SSO-1 (`sso_design`, `design_satellites`) deu lugar aos 6 satélites de `sso_constellation`.
  - Novo `--days` de 1 a 6, com o passo definido em `STEP_BY_DAYS`.
  - `full_force_models` foi mantido, porque o `list_passes.py` o usa.
  - Os testes antigos do build_demo foram adaptados: ids `s1…s6` no lugar de `sso1`, e a
    constelação foi substituída por uma barata (J2 médio) quando o teste só verifica TLE ou
    tamanho.
- `scripts/list_passes.py` (S5-04): nova opção `--constellation`. A saída padrão não mudou.
- `app/static/viewer.js`, `viewer.html`, `viewer.css` (S5-06):
  - legenda agrupada, camadas, painel da constelação e velocidade 5000x;
  - o footprint e a órbita prevista aparecem só para o satélite selecionado quando há
    constelação;
  - arrastar a linha do tempo não redesenha mais a cada evento.
  - Cenas sem `meta.constellation` se comportam como antes.
- `tests/conftest.py`: fixture de sessão `reference_constellation` (Cowell completo, A.1).
- `README.md` (Demo visual, seção "Constelação 91/6") e `docs/adr/0001-arquitetura.md`
  (`Satellite` e `Constellation` implementados, contrato da cena).

## 6. O que a Semana 6 recebe pronto

**De 6 para 18 satélites.** O código já é paramétrico em P: `sso_constellation(epoch, models,
n_sats=18, a=..., i=...)` funciona hoje. O que ainda está fixo em 6 ou 91 são as constantes de
`scripts/build_demo.py` e de `constellation_report.py` (`N_REVS`, `N_DAYS`, nomes), a paleta de 6
cores e os testes.

O Sem4 §8.3 generaliza o teorema: com P = m·D satélites num plano, o padrão continua uniforme.
Para a órbita 91/6, 18 = 3 × 6 e Δu = 20°. Um experimento rápido (Cowell completo, 2026-10-23,
detalhado no Apêndice A do prompt da Semana 6) mostra o que isso significa na prática:

- os 18 depositam **273 trilhas uniformes por dia** (1,31556° a 1,32196°, afastamento 0,007°),
  não 3 passagens sobre as mesmas 91. O vizinho cruza o equador 316,5 s depois, 1/3 de passo a
  oeste. Os satélites s1, s4, …, s16 são exatamente a constelação de 6;
- a cobertura diária completa passa a pedir só **±8,18° de rolagem**, no lugar de 23,2°;
- a revisita quase não muda: 13,79 h contra 13,97 h com 23,2°. As oportunidades extras chegam no
  mesmo trem, minutos depois, porque num só plano todas as passagens têm o mesmo horário solar;
- na estação, são 49 passagens por dia, 36 sobreposições e até 2 satélites simultâneos: seriam
  necessárias duas antenas.

A alternativa, três planos com LTAN diferentes, daria revisita em horários distintos, mas abre mão
do horário solar único. Fica para decisão; o prompt da Semana 6 adota o plano único do Sem4.

**O que falta para o deploy web:**

- A demo já é um HTML autocontido que funciona offline (globe.gl, texturas e cena embutidos) e
  pode ser servida como site estático sem mudanças.
- Faltam a escolha entre 6 e 18 satélites no viewer (dois arquivos irmãos, porque juntos passariam
  de 15 MB), uma página inicial, um script que monte o site e um workflow do GitHub Pages que rode
  os testes antes de publicar. Tudo isso está nos cards S6-02 e S6-03 de
  `docs/PROMPT_SEMANA6_claude_code.html`.
- A ativação do Pages no repositório é manual.

**Lembrete:** o `git push` da branch `semana-5` e a abertura do PR ficam com você.
