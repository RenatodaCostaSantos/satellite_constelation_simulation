# Fechamento da Semana 3 — ISS real, arrasto, estação de solo e passagens

Branch `semana-3` (a partir da `main` com a Semana 2 mergeada, PR #3). Seis cards concluídos e
aprovados (S3-01 a S3-06). Este documento segue a seção 7 do prompt da Semana 3 e acrescenta, no
fim, a **revisão do objetivo do projeto** e as recomendações de planejamento discutidas na
revisão do S3-06.

## 1. O que foi entregue

| Card | Entrega | Commits |
|---|---|---|
| S3-01 · Ingestão de TLE | `satsim.tle`: parse e validação (checksum, comprimento, catálogo), download da CelesTrak, cache em disco com escrita atômica, fallback em camadas (cache < 2 h → rede → cache antigo → TLE embarcado), avisos de idade (> 3 e > 14 dias); `scripts/fetch_tle.py` | `d3c2687`, `430a334` |
| S3-02 · SGP4Propagator | `SGP4Propagator` (biblioteca `sgp4`, WGS72, saída TEME = ECI em SI, vetorizado, erros → `NaN`); `Propagator.propagate_utc` | `70728a1`, `79326b4`, `0a18b50` |
| S3-03 · Arrasto no Cowell | atmosfera exponencial por partes (Vallado), `DragAcceleration` co-rotante (altitude geodésica ou esférica), `decay.py` (taxa, perda por revolução, tempo de vida), reentrada no Cowell, `fit_ballistic_coefficient`; relatório `reports/s3_arrasto/` | `af6f505`, `d53d108` |
| S3-04 · GroundStation | `observer.py`: estação WGS84, ENU, azimute/elevação/distância, `observe`, `footprint_half_angle`; `SAO_BENTO_DO_SAPUCAI` | `fba66ec`, `0aca602` |
| S3-05 · Passagens | `passes.py`: AOS/TCA/LOS por varredura vetorizada + `brentq` + `minimize_scalar`, bordas da janela, vários satélites, tabela UTC/BRT; `scripts/list_passes.py` | `928a00b`, `cc2cdb4`, `9187843`, `27055be` |
| S3-06 · GUI | ISS real (TLE/SGP4) no viewer, painel "ISS × São Bento", lista de passagens clicável, contagem regressiva, linha estação → satélite, footprint no mapa, faixa do TLE; contrato da cena estendido de forma compatível | `b11ef63` |

Os commits "ajustes de revisão" incluem as seções 8 a 11 do `docs/codigo_explicado` (TLE, SGP4,
arrasto e o "V" do coeficiente balístico, sistemas de coordenadas) e comentários linha a linha
pedidos na revisão.

**Validação principal (todas contra o Apêndice A do prompt e o Skyfield 1.55):**

| O quê | Resultado | Referência / critério |
|---|---|---|
| SGP4 × vetores de regressão A.1 | iguais a 1 m e 1e-3 m/s | 1 m / 1e-3 m/s |
| SGP4 → ECEF × Skyfield (3 dias) | máx. 85,6 m, média 70 m | < 150 m (ref. 86 m) |
| Azimute / elevação / distância × Skyfield | 0,026° / 0,009° / 74 m | < 0,05° / 0,05° / 150 m |
| Passagens a 10° (3 dias) | 9, dentro de ±1 s / ±0,05° / ±0,1° | A.2 |
| Passagens × `find_events` do Skyfield | ≤ 0,19 s; 9 (10°) e 11 (5°) | < 1 s |
| Arrasto e decaimento (Sem3) | 1,2e-6 m/s² a 400 km; −11,8 m/rev; Tabela 3 ±5%; vidas 4,64 / 1,16 / 23,2 anos | Sem3 / A.4 |
| Cowell (J2 + arrasto ajustado) × SGP4 | 0,3 / 2,7 / 14,5 km em 24 h / 72 h / 7 d | ≤ 2 / 8 / 40 km |

## 2. Testes e lint

`ruff check .` e `ruff format --check .`:

```text
All checks passed!
59 files already formatted
```

`pytest --cov=satsim`: **389 passed, 1 xfailed** (o xfail é a discrepância conhecida do Sem4 §1,
registrada na Semana 1). Os 2 avisos são um `DeprecationWarning` interno do Skyfield com o
NumPy 2.5 (`searchlib.py`), não do nosso código.

```text
Name                                 Stmts   Miss  Cover
--------------------------------------------------------
src/satsim/__init__.py                   1      0   100%
src/satsim/analysis.py                  43      0   100%
src/satsim/app/__init__.py               0      0   100%
src/satsim/app/scene_export.py         106      0   100%
src/satsim/astrotime.py                 32      0   100%
src/satsim/constants.py                 32      0   100%
src/satsim/decay.py                     39      1    97%
src/satsim/elements.py                 114      1    99%
src/satsim/forces/__init__.py            4      0   100%
src/satsim/forces/atmosphere.py         31      2    94%
src/satsim/forces/drag.py               40      1    98%
src/satsim/forces/j2.py                 26      0   100%
src/satsim/frames.py                    35      0   100%
src/satsim/geodesy.py                   36      1    97%
src/satsim/groundtrack.py               27      0   100%
src/satsim/observer.py                  65      0   100%
src/satsim/orbit_design.py              34      0   100%
src/satsim/passes.py                    95      3    97%
src/satsim/propagators/__init__.py       6      0   100%
src/satsim/propagators/base.py          20      0   100%
src/satsim/propagators/cowell.py       100      3    97%
src/satsim/propagators/kepler.py        34      2    94%
src/satsim/propagators/mean_j2.py       34      0   100%
src/satsim/propagators/sgp4_ref.py      47      0   100%
src/satsim/secular.py                   30      0   100%
src/satsim/tle.py                      190      2    99%
--------------------------------------------------------
TOTAL                                 1221     16    99%
389 passed, 1 xfailed, 2 warnings in 44.31s
```

## 3. Discrepâncias encontradas (prompt × código)

Nenhuma tolerância foi afrouxada. Cada correção está no golden file com a dedução, e quando
possível há um teste de consistência que a teria pegado.

| Card | Prompt dizia | Valor correto | Como foi tratado |
|---|---|---|---|
| S3-02 | subponto com \|lat\| ≤ 51,7° (geodésica) | 51,7° vale para a latitude **geocêntrica** (51,62°); a geodésica é 51,80° (+ f·sin 2φ·R⊕/(R⊕ + h) ≈ 0,18°) | mesma correção da S2-07; teste verifica as duas |
| S3-02 | período 2π/n ≈ 5573,9 s | 5573,82 s | dentro de abs = 1 s; o prompt arredondou |
| S3-03 | altitude geodésica × esférica "até ±21 km" | sempre ≥ 0: de 0 (equador) a 21,4 km (polos); ≈ 13 km para a ISS | docstring de `drag.py` e relatório corrigidos |
| S3-03 | A.4 imprime \|dh/dt\| = 0,16 m/dia a 1000 km | 0,1556 (2 algarismos no prompt) | teste compara na precisão impressa |
| S3-04 | λ = arccos(R⊕/(R⊕ + 400 km)) ≈ 19,7° | 19,78° (truncado no prompt) | golden corrigido |
| S3-05 | AOS da passagem 9 do A.2 = 237 528,4 s | **237 538,4 s** (a própria linha: coluna UTC 10:37:28 e duração 132 s; Skyfield: 237 538,66 s) | golden corrigido + `test_passes_a2_internal_consistency` |

**Resultados que divergem da intuição, mas não do prompt (documentados):**

- **Arrasto padrão piora o Cowell** (C_D = 2,2, A/m = 0,004): 5,8 / 52 / 286 km contra 1,0 / 8,7 /
  47 km só com J2. Dezembro de 2019 foi mínimo solar; a tabela de Vallado (atividade moderada)
  põe ~5× mais ar a 400 km do que o B\* do TLE indica. O erro tem forma de "V" em função de
  C_D·A/m, com mínimo no valor ajustado (0,00163 m²/kg para o RMS em 72 h).
- **O mínimo depende da janela**: o RMS em 7 dias é mínimo em ≈ 0,0014, não em 0,00163. O
  C_D·A/m ajustado é um parâmetro *efetivo*, que absorve também J3, J4 e Sol/Lua, que o SGP4 tem e
  o Cowell não.

**Bugs encontrados e corrigidos durante a semana:**

- `meta.epoch_utc` era truncado ao segundo enquanto as amostras partiam do instante exato
  (0,36 s de defasagem com a época de um TLE). A época da cena agora é arredondada ao
  milissegundo antes da amostragem (S3-06).
- Painéis do viewer não atualizavam imediatamente após um salto do relógio (S3-06).

## 4. Dívida técnica e pontos frágeis

| Ponto | Situação | Impacto |
|---|---|---|
| Idade do TLE | avisos > 3 e > 14 dias; o TLE embarcado envelhece com o tempo | posição real erra km/dia com TLE antigo; atualizar `src/satsim/data/iss_fallback.tle` periodicamente |
| Atmosfera abaixo de 100 km | devolve ρ de 100 km com aviso | irrelevante (reentrada para em 100 km), mas é um patamar artificial |
| Atmosfera estática (Vallado, atividade moderada) | sem F10.7, sem tempestades, sem ciclo diurno | causa principal do erro do arrasto padrão (ver acima) |
| `step_s` na varredura de passagens | passagens mais curtas que o passo escapam | 20 s tem folga para a ISS a 10° (mínima de 132 s); rasantes exigem passo menor |
| Refração | não modelada | até ~0,5° perto do horizonte; desprezível acima de 10° |
| Máscara de horizonte | ε igual em todos os azimutes | relevo da Mantiqueira pode esconder passagens baixas |
| Coordenadas de São Bento | aproximadas (−22,69°, −45,73°, 900 m), **a confirmar** | trocar uma linha em `observer.py` |
| Literais da estação | ainda em `scripts/plot_ground_track.py`, `tests/test_scene_export.py` e `tests/test_sgp4_propagator.py` | S6-01 moverá para configuração |
| Quadros | UT1 ≈ UTC, sem movimento polar | ≈ 86 m contra o Skyfield; aceitável para visibilidade |
| Elevação no painel do viewer | interpolação linear de amostras de 20 s | até ~0,4° abaixo do máximo perto do TCA (a lista usa o valor exato) |
| Formato Alpha-5 do TLE | rejeitado | catálogos ≥ 100 000 não são lidos |
| `SceneSatellite`, `SceneStation` | provisórios | substituir por `Satellite` (S5-01) |
| `codigo_explicado.pdf` | não versionado (`*.pdf` no `.gitignore`) | regenerar do HTML com o comando no cabeçalho |
| Aviso do Skyfield com NumPy 2.5 | `DeprecationWarning` em `searchlib.py` | externo; acompanhar versões do Skyfield |

## 5. Alterações em código de semanas anteriores

| Arquivo | Alteração | Compatibilidade |
|---|---|---|
| `pyproject.toml` | dependência `sgp4>=2.23` | — |
| `propagators/base.py` | método não abstrato `propagate_utc(epoch, t)` | só acrescenta |
| `propagators/cowell.py` | `from_state`, `stop_altitude_m`, evento de reentrada, `reentry_time`, `NaN` após a reentrada | sem arrasto, comportamento idêntico (o evento só existe com `DragAcceleration`) |
| `propagators/__init__.py`, `forces/__init__.py` | novas exportações | só acrescenta |
| `analysis.py` | `cowell_j2_drag`, `position_error`, `fit_ballistic_coefficient` | só acrescenta |
| `app/scene_export.py` | `propagate_utc`; campos opcionais `footprint_deg`, `looks`, `passes`, `meta.el_min_deg`, `meta.tle`; época ao milissegundo; `SceneStation.from_ground_station` | `schema_version` = 1; cenas antigas idênticas |
| `app/static/viewer.*` | painel da estação, passagens, footprint, linha de visada, faixa do TLE | abre cenas antigas sem os painéis novos |
| `scripts/build_demo.py` | ISS real (TLE/SGP4) no lugar das ISS-like da S2-07; `--tle-file`, `--offline`, `--start`, `--cache-dir` | a demo mudou de propósito |
| `tests/test_golden_consistency.py` | teste de consistência interna do A.2 | só acrescenta |
| `docs/adr/0001-arquitetura.md` | TEME ≡ ECI; `propagate_utc`; contrato da cena da S3-06 | — |

## 6. Revisão do objetivo do projeto e recomendações de planejamento

Na revisão do S3-06 o objetivo final foi explicitado:

> O objetivo final é uma **simulação próxima da real para um satélite em SSO**, com um
> propagador **completo e generalizável para qualquer órbita**. A ISS é apenas o **comparador
> real**: serve para verificar, contra o SGP4, que o simulador está próximo de um simulador
> real.

Hoje a demo usa o SGP4 (biblioteca) para a ISS e o `MeanJ2Propagator` (só efeitos seculares
de J2) para a SSO; o Cowell com J2 + arrasto existe, mas só aparece no relatório do S3-03. O
roadmap atual (Semanas 4 a 7) **não cobre** esse objetivo: J3 aparece só como "órbita congelada
na Semana 5" (`secular.py`), a conversão médio → osculador como "volta na Semana 5" (relatório da
S2-06), e J4, arrasto no modelo médio e a generalização não estão previstos.

### Arquitetura recomendada: duas peças e uma ponte

1. **Cowell como propagador fiel e geral.** Já vale para qualquer órbita (integra as forças).
   Para ficar completo, novos modelos de força ao lado de `J2Acceleration` e `DragAcceleration`,
   sem mudar o propagador: **J3 e J4** (acelerações zonais), **Sol e Lua**, **pressão de
   radiação solar** e, depois, atmosfera com atividade solar (F10.7). Este é o "simulador
   próximo do real".
2. **Modelo médio analítico como ferramenta de projeto.** Evoluir o `MeanJ2Propagator` para
   J2 + J3 + J4 + arrasto médio (uma teoria do tipo Brouwer/SGP4). Rápido, ideal para projetar a
   órbita, achar a **órbita congelada** (J3 cancela ω̇ com ω = 90°) e planejar a manutenção da
   SSO, mas sempre uma aproximação.
3. **A ponte: conversão médio ↔ osculador (Brouwer/Kozai).** A peça que falta hoje. Os
   elementos de projeto da SSO são médios; entregues ao Cowell como osculadores, o semieixo médio
   fica 9,5 km fora e o plano atrasa ≈ 2,3°/ano em relação ao Sol (relatório da S2-06). Com a
   conversão, projeta-se no modelo médio e simula-se no Cowell sem perder a heliossincronia.

### Validação contra o SGP4 numa SSO real

A ISS não é SSO (i = 51,6°, 420 km). Satélites reais em SSO têm TLE público, por exemplo o
**Sentinel-2A** (NORAD 40697) e o **Landsat 8** (NORAD 39084). O código atual já suporta:
`get_tle(40697)` e o mesmo `compare_sgp4_cowell.py`. É o comparador real mais próximo do
satélite final.

### Ordem sugerida (para os próximos prompts)

| # | Entrega | Por quê |
|---|---|---|
| 1 | J3 e J4 no Cowell (`forces/zonal.py`) | completa a gravidade zonal; reduz o erro × SGP4 que hoje o C_D·A/m ajustado absorve |
| 2 | Conversão médio ↔ osculador (Brouwer/Kozai, 1ª ordem em J2) | liga projeto (médio) e simulação (Cowell); elimina a deriva de LTAN |
| 3 | Cowell × SGP4 numa SSO real (Sentinel-2A / Landsat 8) | valida o simulador no regime do satélite final |
| 4 | Modelo médio estendido (J2–J4 + arrasto médio) e órbita congelada | projeto e manutenção da SSO; comparação médio × Cowell |
| 5 | Sol/Lua, pressão de radiação e atmosfera com F10.7 no Cowell | fidelidade de longo prazo |
| 6 | Constelação SSO (6 → 18 satélites), cobertura, web | com o propagador fiel já validado |

A visibilidade óptica da ISS (Sol, iluminação; Semana 4 do roadmap original) continua útil e
independente; pode ser mantida, adiada ou feita em paralelo, a critério do planejamento.

## 7. Próximos passos

- `git push` da branch `semana-3` e abertura do PR ficam com o responsável do projeto.
- Confirmar as coordenadas exatas de São Bento do Sapucaí (ponto de observação e altitude).
- Levar a seção 6 ao planejamento das próximas semanas.
