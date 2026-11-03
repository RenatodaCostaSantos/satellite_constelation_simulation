# Fechamento da Semana 4 — Propagador de alta fidelidade e SSO calibrada

Branch `semana-4`. Seis cards concluídos e aprovados (S4-01 a S4-06). Este documento segue a seção
7 do prompt da Semana 4.

## 1. O que foi entregue

| Card | Entrega | Commits |
|---|---|---|
| S4-01 · Gravidade zonal | `forces/zonal.py`: potencial e aceleração J2–J4 por polinômios de Legendre (recorrência de Bonnet), `ZonalAcceleration` para o Cowell | `0b44a89`, `140693a` |
| S4-02 · Sol e Lua | `ephemeris.py`: Sol (Astronomical Almanac) e Lua (Montenbruck & Gill §3.3.2) no ECI da data, α☉/δ☉, `ltan_hours`, `beta_angle`; `forces/third_body.py` (maré de terceiro corpo); constantes `MU_SUN`, `MU_MOON`, `AU`; seção Limitações do README | `87118fd`, `9169443`, `97ab98a`, `74cc642` |
| S4-03 · SRP e eclipse | `forces/srp.py` (bola de canhão com sombra cilíndrica), `in_shadow_cylindrical`, `eclipse_fraction_circular` (Sem4 eq. 30), `P_SUN` | `16d7e79`, `e27c236` |
| S4-04 · Médio ↔ osculador | `mean_osc.py`: médias por janela nodal (vetor excentricidade), `osculating_to_mean`, `mean_to_osculating` (método do tiro, 4 iterações na órbita de projeto) | `d51b662`, `7e2b9e4` |
| S4-05 · Modelo médio e SSO calibrada | `secular.zonal_mean_rates` (J2–J4 de 1ª ordem, A.5), `MeanZonalPropagator` (com arrasto médio), `frozen_eccentricity`, `sso_inclination_calibrated` (secante no Cowell), `raan_for_ltan`; convenção "conversão só com gravidade zonal" (`conversion_models`) | `15aa992`, `c494149` |
| S4-06 · Validação e GUI | `orbit_determination.fit_initial_state` (mínimos quadrados sobre estado e C_D·A/m), `scripts/compare_high_fidelity.py` e `reports/s4_alta_fidelidade/`; TLE do Sentinel-2A congelado e embarcado; demo com a SSO-1 no Cowell completo, o Sentinel-2A real, e LTAN, β e eclipse no painel | `e4a74df`, `e5ef0aa`, `b684cde` |

Os commits "ajustes de revisão" incluem as seções 13 a 16 do `docs/codigo_explicado` (terceiro
corpo e a analogia do peão, eclipse e a corda × arco, elementos médios e o método do tiro, SSO-1 ×
Sentinel-2A na demo), comentários linha a linha pedidos na revisão e duas subseções da seção
Limitações do README: manutenção da SSO e a referência da validação ser o SGP4.

**Validação principal (contra o Apêndice A do prompt):**

| O quê | Resultado | Referência / critério |
|---|---|---|
| Zonal: parte J2 × `j2_acceleration` / gradiente do potencial / energia em 10 d | 1e-15 / 1e-9 / 3e-11 relativos | ≤ 1e-12 / 1e-7 / 1e-9 |
| Sol × DE421 em 2026 (direção) | máx. 0,014° | ≤ 0,02° |
| Lua × DE421 em 2026 (direção / distância) | máx. 0,085° / 439 km | ≤ 0,1° / 500 km |
| SRP em 2026-10-23 | 5,986762e-8 m/s² | rel 1e-4 |
| Fração de eclipse da SSO de demonstração (24 h) | 0,3603 (eq. 30: 0,3650) | A.4: 0,3603 |
| Médio → osculador, órbita de projeto | converge em 4 iterações; voa congelada (ω 90,00…90,10° em 10 d) | ≤ 6; ω ∈ [89,5; 90,5]° |
| Taxas J2–J4 de 1ª ordem | Ω̇ 0,983369 °/d, ω̇ +0,01039 °/d | A.5 |
| Excentricidade congelada J2–J4 | 1,0772101e-3 | A.5: 1,07721e-3 |
| Inclinação calibrada (SSO de projeto / Sentinel-2A) | 97,46330° / 98,56940° (real 98,56361°) | A.7: 97,46325° / a < 0,01° do real |
| LTAN da SSO calibrada no Cowell J2–J4 | −6e-5 min/ano | < 0,1 min/ano |
| Modelo médio × Cowell (30 d) | Ω̇ +0,085%, Tnod +2,3 ppm | ≤ 0,1% / 10 ppm |
| Cowell × SGP4, ajuste de 72 h, forças completas | Sentinel-2A 0,01 / 0,02 / 0,36 km; ISS 0,03 / 0,06 / 0,49 km | ≤ 0,2 / 0,3 / 3 e ≤ 0,5 / 1 / 8 km |
| SSO-1 da demo (Cowell completo, 24 h) | LTAN 22,4997…22,5026 h; eclipse 0,3603 | 22,500 ± 0,01 h; 0,36 ± 0,02 |

## 2. Testes e lint

`ruff check .`:

```text
All checks passed!
```

`ruff format --check .`: `75 files already formatted`.

`pytest --cov=satsim`: **463 passed, 1 xfailed** em 161,6 s com a medição de cobertura (84,8 s sem
ela). O xfail é a discrepância conhecida do Sem4 §1, registrada na Semana 1. Os 2 avisos são um
`DeprecationWarning` interno do Skyfield com o NumPy 2.5 (`searchlib.py`), não do nosso código.

```text
Name                                   Stmts   Miss  Cover
----------------------------------------------------------
src/satsim/__init__.py                     1      0   100%
src/satsim/analysis.py                    43      0   100%
src/satsim/app/__init__.py                 0      0   100%
src/satsim/app/scene_export.py           121      0   100%
src/satsim/astrotime.py                   32      0   100%
src/satsim/constants.py                   40      0   100%
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
src/satsim/groundtrack.py                 27      0   100%
src/satsim/mean_osc.py                   101      0   100%
src/satsim/observer.py                    65      0   100%
src/satsim/orbit_design.py                81      4    95%
src/satsim/orbit_determination.py         65      4    94%
src/satsim/passes.py                      95      3    97%
src/satsim/propagators/__init__.py         7      0   100%
src/satsim/propagators/base.py            20      0   100%
src/satsim/propagators/cowell.py         100      2    98%
src/satsim/propagators/kepler.py          34      2    94%
src/satsim/propagators/mean_j2.py         34      0   100%
src/satsim/propagators/mean_zonal.py      78      4    95%
src/satsim/propagators/sgp4_ref.py        47      0   100%
src/satsim/secular.py                     76      2    97%
src/satsim/tle.py                        195      2    99%
----------------------------------------------------------
TOTAL                                   1769     30    98%
============ 463 passed, 1 xfailed, 2 warnings in 161.62s (0:02:41) ============
```

## 3. Discrepâncias encontradas (prompt × código)

Todas estão registradas no golden (`tests/golden/reference_values.py`) com os valores medidos.

1. **A.3, erro da Lua "|Δr| ≤ 445 km" (S4-02).** O número do A.3 é a diferença de **distância**
   (|r| − |r_DE421|). A diferença vetorial inclui o erro de direção: em 2026-10-23, 439 km de
   distância e 494 km de vetor. Os testes conferem os dois.
2. **Sem4 §10, precessão lunissolar "~1e-3 °/dia em LEO" (S4-02).** O A.3 mede −2,45e-5 °/d (≈ 40×
   menor), e o satsim confirma a ordem de grandeza. O valor aparente depende da janela: com os
   elementos médios usados como osculadores, o ajuste linear de Ω em 30 dias dá −3,52e-4 °/d,
   porque o Sol produz di/dt ≈ −9,8e-5 °/d na SSO e, via J2, ΔΩ ∝ t². Isso motivou a seção
   Limitações do README (manutenção da SSO).
3. **A.6, janela das médias (S4-04).** O A.6 reproduz-se com janelas do período **kepleriano**
   (5 689,30 s), não do nodal (5 696,70 s) que o card pede. Com a janela nodal (implementada):
   osculador a = 6 878 679,4 m (A.6: 6 878 668,405 m) e resíduos 9 483,6 / 74,30 / 0,2460 m (A.6:
   9 495 / 74,46 / 0,2468 m). A janela kepleriana fica 7,4 s curta e enviesa ⟨a⟩ em ≈ 12 m. É o
   mesmo viés que explica o "a − alvo −0,0 … +16,6 m" do A.6: com a janela nodal, o semieixo médio
   fica constante (−0,0 m em 30 dias).
4. **A.7, período nodal medido (S4-05).** O satsim mede Tnod = 5 696,6796 s no Cowell (+2,3 ppm da
   teoria), contra 5 696,6659 s (+4,7 ppm) no A.7. É a consequência do item 3: com a janela
   kepleriana, o estado convertido fica com o semieixo médio real ≈ 11 m baixo, o que encurta o
   período em 2,4 ppm. Equivalente ao longo da órbita em 30 dias: ≈ 45 km (A.7: ≈ 92 km).
   Inclinação calibrada: 97,46330° (A.7: 97,46325°, dentro da tolerância de 0,0005°).
5. **A.5, precisão impressa de ω̇ (S4-05).** O card pede rel 1e-6, mas o A.5 imprime ω̇ congelado
   com 4 algarismos (+0,01039 °/d). O teste compara ω̇ com tolerância de arredondamento (5e-6 °/d).
   O valor medido é +0,0103926 °/d.
6. **A.8, ajuste (B) (S4-06).** O (A) reproduz o A.8 nos 8 casos. Em (B) o satsim chega a erros e
   RMS **menores** em todos os casos, por exemplo:

   | Caso | RMS satsim | RMS A.8 |
   |---|---|---|
   | Sentinel-2A J2–J4 + arrasto | 0,007 km | 0,043 km |
   | ISS J2–J4 + arrasto | 0,008 km | 0,109 km |

   O RMS é o que o ajuste minimiza, então o otimizador do experimento provavelmente parou antes do
   mínimo. Causa provável: o passo das diferenças finitas. Com passo relativo 1e-6 o ajuste
   sintético do satsim também parava a 0,33 m do estado certo; com 1e-4 (o adotado) chega a 1,3 mm.
7. **A.8, C_D·A/m do Sentinel-2A (S4-06).** O ajuste (B) dá 0,0138 m²/kg, contra 0,00208 do A.8
   (na ISS: 0,00119 contra 0,00134). Interpretação: a 786 km o coeficiente efetivo reproduz o
   decaimento que o modelo de densidade do *SGP4* aplica, que não tem por que concordar com a
   tabela de Vallado. Não foi verificado contra o B\* do TLE.

## 4. Dívida técnica e pontos frágeis

| Ponto | Situação | Impacto / o que fazer |
|---|---|---|
| J2² no modelo médio | `zonal_mean_rates` é de 1ª ordem: Ω̇ fica +0,085% acima do Cowell na SSO calibrada e o período nodal +2,3 ppm | ≈ 45 km ao longo da órbita em 30 dias: o `MeanZonalPropagator` é ferramenta de projeto, não de posição. A calibração numérica contorna isso; Brouwer de 2ª ordem fica fora do escopo |
| Singularidade em e → 0 | os termos de J3 em ω̇ e ψ̇ vão com 1/e; `MeanZonalPropagator` recusa e = 0 com J3 | sem efeito em órbitas congeladas; uma órbita que passe perto de e = 0 exigiria variáveis não singulares (e·cos ω, e·sin ω, u) |
| Sombra cilíndrica | sem penumbra: entrada e saída instantâneas (≈ 8 s de transição real em LEO); degrau de ≈ 6e-8 m/s² no Cowell | erro de posição ~1e-4 m por degrau; penumbra cônica ou eventos de sombra no integrador quando a precisão pedir |
| Atmosfera estática | tabela exponencial de Vallado (atividade moderada), sem F10.7/Ap | o C_D·A/m ajustado é efetivo (absorve a atmosfera do dia e, contra o SGP4, o modelo de densidade do próprio SGP4) |
| Referência é o SGP4 | "mais correto" = mais perto do SGP4 (≈ 1 km de erro próprio); Sol/Lua/SRP pioram o acordo | trocar a referência não exige mudar `orbit_determination.py`: basta um `Propagator` sobre efemérides precisas (POD/SP3) (README, Limitações) |
| Manutenção da SSO | o Sol produz di/dt ≈ −0,036°/ano na SSO de projeto; nada no simulador modela manobras | indicador de manutenção no painel proposto (seção 7) |
| Tempo dos testes | suíte em 84,8 s (161,6 s com cobertura), perto do limite de ~2 min. Mais lentos: fixture de taxas do `test_zonal` (15,0 s), ajuste do coeficiente balístico da Semana 3 (11,9 s), Cowell de 30 dias do `test_mean_zonal` (11,3 s), calibração do Sentinel-2A (5,2 s) e ajuste do Sentinel-2A em 24 h (5,1 s) | a Semana 5, com 6 satélites, deve marcar testes longos (`@pytest.mark.slow`) ou reaproveitar fixtures de módulo; a calibração e a conversão custam ≈ 5 s e ≈ 1,2 s por órbita |
| Script de alta fidelidade | 11,9 min; `resultados.json` (325 kB) versionado para refazer o relatório com `--from-json` | rodar de novo só quando mudar o modelo de forças |
| TLEs embarcados | ISS (13/10/2026) e Sentinel-2A (23/10/2026) envelhecem; a demo omite o Sentinel-2A se o TLE estiver a mais de 30 dias do início da cena | atualizar com `scripts/fetch_tle.py --norad <id>` (docs/creditos.md) |
| Calibração usa Ω osculador em 10 dias | com Sol/Lua, o Ω̇ aparente depende da janela (item 2 da seção 3) | o teste calibra só com J2–J4 (opção recomendada no S4-05); calibrar com Sol/Lua dá +0,0002° (A.3) |

## 5. Alterações em código de semanas anteriores

Todas compatíveis (nenhum teste antigo removido ou afrouxado):

- `constants.py`: acréscimo de `MU_SUN`, `MU_MOON`, `AU`, `P_SUN` (S4-02, S4-03).
- `forces/__init__.py` e `propagators/__init__.py`: novas exportações.
- `secular.py` (S4-05): acréscimo de `ZonalMeanRates` e `zonal_mean_rates`; a docstring do módulo
  deixou de dizer que J3 fica "para a Semana 5". `secular_rates` não mudou.
- `orbit_design.py` (S4-05): acréscimo de `frozen_eccentricity`, `raan_for_ltan` e
  `sso_inclination_calibrated` (comentada linha a linha na revisão).
- `propagators/cowell.py`: só a docstring, que dizia que a conversão médio → osculador não existia.
- `tle.py` (S4-06): `BUNDLED_TLES` (catálogo → arquivo), `SENTINEL_2A_NORAD_ID` e
  `load_bundled_tle(norad_id=ISS)`. O padrão continua sendo a ISS.
- `app/scene_export.py` (S4-06): parâmetro `sun_geometry=False`, que acrescenta `ltan_h`,
  `beta_deg` e `eclipse`; `_round_ms` passou a ser `round_ms` (público), usado pelo build_demo.
- `app/static/viewer.js`, `viewer.html`, `viewer.css` (S4-06): linhas LTAN, β e Sol no painel
  Informações e o marcador esmaecido em eclipse; cenas antigas mostram "—".
- `scripts/build_demo.py` (S4-06): a lista `SATELLITES` com o helper `sat()` deu lugar a
  `design_satellites(start)`. A SSO-1 passou de `MeanJ2Propagator` para o Cowell completo, e o
  Sentinel-2A real entrou como terceiro satélite.
- `tests/test_tle.py` e `tests/test_scene_export.py`: testes novos; o teste offline do build_demo
  passou a usar `--start` fixo (data do A.7), para não depender da idade do TLE embarcado.
- `README.md` (Demo visual, Limitações), `docs/adr/0001-arquitetura.md` (contrato da cena),
  `docs/creditos.md` (TLE do Sentinel-2A).

## 6. O que a Semana 5 recebe pronto

**Receita "projeto médio → calibração → conversão → Cowell completo"**, a mesma da SSO-1 da demo
(`scripts/build_demo.py`, `sso_design`), para instanciar os 6 satélites da constelação 91/6:

```python
import math
from satsim.astrotime import datetime_to_jd
from satsim.elements import KeplerianElements
from satsim.forces import (DragAcceleration, SolarRadiationPressure, ThirdBodyAcceleration,
                           ZonalAcceleration)
from satsim.mean_osc import conversion_models, mean_to_osculating
from satsim.orbit_design import frozen_eccentricity, raan_for_ltan, sso_inclination_calibrated
from satsim.propagators import CowellPropagator

a = 6888.089e3                                   # 1. projeto médio (91 rev em 6 dias, Sem4 §7.5)
i0 = math.radians(97.4396)                       #    1ª ordem só J2 (eq. 27), só para o e congelado
e = frozen_eccentricity(a, i0)                   #    J2–J4: 1,07721e-3, ω = 90°
models = (ZonalAcceleration(), ThirdBodyAcceleration("sun", epoch),
          ThirdBodyAcceleration("moon", epoch),
          SolarRadiationPressure(epoch, cr=1.3, area_over_mass=0.01),
          DragAcceleration(cd=2.2, area_over_mass=0.005))
i = sso_inclination_calibrated(a, e, epoch, conversion_models(models))   # 2. calibração J2–J4 (≈ 5 s)
raan = raan_for_ltan(22.5, datetime_to_jd(epoch))                        #    LTAN 22h30

satellites = []
for psi in phases:                               # fases (e/ou planos) do desenho da constelação
    mean = KeplerianElements(a, e, i, raan, math.pi / 2, psi)
    osc = mean_to_osculating(mean, epoch, conversion_models(models)).elements   # 3. conversão
    satellites.append(CowellPropagator(osc, epoch, models))                    # 4. Cowell completo
```

- A inclinação calibrada é a mesma para os 6 satélites (mesmos a e e), então calibra-se uma vez.
- A conversão usa só a gravidade zonal (convenção do S4-05); arrasto, SRP e Sol/Lua entram só na
  propagação.
- Para projeto rápido (taxas, períodos, derivas de meses), há o `MeanZonalPropagator` (J2–J4 de 1ª
  ordem + arrasto médio), a ≤ 0,1% do Cowell em Ω̇.
- Para validar contra satélites reais: `fit_initial_state` (SGP4 como referência) e o relatório
  `reports/s4_alta_fidelidade/`.
- A cena já exporta LTAN, β e eclipse por satélite, e o viewer os mostra.

## 7. Próximos passos

- **Indicador de manutenção da SSO no painel** (proposto na revisão do S4-02): deriva prevista de
  LTAN e de inclinação, data da próxima correção dentro de uma faixa de tolerância e Δv/ano
  (≈ 4,8 m/s por ano na SSO de projeto). Exige um propagador médio com o termo solar em di/dt,
  para escalas de meses a anos.
- Opcional, se útil para a apresentação: "Sentinel-2A (Cowell)" ajustado ao SGP4 na demo, com a
  distância entre os dois em metros no painel (discutido na revisão do S4-06).
- Lembrete: o `git push` da branch `semana-4` e a abertura do PR ficam com você.
