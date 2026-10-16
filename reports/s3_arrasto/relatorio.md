# Relatório do S3-03 — Cowell (J2 + arrasto) × SGP4

Gerado por `python scripts/compare_sgp4_cowell.py` em 19 s (ajuste: 11 s). TLE: ISS (ZARYA) (catálogo 25544), época 2019-12-09 16:38:29 UTC. Todas as configurações partem do estado do SGP4 em t = 0; Cowell DOP853 (rtol = 1e-10, atol = 1e-3); atmosfera exponencial por partes de Vallado (`satsim.forces.atmosphere`), co-rotante. Referência: Apêndice A.5 do prompt da Semana 3.

## Erro de posição 3D contra o SGP4

| Configuração | C_D·A/m [m²/kg] | 24 h [km] | 72 h [km] | 7 dias [km] | A.5 (24 h / 72 h / 7 d) |
|---|---:|---:|---:|---:|---|
| (a) J2 apenas | 0,00000 | 1,0 | 8,7 | 47,2 | 1 / 8,7 / 47 |
| (b) J2 + arrasto padrão, alt. esférica | 0,00880 | 6,6 | 59,6 | 325,2 | 6,6 / 59,6 / 325 |
| (b) J2 + arrasto padrão, alt. geodésica | 0,00880 | 5,8 | 52,4 | 286,1 | 5,8 / 52,4 / 286 |
| (c) J2 + arrasto ajustado (C_D·A/m = 0,00163) | 0,00163 | 0,3 | 2,7 | 14,5 | 0,3 / 2,7 / 14,5 |

**C_D·A/m ajustado:** 0,00163 m²/kg (janela de 72 h, amostras a cada 600 s, mínimos quadrados no RMS do erro de posição). Com C_D = 2,2 equivale a A/m = 0,00074 m²/kg. Referência A.5: 0,00163; faixa aceita: 0,0005 a 0,005.

**Critério do card (ajustado ≤ 2 / 8 / 40 km):** 0,3 / 2,7 / 14,5 km → OK.

![Erro de posição Cowell × SGP4](erro_sgp4_cowell.png)

## Por que o arrasto padrão piora o resultado

Com os parâmetros padrão (C_D = 2,2, A/m = 0,004) o Cowell fica **mais longe** do SGP4 do que só com J2: 5,8 / 52,4 / 286,1 km contra 1,0 / 8,7 / 47,2 km em 24 h / 72 h / 7 dias. Não é erro de código: o arrasto está implementado e validado contra o Sem3 (1,2e-6 m/s² a 400 km, −11,8 m por revolução, Tabela 3 e tempos de vida). O problema é o *valor* do arrasto:

- **Mínimo solar.** Dezembro de 2019 foi o fundo do mínimo do ciclo solar 24/25. Com pouca radiação EUV a termosfera esfria e se contrai, e a densidade a 400 km fica várias vezes abaixo da média. A tabela de Vallado representa atividade solar *moderada* (Sem3 Tabela 3: ρ varia por fatores de 5 a 10 ao longo do ciclo), então superestima o arrasto real daquela semana.
- **A/m efetivo da ISS.** A área frontal depende da atitude e da orientação dos painéis solares, que giram para seguir o Sol; 0,004 m²/kg é um valor nominal.
- **O que o ajuste absorve.** O C_D·A/m ajustado (0,00163, cerca de 19% do padrão) é um parâmetro *efetivo*: compensa a atmosfera do dia, a atitude e também parte do que o Cowell não modela e o SGP4 sim (J3, J4 e, em menor grau, Sol e Lua). Por isso o erro ajustado não vai a zero e cresce com o tempo.

A diferença entre as linhas esférica e geodésica do arrasto padrão vem da altitude usada na densidade: medida sobre o elipsoide WGS84, a altitude da ISS é maior que pela esfera de raio R⊕ (de 0 no equador a ≈ 13 km nas latitudes máximas, 51,6°), e a densidade ali é menor, então o arrasto geodésico é um pouco mais fraco.

**Lição:** um modelo de arrasto correto com parâmetros genéricos pode ser pior do que nenhum arrasto. Na prática, operadores ajustam um coeficiente balístico efetivo contra dados de rastreio (Sem3 §7, item 4), exatamente o que `fit_ballistic_coefficient` faz aqui com o SGP4 no lugar dos radares. A observação sobre atividade solar volta no documento da Semana 3.

## Fora do escopo

- Atmosferas empíricas (NRLMSISE-00, JB2008) e variação com o índice F10.7.
- Tempestades geomagnéticas.
- J3/J4 e atração lunissolar no Cowell.
