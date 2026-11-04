# Relatório da Semana 5 — Constelação 91/6

Gerado por `python scripts/constellation_report.py` em 0,7 min. Época 2026-10-23 00:00 UTC. Constelação de referência do Apêndice A.1 do prompt da Semana 5 (`sso_constellation`: a 6 888,089 km, e congelada, i 97,46325°, ω 90°, LTAN 22h30, ψ0 = 0°, 300°, …, 60°). Cowell completo = J2–J4 + Sol/Lua + SRP (C_R 1,3, A/m 0,01 m²/kg) + arrasto (C_D·A/m 0,011 m²/kg, geodésico). Valores do Apêndice entre parênteses.

## S5-02 · Grade de 91 trilhas e repetição

![Nós ascendentes em 1 dia](grade_nos.png)

### Grade no equador

Nós por `node_crossings` (interpolação linear da latitude, passo 30 s); métricas por `grid_metrics` (espaçamentos com a volta de 360°; afastamento à grade ideal de 360°/91 = 3,95604°, ancorada no 1º nó de S1).

| Caso | Nós | Espaçamento mín. … máx. [°] | Afastamento máx. [°] |
|---|---:|---|---:|
| Modelo médio, P = 1, 6 dias | 91 | 3,95599 … 3,95627 (91; 3,9560 uniforme) | 0,0003  |
| Modelo médio, P = 6, 1 dia | 91 | 3,95600 … 3,95606 (91; 3,9560 uniforme) | 0,0000  |
| Modelo médio, P = 7, 1 dia | 107 | 0,56511 … 3,39090 (107; 0,5651 … 3,3909) | 1,6955  |
| Cowell completo, P = 6, 1 dia | 91 | 3,95173 … 3,96432 (91; 3,95173 … 3,96431) | 0,0064 (0,0064) |
| Cowell completo, P = 1, 6 dias (91 primeiros) | 91 | 3,91683 … 3,96738 (3,9169 … 3,9674) | 0,0414  |
| Cowell completo, P = 7, 1 dia | 106 | 3,38337 … 3,95848 (106; 3,3834 … 3,9585) | 1,6938 (1,694) |

- **Teorema P = D (Sem4 §8.3) no Cowell completo:** seis satélites em um dia depositam as 91 trilhas, com espaçamentos a ±0,006° do ideal; sete satélites na mesma órbita deixam a grade (afastamento de 1,7°, 6/7 de passo), com quase-duplicatas e lacunas mantidas.
- **Modelo médio:** com 1 satélite em 6 dias o passo varia 0,000285°. Não é a teoria: o ECEF gira com o GMST real, 8,554e-12 rad/s mais rápido que o ω⊕ = 7,292115e-5 rad/s do Sem4, o que desloca a grade em 0,000254° (≈ 28 m) em 6 dias.
- **Afastamento no Cowell (P = 6):** os desvios à grade têm padrão de período 3 (S1/S4, S2/S5, S3/S6, ±0,0045°): é o termo de curto período de J2 em Ω e u, que a conversão médio → osculador (S4-04) não corrige (ela ajusta a, o vetor e e i). Ancorado no 1º nó no tempo (o de S3), o afastamento máximo seria 0,0107°.
- **Precisão do nó:** contra nós refinados por `brentq` em z_ECI = 0 (S1, 15 nós em 24 h), a interpolação linear erra 7,3e-5° / 1,1e-3 s com passo de 60 s, 9,3e-6° / 2,7e-4 s com 30 s e 5,2e-7° com 10 s (`tests/test_repeat_grid.py`; A.3: 7,3e-5°, 9,3e-6°, 4,8e-7°). O refinamento por `brentq` não é necessário.

### Repetição de 6 dias e semieixo calibrado

`repeat_semi_major_axis_calibrated` (Cowell J2–J4, conversão só zonal): secante em a até N·Tnod = D·Tdia,nod, Tnod pelo ajuste linear dos nós e Tdia,nod = 2π/(ω⊕ − Ω̇) com Ω̇ por `fit_rate`. Resultado: **a = 6 888 108,37 m, Δa = +19,37 m** (A.4: 6 888 108,37 m, +19,37 m), em 3,6 s.

| Forças | a | Fechamento após 91 voltas | N·Tnod − D·Tdia,nod [s] | Tnod [s] | Ω̇ [°/d] |
|---|---|---|---:|---:|---:|
| J2–J4 | projeto | 0,00895° = 0,996 km (+0,00895° = 1,00 km) | −2,1973 | 5 696,6796 | 0,985632 |
| J2–J4 | calibrado | −0,00020° = −0,022 km (−0,00020° = −0,022 km) | 0,0000 | 5 696,7036 | 0,985622 |
| completo (com arrasto) | projeto | 0,04223° = 4,701 km (+0,0422° = 4,69 km) | −10,2042 | 5 696,5900 | 0,985531 |
| completo (com arrasto) | calibrado | 0,03306° = 3,680 km (+0,0330° = 3,68 km) | −8,0043 | 5 696,6140 | 0,985521 |

- **Fechamento com o a calibrado (J2–J4):** o resíduo de −22 m é, na maior parte, a mesma diferença ω⊕ × GMST real do modelo médio (−28 m previstos em 6 dias): a calibração usa o ω⊕ do Sem4 em Tdia,nod, como pede o card.
- **Com arrasto o ciclo não fecha** mesmo com o a calibrado (≈ 3,7 km por ciclo): é a deriva quadrática do Sem4 §7.6, assunto do S5-05. A calibração é feita sem arrasto, como a da inclinação.
- **Efeito no LTAN:** Δa = +19 m muda Ω̇ por Ω̇ ∝ a^(−7/2) (Sem4 eq. 32): −0,0142 min/ano (Sem4: −0,014). Medido no Cowell J2–J4 (diferença dos Ω̇ ajustados em 6 dias): −0,0142 min/ano. Desprezível diante da banda de ±5 min.

