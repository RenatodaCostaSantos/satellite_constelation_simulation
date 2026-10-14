# Relatório da Semana 2 — Kepler × MeanJ2 × Cowell

Gerado por `python scripts/compare_propagators.py` em 5 s. Época 2026-09-29 12:00 UTC; 10 dias a 60 s; Cowell DOP853 (rtol = 1e-11, atol = 1e-6), só J2. Ω0 = 30° e ψ0 = 0 em todas as órbitas; ω0 = 90° (SSO), 0° (ISS-like, ver nota) e 30° (e = 0,02). Faixas de referência: Apêndice A.4 do prompt da Semana 2.

## Números

| Órbita | a_médio − a_nominal [km] | Ω̇ ajustado [°/dia] | Ω̇ teoria@médio (erro) | Ω̇ teoria@nominal (erro) | erro (a) 24 h [km] | erro (b) 24 h [km] | \|Kepler − Cowell\| 10 d [km] | faixas A.4 |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| SSO de projeto | +9,50 | +0,9794 | +0,9803 (-0,085%) | +0,9857 (-0,633%) | 1345 | 6,6 | 13798 | OK |
| ISS-like | -5,97 | -4,9680 | -4,9641 (+0,080%) | -4,9466 (+0,432%) | 874 | 14,5 | 9749 | OK |
| órbita de teste e = 0,02 | -2,52 | -4,6251 | -4,6197 (+0,118%) | -4,6135 (+0,252%) | 342 | 6,6 | 7871 | OK |

**ω̇ na órbita e = 0,02:** ajustado +4,9095 °/dia vs teoria no semieixo médio +4,9018 °/dia (+0,157%; critério < 1%: OK).

Erro (a): MeanJ2 com os elementos nominais tratados como osculadores. Erro (b): MeanJ2 com a, e, i médios medidos no Cowell (ângulos nominais). Critérios: Ω̇ teoria@médio < 0,5%; (b) < 20 km; (a) > 500 km (SSO e ISS-like).

**Nota (ISS-like, ω0 = 0°):** o semieixo osculador oscila como 1,5·J2·R⊕²/a·sin²i·cos 2u (u = ω + ν), então o semieixo médio depende de onde a órbita começa. Com ω0 = 0° (início no nó) reproduz-se a linha ISS-like do Apêndice A.4; com ω0 = 90° o sinal de a_médio − a_nominal se inverte (ver o golden file, S2-02).

## Figuras

1. [Ω(t): Kepler, MeanJ2 e Cowell](fig1_raan.png)
2. [ω(t) na órbita e = 0,02](fig2_argp_e002.png)
3. [Semieixo osculador e médio](fig3_semieixo.png)
4. [\|MeanJ2 − Cowell\| em 24 h, (a) e (b)](fig4_erro_meanj2.png)
5. [\|Kepler − Cowell\| em 10 dias](fig5_erro_kepler.png)

## Conclusões

J2 faz o plano da órbita precessar: Ω̇ = +0,979 °/dia na SSO (a teoria nominal dá a taxa do Sol, 0,9856 °/dia; a diferença de -0,6% vem do semieixo médio, ver abaixo) e -4,968 °/dia na ISS-like; na órbita e = 0,02 ele também gira a linha dos apsides (ω̇ = +4,909 °/dia). Ignorar J2 (Kepler) leva a erros de 13798 km (SSO) e 9749 km (ISS-like) em 10 dias. O MeanJ2 reproduz a taxa medida no Cowell dentro de 0,12%, *desde que* a teoria use o semieixo médio medido: elementos médios e osculadores diferem (a_médio − a_nominal = +9,5 km na SSO, -6,0 km na ISS-like). Tratar os elementos nominais como osculadores dá 1345 km de erro em 24 h na SSO, contra 6,6 km com os médios medidos; o resíduo são as oscilações de curto período (Figura 3), que o MeanJ2 não representa por construção.
