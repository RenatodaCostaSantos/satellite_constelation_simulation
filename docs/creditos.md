# Créditos e licenças dos assets do viewer

Os arquivos abaixo estão versionados em `src/satsim/app/static/` para que o viewer funcione
offline (ADR 0001, requisito R10). Nenhum deles é baixado em tempo de execução.

## Bibliotecas JavaScript

### globe.gl 2.46.2 — `static/vendor/globe.gl.min.js`

- Fonte: <https://github.com/vasturiano/globe.gl>; arquivo obtido de
  `https://cdn.jsdelivr.net/npm/globe.gl@2.46.2/dist/globe.gl.min.js` em 02/10/2026.
- SHA-256: `2c3e445c04d121215910a89688b96091c8a72071c122a4f830081a39b636c94c`
- Licença: MIT, © 2019 Vasco Asturiano. Texto em `static/vendor/LICENSE-globe.gl.txt`.
- O bundle UMD inclui as dependências abaixo, com suas licenças:

| Componente embutido | Licença | Observação |
|---|---|---|
| three.js r185 | MIT, © 2010–2026 three.js authors | texto em `static/vendor/LICENSE-three.txt` |
| three-globe, three-render-objects, kapsule, accessor-fn | MIT, © Vasco Asturiano | mesmo autor do globe.gl |
| @tweenjs/tween.js | MIT | animações de câmera |
| h3-js | Apache-2.0, © 2018–2022 Uber Technologies, Inc. | camada de hexágonos do three-globe (não usada pelo viewer) |
| Fonte Helvetiker (typeface) | licença MgOpen, © 2004 MAGENTA Ltd. (permissiva) | rótulos 3D do three-globe (não usada pelo viewer) |

## Texturas

### Terra — `static/textures/earth_2048.jpg`

- "Blue Marble: Land Surface, Shallow Water, and Shaded Topography", NASA Visible Earth,
  imagem 57752, arquivo `land_shallow_topo_2048.jpg` (2048 × 1024, projeção equiretangular).
- Fonte: <https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57752/land_shallow_topo_2048.jpg>,
  obtido em 02/10/2026.
- SHA-256: `5b54cc586c6cbf2b28762ef4d4011f6cf4227a8b93a637b818a0c54090ce6c2c`
- Crédito: NASA Goddard Space Flight Center / NASA Earth Observatory (projeto Blue Marble).
- Licença: domínio público. Imagens da NASA não são protegidas por direitos autorais nos EUA
  (diretrizes de uso de mídia da NASA); pede-se o crédito acima.

### Céu estrelado

Gerado em tempo de execução por `viewer.js` (canvas com gerador pseudoaleatório de semente fixa).
Não há arquivo de imagem nem licença de terceiros.

## Verificação de uso offline

`grep -n "http"` nos assets do viewer (S1-08):

- `viewer.html`, `viewer.js` e `viewer.css`: nenhuma ocorrência.
- `globe.gl.min.js`: só namespaces XML/SVG (`http://www.w3.org/...`), links em comentários, em
  mensagens de erro e em avisos de licença (three.js, h3-js, Helvetiker). Nenhum é um recurso
  carregado em tempo de execução.
- Teste no navegador (Chrome, `file://`, com todas as requisições que não fossem `file:`,
  `data:` ou `blob:` registradas e bloqueadas): **zero requisições externas**.
