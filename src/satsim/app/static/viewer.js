/*
 * Viewer do satsim (ADR 0001, seção GUI).
 *
 * Regras: a física fica em Python. Este arquivo só lê a cena (window.SATSIM_SCENE), interpola
 * no tempo entre amostras e projeta para desenhar. Nenhuma equação orbital aqui.
 * Uma única variável, simTime (segundos desde meta.epoch_utc), comanda o globo 3D e o mapa 2D.
 * Qualquer número de satélites é lido do JSON (R7).
 *
 * S3-06 (campos opcionais da cena): estação com az/el/distância por satélite (``looks``),
 * passagens, círculo de visibilidade (``footprint_deg``) e o TLE da ISS real (``meta.tle``).
 * O círculo do footprint é desenhado com geometria de esfera (destino a uma distância angular
 * dada), não com equações orbitais: o raio vem pronto do Python.
 */
(function () {
  "use strict";

  // --- parâmetros de exibição (não são física) ------------------------------------------------
  var SPEEDS = [1, 10, 60, 300, 1000];
  var DEFAULT_SPEED = 60;
  var TRAIL_3D_S = 15 * 60; // trilha recente no globo
  var AHEAD_3D_S = 100 * 60; // órbita prevista no globo (~1 volta em LEO)
  var PAST_2D_S = 100 * 60; // traço passado no mapa
  var AHEAD_2D_S = 100 * 60; // traço futuro no mapa
  var EARTH_RADIUS_KM = 6378.137; // só para converter alt_km em raios do globo (escala de desenho)
  var STRIDE_3D_AHEAD = 3; // órbita prevista no globo: 1 amostra a cada 3 passos
  var STRIDE_2D = 2; // traços no mapa: 1 amostra a cada 2 passos
  var TEXTURE = window.SATSIM_TEXTURE || "textures/earth_2048.jpg";
  var D2R = Math.PI / 180;
  var BRT_OFFSET_MS = -3 * 3600 * 1000; // horário de Brasília: UTC−03:00 fixo (sem horário de verão)
  var TLE_STALE_DAYS = 3; // idade a partir da qual a faixa do TLE fica em destaque
  var FOOTPRINT_POINTS = 96; // vértices do círculo de visibilidade no mapa
  var LOS_POINTS = 24; // vértices da linha estação → satélite no globo
  var VIS_COLORS = { up: "#3ddc84", low: "#ffd166", down: "#8a97ad" };

  var scene = window.SATSIM_SCENE;
  if (!scene || typeof Globe === "undefined") {
    showError(
      "Cena ou biblioteca do globo não encontrada. Gere a demo com: python scripts/build_demo.py"
    );
    return;
  }
  if (scene.meta.schema_version !== 1) {
    showError("Versão de cena não suportada: " + scene.meta.schema_version);
    return;
  }

  var STEP = scene.meta.step_s;
  var DURATION = scene.meta.duration_s;
  var EPOCH_MS = Date.parse(scene.meta.epoch_utc);
  var sats = scene.satellites.map(function (s) {
    return {
      id: s.id,
      name: s.name,
      color: s.color,
      lla: s.lla,
      footprint: s.footprint_deg || null, // S3-06 (opcional)
      visible: true,
      now: null,
    };
  });
  var stations = scene.stations;
  // Cenas da S2-07 em diante trazem meta.geodesy = "WGS84"; as anteriores, subponto esférico.
  var GEODETIC = scene.meta.geodesy === "WGS84";

  // --- S3-06: estação, passagens e TLE (todos opcionais; cenas antigas não os têm) ------------
  var EL_MIN = typeof scene.meta.el_min_deg === "number" ? scene.meta.el_min_deg : null;
  var STATION = stations.length ? stations[0] : null;
  var TLE = scene.meta.tle || null;
  var PASSES = (scene.passes || []).map(function (p) {
    return {
      sat_id: p.sat_id,
      station_id: p.station_id,
      aos: (Date.parse(p.aos_utc) - EPOCH_MS) / 1000, // segundos desde o início da cena
      tca: (Date.parse(p.tca_utc) - EPOCH_MS) / 1000,
      los: (Date.parse(p.los_utc) - EPOCH_MS) / 1000,
      aosMs: Date.parse(p.aos_utc),
      maxEl: p.max_el_deg,
      duration: p.duration_s,
    };
  });
  sats.forEach(function (s) {
    s.looks = STATION && STATION.looks && STATION.looks[s.id] ? STATION.looks[s.id] : null;
    s.passes = PASSES.filter(function (p) {
      return p.sat_id === s.id && STATION && p.station_id === STATION.id;
    });
  });
  var HAS_STATION_VIEW = EL_MIN !== null && sats.some(function (s) {
    return s.looks;
  });

  var state = {
    simTime: 0,
    playing: true,
    speed: DEFAULT_SPEED,
    dragging: false,
    selected: sats.length ? sats[0] : null,
  };

  // --- interpolação temporal (trata o salto de ±180° na longitude) -----------------------------
  function wrapLon(lon) {
    if (lon >= 180) return lon - 360;
    if (lon < -180) return lon + 360;
    return lon;
  }

  function sampleAt(lla, t) {
    var n = lla.length;
    var f = t / STEP;
    if (f <= 0) return lla[0].slice();
    if (f >= n - 1) return lla[n - 1].slice();
    var k = Math.floor(f);
    var u = f - k;
    var a = lla[k];
    var b = lla[k + 1];
    var dlon = b[1] - a[1];
    if (dlon > 180) dlon -= 360;
    else if (dlon < -180) dlon += 360;
    return [a[0] + u * (b[0] - a[0]), wrapLon(a[1] + u * dlon), a[2] + u * (b[2] - a[2])];
  }

  // Interpolação linear de um array amostrado a cada STEP (footprint, el, distância).
  function sampleArr(arr, t) {
    var n = arr.length;
    var f = t / STEP;
    if (f <= 0) return arr[0];
    if (f >= n - 1) return arr[n - 1];
    var k = Math.floor(f);
    return arr[k] + (f - k) * (arr[k + 1] - arr[k]);
  }

  // Idem para azimute em [0°, 360°): interpola pelo caminho curto (359° → 1° passa por 0°).
  function sampleAz(arr, t) {
    var n = arr.length;
    var f = t / STEP;
    if (f <= 0) return arr[0];
    if (f >= n - 1) return arr[n - 1];
    var k = Math.floor(f);
    var d = arr[k + 1] - arr[k];
    if (d > 180) d -= 360;
    else if (d < -180) d += 360;
    var az = arr[k] + (f - k) * d;
    return az < 0 ? az + 360 : az >= 360 ? az - 360 : az;
  }

  // O satélite visto da estação no instante t, e o estado de visibilidade.
  function lookAt(s, t) {
    if (!s.looks) return null;
    var el = sampleArr(s.looks.el_deg, t);
    return {
      az: sampleAz(s.looks.az_deg, t),
      el: el,
      range: sampleArr(s.looks.range_km, t),
      state: el >= EL_MIN ? "up" : el >= 0 ? "low" : "down",
    };
  }

  // Pontos entre t0 e t1: extremos interpolados + amostras intermediárias (uma a cada
  // ``stride`` amostras, para aliviar o desenho de caminhos longos).
  function track(lla, t0, t1, stride) {
    stride = stride || 1;
    t0 = Math.max(0, t0);
    t1 = Math.min(DURATION, t1);
    if (t1 <= t0) return [sampleAt(lla, t0)];
    var pts = [sampleAt(lla, t0)];
    var k0 = Math.floor(t0 / STEP) + 1;
    var k1 = Math.ceil(t1 / STEP) - 1;
    for (var k = k0; k <= k1; k++) if (k % stride === 0) pts.push(lla[k]);
    pts.push(sampleAt(lla, t1));
    return pts;
  }

  // --- utilidades -------------------------------------------------------------------------------
  function showError(msg) {
    var el = document.getElementById("error");
    if (el) {
      el.textContent = msg;
      el.hidden = false;
    }
  }

  function rgba(hex, alpha) {
    var h = hex.replace("#", "");
    if (h.length === 3) h = h[0] + h[0] + h[1] + h[1] + h[2] + h[2];
    var n = parseInt(h, 16);
    return "rgba(" + ((n >> 16) & 255) + "," + ((n >> 8) & 255) + "," + (n & 255) + "," + alpha + ")";
  }

  function fmtUtc(ms) {
    return new Date(ms).toISOString().replace("T", " ").slice(0, 19) + " UTC";
  }

  function fmtNum(v, digits) {
    return v.toFixed(digits).replace(".", ",");
  }

  // "dd/mm HH:MM:SS" de um instante (ms) deslocado por offsetMs (0 = UTC, BRT_OFFSET_MS = BRT).
  function fmtDayTime(ms, offsetMs) {
    var iso = new Date(ms + offsetMs).toISOString();
    return iso.slice(8, 10) + "/" + iso.slice(5, 7) + " " + iso.slice(11, 19);
  }

  function fmtDuration(seconds) {
    var s = Math.max(0, Math.round(seconds));
    var h = Math.floor(s / 3600);
    var m = Math.floor((s % 3600) / 60);
    var sec = s % 60;
    if (h > 0) return h + " h " + (m < 10 ? "0" : "") + m + " min";
    if (m > 0) return m + " min " + (sec < 10 ? "0" : "") + sec + " s";
    return sec + " s";
  }

  function fmtLat(v) {
    return Math.abs(v).toFixed(2) + "°" + (v >= 0 ? "N" : "S");
  }

  function fmtLon(v) {
    return Math.abs(v).toFixed(2) + "°" + (v >= 0 ? "L" : "O");
  }

  // Céu estrelado gerado em tempo de execução (sem arquivo externo); PRNG fixo = mesmo céu.
  function starfieldDataUrl() {
    var c = document.createElement("canvas");
    c.width = 4096;
    c.height = 2048;
    var g = c.getContext("2d");
    g.fillStyle = "#01030a";
    g.fillRect(0, 0, c.width, c.height);
    var seed = 20260929;
    function rnd() {
      seed = (seed * 1664525 + 1013904223) >>> 0;
      return seed / 4294967296;
    }
    for (var i = 0; i < 9000; i++) {
      var x = rnd() * c.width;
      var y = rnd() * c.height;
      var m = Math.pow(rnd(), 6);
      var r = 0.4 + 1.6 * m;
      var tint = rnd();
      var col = tint < 0.15 ? "170,200,255" : tint > 0.92 ? "255,220,180" : "255,255,255";
      g.fillStyle = "rgba(" + col + "," + (0.35 + 0.65 * m) + ")";
      g.beginPath();
      g.arc(x, y, r, 0, 2 * Math.PI);
      g.fill();
      if (m > 0.55) {
        var glow = g.createRadialGradient(x, y, 0, x, y, r * 5);
        glow.addColorStop(0, "rgba(" + col + ",0.35)");
        glow.addColorStop(1, "rgba(" + col + ",0)");
        g.fillStyle = glow;
        g.fillRect(x - r * 5, y - r * 5, r * 10, r * 10);
      }
    }
    return c.toDataURL("image/png");
  }

  // --- globo 3D ---------------------------------------------------------------------------------
  var globeEl = document.getElementById("globe");

  sats.forEach(function (s) {
    var el = document.createElement("div");
    el.className = "sat-marker";
    el.style.color = s.color;
    el.innerHTML = '<div class="dot"></div><div class="label"></div>';
    el.querySelector(".label").textContent = s.name;
    s.marker = { kind: "sat", sat: s, el: el, lat: 0, lng: 0, alt: 0 };
  });
  var stationMarkers = stations.map(function (st) {
    var el = document.createElement("div");
    el.className = "station-marker";
    el.innerHTML = '<div class="label"></div>';
    el.querySelector(".label").textContent = "◉ " + st.name;
    return { kind: "station", el: el, lat: st.lat_deg, lng: st.lon_deg, alt: 0.002 };
  });

  var globe = Globe()(globeEl)
    .globeImageUrl(TEXTURE)
    .backgroundImageUrl(starfieldDataUrl())
    .showAtmosphere(true)
    .atmosphereColor("#4ea8ff")
    .atmosphereAltitude(0.2)
    .htmlElementsData([])
    .htmlLat("lat")
    .htmlLng("lng")
    .htmlAltitude("alt")
    .htmlElement(function (d) {
      return d.el;
    })
    .htmlElementVisibilityModifier(function (el, isVisible) {
      el.style.opacity = isVisible ? "1" : "0";
    })
    .htmlTransitionDuration(0)
    .pathsData([])
    .pathPoints("pts")
    .pathPointLat(function (p) {
      return p[0];
    })
    .pathPointLng(function (p) {
      return p[1];
    })
    .pathPointAlt(function (p) {
      return p[2] / EARTH_RADIUS_KM;
    })
    .pathColor(function (d) {
      return d.colors;
    })
    .pathStroke(function (d) {
      return d.stroke;
    })
    .pathDashLength(function (d) {
      return d.dash;
    })
    .pathDashGap(function (d) {
      return d.gap;
    })
    .pathTransitionDuration(0)
    .ringsData(stations)
    .ringLat("lat_deg")
    .ringLng("lon_deg")
    .ringColor(function () {
      return function (t) {
        return "rgba(255,77,141," + (1 - t) + ")";
      };
    })
    .ringMaxRadius(4)
    .ringPropagationSpeed(2)
    .ringRepeatPeriod(1400);

  if (stations.length) {
    globe.pointOfView({ lat: stations[0].lat_deg + 8, lng: stations[0].lon_deg, altitude: 2.4 });
  }

  function update3D(t) {
    var markers = [];
    var paths = [];
    // Caminhos com menos de 2 pontos (ex.: trilha em t = 0) quebram o gradiente do globe.gl.
    function addPath(path) {
      if (path.pts.length >= 2) paths.push(path);
    }
    sats.forEach(function (s) {
      if (!s.visible) return;
      var p = s.now;
      s.marker.lat = p[0];
      s.marker.lng = p[1];
      s.marker.alt = p[2] / EARTH_RADIUS_KM;
      markers.push(s.marker);
      // Cor única por caminho (gradiente por vértice custa caro a cada quadro no globe.gl).
      addPath({
        pts: track(s.lla, t - TRAIL_3D_S, t),
        colors: s.color,
        stroke: 2.2,
        dash: 1,
        gap: 0,
      });
      addPath({
        pts: track(s.lla, t, t + AHEAD_3D_S, STRIDE_3D_AHEAD),
        colors: rgba(s.color, 0.55),
        stroke: null,
        dash: 0.012,
        gap: 0.008,
      });
    });
    // S3-06: linha de visada estação → satélite selecionado, colorida pelo estado; só quando o
    // satélite está acima do horizonte (abaixo dele a reta atravessaria a Terra).
    var sel = state.selected;
    if (STATION && sel && sel.visible && sel.looks) {
      var look = lookAt(sel, t);
      if (look.state !== "down") {
        addPath({
          pts: chordPoints(STATION, sel.now),
          colors: VIS_COLORS[look.state],
          stroke: look.state === "up" ? 1.6 : 0.9,
          dash: 1,
          gap: 0,
        });
      }
    }
    globe.htmlElementsData(markers.concat(stationMarkers));
    globe.pathsData(paths);
  }

  // Pontos da reta (corda) entre a estação e o satélite, em (lat, lon, alt_km) para o globe.gl.
  // Geometria de desenho no globo esférico do viewer: interpola em 3D e volta para lat/lon/alt.
  function chordPoints(st, satNow) {
    function xyz(lat, lon, r) {
      var la = lat * D2R;
      var lo = lon * D2R;
      return [r * Math.cos(la) * Math.cos(lo), r * Math.cos(la) * Math.sin(lo), r * Math.sin(la)];
    }
    var a = xyz(st.lat_deg, st.lon_deg, 1 + st.alt_m / 1000 / EARTH_RADIUS_KM);
    var b = xyz(satNow[0], satNow[1], 1 + satNow[2] / EARTH_RADIUS_KM);
    var pts = [];
    for (var i = 0; i <= LOS_POINTS; i++) {
      var u = i / LOS_POINTS;
      var x = a[0] + u * (b[0] - a[0]);
      var y = a[1] + u * (b[1] - a[1]);
      var z = a[2] + u * (b[2] - a[2]);
      var r = Math.sqrt(x * x + y * y + z * z);
      pts.push([Math.asin(z / r) / D2R, Math.atan2(y, x) / D2R, Math.max(0, (r - 1) * EARTH_RADIUS_KM)]);
    }
    return pts;
  }

  // --- mapa 2D (equiretangular, mesma textura) ------------------------------------------------
  var mapCanvas = document.getElementById("map");
  var ctx = mapCanvas.getContext("2d");
  var bgCanvas = document.createElement("canvas");
  var texture = new Image();
  texture.onload = function () {
    renderMapBackground();
  };
  texture.src = TEXTURE;

  function mapX(lon, w) {
    return ((lon + 180) / 360) * w;
  }

  function mapY(lat, h) {
    return ((90 - lat) / 180) * h;
  }

  function resizeMap() {
    var dpr = window.devicePixelRatio || 1;
    var w = Math.max(1, Math.round(mapCanvas.clientWidth * dpr));
    var h = Math.max(1, Math.round(mapCanvas.clientHeight * dpr));
    if (mapCanvas.width !== w || mapCanvas.height !== h) {
      mapCanvas.width = w;
      mapCanvas.height = h;
      renderMapBackground();
    }
  }

  // Fundo pré-renderizado: textura escurecida + grade de lat/lon a cada 30°.
  function renderMapBackground() {
    var w = mapCanvas.width;
    var h = mapCanvas.height;
    bgCanvas.width = w;
    bgCanvas.height = h;
    var g = bgCanvas.getContext("2d");
    g.fillStyle = "#06101d";
    g.fillRect(0, 0, w, h);
    if (texture.complete && texture.naturalWidth) g.drawImage(texture, 0, 0, w, h);
    g.fillStyle = "rgba(3,8,18,0.38)";
    g.fillRect(0, 0, w, h);
    var dpr = window.devicePixelRatio || 1;
    g.lineWidth = 1 * dpr;
    g.font = 10 * dpr + "px DejaVu Sans, sans-serif";
    g.fillStyle = "rgba(200,215,235,0.55)";
    for (var lon = -180; lon <= 180; lon += 30) {
      g.strokeStyle = lon === 0 ? "rgba(200,215,235,0.35)" : "rgba(200,215,235,0.16)";
      g.beginPath();
      g.moveTo(mapX(lon, w), 0);
      g.lineTo(mapX(lon, w), h);
      g.stroke();
      if (lon > -180 && lon < 180) g.fillText(lon + "°", mapX(lon, w) + 3 * dpr, h - 4 * dpr);
    }
    for (var lat = -60; lat <= 60; lat += 30) {
      g.strokeStyle = lat === 0 ? "rgba(200,215,235,0.35)" : "rgba(200,215,235,0.16)";
      g.beginPath();
      g.moveTo(0, mapY(lat, h));
      g.lineTo(w, mapY(lat, h));
      g.stroke();
      g.fillText(lat + "°", 3 * dpr, mapY(lat, h) - 3 * dpr);
    }
  }

  // Polilinha com quebra no antimeridiano: ao cruzar ±180°, termina o trecho na borda e
  // recomeça na borda oposta, na latitude do cruzamento (sem linha atravessando o mapa).
  function strokeTrack(pts, w, h) {
    ctx.beginPath();
    ctx.moveTo(mapX(pts[0][1], w), mapY(pts[0][0], h));
    for (var i = 1; i < pts.length; i++) {
      var a = pts[i - 1];
      var b = pts[i];
      var d = b[1] - a[1];
      if (Math.abs(d) > 180) {
        var edge = d > 0 ? -180 : 180;
        var bUnwrapped = b[1] + (d > 0 ? -360 : 360);
        var frac = (edge - a[1]) / (bUnwrapped - a[1]);
        var latC = a[0] + frac * (b[0] - a[0]);
        ctx.lineTo(mapX(edge, w), mapY(latC, h));
        ctx.moveTo(mapX(-edge, w), mapY(latC, h));
      }
      ctx.lineTo(mapX(b[1], w), mapY(b[0], h));
    }
    ctx.stroke();
  }

  // Rótulo à direita do ponto; se não couber no mapa, vai para a esquerda. ``placed`` guarda os
  // rótulos já desenhados no quadro: um rótulo que colidiria com outro desce uma linha (ex.:
  // dois satélites na mesma posição, como a órbita com e sem J2 no início da cena).
  function drawLabel(text, x, y, offset, w, placed) {
    var width = ctx.measureText(text).width;
    var lineH = parseFloat(ctx.font.match(/([\d.]+)px/)[1]) * 1.25;
    var tx = x + offset + width > w ? x - offset - width : x + offset;
    if (placed) {
      var moved = true;
      while (moved) {
        moved = false;
        for (var i = 0; i < placed.length; i++) {
          var p = placed[i];
          // folga de 1 %: depois de mover, |y − p.y| = lineH exato pode sair um ulp abaixo
          if (tx < p.x + p.w && tx + width > p.x && Math.abs(y - p.y) < 0.99 * lineH) {
            y = p.y + lineH;
            moved = true;
          }
        }
      }
      placed.push({ x: tx, y: y, w: width });
    }
    ctx.fillText(text, tx, y);
  }

  // S3-06: círculo de visibilidade (footprint) em torno do subponto. Cada vértice é o ponto à
  // distância angular λ do centro na direção θ (geometria da esfera). As longitudes são
  // desembrulhadas em torno do centro e o polígono é desenhado três vezes (deslocado em −360°,
  // 0 e +360°): o canvas recorta o que sai do mapa, e o círculo que cruza o antimeridiano
  // aparece inteiro, metade em cada borda. Se o círculo contém um polo, o contorno percorre todas
  // as longitudes e o preenchimento fecha pela borda do polo.
  function footprintRing(lat0, lon0, lamDeg) {
    var p1 = lat0 * D2R;
    var d = lamDeg * D2R;
    var pts = [];
    for (var k = 0; k <= FOOTPRINT_POINTS; k++) {
      var th = (2 * Math.PI * k) / FOOTPRINT_POINTS;
      var p2 = Math.asin(Math.sin(p1) * Math.cos(d) + Math.cos(p1) * Math.sin(d) * Math.cos(th));
      var dl = Math.atan2(
        Math.sin(th) * Math.sin(d) * Math.cos(p1),
        Math.cos(d) - Math.sin(p1) * Math.sin(p2)
      );
      var lon = lon0 + dl / D2R;
      if (k > 0) {
        var prev = pts[k - 1][1];
        while (lon - prev > 180) lon -= 360;
        while (lon - prev < -180) lon += 360;
      }
      pts.push([p2 / D2R, lon]);
    }
    return pts;
  }

  function drawFootprint(s, t, w, h, dpr) {
    var lam = sampleArr(s.footprint, t);
    var ring = footprintRing(s.now[0], s.now[1], lam);
    var polar = Math.abs(s.now[0]) + lam > 90;
    var poleLat = s.now[0] >= 0 ? 90 : -90;
    var strong = s === state.selected;
    [-360, 0, 360].forEach(function (shift) {
      ctx.beginPath();
      ring.forEach(function (p, i) {
        var x = mapX(p[1] + shift, w);
        var y = mapY(p[0], h);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      });
      ctx.strokeStyle = rgba(s.color, strong ? 0.85 : 0.45);
      ctx.lineWidth = (strong ? 1.4 : 1) * dpr;
      ctx.stroke();
      if (polar) {
        var last = ring[ring.length - 1];
        ctx.lineTo(mapX(last[1] + shift, w), mapY(poleLat, h));
        ctx.lineTo(mapX(ring[0][1] + shift, w), mapY(poleLat, h));
      }
      ctx.closePath();
      ctx.fillStyle = rgba(s.color, strong ? 0.12 : 0.06);
      ctx.fill();
    });
  }

  function drawMap(t, wallMs) {
    resizeMap();
    var w = mapCanvas.width;
    var h = mapCanvas.height;
    var dpr = window.devicePixelRatio || 1;
    ctx.drawImage(bgCanvas, 0, 0);
    ctx.lineJoin = "round";
    ctx.lineCap = "round";

    stations.forEach(function (st) {
      var x = mapX(st.lon_deg, w);
      var y = mapY(st.lat_deg, h);
      var pulse = (wallMs % 1400) / 1400;
      ctx.strokeStyle = "rgba(255,77,141," + (1 - pulse) + ")";
      ctx.lineWidth = 1.5 * dpr;
      ctx.beginPath();
      ctx.arc(x, y, (4 + 10 * pulse) * dpr, 0, 2 * Math.PI);
      ctx.stroke();
      ctx.fillStyle = "#ff4d8d";
      ctx.beginPath();
      ctx.arc(x, y, 3.5 * dpr, 0, 2 * Math.PI);
      ctx.fill();
      ctx.font = 11 * dpr + "px DejaVu Sans, sans-serif";
      ctx.fillStyle = "#ffb3cf";
      drawLabel(st.name, x, y + 4 * dpr, 7 * dpr, w);
    });

    var placedLabels = [];
    sats.forEach(function (s) {
      if (s.visible && s.footprint) drawFootprint(s, t, w, h, dpr);
    });
    sats.forEach(function (s) {
      if (!s.visible) return;
      ctx.lineWidth = 2 * dpr;
      ctx.strokeStyle = rgba(s.color, 0.95);
      ctx.setLineDash([]);
      strokeTrack(track(s.lla, t - PAST_2D_S, t, STRIDE_2D), w, h);
      ctx.strokeStyle = rgba(s.color, 0.7);
      ctx.lineWidth = 1.6 * dpr;
      ctx.setLineDash([6 * dpr, 5 * dpr]);
      strokeTrack(track(s.lla, t, t + AHEAD_2D_S, STRIDE_2D), w, h);
      ctx.setLineDash([]);

      var x = mapX(s.now[1], w);
      var y = mapY(s.now[0], h);
      var glow = ctx.createRadialGradient(x, y, 0, x, y, 14 * dpr);
      glow.addColorStop(0, rgba(s.color, 0.9));
      glow.addColorStop(1, rgba(s.color, 0));
      ctx.fillStyle = glow;
      ctx.beginPath();
      ctx.arc(x, y, 14 * dpr, 0, 2 * Math.PI);
      ctx.fill();
      ctx.fillStyle = "#ffffff";
      ctx.beginPath();
      ctx.arc(x, y, 3.5 * dpr, 0, 2 * Math.PI);
      ctx.fill();
      if (s === state.selected) {
        ctx.strokeStyle = "#ffffff";
        ctx.lineWidth = 1.5 * dpr;
        ctx.beginPath();
        ctx.arc(x, y, 8 * dpr, 0, 2 * Math.PI);
        ctx.stroke();
      }
      ctx.font = "600 " + 11.5 * dpr + "px DejaVu Sans, sans-serif";
      ctx.fillStyle = s.color;
      drawLabel(s.name, x, y - 7 * dpr, 9 * dpr, w, placedLabels);
    });
  }

  // --- legenda e controles ----------------------------------------------------------------------
  var legend = document.getElementById("legend");
  sats.forEach(function (s) {
    var li = document.createElement("li");
    li.innerHTML =
      '<input type="checkbox" checked><span class="swatch"></span><span class="name"></span>' +
      '<span class="readout"></span>';
    var cb = li.querySelector("input");
    cb.title = "Mostrar / ocultar " + s.name;
    cb.addEventListener("change", function () {
      s.visible = cb.checked;
      render();
    });
    li.querySelector(".swatch").style.color = s.color;
    li.querySelector(".swatch").style.background = s.color;
    li.querySelector(".name").textContent = s.name;
    li.querySelector(".name").style.color = s.color;
    li.querySelector(".name").title = "Clique para ver as informações de " + s.name;
    s.readoutEl = li.querySelector(".readout");
    s.legendEl = li;
    li.addEventListener("click", function (e) {
      if (e.target === cb) return; // a caixa só mostra/oculta
      selectSat(s);
    });
    legend.appendChild(li);
  });

  // --- painel de informações do satélite selecionado (clique na legenda) ---------------------
  var info = {
    name: document.getElementById("info-name"),
    lat: document.getElementById("info-lat"),
    lon: document.getElementById("info-lon"),
    alt: document.getElementById("info-alt"),
  };
  document.getElementById("info-note").textContent = GEODETIC
    ? "Coordenadas geodésicas WGS84: latitude da vertical local, altitude acima do elipsoide. " +
      "Clique em outro satélite na lista para trocar."
    : "Subponto esférico (cena antiga). Clique em outro satélite na lista para trocar.";

  function selectSat(s) {
    state.selected = s;
    lastReadout = 0; // atualiza o painel já no próximo quadro
    sats.forEach(function (o) {
      o.legendEl.classList.toggle("selected", o === s);
    });
    buildPassList();
    render();
  }

  function updateInfo() {
    var s = state.selected;
    if (!s || !s.now) return;
    info.name.textContent = s.name;
    info.name.style.color = s.color;
    info.lat.textContent = Math.abs(s.now[0]).toFixed(3) + "° " + (s.now[0] >= 0 ? "N" : "S");
    info.lon.textContent = Math.abs(s.now[1]).toFixed(3) + "° " + (s.now[1] >= 0 ? "L" : "O");
    info.alt.textContent = s.now[2].toFixed(2) + " km";
  }

  // --- S3-06: painel "satélite × estação", lista de passagens e faixa do TLE ------------------
  var stationPanel = document.getElementById("station-panel");
  var stEls = {
    title: document.getElementById("station-title"),
    state: document.getElementById("st-state"),
    el: document.getElementById("st-el"),
    az: document.getElementById("st-az"),
    range: document.getElementById("st-range"),
    next: document.getElementById("st-next"),
    passesTitle: document.getElementById("passes-title"),
    passes: document.getElementById("passes"),
    marks: document.getElementById("tl-marks"),
    banner: document.getElementById("tle-banner"),
  };
  var STATE_TEXT = {
    up: function () {
      return "visível · acima de " + fmtNum(EL_MIN, 0) + "°";
    },
    low: function () {
      return "baixo · entre 0° e " + fmtNum(EL_MIN, 0) + "°";
    },
    down: function () {
      return "abaixo do horizonte";
    },
  };
  if (HAS_STATION_VIEW) {
    stationPanel.hidden = false;
    Array.prototype.forEach.call(document.querySelectorAll(".k-optional"), function (el) {
      el.hidden = false;
    });
  }

  // Lista de passagens do satélite selecionado sobre a estação; clicar leva o relógio a AOS − 60 s.
  function buildPassList() {
    if (!HAS_STATION_VIEW) return;
    var s = state.selected;
    stEls.title.textContent = (s ? s.name : "—") + " × " + STATION.name;
    stEls.passes.innerHTML = "";
    stEls.marks.innerHTML = "";
    var list = s ? s.passes : [];
    stEls.passesTitle.textContent =
      list.length +
      " passage" +
      (list.length === 1 ? "m" : "ns") +
      " acima de " +
      fmtNum(EL_MIN, 0) +
      "° na janela · clique para ir ao AOS − 60 s · horário UTC e de Brasília (BRT)";
    list.forEach(function (p, i) {
      var li = document.createElement("li");
      li.innerHTML =
        '<span class="n"></span><span class="when"></span><span class="el"></span>' +
        '<span class="dur"></span>';
      li.querySelector(".n").textContent = i + 1;
      li.querySelector(".when").textContent =
        fmtDayTime(p.aosMs, 0) + " UTC · " + fmtDayTime(p.aosMs, BRT_OFFSET_MS).slice(6, 11) + " BRT";
      li.querySelector(".el").textContent = "máx " + fmtNum(p.maxEl, 1) + "°";
      li.querySelector(".dur").textContent = fmtNum(p.duration / 60, 1) + " min";
      li.title = "AOS → LOS: " + fmtDayTime(p.aosMs, 0) + " → " +
        fmtDayTime(p.aosMs + p.duration * 1000, 0) + " UTC";
      li.addEventListener("click", function () {
        state.simTime = Math.max(0, Math.min(DURATION, p.aos - 60));
        jumpRender();
      });
      p.li = li;
      stEls.passes.appendChild(li);
      var mark = document.createElement("span");
      mark.style.left = (100 * Math.max(0, p.aos)) / DURATION + "%";
      mark.style.width = (100 * (Math.min(DURATION, p.los) - Math.max(0, p.aos))) / DURATION + "%";
      mark.title = "Passagem " + (i + 1) + " · máx " + fmtNum(p.maxEl, 1) + "°";
      stEls.marks.appendChild(mark);
    });
  }

  function updateStationPanel(t) {
    if (!HAS_STATION_VIEW) return;
    var s = state.selected;
    var look = s ? lookAt(s, t) : null;
    if (!look) {
      stEls.state.className = "st-state down";
      stEls.state.textContent = "sem dados de visibilidade";
      return;
    }
    stEls.state.className = "st-state " + look.state;
    stEls.state.textContent = STATE_TEXT[look.state]();
    stEls.el.textContent = fmtNum(look.el, 1) + "°";
    stEls.az.textContent = fmtNum(look.az, 1) + "°";
    stEls.range.textContent = Math.round(look.range) + " km";
    var current = null;
    var next = null;
    s.passes.forEach(function (p) {
      if (p.aos <= t && t <= p.los) current = p;
      else if (p.aos > t && (!next || p.aos < next.aos)) next = p;
      if (p.li) {
        p.li.classList.toggle("now", p.aos <= t && t <= p.los);
        p.li.classList.toggle("past", p.los < t);
      }
    });
    s.passes.forEach(function (p) {
      if (p.li) p.li.classList.toggle("next", p === next && !current);
    });
    if (current) {
      stEls.next.innerHTML = "Passagem em andamento · termina em <b></b>";
      stEls.next.querySelector("b").textContent = fmtDuration(current.los - t);
    } else if (next) {
      stEls.next.innerHTML = "Próxima passagem em <b></b> · ";
      stEls.next.querySelector("b").textContent = fmtDuration(next.aos - t);
      stEls.next.appendChild(
        document.createTextNode(
          "AOS " + fmtDayTime(next.aosMs, BRT_OFFSET_MS).slice(6, 11) + " BRT, máx " +
            fmtNum(next.maxEl, 1) + "°"
        )
      );
    } else {
      stEls.next.textContent = "Sem mais passagens nesta janela.";
    }
  }

  var SOURCE_TEXT = {
    network: "CelesTrak (rede)",
    cache: "cache local",
    bundled: "TLE embarcado no pacote",
  };
  function updateTleBanner(t) {
    if (!TLE) return;
    var age = TLE.age_days_at_start + t / 86400;
    var stale = age > TLE_STALE_DAYS;
    stEls.banner.hidden = false;
    stEls.banner.classList.toggle("stale", stale);
    stEls.banner.innerHTML =
      "<strong></strong> · época do TLE <strong></strong> · idade <strong></strong> · fonte: " +
      "<span></span><span></span>";
    var parts = stEls.banner.querySelectorAll("strong");
    parts[0].textContent = (TLE.name || "satélite " + TLE.norad_id) + " real (TLE/SGP4)";
    parts[1].textContent = fmtDayTime(Date.parse(TLE.epoch_utc), 0) + " UTC";
    parts[2].textContent = fmtNum(age, 1) + " dias";
    var spans = stEls.banner.querySelectorAll("span");
    spans[0].textContent = SOURCE_TEXT[TLE.source] || TLE.source || "—";
    spans[1].textContent = stale
      ? " · ⚠ TLE com mais de " + TLE_STALE_DAYS + " dias: a posição real pode diferir de dezenas de km"
      : "";
  }

  var clockEl = document.getElementById("clock");
  var playBtn = document.getElementById("btn-play");
  var timeline = document.getElementById("timeline");
  timeline.max = String(DURATION);
  timeline.step = String(Math.min(1, STEP));
  document.getElementById("tl-start").textContent = fmtUtc(EPOCH_MS).slice(0, 16);
  document.getElementById("tl-end").textContent = fmtUtc(EPOCH_MS + DURATION * 1000).slice(0, 16);
  document.getElementById("frame-note").textContent =
    "Vista com a Terra fixa (ECEF) · subponto " +
    (GEODETIC ? "geodésico WGS84" : "esférico") +
    " · janela de " +
    (DURATION / 3600).toFixed(1).replace(".", ",") +
    " h com passo de " +
    STEP +
    " s · propagação em Python (" +
    sats.length +
    " satélite" +
    (sats.length === 1 ? "" : "s") +
    ")";

  var speedsEl = document.getElementById("speeds");
  SPEEDS.forEach(function (v) {
    var b = document.createElement("button");
    b.className = "btn" + (v === state.speed ? " active" : "");
    b.textContent = v + "x";
    b.title = "Velocidade " + v + "x";
    b.addEventListener("click", function () {
      state.speed = v;
      Array.prototype.forEach.call(speedsEl.children, function (c) {
        c.classList.toggle("active", c === b);
      });
    });
    speedsEl.appendChild(b);
  });

  function setPlaying(p) {
    state.playing = p;
    playBtn.textContent = p ? "❚❚" : "▶";
  }
  playBtn.addEventListener("click", function () {
    setPlaying(!state.playing);
  });
  document.getElementById("btn-start").addEventListener("click", function () {
    state.simTime = 0;
    jumpRender();
  });
  timeline.addEventListener("pointerdown", function () {
    state.dragging = true;
  });
  window.addEventListener("pointerup", function () {
    state.dragging = false;
  });
  timeline.addEventListener("input", function () {
    state.simTime = Number(timeline.value);
    jumpRender();
  });
  window.addEventListener("keydown", function (e) {
    if (e.code === "Space" && e.target === document.body) {
      e.preventDefault();
      setPlaying(!state.playing);
    }
  });

  // --- laço de animação: uma única fonte de tempo (simTime) -----------------------------------
  var lastReadout = 0;

  // Depois de um salto do relógio (passagem, linha do tempo, início), os painéis de texto, que
  // normalmente atualizam a cada 100 ms, são refeitos já neste quadro.
  function jumpRender() {
    lastReadout = 0;
    render();
  }

  function render(wallMs) {
    wallMs = wallMs || performance.now();
    var t = state.simTime;
    sats.forEach(function (s) {
      s.now = sampleAt(s.lla, t);
    });
    update3D(t);
    drawMap(t, wallMs);
    clockEl.textContent = fmtUtc(EPOCH_MS + t * 1000);
    if (!state.dragging) timeline.value = String(t);
    if (wallMs - lastReadout > 100) {
      lastReadout = wallMs;
      sats.forEach(function (s) {
        s.readoutEl.textContent =
          fmtLat(s.now[0]) + "  " + fmtLon(s.now[1]) + "  " + s.now[2].toFixed(1) + " km";
      });
      updateInfo();
      updateStationPanel(t);
      updateTleBanner(t);
    }
  }

  function resizeGlobe() {
    globe.width(globeEl.clientWidth).height(globeEl.clientHeight);
  }
  window.addEventListener("resize", function () {
    resizeGlobe();
    resizeMap();
  });
  resizeGlobe();
  if (state.selected) state.selected.legendEl.classList.add("selected");
  buildPassList();

  var lastWall = null;
  function frame(wallMs) {
    if (lastWall !== null && state.playing && !state.dragging) {
      var dt = Math.min(0.25, (wallMs - lastWall) / 1000);
      state.simTime += dt * state.speed;
      if (state.simTime > DURATION) state.simTime = 0; // volta ao início (demo em laço)
    }
    lastWall = wallMs;
    render(wallMs);
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);

  window.SATSIM_VIEWER = { state: state, globe: globe, render: render }; // inspeção/depuração
})();
