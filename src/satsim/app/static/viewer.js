/*
 * Viewer do satsim (ADR 0001, seção GUI).
 *
 * Regras: a física fica em Python. Este arquivo só lê a cena (window.SATSIM_SCENE), interpola
 * no tempo entre amostras e projeta para desenhar. Nenhuma equação orbital aqui.
 * Uma única variável, simTime (segundos desde meta.epoch_utc), comanda o globo 3D e o mapa 2D.
 * Qualquer número de satélites é lido do JSON (R7).
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
  var STRIDE_3D_AHEAD = 3; // órbita prevista no globo: 1 amostra a cada 3 passos (30 s)
  var STRIDE_2D = 2; // traços no mapa: 1 amostra a cada 2 passos (20 s ≈ 1,2° na órbita)
  var TEXTURE = window.SATSIM_TEXTURE || "textures/earth_2048.jpg";

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
    return { id: s.id, name: s.name, color: s.color, lla: s.lla, visible: true, now: null };
  });
  var stations = scene.stations;

  var state = { simTime: 0, playing: true, speed: DEFAULT_SPEED, dragging: false };

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
    globe.htmlElementsData(markers.concat(stationMarkers));
    globe.pathsData(paths);
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

  // Rótulo à direita do ponto; se não couber no mapa, vai para a esquerda.
  function drawLabel(text, x, y, offset, w) {
    var width = ctx.measureText(text).width;
    var tx = x + offset + width > w ? x - offset - width : x + offset;
    ctx.fillText(text, tx, y);
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
      ctx.font = "600 " + 11.5 * dpr + "px DejaVu Sans, sans-serif";
      ctx.fillStyle = s.color;
      drawLabel(s.name, x, y - 7 * dpr, 9 * dpr, w);
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
    s.readoutEl = li.querySelector(".readout");
    legend.appendChild(li);
  });

  var clockEl = document.getElementById("clock");
  var playBtn = document.getElementById("btn-play");
  var timeline = document.getElementById("timeline");
  timeline.max = String(DURATION);
  timeline.step = String(Math.min(1, STEP));
  document.getElementById("tl-start").textContent = fmtUtc(EPOCH_MS).slice(0, 16);
  document.getElementById("tl-end").textContent = fmtUtc(EPOCH_MS + DURATION * 1000).slice(0, 16);
  document.getElementById("frame-note").textContent =
    "Vista com a Terra fixa (ECEF) · subponto esférico · janela de " +
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
    render();
  });
  timeline.addEventListener("pointerdown", function () {
    state.dragging = true;
  });
  window.addEventListener("pointerup", function () {
    state.dragging = false;
  });
  timeline.addEventListener("input", function () {
    state.simTime = Number(timeline.value);
    render();
  });
  window.addEventListener("keydown", function (e) {
    if (e.code === "Space" && e.target === document.body) {
      e.preventDefault();
      setPlaying(!state.playing);
    }
  });

  // --- laço de animação: uma única fonte de tempo (simTime) -----------------------------------
  var lastReadout = 0;

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
