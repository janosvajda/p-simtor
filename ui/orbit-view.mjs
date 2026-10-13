// Orbit view: renders point masses, their computed trajectories and a reference
// grid with an orthographic 3D camera on a 2D canvas.
//
// Scientific conventions (disclosed in the page legend):
//   - Trails are polylines through integrator states only (no splines, no smoothing).
//   - Marker radius ∝ m^(1/3) for legibility; bodies are point masses in the model.
//   - Orthographic projection: lengths in the orbital plane are true to scale.

import { cssVar } from './instrument.mjs';

const TRAIL_CAP = 24000;

export function createOrbitView(canvas, { n: nBodies, colors: initialColors, masses, editable }) {
  let n = nBodies, colors = initialColors;
  const ctx = canvas.getContext('2d');
  const cam = { yaw: 0, pitch: 0, scale: 100, cx: 0, cy: 0, cz: 0 };
  let W = 0, H = 0, dpr = 1;
  let trails = [];
  let twinTrails = [];
  let trailTime = 20;
  let minSeg = 1e-3;
  let massesNow = masses;
  const opts = { grid: true, vectors: false, barycentre: true, labels: true };
  const view = { cam, opts };

  function resize() {
    dpr = devicePixelRatio || 1;
    W = canvas.clientWidth; H = canvas.clientHeight;
    canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
  }
  new ResizeObserver(resize).observe(canvas);
  resize();

  function makeTrails() {
    return Array.from({ length: n }, () => ({ buf: new Float64Array(TRAIL_CAP * 4), head: 0, len: 0 }));
  }
  view.resetTrails = (count = n) => { n = count; trails = makeTrails(); twinTrails = makeTrails(); };
  view.resetTrails();
  view.setTrailTime = t => { trailTime = t; };
  view.setMinSegment = d => { minSeg = d; };
  view.setMasses = m => { massesNow = m; };
  view.setColors = c => { colors = c; };

  function pushTo(set, x, t) {
    for (let i = 0; i < n; i++) {
      const tr = set[i];
      const px = x[3 * i], py = x[3 * i + 1], pz = x[3 * i + 2];
      if (tr.len) {
        const j = ((tr.head - 1 + TRAIL_CAP) % TRAIL_CAP) * 4;
        const d = Math.hypot(px - tr.buf[j], py - tr.buf[j + 1], pz - tr.buf[j + 2]);
        if (d < minSeg) continue;
      }
      const k = tr.head * 4;
      tr.buf[k] = px; tr.buf[k + 1] = py; tr.buf[k + 2] = pz; tr.buf[k + 3] = t;
      tr.head = (tr.head + 1) % TRAIL_CAP;
      tr.len = Math.min(TRAIL_CAP, tr.len + 1);
    }
  }
  view.pushTrail = (x, t) => pushTo(trails, x, t);
  view.pushTwinTrail = (x, t) => pushTo(twinTrails, x, t);

  // ---------- camera ----------
  function project(wx, wy, wz, out) {
    const x = wx - cam.cx, y = wy - cam.cy, z = wz - cam.cz;
    const cy = Math.cos(cam.yaw), sy = Math.sin(cam.yaw);
    const x1 = cy * x - sy * y, y1 = sy * x + cy * y;
    const cp = Math.cos(cam.pitch), sp = Math.sin(cam.pitch);
    const y2 = cp * y1 - sp * z, z2 = sp * y1 + cp * z;
    out[0] = W / 2 + x1 * cam.scale;
    out[1] = H / 2 - y2 * cam.scale;
    out[2] = z2;
    return out;
  }
  view.project = project;

  /** Inverse projection onto the z = cam.cz plane (valid for any yaw when pitch = 0). */
  view.unproject = (sx, sy) => {
    const x1 = (sx - W / 2) / cam.scale, y1 = -(sy - H / 2) / cam.scale;
    const c = Math.cos(cam.yaw), s = Math.sin(cam.yaw);
    return [c * x1 + s * y1 + cam.cx, -s * x1 + c * y1 + cam.cy, cam.cz];
  };

  view.fit = (x, { animate = false } = {}) => {
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity, z0 = Infinity, z1 = -Infinity;
    for (let i = 0; i < n; i++) {
      x0 = Math.min(x0, x[3 * i]); x1 = Math.max(x1, x[3 * i]);
      y0 = Math.min(y0, x[3 * i + 1]); y1 = Math.max(y1, x[3 * i + 1]);
      z0 = Math.min(z0, x[3 * i + 2]); z1 = Math.max(z1, x[3 * i + 2]);
    }
    const span = Math.max(x1 - x0, y1 - y0, z1 - z0, 1e-9);
    const target = {
      cx: (x0 + x1) / 2, cy: (y0 + y1) / 2, cz: (z0 + z1) / 2,
      scale: 0.62 * Math.min(W, H) / span,
    };
    if (!animate) Object.assign(cam, target);
    else view._fitTarget = target;
  };

  // ---------- drawing ----------
  const P = new Float64Array(3), Q = new Float64Array(3);

  function niceStep(raw) {
    const e = Math.pow(10, Math.floor(Math.log10(raw)));
    const f = raw / e;
    return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * e;
  }
  view.niceStep = niceStep;

  function drawGrid() {
    const step = niceStep(70 / cam.scale);
    const reach = Math.hypot(W, H) / cam.scale;
    const i0 = Math.floor((cam.cx - reach) / step), i1 = Math.ceil((cam.cx + reach) / step);
    const j0 = Math.floor((cam.cy - reach) / step), j1 = Math.ceil((cam.cy + reach) / step);
    const minor = cssVar('--grid'), major = cssVar('--grid-strong');
    ctx.lineWidth = 1;
    const line = (ax, ay, bx, by, strong) => {
      project(ax, ay, 0, P); project(bx, by, 0, Q);
      ctx.strokeStyle = strong ? major : minor;
      ctx.beginPath(); ctx.moveTo(P[0], P[1]); ctx.lineTo(Q[0], Q[1]); ctx.stroke();
    };
    if (i1 - i0 > 400 || j1 - j0 > 400) return step;
    for (let i = i0; i <= i1; i++) line(i * step, j0 * step, i * step, j1 * step, i === 0 || i % 5 === 0);
    for (let j = j0; j <= j1; j++) line(i0 * step, j * step, i1 * step, j * step, j === 0 || j % 5 === 0);
    return step;
  }

  function drawTrail(tr, color, tNow, { dashed = false, alpha = 1 } = {}) {
    if (tr.len < 2) return;
    const tMin = tNow - trailTime;
    const CH = 28;
    // Collect visible points oldest -> newest that are not in the future of tNow.
    const pts = [];
    for (let k = tr.len - 1; k >= 0; k--) {
      const idx = ((tr.head - 1 - k + TRAIL_CAP * 2) % TRAIL_CAP) * 4;
      const t = tr.buf[idx + 3];
      if (t < tMin || t > tNow) continue;
      pts.push(idx);
    }
    if (pts.length < 2) return;
    ctx.strokeStyle = color;
    ctx.lineCap = 'round'; ctx.lineJoin = 'round';
    if (dashed) ctx.setLineDash([2, 4]);
    const per = Math.ceil(pts.length / CH);
    for (let c = 0; c < CH; c++) {
      const a = c * per, b = Math.min(pts.length - 1, (c + 1) * per);
      if (a >= b) break;
      const age = 1 - (c + 0.5) / CH;
      ctx.globalAlpha = alpha * Math.pow(1 - age, 1.6) * 0.95 + 0.03;
      ctx.lineWidth = 0.6 + 1.1 * (1 - age);
      ctx.beginPath();
      for (let k = a; k <= b; k++) {
        const idx = pts[k];
        project(tr.buf[idx], tr.buf[idx + 1], tr.buf[idx + 2], P);
        k === a ? ctx.moveTo(P[0], P[1]) : ctx.lineTo(P[0], P[1]);
      }
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
    ctx.setLineDash([]);
  }

  function radiusFor(m) {
    const mMax = Math.max(...massesNow);
    return 3.2 + 6.5 * Math.cbrt(m / mMax);
  }

  view.draw = ({ x, v, t, twinX, labels, hover }) => {
    const glow = cssVar('--glow') === '1';
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    if (view._fitTarget) {
      const f = view._fitTarget, k = 0.12;
      cam.cx += (f.cx - cam.cx) * k; cam.cy += (f.cy - cam.cy) * k; cam.cz += (f.cz - cam.cz) * k;
      cam.scale *= Math.pow(f.scale / cam.scale, k);
      if (Math.abs(Math.log(f.scale / cam.scale)) < 1e-3) view._fitTarget = null;
    }
    const gridStep = opts.grid ? drawGrid() : niceStep(70 / cam.scale);
    view.gridStep = gridStep;

    const col = i => cssVar(colors[i % colors.length]);

    if (twinX) for (let i = 0; i < n; i++) drawTrail(twinTrails[i], col(i), t, { dashed: true, alpha: 0.55 });
    for (let i = 0; i < n; i++) drawTrail(trails[i], col(i), t);

    if (opts.barycentre) {
      project(0, 0, 0, P);
      ctx.strokeStyle = cssVar('--muted'); ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(P[0] - 5, P[1]); ctx.lineTo(P[0] + 5, P[1]); ctx.moveTo(P[0], P[1] - 5); ctx.lineTo(P[0], P[1] + 5); ctx.stroke();
    }

    // Depth-sort bodies so nearer ones paint last.
    const order = [...Array(n).keys()].map(i => { project(x[3 * i], x[3 * i + 1], x[3 * i + 2], P); return [i, P[0], P[1], P[2]]; })
      .sort((a, b) => a[3] - b[3]);

    if (twinX) {
      for (let i = 0; i < n; i++) {
        project(twinX[3 * i], twinX[3 * i + 1], twinX[3 * i + 2], Q);
        ctx.strokeStyle = col(i); ctx.lineWidth = 1.2; ctx.globalAlpha = 0.8;
        ctx.beginPath(); ctx.arc(Q[0], Q[1], radiusFor(massesNow[i]) + 3, 0, Math.PI * 2); ctx.stroke();
        ctx.globalAlpha = 1;
      }
    }

    for (const [i, sx, sy] of order) {
      if (!Number.isFinite(sx) || !Number.isFinite(sy) || Math.abs(sx) > 1e6 || Math.abs(sy) > 1e6) continue;   // off-screen or diverged
      const r = radiusFor(massesNow[i]);
      const c = col(i);
      if (glow) {
        ctx.globalCompositeOperation = 'lighter';
        const g = ctx.createRadialGradient(sx, sy, 0, sx, sy, r * 6);
        g.addColorStop(0, withAlpha(c, 0.55));
        g.addColorStop(0.25, withAlpha(c, 0.18));
        g.addColorStop(1, withAlpha(c, 0));
        ctx.fillStyle = g;
        ctx.beginPath(); ctx.arc(sx, sy, r * 6, 0, Math.PI * 2); ctx.fill();
        ctx.globalCompositeOperation = 'source-over';
        ctx.fillStyle = c;
        ctx.beginPath(); ctx.arc(sx, sy, r, 0, Math.PI * 2); ctx.fill();
        ctx.fillStyle = 'rgba(255,255,255,0.85)';
        ctx.beginPath(); ctx.arc(sx - r * 0.25, sy - r * 0.25, r * 0.42, 0, Math.PI * 2); ctx.fill();
      } else {
        ctx.fillStyle = withAlpha(c, 0.16);
        ctx.beginPath(); ctx.arc(sx, sy, r + 4, 0, Math.PI * 2); ctx.fill();
        ctx.fillStyle = c;
        ctx.beginPath(); ctx.arc(sx, sy, r, 0, Math.PI * 2); ctx.fill();
        ctx.strokeStyle = cssVar('--ink'); ctx.lineWidth = 1; ctx.globalAlpha = 0.35;
        ctx.stroke(); ctx.globalAlpha = 1;
      }
      if (opts.vectors && v) {
        const k = 0.25 * Math.min(W, H) / Math.max(1e-12, maxSpeed(v)) / cam.scale;
        project(x[3 * i] + v[3 * i] * k, x[3 * i + 1] + v[3 * i + 1] * k, x[3 * i + 2] + v[3 * i + 2] * k, Q);
        arrow(sx, sy, Q[0], Q[1], c);
      }
      if (opts.labels && labels) {
        ctx.font = '500 11px ui-sans-serif, -apple-system, system-ui, sans-serif';
        ctx.fillStyle = cssVar(hover === i ? '--ink' : '--ink-2');
        ctx.fillText(labels[i], sx + r + 6, sy - r - 2);
      }
    }
  };

  function maxSpeed(v) {
    let s = 0;
    for (let i = 0; i < n; i++) s = Math.max(s, Math.hypot(v[3 * i], v[3 * i + 1], v[3 * i + 2]));
    return s;
  }

  function arrow(x0, y0, x1, y1, color) {
    const a = Math.atan2(y1 - y0, x1 - x0), L = Math.hypot(x1 - x0, y1 - y0);
    if (L < 4) return;
    ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = 1.3;
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(x1, y1);
    ctx.lineTo(x1 - 7 * Math.cos(a - 0.4), y1 - 7 * Math.sin(a - 0.4));
    ctx.lineTo(x1 - 7 * Math.cos(a + 0.4), y1 - 7 * Math.sin(a + 0.4));
    ctx.closePath(); ctx.fill();
  }

  /** Screen-space hit test; returns body index or -1. */
  view.pick = (sx, sy, x) => {
    let best = -1, bd = 14;
    for (let i = 0; i < n; i++) {
      project(x[3 * i], x[3 * i + 1], x[3 * i + 2], P);
      const d = Math.hypot(P[0] - sx, P[1] - sy);
      if (d < bd) { bd = d; best = i; }
    }
    return best;
  };

  // ---------- interaction ----------
  // drag: pan · shift/right-drag: rotate · wheel: zoom about cursor · drag body (when editable): move it
  let drag = null;
  canvas.classList.add('grab');
  canvas.addEventListener('contextmenu', e => e.preventDefault());
  canvas.addEventListener('pointerdown', e => {
    canvas.setPointerCapture(e.pointerId);
    const r = canvas.getBoundingClientRect();
    const sx = e.clientX - r.left, sy = e.clientY - r.top;
    const body = editable?.canEdit() && Math.abs(cam.pitch) < 1e-9 ? view.pick(sx, sy, editable.positions()) : -1;
    drag = { sx, sy, mode: body >= 0 ? 'body' : (e.shiftKey || e.button === 2 ? 'rotate' : 'pan'), body, cam: { ...cam } };
    canvas.classList.replace('grab', drag.mode === 'body' ? 'move' : 'grabbing');
  });
  canvas.addEventListener('pointermove', e => {
    const r = canvas.getBoundingClientRect();
    const sx = e.clientX - r.left, sy = e.clientY - r.top;
    if (!drag) {
      if (editable?.canEdit() && Math.abs(cam.pitch) < 1e-9) {
        const hit = view.pick(sx, sy, editable.positions());
        canvas.classList.toggle('move', hit >= 0);
        canvas.classList.toggle('grab', hit < 0);
      }
      return;
    }
    const dx = sx - drag.sx, dy = sy - drag.sy;
    if (drag.mode === 'pan') {
      const c = Math.cos(cam.yaw), s = Math.sin(cam.yaw);
      const wx = -dx / cam.scale, wy = dy / cam.scale;
      cam.cx = drag.cam.cx + c * wx + s * wy;
      cam.cy = drag.cam.cy - s * wx + c * wy;
      view._fitTarget = null;
    } else if (drag.mode === 'rotate') {
      cam.yaw = drag.cam.yaw + dx * 0.006;
      cam.pitch = Math.max(-Math.PI / 2, Math.min(Math.PI / 2, drag.cam.pitch - dy * 0.006));
    } else if (drag.mode === 'body') {
      const p = view.unproject(sx, sy);
      editable.move(drag.body, p[0], p[1]);
    }
  });
  const end = () => {
    if (drag?.mode === 'body') editable.commit();
    drag = null;
    canvas.classList.remove('grabbing', 'move');
    canvas.classList.add('grab');
  };
  canvas.addEventListener('pointerup', end);
  canvas.addEventListener('pointercancel', end);
  canvas.addEventListener('wheel', e => {
    e.preventDefault();
    const r = canvas.getBoundingClientRect();
    const sx = e.clientX - r.left, sy = e.clientY - r.top;
    const before = view.unproject(sx, sy);
    cam.scale *= Math.exp(-e.deltaY * 0.0015);
    const after = view.unproject(sx, sy);
    if (Math.abs(cam.pitch) < 1e-9) { cam.cx += before[0] - after[0]; cam.cy += before[1] - after[1]; }
    view._fitTarget = null;
  }, { passive: false });
  canvas.addEventListener('dblclick', () => { cam.yaw = 0; cam.pitch = 0; view.onReset?.(); });

  return view;
}

const alphaCache = new Map();
function withAlpha(color, a) {
  const key = color + a;
  let out = alphaCache.get(key);
  if (!out) {
    out = `color-mix(in srgb, ${color} ${Math.round(a * 100)}%, transparent)`;
    if (color.startsWith('#') && color.length === 7) {
      const r = parseInt(color.slice(1, 3), 16), g = parseInt(color.slice(3, 5), 16), b = parseInt(color.slice(5, 7), 16);
      out = `rgba(${r},${g},${b},${a})`;
    }
    alphaCache.set(key, out);
  }
  return out;
}
