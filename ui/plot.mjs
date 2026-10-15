// Scientific 2D plot: linear or log axes with 1-2-5 ticks, function and data
// series, shaded validity bands, reference lines, markers, a draggable cursor,
// and highlighting linked to equation terms. Theme-aware (reads CSS tokens).
//
//   const p = createPlot({ title: 'Spectral radiance', x: { label: 'λ', unit: 'nm', min: 100, max: 3000 },
//                          y: { label: 'B_λ', unit: 'W sr⁻¹ m⁻³', log: false } });
//   p.fn('planck', λ => B(λ, T), { color: '--body-0', label: 'T = 5772 K' });
//   p.vline('peak', λmax, { label: 'λ_max' });  p.draw();

import { el, cssVar, sci } from './instrument.mjs';

function niceTicks(min, max, target = 5) {
  const span = max - min;
  if (!(span > 0)) return [min];
  const raw = span / target, e = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / e;
  const step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * e;
  const out = [];
  for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return out;
}

const SUPS = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
/** Powers of ten on log axes read as 10ⁿ. */
function fmtPow10(v) {
  const e = Math.round(Math.log10(v));
  if (e >= 0 && e <= 3) return String(10 ** e);
  return `10${String(e).split('').map(ch => SUPS[ch]).join('')}`;
}

function fmtTick(v) {
  const a = Math.abs(v);
  if (a !== 0 && (a < 1e-3 || a >= 1e5)) return sci(v, 2);
  return String(+v.toPrecision(6)).replace('-', '−');
}

export function createPlot({ title = '', x = {}, y = {}, height = 220, legend = true, onCursor } = {}) {
  const canvas = el('canvas');
  const head = el('div', { class: 'psim-plot-head' }, el('span', { class: 't' }, title), el('span', { class: 'l' }));
  const box = el('div', { class: 'psim-plot2', style: `height:${height}px` }, canvas);
  const root = el('figure', { class: 'psim-figure' }, head, box);
  const ctx = canvas.getContext('2d');
  const items = new Map();
  let order = [], focus = null, cursor = null, dirty = true;
  const ax = { x: { min: 0, max: 1, log: false, ...x }, y: { min: 0, max: 1, log: false, auto: y.min == null || y.max == null, ...y } };
  const pad = { l: 60, r: 14, t: 10, b: 34 };
  let W = 0, H = 0;

  const tr = (a, v) => a.log ? Math.log10(v) : v;
  const X = v => pad.l + (tr(ax.x, v) - tr(ax.x, ax.x.min)) / (tr(ax.x, ax.x.max) - tr(ax.x, ax.x.min)) * (W - pad.l - pad.r);
  const Y = v => H - pad.b - (tr(ax.y, v) - tr(ax.y, ax.y.min)) / (tr(ax.y, ax.y.max) - tr(ax.y, ax.y.min)) * (H - pad.t - pad.b);
  const invX = px => {
    const f = (px - pad.l) / (W - pad.l - pad.r);
    const v = tr(ax.x, ax.x.min) + f * (tr(ax.x, ax.x.max) - tr(ax.x, ax.x.min));
    return ax.x.log ? Math.pow(10, v) : v;
  };

  const add = (id, item) => { if (!items.has(id)) order.push(id); items.set(id, item); dirty = true; return api; };

  const api = {
    el: root,
    /** Plot a function y = f(x), sampled per pixel column. */
    fn(id, f, opts = {}) { return add(id, { kind: 'fn', f, ...opts }); },
    /** Plot data points [[x, y], ...] as a line (or dots with { dots: true }). */
    data(id, pts, opts = {}) { return add(id, { kind: 'data', pts, ...opts }); },
    /** Shaded x-range, e.g. where a formula is invalid. */
    band(id, x0, x1, opts = {}) { return add(id, { kind: 'band', x0, x1, ...opts }); },
    hline(id, yv, opts = {}) { return add(id, { kind: 'hline', y: yv, ...opts }); },
    vline(id, xv, opts = {}) { return add(id, { kind: 'vline', x: xv, ...opts }); },
    marker(id, xv, yv, opts = {}) { return add(id, { kind: 'marker', x: xv, y: yv, ...opts }); },
    remove(id) { items.delete(id); order = order.filter(o => o !== id); dirty = true; return api; },
    setX(o) { Object.assign(ax.x, o); dirty = true; return api; },
    setY(o) { Object.assign(ax.y, o); dirty = true; return api; },
    /** Highlight one series (dims the others); null clears. Used by equation terms. */
    highlight(id) { focus = id; dirty = true; api.draw(); },
    /** Show a draggable vertical cursor at x; drag calls onCursor(x). */
    cursor(xv) { cursor = xv; dirty = true; return api; },
    invalidate() { dirty = true; },
    draw,
  };

  function autoY() {
    let lo = Infinity, hi = -Infinity;
    for (const it of items.values()) {
      if (it.kind === 'fn' && !it.ignoreScale) {
        for (let i = 0; i <= 200; i++) {
          const xv = ax.x.log ? Math.pow(10, Math.log10(ax.x.min) + i / 200 * (Math.log10(ax.x.max) - Math.log10(ax.x.min))) : ax.x.min + i / 200 * (ax.x.max - ax.x.min);
          const v = it.f(xv);
          if (Number.isFinite(v) && (!ax.y.log || v > 0)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
        }
      } else if (it.kind === 'data' && !it.ignoreScale) {
        for (const [, v] of it.pts) if (Number.isFinite(v) && (!ax.y.log || v > 0)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
      }
    }
    if (!Number.isFinite(lo)) return;
    if (ax.y.log) { ax.y.min = Math.pow(10, Math.floor(Math.log10(lo))); ax.y.max = Math.pow(10, Math.ceil(Math.log10(hi))); }
    else {
      const span = hi - lo || Math.abs(hi) || 1;
      ax.y.min = ax.y.zero !== false && lo >= 0 ? 0 : lo - 0.06 * span;
      ax.y.max = hi + 0.08 * span;
    }
  }

  function draw() {
    if (!dirty) return;
    dirty = false;
    const dpr = devicePixelRatio || 1;
    W = box.clientWidth; H = box.clientHeight;
    if (!W) return;
    if (canvas.width !== Math.round(W * dpr) || canvas.height !== Math.round(H * dpr)) { canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr); }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    if (ax.y.auto) autoY();
    const C = cssVar;
    ctx.font = '10.5px ui-monospace, Menlo, monospace';

    const labelBoxes = [];                              // band and vertical-line labels placed this frame
    // bands (behind everything)
    for (const id of order) {
      const it = items.get(id);
      if (it.kind !== 'band') continue;
      const x0 = Math.max(pad.l, X(it.x0)), x1 = Math.min(W - pad.r, X(it.x1));
      ctx.fillStyle = C(it.color ?? '--grid-strong');
      ctx.globalAlpha = it.alpha ?? 0.5;
      ctx.fillRect(x0, pad.t, x1 - x0, H - pad.t - pad.b);
      ctx.globalAlpha = 1;
      if (it.label) { ctx.fillStyle = C('--muted'); ctx.fillText(it.label, x0 + 6, pad.t + 12); labelBoxes.push({ x: x0 + 6, w: ctx.measureText(it.label).width, y: pad.t + 12 }); }
    }

    // grid + ticks
    ctx.strokeStyle = C('--grid'); ctx.fillStyle = C('--faint'); ctx.lineWidth = 1;
    const xt = ax.x.log ? logTicks(ax.x.min, ax.x.max) : niceTicks(ax.x.min, ax.x.max, Math.max(3, Math.floor(W / 90)));
    const yt = ax.y.log ? logTicks(ax.y.min, ax.y.max) : niceTicks(ax.y.min, ax.y.max, Math.max(3, Math.floor(H / 45)));
    ctx.textAlign = 'center';
    for (const v of xt) { const px = Math.round(X(v)) + 0.5; ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, H - pad.b); ctx.stroke(); ctx.fillText(ax.x.log ? fmtPow10(v) : fmtTick(v), px, H - pad.b + 14); }
    ctx.textAlign = 'right';
    for (const v of yt) { const py = Math.round(Y(v)) + 0.5; ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(W - pad.r, py); ctx.stroke(); ctx.fillText(ax.y.log ? fmtPow10(v) : fmtTick(v), pad.l - 6, py + 3.5); }
    ctx.textAlign = 'start';
    ctx.strokeStyle = C('--line-strong');
    ctx.beginPath(); ctx.moveTo(pad.l + 0.5, pad.t); ctx.lineTo(pad.l + 0.5, H - pad.b + 0.5); ctx.lineTo(W - pad.r, H - pad.b + 0.5); ctx.stroke();
    // axis labels
    ctx.fillStyle = C('--muted'); ctx.font = '11px ui-sans-serif, -apple-system, system-ui, sans-serif';
    const xl = `${ax.x.label ?? ''}${ax.x.unit ? ` (${ax.x.unit})` : ''}`;
    ctx.textAlign = 'right'; ctx.fillText(xl, W - pad.r, H - 4); ctx.textAlign = 'start';
    const yl = `${ax.y.label ?? ''}${ax.y.unit ? ` (${ax.y.unit})` : ''}`;
    ctx.save(); ctx.translate(11, pad.t); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'right'; ctx.fillText(yl, 0, 0); ctx.restore();

    // clip to plot area
    ctx.save();
    ctx.beginPath(); ctx.rect(pad.l, pad.t - 2, W - pad.l - pad.r, H - pad.t - pad.b + 4); ctx.clip();
    const legendItems = [];
    for (const id of order) {
      const it = items.get(id);
      const dim = focus && focus !== id && it.kind !== 'band' ? 0.18 : 1;
      const col = C(it.color ?? '--accent');
      ctx.globalAlpha = dim * (it.alpha ?? 1);
      ctx.strokeStyle = col; ctx.fillStyle = col;
      ctx.lineWidth = (it.width ?? 1.8) * (focus === id ? 1.6 : 1);
      ctx.setLineDash(it.dash ? [5, 4] : it.dot ? [1.5, 3.5] : []);
      if (it.kind === 'fn') {
        ctx.beginPath();
        let pen = false;
        const n = Math.max(2, Math.round((W - pad.l - pad.r) * (it.oversample ?? 1)));
        for (let i = 0; i <= n; i++) {
          const px = pad.l + i / n * (W - pad.l - pad.r);
          const v = it.f(invX(px));
          if (!Number.isFinite(v) || (ax.y.log && v <= 0)) { pen = false; continue; }
          const py = Y(v);
          pen ? ctx.lineTo(px, py) : ctx.moveTo(px, py);
          pen = true;
        }
        ctx.stroke();
        if (it.fill) { ctx.globalAlpha = 0.12 * dim; ctx.lineTo(W - pad.r, Y(Math.max(ax.y.min, 0))); ctx.lineTo(pad.l, Y(Math.max(ax.y.min, 0))); ctx.fill(); }
      } else if (it.kind === 'data') {
        if (it.dots) {
          const shown = ([a, b]) => Number.isFinite(a) && Number.isFinite(b) && (!ax.y.log || b > 0) && (!ax.x.log || a > 0);
          for (const p of it.pts) { if (!shown(p)) continue; ctx.beginPath(); ctx.arc(X(p[0]), Y(p[1]), it.r ?? 2.4, 0, 7); ctx.fill(); }
          if (it.errors) {for (let i = 0; i < it.pts.length; i++) {
            const [a, b] = it.pts[i], e = Math.abs(it.errors[i]);
            if (!shown(it.pts[i]) || !Number.isFinite(e)) continue;
            // On a log axis a bar reaching zero or below has no finite end: run it to the axis floor.
            const lo = ax.y.log && b - e <= ax.y.min ? ax.y.min : b - e;
            ctx.beginPath(); ctx.moveTo(X(a), Y(lo)); ctx.lineTo(X(a), Y(b + e)); ctx.stroke();
          }}
        } else {
          ctx.beginPath();
          it.pts.forEach(([a, b], i) => i ? ctx.lineTo(X(a), Y(b)) : ctx.moveTo(X(a), Y(b)));
          ctx.stroke();
        }
      } else if (it.kind === 'hline' || it.kind === 'vline') {
        ctx.lineWidth = 1; ctx.setLineDash(it.solid ? [] : [3, 4]);
        ctx.beginPath();
        if (it.kind === 'hline') { ctx.moveTo(pad.l, Y(it.y)); ctx.lineTo(W - pad.r, Y(it.y)); }
        else { ctx.moveTo(X(it.x), pad.t); ctx.lineTo(X(it.x), H - pad.b); }
        ctx.stroke();
        if (it.label) {
          ctx.setLineDash([]); ctx.font = '10.5px ui-monospace, Menlo, monospace';
          if (it.kind === 'hline') ctx.fillText(it.label, W - pad.r - ctx.measureText(it.label).width - 4, Y(it.y) - 5 < pad.t + 10 ? Y(it.y) + 13 : Y(it.y) - 5);   // below the line near the top edge
          else {
            // stack labels of nearby vertical lines instead of drawing them on top of each other
            const x0 = X(it.x) + 5, w = ctx.measureText(it.label).width;
            let ly = pad.t + 12;
            while (labelBoxes.some(b => x0 < b.x + b.w + 4 && b.x < x0 + w + 4 && Math.abs(b.y - ly) < 12)) ly += 13;
            labelBoxes.push({ x: x0, w, y: ly });
            ctx.fillText(it.label, x0, ly);
          }
        }
      } else if (it.kind === 'marker') {
        ctx.setLineDash([]);
        ctx.beginPath(); ctx.arc(X(it.x), Y(it.y), it.r ?? 4.5, 0, 7); ctx.fill();
        ctx.strokeStyle = C('--panel'); ctx.lineWidth = 1.5; ctx.stroke();
        if (it.label) { ctx.fillStyle = C('--ink-2'); ctx.font = '10.5px ui-monospace, Menlo, monospace'; ctx.fillText(it.label, X(it.x) + 8, Y(it.y) - 7); }
      }
      ctx.setLineDash([]); ctx.globalAlpha = 1;
      if (it.label && (it.kind === 'fn' || it.kind === 'data') && legend) legendItems.push([col, it.label, it.dash, id]);
    }
    if (cursor != null) {
      ctx.strokeStyle = C('--ink'); ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
      ctx.beginPath(); ctx.moveTo(X(cursor), pad.t); ctx.lineTo(X(cursor), H - pad.b); ctx.stroke(); ctx.setLineDash([]);
    }
    ctx.restore();
    const l = head.querySelector('.l');
    l.replaceChildren(...legendItems.map(([c, t, dash, id]) => el('span', {
      class: `k${focus === id ? ' on' : ''}`, onmouseenter: () => api.highlight(id), onmouseleave: () => api.highlight(null),
    }, el('i', { style: `background:${c};${dash ? 'opacity:.6' : ''}` }), t)));
  }

  function logTicks(min, max) {
    const out = [], e0 = Math.ceil(Math.log10(min) - 1e-9), e1 = Math.floor(Math.log10(max) + 1e-9);
    const stride = Math.max(1, Math.ceil((e1 - e0 + 1) / 6));
    for (let e = e0; e <= e1; e += stride) out.push(Math.pow(10, e));
    return out;
  }

  // cursor dragging
  let dragging = false;
  const setFromEvent = e => {
    if (cursor == null || !onCursor) return;
    const r = canvas.getBoundingClientRect();
    const xv = Math.min(ax.x.max, Math.max(ax.x.min, invX(e.clientX - r.left)));
    cursor = xv; dirty = true; onCursor(xv); draw();
  };
  canvas.addEventListener('pointerdown', e => { if (cursor == null || !onCursor) return; dragging = true; canvas.setPointerCapture(e.pointerId); setFromEvent(e); });
  canvas.addEventListener('pointermove', e => { if (dragging) setFromEvent(e); });
  canvas.addEventListener('pointerup', () => { dragging = false; });
  if (onCursor) canvas.style.cursor = 'ew-resize';

  new ResizeObserver(() => { dirty = true; draw(); }).observe(box);
  document.documentElement.addEventListener('psim-theme', () => { dirty = true; draw(); });
  return api;
}
