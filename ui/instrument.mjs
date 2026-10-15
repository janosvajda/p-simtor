// Reusable instrument UI: DOM helpers, controls, log-scale plots, the simulation
// clock, number formatting, theme handling and data export.
// Domain-agnostic: orbit, field and density views all build on these pieces.

export function el(tag, attrs = {}, ...children) {
  const node = document.createElementNS(
    tag === 'svg' || attrs.xmlns === 'svg' ? 'http://www.w3.org/2000/svg' : 'http://www.w3.org/1999/xhtml', tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v == null || v === false || k === 'xmlns') continue;
    if (k.startsWith('on') && typeof v === 'function') node.addEventListener(k.slice(2), v);
    else if (k === 'class') node.className = v;
    else if (k === 'html') node.innerHTML = v;
    else if (k in node && typeof v !== 'string') node[k] = v;
    else node.setAttribute(k, v === true ? '' : v);
  }
  for (const c of children.flat()) if (c != null && c !== false) node.append(c.nodeType ? c : document.createTextNode(String(c)));
  return node;
}

export const ICONS = {
  play: '<svg viewBox="0 0 16 16"><path d="M4 2.5v11l9-5.5z" fill="currentColor"/></svg>',
  pause: '<svg viewBox="0 0 16 16"><path d="M4 2.5h3v11H4zM9 2.5h3v11H9z" fill="currentColor"/></svg>',
  step: '<svg viewBox="0 0 16 16"><path d="M3 2.5v11l7-5.5zM11 2.5h2v11h-2z" fill="currentColor"/></svg>',
  restart: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.6"><path d="M3 8a5 5 0 1 0 1.6-3.7"/><path d="M3 2.5v3h3"/></svg>',
  fit: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="8" cy="8" r="4.5"/><path d="M8 1v3M8 12v3M1 8h3M12 8h3"/></svg>',
  fullscreen: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M2 6V2h4M14 6V2h-4M2 10v4h4M14 10v4h-4"/></svg>',
  exitFullscreen: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M6 2v4H2M10 2v4h4M6 14v-4H2M10 14v-4h4"/></svg>',
  sun: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="8" cy="8" r="3"/><path d="M8 1v2M8 13v2M1 8h2M13 8h2M3 3l1.4 1.4M11.6 11.6L13 13M3 13l1.4-1.4M11.6 4.4L13 3"/></svg>',
  moon: '<svg viewBox="0 0 16 16"><path d="M13.5 10A6 6 0 0 1 6 2.5a6 6 0 1 0 7.5 7.5z" fill="currentColor"/></svg>',
  info: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="8" cy="8" r="6.5"/><path d="M8 7v4.5M8 4.6v.1"/></svg>',
};

// ---------- formatting ----------

const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };

/** 2.44e-11 -> "2.4×10⁻¹¹". */
export function sci(x, digits = 2) {
  if (x == null || Number.isNaN(x)) return '—';
  if (x === 0) return '0';
  if (!Number.isFinite(x)) return x > 0 ? '∞' : '−∞';
  const e = Math.floor(Math.log10(Math.abs(x)));
  if (e >= -2 && e < 4) return fix(x, Math.max(0, digits - 1 - e));
  const m = x / Math.pow(10, e);
  return `${m.toFixed(digits - 1)}×10${String(e).split('').map(c => SUP[c]).join('')}`;
}

export function fix(x, d = 3) {
  if (x == null || Number.isNaN(x)) return '—';
  return x.toFixed(d).replace('-', '−');
}

/** Classify a relative error against thresholds -> 'good' | 'warn' | 'bad'. */
export function grade(x, good, warn) {
  if (!(x <= warn)) return 'bad';
  return x <= good ? 'good' : 'warn';
}

// ---------- theme ----------

export function initTheme(toggleButton) {
  const root = document.documentElement;
  let stored = null;
  try { stored = localStorage.getItem('psim-theme'); } catch { /* storage unavailable */ }
  const query = new URLSearchParams(location.search).get('theme');
  const sys = matchMedia('(prefers-color-scheme: light)').matches ? 'light' : 'dark';
  const apply = t => {
    root.dataset.theme = t;
    if (toggleButton) toggleButton.innerHTML = t === 'dark' ? ICONS.sun : ICONS.moon;
    if (toggleButton) toggleButton.title = t === 'dark' ? 'Light theme' : 'Dark theme';
    root.dispatchEvent(new CustomEvent('psim-theme', { detail: t }));
  };
  apply(query || stored || root.dataset.theme || sys);
  toggleButton?.addEventListener('click', () => {
    const t = root.dataset.theme === 'dark' ? 'light' : 'dark';
    try { localStorage.setItem('psim-theme', t); } catch { /* ignore */ }
    apply(t);
  });
}

/** A theme token ('--accent') → its value; any other colour ('#ff8800', 'rgb(…)', 'tomato') passes through. */
export function cssVar(name) {
  if (typeof name !== 'string' || !name.startsWith('--')) return name;
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

// ---------- layout pieces ----------

export function section(title, { open = true, aside } = {}) {
  const body = el('div', { class: 'psim-body' });
  const asideEl = aside != null ? el('span', { class: 'aside' }, aside) : null;
  const sec = el('details', { class: 'psim-sec', open }, el('summary', {}, title, asideEl), body);
  sec.body = body;
  sec.aside = asideEl;
  return sec;
}

export function row(label, value = '—', cls = '') {
  const v = el('span', { class: `v ${cls}` }, value);
  const r = el('div', { class: 'psim-row' }, el('span', { class: 'k' }, label), v);
  r.set = (text, c) => { v.textContent = text; v.className = `v ${c ?? ''}`; };
  return r;
}

/** Value for a numeric input: short, exact enough, scientific notation for small or huge numbers. */
export const inputValue = v => {
  const a = Math.abs(v);
  return a !== 0 && (a < 1e-3 || a >= 1e7) ? v.toExponential(4) : String(+Number(v).toPrecision(6));
};

/**
 * Slider + exact numeric entry, linear or logarithmic. onInput receives the number.
 * Typed values are clamped to `limits` (default: the slider's [min, max]) and rounded when
 * `integer` is set, so a typo like 100000000 can never freeze the page. Pass wider `limits`
 * when a quantity may legitimately exceed the slider's range.
 */
export function slider({ label, min, max, step = 'any', value, log = false, unit = '', digits = 3, integer = false, limits = [min, max], onInput }) {
  const toS = v => log ? Math.log10(v) : v;
  const fromS = s => log ? Math.pow(10, s) : s;
  const [lo, hi] = limits;
  const tidy = v => {
    const c = Math.min(hi, Math.max(lo, v));
    return integer ? Math.round(c) : c;
  };
  const range = el('input', { type: 'range', min: toS(min), max: toS(max), step: log ? 'any' : step, value: toS(value), 'aria-label': label });
  const num = el('input', { type: 'number', step: integer ? 1 : 'any', min: lo, max: hi, value: inputValue(value), 'aria-label': `${label} value` });
  const unitEl = unit ? el('span', { style: 'color:var(--muted);font:11px var(--mono)' }, unit) : null;
  const note = el('span', { class: 'psim-limit', role: 'status' });
  const field = el('label', { class: 'psim-field' },
    el('div', { class: 'top' }, el('span', {}, label), el('div', { style: 'display:flex;align-items:baseline;gap:4px' }, num, unitEl)),
    range, note);
  range.addEventListener('input', () => {
    const v = tidy(+(+fromS(+range.value)).toPrecision(digits + 1));
    num.value = inputValue(v);
    note.textContent = '';
    onInput?.(v);
  });
  num.addEventListener('change', () => {
    const raw = +num.value;
    if (!Number.isFinite(raw)) { num.value = inputValue(tidy(fromS(+range.value))); return; }
    const v = tidy(raw);
    note.textContent = v !== raw ? `limited to ${inputValue(v)} (allowed ${inputValue(lo)}–${inputValue(hi)})` : '';
    num.value = inputValue(v);
    range.value = toS(v);
    onInput?.(v);
  });
  field.setValue = v => { range.value = toS(v); num.value = inputValue(v); };
  return field;
}

export function toggle(label, checked, onChange) {
  const input = el('input', { type: 'checkbox', checked });
  input.addEventListener('change', () => onChange(input.checked));
  const t = el('label', { class: 'psim-toggle' }, el('span', {}, label), input);
  t.input = input;
  return t;
}

export function select(options, value, onChange) {
  const s = el('select', { class: 'psim-select' },
    options.map(o => el('option', { value: o.value, selected: o.value === value }, o.label)));
  s.addEventListener('change', () => onChange(s.value));
  return s;
}

export function button(label, onClick, cls = '') {
  return el('button', { class: `psim-btn ${cls}`, type: 'button', onclick: onClick }, label);
}

export function toast(stage, text) {
  let t = stage.querySelector('.psim-toast');
  if (!t) { t = el('div', { class: 'psim-toast', role: 'status' }); stage.append(t); }
  t.textContent = text;
  t.classList.add('on');
  clearTimeout(t._h);
  t._h = setTimeout(() => t.classList.remove('on'), 1800);
}

// ---------- transport ----------

export function transport({ onToggle, onStep, onRestart, onFit, onFullscreen, speed, speedLabel, onSpeed, range: [lo, hi] = [speed / 100, speed * 100] }) {
  const play = el('button', { class: 'primary', title: 'Play / pause (Space)', 'aria-label': 'Play or pause', html: ICONS.play, onclick: onToggle });
  const out = el('output', {}, speedLabel(speed));
  const range = el('input', { type: 'range', min: Math.log10(lo), max: Math.log10(hi), step: 'any', value: Math.log10(speed), 'aria-label': 'Playback speed' });
  range.addEventListener('input', () => { const s = Math.pow(10, +range.value); out.textContent = speedLabel(s); onSpeed(s); });
  const bar = el('div', { class: 'psim-transport' },
    el('button', { title: 'Restart from initial conditions (R)', 'aria-label': 'Restart', html: ICONS.restart, onclick: onRestart }),
    play,
    el('button', { title: 'Advance one integrator step (.)', 'aria-label': 'Step', html: ICONS.step, onclick: onStep }),
    el('span', { class: 'sep' }),
    el('label', { class: 'speed' }, 'speed', range, out),
    el('span', { class: 'sep' }),
    el('button', { title: 'Fit the view', 'aria-label': 'Fit view', html: ICONS.fit, onclick: onFit }),
    onFullscreen ? el('button', { class: 'fs', title: 'Full screen (F)', 'aria-label': 'Full screen', html: ICONS.fullscreen, onclick: onFullscreen }) : null);
  bar.setPlaying = p => { play.innerHTML = p ? ICONS.pause : ICONS.play; };
  bar.setFullscreen = on => {
    const b = bar.querySelector('button.fs');
    if (b) { b.innerHTML = on ? ICONS.exitFullscreen : ICONS.fullscreen; b.title = on ? 'Leave full screen (F or Esc)' : 'Full screen (F)'; }
  };
  return bar;
}

/**
 * Simulation clock. Converts wall-clock frames into a target *physical* time.
 * The integrator then advances to that target with its own step sequence, so
 * the physics never depends on frame rate. Long frames (background tabs) are capped.
 */
export function createClock({ speed = 1 } = {}) {
  return {
    speed, playing: false, last: null,
    tick(now, t) {
      if (this.last == null || !this.playing) { this.last = now; return t; }
      const dt = Math.min(0.1, (now - this.last) / 1000);
      this.last = now;
      return t + dt * this.speed;
    },
  };
}

// ---------- plots ----------

/**
 * Lightweight log-y line plot for diagnostics (energy error, separations).
 * series: [{ color: cssVarName, data: [[x, y], ...] }]
 */
export function logPlot({ caption, yMin = 1e-16, yMax = 1, height = 112 }) {
  const canvas = el('canvas');
  const box = el('div', { class: 'psim-plot', style: `height:${height}px` }, canvas, el('span', { class: 'cap' }, caption));
  const series = [];
  let xMax = 1, dirty = true, refLines = [];
  const ctx = canvas.getContext('2d');
  const plot = {
    el: box,
    addSeries(color, opts = {}) { const s = { color, data: [], ...opts }; series.push(s); return s; },
    push(s, x, y) { s.data.push([x, y]); if (s.data.length > 4000) decimate(s); xMax = Math.max(xMax, x); dirty = true; },
    clear() { for (const s of series) s.data.length = 0; xMax = 1; dirty = true; },
    setXMax(x) { xMax = x; dirty = true; },
    setRefLines(lines) { refLines = lines; dirty = true; },
    invalidate() { dirty = true; },
    draw() {
      if (!dirty) return;
      dirty = false;
      const dpr = devicePixelRatio || 1;
      const w = box.clientWidth, h = box.clientHeight;
      if (!w) return;
      if (canvas.width !== Math.round(w * dpr)) { canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); }
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, w, h);
      const l0 = Math.log10(yMin), l1 = Math.log10(yMax);
      const pad = { l: 34, r: 8, t: 20, b: 14 };
      const X = x => pad.l + (x / xMax) * (w - pad.l - pad.r);
      const Y = y => {
        const ly = Math.min(l1, Math.max(l0, Math.log10(Math.max(y, 1e-300))));
        return pad.t + (1 - (ly - l0) / (l1 - l0)) * (h - pad.t - pad.b);
      };
      ctx.font = '9.5px ui-monospace, Menlo, monospace';
      ctx.fillStyle = cssVar('--faint');
      ctx.strokeStyle = cssVar('--grid-strong');
      ctx.lineWidth = 1;
      const stepDec = Math.max(1, Math.round((l1 - l0) / 4));
      for (let e = Math.ceil(l0); e <= l1; e += stepDec) {
        const y = Math.round(Y(Math.pow(10, e))) + 0.5;
        ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();
        ctx.fillText(`1e${e}`, 2, y + 3);
      }
      for (const r of refLines) {
        ctx.setLineDash([3, 3]); ctx.strokeStyle = cssVar(r.color);
        const y = Y(r.y); ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();
        ctx.setLineDash([]);
      }
      for (const s of series) {
        if (!s.data.length) continue;
        ctx.strokeStyle = cssVar(s.color);
        ctx.lineWidth = s.width ?? 1.4;
        ctx.globalAlpha = s.alpha ?? 1;
        ctx.beginPath();
        s.data.forEach(([x, y], i) => i ? ctx.lineTo(X(x), Y(y)) : ctx.moveTo(X(x), Y(y)));
        ctx.stroke();
        ctx.globalAlpha = 1;
      }
    },
  };
  new ResizeObserver(() => { dirty = true; plot.draw(); }).observe(box);
  document.documentElement.addEventListener('psim-theme', () => { dirty = true; plot.draw(); });
  return plot;
}

function decimate(s) {
  // Keep the envelope: max of each pair (errors are what matter, so keep the worse one).
  const out = [];
  for (let i = 0; i < s.data.length; i += 2) {
    const a = s.data[i], b = s.data[i + 1];
    out.push(!b || a[1] >= b[1] ? a : b);
  }
  s.data = out;
}

// ---------- export ----------

export function download(filename, text, mime = 'text/plain') {
  const url = URL.createObjectURL(new Blob([text], { type: mime }));
  const a = el('a', { href: url, download: filename });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function readJsonScript(id) {
  const node = document.getElementById(id);
  if (!node) return null;
  const text = node.textContent.trim();
  return text ? JSON.parse(text) : null;
}

/**
 * Status beacon read by scripts/psim.mjs check (headless Chrome).
 * Pages set fields on it; it is mirrored into <html data-psim-*> attributes.
 */
export function beacon() {
  const root = document.documentElement;
  const errors = [];
  addEventListener('error', e => { errors.push(String(e.message)); root.dataset.psimErrors = JSON.stringify(errors); });
  addEventListener('unhandledrejection', e => { errors.push(String(e.reason)); root.dataset.psimErrors = JSON.stringify(errors); });
  root.dataset.psimErrors = '[]';
  return {
    set(key, value) { root.dataset[`psim${key[0].toUpperCase()}${key.slice(1)}`] = typeof value === 'string' ? value : JSON.stringify(value); },
  };
}
