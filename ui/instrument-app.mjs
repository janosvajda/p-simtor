// Generic instrument shell for composed simulations (any domain).
//
// The page supplies a physics module that follows the simulation contract
// (see references/compose.md) and a draw function; this shell provides the
// header, trust tier, theme, transport, frame-rate-independent clock,
// parameter controls, live diagnostics, verification report, model panel,
// export and the window.psim hook used by automated checks.
//
// Simulation contract (returned by create(params)):
//   sim.t                    physical time of the latest computed state
//   sim.step()               advance by one integrator step (or one analytic tick)
//   sim.sample?(t)           state at t between steps, for drawing (default: current state)
//   sim.diagnostics()        [{ key, label, value, kind: 'error'|'value', unit?, good?, warn? }]
//   sim.state()              JSON-serialisable state (export, tests)
//   sim.record?()            { header: [...], rows: [[...], ...] } for CSV export
//   sim.timeDependent        false for static analytic visualisations (no transport)

import {
  el, sci, fix, grade, initTheme, section, row, slider, toggle, select, button, toast, ICONS,
  transport, createClock, logPlot, download, readJsonScript, beacon, cssVar,
} from './instrument.mjs';
import { createPlot } from './plot.mjs';
import { createGLLayer } from './gl.mjs';
import { createEquation, createDerivation, M } from './equation.mjs';
import { validateModuleReport } from '../core/verification.mjs';

const TIERS = {
  validated: 'Validated numerical',
  approximate: 'Approximate numerical',
  analytical: 'Analytical',
  conceptual: 'Conceptual illustration',
};

export function createInstrument(spec) {
  const {
    title, subtitle = '', tier, tierReason, create, draw, params = [], model = {},
    timeUnit = 'T', timeLabel = 'simulated time', speed: speed0 = 1, gl: wantGL = false, drawGL, controls, onPointer, layout = 'stage', bench,
    presets = [], texts = {}, hudError = true, formatTime, timeDigits: digitsGiven, frameBudget = 9, speedRange, viewControls, testSteps = 300,
    glResolution = 'auto', timeline: timelineSpec,
  } = spec;
  // counted progress (sweeps, steps, generations, events) shows whole numbers unless told otherwise
  const timeDigits = digitsGiven ?? (/sweep|step|generation|event|iteration|cycle|round/i.test(timeUnit) ? 0 : 3);
  if (!TIERS[tier]) throw new Error(`tier must be one of ${Object.keys(TIERS).join(', ')}`);
  const report = readJsonScript('psim-verification');
  if (report) validateModuleReport(report);
  const verified = report?.checks.every(c => c.pass === true) === true;
  const trustLabel = verified ? TIERS[tier] : 'Verification needed';
  const status = beacon();
  const query = new URLSearchParams(location.search);

  document.title = title;
  document.body.innerHTML = `
  <div class="psim-app ${layout === 'workbench' ? 'workbench' : ''}">
    <header class="psim-head">
      <div class="psim-mark" aria-hidden="true"></div>
      <div class="psim-title"><h1></h1><p></p></div>
      <div class="spacer"></div>
      <button class="psim-tier" type="button"><i></i><span></span></button>
      <button class="psim-iconbtn" type="button" aria-label="Toggle theme"></button>
    </header>
    <main class="psim-stage"><canvas></canvas>
      <div class="psim-hud"><span class="dim"></span><b>0</b><span class="h1"></span><span class="dim h2"></span><span class="flag"></span></div>
      <div class="psim-timeline" hidden><div class="ticks"></div><div class="bar"><i class="fill"></i></div><div class="ends"><span></span><span></span></div></div>
      <div class="psim-legend"></div>
    </main>
    <aside class="psim-side"></aside>
  </div>`;
  const $ = s => document.querySelector(s);
  $('.psim-title h1').textContent = title;
  $('.psim-title p').textContent = subtitle;
  const tierBtn = $('.psim-tier');
  tierBtn.dataset.tier = verified ? tier : 'unverified';
  tierBtn.querySelector('span').textContent = trustLabel;
  initTheme($('.psim-iconbtn'));
  const stage = $('.psim-stage'), side = $('.psim-side');
  // Workbench layout: the live canvas on top, then term-split equations and figures that
  // update with the model (spec.bench builds them; see references/compose.md).
  let benchParts = null;
  if (layout === 'workbench') {
    const wrap = el('div', { class: 'psim-canvaswrap' }, stage.querySelector('canvas'), stage.querySelector('.psim-hud'), stage.querySelector('.psim-legend'));
    const eqs = el('div', { class: 'psim-eqs' }), figs = el('div', { class: 'psim-figs' });
    stage.prepend(el('div', { class: 'psim-bench' }, wrap, eqs, figs));
    benchParts = { wrap, eqs, figs };
  }
  const canvas = stage.querySelector('canvas');
  const ctx = canvas.getContext('2d');
  // optional GPU layer under the 2D canvas (ui/gl.mjs): drawGL renders the scene, draw adds labels
  const glLayer = wantGL ? createGLLayer(canvas) : null;
  if (glLayer) $('.psim-app').classList.add('has-gl');      // HUD and legend get a backing over bright GPU scenes
  $('.psim-hud .dim').textContent = timeLabel;
  const hud = { t: $('.psim-hud b'), a: $('.psim-hud .h1'), b: $('.psim-hud .h2'), flag: $('.psim-hud .flag') };
  // the legend beside the view; view.setLegend(items) replaces it when the representation changes
  const setLegend = items => $('.psim-legend').replaceChildren(...(items ?? []).map(t => el('span', {}, t)));
  setLegend(model.legend);

  // ------------------------------------------------------------ state
  const values = Object.fromEntries(params.map(p => [p.key, p.value]));
  let sim, tDisplay = 0, limited = false;
  const clock = createClock({ speed: speed0 });
  const view = {
    W: 0, H: 0, dpr: 1, color: cssVar, opts: {}, pointer: null, gl: glLayer?.gl ?? null, glError: glLayer && !glLayer.ok ? glLayer.reason : '',
    /** Current playback speed (simulated time per second) and state, e.g. to disclose time compression. */
    get speed() { return clock.speed; }, get playing() { return clock.playing; },
    /** GL layer resolution as a fraction of the device pixel ratio (adaptive unless glResolution is a number). */
    get glScale() { return glLayer?.scale ?? 1; },
    /** 'dark' or 'light': the theme now shown (e.g. additive vs alpha blending, mesh colours). */
    get theme() { return document.documentElement.dataset.theme === 'light' ? 'light' : 'dark'; },
    setLegend,
    /** The embedded build-time report (verify()'s return value, including any extra data fields), or null. */
    report,
  };

  let benchApi = null;
  // ------------------------------------------------------------ timeline: the run's span and its key times
  // spec.timeline = { duration, markers: [{ t, label, kind: 'computed' | 'measured' | 'requested' }], note }
  // (each may also be a function of the sim). The run stops exactly at `duration`; markers are drawn
  // on a bar under the time readout; clicking the bar jumps to that time (recomputing if earlier).
  let tl = null;
  const tlEl = $('.psim-timeline');
  if (benchParts) benchParts.wrap.append(tlEl);              // inside the scene card in the workbench layout
  const secTimeline = timelineSpec ? section(texts.timeline ?? 'Timeline', { aside: 'key times' }) : null;
  const shortLabel = m => m.short ?? m.label.replace(/\s*\(.*$/, '').slice(0, 28);
  function resolveTimeline() {
    if (!timelineSpec) { tl = null; tlEl.hidden = true; return; }
    const val = v => (typeof v === 'function' ? v(sim, values) : v);
    // the run spans [start, duration] on the page's own clock: start = the sim's initial time (0 usually;
    // e.g. −1 for a log₁₀ clock); duration is the END time on that clock
    const start = Number.isFinite(+val(timelineSpec.start)) ? +val(timelineSpec.start) : (sim.t ?? 0);
    const duration = +val(timelineSpec.duration);
    if (!(duration > start)) throw new Error(`timeline.duration (the end time) must be later than the start ${start}`);
    const span = duration - start, tol = 1e-9 * Math.max(Math.abs(duration), span);
    const markers = (val(timelineSpec.markers) ?? []).filter(m => Number.isFinite(m.t));
    for (const m of markers) if (m.t < start - tol || m.t > duration + tol) throw new Error(`timeline marker "${m.label}" at ${m.t} lies outside the run ${start}…${duration}`);
    tl = { start, duration, markers, note: val(timelineSpec.note) ?? '' };
    tlEl.hidden = false;
    const fmt = t => (formatTime ? formatTime(t) : `${+t.toPrecision(4)} ${timeUnit}`);
    // the bar: one tick per key time, short labels above (staggered when close), full text as a tooltip
    const bar = tlEl.querySelector('.bar'), ticks = tlEl.querySelector('.ticks');
    bar.querySelectorAll('.m').forEach(e => e.remove());
    let lastX = -1, level = 0;
    ticks.replaceChildren(...markers.slice().sort((a, b) => a.t - b.t).map(m => {
      const x = 100 * (m.t - start) / span;
      level = x - lastX < 22 ? 1 - level : 0; lastX = x;
      bar.append(el('span', { class: `m ${m.kind ?? 'computed'}`, style: `left:${x.toFixed(3)}%`, title: `${m.label}: ${fmt(m.t)}` }));
      const sl = shortLabel(m), withTime = /\d/.test(sl) ? sl : `${sl} · ${fmt(m.t)}`;   // no repeated time
      // near the ends, anchor the label at its tick instead of centring it (never clipped by the bar's edges)
      const anchor = x < 15 ? 'translateX(0)' : x > 85 ? 'translateX(-100%)' : 'translateX(-50%)';
      return el('span', { class: `lab ${m.kind ?? 'computed'} r${level}`, style: `left:${x.toFixed(3)}%;transform:${anchor}`, title: `${m.label}: ${fmt(m.t)}` }, withTime);
    }));
    const [a, b] = tlEl.querySelectorAll('.ends span');
    a.textContent = fmt(start); b.textContent = `${fmt(duration)} (end)`;
    // the side panel lists every key time in full, with its kind, and the note
    if (secTimeline) {
      const kinds = { computed: 'computed in this model', measured: 'measured (literature)', requested: 'named in the request' };
      secTimeline.body.replaceChildren(
        row('run', `${fmt(start)} – ${fmt(duration)}`),
        ...markers.map(m => el('div', { class: 'psim-tlrow' }, el('i', { class: m.kind ?? 'computed' }), el('span', {}, m.label), el('code', {}, fmt(m.t)), el('small', {}, kinds[m.kind ?? 'computed']))),
        tl.note ? el('p', { class: 'psim-hint' }, tl.note) : null);
    }
  }
  tlEl.querySelector('.bar').addEventListener('click', e => {
    if (!tl) return;
    const r = e.currentTarget.getBoundingClientRect(), T = tl.start + Math.max(0, Math.min(1, (e.clientX - r.left) / r.width)) * (tl.duration - tl.start);
    if (T < tDisplay) rebuild();
    window.psim.advanceTo(T);
  });

  function rebuild() {
    sim = create({ ...values });
    tDisplay = sim.t ?? 0;
    resolveTimeline();
    diagPlot.clear();
    renderDiag(true);
    benchApi?.update?.(sim, tDisplay, { rebuilt: true });
    status.set('params', values);
  }

  // ------------------------------------------------------------ full screen
  // The scene (canvas, HUD, legend) and the playback bar fill the screen. Falls back to filling the
  // window with CSS where the Fullscreen API is missing or refused (iPhone Safari, sandboxed frames).
  const scene = benchParts?.wrap ?? stage;
  let pseudoFs = false, fsBtn = null;
  const fsElement = () => document.fullscreenElement ?? document.webkitFullscreenElement ?? null;
  const isFullscreen = () => pseudoFs || fsElement() === scene;
  function afterFullscreen() {
    const on = isFullscreen();
    scene.classList.toggle('psim-fs', on);
    if (tp) { if (on && tp.parentNode !== scene) scene.append(tp); else if (!on && tp.parentNode !== stage) stage.append(tp); tp.setFullscreen(on); }
    if (fsBtn) fsBtn.innerHTML = on ? ICONS.exitFullscreen : ICONS.fullscreen;
    requestAnimationFrame(() => spec.onFit?.(view));
  }
  function setPseudo(on) { pseudoFs = on; scene.classList.toggle('psim-pseudo-fs', on); afterFullscreen(); }
  function toggleFullscreen() {
    if (isFullscreen()) {
      if (pseudoFs) setPseudo(false);
      else (document.exitFullscreen ?? document.webkitExitFullscreen)?.call(document);
      return;
    }
    const request = scene.requestFullscreen ?? scene.webkitRequestFullscreen;
    if (!request) { setPseudo(true); return; }
    try { Promise.resolve(request.call(scene)).catch(() => setPseudo(true)); } catch { setPseudo(true); }
  }
  addEventListener('fullscreenchange', afterFullscreen);
  addEventListener('webkitfullscreenchange', afterFullscreen);

  // ------------------------------------------------------------ transport
  const timeDependent = spec.timeDependent !== false;
  let tp = null;
  if (timeDependent) {
    tp = transport({
      speed: speed0, speedLabel: s => `${sci(s, 2)} ${timeUnit}/s`, range: speedRange,
      onToggle: () => togglePlay(),
      onStep: () => { setPlaying(false); if (!isDone()) { sim.step(); tDisplay = sim.t; } },
      onRestart: () => rebuild(),
      onFit: () => spec.onFit?.(view),
      onFullscreen: () => toggleFullscreen(),
      onSpeed: s => { clock.speed = s; },
    });
    stage.append(tp);
  } else {
    fsBtn = el('button', { class: 'psim-fsbtn', title: 'Full screen (F)', 'aria-label': 'Full screen', html: ICONS.fullscreen, onclick: () => toggleFullscreen() });
    scene.append(fsBtn);
    stage.querySelector('.psim-hud').querySelector('.dim').textContent = '';
  }
  function setPlaying(p) { clock.playing = p; tp?.setPlaying(p); }
  // A model that ends (a merger, everything decayed, the last particle left) sets sim.done (or
  // sim.ended): the clock stops there instead of counting time in which nothing is computed,
  // and play starts the run again from the beginning.
  // the end shown on the clock: exactly the declared duration when the steps' float sum lands within round-off of it
  const endTime = () => (tl && Math.abs(sim.t - tl.duration) <= 1e-9 * Math.max(Math.abs(tl.duration), tl.duration - tl.start) ? tl.duration : tl ? Math.min(sim.t, tl.duration) : sim.t);
  const isDone = () => sim?.done === true || sim?.ended === true || (tl != null && (sim?.t ?? 0) >= tl.duration - 1e-12 * Math.max(Math.abs(tl.duration), tl.duration - tl.start));
  function togglePlay() {
    if (!clock.playing && isDone()) { rebuild(); setPlaying(true); return; }
    setPlaying(!clock.playing);
  }


  // ------------------------------------------------------------ side panel
  const secParams = section(texts.parameters ?? 'Parameters');
  const paramInputs = new Map();
  // set one or several parameters as if typed: every control shows the new value, one rebuild
  function setParams(next, { play } = {}) {
    for (const [k, v] of Object.entries(next)) { values[k] = v; paramInputs.get(k)?.setValue?.(v); }
    rebuild();
    if (play != null && timeDependent) setPlaying(play);
  }
  const setParam = (k, v) => setParams({ [k]: v });
  if (presets.length) {
    secParams.body.append(el('div', { class: 'psim-presets' },
      el('span', { class: 'psim-hint' }, texts.presets ?? 'Try this'),
      ...presets.map(pr => button(pr.label, () => setParams(pr.values ?? {}, { play: pr.play ?? true })))));
  }
  if (params.length) {
    const note = texts.paramsNote ?? 'Changing a parameter recomputes the model from its initial state.';
    if (note) secParams.body.append(el('p', { class: 'psim-hint' }, note));
    for (const p of params) {
      if (p.options) {
        const sel = select(p.options, p.value, v => { values[p.key] = isNaN(+v) ? v : +v; rebuild(); });
        sel.setValue = v => { sel.value = String(v); };
        paramInputs.set(p.key, sel);
        secParams.body.append(el('div', { class: 'psim-field' }, el('div', { class: 'top' }, el('span', {}, p.label)), sel));
      } else {
        const sl = slider({
          label: p.label, min: p.min, max: p.max, value: p.value, log: p.log, unit: p.unit ?? '', digits: p.digits ?? 3,
          integer: p.integer, limits: p.limits,
          onInput: v => { values[p.key] = v; rebuild(); },
        });
        paramInputs.set(p.key, sl);
        secParams.body.append(sl);
      }
      if (p.hint) secParams.body.append(el('p', { class: 'psim-hint', style: 'margin-top:-4px;font-size:11.5px' }, p.hint));
    }
  }
  if (controls) controls(secParams.body, { values, rebuild, setParam, setParams, play: () => setPlaying(true), pause: () => setPlaying(false), getSim: () => sim });

  // presentation-only controls (projection, colour map, slice): never rebuild the model
  const secView = viewControls ? section(texts.view ?? 'View') : null;
  if (secView) secView.body.append(el('p', { class: 'psim-hint' }, texts.viewNote ?? 'These change only how the results are shown, not the calculation.'));

  const secDiag = section(texts.diagnostics ?? 'Diagnostics', { aside: 'live' });
  const diagPlot = logPlot({ caption: '', yMin: 1e-16, yMax: 1 });
  const diagSeries = diagPlot.addSeries('--accent');
  const diagRows = new Map();
  secDiag.body.append(diagPlot.el);
  let nextPlot = 0;
  function renderDiag(reset = false) {
    const d = sim.diagnostics();
    if (reset) {
      diagRows.forEach(r => r.remove());
      diagRows.clear();
      nextPlot = 0;
      const firstErr = d.find(x => x.kind === 'error');
      diagPlot.el.style.display = firstErr ? '' : 'none';
      diagPlot.el.querySelector('.cap').textContent = firstErr ? `${firstErr.label}   vs   t` : '';
      if (spec.duration) diagPlot.setXMax(spec.duration);
    }
    for (const x of d) {
      let r = diagRows.get(x.key);
      if (!r) { r = row(x.label); diagRows.set(x.key, r); secDiag.body.append(r); }
      // pending: a quantity that is not measured yet (e.g. before enough Monte Carlo samples)
      const text = x.pending ? (x.text ?? texts.pending ?? 'measuring…') : x.text ?? (x.kind === 'error' ? sci(x.value) : `${sci(x.value, 4)}${x.unit ? ` ${x.unit}` : ''}`);
      r.set(text, !x.pending && x.kind === 'error' && x.good != null ? grade(x.value, x.good, x.warn ?? x.good * 1000) : '');
    }
    const firstErr = d.find(x => x.kind === 'error' && !x.pending);
    if (firstErr && timeDependent && sim.t >= nextPlot) {
      diagPlot.push(diagSeries, sim.t, Math.max(firstErr.value, 1e-17));
      nextPlot = sim.t + (spec.duration ?? 10) / 800;
    }
    diagPlot.draw();
    return d;
  }

  const secVerify = section('Verification', { aside: report ? 'report embedded' : 'no report' });
  if (report?.checks) {
    secVerify.body.append(el('p', { class: 'psim-hint' }, report.summary ?? `Measured when this page was built (${(report.generated ?? '').slice(0, 10)}):`));
    for (const c of report.checks) {
      secVerify.body.append(el('div', { class: 'psim-check' },
        el('i', { class: c.pass !== true ? 'bad' : c.state === 'warn' ? 'warn' : '' }),
        el('span', {}, c.label), el('code', {}, c.value)));
    }
    if (report.notes) for (const n of report.notes) secVerify.body.append(el('p', { class: 'psim-hint' }, n));
  } else {
    secVerify.body.append(el('p', { class: 'psim-hint' }, 'No build-time verification report is embedded in this page.'));
  }

  const secModel = section('Model, assumptions & limits', { open: false });
  secModel.body.append(
    el('div', { class: 'psim-prose' }, el('p', {}, el('b', {}, `${trustLabel}. `), tierReason ?? '')),
    ...(model.equations ?? []).flatMap(eq => [
      eq.caption ? el('p', { class: 'psim-hint' }, eq.caption) : null,
      el('div', { class: 'psim-math', html: eq.mathml }),
    ]),
    model.assumptions ? el('div', { class: 'psim-prose' }, el('p', {}, el('b', {}, 'Assumptions and limits')), el('ul', {}, model.assumptions.map(a => el('li', {}, a)))) : null,
    model.conventions ? el('div', { class: 'psim-prose' }, el('p', {}, el('b', {}, 'Visual conventions')), el('ul', {}, model.conventions.map(a => el('li', {}, a)))) : null,
    model.references ? el('div', { class: 'psim-prose' }, el('p', {}, el('b', {}, 'References')), el('ul', {}, model.references.map(a => el('li', {}, el('cite', {}, a))))) : null,
  );
  tierBtn.addEventListener('click', () => { secModel.open = true; secModel.scrollIntoView({ behavior: 'smooth' }); });

  const secExport = section('Export', { open: false });
  const slug = title.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
  secExport.body.append(el('div', { class: 'psim-btns' },
    button('State JSON', () => download(`${slug}-state.json`, JSON.stringify({
      generator: 'p-simtor', title, tier: tierBtn.dataset.tier, params: values, t: sim.t, state: sim.state(), diagnostics: sim.diagnostics(), verification: report,
    }, null, 2), 'application/json')),
    button('Data CSV', () => {
      const r = sim.record?.();
      if (!r) { toast(stage, 'This model records no tabular data'); return; }
      download(`${slug}-data.csv`, [`# p-simtor | ${title}`, r.header.join(','), ...r.rows.map(x => x.join(','))].join('\n'), 'text/csv');
    }),
    button('PNG image', () => {
      // composite the GPU layer (if any) under the 2D overlay
      const out = document.createElement('canvas');
      out.width = canvas.width; out.height = canvas.height;
      const c = out.getContext('2d');
      if (glLayer?.ok) c.drawImage(glLayer.canvas, 0, 0);
      c.drawImage(canvas, 0, 0);
      out.toBlob(b => el('a', { href: URL.createObjectURL(b), download: `${slug}.png` }).click());
    }),
  ));

  side.append(...[secParams, secView, secTimeline, secDiag, secVerify, secModel, secExport].filter(s => s && (s !== secParams || params.length || controls)));
  if (viewControls) viewControls(secView.body, { el, slider, select, toggle, button, view, getSim: () => sim });
  if (spec.extraSections) spec.extraSections(side, { section, el, row, slider, toggle, button, values, rebuild, view });

  // ------------------------------------------------------------ canvas & loop
  function resize() {
    view.dpr = devicePixelRatio || 1;
    view.W = canvas.clientWidth; view.H = canvas.clientHeight;
    canvas.width = Math.round(view.W * view.dpr); canvas.height = Math.round(view.H * view.dpr);
  }
  new ResizeObserver(resize).observe(canvas);
  resize();
  if (onPointer) onPointer(canvas, view, () => sim);

  // Adaptive GL resolution: when frames take longer than the display refresh while the JavaScript
  // part is small, the GPU is the bottleneck (large or fullscreen windows, ray marching): render the
  // GL layer at a lower pixel ratio, and raise it again when there is headroom. Display only.
  if (glLayer && typeof glResolution === 'number') glLayer.scale = Math.min(1, Math.max(0.25, glResolution));
  const perf = { last: null, interval: null, js: 0, refresh: Infinity, lastAdjust: 0, calm: 0 };
  function adaptResolution(now, jsMs) {
    if (!glLayer?.ok || glResolution !== 'auto') return;
    if (perf.last != null) {
      const dt = now - perf.last;
      if (dt > 0 && dt < 250) {
        perf.refresh = Math.min(perf.refresh * 1.0005, Math.max(dt, 4));
        perf.interval = perf.interval == null ? dt : 0.9 * perf.interval + 0.1 * dt;
        perf.js = 0.9 * perf.js + 0.1 * jsMs;
        const floor = Math.max(0.35, 0.75 / (view.dpr || 1)), slow = perf.interval > 1.3 * perf.refresh;
        if (slow && perf.interval - perf.js > 0.5 * perf.refresh && now - perf.lastAdjust > 500 && glLayer.scale > floor) {
          glLayer.scale = Math.max(floor, glLayer.scale * 0.8); perf.lastAdjust = now; perf.calm = 0;
        } else if (!slow && perf.interval < 1.1 * perf.refresh) {
          if (++perf.calm > 120 && glLayer.scale < 1 && now - perf.lastAdjust > 2000) { glLayer.scale = Math.min(1, glLayer.scale * 1.15); perf.lastAdjust = now; perf.calm = 0; }
        } else perf.calm = 0;
      }
    }
    perf.last = now;
  }

  let frames = 0, lastUi = 0;
  const frameWaiters = [];   // resolved right after a frame has been drawn (GPU self-tests read pixels then)
  function frame(now) {
    const frameStart = performance.now();
    if (timeDependent) {
      const target = tl ? Math.min(tl.duration, clock.tick(now, tDisplay)) : clock.tick(now, tDisplay);
      if (target > tDisplay) {
        // frameBudget ms of computing per frame, checked after every step, so a slow step
        // (a Monte Carlo sweep) runs at most once past the budget; fast steps run in thousands
        const deadline = performance.now() + frameBudget;
        while (sim.t < target && !isDone()) {
          sim.step();
          if (performance.now() > deadline) break;
        }
        limited = sim.t < target && !isDone();
        tDisplay = Math.min(target, sim.t);
        if (isDone()) { tDisplay = endTime(); setPlaying(false); }
      }
    }
    if (glLayer?.ok) { glLayer.sync(view.dpr); drawGL?.(glLayer.gl, sim, view, tDisplay); }
    ctx.setTransform(view.dpr, 0, 0, view.dpr, 0, 0);
    ctx.clearRect(0, 0, view.W, view.H);
    if (glLayer && !glLayer.ok) { ctx.fillStyle = cssVar('--muted'); ctx.font = '13px system-ui, sans-serif'; ctx.fillText(`${glLayer.reason}: the 3D view cannot be shown.`, 16, view.H - 20); }
    draw(ctx, sim, view, tDisplay);
    if (now - lastUi > 100) {
      lastUi = now;
      const d = renderDiag();
      hud.t.textContent = timeDependent ? (formatTime ? formatTime(tDisplay) : `${fix(tDisplay, timeDigits)} ${timeUnit}`) : '';
      if (tl) tlEl.querySelector('.fill').style.width = `${(100 * Math.max(0, Math.min(1, (tDisplay - tl.start) / (tl.duration - tl.start)))).toFixed(2)}%`;
      const firstErr = d.find(x => x.kind === 'error' && !x.pending);
      hud.a.textContent = firstErr && hudError ? `${firstErr.label} ${sci(firstErr.value)}` : '';
      hud.b.textContent = spec.hudLine?.(sim) ?? '';
      hud.flag.textContent = isDone() ? (texts.finished ?? 'finished · play restarts')
        : limited ? 'compute-limited: playback slowed to keep the calculation exact' : '';
      benchApi?.update?.(sim, tDisplay, { rebuilt: false });
      status.set('t', timeDependent ? tDisplay : 1);
      if (firstErr) { status.set('drift', firstErr.value); status.set('driftLabel', firstErr.label); }
    }
    if (++frames === 1 || frames % 30 === 0) status.set('frames', frames);   // first frame at once: slow pages count as live
    adaptResolution(now, performance.now() - frameStart);
    if (frameWaiters.length) { const w = frameWaiters.splice(0); w.forEach(f => f()); }
    requestAnimationFrame(frame);
  }

  window.psim = {
    model: 'generic',
    get t() { return tDisplay; },
    state: () => ({ t: sim.t, ...sim.state() }),
    params: () => ({ ...values }),
    advanceTo: t => {
      const T = tl ? Math.min(t, tl.duration) : t;
      while (sim.t < T && !isDone()) sim.step();
      tDisplay = isDone() ? endTime() : Math.max(tDisplay, Math.min(T, sim.t));
      return window.psim.state();
    },
    /** The run's declared span and key times ({ duration, markers, note }), or null. */
    get timeline() { return tl ? { start: tl.start, duration: tl.duration, markers: tl.markers.map(m => ({ ...m })), note: tl.note } : null; },
    get done() { return isDone(); },
    restart: () => rebuild(),
    setParams: (next, opts) => setParams(next, opts),
    setParam: (k, v) => setParam(k, v),
    getSim: () => sim,
    toggleFullscreen: () => toggleFullscreen(),
    /** Resolves after the next frame has been drawn: await it in a self-test before reading GPU pixels. */
    nextFrame: () => new Promise(res => frameWaiters.push(res)),
    report,
    get fullscreen() { return isFullscreen(); },
    testSteps,
    play: () => setPlaying(true), pause: () => setPlaying(false),
    steps: n => { if (typeof sim.step === "function") for (let i = 0; i < n; i++) sim.step(); tDisplay = sim.t ?? 0; return window.psim.state(); },
  };

  if (bench && benchParts) {
    benchApi = bench({
      ...benchParts, createPlot, createEquation, createDerivation, M, values, report,
      setParam, setParams, play: () => setPlaying(true), pause: () => setPlaying(false),
      getSim: () => sim,
    }) ?? null;
  }
  rebuild();
  const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
  setPlaying(timeDependent && (query.get('autoplay') ? query.get('autoplay') !== '0' : !reduced));
  addEventListener('keydown', e => {
    if (e.target.closest('input, select, textarea')) return;
    if ((e.key === 'f' || e.key === 'F') && !e.metaKey && !e.ctrlKey) { e.preventDefault(); toggleFullscreen(); return; }
    if (e.key === 'Escape' && pseudoFs) { setPseudo(false); return; }
    if (!timeDependent) return;
    if (e.code === 'Space') { e.preventDefault(); togglePlay(); }
    else if (e.key === '.') { setPlaying(false); if (!isDone()) { sim.step(); tDisplay = sim.t; } }
    else if (e.key === 'r' || e.key === 'R') rebuild();
  });
  status.set('model', 'generic');
  status.set('ready', true);
  requestAnimationFrame(frame);
  return { values, rebuild, view, get sim() { return sim; } };
}
