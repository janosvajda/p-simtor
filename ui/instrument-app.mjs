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
  el, sci, fix, grade, initTheme, section, row, slider, toggle, select, button, toast,
  transport, createClock, logPlot, download, readJsonScript, beacon, cssVar,
} from './instrument.mjs';
import { createPlot } from './plot.mjs';
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
    timeUnit = 'T', timeLabel = 'simulated time', speed: speed0 = 1, controls, onPointer, layout = 'stage', bench,
  } = spec;
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
  $('.psim-hud .dim').textContent = timeLabel;
  const hud = { t: $('.psim-hud b'), a: $('.psim-hud .h1'), b: $('.psim-hud .h2'), flag: $('.psim-hud .flag') };
  if (model.legend) $('.psim-legend').append(...model.legend.map(t => el('span', {}, t)));

  // ------------------------------------------------------------ state
  const values = Object.fromEntries(params.map(p => [p.key, p.value]));
  let sim, tDisplay = 0, limited = false;
  const clock = createClock({ speed: speed0 });
  const view = { W: 0, H: 0, dpr: 1, color: cssVar, opts: {}, pointer: null };

  let benchApi = null;
  function rebuild() {
    sim = create({ ...values });
    tDisplay = sim.t ?? 0;
    diagPlot.clear();
    renderDiag(true);
    benchApi?.update?.(sim, tDisplay, { rebuilt: true });
    status.set('params', values);
  }

  // ------------------------------------------------------------ transport
  const timeDependent = spec.timeDependent !== false;
  let tp = null;
  if (timeDependent) {
    tp = transport({
      speed: speed0, speedLabel: s => `${sci(s, 2)} ${timeUnit}/s`,
      onToggle: () => setPlaying(!clock.playing),
      onStep: () => { setPlaying(false); sim.step(); tDisplay = sim.t; },
      onRestart: () => rebuild(),
      onFit: () => spec.onFit?.(view),
      onSpeed: s => { clock.speed = s; },
    });
    stage.append(tp);
  } else {
    stage.querySelector('.psim-hud').querySelector('.dim').textContent = '';
  }
  function setPlaying(p) { clock.playing = p; tp?.setPlaying(p); }

  // ------------------------------------------------------------ side panel
  const secParams = section('Parameters');
  const paramInputs = new Map();
  if (params.length) {
    secParams.body.append(el('p', { class: 'psim-hint' }, 'Changing a parameter recomputes the model from its initial state.'));
    for (const p of params) {
      if (p.options) {
        secParams.body.append(el('div', { class: 'psim-field' }, el('div', { class: 'top' }, el('span', {}, p.label)),
          select(p.options, p.value, v => { values[p.key] = isNaN(+v) ? v : +v; rebuild(); })));
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
  if (controls) controls(secParams.body, { values, rebuild });

  const secDiag = section('Diagnostics', { aside: 'live' });
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
      const text = x.text ?? (x.kind === 'error' ? sci(x.value) : `${sci(x.value, 4)}${x.unit ? ` ${x.unit}` : ''}`);
      r.set(text, x.kind === 'error' && x.good != null ? grade(x.value, x.good, x.warn ?? x.good * 1000) : '');
    }
    const firstErr = d.find(x => x.kind === 'error');
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
    button('PNG image', () => canvas.toBlob(b => el('a', { href: URL.createObjectURL(b), download: `${slug}.png` }).click())),
  ));

  side.append(...[secParams, secDiag, secVerify, secModel, secExport].filter(s => s !== secParams || params.length || controls));
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

  let frames = 0, lastUi = 0;
  function frame(now) {
    if (timeDependent) {
      const target = clock.tick(now, tDisplay);
      if (target > tDisplay) {
        const deadline = performance.now() + 9;
        let k = 0;
        while (sim.t < target) {
          sim.step();
          if ((++k & 63) === 0 && performance.now() > deadline) break;
        }
        limited = sim.t < target;
        tDisplay = Math.min(target, sim.t);
      }
    }
    ctx.setTransform(view.dpr, 0, 0, view.dpr, 0, 0);
    ctx.clearRect(0, 0, view.W, view.H);
    draw(ctx, sim, view, tDisplay);
    if (now - lastUi > 100) {
      lastUi = now;
      const d = renderDiag();
      hud.t.textContent = timeDependent ? `${fix(tDisplay, 3)} ${timeUnit}` : '';
      const firstErr = d.find(x => x.kind === 'error');
      hud.a.textContent = firstErr ? `${firstErr.label} ${sci(firstErr.value)}` : '';
      hud.b.textContent = spec.hudLine?.(sim) ?? '';
      hud.flag.textContent = limited ? 'compute-limited: playback slowed to keep the calculation exact' : '';
      benchApi?.update?.(sim, tDisplay, { rebuilt: false });
      status.set('t', timeDependent ? tDisplay : 1);
      if (firstErr) status.set('drift', firstErr.value);
    }
    if (++frames === 1 || frames % 30 === 0) status.set('frames', frames);   // first frame at once: slow pages count as live
    requestAnimationFrame(frame);
  }

  window.psim = {
    model: 'generic',
    get t() { return tDisplay; },
    state: () => ({ t: sim.t, ...sim.state() }),
    params: () => ({ ...values }),
    advanceTo: t => { while (sim.t < t) sim.step(); tDisplay = Math.max(tDisplay, Math.min(t, sim.t)); return window.psim.state(); },
    restart: () => rebuild(),
    steps: n => { if (typeof sim.step === "function") for (let i = 0; i < n; i++) sim.step(); tDisplay = sim.t ?? 0; return window.psim.state(); },
  };

  if (bench && benchParts) {
    benchApi = bench({
      ...benchParts, createPlot, createEquation, createDerivation, M, values,
      setParam: (k, v) => { values[k] = v; paramInputs.get(k)?.setValue?.(v); rebuild(); },
      getSim: () => sim,
    }) ?? null;
  }
  rebuild();
  const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
  setPlaying(timeDependent && (query.get('autoplay') ? query.get('autoplay') !== '0' : !reduced));
  addEventListener('keydown', e => {
    if (e.target.closest('input, select, textarea') || !timeDependent) return;
    if (e.code === 'Space') { e.preventDefault(); setPlaying(!clock.playing); }
    else if (e.key === '.') { setPlaying(false); sim.step(); tDisplay = sim.t; }
    else if (e.key === 'r' || e.key === 'R') rebuild();
  });
  status.set('model', 'generic');
  status.set('ready', true);
  requestAnimationFrame(frame);
  return { values, rebuild, view, get sim() { return sim; } };
}
