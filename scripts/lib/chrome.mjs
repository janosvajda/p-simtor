// Headless Chrome checks over the DevTools protocol, with no npm dependencies
// (Node >= 22 provides a global WebSocket). Pages report state through
// <html data-psim-*> attributes (ui/instrument.mjs beacon) and, for interaction
// tests, a window.psim scripting hook.

import { spawn } from 'node:child_process';
import { existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, basename } from 'node:path';
import { pathToFileURL } from 'node:url';

const CANDIDATES = [
  process.env.CHROME_PATH,
  '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
  '/Applications/Chromium.app/Contents/MacOS/Chromium',
  '/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge',
  '/usr/bin/google-chrome', '/usr/bin/google-chrome-stable', '/usr/bin/chromium', '/usr/bin/chromium-browser',
  'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe',
  'C:\\Program Files (x86)\\Google\\Chrome\\Application\\chrome.exe',
].filter(Boolean);

export function findChrome() {
  return CANDIDATES.find(p => existsSync(p)) ?? null;
}

const sleep = ms => new Promise(r => setTimeout(r, ms));

/** Launch headless Chrome and return a minimal CDP session on its first page. */
export async function launch({ args = [] } = {}) {
  const bin = findChrome();
  if (!bin) throw new Error('Chrome/Chromium not found. Set CHROME_PATH to the browser executable.');
  if (typeof WebSocket === 'undefined') throw new Error('Node >= 22 is required for browser checks (global WebSocket).');
  const profile = mkdtempSync(join(tmpdir(), 'psim-chrome-'));
  const proc = spawn(bin, [
    '--headless=new', '--remote-debugging-port=0', `--user-data-dir=${profile}`,
    '--no-first-run', '--no-default-browser-check', '--disable-extensions', '--mute-audio',
    '--hide-scrollbars', ...args, 'about:blank',
  ], { stdio: 'ignore' });
  const portFile = join(profile, 'DevToolsActivePort');
  for (let i = 0; i < 200 && !existsSync(portFile); i++) await sleep(50);
  if (!existsSync(portFile)) { proc.kill(); throw new Error('Chrome did not start (no DevToolsActivePort)'); }
  await sleep(100);
  const port = readFileSync(portFile, 'utf8').split('\n')[0].trim();
  let targets = [];
  for (let i = 0; i < 50 && !targets.some(t => t.type === 'page'); i++) {
    targets = await fetch(`http://127.0.0.1:${port}/json/list`).then(r => r.json()).catch(() => []);
    if (!targets.some(t => t.type === 'page')) await sleep(100);
  }
  const page = targets.find(t => t.type === 'page');
  const ws = new WebSocket(page.webSocketDebuggerUrl);
  await new Promise((res, rej) => { ws.onopen = res; ws.onerror = rej; });
  let id = 0;
  const pending = new Map(), listeners = [];
  ws.onmessage = ev => {
    const msg = JSON.parse(ev.data);
    if (msg.id && pending.has(msg.id)) {
      const { res, rej, timer } = pending.get(msg.id);
      clearTimeout(timer);
      pending.delete(msg.id);
      msg.error ? rej(new Error(msg.error.message)) : res(msg.result);
    } else if (msg.method) listeners.forEach(l => l(msg));
  };
  ws.onclose = () => {
    for (const { rej, timer } of pending.values()) { clearTimeout(timer); rej(new Error('Chrome connection closed')); }
    pending.clear();
  };
  const send = (method, params = {}, timeoutMs = 10000) => new Promise((res, rej) => {
    const mid = ++id;
    const timer = setTimeout(() => { pending.delete(mid); rej(new Error(`Chrome timed out: ${method}`)); }, timeoutMs);
    pending.set(mid, { res, rej, timer });
    ws.send(JSON.stringify({ id: mid, method, params }));
  });
  const once = method => new Promise(res => {
    const l = msg => { if (msg.method === method) { listeners.splice(listeners.indexOf(l), 1); res(msg.params); } };
    listeners.push(l);
  });
  const consoleErrors = [];
  listeners.push(msg => {
    if (msg.method === 'Runtime.exceptionThrown') consoleErrors.push(msg.params.exceptionDetails?.exception?.description ?? msg.params.exceptionDetails?.text);
    if (msg.method === 'Runtime.consoleAPICalled' && msg.params.type === 'error') consoleErrors.push(msg.params.args.map(a => a.value ?? a.description).join(' '));
  });
  await send('Page.enable');
  await send('Runtime.enable');
  return {
    send, consoleErrors,
    async goto(url, { width, height, mobile = false, scale = 2 }) {
      await send('Emulation.setDeviceMetricsOverride', { width, height, deviceScaleFactor: scale, mobile });
      const loaded = once('Page.loadEventFired');
      await send('Page.navigate', { url });
      let timer;
      try {
        await Promise.race([loaded, new Promise((_res, rej) => { timer = setTimeout(() => rej(new Error('Chrome page load timed out')), 10000); })]);
      } finally { clearTimeout(timer); }
    },
    async eval(expression, timeoutMs = 10000) {
      const r = await send('Runtime.evaluate', { expression, returnByValue: true, awaitPromise: true }, timeoutMs)
        .catch(e => { throw new Error(`${e.message} (${timeoutMs / 1000} s) while evaluating: ${expression.replace(/\s+/g, ' ').slice(0, 70)}…`); });
      if (r.exceptionDetails) throw new Error(r.exceptionDetails.exception?.description ?? r.exceptionDetails.text);
      return r.result.value;
    },
    async screenshot(file) {
      const r = await send('Page.captureScreenshot', { format: 'png' });
      writeFileSync(file, Buffer.from(r.data, 'base64'));
      return file;
    },
    async close() {
      try { await send('Browser.close', {}, 1000); } catch { /* closed or unresponsive: kill below */ }
      ws.close();
      await sleep(200);
      proc.kill();
      try { rmSync(profile, { recursive: true, force: true }); } catch { /* best effort */ }
    },
  };
}

const READ_ATTRS = `(() => {
  const out = {};
  for (const [k, v] of Object.entries(document.documentElement.dataset)) {
    if (!k.startsWith('psim')) continue;
    const key = k.slice(4, 5).toLowerCase() + k.slice(5);
    try { out[key] = JSON.parse(v); } catch { out[key] = v; }
  }
  return out;
})()`;

/**
 * Load the page in dark/light themes at desktop and phone sizes: collect runtime
 * errors, verify the simulation is live, check for horizontal overflow, save
 * screenshots. Then run model-specific interaction tests (see interactionTests).
 */
export async function checkPage(file, { shotsDir, waitMs = 3000, advance, params: given, preset } = {}) {
  let params = given;
  if (!findChrome()) {
    return { runs: [], problems: ['Chrome/Chromium not found. Set CHROME_PATH, or open the page manually and check the browser console.'] };
  }
  mkdirSync(shotsDir, { recursive: true });
  const name = basename(file, '.html');
  const variants = [
    { label: 'dark · desktop', q: 'theme=dark', width: 1440, height: 900, scroll: true },
    { label: 'light · desktop', q: 'theme=light', width: 1440, height: 900 },
    { label: 'dark · phone', q: 'theme=dark', width: 390, height: 844, mobile: true, scroll: true },
  ];
  const runs = [], problems = [];
  const b = await launch();
  try {
    for (const v of variants) {
      b.consoleErrors.length = 0;
      await b.goto(`${pathToFileURL(file).href}?${v.q}&autoplay=1`, v);
      if (preset && !params) {
        // a "Try this" preset by its button label: read the parameters it sets, then test that mode
        const r = await evalWithin(b, `(() => { const btn = [...document.querySelectorAll('.psim-presets button')].find(x => x.textContent.trim() === ${JSON.stringify(preset)}); if (!btn) return 'none'; btn.click(); return JSON.stringify(window.psim.params()); })()`, 60000);
        if (r.value === 'none' || r.error || r.frozen) {
          const summary = r.error ?? (r.frozen ? 'the preset froze the page' : `no preset button labelled "${preset}"`);
          problems.push(`${v.label}: ${summary}`); runs.push({ label: v.label, ok: false, summary });
          return { runs, problems };
        }
        params = JSON.parse(r.value);
      }
      if (params) {
        // a non-default mode or setting: set through the page, exactly as a preset would
        const r = await evalWithin(b, `window.psim?.setParams ? (window.psim.setParams(${JSON.stringify(params)}, { play: true }), JSON.stringify(window.psim.params())) : 'no setParams'`, 60000);
        if (r.frozen || r.error || r.value === 'no setParams') {
          const summary = r.error ?? (r.frozen ? 'setting --params froze the page' : 'the page has no window.psim.setParams');
          problems.push(`${v.label}: ${summary}`);
          runs.push({ label: v.label, ok: false, summary });
          return { runs, problems };
        }
      }
      await sleep(waitMs);
      if (advance != null) {
        // fast-forward the simulation (in the page's own time unit) and hold it at T, so the
        // screenshots show exactly that moment
        const r = await evalWithin(b, `window.psim?.advanceTo ? (window.psim.advanceTo(${+advance}), window.psim.pause?.(), window.psim.fit?.(), 'ok') : 'no advanceTo'`, 600000);
        if (r.frozen || r.error) {
          const summary = r.error ?? `advancing to t = ${advance} took longer than 10 min`;
          problems.push(`${v.label}: ${summary}`);
          runs.push({ label: v.label, ok: false, summary });
          return { runs, problems };
        }
        await sleep(600);
      }
      // A heavy page may need longer than waitMs for its first frames: wait for it to go live.
      // a run may start before zero (timeline.start < 0): live means the clock has moved past its start
      const t0 = await b.eval('window.psim?.timeline?.start ?? 0').catch(() => 0);
      let a = await b.eval(READ_ATTRS);
      for (let k = 0; k < 60 && a.ready === true && !((a.frames ?? 0) > 0 && (a.t ?? t0) > t0); k++) { await sleep(250); a = await b.eval(READ_ATTRS); }
      const overflow = await b.eval('document.documentElement.scrollWidth - document.documentElement.clientWidth');
      // name the widest culprits, so a phone overflow can be fixed without guessing
      const wide = overflow > 1 ? await b.eval(`(() => { const W = document.documentElement.clientWidth; return [...document.querySelectorAll('body *')].filter(e => e.getBoundingClientRect().right > W + 1 && ![...e.children].some(c => c.getBoundingClientRect().right > W + 1)).slice(0, 3).map(e => e.tagName.toLowerCase() + (e.className && typeof e.className === 'string' ? '.' + e.className.trim().split(/\\s+/).join('.') : '') + ' "' + (e.textContent || '').trim().slice(0, 40) + '"'); })()`) : [];
      // Every diagnostic must show a value: "—" or NaN means a quantity failed to compute in the page.
      const blank = await b.eval(`[...document.querySelectorAll('.psim-row')].filter(r => /^(—|NaN|[−-]?∞|[+-]?Infinity)(?:$|\\s)|NaN/.test(r.querySelector('.v')?.textContent.trim() ?? '')).map(r => r.querySelector('.k')?.textContent)`);
      const base = join(shotsDir, `${name}-${v.label.replace(/[^a-z]+/g, '-')}`);
      const shot = await b.screenshot(`${base}.png`);
      const more = v.scroll ? await scrollShots(b, base) : [];
      // a mid-run picture too (60 % through the timeline): the start alone is often the least telling moment
      if (v.label === 'dark · desktop' && advance == null) {
        const mid = await b.eval(`(() => { const tl = window.psim?.timeline; if (!tl) return false; window.psim.advanceTo(tl.start + 0.6 * (tl.duration - tl.start)); window.psim.pause?.(); return true; })()`, 120000).catch(() => false);
        if (mid) { await sleep(400); more.push(await b.screenshot(`${base}-midrun.png`)); }
      }
      const errors = [...(Array.isArray(a.errors) ? a.errors : []), ...b.consoleErrors];
      const live = (a.frames ?? 0) > 0 && (a.t ?? t0) > t0;
      const ok = a.ready === true && errors.length === 0 && live && overflow <= 1 && blank.length === 0;
      const summary = a.ready === undefined
        ? 'page did not report ready (script failed to start?)'
        : `frames ${a.frames ?? 0}, t = ${typeof a.t === 'number' ? a.t.toFixed(3) : '—'}${a.steps != null ? `, steps ${a.steps}` : ''}`
          + `${typeof a.drift === 'number' ? `, ${a.driftLabel ?? 'error'} ${a.drift.toExponential(2)}` : ''}`
          + `${a.ready === true && !live ? `, not running: the clock did not move past its start (${t0}) after ${waitMs} ms` : ''}`
          + `${overflow > 1 ? `, horizontal overflow ${overflow}px (widest: ${wide.join('; ')})` : ''}`
          + `${blank.length ? `, diagnostics without a value: ${blank.join(', ')}` : ''}`
          + `${errors.length ? `, errors: ${errors.join(' | ')}` : ''}`;
      if (!ok) problems.push(`${v.label}: ${summary}`);
      runs.push({ label: v.label, ok, summary, screenshot: shot, more });
    }
    for (const t of await stressInputs(b, file)) {
      runs.push({ label: t.label, ok: t.ok, summary: t.summary });
      if (!t.ok) problems.push(`${t.label}: ${t.summary}`);
    }
    // A failed stress test may have left the renderer frozen. Close it in finally,
    // without issuing further page commands that would wait on that renderer.
    if (problems.length) {
      runs.push({ label: 'later tests', ok: false, summary: 'skipped: the self-test and the UI→physics test run only when every view above passes' });
      return { runs, problems };
    }
    const model = await b.eval('window.psim?.model ?? null').catch(() => null);
    if (interactionTests[model]) {
      for (const t of await interactionTests[model](b, file, { params })) {
        runs.push({ label: t.label, ok: t.ok, summary: t.summary });
        if (!t.ok) problems.push(`${t.label}: ${t.summary}`);
      }
    }
    for (const t of await actionTests(b, file, { params, base: join(shotsDir, `${name}-dark-desktop`) })) {
      runs.push(t);
      if (!t.ok) problems.push(`${t.label}: ${t.summary}`);
    }
  } finally {
    await b.close();
  }
  return { runs, problems };
}

/**
 * Press every main action button (spec.actions) the way a visitor would: from 40 % through the
 * run, then play for two seconds and save a screenshot (<page>-dark-desktop-action-N.png), so what
 * the button does is seen and judged, not assumed. Fails on an error or a frozen page; a press
 * that changes neither the parameters nor the next drawn frame is reported for a look.
 */
export async function actionTests(b, file, { params, base, playMs = 2000 } = {}) {
  const out = [];
  const url = `${pathToFileURL(file).href}?theme=dark&autoplay=0`;
  await b.goto(url, { width: 1440, height: 900 });
  const labels = await b.eval(`[...document.querySelectorAll('.psim-action')].map(x => x.textContent.trim())`).catch(() => []);
  for (let i = 0; i < Math.min(labels.length, 6); i++) {
    if (i) await b.goto(url, { width: 1440, height: 900 });
    if (params) await b.eval(`window.psim?.setParams?.(${JSON.stringify(params)}, { play: false })`, 60000);
    await sleep(500);
    b.consoleErrors.length = 0;
    const frame = `(window.psim?.nextFrame?.() ?? Promise.resolve()).then(() => { try { return [...document.querySelectorAll('.psim-view canvas')].map(c => c.toDataURL()).join('|'); } catch { return ''; } })`;
    const before = await evalWithin(b, `(() => { const tl = window.psim?.timeline; window.psim?.restart?.(); if (tl) window.psim.advanceTo(tl.start + 0.4 * (tl.duration - tl.start)); window.psim?.pause?.(); return (${frame}).then(img => JSON.stringify({ img, p: window.psim?.params?.() ?? null, t: window.psim?.t ?? null })); })()`, 120000);
    const pressed = await evalWithin(b, `(() => { document.querySelectorAll('.psim-action')[${i}].click(); window.psim?.pause?.(); return (${frame}).then(img => JSON.stringify({ img, p: window.psim?.params?.() ?? null, t: window.psim?.t ?? null })); })()`, 60000);
    const label = `action "${labels[i]}" pressed`;
    if (before.frozen || before.error || pressed.frozen || pressed.error) {
      out.push({ label, ok: false, summary: pressed.error ?? before.error ?? 'the page froze when the button was pressed' });
      continue;
    }
    await b.eval('window.psim?.play?.()');
    await sleep(playMs);
    await b.eval('window.psim?.pause?.()');
    await sleep(300);
    const shot = await b.screenshot(`${base}-action-${i + 1}.png`);
    // and once its result should be plainly visible (spec.actions[i].reveal seconds of playback)
    const reveal = Math.min(30, Math.max(playMs / 1000, +(await b.eval(`window.psim?.actions?.[${i}]?.reveal ?? 8`)) || 8));
    let result = null;
    if (reveal * 1000 > playMs) {
      await b.eval('window.psim?.play?.()');
      await sleep(reveal * 1000 - playMs);
      await b.eval('window.psim?.pause?.()');
      await sleep(300);
      result = await b.screenshot(`${base}-action-${i + 1}-result.png`);
    }
    const a = JSON.parse(before.value ?? '{}'), z = JSON.parse(pressed.value ?? '{}');
    const errors = [...JSON.parse(await b.eval(`document.documentElement.dataset.psimErrors || '[]'`)), ...b.consoleErrors];
    const changed = JSON.stringify(a.p) !== JSON.stringify(z.p) || a.t !== z.t || (a.img && z.img && a.img !== z.img);
    const compared = changed || !!(a.img && z.img);
    out.push({ label, ok: errors.length === 0, screenshot: shot, more: result ? [result] : [],
      summary: errors.length ? `errors: ${errors.join(' | ')}`
        : `${changed ? 'changed the picture' : compared ? 'warning: the next frame, the time and the parameters are unchanged (fine if its effect grows while playing; look at the screenshot)' : 'pressed without errors (picture not comparable)'}; screenshots after ${playMs / 1000} s${result ? ` and ${reveal} s (result)` : ''} of play` });
  }
  return out;
}

/**
 * Screenshots of everything below the first screen: figures under the stage (workbench), the
 * control column, and the scrolling page on phones. Each scrollable area is paged through and
 * returned to the top afterwards.
 */
async function scrollShots(b, base, maxPages = 6) {
  const areas = await b.eval(`(() => {
    const named = [['figures', document.querySelector('.psim-app.workbench .psim-stage')], ['controls', document.querySelector('.psim-side')], ['page', document.scrollingElement]];
    return named.filter(([, e]) => e && e.scrollHeight - e.clientHeight > 40).map(([label, e], i) => { e.dataset.psimShot = i; return { label, i, pages: Math.ceil((e.scrollHeight - e.clientHeight) / (0.9 * e.clientHeight)) }; });
  })()`);
  const files = [];
  for (const { label, i, pages } of areas) {
    for (let k = 1; k <= Math.min(pages, maxPages); k++) {
      await b.eval(`(() => { const e = document.querySelector('[data-psim-shot="${i}"]'); e.scrollTop = ${k} * 0.9 * e.clientHeight; })()`);
      await sleep(350);
      files.push(await b.screenshot(`${base}-${label}-${k}.png`));
    }
    await b.eval(`document.querySelector('[data-psim-shot="${i}"]').scrollTop = 0`);
  }
  return files;
}

/** Evaluate with a deadline: a frozen page never answers. */
async function evalWithin(b, expr, ms) {
  let timer;
  const timeout = new Promise(res => { timer = setTimeout(() => res({ frozen: true }), ms); });
  const r = await Promise.race([b.eval(expr, ms + 1000).then(v => ({ value: v })).catch(e => ({ error: e.message })), timeout]);
  clearTimeout(timer);
  return r;
}

/**
 * Abuse every numeric input in the side panel (huge, negative, zero, nonsense): the page must
 * stay responsive and error-free. Catches unbounded parameters (e.g. 10⁸ stars) and inputs that
 * crash the model. Each input is tested on a fresh page load.
 */
export async function stressInputs(b, file, { timeoutMs = 3000, settleMs = 400, initialWaitMs = 300 } = {}) {
  const url = `${pathToFileURL(file).href}?autoplay=1`;
  await b.goto(url, { width: 1440, height: 900 });
  await sleep(initialWaitMs);
  const count = await b.eval(`document.querySelectorAll('.psim-side input[type=number]').length`);
  const bad = [];
  const result = () => [{
    label: 'extreme inputs',
    ok: bad.length === 0,
    summary: bad.length ? bad.slice(0, 4).join(' · ') : `${Math.min(count, 24)} numeric inputs × {10¹⁵, −10¹⁵, 0}: page stayed responsive, no errors`,
  }];
  for (let i = 0; i < Math.min(count, 24); i++) {
    for (const value of ['1e15', '-1e15', '0']) {
      b.consoleErrors.length = 0;
      const edited = await evalWithin(b, `(() => { const inp = document.querySelectorAll('.psim-side input[type=number]')[${i}]; inp.value = '${value}'; inp.dispatchEvent(new Event('change')); })()`, timeoutMs);
      if (edited.frozen || edited.error) {
        bad.push(`input #${i} = ${value}: ${edited.error ?? 'page froze during input change'}`);
        return result();
      }
      await sleep(settleMs);
      const alive = await evalWithin(b, `JSON.stringify({ label: document.querySelectorAll('.psim-side input[type=number]')[${i}].getAttribute('aria-label'), errors: JSON.parse(document.documentElement.dataset.psimErrors || '[]'), invalid: [...document.querySelectorAll('.psim-row')].filter(r => /^(—|NaN|[−-]?∞|[+-]?Infinity)(?:$|\\s)|NaN/.test(r.querySelector('.v')?.textContent.trim() ?? '')).map(r => r.querySelector('.k')?.textContent + ' = ' + r.querySelector('.v')?.textContent.trim()) })`, timeoutMs);
      if (alive.frozen || alive.error) {
        bad.push(`input #${i} = ${value}: page froze`);
        return result();
      }
      const info = JSON.parse(alive.value ?? '{}');
      const errs = [...(info.errors ?? []), ...b.consoleErrors];
      if (errs.length) { bad.push(`"${info.label}" = ${value}: ${errs[0]}`); break; }
      if (info.invalid?.length) { bad.push(`"${info.label}" = ${value}: diagnostic without a finite value (${info.invalid.join(', ')}; use { pending: true } for a quantity not measured yet)`); break; }
    }
    await b.goto(url, { width: 1440, height: 900 });
    await sleep(150);
  }
  return result();
}

/**
 * Model-specific end-to-end tests: drive the real UI, then compare the page's
 * computed state with an independent Node integration of the same setup.
 */
export const interactionTests = {
  async generic(b, file, { params } = {}) {
    await b.goto(`${pathToFileURL(file).href}?autoplay=0`, { width: 1440, height: 900 });
    // the mode under test (--params) applies to the self-test and the parameter-edit test too
    if (params) await b.eval(`window.psim?.setParams?.(${JSON.stringify(params)}, { play: false })`, 60000);
    await sleep(600);
    const out = [];
    // Optional page self-test (e.g. GPU output vs an exact formula): window.psim.selfTest()
    // returns { label, pass, ...measurements } or an array of them.
    if (await b.eval('typeof window.psim?.selfTest === "function"')) {
      // after a drawn frame; a failed GPU read is repeated once after two more frames (headless GPUs
      // occasionally return a blank buffer) and the first reading is kept in the summary
      const run = () => b.eval('(window.psim.nextFrame?.() ?? Promise.resolve()).then(() => window.psim.selfTest())', 120000)
        .then(r => [].concat(r)).catch(e => [{ label: 'page self-test', pass: false, error: e.message }]);
      const first = await run();
      let r = first, retried = false;
      if (first.some(t => t.pass !== true)) {
        await b.eval('window.psim.nextFrame?.().then(() => window.psim.nextFrame())', 10000).catch(() => {});
        r = await run(); retried = true;
      }
      r.forEach((t, i) => {
        const { label, pass, ...rest } = t;
        const note = retried && first[i]?.pass !== true && pass === true ? ` (passed on the 2nd attempt; 1st read ${JSON.stringify((({ label: _l, pass: _p, ...x }) => x)(first[i]))})` : '';
        out.push({ label: label ?? 'page self-test', ok: pass === true, summary: JSON.stringify(rest) + note });
      });
    }
    // Timer precision: a declared run span must end exactly there, with every key time inside it.
    const tlJson = await b.eval('JSON.stringify(window.psim?.timeline ?? null)');
    const tlv = JSON.parse(tlJson ?? 'null');
    const timeDependent = await b.eval('!!document.querySelector(".psim-transport")');
    if (tlv) {
      await b.eval('window.psim.restart()');
      const start = tlv.start ?? 0, span = tlv.duration - start, far = tlv.duration + 2 * span;
      const r = await evalWithin(b, `(window.psim.advanceTo(${far}), JSON.stringify({ t: window.psim.t, done: window.psim.done }))`, 600000);
      const end = JSON.parse(r.value ?? '{}'), exact = Math.abs((end.t ?? NaN) - tlv.duration) <= 1e-9 * Math.max(Math.abs(tlv.duration), span);
      const inside = tlv.markers.every(m => m.t >= start && m.t <= tlv.duration);
      out.push({ label: 'timeline: run ends exactly at its declared end', ok: !!end.done && exact && inside,
        summary: `${r.error ?? ''}run ${start}–${tlv.duration}; after advancing to the end + 2× the span: t = ${end.t}, done = ${end.done}; key times: ${tlv.markers.map(m => `${m.label} @ ${m.t}`).join(', ') || 'none'}` });
    } else if (timeDependent) {
      out.push({ label: 'timeline', ok: true, summary: 'warning: no spec.timeline declared. Declare the run span and its key times (computed and requested) so the clock shows when things happen and stops at the end' });
    }
    // A visitor without a physics degree must learn what they are looking at: a plain intro.
    const intro = JSON.parse(await b.eval('JSON.stringify(window.psim?.intro ?? null)') ?? 'null');
    if (intro !== null || await b.eval('!!window.psim?.defaultSpeed')) {
      const words = intro?.text ? intro.text.trim().split(/\s+/).length : 0;
      out.push({ label: 'plain-language intro (spec.intro)', ok: words >= 12 && words <= 90 && (intro?.watch ?? []).length <= 4,
        summary: intro ? `${words} words, ${(intro.watch ?? []).length} "look for" items${words > 90 ? ': too long, keep it under 90 words' : words < 12 ? ': too short to explain what is shown' : ''}` : 'missing: add intro: { text, watch } in everyday words (what is shown, why it matters, what to look for)' });
    }
    // The page's main action (spec.actions) must be on the first screen, under the view.
    const act = await b.eval(`(() => { const a = [...document.querySelectorAll('.psim-action')]; if (!a.length) return null; window.scrollTo(0, 0); return { n: a.length, visible: a.every(x => { const r = x.getBoundingClientRect(); return r.bottom <= innerHeight && r.top >= 0; }) }; })()`);
    if (act) out.push({ label: 'main action visible on the first screen (1440×900)', ok: act.visible, summary: `${act.n} action button(s)${act.visible ? '' : ': below the fold, make the view shorter or the page will hide its main action'}` });
    // The run must be watchable: at the default speed it should take long enough to follow.
    if (tlv) {
      const speed = await b.eval('window.psim.defaultSpeed ?? 0');
      const secs = speed > 0 ? (tlv.duration - (tlv.start ?? 0)) / speed : Infinity;
      out.push({ label: 'playback length at default speed', ok: secs >= 30,
        summary: `${Number.isFinite(secs) ? secs.toFixed(1) : '∞'} s for the whole run${secs < 30 ? ': too short to follow, aim for 45–120 s (lower the default speed)' : secs < 45 ? ' (short; 45–120 s reads better)' : ''}` });
    }
    // Smoothness: frames per second while playing, at a laptop window size.
    if (timeDependent) {
      await b.eval('window.psim.restart(); window.psim.play?.()');
      await sleep(1500);
      const fps = await b.eval(`new Promise(res => { let n = 0; const t0 = performance.now(); const f = () => { n++; const dt = performance.now() - t0; if (dt < 3000) requestAnimationFrame(f); else res(n / (dt / 1000)); }; requestAnimationFrame(f); })`, 20000);
      const limited = await b.eval(`/compute-limited/.test(document.querySelector('.psim-hud .flag')?.textContent ?? '')`);
      out.push({ label: 'smooth playback (frames per second, 1440×900)', ok: fps >= 24,
        summary: `${fps.toFixed(0)} fps${limited ? ', compute-limited (the physics cannot keep up with the default speed)' : ''}${fps < 24 ? ': too slow, reduce the per-frame work or the resolution' : fps < 45 ? ' (below 45: check on a laptop)' : ''}` });
      await b.eval('window.psim.pause?.()');
    }
    const input = '.psim-side .psim-field input[type=number]';
    if (!(await b.eval(`!!document.querySelector('${input}')`))) {
      return [...out, { label: 'parameter edit → model', ok: true, summary: 'no numeric parameters to test' }];
    }
    await b.eval('window.psim.restart()');
    const n = await b.eval('window.psim.testSteps ?? 300');
    const s1 = await b.eval(`window.psim.steps(${n})`, 120000);
    const p0 = await b.eval('window.psim.params()');
    // a new value inside [min, max]: ×1.3, or towards the middle of the range when that would leave it
    await b.eval(`(() => {
      const i = document.querySelector('${input}'), v = +i.value;
      const lo = i.min === '' ? -Infinity : +i.min, hi = i.max === '' ? Infinity : +i.max;
      let next = v * 1.3 || 1;
      if (!(next >= lo && next <= hi) || next === v) next = Number.isFinite(lo) && Number.isFinite(hi) ? (v - lo > hi - v ? v - 0.3 * (v - lo) : v + 0.3 * (hi - v)) : v * 0.7 || -1;
      const st = +i.step > 0 ? +i.step : 0;
      if (st) { next = Math.round(next / st) * st; if (next === v) next = v + st <= hi ? v + st : v - st; }
      i.value = String(next); i.dispatchEvent(new Event('change'));
    })()`);
    const p1 = await b.eval('window.psim.params()');
    const s2 = await b.eval(`window.psim.steps(${n})`, 120000);
    const changedKey = Object.keys(p0).find(k => p0[k] !== p1[k]);
    const differs = JSON.stringify(s1) !== JSON.stringify(s2);
    // Blackbody's scale selector is a display control: changing only the scale
    // must redraw its spectrum even though the temperature stays the same.
    // Only a lin/log *parameter* (Parameters section) is expected to redraw the first figure;
    // View-section controls may act on the 3D scene instead.
    const PARAM_SELECTS = `[...(document.querySelector('.psim-side details')?.querySelectorAll('select') ?? [])]`;
    const hasScale = await b.eval(`${PARAM_SELECTS}.some(s => [...s.options].some(o => o.value === 'lin') && [...s.options].some(o => o.value === 'log'))`);
    if (hasScale) {
      const before = await b.eval(`document.querySelector('.psim-figure canvas').toDataURL()`);
      await b.eval(`(() => { const s = ${PARAM_SELECTS}.find(q => [...q.options].some(o => o.value === 'lin')); s.value = s.value === 'lin' ? 'log' : 'lin'; s.dispatchEvent(new Event('change')); })()`);
      await sleep(200);
      const after = await b.eval(`document.querySelector('.psim-figure canvas').toDataURL()`);
      out.push({ label: 'spectrum axis control', ok: before !== after, summary: before !== after ? 'changing only the axis scale redrew the spectrum' : 'axis selector left the plot unchanged' });
    }
    // in a non-default mode (--params / --preset) a parameter may legitimately do nothing (a field
    // that only acts on polarised particles, with polarisation off): a warning there, a failure by default
    const modeOnly = !!params && !!changedKey && !differs;
    return [...out, {
      label: 'parameter edit → model',
      ok: (!!changedKey && differs) || modeOnly,
      summary: changedKey
        ? `${changedKey}: ${p0[changedKey]} → ${p1[changedKey]} via the UI; recomputed state after ${n} steps ${differs ? 'differs' : modeOnly ? 'is identical in this mode (warning: confirm it has no effect here; the default mode must pass)' : 'is IDENTICAL (parameter not wired to the physics?)'}`
        : 'editing the first parameter input did not change the model parameters',
    }];
  },
  async nbody(b, file) {
    const { integrateTo } = await import('../../core/nbody.mjs');
    const results = [];
    await b.goto(`${pathToFileURL(file).href}?autoplay=0`, { width: 1440, height: 900 });
    await sleep(300);
    const compare = async tEnd => {
      const setup = await b.eval('window.psim.setup()');
      const page = await b.eval(`window.psim.advanceTo(${tEnd})`);
      // Same initial value problem and step parameter, integrated independently in Node.
      const ref = integrateTo(setup.bodies,
        { integrator: setup.integrator, h: setup.h, G: setup.G, softening: setup.softening, barycentric: false }, page.t);
      let d = 0;
      for (let i = 0; i < page.x.length; i++) d = Math.max(d, Math.abs(page.x[i] - ref.sys.x[i]));
      return { setup, page, d };
    };
    const base = await compare(2);
    results.push({
      label: 'page ≡ Node integrator',
      ok: base.d === 0,
      summary: `max |x_page − x_node| = ${base.d.toExponential(2)} at t = ${base.page.t.toFixed(4)} (${base.page.steps} steps)`,
    });
    // Editing either solver control must withdraw the old report and trust label.
    for (const control of ['integrator', 'accuracy']) {
      const selector = control === 'integrator'
        ? "[...document.querySelectorAll('select')].find(s => [...s.options].some(o => o.value === 'leapfrog'))"
        : "document.querySelector('input[aria-label=\"Steps per crossing time value\"]')";
      await b.eval(`(() => {
        const input = ${selector};
        input.value = ${control === 'integrator' ? JSON.stringify(base.setup.integrator === 'leapfrog' ? 'logh4' : 'leapfrog') : String(base.setup.stepsPerCrossing * 2)};
        input.dispatchEvent(new Event('change'));
      })()`);
      const withdrawn = await b.eval(`document.getElementById('tier').dataset.tier === 'unverified' && document.querySelector('.psim-side').textContent.includes('build-time report does not apply')`);
      results.push({ label: `${control} → report`, ok: withdrawn, summary: withdrawn ? 'changed solver settings withdrew the old verification and trust label' : 'stale verification remained visible' });
      await b.goto(`${pathToFileURL(file).href}?autoplay=0`, { width: 1440, height: 900 });
      await sleep(300);
    }
    // Change m₁ through the real number input, as a user would.
    await b.eval('window.psim.load(window.psim.scenarios()[0])');
    const m0 = (await b.eval('window.psim.setup()')).bodies[0].m;
    const before = await b.eval('window.psim.advanceTo(2)');
    await b.eval(`(() => {
      const input = document.querySelector('.psim-bodycard .psim-field input[type=number]');
      input.value = String(+input.value * 1.5);
      input.dispatchEvent(new Event('change'));
    })()`);
    const changed = await compare(2);
    let moved = 0;
    for (let i = 0; i < before.x.length; i++) moved = Math.max(moved, Math.abs(before.x[i] - changed.page.x[i]));
    const m1 = changed.setup.bodies[0].m;
    results.push({
      label: 'mass edit → physics',
      ok: Math.abs(m1 / m0 - 1.5) < 1e-9 && changed.page.t > 0 && moved > 1e-3 && changed.d === 0,
      summary: `m₁ ${m0} → ${m1} via the UI moved positions at t = 2 by ${moved.toExponential(2)}; page vs Node with the new mass: ${changed.d.toExponential(2)}`,
    });
    return results;
  },
};
