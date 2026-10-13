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
export async function launch() {
  const bin = findChrome();
  if (!bin) throw new Error('Chrome/Chromium not found. Set CHROME_PATH to the browser executable.');
  if (typeof WebSocket === 'undefined') throw new Error('Node >= 22 is required for browser checks (global WebSocket).');
  const profile = mkdtempSync(join(tmpdir(), 'psim-chrome-'));
  const proc = spawn(bin, [
    '--headless=new', '--remote-debugging-port=0', `--user-data-dir=${profile}`,
    '--no-first-run', '--no-default-browser-check', '--disable-extensions', '--mute-audio',
    '--hide-scrollbars', 'about:blank',
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
      const { res, rej } = pending.get(msg.id);
      pending.delete(msg.id);
      msg.error ? rej(new Error(msg.error.message)) : res(msg.result);
    } else if (msg.method) listeners.forEach(l => l(msg));
  };
  const send = (method, params = {}) => new Promise((res, rej) => {
    const mid = ++id;
    pending.set(mid, { res, rej });
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
      await loaded;
    },
    async eval(expression) {
      const r = await send('Runtime.evaluate', { expression, returnByValue: true, awaitPromise: true });
      if (r.exceptionDetails) throw new Error(r.exceptionDetails.exception?.description ?? r.exceptionDetails.text);
      return r.result.value;
    },
    async screenshot(file) {
      const r = await send('Page.captureScreenshot', { format: 'png' });
      writeFileSync(file, Buffer.from(r.data, 'base64'));
      return file;
    },
    async close() {
      try { await send('Browser.close'); } catch { /* already closing */ }
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
export async function checkPage(file, { shotsDir, waitMs = 3000 } = {}) {
  if (!findChrome()) {
    return { runs: [], problems: ['Chrome/Chromium not found. Set CHROME_PATH, or open the page manually and check the browser console.'] };
  }
  mkdirSync(shotsDir, { recursive: true });
  const name = basename(file, '.html');
  const variants = [
    { label: 'dark · desktop', q: 'theme=dark', width: 1440, height: 900 },
    { label: 'light · desktop', q: 'theme=light', width: 1440, height: 900 },
    { label: 'dark · phone', q: 'theme=dark', width: 390, height: 844, mobile: true },
  ];
  const runs = [], problems = [];
  const b = await launch();
  try {
    for (const v of variants) {
      b.consoleErrors.length = 0;
      await b.goto(`${pathToFileURL(file).href}?${v.q}&autoplay=1`, v);
      await sleep(waitMs);
      const a = await b.eval(READ_ATTRS);
      const overflow = await b.eval('document.documentElement.scrollWidth - document.documentElement.clientWidth');
      const shot = await b.screenshot(join(shotsDir, `${name}-${v.label.replace(/[^a-z]+/g, '-')}.png`));
      const errors = [...(Array.isArray(a.errors) ? a.errors : []), ...b.consoleErrors];
      const live = (a.frames ?? 0) > 0 && (a.t ?? 0) > 0;
      const ok = a.ready === true && errors.length === 0 && live && overflow <= 1;
      const summary = a.ready === undefined
        ? 'page did not report ready (script failed to start?)'
        : `frames ${a.frames ?? 0}, t = ${typeof a.t === 'number' ? a.t.toFixed(3) : '—'}, steps ${a.steps ?? '—'}`
          + `${typeof a.drift === 'number' ? `, |ΔE/E0| ${a.drift.toExponential(2)}` : ''}`
          + `${overflow > 1 ? `, horizontal overflow ${overflow}px` : ''}`
          + `${errors.length ? `, errors: ${errors.join(' | ')}` : ''}`;
      if (!ok) problems.push(`${v.label}: ${summary}`);
      runs.push({ label: v.label, ok, summary, screenshot: shot });
    }
    const model = await b.eval('window.psim?.model ?? null').catch(() => null);
    if (interactionTests[model]) {
      for (const t of await interactionTests[model](b, file)) {
        runs.push({ label: t.label, ok: t.ok, summary: t.summary });
        if (!t.ok) problems.push(`${t.label}: ${t.summary}`);
      }
    }
  } finally {
    await b.close();
  }
  return { runs, problems };
}

/**
 * Model-specific end-to-end tests: drive the real UI, then compare the page's
 * computed state with an independent Node integration of the same setup.
 */
export const interactionTests = {
  async generic(b, file) {
    await b.goto(`${pathToFileURL(file).href}?autoplay=0`, { width: 1440, height: 900 });
    await sleep(300);
    const input = '.psim-side .psim-field input[type=number]';
    if (!(await b.eval(`!!document.querySelector('${input}')`))) {
      return [{ label: 'parameter edit → model', ok: true, summary: 'no numeric parameters to test' }];
    }
    await b.eval('window.psim.restart()');
    const s1 = await b.eval('window.psim.steps(300)');
    const p0 = await b.eval('window.psim.params()');
    await b.eval(`(() => { const i = document.querySelector('${input}'); i.value = String(+i.value * 1.3 || 1); i.dispatchEvent(new Event('change')); })()`);
    const p1 = await b.eval('window.psim.params()');
    const s2 = await b.eval('window.psim.steps(300)');
    const changedKey = Object.keys(p0).find(k => p0[k] !== p1[k]);
    const differs = JSON.stringify(s1) !== JSON.stringify(s2);
    return [{
      label: 'parameter edit → model',
      ok: !!changedKey && differs,
      summary: changedKey
        ? `${changedKey}: ${p0[changedKey]} → ${p1[changedKey]} via the UI; recomputed state after 300 steps ${differs ? 'differs' : 'is IDENTICAL (parameter not wired to the physics?)'}`
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
