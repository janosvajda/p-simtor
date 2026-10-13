// Regressions from the read-only review: test the original failing inputs and
// observable results, including CLI failures and a frozen browser input handler.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, mkdtempSync, writeFileSync, rmSync, existsSync } from 'node:fs';
import { spawnSync, execFileSync } from 'node:child_process';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import vm from 'node:vm';
import '../scripts/register.mjs';
import { dopri5 } from '../core/numerics.mjs';
import { deflectionExact, deflectionWeak, relativisticOrbit, precessionPerOrbit } from '../core/schwarzschild.mjs';
import { circularSpeed, acceleration, MW, G, NFW_SCALE } from '../core/galaxy.mjs';
import { createLattice, sweep } from '../core/lattice-gauge.mjs';
import { keplerBodies, keplerRelative } from '../core/kepler.mjs';
import { integrateTo, createSystem, createStepper, toBarycentric, pairCrossingTime } from '../core/nbody.mjs';
import { solarSystemScenarios } from '../core/solar-system.mjs';
import { THREE_BODY_PRESETS } from '../core/presets.mjs';
import { verifyScenario } from '../core/nbody-verify.mjs';
import { validateModuleReport, gravityReportProblems, gravityReportMatches } from '../core/verification.mjs';
import { stressInputs } from '../scripts/lib/chrome.mjs';
import { M } from '../ui/equation.mjs';

const ROOT = fileURLToPath(new URL('../', import.meta.url));

test('Dormand–Prince enforces tolerance under changes of time units', () => {
  for (const scale of [1, 1e16, 1e-16]) {
    const r = dopri5((_t, y) => [scale * y[0]], 0, [1], 10 / scale, { h0: 10 / scale, rtol: 1e-9, atol: 1e-12 });
    assert.equal(r.t, 10 / scale);
    assert.ok(r.rejected > 0 && r.steps > 100);
    assert.ok(Math.abs(r.y[0] / Math.exp(10) - 1) < 3e-9, `scale=${scale}: ${r.y[0]}`);
  }
});

test('Dormand–Prince distinguishes failure from an intentional early stop', () => {
  const f = (_t, y) => [y[0]];
  assert.throws(() => dopri5(f, 0, [1], 1, { maxSteps: 1 }), /budget exhausted/);
  for (const opts of [{ h0: 0 }, { hmax: 0 }, { rtol: NaN }, { atol: -1 }]) assert.throws(() => dopri5(f, 0, [1], 1, opts));
  assert.throws(() => dopri5(f, 1, [1], 2, { h0: Number.MIN_VALUE }), /underflow/);
  const r = dopri5(f, 0, [1], 1, { onStep: t => t < 0.1 });
  assert.equal(r.stopped, true);
  assert.ok(r.t >= 0.1 && r.t < 1);
  const back = dopri5(f, 1, [Math.E], 0);
  assert.equal(back.t, 0);
  assert.ok(Math.abs(back.y[0] - 1) < 1e-9);
});

test('exact deflection stays positive and accurate in arbitrarily weak fields', () => {
  // The first omitted weak-field term is O(b^-4); beyond 10^5 M it is negligible.
  for (const b of [1e5, 471127.6, 1e7, 1e9, 1e16]) {
    const ref = deflectionWeak(b);
    for (const n of [32, 64, 256]) {
      const value = deflectionExact(b, n);
      assert.ok(value > 0 && Math.abs(value / ref - 1) < 1e-12, `b=${b}, n=${n}: ${value} vs ${ref}`);
    }
  }
  assert.equal(deflectionExact(Infinity), 0);
});

test('strong-field orbits complete every requested radial period and preserve coordinate time', () => {
  for (const [p, e] of [[10, 0.5], [7.1, 0.5], [6.9, 0.4]]) {
    const r = relativisticOrbit(p, e, { orbits: 3 });
    const period = 2 * Math.PI + precessionPerOrbit(p, e);
    assert.equal(r.periapsides.length, 4);
    assert.equal(r.periapsisTimes.length, 4);
    for (let k = 1; k <= 3; k++) {
      assert.ok(Math.abs(r.periapsides[k] - k * period) < 1e-7, `p=${p}, passage ${k}`);
      assert.ok(r.periapsisTimes[k] > r.periapsisTimes[k - 1]);
      assert.ok(Math.abs(r.at(r.periapsisTimes[k]).r / (p / (1 + e)) - 1) < 1e-5);
    }
  }
});

test('circular relativistic orbits count revolutions without spurious periapsis detection', () => {
  const r = relativisticOrbit(10, 0, { orbits: 3 });
  assert.equal(r.periapsides.length, 4);
  assert.ok(Math.abs(r.duration / (6 * Math.PI * 10 ** 1.5) - 1) < 1e-15);
  assert.ok(r.points.every(q => q[1] === 10));
  assert.throws(() => relativisticOrbit(7, 0.5), /separatrix/);
  assert.throws(() => relativisticOrbit(10, 1));
});

test('zero and tiny lattice coupling sample Haar measure without hanging', () => {
  // A subprocess deadline also makes a recurrence of the old infinite loop fail.
  const code = `
    import { createLattice, sweep, plaquette } from './core/lattice-gauge.mjs';
    const out = [];
    for (const beta of [0, 1e-12]) {
      const lat = createLattice({ L: 8, d: 2, beta, seed: 19 });
      let mean = 0, second = 0, p = 0, count = 0, unitError = 0;
      for (let k = 0; k < 64; k++) {
        sweep(lat, { overrelax: 0 }); p += plaquette(lat) / 64;
        for (let i = 0; i < lat.U.length; i += 4) {
          const a = lat.U[i]; mean += a; second += a*a; count++;
          unitError = Math.max(unitError, Math.abs(Math.hypot(...lat.U.subarray(i,i+4))-1));
        }
      }
      out.push({ beta, mean: mean/count, second: second/count, p, unitError });
    }
    console.log(JSON.stringify(out));`;
  const values = JSON.parse(execFileSync(process.execPath, ['--input-type=module', '-e', code], { cwd: ROOT, encoding: 'utf8', timeout: 5000 }));
  for (const r of values) {
    assert.ok(Math.abs(r.mean) < 0.02);
    assert.ok(Math.abs(r.second - 0.25) < 0.02, JSON.stringify(r));
    assert.ok(Math.abs(r.p) < 0.03);
    assert.ok(r.unitError < 1e-14);
  }
  for (const beta of [-1, NaN, Infinity]) assert.throws(() => createLattice({ beta }));
  const lat = createLattice(); lat.beta = NaN;
  assert.throws(() => sweep(lat));
});

test('galaxy circular speed has the correct central limit', () => {
  assert.deepEqual(circularSpeed(0), { total: 0, bulge: 0, disk: 0, halo: 0 });
  const halo = { bulge: false, disk: false, halo: true };
  const centralForce = G * MW.Mh / (2 * NFW_SCALE ** 2);
  for (const r of [1e-9, 1e-6, 1e-3]) {
    const a = acceleration(r, 0, 0, halo, MW, new Float64Array(3));
    assert.ok(Math.abs(circularSpeed(r, halo).halo ** 2 / (-r * a[0]) - 1) < 1e-12);
    assert.ok(Math.abs(-a[0] / centralForce - 1) < 1e-4);
  }
});

test('conservation tracking catches an error peak that an endpoint check misses', () => {
  const bodies = keplerBodies(1, 0.5, 1, 0.6);
  const period = keplerRelative(1, 0.6, 1.5, 0).period;
  const r = integrateTo(bodies, { integrator: 'leapfrog', h: 0.01, trackConservation: true }, period);
  assert.ok(r.maxDrift.dE > 2 * r.dE, `${r.maxDrift.dE} vs endpoint ${r.dE}`);
  assert.ok(r.maxDrift.dP >= r.dP && r.maxDrift.dL >= r.dL);
});

test('N-body step selection resolves extreme Solar System velocity edits', () => {
  const sc = solarSystemScenarios()[0], Gsolar = sc.units.G;
  const initial = toBarycentric(createSystem(sc.bodies, { G: Gsolar }));
  const bodies = sc.bodies.map((b, i) => ({ ...b, x: Array.from(initial.x.subarray(3 * i, 3 * i + 3)), v: Array.from(initial.v.subarray(3 * i, 3 * i + 3)) }));
  const vmax = Math.max(...sc.bodies.map(b => Math.hypot(...b.v)));
  bodies[2].v[0] = 100 * vmax;
  const sys = toBarycentric(createSystem(bodies, { G: Gsolar }));
  assert.ok(pairCrossingTime(sys) < 0.001);
  const st = createStepper(sys, { integrator: 'logh4', stepsPerCrossing: 400 });
  while (sys.t < 0.2) { const before = sys.t; st.step(); assert.ok(sys.t > before); }
  assert.ok(sys.x.every(Number.isFinite) && sys.v.every(Number.isFinite));
  assert.equal(pairCrossingTime(createSystem([{ m: 1, x: [0, 0, 0], v: [0, 0, 0] }, { m: 1, x: [1, 0, 0], v: [0, 0, 0] }])), Infinity);
});

test('a singular N-body step restores the last valid state and reports failure', () => {
  const sys = createSystem([{ m: 1, x: [-1, 0, 0], v: [1, 0, 0] }, { m: 1, x: [1, 0, 0], v: [-1, 0, 0] }], { G: 1e-30 });
  const before = { x: sys.x.slice(), v: sys.v.slice(), t: sys.t };
  const st = createStepper(sys, { integrator: 'leapfrog', h: 1 });
  assert.throws(() => st.step(), /numerical step failed/);
  assert.deepEqual(sys.x, before.x); assert.deepEqual(sys.v, before.v);
  assert.equal(sys.t, before.t); assert.equal(sys.steps, 0);
  assert.ok(sys.a.every(Number.isFinite));
});

test('chaotic gravity recipes resolve conservation error throughout their runs', () => {
  for (const name of ['pythagorean', 'lagrange']) {
    const p = THREE_BODY_PRESETS[name];
    const report = verifyScenario({ ...p, integrator: 'logh4' });
    assert.ok(report.conservation.dE < 1e-8, `${name}: ${report.conservation.dE}`);
    assert.deepEqual(gravityReportProblems(report), []);
    assert.ok(report.convergence.orderH > report.conservation.h);
    assert.ok(Math.abs(report.convergence.observedOrder - 4) < 0.6);
  }
});

test('module reports require explicit boolean results and finite measured values', () => {
  const valid = { checks: [{ label: 'measured', value: '1e-9', pass: true }] };
  assert.equal(validateModuleReport(valid), valid);
  assert.doesNotThrow(() => validateModuleReport({ checks: [{ label: 'measured', value: 2, pass: false }] }));
  for (const pass of [undefined, null, 'true', 1]) assert.throws(() => validateModuleReport({ checks: [{ label: 'measured', value: 1, pass }] }), /explicitly/);
  for (const value of [NaN, Infinity, null, '', 'NaN', 'rel. error = NaN', '∞', '—']) assert.throws(() => validateModuleReport({ checks: [{ label: 'measured', value, pass: true }] }), /finite/);
});

test('gravity validation rejects NaN and failed quantitative benchmarks; settings changes invalidate evidence', () => {
  const report = verifyScenario({ bodies: THREE_BODY_PRESETS['figure-eight'].bodies, duration: 2 });
  assert.deepEqual(gravityReportProblems(report), []);
  const settings = { ...report.settings, h: report.conservation.h };
  assert.equal(gravityReportMatches(report, settings), true);
  for (const changed of [{ integrator: 'leapfrog' }, { stepsPerCrossing: 100 }, { h: settings.h / 2 }, { G: 2 }, { softening: 0.01 }]) {
    assert.equal(gravityReportMatches(report, { ...settings, ...changed }), false);
  }
  for (const key of ['dE', 'dP', 'dL']) {
    const bad = structuredClone(report); bad.conservation[key] = NaN;
    assert.ok(gravityReportProblems(bad).length > 0);
  }
  for (const [index, key] of [[0, 'positionError'], [1, 'returnError']]) {
    const bad = structuredClone(report); bad.benchmarks[index][key] = 1;
    assert.ok(gravityReportProblems(bad).some(p => p.includes(bad.benchmarks[index].name)));
  }
  assert.ok(report.conservation.dE >= report.conservation.final.dE);
});

test('CLI rejects malformed and failed verification, including imported report files', t => {
  const dir = mkdtempSync(join(tmpdir(), 'psim-verification-'));
  t.after(() => rmSync(dir, { recursive: true, force: true }));
  const moduleFile = join(dir, 'physics.mjs');
  const run = (...args) => spawnSync(process.execPath, [join(ROOT, 'scripts/psim.mjs'), ...args], { encoding: 'utf8', timeout: 10000 });
  for (const check of [{ label: 'missing result', value: 'NaN' }, { label: 'null result', value: 1, pass: null }, { label: 'failed', value: 1, pass: false }]) {
    writeFileSync(moduleFile, `export function verify() { return ${JSON.stringify({ checks: [check] })}; }`);
    const r = run('verify-module', moduleFile);
    assert.equal(r.error, undefined);
    assert.notEqual(r.status, 0, r.stdout);
    if (check.pass === false) assert.equal(r.status, 2);
  }
  const page = join(dir, 'page.html'), report = join(dir, 'report.json'), out = join(dir, 'out.html');
  writeFileSync(page, '<script type="application/json" id="psim-verification">null</script>');
  writeFileSync(report, JSON.stringify({ checks: [{ label: 'missing result', value: 1 }] }));
  assert.notEqual(run('build', page, '--verification', report, '--out', out).status, 0);
  assert.equal(existsSync(out), false);
});

test('browser stress check times out a frozen input handler and avoids further renderer commands', async () => {
  let calls = 0;
  const browser = {
    consoleErrors: [], goto: async () => {},
    eval: async () => { calls++; return calls === 1 ? 1 : new Promise(() => {}); },
  };
  let timer;
  const deadline = new Promise((_resolve, reject) => { timer = setTimeout(() => reject(new Error('stress check hung')), 1000); });
  try {
    const [r] = await Promise.race([stressInputs(browser, '/memory.html', { timeoutMs: 20, initialWaitMs: 0, settleMs: 0 }), deadline]);
    assert.equal(r.ok, false);
    assert.match(r.summary, /froze during input change/);
    assert.equal(calls, 2);
  } finally { clearTimeout(timer); }
});

test('browser stress check treats evaluation errors as failures', async () => {
  let calls = 0;
  const browser = { consoleErrors: [], goto: async () => {}, eval: async () => { if (++calls === 1) return 1; throw new Error('handler failed'); } };
  const [r] = await stressInputs(browser, '/memory.html', { timeoutMs: 20, initialWaitMs: 0, settleMs: 0 });
  assert.equal(r.ok, false);
  assert.match(r.summary, /handler failed/);
});

test('blackbody axes respond to a scale-only edit and contain cold and hot spectral peaks', async () => {
  const physics = await import('../examples/blackbody/physics.mjs');
  const text = readFileSync(join(ROOT, 'examples/blackbody/page.html'), 'utf8');
  const code = text.slice(text.indexOf('function bench('), text.indexOf('\ncreateInstrument({'));
  const plots = [];
  const createPlot = opts => {
    const plot = { el: {}, x: { ...opts.x }, y: { ...opts.y }, setX(x) { Object.assign(this.x, x); return this; }, setY(y) { Object.assign(this.y, y); return this; } };
    for (const key of ['fn', 'band', 'marker', 'vline', 'cursor', 'invalidate', 'draw']) plot[key] = () => plot;
    plots.push(plot); return plot;
  };
  const element = () => ({ append() {}, prepend() {} });
  const args = { eqs: element(), figs: element(), createPlot, createEquation: () => ({ el: {}, update() {} }), createDerivation: element, M, setParam() {} };
  const ctx = { ...physics, args, result: null, sci: String, document: { createElement: element } };
  vm.runInNewContext(`${code}\nresult = bench(args);`, ctx);
  ctx.result.update(physics.createSim({ T: 5772, scale: 'log' }));
  assert.equal(plots[0].x.log, true);
  ctx.result.update(physics.createSim({ T: 5772, scale: 'lin' }));
  assert.equal(plots[0].x.log, false);
  for (const T of [1, 2.7255, 100000]) {for (const scale of ['log', 'lin']) {
    const sim = physics.createSim({ T, scale });
    ctx.result.update(sim);
    assert.ok(plots[0].x.min < sim.lmax * 1e9 && plots[0].x.max > sim.lmax * 1e9);
  }}
});
