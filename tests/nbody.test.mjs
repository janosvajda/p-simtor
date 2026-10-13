import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  createSystem, createStepper, integrateTo, totalEnergy, toBarycentric, hermite, stateDistance, INTEGRATORS,
} from '../core/nbody.mjs';
import { THREE_BODY_PRESETS } from '../core/presets.mjs';
import {
  convergenceStudy, determinismCheck, frameIndependenceCheck,
  keplerBenchmark, figureEightBenchmark, pythagoreanBenchmark,
} from '../core/nbody-verify.mjs';
import { solveKepler, keplerRelative } from '../core/kepler.mjs';

test('Kepler solver satisfies E − e sin E = M', () => {
  for (const e of [0, 0.3, 0.9, 0.99]) {for (const M of [0.1, 1, 3, 6]) {
    const E = solveKepler(M, e);
    assert.ok(Math.abs(E - e * Math.sin(E) - M) < 1e-12, `e=${e} M=${M}`);
  }}
});

test('analytic Kepler orbit conserves energy', () => {
  const mu = 1.5, a = 1, e = 0.7;
  const eps = t => { const s = keplerRelative(a, e, mu, t); return 0.5 * (s.vx ** 2 + s.vy ** 2) - mu / Math.hypot(s.x, s.y); };
  assert.ok(Math.abs(eps(0.3) - (-mu / (2 * a))) < 1e-12);
  assert.ok(Math.abs(eps(2.9) - (-mu / (2 * a))) < 1e-12);
});

test('every preset conserves energy, momentum and angular momentum with logh4', () => {
  for (const [k, p] of Object.entries(THREE_BODY_PRESETS)) {
    const r = integrateTo(p.bodies, { integrator: 'logh4' }, p.duration);
    assert.ok(r.dE < 1e-9, `${k}: dE ${r.dE}`);
    assert.ok(r.dP < 1e-12, `${k}: dP ${r.dP}`);
    assert.ok(r.dL < 1e-12, `${k}: dL ${r.dL}`);
  }
});

test('observed convergence order matches theory for all integrators', () => {
  const p = THREE_BODY_PRESETS['figure-eight'];
  for (const [name, spec] of Object.entries(INTEGRATORS)) {
    const r = convergenceStudy(p.bodies, { integrator: name }, p.period, { samples: 40 });
    assert.ok(r.observedOrder != null, `${name}: order unresolved`);
    assert.ok(Math.abs(r.observedOrder - spec.order) < 0.3, `${name}: observed ${r.observedOrder}, expected ${spec.order}`);
  }
});

test('eccentric Kepler orbit matches the analytic solution', () => {
  const r = keplerBenchmark({ e: 0.9, orbits: 10 });
  assert.ok(r.positionError < 1e-3, `position error ${r.positionError}`);
  assert.ok(r.dE < 1e-10);
});

test('figure-eight returns to its initial state after one period', () => {
  const r = figureEightBenchmark();
  assert.ok(r.returnError < 1e-6, `return error ${r.returnError}`);
});

test('Pythagorean problem reproduces the published outcome', () => {
  const r = pythagoreanBenchmark();
  assert.equal(r.matchesLiterature, true, r.outcome);
  assert.ok(r.dE < 1e-9);
});

test('fixed-step integrators fail the Pythagorean close encounters (why logh4 is the default)', () => {
  const r = integrateTo(THREE_BODY_PRESETS.pythagorean.bodies, { integrator: 'yoshida4' }, 80);
  assert.ok(r.dE > 1e-3, `expected large energy error, got ${r.dE}`);
});

test('integration is bit-for-bit deterministic', () => {
  assert.equal(determinismCheck(THREE_BODY_PRESETS.pythagorean.bodies, { integrator: 'logh4' }, 20).identical, true);
});

test('state at time t does not depend on the frame schedule', () => {
  for (const k of ['pythagorean', 'figure-eight']) {
    const p = THREE_BODY_PRESETS[k];
    assert.equal(frameIndependenceCheck(p.bodies, { integrator: 'logh4' }, 10).identical, true, k);
  }
});

test('Hermite interpolation is 4th-order accurate between steps', () => {
  // Circular orbit: compare interpolated position at mid-step with the exact one.
  const err = dt => {
    const pos = t => [Math.cos(t), Math.sin(t), 0];
    const vel = t => [-Math.sin(t), Math.cos(t), 0];
    const out = new Float64Array(3);
    hermite(0, Float64Array.from(pos(0)), Float64Array.from(vel(0)), dt, Float64Array.from(pos(dt)), Float64Array.from(vel(dt)), dt / 2, out);
    const ex = pos(dt / 2);
    return Math.hypot(out[0] - ex[0], out[1] - ex[1]);
  };
  const ratio = err(0.1) / err(0.05);
  assert.ok(ratio > 14 && ratio < 18, `error ratio ${ratio} (expect ≈16)`);
});

test('changing a mass changes the computed evolution', () => {
  const bodies = THREE_BODY_PRESETS.pythagorean.bodies;
  const heavier = bodies.map((b, i) => ({ ...b, m: i === 0 ? b.m * 1.5 : b.m }));
  const a = integrateTo(bodies, { integrator: 'logh4' }, 2);
  const b = integrateTo(heavier, { integrator: 'logh4' }, 2);
  assert.ok(stateDistance(a.x, b.x) > 1e-2);
});

test('barycentric shift zeroes total momentum and centre of mass', () => {
  const s = toBarycentric(createSystem([{ m: 1, x: [1, 2, 3], v: [1, 0, 0] }, { m: 3, x: [0, 0, 0], v: [0, 1, 0] }]));
  let px = 0, py = 0, cx = 0;
  for (let i = 0; i < 2; i++) { px += s.m[i] * s.v[3 * i]; py += s.m[i] * s.v[3 * i + 1]; cx += s.m[i] * s.x[3 * i]; }
  assert.ok(Math.abs(px) < 1e-15 && Math.abs(py) < 1e-15 && Math.abs(cx) < 1e-15);
});

test('LogH physical step shrinks during close encounters', () => {
  const sys = createSystem(THREE_BODY_PRESETS.pythagorean.bodies);
  toBarycentric(sys);
  const st = createStepper(sys, { integrator: 'logh4' });
  let dtMin = Infinity, dtMax = 0;
  while (sys.t < 30) {
    const t0 = sys.t;
    st.step();
    const dt = sys.t - t0;
    dtMin = Math.min(dtMin, dt); dtMax = Math.max(dtMax, dt);
  }
  assert.ok(dtMax / dtMin > 100, `dt range ${dtMin}..${dtMax}`);
  assert.ok(Number.isFinite(totalEnergy(sys)));
});
