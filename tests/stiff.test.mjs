import { test } from 'node:test';
import assert from 'node:assert/strict';
import { stiffStepper, stiffnessRatio, odeAuto } from '../core/stiff.mjs';
import { refinementOrder } from '../core/ode.mjs';
import { solveLinear } from '../core/linalg.mjs';

const close = (a, b, tol, msg = '') => assert.ok(Math.abs(a - b) <= tol, `${msg} ${a} vs ${b} (tol ${tol})`);
const rel = (a, b, tol, msg = '') => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg} ${a} vs ${b} (rel tol ${tol})`);
const METHODS = ['radau5', 'rosenbrock', 'bdf'];

// Robertson chemical kinetics (Hairer & Wanner II, §IV.1 and §IV.10)
const rob = (_t, y) => [-0.04 * y[0] + 1e4 * y[1] * y[2], 0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1], 3e7 * y[1] * y[1]];
const robJ = (_t, y) => [[-0.04, 1e4 * y[2], 1e4 * y[1]], [0.04, -1e4 * y[2] - 6e7 * y[1], -1e4 * y[1]], [0, 6e7 * y[1], 0]];

test('Robertson kinetics: Hairer–Wanner reference values at t = 40 and t = 1e11, mass conserved', () => {
  const ref40 = [0.7158270687e0, 0.9185534764e-5, 0.2841637457e0];               // Hairer & Wanner II, Table IV.10 / ROBER
  const ref11 = [0.2083340149701255e-7, 0.8333360770334713e-13, 0.9999999791665050];  // IVP test set (Hairer–Wanner lineage), t = 1e11
  for (const method of METHODS) {
    const rtol = method === 'rosenbrock' ? 1e-8 : 1e-10;
    const s = stiffStepper({ f: rob, jac: robJ, y0: [1, 0, 0], method, rtol, atol: 1e-16 });
    s.advanceTo(40, 1e6, { exact: true });
    assert.equal(s.t, 40);
    ref40.forEach((v, i) => rel(s.y[i], v, 1e-7, `${method} y${i + 1}(40)`));
    close(s.y[0] + s.y[1] + s.y[2], 1, 1e-13, `${method} mass at 40`);
    s.advanceTo(1e11, 1e6, { exact: true });
    rel(s.y[0], ref11[0], 2e-7, `${method} y1(1e11)`);
    rel(s.y[1], ref11[1], 2e-7, `${method} y2(1e11)`);
    rel(s.y[2], ref11[2], 1e-12, `${method} y3(1e11)`);
    close(s.y[0] + s.y[1] + s.y[2], 1, 1e-13, `${method} mass at 1e11`);
    assert.ok(s.steps < 1e5, `${method} steps ${s.steps}`);
  }
});

test('Robertson at t = 4e5: the three methods agree (y1 = 4.9382745e-3)', () => {
  const y = METHODS.map(method => {
    const s = stiffStepper({ f: rob, jac: robJ, y0: [1, 0, 0], method, rtol: method === 'rosenbrock' ? 1e-9 : 1e-10, atol: 1e-16 });
    s.advanceTo(4e5, 1e6, { exact: true });
    return s.y;
  });
  for (const v of y) {
    rel(v[0], 4.938274521e-3, 5e-9, 'y1');
    rel(v[1], 1.984994088e-8, 5e-9, 'y2');
    rel(v[2], 0.9950617056, 1e-10, 'y3');
  }
});

test('Van der Pol μ = 1000: relaxation period matches the asymptotic (3 − 2 ln 2)μ + 7.0143 μ^{−1/3}', () => {
  const mu = 1000;
  const f = (_t, y) => [y[1], mu * (1 - y[0] * y[0]) * y[1] - y[0]];
  const jac = (_t, y) => [[0, 1], [-2 * mu * y[0] * y[1] - 1, mu * (1 - y[0] * y[0])]];
  const asymptotic = (3 - 2 * Math.log(2)) * mu + 7.0143 * Math.cbrt(1 / mu);   // Dorodnitsyn (1947); Hairer & Wanner II, §IV.1
  for (const method of ['radau5', 'bdf']) {
    const s = stiffStepper({ f, jac, y0: [2, 0], method, rtol: 1e-8, atol: 1e-10, events: [{ g: (_t, y) => y[0], direction: -1 }] });
    s.advanceTo(2.6 * asymptotic);
    assert.equal(s.events.length, 3, method);
    const T1 = s.events[1].t - s.events[0].t, T2 = s.events[2].t - s.events[1].t;
    close(T1, T2, 1e-3, `${method} periodic`);
    close(T1, asymptotic, 0.02, `${method} period`);
    rel(T1, (3 - 2 * Math.log(2)) * mu, 1e-3, `${method} leading term`);
    // an explicit method would need ~μ·T/3 ≈ 5e5 stability-limited steps; implicit needs a few thousand
    assert.ok(s.steps < 1e4, `${method} steps ${s.steps}`);
  }
});

test('HIRES and Oregonator: Hairer–Wanner test-set reference solutions', () => {
  const hires = (_t, y) => [
    -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007, 1.71 * y[0] - 8.75 * y[1],
    -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4], 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3],
    -1.745 * y[4] + 0.43 * y[5] + 0.43 * y[6], -280 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6],
    280 * y[5] * y[7] - 1.81 * y[6], -280 * y[5] * y[7] + 1.81 * y[6]];
  const href = [0.7371312573325668e-3, 0.1442485726316185e-3, 0.5888729740967575e-4, 0.1175651343283149e-2,
    0.2386356198831331e-2, 0.6238968252742796e-2, 0.2849998395185769e-2, 0.2850001604814231e-2];
  const orego = (_t, y) => [77.27 * (y[1] + y[0] * (1 - 8.375e-6 * y[0] - y[1])), (y[2] - (1 + y[0]) * y[1]) / 77.27, 0.161 * (y[0] - y[2])];
  const oref = [0.1000814870318523e1, 0.1228178521549917e4, 0.1320554942846706e3];
  for (const method of ['radau5', 'bdf']) {
    const s = stiffStepper({ f: hires, y0: [1, 0, 0, 0, 0, 0, 0, 0.0057], method, rtol: 1e-10, atol: 1e-14 });
    s.advanceTo(321.8122, 1e6, { exact: true });
    href.forEach((v, i) => rel(s.y[i], v, 1e-7, `${method} HIRES y${i + 1}`));
    const o = stiffStepper({ f: orego, y0: [1, 2, 3], method, rtol: 1e-10, atol: 1e-10 });
    o.advanceTo(360, 1e6, { exact: true });
    oref.forEach((v, i) => rel(o.y[i], v, 1e-6, `${method} OREGO y${i + 1}`));
  }
  const r = stiffStepper({ f: orego, y0: [1, 2, 3], method: 'rosenbrock', rtol: 1e-7, atol: 1e-8 });
  r.advanceTo(360, 1e6, { exact: true });
  oref.forEach((v, i) => rel(r.y[i], v, 1e-4, `rosenbrock OREGO y${i + 1}`));
});

test('linear test y′ = λy, λ = −1e6: L-stability (one step h = 1 kills the mode) and A-stability', () => {
  for (const method of METHODS) {
    const s = stiffStepper({ f: (_t, y) => [-1e6 * y[0]], jac: () => [[-1e6]], y0: [1], method, adaptive: false, h: 1 });
    s.step();
    assert.ok(Math.abs(s.y[0]) < 1e-5, `${method} R(−1e6) = ${s.y[0]}`);
    // oscillatory decaying mode λ = −1 ± 1e3 i with hλ far outside any explicit stability region
    const A = [[-1, -1e3], [1e3, -1]];
    const o = stiffStepper({ f: (_t, y) => [A[0][0] * y[0] + A[0][1] * y[1], A[1][0] * y[0] + A[1][1] * y[1]], jac: () => A, y0: [1, 0], method, adaptive: false, h: 0.5, ...(method === 'bdf' ? { maxOrder: 2 } : {}) });
    let norm = 1;
    for (let i = 0; i < 10; i++) { o.step(); const nn = Math.hypot(o.y[0], o.y[1]); assert.ok(nn <= norm * (1 + 1e-12), `${method} amplification`); norm = nn; }
  }
});

test('stiff tracking y′ = −1e6 (y − cos t): few adaptive steps, solution on the slow manifold', () => {
  for (const method of METHODS) {
    const s = stiffStepper({ f: (t, y) => [-1e6 * (y[0] - Math.cos(t))], y0: [0], method, rtol: 1e-6, atol: 1e-9 });
    s.advanceTo(10, 1e5, { exact: true });
    close(s.y[0], Math.cos(10) - 1e-6 * Math.sin(10), 1e-5, method);
    assert.ok(s.steps < 2000, `${method} steps ${s.steps}`);
  }
});

test('measured order: Radau IIA 5, Rodas3 3, BDF1 1, BDF2 2 (fixed steps, nonlinear oscillator)', () => {
  const f = (_t, y) => [y[1], -y[0] - 0.5 * y[1] * y[1] * y[0]];
  const jac = (_t, y) => [[0, 1], [-1 - 0.5 * y[1] * y[1], -y[1] * y[0]]];
  for (const [method, expected, extra] of [['radau5', 5, {}], ['rosenbrock', 3, {}], ['bdf', 1, { maxOrder: 1 }], ['bdf', 2, { maxOrder: 2 }]]) {
    const run = h => { const s = stiffStepper({ f, jac, y0: [0.3, 1], method, adaptive: false, h, rtol: 1e-13, atol: 1e-15, ...extra }); s.advanceTo(1 - 1e-12); return s.y; };
    close(refinementOrder(run, 0.1).order, expected, 0.15, `${method} order`);
  }
});

test('BDF is variable order (reaches 5) and tolerance-proportional', () => {
  const f = (_t, y) => [y[1], -y[0]];
  const err = rtol => {
    const s = stiffStepper({ f, y0: [1, 0], method: 'bdf', rtol, atol: rtol * 1e-3 });
    s.advanceTo(10, 1e6, { exact: true });
    return { e: Math.hypot(s.y[0] - Math.cos(10), s.y[1] + Math.sin(10)), order: s.stats.order };
  };
  const a = err(1e-5), b = err(1e-8);
  assert.equal(b.order, 5);
  assert.ok(a.e / b.e > 100 && a.e / b.e < 1e5, `ratio ${a.e / b.e}`);
  assert.ok(b.e < 1e-5);
});

test('radau5 with a singular mass matrix: index-1 DAE against its exact solution', () => {
  // y1′ = y2 − 2y1, 0 = y2 − y1²  ⇒  y1 = 2/(1 + e^{2t})
  const s = stiffStepper({ f: (_t, y) => [y[1] - 2 * y[0], y[1] - y[0] * y[0]], y0: [1, 1], mass: [[1, 0], [0, 0]], rtol: 1e-10, atol: 1e-12 });
  s.advanceTo(2, 1e5, { exact: true });
  const y1 = 2 / (1 + Math.exp(4));
  close(s.y[0], y1, 1e-10); close(s.y[1], y1 * y1, 1e-10);
  const mid = s.sampleAt(1.99);
  close(mid[0], 2 / (1 + Math.exp(3.98)), 1e-7, 'dense');
});

test('events (odeStepper semantics) and dense output for every method', () => {
  for (const method of METHODS) {
    const s = stiffStepper({ f: (_t, y) => [-y[0]], y0: [1], method, rtol: 1e-10, atol: 1e-12, events: [{ g: (_t, y) => y[0] - 0.75, direction: -1 }, { g: (_t, y) => y[0] - 0.5, terminal: true }] });
    s.advanceTo(5);
    assert.ok(s.stopped, method);
    assert.equal(s.events.length, 2);
    close(s.events[0].t, Math.log(4 / 3), 1e-7, `${method} event 1`);
    close(s.t, Math.LN2, 1e-7, `${method} terminal`);
    close(s.y[0], 0.5, 1e-7);
    const steps = s.steps;
    s.step();
    assert.equal(s.steps, steps, 'no steps after a terminal event');
    // dense output between steps
    const d = stiffStepper({ f: (_t, y) => [-y[0]], y0: [1], method, rtol: 1e-10, atol: 1e-12 });
    d.advanceTo(1);
    const before = d.t;
    d.step();
    const tm = before + 0.37 * (d.t - before);
    close(d.sampleAt(tm)[0], Math.exp(-tm), 1e-7, `${method} dense`);
  }
});

test('finite-difference Jacobian gives the same answer as the analytic one', () => {
  for (const method of METHODS) {
    const a = stiffStepper({ f: rob, jac: robJ, y0: [1, 0, 0], method, rtol: 1e-7, atol: 1e-14 });
    const b = stiffStepper({ f: rob, y0: [1, 0, 0], method, rtol: 1e-7, atol: 1e-14 });
    a.advanceTo(100, 1e6, { exact: true }); b.advanceTo(100, 1e6, { exact: true });
    for (let i = 0; i < 3; i++) close(a.y[i], b.y[i], 1e-5 * Math.abs(a.y[i]) + 1e-14, `${method} y${i + 1}`);
    assert.ok(b.stats.fevals > a.stats.fevals && a.stats.jacobians > 0 && a.stats.decompositions > 0 && a.stats.solves > 0);
  }
});

test('stiffnessRatio: eigenvalues of diagonal, similar and rotation matrices; Robertson becomes stiff', () => {
  const d = stiffnessRatio([[-1, 0], [0, -1e6]]);
  close(d.ratio, 1e6, 1e-6); close(d.spectralRadius, 1e6, 1e-6);
  const rot = stiffnessRatio([[0, 1], [-1, 0]]).eigenvalues;
  close(rot[0].re, 0, 1e-15); close(Math.abs(rot[0].im), 1, 1e-15);
  // S · diag-block(−2 ± 3i, −5, −1e4) · S⁻¹ with a non-orthogonal S: eigenvalues are preserved
  const B = [[-2, 3, 0, 0], [-3, -2, 0, 0], [0, 0, -5, 0], [0, 0, 0, -1e4]];
  const S = [[1, 2, 0, 1], [0, 1, 3, 0], [1, 0, 1, 2], [2, 1, 0, 1]];
  const SB = S.map(row => B[0].map((_, j) => row.reduce((acc, v, k) => acc + v * B[k][j], 0)));
  const inv = [0, 1, 2, 3].map(j => solveLinear(S, [0, 1, 2, 3].map(i => (i === j ? 1 : 0))));     // columns of S⁻¹
  const A = SB.map(row => [0, 1, 2, 3].map(j => row.reduce((acc, v, k) => acc + v * inv[j][k], 0)));
  const ev = stiffnessRatio(A).eigenvalues;
  const expect = [{ re: -1e4, im: 0 }, { re: -5, im: 0 }, { re: -2, im: -3 }, { re: -2, im: 3 }];
  ev.forEach((e, i) => { close(e.re, expect[i].re, 1e-9 * 1e4); close(e.im, expect[i].im, 1e-9 * 1e4); });
  // Robertson: not stiff at t = 0 (eigenvalues −0.04, 0, 0), stiff on the slow manifold
  assert.ok(stiffnessRatio({ f: rob, jac: robJ, t: 0, y: [1, 0, 0] }).fastest < 0.05);
  const s = stiffStepper({ f: rob, jac: robJ, y0: [1, 0, 0] });
  s.advanceTo(1);
  const late = stiffnessRatio({ f: rob, t: s.t, y: s.y });
  assert.ok(late.fastest > 1e3 && late.ratio > 1e3, `fastest ${late.fastest}`);
  assert.ok(late.eigenvalues.some(e => Math.abs(e.re) < 1e-9 * late.spectralRadius), 'conservation law → zero eigenvalue');
});

test('odeAuto recommends radau5 for Robertson and Van der Pol, dopri5 for a harmonic oscillator', () => {
  const a = odeAuto({ f: rob, y0: [1, 0, 0], tEnd: 1e5 });
  assert.equal(a.stiff, true); assert.equal(a.method, 'radau5');
  a.stepper.advanceTo(1e5);
  close(a.stepper.y[0] + a.stepper.y[1] + a.stepper.y[2], 1, 1e-12);
  const v = odeAuto({ f: (_t, y) => [y[1], 1000 * (1 - y[0] * y[0]) * y[1] - y[0]], y0: [2, 0], tEnd: 3000 });
  assert.equal(v.stiff, true);
  const o = odeAuto({ f: (_t, y) => [y[1], -y[0]], y0: [1, 0], tEnd: 20 });
  assert.equal(o.stiff, false); assert.equal(o.method, 'dopri5');
  o.stepper.advanceTo(20);
  close(o.stepper.sampleAt(20)[0], Math.cos(20), 1e-5);
});

test('stiffStepper rejects bad input and stops on non-finite derivatives', () => {
  assert.throws(() => stiffStepper({ f: rob, y0: [1, 0, 0], method: 'euler' }), /unknown method/);
  assert.throws(() => stiffStepper({ f: rob, y0: [1, NaN, 0] }), /finite/);
  assert.throws(() => stiffStepper({ f: rob, y0: [1, 0, 0], method: 'bdf', mass: [[1, 0, 0], [0, 1, 0], [0, 0, 1]] }), /supported by method 'radau5' only/);
  assert.throws(() => stiffStepper({ f: rob, y0: [1, 0, 0], adaptive: false }), /positive step/);
  assert.throws(() => stiffStepper({ f: rob, jac: () => [[1]], y0: [1, 0, 0] }).step(), /rows/);
  for (const method of METHODS) {
    const s = stiffStepper({ f: (t, y) => [t > 0.5 ? NaN : -y[0]], y0: [1], method, h: 0.1 });
    assert.throws(() => s.advanceTo(1), /non-finite|cannot meet|underflow|singular/, method);
  }
  assert.throws(() => stiffnessRatio([[1, 2]]), /square/);
  assert.throws(() => odeAuto({ f: rob, y0: [1, 0, 0] }), /tEnd/);
});
