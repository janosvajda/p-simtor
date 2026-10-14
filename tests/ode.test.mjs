import { test } from 'node:test';
import assert from 'node:assert/strict';
import { odeStepper, mechanicsStepper, refinementOrder } from '../core/ode.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const oscillator = (_t, [x, v]) => [v, -x];

test('odeStepper RK4: exact exponential decay and 4th-order convergence', () => {
  const run = h => { const s = odeStepper({ f: (_t, y) => [-y[0]], y0: [1], method: 'rk4', h }); s.advanceTo(1 - 1e-12); return s.y; };
  close(run(0.01)[0], Math.exp(-1), 1e-10);
  close(refinementOrder(run, 0.1).order, 4, 0.1, 'order');
});

test('odeStepper Dormand–Prince: harmonic oscillator to tolerance over 10 periods', () => {
  const s = odeStepper({ f: oscillator, y0: [1, 0], h: 0.1, rtol: 1e-11, atol: 1e-13 });
  s.advanceTo(20 * Math.PI);
  const out = s.sampleAt(20 * Math.PI);                   // exact time between steps
  close(out[0], 1, 1e-8); close(out[1], 0, 1e-8);
  assert.ok(s.steps > 20 && s.steps < 5000);
});

test('odeStepper dense output is 4th order between steps and frame-independent', () => {
  const s = odeStepper({ f: oscillator, y0: [1, 0], method: 'rk4', h: 0.05 });
  s.advanceTo(1.03);
  close(s.sampleAt(1.03)[0], Math.cos(1.03), 1e-6);
  // advancing in irregular chunks gives exactly the same states as one call
  const a = odeStepper({ f: oscillator, y0: [1, 0], h: 0.1 }), b = odeStepper({ f: oscillator, y0: [1, 0], h: 0.1 });
  a.advanceTo(5);
  for (let t = 0; t < 5; t += 0.0137 + (t % 0.05)) b.advanceTo(t);
  b.advanceTo(5);
  assert.deepEqual(Array.from(a.y), Array.from(b.y));
  assert.equal(a.t, b.t);
});

test('odeStepper rejects bad input and stops on non-finite derivatives', () => {
  assert.throws(() => odeStepper({ f: oscillator, y0: [NaN], h: 0.1 }), /finite/);
  assert.throws(() => odeStepper({ f: oscillator, y0: [1, 0], h: 0 }), /step/);
  const s = odeStepper({ f: (t, y) => [t > 0.5 ? NaN : y[0]], y0: [1], method: 'rk4', h: 0.1 });
  assert.throws(() => s.advanceTo(1), /non-finite/);
});

test('mechanicsStepper: symplectic energy bound, orders 2 and 4, time reversibility', () => {
  const accel = (_t, q, a) => { a[0] = -q[0]; };
  // energy error stays bounded (no drift) over 1000 periods
  const s = mechanicsStepper({ accel, q0: [1], v0: [0], h: 0.05, method: 'leapfrog' });
  let worst = 0;
  for (let i = 0; i < 1000 * 126; i++) { s.step(); worst = Math.max(worst, Math.abs(0.5 * (s.q[0] ** 2 + s.v[0] ** 2) - 0.5)); }
  assert.ok(worst < 1e-3, `${worst}`);
  for (const [method, order] of [['leapfrog', 2], ['yoshida4', 4]]) {
    const run = h => { const m = mechanicsStepper({ accel, q0: [1], v0: [0], h, method }); for (let i = 0; i < Math.round(2 / h); i++) m.step(); return m.q; };
    close(refinementOrder(run, 0.1).order, order, 0.15, method);
  }
  const r = mechanicsStepper({ accel, q0: [0.3], v0: [1.1], h: 0.07 });
  for (let i = 0; i < 300; i++) r.step();
  r.h = -0.07;
  for (let i = 0; i < 300; i++) r.step();
  close(r.q[0], 0.3, 1e-12); close(r.v[0], 1.1, 1e-12);
});

test('mechanicsStepper sampleAt interpolates positions to O(h⁴)', () => {
  const accel = (_t, q, a) => { a[0] = -q[0]; };
  const err = h => {
    const m = mechanicsStepper({ accel, q0: [1], v0: [0], h, method: 'yoshida4' });
    m.step();
    return Math.abs(m.sampleAt(h / 2)[0] - Math.cos(h / 2));
  };
  // asymptotic range: at larger h the Hermite (h⁴) and step (h⁵) errors partly cancel
  assert.ok(Math.log2(err(0.1) / err(0.05)) > 3.5, `${err(0.1)} ${err(0.05)}`);
});

test('lyapunovExponent: exact for a linear system, Lorenz λ₁ ≈ 0.906, regular motion → 0', async () => {
  const { lyapunovExponent } = await import('../core/ode.mjs');
  // y' = A y, A = [[0.3, 1], [0, −0.5]]: largest exponent = 0.3 exactly
  const A = (y) => [0.3 * y[0] + y[1], -0.5 * y[1]];
  for (const jvp of [(t, y, w) => A(w), undefined]) {
    const { lambda } = lyapunovExponent({ f: (t, y) => A(y), jvp, y0: [1, 1], duration: 200, interval: 2 });
    assert.ok(Math.abs(lambda - 0.3) < 0.01, `linear ${jvp ? 'variational' : 'twin'}: ${lambda}`);
  }
  // Lorenz σ = 10, ρ = 28, β = 8/3: λ₁ = 0.9056 (Sprott, Chaos and Time-Series Analysis, 2003)
  const lorenz = (t, y) => [10 * (y[1] - y[0]), y[0] * (28 - y[2]) - y[1], y[0] * y[1] - (8 / 3) * y[2]];
  const lorenzJ = (t, y, w) => [10 * (w[1] - w[0]), (28 - y[2]) * w[0] - w[1] - y[0] * w[2], y[1] * w[0] + y[0] * w[1] - (8 / 3) * w[2]];
  const v = lyapunovExponent({ f: lorenz, jvp: lorenzJ, y0: [1, 1, 20], duration: 1000, transient: 20, interval: 0.5, rtol: 1e-9 });
  const b = lyapunovExponent({ f: lorenz, y0: [1, 1, 20], duration: 1000, transient: 20, interval: 0.5, rtol: 1e-9 });
  assert.ok(Math.abs(v.lambda - 0.9056) < 0.03, `Lorenz variational ${v.lambda}`);
  assert.ok(Math.abs(b.lambda - v.lambda) < 0.03, `twin-trajectory agrees: ${b.lambda}`);
  // harmonic oscillator: perturbations grow at most linearly, so λ(t) → 0
  const osc = lyapunovExponent({ f: (t, y) => [y[1], -y[0]], y0: [1, 0], duration: 500, interval: 1 });
  assert.ok(Math.abs(osc.lambda) < 0.01, `regular ${osc.lambda}`);
});
