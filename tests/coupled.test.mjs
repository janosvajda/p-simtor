import { test } from 'node:test';
import assert from 'node:assert/strict';
import { coupledSchrodinger2, schrodingerND } from '../core/fields.mjs';

const N = 128, dx = 0.25, xs = Array.from({ length: N }, (_, i) => (i - N / 2) * dx);
const packet = (x0, k0, s = 1) => { const re = new Float64Array(N), im = new Float64Array(N); xs.forEach((x, i) => { const g = Math.exp(-((x - x0) ** 2) / (4 * s * s)) / Math.pow(2 * Math.PI * s * s, 0.25); re[i] = g * Math.cos(k0 * x); im[i] = g * Math.sin(k0 * x); }); return { re, im }; };

test('coupledSchrodinger2: constant coupling gives Rabi oscillation P2 = sin²(c t/ħ); norm exact', () => {
  const c = 0.3, zero = new Float64Array(N), V12 = new Float64Array(N).fill(c);
  const s = coupledSchrodinger2({ psi0: [packet(0, 0), { re: zero, im: zero }], V11: zero, V22: zero, V12, shape: [N], spacing: dx, dt: 0.05 });
  for (const T of [1, 2.5, 4]) { s.advanceTo(T); const [, P2] = s.populations(); assert.ok(Math.abs(P2 - Math.sin(c * s.t) ** 2) < 1e-10, `t = ${s.t}: ${P2}`); }
  assert.ok(Math.abs(s.norm() - 1) < 1e-12);
});

test('coupledSchrodinger2 without coupling equals schrodingerND on each surface; per-axis mass and energy()', () => {
  const zero = new Float64Array(N), V1 = Float64Array.from(xs, x => 0.5 * x * x), V2 = Float64Array.from(xs, x => 0.2 * x * x + 1);
  const p0 = packet(-1, 1.5);
  const s = coupledSchrodinger2({ psi0: [p0, { re: zero, im: zero }], V11: V1, V22: V2, V12: zero, shape: [N], spacing: dx, dt: 0.01, mass: 1.7 });
  const r = schrodingerND({ re0: p0.re, im0: p0.im, V: V1, shape: [N], spacing: dx, dt: 0.01, mass: [1.7] });
  const E0 = r.energy();
  s.advanceTo(2); r.advanceTo(2);
  let d = 0; for (let i = 0; i < N; i++) d = Math.max(d, Math.hypot(s.psi[0].re[i] - r.re[i], s.psi[0].im[i] - r.im[i]));
  assert.ok(d < 1e-12, `${d}`);
  assert.ok(Math.abs(s.energy() - r.energy()) < 1e-12);
  assert.ok(Math.abs(r.energy() - E0) < 1e-3 * Math.abs(E0), 'energy bounded');
});
