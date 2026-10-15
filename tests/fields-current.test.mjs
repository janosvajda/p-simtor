import { test } from 'node:test';
import assert from 'node:assert/strict';
import { probabilityCurrent, schrodingerND } from '../core/fields.mjs';

test('probabilityCurrent: a plane-wave packet carries j = (ħk/m)|ψ|²; flux out of a box balances the norm lost to absorbers', () => {
  const n = 128, dx = 0.25, k = 3, re = new Float64Array(n * n), im = new Float64Array(n * n);
  for (let j = 0; j < n; j++) {for (let i = 0; i < n; i++) {
    const x = (i - 64) * dx, y = (j - 64) * dx, g = Math.exp(-(x * x + y * y) / 8);
    re[i + n * j] = g * Math.cos(k * x); im[i + n * j] = g * Math.sin(k * x);
  }}
  const [jx, jy] = probabilityCurrent(re, im, { shape: [n, n], spacing: dx });
  const c = 64 + 64 * n, rho = re[c] ** 2 + im[c] ** 2;
  assert.ok(Math.abs(jx[c] / rho - k) < 1e-6, `${jx[c] / rho}`);
  assert.ok(Math.abs(jy[c]) < 1e-9);
  // time-integrated flux through a line beyond which the absorber sits equals the absorbed probability
  const s = schrodingerND({ re0: re, im0: im, shape: [n, n], spacing: dx, dt: 0.005, absorb: { width: 16, strength: 20 } });
  const line = 100, steps = 1800, norm0 = s.norm();
  let flux = 0;
  const lineFlux = () => { const [fx] = probabilityCurrent(s.re, s.im, { shape: [n, n], spacing: dx }); let f = 0; for (let jj = 0; jj < n; jj++) f += fx[line + n * jj] * dx; return f; };
  let prev = lineFlux();
  for (let t = 0; t < steps; t++) { s.step(); const cur = lineFlux(); flux += 0.5 * (prev + cur) * 0.005; prev = cur; }
  const lost = s.absorbed();
  assert.ok(lost > 0.5 * norm0, `absorbed ${lost} of ${norm0}`);
  assert.ok(Math.abs(flux / lost - 1) < 0.02, `flux ${flux} vs absorbed ${lost}`);
});
