import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  rng, fft, fftFreq, rk4Step, dopri5, gaussLegendre, adaptiveSimpson, integrateToInfinity, brent,
  lgamma, gamma, erf, besselI, laguerre, legendreP, realYlm, solveTridiagonal, linearFit, meanWithError,
} from '../core/numerics.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);

test('rng is deterministic and roughly uniform / normal', () => {
  const a = rng(42), b = rng(42);
  for (let i = 0; i < 10; i++) assert.equal(a(), b());
  const r = rng(7); let s = 0, s2 = 0; const N = 200000;
  for (let i = 0; i < N; i++) { const x = r.normal(); s += x; s2 += x * x; }
  close(s / N, 0, 0.01, 'normal mean'); close(s2 / N, 1, 0.01, 'normal variance');
});

test('FFT matches a direct DFT and inverts exactly', () => {
  const n = 64, r = rng(3);
  const re = Float64Array.from({ length: n }, () => r()), im = Float64Array.from({ length: n }, () => r());
  const R = re.slice(), I = im.slice();
  fft(R, I);
  for (const k of [0, 1, 7, 31, 63]) {
    let dr = 0, di = 0;
    for (let j = 0; j < n; j++) { const a = -2 * Math.PI * k * j / n; dr += re[j] * Math.cos(a) - im[j] * Math.sin(a); di += re[j] * Math.sin(a) + im[j] * Math.cos(a); }
    close(R[k], dr, 1e-12); close(I[k], di, 1e-12);
  }
  fft(R, I, true);
  for (let j = 0; j < n; j++) { close(R[j], re[j], 1e-14); close(I[j], im[j], 1e-14); }
  const k = fftFreq(8, 0.5);
  close(k[1], 2 * Math.PI / 4, 1e-15); close(k[4], -Math.PI / 0.5, 1e-15);
});

test('RK4 converges at 4th order on y\' = y', () => {
  const err = h => { let y = [1]; for (let t = 0; t < 1 - 1e-12; t += h) y = rk4Step((_t, v) => [v[0]], t, y, h); return Math.abs(y[0] - Math.E); };
  close(Math.log2(err(0.1) / err(0.05)), 4, 0.1);
});

test('Dormand–Prince solves the harmonic oscillator to tolerance', () => {
  const r = dopri5((t, y) => [y[1], -y[0]], 0, [1, 0], 20 * Math.PI, { rtol: 1e-11, atol: 1e-13 });
  close(r.y[0], 1, 1e-8); close(r.y[1], 0, 1e-8);
  assert.ok(r.steps > 50 && r.steps < 5000);
});

test('quadrature: Gauss–Legendre, adaptive Simpson, semi-infinite', () => {
  close(gaussLegendre(x => x ** 9, -1, 2, 5), (2 ** 10 - 1) / 10, 1e-12, 'GL exact for degree 9 with n=5');
  close(gaussLegendre(Math.sin, 0, Math.PI, 20), 2, 1e-14);
  close(adaptiveSimpson(x => Math.exp(-x * x), 0, 3, 1e-13), Math.sqrt(Math.PI) / 2 * erf(3), 1e-12);
  // ∫0^∞ x³/(eˣ−1) dx = π⁴/15 (the Stefan–Boltzmann integral)
  close(integrateToInfinity(x => x === 0 ? 0 : x ** 3 / Math.expm1(x)), Math.PI ** 4 / 15, 1e-10);
});

test('Brent finds the Wien displacement root', () => {
  const root = brent(x => x - 5 * (1 - Math.exp(-x)), 1, 10);
  close(root, 4.965114231744276, 1e-13);
});

test('special functions match tabulated values', () => {
  close(gamma(5), 24, 1e-12); close(gamma(0.5), Math.sqrt(Math.PI), 1e-14);
  close(lgamma(100), 359.13420536957540, 1e-10);
  close(erf(0.5), 0.5204998778130465, 1e-15); close(erf(3), 0.9999779095030014, 1e-15); close(erf(-1.2), -0.9103139782296353, 1e-15);
  close(besselI(1, 2), 1.590636854637329, 1e-14); close(besselI(2, 2), 0.688948447698738, 1e-14);
  close(laguerre(3, 1, 0.7), (-(0.7 ** 3) + 12 * 0.7 ** 2 - 36 * 0.7 + 24) / 6, 1e-14);
  close(legendreP(3, 0, 0.3), 0.5 * (5 * 0.027 - 3 * 0.3), 1e-15);
  close(legendreP(2, 1, 0.5), -3 * 0.5 * Math.sqrt(0.75), 1e-15);
});

test('real spherical harmonics are orthonormal', () => {
  const lm = [[0, 0], [1, -1], [1, 0], [1, 1], [2, -2], [2, 1], [3, 2]];
  const inner = (a, b) => gaussLegendre(th => Math.sin(th) * gaussLegendre(ph => realYlm(a[0], a[1], th, ph) * realYlm(b[0], b[1], th, ph), 0, 2 * Math.PI, 24), 0, Math.PI, 24);
  for (const a of lm) for (const b of lm) close(inner(a, b), a === b ? 1 : 0, 1e-12, `<${a}|${b}>`);
  assert.ok(realYlm(1, 1, Math.PI / 2, 0) > 0, 'p_x positive along +x');
});

test('tridiagonal solver, linear fit, blocking error', () => {
  const x = solveTridiagonal([0, 1, 1], [4, 4, 4], [1, 1, 0], [5, 6, 5]);
  close(x[0], 1, 1e-15); close(x[1], 1, 1e-15); close(x[2], 1, 1e-15);
  const f = linearFit([0, 1, 2, 3], [1, 3, 5, 7]);
  close(f.a, 1, 1e-14); close(f.b, 2, 1e-14); close(f.r2, 1, 1e-14);
  const r = rng(1), xs = Array.from({ length: 4096 }, () => r());
  const m = meanWithError(xs);
  close(m.error, Math.sqrt(1 / 12 / 4096), 0.002);
});

test('complete elliptic integral K(m) (AGM) matches tabulated values', async () => {
  const { ellipticK } = await import('../core/numerics.mjs');
  close(ellipticK(0), Math.PI / 2, 1e-15);
  close(ellipticK(0.5), 1.854074677301372, 1e-14);
  close(ellipticK(0.9), 2.578092113348173, 1e-14);
  close(ellipticK(-1), 1.311028777146060, 1e-14);
  close(ellipticK(0.3), gaussLegendre(t => 1 / Math.sqrt(1 - 0.3 * Math.sin(t) ** 2), 0, Math.PI / 2, 32), 1e-14);
  assert.equal(ellipticK(1), Infinity);
});
