// Quasi-Monte Carlo checked against published Sobol points, the algebraic properties of the
// Joe–Kuo (2008) direction numbers, exact (t, m, s)-net stratification, exact discrepancies, and
// measured integration convergence (~N⁻¹ or better for RQMC vs N⁻¹ᐟ² for MC) on Genz (1984)
// test integrands with closed-form values.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { sobol, halton, qmcIntegrate, starDiscrepancy, l2StarDiscrepancy, SOBOL_DIRECTIONS, SOBOL_MAX_DIM } from '../core/qmc.mjs';
import { rng, linearFit, erf } from '../core/numerics.mjs';

const close = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol, `${a} vs ${b} (tolerance ${tol})`);

test('Sobol reproduces the published first points (Gray-code order) and the Joe–Kuo table is valid', () => {
  const expected = [[0, 0, 0], [0.5, 0.5, 0.5], [0.75, 0.25, 0.25], [0.25, 0.75, 0.75], [0.375, 0.375, 0.625], [0.875, 0.875, 0.125], [0.625, 0.125, 0.875], [0.125, 0.625, 0.375]];
  const s = sobol(3);
  for (const e of expected) assert.deepEqual(Array.from(s.next()), e);
  assert.deepEqual(Array.from(sobol(3).point(6)), expected[6]);
  assert.deepEqual(Array.from(sobol(3, { skip: 5 }).next()), expected[5]);
  assert.equal(SOBOL_MAX_DIM, 21);
  // each row: x^s + a₁x^{s−1} + … + 1 is primitive over GF(2); m_k odd with m_k < 2^k
  for (const [s0, a, ...m] of SOBOL_DIRECTIONS) {
    assert.equal(m.length, s0);
    m.forEach((mk, k) => { assert.equal(mk % 2, 1); assert.ok(mk < 2 ** (k + 1)); });
    const P = (1 << s0) | (a << 1) | 1, period = (1 << s0) - 1;
    let x = 1, order = 0;
    do { x <<= 1; if (x & (1 << s0)) x ^= P; order++; } while (x !== 1 && order <= period);
    assert.equal(order, period, `polynomial ${P.toString(2)} is not primitive`);
  }
  assert.throws(() => sobol(22), /1…21/);
});

test('every Sobol coordinate is stratified and the first two dimensions form a (0, m, 2)-net, scrambled or not', () => {
  const m = 10, N = 2 ** m;
  for (const scramble of ['none', 'shift', 'owen']) {
    const pts = sobol(21, { scramble, seed: 4 }).points(N);
    for (let j = 0; j < 21; j++) {
      const seen = new Uint8Array(N);
      for (const p of pts) seen[Math.floor(p[j] * N)]++;
      assert.ok(seen.every(c => c === 1), `${scramble}: coordinate ${j} not stratified`);
    }
    for (let a = 0; a <= m; a++) {           // boxes 2^-a × 2^-(m−a) hold exactly one point each
      const count = new Uint16Array(N);
      for (const p of pts) count[Math.floor(p[0] * 2 ** a) * 2 ** (m - a) + Math.floor(p[1] * 2 ** (m - a))]++;
      assert.ok(count.every(c => c === 1), `${scramble}: not a (0,m,2)-net at a = ${a}`);
    }
  }
});

test('Halton radical inverses are exact; digit-permutation scrambling keeps b^k stratification', () => {
  const h = halton(3).points(5);
  const exact = [[0, 0, 0], [1 / 2, 1 / 3, 1 / 5], [1 / 4, 2 / 3, 2 / 5], [3 / 4, 1 / 9, 3 / 5], [1 / 8, 4 / 9, 4 / 5]];
  h.forEach((p, i) => p.forEach((v, j) => close(v, exact[i][j], 1e-15)));
  const pts = halton(3, { scramble: 'permute', seed: 9 }).points(125);
  const cells = new Uint8Array(125);
  for (const p of pts) cells[Math.floor(p[2] * 125)]++;           // base 5, k = 3
  assert.ok(cells.every(c => c === 1));
  const cells3 = new Uint8Array(81);
  for (const p of halton(2, { scramble: 'permute', seed: 9 }).points(81)) cells3[Math.floor(p[1] * 81)]++;
  assert.ok(cells3.every(c => c === 1));
});

test('exact star discrepancy: 1D closed form, single-point formula, and L2 by Warnock vs quadrature', () => {
  const N = 20, centred = Array.from({ length: N }, (_, i) => [(2 * i + 1) / (2 * N)]);
  close(starDiscrepancy(centred), 1 / (2 * N), 1e-15);
  const rand = rng(8), xs = Array.from({ length: 30 }, () => rand()).sort((a, b) => a - b);
  const closed1D = 1 / (2 * xs.length) + Math.max(...xs.map((x, i) => Math.abs(x - (2 * i + 1) / (2 * xs.length))));
  close(starDiscrepancy(xs.map(x => [x])), closed1D, 1e-15);
  for (const a of [[0.3, 0.6], [0.9, 0.2], [0.5, 0.5]]) close(starDiscrepancy([a]), Math.max(a[0], a[1], 1 - a[0] * a[1]), 1e-15);
  // Warnock's formula against midpoint quadrature of the defining integral
  const P = Array.from({ length: 6 }, () => [rand(), rand()]), G = 600;
  let integral = 0;
  for (let i = 0; i < G; i++) {for (let j = 0; j < G; j++) {
    const y = [(i + 0.5) / G, (j + 0.5) / G];
    const frac = P.filter(p => p[0] < y[0] && p[1] < y[1]).length / P.length;
    integral += (frac - y[0] * y[1]) ** 2 / (G * G);
  }}
  close(l2StarDiscrepancy(P), Math.sqrt(integral), 2e-3);
  // Sobol beats random points in both measures
  const sob = sobol(2).points(256), ran = Array.from({ length: 256 }, () => [rand(), rand()]);
  assert.ok(starDiscrepancy(sob) < starDiscrepancy(ran) / 3);
  assert.ok(l2StarDiscrepancy(sob) < l2StarDiscrepancy(ran) / 3);
  assert.throws(() => starDiscrepancy(Array.from({ length: 500 }, () => [rand(), rand(), rand()])), /maxWork/);
});

// Genz test integrands on [0, 1]^5 with closed-form integrals
const A = [1.1, 0.7, 1.6, 0.9, 1.3], U = [0.3, 0.6, 0.2, 0.8, 0.5];
const oscillatory = x => Math.cos(2 * Math.PI * U[0] + A.reduce((s, a, k) => s + a * x[k], 0));
const oscillatoryExact = (() => {
  let re = Math.cos(2 * Math.PI * U[0]), im = Math.sin(2 * Math.PI * U[0]);
  for (const a of A) {           // multiply by (e^{ia} − 1)/(ia) = (sin a + i(1 − cos a))/a
    const fr = Math.sin(a) / a, fi = (1 - Math.cos(a)) / a;
    [re, im] = [re * fr - im * fi, re * fi + im * fr];
  }
  return re;
})();
const productPeak = x => A.reduce((p, a, k) => p / (a ** -2 + (x[k] - U[k]) ** 2), 1);
const productPeakExact = A.reduce((p, a, k) => p * a * (Math.atan(a * (1 - U[k])) + Math.atan(a * U[k])), 1);
const gaussian = x => Math.exp(-A.reduce((s, a, k) => s + (a * (x[k] - U[k])) ** 2, 0));
const gaussianExact = A.reduce((p, a, k) => p * Math.sqrt(Math.PI) / (2 * a) * (erf(a * (1 - U[k])) + erf(a * U[k])), 1);

const measuredRate = (f, exact, sequence) => {
  const ms = [7, 8, 9, 10, 11, 12], logErr = [];
  for (const m of ms) {
    const r = qmcIntegrate(f, 5, { n: 2 ** m, replicates: 16, sequence, seed: m });
    const rmse = Math.sqrt(r.estimates.reduce((a, v) => a + (v - exact) ** 2, 0) / r.replicates);
    logErr.push(Math.log2(rmse));
    assert.ok(Math.abs(r.value - exact) <= 4.5 * r.error, `${sequence} n=2^${m}: ${r.value} vs ${exact} ± ${r.error}`);
  }
  return linearFit(ms, logErr).b;
};

test('randomised QMC converges like N⁻¹ or faster on Genz integrands while MC converges like N⁻¹ᐟ²', () => {
  for (const [name, f, exact] of [['oscillatory', oscillatory, oscillatoryExact], ['product peak', productPeak, productPeakExact], ['gaussian', gaussian, gaussianExact]]) {
    const sobolRate = measuredRate(f, exact, 'sobol'), haltonRate = measuredRate(f, exact, 'halton'), mcRate = measuredRate(f, exact, 'random');
    assert.ok(sobolRate < -1.0, `${name}: Owen-scrambled Sobol slope ${sobolRate}`);
    assert.ok(haltonRate < -0.8, `${name}: scrambled Halton slope ${haltonRate}`);
    assert.ok(mcRate > -0.75 && mcRate < -0.3, `${name}: Monte Carlo slope ${mcRate}`);
  }
});

test('all 21 Sobol dimensions integrate a 21-dimensional product exactly enough, box mapping works', () => {
  const f = x => x.reduce((p, v, k) => p * (1 + (v - 0.5) / (k + 1)), 1);  // exact integral 1
  const q = qmcIntegrate(f, 21, { n: 4096, replicates: 10 }), mc = qmcIntegrate(f, 21, { n: 4096, replicates: 10, sequence: 'random' });
  assert.ok(Math.abs(q.value - 1) <= 4 * q.error);
  assert.ok(q.error < mc.error / 10, `${q.error} vs MC ${mc.error}`);
  const box = qmcIntegrate(x => x[0] * x[1], 2, { n: 1024, lower: [1, -1], upper: [3, 2] });   // ∫∫ xy = 4 · 1.5
  assert.ok(Math.abs(box.value - 6) <= 4 * box.error && box.error < 1e-3);
  assert.throws(() => qmcIntegrate(f, 2, { n: 1000 }), /2\^m/);
});
