import { test } from 'node:test';
import assert from 'node:assert/strict';
import { jackknife, bootstrap, createHistogram, chiSquareTest, meanSem, pull } from '../core/data.mjs';
import { rng, erf, gammaP, gammaQ } from '../core/numerics.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);

// AR(1) series x_t = ρ x_{t−1} + √(1−ρ²) ξ_t: unit variance, τ_int = (1+ρ)/(2(1−ρ)).
const ar1 = (n, rho, rand) => {
  const xs = new Float64Array(n), c = Math.sqrt(1 - rho * rho);
  let x = rand.normal();
  for (let i = 0; i < n; i++) { x = rho * x + c * rand.normal(); xs[i] = x; }
  return xs;
};
const sd = xs => meanSem(xs).sd;

test('gammaQ / gammaP: exact special cases', () => {
  for (const x of [0, 0.1, 0.7, 1, 2.5, 6, 20]) {
    close(gammaQ(1, x), Math.exp(-x), 1e-14, `Q(1, ${x})`);
    close(gammaQ(0.5, x), 1 - erf(Math.sqrt(x)), 1e-13, `Q(½, ${x})`);
    close(gammaP(1, x) + gammaQ(1, x), 1, 1e-15, 'P + Q');
    // χ² with 2 dof: p = e^{−χ²/2}; with 4 dof: (1 + χ²/2) e^{−χ²/2}
    close(gammaQ(1, x / 2), Math.exp(-x / 2), 1e-14, 'χ² dof 2');
    close(gammaQ(2, x / 2), (1 + x / 2) * Math.exp(-x / 2), 1e-13, 'χ² dof 4');
  }
  close(gammaP(3, 3), 1 - 8.5 * Math.exp(-3), 1e-14, 'P(3,3) = 1 − (1+3+9/2)e^{−3}');
  assert.equal(gammaQ(2, Infinity), 0);
  assert.throws(() => gammaQ(0, 1), /gammaQ/);
});

test('jackknife of the mean equals the standard error of the block means', () => {
  const xs = ar1(1000, 0.5, rng(2));
  for (const bins of [10, 50]) {
    const size = 1000 / bins, means = Array.from({ length: bins }, (_, b) => xs.slice(b * size, (b + 1) * size).reduce((a, v) => a + v, 0) / size);
    const jk = jackknife(xs, m => m, { bins }), ref = meanSem(means);
    close(jk.value, ref.mean, 1e-12, 'value');
    close(jk.error, ref.sem, 1e-12, 'error');
    close(jk.bias, 0, 1e-12, 'bias of a linear estimator');
    assert.equal(jk.blocks, bins);
  }
});

test('jackknife: blocks drop the remainder at the start; arrays and typed arrays both work', () => {
  const xs = [100, 1, 2, 3, 4];                 // n = 5, 2 blocks of 2: the leading 100 is dropped
  assert.equal(jackknife(xs, m => m, { bins: 2 }).value, 2.5);
  const v = jackknife([[1, 2], [3, 4], Float64Array.of(5, 6), [7, 8]], m => [m[0], m[1] / m[0]], { bins: 4 });
  assert.deepEqual(v.value, [4, 5 / 4]);
  assert.throws(() => jackknife([1], m => m), /jackknife/);
  assert.throws(() => jackknife([1, 2, 3], m => m, { bins: 4 }), /jackknife/);
});

test('jackknife of nonlinear estimators matches the spread over independent repeats (correlated data)', () => {
  const rand = rng(7), repeats = 200, n = 8000, rho = 0.8;
  const logs = [], ratios = [], jkLog = [], jkRatio = [];
  for (let r = 0; r < repeats; r++) {
    const a = ar1(n, rho, rand), b = ar1(n, rho, rand);
    const samples = Array.from({ length: n }, (_, i) => [3 + a[i], 2 + 0.5 * a[i] + b[i]]);
    const L = jackknife(samples, m => Math.log(m[0]), { bins: 40 });
    const R = jackknife(samples, m => [m[0] / m[1], m[0] * m[1]], { bins: 40 });   // vector estimator
    logs.push(L.value); jkLog.push(L.error); ratios.push(R.value[0]); jkRatio.push(R.error[0]);
    assert.equal(R.value.length, 2);
  }
  const mean = xs => xs.reduce((s, v) => s + v, 0) / xs.length;
  close(mean(jkLog) / sd(logs), 1, 0.2, 'log ⟨x⟩: jackknife / true spread');
  close(mean(jkRatio) / sd(ratios), 1, 0.2, '⟨a⟩/⟨b⟩: jackknife / true spread');
});

test('jackknife: NaN from any block propagates to the error', () => {
  const xs = Array.from({ length: 100 }, (_, i) => (i < 10 ? 0.5 : 1));   // block 0 holds the 0.5s
  const jk = jackknife(xs, m => (m >= 1 ? NaN : Math.log(m)), { bins: 10 });
  assert.ok(Number.isFinite(jk.value));
  assert.ok(Number.isNaN(jk.error));
  const v = jackknife(xs, m => [m, m >= 1 ? NaN : m], { bins: 10 });
  assert.ok(Number.isFinite(v.error[0]) && Number.isNaN(v.error[1]));
});

test('bootstrap error agrees with the jackknife; quantiles bracket the value; reproducible', () => {
  const a = ar1(4000, 0.5, rng(3)), samples = Array.from(a, x => [2 + x, 1 + 0.3 * x]);
  const est = m => m[0] / m[1];
  const jk = jackknife(samples, est, { bins: 40 }), bs = bootstrap(samples, est, { bins: 40, resamples: 2000, rand: rng(4) });
  close(bs.value, jk.value, 1e-12, 'value');
  close(bs.error / jk.error, 1, 0.12, 'bootstrap / jackknife error');
  const [lo, hi] = bs.quantiles([0.16, 0.84]);
  assert.ok(lo < bs.value && bs.value < hi);
  close((hi - lo) / 2 / bs.error, 1, 0.15, '68 % half-width ≈ σ');
  assert.equal(bootstrap(samples, est).error, bootstrap(samples, est).error, 'default seed is fixed');
  const v = bootstrap(samples, m => [m[0], NaN], { resamples: 50 });
  assert.ok(Number.isFinite(v.error[0]) && Number.isNaN(v.error[1]) && Number.isNaN(v.quantiles(0.5)[1]));
});

test('createHistogram: counts, under/overflow, weights, density integrates to the in-range fraction', () => {
  const h = createHistogram({ bins: 4, range: [0, 2] });
  for (const x of [-1, 0, 0.2, 0.5, 1.99, 2, 2.5, NaN]) h.add(x);
  h.add(1.2, 3);
  assert.deepEqual(Array.from(h.counts), [2, 1, 3, 2]);
  assert.equal(h.underflow, 1); assert.equal(h.overflow, 1); assert.equal(h.nan, 1); assert.equal(h.total, 11);
  assert.deepEqual(h.errors(), [Math.SQRT2, 1, 3, Math.SQRT2]);
  assert.deepEqual(h.edges, [0, 0.5, 1, 1.5, 2]);
  assert.deepEqual(h.centers, [0.25, 0.75, 1.25, 1.75]);
  const rand = rng(5), g = createHistogram({ bins: 40, range: [-1, 1] });
  for (let i = 0; i < 20000; i++) g.add(rand.normal());
  const d = g.density(), integral = d.values.reduce((s, v) => s + v * g.width, 0);
  close(integral, (g.total - g.underflow - g.overflow) / g.total, 1e-12, '∫ density');
  close(integral, erf(1 / Math.SQRT2), 0.015, 'in-range fraction of N(0,1)');
  close(d.errors[20], Math.sqrt(g.counts[20]) / (g.total * g.width), 1e-15, 'density error');
  assert.throws(() => createHistogram({ bins: 4 }), /createHistogram/);
});

test('chiSquareTest: exact small case; p-values on Poisson data are uniform', () => {
  const r = chiSquareTest([12, 8, 0], [10, 10, 0]);
  close(r.chi2, 0.8, 1e-12); assert.equal(r.dof, 2); assert.equal(r.skipped, 1);
  close(r.p, Math.exp(-0.4), 1e-14, 'p for dof 2');
  close(chiSquareTest([1, 3], [2, 2], { errors: [0.5, 1] }).chi2, 5, 1e-12, 'with errors');
  const rand = rng(6), poisson = mu => { let k = 0, p = Math.exp(-mu), s = p; const u = rand(); while (u > s) { k++; p *= mu / k; s += p; } return k; };
  const expected = [20, 35, 50, 40, 25, 30], ps = [];
  for (let rep = 0; rep < 200; rep++) ps.push(chiSquareTest(expected.map(poisson), expected).p);
  const m = ps.reduce((s, v) => s + v, 0) / ps.length;
  assert.ok(m > 0.4 && m < 0.6, `mean p ${m}`);
  assert.ok(ps.filter(p => p < 0.1).length < 40, 'not too many small p');
});

test('meanSem and pull', () => {
  const r = meanSem([1, 2, 3, 4]);
  assert.equal(r.mean, 2.5); assert.equal(r.n, 4);
  close(r.sd, Math.sqrt(5 / 3), 1e-15); close(r.sem, Math.sqrt(5 / 3) / 2, 1e-15);
  assert.equal(pull(11, 0.5, 10), 2);
});
