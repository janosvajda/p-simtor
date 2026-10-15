import { test } from 'node:test';
import assert from 'node:assert/strict';
import { curveFit, jackknifeRefit } from '../core/data.mjs';
import { rng } from '../core/numerics.mjs';

const damped = (t, [A, g, w, ph, c]) => A * Math.exp(-g * t) * Math.cos(w * t + ph) + c;

test('curveFit recovers a damped sine exactly from noise-free data', () => {
  const truth = [2, 0.05, 1.3, 0.4, 0.7], xs = Array.from({ length: 200 }, (_, i) => i * 0.2), ys = xs.map(x => damped(x, truth));
  const f = curveFit(damped, xs, ys, [1.8, 0.04, 1.2, 0.3, 0.6]);
  assert.ok(f.converged);
  f.p.forEach((v, k) => assert.ok(Math.abs(v - truth[k]) < 1e-8, `p${k} ${v} vs ${truth[k]}`));
});

test('curveFit errors are honest: pulls over 200 noisy datasets have rms ≈ 1; χ² p-value is uniform-ish', () => {
  const truth = [2, 0.05, 1.3, 0.4, 0.7], xs = Array.from({ length: 120 }, (_, i) => i * 0.25), r = rng(11), s = 0.05;
  const pulls = [], ps = [];
  for (let k = 0; k < 200; k++) {
    const ys = xs.map(x => damped(x, truth) + s * r.normal());
    const f = curveFit(damped, xs, ys, [1.9, 0.045, 1.28, 0.35, 0.68], { sigma: s });
    pulls.push((f.p[2] - truth[2]) / f.errors[2]); ps.push(f.p_value);
  }
  const rms = Math.sqrt(pulls.reduce((a, v) => a + v * v, 0) / pulls.length);
  assert.ok(rms > 0.85 && rms < 1.15, `pull rms ${rms}`);
  const meanP = ps.reduce((a, v) => a + v, 0) / ps.length;
  assert.ok(meanP > 0.4 && meanP < 0.6, `mean p ${meanP}`);
});

test('curveFit without sigma estimates errors from the scatter; a straight line matches the closed form', () => {
  const xs = [0, 1, 2, 3, 4, 5], ys = [1.1, 2.9, 5.2, 7.1, 8.8, 11.2];
  const f = curveFit((x, [a, b]) => a + b * x, xs, ys, [0, 1]);
  const n = xs.length, mx = xs.reduce((a, v) => a + v) / n, my = ys.reduce((a, v) => a + v) / n;
  const sxx = xs.reduce((a, v) => a + (v - mx) ** 2, 0), b = xs.reduce((a, v, i) => a + (v - mx) * (ys[i] - my), 0) / sxx;
  const res = ys.map((y, i) => y - (my - b * mx) - b * xs[i]), s2 = res.reduce((a, v) => a + v * v, 0) / (n - 2);
  assert.ok(Math.abs(f.p[1] - b) < 1e-9 && Math.abs(f.errors[1] - Math.sqrt(s2 / sxx)) < 1e-9);
});

test('jackknifeRefit: refitting the mean gives the standard error; vector results; too few samples throw', () => {
  const r = rng(3), xs = Array.from({ length: 400 }, () => r.normal() + 5);
  const jk = jackknifeRefit(xs, sub => sub.reduce((a, v) => a + v, 0) / sub.length, { bins: 40 });
  const sd = Math.sqrt(xs.reduce((a, v) => a + (v - 5) ** 2, 0) / xs.length);
  assert.ok(Math.abs(jk.error / (sd / Math.sqrt(400)) - 1) < 0.3, `${jk.error}`);
  const vec = jackknifeRefit(xs, sub => [Math.min(...sub), sub.length]);
  assert.equal(vec.value.length, 2);
  assert.throws(() => jackknifeRefit([1], x => x[0]), /2 samples/);
});
