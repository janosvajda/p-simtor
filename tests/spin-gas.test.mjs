import { test } from 'node:test';
import assert from 'node:assert/strict';
import { ketAlong, rotate, bloch, probability, measure } from '../core/spin.mjs';
import { maxwellSpeed, effusiveSpeed } from '../core/probability.mjs';
import { rng, rngFor, gaussLegendre } from '../core/numerics.mjs';
import { ksTest } from '../core/data.mjs';
import { u, ATOMIC_WEIGHT, muB } from '../core/constants.mjs';

const close = (a, b, t = 1e-12) => a.every((v, i) => Math.abs(v - b[i]) < t);

test('ketAlong is the ±1 eigenstate: Bloch vector ±n̂, unit norm', () => {
  for (const n of [[0, 0, 1], [1, 0, 0], [0.3, -0.5, 0.81], [0, 0, -1]]) {
    const l = Math.hypot(...n), nh = n.map(v => v / l);
    for (const s of [1, -1]) {
      const k = ketAlong(n, s);
      assert.ok(Math.abs(k[0] ** 2 + k[1] ** 2 + k[2] ** 2 + k[3] ** 2 - 1) < 1e-14);
      assert.ok(close(bloch(k), nh.map(v => s * v)), `${n} ${s}`);
    }
  }
});

test('rotate precesses the Bloch vector rigidly about n̂ (Larmor) and preserves the norm', () => {
  const psi = ketAlong([1, 0, 0]), r = rotate(psi, [0, 0, 1], Math.PI / 2);
  assert.ok(close(bloch(r), [0, 1, 0]));
  const full = rotate(psi, [0.2, 0.7, -0.4], 2 * Math.PI);           // 2π: −ψ (spinor sign), same Bloch vector
  assert.ok(close(Array.from(full), Array.from(psi).map(v => -v)));
});

test('Born rule: P(+ along n̂) = cos²(θ/2) for a state prepared along ẑ; sequential z→x→z gives ½', () => {
  for (const th of [0, 0.5, Math.PI / 3, 2, Math.PI]) {
    assert.ok(Math.abs(probability(ketAlong([0, 0, 1]), [Math.sin(th), 0, Math.cos(th)], 1) - Math.cos(th / 2) ** 2) < 1e-14);
  }
  const r = rng(5);
  let up = 0; const N = 20000;
  for (let i = 0; i < N; i++) {
    const a = measure(ketAlong([0, 0, 1]), [1, 0, 0], r);
    if (measure(a.state, [0, 0, 1], r).s > 0) up++;
  }
  assert.ok(Math.abs(up / N - 0.5) < 4 * Math.sqrt(0.25 / N));
});

test('Maxwell and effusive speed laws: pdfs normalised, cdfs consistent, exact samplers pass KS, known means', () => {
  const m = ATOMIC_WEIGHT.Ag * u, T = 1273;
  for (const D of [maxwellSpeed(T, m), effusiveSpeed(T, m)]) {
    const total = gaussLegendre(D.pdf, 0, 12 * D.sigma, 64);
    assert.ok(Math.abs(total - 1) < 1e-10);
    assert.ok(Math.abs(gaussLegendre(D.pdf, 0, D.sigma, 32) - D.cdf(D.sigma)) < 1e-12);
    assert.ok(Math.abs(gaussLegendre(v => v * D.pdf(v), 0, 12 * D.sigma, 64) - D.mean) < 1e-9 * D.mean);
    const r = rng(9), xs = Array.from({ length: 5000 }, () => D.sample(r));
    assert.ok(ksTest(xs, D.cdf).p > 0.01);
  }
  assert.ok(Math.abs(muB / 9.2740100783e-24 - 1) < 1e-10);
});

test('rngFor gives independent, reproducible streams per item', () => {
  const a = rngFor(7, 3), b = rngFor(7, 3), c = rngFor(7, 4);
  const sa = [a(), a(), a()], sb = [b(), b(), b()], sc = [c(), c(), c()];
  assert.deepEqual(sa, sb);
  assert.notDeepEqual(sa, sc);
});
