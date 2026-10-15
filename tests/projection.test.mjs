import { test } from 'node:test';
import assert from 'node:assert/strict';
import { rotation4, apply4, project4, stereographic, inverseStereographic, axisAngle, hopfMap, hopfFibre } from '../core/projection.mjs';
import { qexp, qmul, qdag } from '../core/lattice-gauge.mjs';
import { rng } from '../core/numerics.mjs';

const close = (a, b, tol = 1e-12) => a.every((v, i) => Math.abs(v - b[i]) < tol);

test('rotation4 is orthogonal with determinant +1; a single plane rotates only that plane', () => {
  const R = rotation4({ xy: 0.3, xw: 1.1, yz: -0.4, zw: 2 });
  for (let i = 0; i < 4; i++) {for (let j = 0; j < 4; j++) {
    let d = 0; for (let k = 0; k < 4; k++) d += R[4 * i + k] * R[4 * j + k];
    assert.ok(Math.abs(d - (i === j ? 1 : 0)) < 1e-14);
  }}
  const Rw = rotation4({ xw: Math.PI / 2 });
  assert.ok(close(apply4(Rw, [1, 0, 0, 0]), [0, 0, 0, 1]));
  assert.ok(close(apply4(Rw, [0, 1, 0, 0]), [0, 1, 0, 0]));
});

test('project4: orthographic drops w; perspective scales by eye/(eye − w) and refuses points behind the eye', () => {
  assert.deepEqual(project4([1, 2, 3, 9], { mode: 'orthographic' }), [1, 2, 3]);
  assert.ok(close(project4([1, 2, 3, 1], { eye: 3 }), [1.5, 3, 4.5]));
  assert.throws(() => project4([0, 0, 0, 3], { eye: 3 }), /behind/);
});

test('stereographic and its inverse are inverse; the identity maps to the origin', () => {
  const r = rng(3);
  for (let k = 0; k < 50; k++) {
    const q = qexp([r() * 2 - 1, r() * 2 - 1, r() * 2 - 1]);
    assert.ok(close(inverseStereographic(stereographic(q)), Array.from(q), 1e-12));
  }
  assert.ok(close(stereographic([1, 0, 0, 0]), [0, 0, 0]));
});

test('axisAngle: radius = angle from the identity; conjugation q → g q g† rotates it rigidly', () => {
  const v = [0.4, -0.2, 0.9], q = qexp(v);
  assert.ok(close(axisAngle(q), v));
  const g = qexp([0.3, 0.7, -0.1]), c = qmul(qmul(g, q), qdag(g));
  assert.ok(Math.abs(Math.hypot(...axisAngle(c)) - Math.hypot(...v)) < 1e-14);
});

test('hopfMap sends every point of hopfFibre(p) to p, including the south pole; fibres are unit quaternions', () => {
  for (const p of [[0, 0, 1], [0.6, 0, 0.8], [0.36, -0.48, 0.8], [1, 0, 0], [0, 0, -1]]) {
    for (const q of hopfFibre(p, 12)) {
      assert.ok(Math.abs(Math.hypot(...q) - 1) < 1e-14);
      assert.ok(close(hopfMap(q), p, 1e-12), `${p}: ${hopfMap(q)}`);
    }
  }
});

test('hopfMap is unchanged by left multiplication with exp(i t σ₃): fibres are U(1) orbits', () => {
  const q = qexp([0.3, -0.8, 0.5]), h = qmul(qexp([0, 0, 0.9]), q);
  assert.ok(close(hopfMap(q), hopfMap(h), 1e-14));
});
