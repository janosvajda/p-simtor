import { test } from 'node:test';
import assert from 'node:assert/strict';
import { cholesky, eigenGeneralized } from '../core/linalg.mjs';

const mul = (A, B) => A.map((r, i) => B[0].map((_, j) => r.reduce((s, _v, k) => s + A[i][k] * B[k][j], 0)));
const T = A => A[0].map((_, j) => A.map(r => r[j]));

test('cholesky reproduces A = LLᵀ and rejects indefinite matrices', () => {
  const A = [[4, 12, -16], [12, 37, -43], [-16, -43, 98]];
  const L = cholesky(A);
  assert.deepEqual(L.map(r => Array.from(r)), [[2, 0, 0], [6, 1, 0], [-8, 5, 3]]);
  const Ls = L.map(r => Array.from(r)), back = mul(Ls, T(Ls));
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) assert.ok(Math.abs(back[i][j] - A[i][j]) < 1e-12);
  assert.throws(() => cholesky([[1, 2], [2, 1]]), /positive definite/);
  assert.throws(() => cholesky([[1, 2], [0, 1]]), /symmetric/);
});

test('eigenGeneralized: A v = λ B v with B-normalised vectors', () => {
  const A = [[2, 1, 0], [1, 3, 1], [0, 1, 4]], B = [[2, 0.5, 0], [0.5, 1, 0.2], [0, 0.2, 1.5]];
  const { values, vectors } = eigenGeneralized(A, B);
  for (let i = 0; i < 3; i++) {
    const v = vectors[i];
    const Av = A.map(r => r.reduce((s, x, k) => s + x * v[k], 0)), Bv = B.map(r => r.reduce((s, x, k) => s + x * v[k], 0));
    for (let k = 0; k < 3; k++) assert.ok(Math.abs(Av[k] - values[i] * Bv[k]) < 1e-12);
    assert.ok(Math.abs(v.reduce((s, x, k) => s + x * Bv[k], 0) - 1) < 1e-12);
  }
  assert.ok(values[0] < values[1] && values[1] < values[2]);
});

test('GEVP isolates the ground state of a two-state correlation matrix', () => {
  // C_ij(t) = Σ_n Z_in Z_jn e^{−E_n t} with E = 0.8, 1.9: a single operator is contaminated
  const E = [0.8, 1.9], Z = [[1, 0.6], [0.7, -0.9]];
  const C = t => Z.map(zi => Z.map(zj => zi[0] * zj[0] * Math.exp(-E[0] * t) + zi[1] * zj[1] * Math.exp(-E[1] * t)));
  const { values } = eigenGeneralized(C(1), C(0));
  assert.ok(Math.abs(-Math.log(values[1]) - 0.8) < 1e-12);
  assert.ok(Math.abs(-Math.log(values[0]) - 1.9) < 1e-12);
  const single = Math.log(C(0)[0][0] / C(1)[0][0]);
  assert.ok(single > 0.85, 'a single operator overestimates E₀');
});
