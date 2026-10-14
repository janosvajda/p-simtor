import { test } from 'node:test';
import assert from 'node:assert/strict';
import { matrixVector, solveLinear, leastSquares, eigenSymmetric, conjugateGradient, newtonSolve } from '../core/linalg.mjs';

const close = (a, b, tol = 1e-12) => assert.ok(Math.abs(a - b) <= tol, `${a} vs ${b} (tolerance ${tol})`);
const vectorClose = (a, b, tol = 1e-12) => { assert.equal(a.length, b.length); a.forEach((value, i) => close(value, b[i], tol)); };

test('dense solve pivots, preserves inputs, and is invariant under dimensional rescaling', () => {
  const A = [[0, 2, 1], [1, -2, -3], [2, 3, 1]], b = [1, -5, 4];
  const original = A.map(row => row.slice());
  for (const scale of [1e-180, 1, 1e180]) {
    vectorClose(solveLinear(A.map(row => row.map(value => value * scale)), b.map(value => value * scale)), [2, -1, 3]);
  }
  assert.deepEqual(A, original); assert.deepEqual(b, [1, -5, 4]);
  vectorClose(solveLinear([[1e-200, 0], [0, 1e200]], [3e-200, 4e200]), [3, 4]);
  vectorClose(matrixVector([[1, 2, 3], [-1, 0, 2]], [2, -1, 4]), [12, 6]);
  // A compensated dot product retains a term lost by ordinary left-to-right addition.
  vectorClose(matrixVector([[1e16, 1, -1e16]], [1, 1, 1]), [1], 0);
});

test('dense solve rejects singular, ragged and non-finite problems', () => {
  assert.throws(() => solveLinear([[1, 2], [2, 4]], [1, 2]), /rank-deficient/);
  assert.throws(() => solveLinear([[0]], [1]), /singular/);
  assert.throws(() => solveLinear([[1, 2], [3]], [1, 2]), /finite values/);
  assert.throws(() => solveLinear([[NaN]], [1]), /finite values/);
  assert.throws(() => solveLinear([[1]], [Infinity]), /finite values/);
  assert.throws(() => solveLinear([[1]], [1], { pivotTolerance: -1 }), /tolerance/);
});

test('QR least squares reproduces an independent fit and survives nearly dependent columns', () => {
  const fit = leastSquares([[1, -1], [1, 0], [1, 1]], [1, 0, 3]);
  vectorClose(fit.x, [4 / 3, 1]); close(fit.residualNorm, Math.sqrt(8 / 3));
  const xs = [-2, -1, 0, 1, 2], perturbation = [1, -2, 1, 0, 3];
  const A = xs.map((x, i) => [1 + x, 1 + x + 1e-8 * perturbation[i]]);
  const b = A.map(row => row[0] - 3 * row[1]);
  vectorClose(leastSquares(A, b).x, [1, -3], 2e-7);
  const units = leastSquares([[1e-100, 0], [0, 1e100], [1e-100, 1e100]], [1, 2, 3]);
  close(units.x[0] / 1e100, 1); close(units.x[1] / 1e-100, 2);
  const pivot = leastSquares([[1, 1], [0, 1], [0, 1], [0, 1]], [3, 2, 2, 2]);
  assert.equal(pivot.permutation[0], 1); vectorClose(pivot.x, [1, 2]);
});

test('QR reports rank deficiency instead of manufacturing a unique fit', () => {
  assert.throws(() => leastSquares([[1, 2], [2, 4], [3, 6]], [1, 2, 3]), /rank-deficient/);
  assert.throws(() => leastSquares([[1, 2]], [1]), /rows as columns/);
  vectorClose(leastSquares([[1, 0], [0, 1]], [0, 0]).x, [0, 0], 0);
});

test('dense symmetric eigensystem matches exact Laplacian modes and verifies residuals/orthogonality', () => {
  const n = 9;
  const A = Array.from({ length: n }, (_row, i) => Array.from({ length: n }, (_column, j) => i === j ? 2 : Math.abs(i - j) === 1 ? -1 : 0));
  for (const scale of [1e-180, 1, 1e180]) {
    const scaled = A.map(row => row.map(value => value * scale));
    const { values, vectors } = eigenSymmetric(scaled);
    for (let k = 0; k < n; k++) {
      close(values[k] / scale, 2 - 2 * Math.cos((k + 1) * Math.PI / (n + 1)), 2e-12);
      const Av = matrixVector(scaled, vectors[k]);
      for (let i = 0; i < n; i++) close((Av[i] - values[k] * vectors[k][i]) / scale, 0, 5e-12);
      for (let j = 0; j < n; j++) close(vectors[k].reduce((sum, value, i) => sum + value * vectors[j][i], 0), k === j ? 1 : 0, 2e-14);
    }
  }
});

test('Jacobi handles degeneracy and rejects invalid symmetry or exhausted sweeps', () => {
  vectorClose(eigenSymmetric([[3]]).values, [3], 0);
  vectorClose(eigenSymmetric([[0, 0], [0, 0]]).values, [0, 0], 0);
  vectorClose(eigenSymmetric([[2, 1], [1, 2]]).values, [1, 3], 1e-14);
  assert.throws(() => eigenSymmetric([[1, 2], [0, 1]]), /symmetric/);
  assert.throws(() => eigenSymmetric([[2, 1, 0.4], [1, 3, 1], [0.4, 1, 4]], { maxSweeps: 1, rtol: 1e-15 }), /budget/);
});

test('matrix-free CG solves a manufactured boundary problem, including tiny-dimensional rhs', () => {
  const n = 40;
  const apply = (x, out) => { for (let i = 0; i < n; i++) out[i] = 2 * x[i] - (x[i - 1] ?? 0) - (x[i + 1] ?? 0); };
  const exact = Float64Array.from({ length: n }, (_, i) => Math.sin(Math.PI * (i + 1) / (n + 1)) + 0.3 * Math.sin(4 * Math.PI * (i + 1) / (n + 1)));
  for (const scale of [1, 1e-200]) {
    const b = Float64Array.from({ length: n }, (_, i) => scale * (2 * exact[i] - (exact[i - 1] ?? 0) - (exact[i + 1] ?? 0)));
    const result = conjugateGradient(apply, b, { rtol: 1e-10 });
    assert.equal(result.converged, true);
    vectorClose(Float64Array.from(result.x, value => value / scale), exact, 1e-10);
    assert.ok(result.residualNorm <= 1e-10 * Math.hypot(...b));
  }
});

test('CG preconditioning, iteration failure and operator breakdown are explicit', () => {
  const A = [[1, 0, 0], [0, 100, 0], [0, 0, 10000]], b = [1, 2, 3];
  const M = (r, out) => { out[0] = r[0]; out[1] = r[1] / 100; out[2] = r[2] / 10000; };
  const solved = conjugateGradient(A, b, { preconditioner: M });
  assert.equal(solved.converged, true); assert.equal(solved.iterations, 1);
  vectorClose(solved.x, [1, 0.02, 0.0003]);
  const failed = conjugateGradient(A, b, { maxIterations: 1 });
  assert.equal(failed.converged, false); assert.ok(failed.residualNorm > 1e-10 * Math.hypot(...b));
  assert.throws(() => conjugateGradient([[-1]], [1]), /SPD operator/);
  assert.throws(() => conjugateGradient([[1]], [1], { preconditioner: [[-1]] }), /SPD preconditioner/);
  assert.throws(() => conjugateGradient(() => [NaN], [1]), /finite values/);
  assert.throws(() => conjugateGradient((_x, _out) => undefined, [1]), /finite values/);
  assert.equal(conjugateGradient([[1]], [0]).converged, true);
});

test('damped Newton solves coupled equations with analytic and numerical Jacobians', () => {
  const f = ([x, y]) => [x * x + y * y - 1, x - y];
  for (const jacobian of [undefined, ([x, y]) => [[2 * x, 2 * y], [1, -1]]]) {
    const result = newtonSolve(f, [0.8, 0.3], { jacobian });
    assert.equal(result.converged, true); vectorClose(result.x, [Math.SQRT1_2, Math.SQRT1_2], 1e-10);
  }
  const damped = newtonSolve(([x]) => [x ** 3 - 1], [0.1]);
  assert.equal(damped.converged, true); close(damped.x[0], 1, 1e-10);
  const domain = newtonSolve(([x]) => [Math.log(x)], [10]);
  assert.equal(domain.converged, true); close(domain.x[0], 1, 1e-10);
});

test('Newton distinguishes initial roots, insufficient work and invalid functions', () => {
  assert.equal(newtonSolve(([x]) => [x], [0]).iterations, 0);
  const unfinished = newtonSolve(([x]) => [x * x - 2], [1], { maxIterations: 1 });
  assert.equal(unfinished.converged, false); assert.equal(unfinished.reason, 'iteration budget exhausted');
  const stalled = newtonSolve(([x]) => [x ** 3 - 1], [0.1], { maxBacktracks: 1 });
  assert.equal(stalled.converged, false); vectorClose(stalled.x, [0.1], 0);
  assert.throws(() => newtonSolve(() => [NaN], [1]), /finite values/);
  assert.throws(() => newtonSolve(() => [1], [1]), /singular|rank-deficient/);
  assert.throws(() => newtonSolve(([x]) => [x], [1], { xScale: [0] }), /positive/);
  const reused = new Float64Array(1);
  assert.equal(newtonSolve(([x]) => { reused[0] = x * x - 2; return reused; }, [1]).converged, true);
});
