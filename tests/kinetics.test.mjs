import { test } from 'node:test';
import assert from 'node:assert/strict';
import { expm } from '../core/linalg.mjs';
import { linearChain } from '../core/ode.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);

test('expm: rotation, nilpotent, diagonal, e^A e^−A = I, large norm', () => {
  const t = 2.3, R = expm([[0, t], [-t, 0]]);
  close(R[0][0], Math.cos(t), 1e-15); close(R[0][1], Math.sin(t), 1e-15); close(R[1][0], -Math.sin(t), 1e-15);
  const N = expm([[0, 1, 0], [0, 0, 1], [0, 0, 0]]);
  [[1, 1, 0.5], [0, 1, 1], [0, 0, 1]].forEach((row, i) => row.forEach((v, j) => close(N[i][j], v, 1e-15)));
  const D = expm([[-40, 0], [0, 3]]);
  close(D[0][0] / Math.exp(-40), 1, 1e-13); close(D[1][1], Math.exp(3), 1e-13 * Math.exp(3));
  const M = [[0.3, -1.2, 0.5], [2, 0.1, -0.7], [-0.4, 0.9, -0.2]].map(r => r.map(v => 6 * v));
  const P = expm(M), Q = expm(M.map(r => r.map(v => -v)));
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) close(P[i].reduce((s, v, k) => s + v * Q[k][j], 0), i === j ? 1 : 0, 1e-10);
});

test('linearChain: Bateman two-member formula, conservation, branching', () => {
  const l1 = 0.7, l2 = 0.2, ch = linearChain({ rates: [l1, l2, 0], links: [[0, 1, 1], [1, 2, 1]] });
  for (const t of [0, 0.5, 3, 20]) {
    const x = ch.at(t, [1, 0, 0]);
    close(x[0], Math.exp(-l1 * t), 1e-15);
    close(x[1], l1 / (l2 - l1) * (Math.exp(-l1 * t) - Math.exp(-l2 * t)), 1e-15);
    close(x[0] + x[1] + x[2], 1, 1e-14, 'nothing is lost');
  }
  // Bi-212 → Po-212 (64.06 %) or Tl-208 (35.94 %) → Pb-208: at long times all is lead, split by branch along the way
  const lam = h => Math.LN2 / h;
  const bi = linearChain({ rates: [lam(60.55 * 60), lam(294.4e-9), lam(3.053 * 60), 0], links: [[0, 1, 0.6406], [0, 2, 0.3594], [1, 3, 1], [2, 3, 1]] });
  const x = bi.at(1e6, [1, 0, 0, 0]);
  close(x[3], 1, 1e-12);
  assert.throws(() => linearChain({ rates: [1, 1], links: [[0, 1, 1]] }), /expm/);
});

test('linearChain agrees with expm on the stiff radon-222 progeny chain (µs to decades)', () => {
  const lam = h => Math.LN2 / h, min = 60, yr = 365.25 * 86400;
  const rates = [lam(3.098 * min), lam(27.06 * min), lam(19.9 * min), lam(164.3e-6), lam(22.2 * yr)];
  const ch = linearChain({ rates, links: [[0, 1, 1], [1, 2, 1], [2, 3, 1], [3, 4, 1]] });
  for (const t of [60, 1800, 7200]) {
    const exact = ch.at(t, [1, 0, 0, 0, 0]);
    const E = expm(ch.matrix.map(r => Array.from(r, v => v * t)));
    const viaExpm = E.map(row => row[0]);
    exact.forEach((v, i) => close(viaExpm[i], v, 1e-12 + 1e-8 * Math.abs(v), `species ${i} at ${t} s`));
  }
});
