import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createParticles, forEachPairWithin, lennardJones, gravity } from '../core/particles.mjs';
import { rng } from '../core/numerics.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);

/** Brute force: all pairs, minimum image. */
function bruteForce(x, n, rc, L) {
  const out = new Set();
  for (let i = 0; i < n; i++) {for (let j = i + 1; j < n; j++) {
    let r2 = 0;
    for (let k = 0; k < 3; k++) { let d = x[3 * j + k] - x[3 * i + k]; if (L[k]) d -= L[k] * Math.round(d / L[k]); r2 += d * d; }
    if (r2 < rc * rc) out.add(`${i},${j}`);
  }}
  return out;
}

test('cell-list neighbour search finds exactly the brute-force pairs (3D, 2D, open, tiny boxes)', () => {
  const rand = rng(8);
  for (const [n, L, rc] of [[300, [10, 10, 10], 2.5], [200, [12, 7, 0], 1.6], [150, [0, 0, 0], 1.2], [40, [5, 5, 5], 2.4], [30, [4.9, 6, 0], 2.0]]) {
    const x = new Float64Array(3 * n);
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) x[3 * i + k] = L[k] ? (rand() * 3 - 1) * L[k] : (k === 2 && L[0] && !L[2] ? 0 : rand() * 8);
    const seen = new Set();
    forEachPairWithin(x, n, rc, L, (i, j, dx, dy, dz, r2) => {
      assert.ok(i < j && !seen.has(`${i},${j}`), 'each pair once');
      seen.add(`${i},${j}`);
      close(dx * dx + dy * dy + dz * dz, r2, 1e-12);
    });
    assert.deepEqual([...seen].sort(), [...bruteForce(x, n, rc, L)].sort(), `n=${n} L=${L}`);
  }
  assert.throws(() => forEachPairWithin(new Float64Array(6), 2, 3, [5, 5, 5], () => {}), /half the box/);
});

test('Lennard-Jones dimer: minimum at 2^(1/6)σ and vibration ω² = 72·2^(−1/3)ε/(σ²μ)', () => {
  const r0 = 2 ** (1 / 6), m = 1, mu = m / 2, A = 1e-6;                // small: the cubic (anharmonic) term shifts the turning point by ~A²
  const sys = createParticles({
    bodies: [{ m, x: [0, 0, 0], v: [0, 0, 0] }, { m, x: [r0 + A, 0, 0], v: [0, 0, 0] }],
    forces: [lennardJones({ shift: 'none', cutoff: 1e6 })], h: 0.001,
  });
  const at = sep => createParticles({ bodies: [{ m, x: [0, 0, 0], v: [0, 0, 0] }, { m, x: [sep, 0, 0], v: [0, 0, 0] }], forces: [lennardJones({ shift: 'none', cutoff: 1e6 })], h: 1 });
  close(at(r0).energy().U, -1, 1e-14, 'well depth −ε');
  const omega = Math.sqrt(72 * 2 ** (-1 / 3) / mu), T = 2 * Math.PI / omega;
  sys.advanceTo(T / 2);                                   // half a period: compressed by A on the other side
  close((sys.x[3] - sys.x[0] - r0) / A, -1, 2e-3, 'displacement after half a period, relative to the amplitude');
});

test('periodic LJ fluid: cell lists = all pairs, momentum to round-off, energy error O(h²), wrap and min image', () => {
  const rand = rng(4), n = 108, rho = 0.6, Lb = Math.cbrt(n / rho), cellsPerSide = 3, a = Lb / cellsPerSide;
  const bodies = [];
  // fcc lattice (4 per cell × 27 cells), small random velocities with zero total momentum
  for (let i = 0; i < 3; i++) {for (let j = 0; j < 3; j++) {for (let k = 0; k < 3; k++) {
    for (const [p, q, r] of [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) {
      bodies.push({ m: 1, x: [(i + p) * a, (j + q) * a, (k + r) * a], v: [rand() - 0.5, rand() - 0.5, rand() - 0.5] });
    }
  }}}
  const mean = [0, 1, 2].map(k => bodies.reduce((s, b) => s + b.v[k], 0) / n);
  bodies.forEach(b => { for (let k = 0; k < 3; k++) b.v[k] -= mean[k]; });
  const run = h => {
    const sys = createParticles({ bodies, forces: [lennardJones()], h, method: 'leapfrog', box: [Lb, Lb, Lb] });
    let worst = 0;
    while (sys.t < 2 - 1e-9) { sys.step(); worst = Math.max(worst, sys.energyBalanceError()); }
    return { sys, worst };
  };
  const coarse = run(0.004), fine = run(0.002);
  close(Math.log2(coarse.worst / fine.worst), 2, 0.25, 'energy error order');
  assert.ok(fine.worst < 1e-4, `energy error ${fine.worst}`);
  assert.ok(fine.sys.momentumDrift() < 1e-13, `momentum ${fine.sys.momentumDrift()}`);
  const w = fine.sys.wrapped();
  assert.ok(w.every(v => v >= 0 && v < Lb));

  // two bodies straddling the boundary interact through the minimum image
  const edge = createParticles({ bodies: [{ m: 1, x: [0.1, 1, 1], v: [0, 0, 0] }, { m: 1, x: [Lb - 1.0, 1, 1], v: [0, 0, 0] }], forces: [lennardJones({ shift: 'none' })], h: 0.001, box: [Lb, Lb, Lb] });
  close(edge.energy().U, 4 * (1.1 ** -12 - 1.1 ** -6), 1e-14, 'distance 1.1 across the wall');
});

test('periodic box refuses long-range forces', () => {
  assert.throws(() => createParticles({ bodies: [{ m: 1, x: [0, 0, 0], v: [0, 0, 0] }], forces: [gravity()], h: 0.1, box: [10, 10, 10] }), /Ewald/);
});
