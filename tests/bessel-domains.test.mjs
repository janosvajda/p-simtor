import { test } from 'node:test';
import assert from 'node:assert/strict';
import { besselJ, besselJZeros, factorial } from '../core/numerics.mjs';
import { laplacianND } from '../core/fields.mjs';
import { conjugateGradient } from '../core/linalg.mjs';
import { waveND } from '../core/pde.mjs';
import { rng } from '../core/numerics.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
// Independent reference: J_n(x) = (1/2π)∫₀^{2π} cos(nτ − x sin τ) dτ; the trapezoidal rule is
// spectrally accurate for this periodic integrand (absolute accuracy ~1e-16).
const jIntegral = (n, x) => {
  const N = Math.ceil(2 * (Math.abs(x) + n) + 80);
  let s = 0;
  for (let i = 0; i < N; i++) { const t = 2 * Math.PI * i / N; s += Math.cos(n * t - x * Math.sin(t)); }
  return s / N;
};

test('besselJ: published values, independent integral representation, tiny values, parity', () => {
  // Abramowitz & Stegun Table 9.1 / standard references
  for (const [n, x, ref] of [[0, 1, 0.7651976865579666], [1, 1, 0.44005058574493355], [0, 10, -0.2459357644513483],
    [1, 10, 0.04347274616886144], [10, 10, 0.2074861066333589], [0, 100, 0.019985850304223122]]) close(besselJ(n, x), ref, 5e-15);
  for (const n of [0, 1, 2, 5, 13, 40]) for (const x of [0.3, 2, 7.7, 25, 60, 300]) close(besselJ(n, x), jIntegral(n, x), 2e-14, `J_${n}(${x})`);
  // where J_n is tiny, the power series Σ(−1)ᵏ(x/2)^{2k+n}/(k!(n+k)!) has no cancellation: relative accuracy must survive
  for (const [n, x] of [[30, 1], [20, 0.01], [8, 0.5]]) {
    let term = (x / 2) ** n / factorial(n), sum = term;
    for (let k = 1; k < 30; k++) { term *= -(x * x / 4) / (k * (n + k)); sum += term; }
    close(besselJ(n, x) / sum, 1, 1e-13, `relative J_${n}(${x})`);
  }
  close(besselJ(3, -2), -besselJ(3, 2), 0);
  assert.equal(besselJ(0, 0), 1); assert.equal(besselJ(4, 0), 0);
  assert.throws(() => besselJ(1.5, 1), /integer/);
});

test('besselJZeros: published zeros, membrane overtone ratios, McMahon asymptotics', () => {
  const z0 = besselJZeros(0, 3);
  [2.404825557695773, 5.520078110286311, 8.653727912911013].forEach((z, i) => close(z0[i], z, 1e-13));
  close(besselJZeros(1, 1)[0], 3.831705970207512, 1e-13);
  close(besselJZeros(2, 1)[0], 5.135622301840683, 1e-13);
  close(besselJZeros(1, 1)[0] / z0[0], 1.5933, 1e-4, 'first overtone of an ideal drum is inharmonic');
  for (const n of [0, 3]) {
    const s = 60, beta = (s + n / 2 - 0.25) * Math.PI, mu = 4 * n * n;
    close(besselJZeros(n, s)[s - 1], beta - (mu - 1) / (8 * beta), 1e-5, `j_${n},${s}`);
  }
  for (const z of besselJZeros(7, 10)) close(besselJ(7, z), 0, 1e-15);
});

/** Disk of radius R on a cell-centred square grid; returns grid, mask and cell radii. */
function disk(R, cells) {
  const h = 2 * R / cells, n = cells + 2, mask = new Uint8Array(n * n), r = new Float64Array(n * n);
  for (let j = 0; j < n; j++) {for (let i = 0; i < n; i++) {
    const k = i + n * j, x = (i - n / 2 + 0.5) * h, y = (j - n / 2 + 0.5) * h;
    r[k] = Math.hypot(x, y); mask[k] = r[k] < R ? 1 : 0;
  }}
  return { grid: { shape: [n, n], spacing: h, bc: 'dirichlet', mask }, h, mask, r };
}

test('mask: all-inside mask changes nothing; masked Neumann has zero net flux', () => {
  const rand = rng(2), shape = [12, 9], u = Float64Array.from({ length: 108 }, () => rand() - 0.5);
  const plain = laplacianND(u, { shape, spacing: 0.3, bc: 'neumann' });
  const masked = laplacianND(u, { shape, spacing: 0.3, bc: 'neumann', mask: new Uint8Array(108).fill(1) });
  plain.forEach((v, i) => assert.equal(masked[i], v));
  const { grid, mask } = disk(1, 30);
  const f = Float64Array.from(mask, m => (m ? rand() : 0));
  const lap = laplacianND(f, { ...grid, maskBc: 'neumann' });
  close(lap.reduce((s, v) => s + v, 0), 0, 1e-9, 'Σ∇²u over the domain = boundary flux = 0');
  assert.ok(lap.every((v, i) => mask[i] || v === 0), 'outside cells stay zero');
});

test('mask: wave energy on a circular drum is conserved to round-off, clamped or free edge', () => {
  const { grid, h, mask, r } = disk(1, 40), rand = rng(3);
  for (const maskBc of ['dirichlet', 'neumann']) {
    const u0 = Float64Array.from(r, (rr, i) => (mask[i] ? Math.exp(-((rr - 0.3) ** 2) / 0.02) : 0));
    const v0 = Float64Array.from(mask, m => (m ? 0.1 * rand() : 0));
    const sim = waveND({ u0, v0, ...grid, maskBc, c: 1, dt: 0.7 * h / Math.SQRT2 });
    const E0 = sim.energy();
    for (let k = 0; k < 400; k++) sim.step();
    close(sim.energy(), E0, 1e-12 * E0, maskBc);
    assert.ok(sim.u.every((v, i) => mask[i] || v === 0));
  }
});

test('mask: clamped circular membrane fundamental → c·j₀₁/R; staircase boundary converges at first order', () => {
  const j01 = besselJZeros(0, 1)[0];
  const fundamental = cells => {
    const { grid, mask, r } = disk(1, cells);
    const A = (x, out) => { laplacianND(x, grid, out); for (let i = 0; i < out.length; i++) out[i] = -out[i]; };
    let x = Float64Array.from(r, (rr, i) => (mask[i] ? besselJ(0, j01 * rr) : 0));
    for (let it = 0; it < 6; it++) {                    // inverse iteration from the sampled exact mode
      const res = conjugateGradient(A, x, { rtol: 1e-12, maxIterations: 5000 });
      assert.ok(res.converged);
      const norm = Math.hypot(...res.x); x = res.x.map(v => v / norm);
    }
    const Ax = new Float64Array(x.length); A(x, Ax);
    return Math.sqrt(x.reduce((s, v, i) => s + v * Ax[i], 0));   // ω/c for unit c
  };
  const errors = [40, 80].map(cells => Math.abs(fundamental(cells) - j01) / j01);
  assert.ok(errors[1] < 0.01, `1% at 80 cells across (${errors[1]})`);
  assert.ok(errors[0] / errors[1] > 1.5, `error falls with resolution (${errors[0]} → ${errors[1]})`);
});

test('contourLines: a circle level set is one closed loop of radius R (O(h²)); two blobs give two loops', async () => {
  const { contourLines, densityLevel } = await import('../core/fields.mjs');
  const n = 101, h = 4 / (n - 1), origin = [-2, -2];
  const f = Float64Array.from({ length: n * n }, (_, k) => { const x = -2 + (k % n) * h, y = -2 + Math.floor(k / n) * h; return x * x + y * y; });
  const lines = contourLines(f, { shape: [n, n], spacing: h, origin }, 1.3 ** 2);
  assert.equal(lines.length, 1); assert.ok(lines[0].closed);
  for (const [x, y] of lines[0].points) close(Math.hypot(x, y), 1.3, 2 * h * h);
  const len = lines[0].points.slice(1).reduce((s, p, i) => s + Math.hypot(p[0] - lines[0].points[i][0], p[1] - lines[0].points[i][1]), 0);
  close(len, 2 * Math.PI * 1.3, 0.01);
  const g = Float64Array.from({ length: n * n }, (_, k) => { const x = -2 + (k % n) * h, y = -2 + Math.floor(k / n) * h; return Math.exp(-((x - 1) ** 2 + y * y) * 8) + Math.exp(-((x + 1) ** 2 + y * y) * 8); });
  assert.equal(contourLines(g, { shape: [n, n], spacing: h, origin }, 0.5).length, 2);

  // 2D Gaussian: the region holding 90 % has density exp(−r²/2σ²) = 0.1 of the peak
  const s = 0.4, G = Float64Array.from({ length: n * n }, (_, k) => { const x = -2 + (k % n) * h, y = -2 + Math.floor(k / n) * h; return Math.exp(-(x * x + y * y) / (2 * s * s)); });
  close(densityLevel(G, 0.9), 0.1, 0.005);
});
