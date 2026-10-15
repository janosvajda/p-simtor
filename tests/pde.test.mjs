import { test } from 'node:test';
import assert from 'node:assert/strict';
import { diffusion1D, waveND, stationaryStates1D } from '../core/pde.mjs';
import { schrodingerND } from '../core/fields.mjs';
import { eigenTridiagonal, rng } from '../core/numerics.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const grid = (n, x0, dx) => Float64Array.from({ length: n }, (_, i) => x0 + i * dx);

test('eigenTridiagonal: exact eigenvalues of the discrete Laplacian, orthonormal vectors', () => {
  const n = 40, { values, vectors } = eigenTridiagonal(new Float64Array(n).fill(2), new Float64Array(n - 1).fill(-1));
  const exact = [...Array(n)].map((_, k) => 2 - 2 * Math.cos((k + 1) * Math.PI / (n + 1))).sort((a, b) => a - b);
  values.slice().sort((a, b) => a - b).forEach((v, i) => close(v, exact[i], 1e-13));
  for (const [i, j] of [[0, 0], [0, 5], [7, 7], [3, 30]]) {
    const dot = vectors[i].reduce((s, a, k) => s + a * vectors[j][k], 0);
    close(dot, i === j ? 1 : 0, 1e-13);
  }
});

test('diffusion: Gaussian spreads as the exact heat kernel; total conserved with zero-flux walls', () => {
  const D = 0.5, n = 401, dx = 0.05, x = grid(n, -10, dx), s0 = 0.5;
  const u0 = Float64Array.from(x, xi => Math.exp(-xi * xi / (2 * s0 * s0)));
  const sim = diffusion1D({ u0, D, dx, dt: 0.005, bc: 'neumann' });
  const total0 = sim.total();
  sim.advanceTo(1);
  const s2 = s0 * s0 + 2 * D * sim.t;
  x.forEach((xi, i) => close(sim.u[i], s0 / Math.sqrt(s2) * Math.exp(-xi * xi / (2 * s2)), 2e-4));
  close(sim.total(), total0, 1e-12 * total0, 'conservation');
});

test('diffusion: 2nd-order convergence in dt; periodic boundaries conserve the total', () => {
  const n = 64, dx = 1 / n, x = grid(n, 0, dx);
  const u0 = Float64Array.from(x, xi => Math.sin(2 * Math.PI * xi) + 2);
  const at = dt => { const s = diffusion1D({ u0, D: 0.1, dx, dt, bc: 'periodic' }); s.advanceTo(0.2); return s; };
  const a = at(0.02), b = at(0.01), c = at(0.005);
  const e1 = Math.max(...a.u.map((v, i) => Math.abs(v - b.u[i]))), e2 = Math.max(...b.u.map((v, i) => Math.abs(v - c.u[i])));
  close(Math.log2(e1 / e2), 2, 0.1, 'order in dt');
  close(c.total(), 2, 1e-12);
});

test('wave: standing mode oscillates at ω = cπ/L; CFL is enforced', () => {
  const L = 1, n = 201, dx = L / (n + 1), c = 2, x = grid(n, dx, dx);
  const u0 = Float64Array.from(x, xi => Math.sin(Math.PI * xi / L));
  const sim = waveND({ u0, shape: [n], spacing: dx, c, dt: 0.5 * dx / c, bc: 'dirichlet' });
  sim.advanceTo(L / c);                                   // half the fundamental period
  x.forEach((xi, i) => close(sim.u[i], Math.sin(Math.PI * xi / L) * Math.cos(Math.PI * c * sim.t / L), 2e-4));
  assert.throws(() => waveND({ u0, shape: [n], spacing: dx, c, dt: 2 * dx / c }), /Courant/);
  assert.throws(() => waveND({ u0: new Float64Array(n * n), shape: [n, n], spacing: dx, c, dt: 0.8 * dx / c }), /Courant/, '2D limit is c·dt·√2/dx ≤ 1');
});

test('wave: the (1,1) mode of a square membrane has ω = cπ√2/L', () => {
  const L = 1, n = 61, d = L / (n + 1), c = 1;
  const u0 = new Float64Array(n * n);
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) u0[j * n + i] = Math.sin(Math.PI * (i + 1) * d) * Math.sin(Math.PI * (j + 1) * d);
  const sim = waveND({ u0, shape: [n, n], spacing: d, c, dt: 0.5 * d / c });
  sim.advanceTo(1);
  const mid = Math.floor(n / 2) * n + Math.floor(n / 2);
  close(sim.u[mid], u0[mid] * Math.cos(c * Math.PI * Math.SQRT2 / L * sim.t), 2e-3);
});

test('wave: discrete energy is conserved to round-off in 1D/2D/3D for every boundary type', () => {
  const rand = rng(7);
  for (const shape of [[64], [24, 16], [10, 8, 6]]) {
    for (const bc of ['dirichlet', 'neumann', 'periodic']) {
      const n = shape.reduce((a, b) => a * b, 1), spacing = shape.map((_, k) => 0.1 + 0.02 * k);
      const u0 = Float64Array.from({ length: n }, () => rand() - 0.5), v0 = Float64Array.from({ length: n }, () => rand() - 0.5);
      const dt = 0.9 / Math.sqrt(spacing.reduce((s, h) => s + 1 / (h * h), 0));
      const sim = waveND({ u0, v0, shape, spacing, bc, c: 1, dt });
      const E0 = sim.energy();
      assert.ok(E0 > 0, 'positive below the CFL limit');
      for (let k = 0; k < 300; k++) sim.step();
      close(sim.energy(), E0, 1e-12 * E0, `${shape.length}D ${bc}`);
    }
  }
});

test('Schrödinger (fields.mjs): free Gaussian spreads exactly; coherent state oscillates at ω', () => {
  const n = 1024, dx = 0.1, x = grid(n, -51.2, dx), s0 = 1, k0 = 2;
  const moments = sim => {
    const p = sim.probability(), N = sim.norm();
    let m = 0, m2 = 0;
    for (let i = 0; i < n; i++) { m += x[i] * p[i] * dx; m2 += x[i] * x[i] * p[i] * dx; }
    return { mean: m / N, sd: Math.sqrt(m2 / N - (m / N) ** 2) };
  };
  const re = Float64Array.from(x, xi => Math.exp(-xi * xi / (4 * s0 * s0)) * Math.cos(k0 * xi));
  const im = Float64Array.from(x, xi => Math.exp(-xi * xi / (4 * s0 * s0)) * Math.sin(k0 * xi));
  const free = schrodingerND({ re0: re, im0: im, shape: [n], spacing: dx, dt: 0.01 });
  free.advanceTo(5);
  const f = moments(free);
  close(f.mean, k0 * free.t, 1e-6, '⟨x⟩ = ħk₀t/m');
  close(f.sd, s0 * Math.sqrt(1 + (free.t / (2 * s0 * s0)) ** 2), 1e-6, 'σ(t) = σ₀√(1 + (ħt/2mσ₀²)²)');

  // harmonic trap ω = 1: a displaced ground state oscillates rigidly, ⟨x⟩ = x₀ cos t
  const V = Float64Array.from(x, xi => 0.5 * xi * xi);
  const g = Float64Array.from(x, xi => Math.exp(-((xi - 3) ** 2) / 2));
  const osc = schrodingerND({ re0: g, V, shape: [n], spacing: dx, dt: 0.005 });
  osc.advanceTo(Math.PI);
  close(moments(osc).mean, 3 * Math.cos(osc.t), 1e-4, 'half a period later on the other side');
});

test('stationary states: harmonic oscillator levels (n + ½)ħω converge at O(dx²)', () => {
  const levels = dx => {
    const n = Math.round(20 / dx), x = grid(n, -10 + dx, dx);
    return stationaryStates1D({ V: Float64Array.from(x, xi => 0.5 * xi * xi), dx, count: 4 });
  };
  const a = levels(0.1), b = levels(0.05);
  b.energies.forEach((E, i) => close(E, i + 0.5, 2e-3 * (i + 1), `level ${i}`));
  const errA = Math.abs(a.energies[0] - 0.5), errB = Math.abs(b.energies[0] - 0.5);
  close(Math.log2(errA / errB), 2, 0.1, 'order');
  // states are normalised and orthogonal
  const [s0, s1] = b.states;
  close(s0.reduce((s, v) => s + v * v, 0) * 0.05, 1, 1e-12);
  close(s0.reduce((s, v, i) => s + v * s1[i], 0) * 0.05, 0, 1e-10);
});

test('diffusion walls: held values give the exact linear steady state; erfc for a held wall converges at 2nd order', async () => {
  const { erf } = await import('../core/numerics.mjs');
  const n = 50, dx = 1 / n;
  const st = diffusion1D({ u0: new Float64Array(n), D: 1, dx, dt: 0.002, left: { type: 'value', value: 2 }, right: { type: 'value', value: -1 } });
  st.advanceTo(5);
  st.u.forEach((v, i) => close(v, 2 - 3 * (i + 0.5) * dx, 1e-10, `cell ${i}`));

  // semi-infinite medium, wall held at u = 1 from t = 0: u = erfc(x / 2√(Dt)) with x from the wall face
  const D = 0.5, T = 1;
  const err = h => {
    const m = Math.round(8 / h), s = diffusion1D({ u0: new Float64Array(m), D, dx: h, dt: h / 4, left: { type: 'value', value: 1 }, right: 'neumann' });
    s.advanceTo(T);
    let e = 0;
    for (let i = 0; i < m; i++) e = Math.max(e, Math.abs(s.u[i] - (1 - erf((i + 0.5) * h / (2 * Math.sqrt(D * s.t))))));
    return e;
  };
  const e1 = err(0.05), e2 = err(0.025);
  assert.ok(e2 < 2e-3, `max error ${e2}`);
  close(Math.log2(e1 / e2), 2, 0.3, 'order');
});

test('diffusion walls: inflow accounts for every change of the total; a wall can switch mid-run', () => {
  const rand = rng(3), n = 40, dx = 0.1;
  const u0 = Float64Array.from({ length: n }, () => rand());
  const s = diffusion1D({ u0, D: 0.3, dx, dt: 0.02, left: { type: 'value', value: t => 1 + Math.sin(t) }, right: { type: 'flux', value: t => 0.2 * Math.cos(3 * t) } });
  const total0 = s.total();
  s.advanceTo(3);
  close(s.total() - total0, s.inflow, 1e-12, 'Σu dx − Σu₀ dx = inflow');
  assert.ok(Math.abs(s.inflow) > 0.1, 'walls really exchange material');
  s.setBoundary('left', 'neumann'); s.setBoundary('right', { type: 'flux', value: 0 });
  const inflowAtSwitch = s.inflow, totalAtSwitch = s.total();
  s.advanceTo(6);
  close(s.total(), totalAtSwitch, 1e-12, 'insulated after the switch');
  close(s.inflow, inflowAtSwitch, 0);
  assert.throws(() => s.setBoundary('top', 'neumann'), /side/);
});

test('Schrödinger absorbing layers: outgoing packet leaves (≤ 10⁻³ reflected), nothing wraps, interior untouched before arrival', () => {
  const n = 1024, dx = 0.1, x = grid(n, -51.2, dx), k0 = 3;
  const re = Float64Array.from(x, xi => Math.exp(-xi * xi / 4) * Math.cos(k0 * xi));
  const im = Float64Array.from(x, xi => Math.exp(-xi * xi / 4) * Math.sin(k0 * xi));
  const plain = schrodingerND({ re0: re, im0: im, shape: [n], spacing: dx, dt: 0.005 });
  const absorbing = schrodingerND({ re0: re, im0: im, shape: [n], spacing: dx, dt: 0.005, absorb: { width: 120, strength: 5 } });
  const N0 = plain.norm();
  plain.advanceTo(5); absorbing.advanceTo(5);                // packet centre at x = 15, layer starts at x = 39
  let diff = 0;
  for (let i = 0; i < n; i++) diff = Math.max(diff, Math.hypot(plain.re[i] - absorbing.re[i], plain.im[i] - absorbing.im[i]));
  assert.ok(diff < 1e-9, `interior unchanged before arrival (${diff})`);
  absorbing.advanceTo(30); plain.advanceTo(30);              // the plain packet has wrapped round to the left half
  const left = sim => { const p = sim.probability(); let s = 0; for (let i = 0; i < n / 2; i++) s += p[i] * dx; return s; };
  assert.ok(left(plain) > 0.3 * N0, 'without absorption the packet wraps round');
  assert.ok(absorbing.norm() < 1e-3 * N0, `left in the box: ${absorbing.norm() / N0}`);
  close(absorbing.absorbed() + absorbing.norm(), N0, 1e-12 * N0, 'absorbed + remaining = initial');
  assert.throws(() => schrodingerND({ re0: re, shape: [n], spacing: dx, dt: 0.01, absorb: { width: 600, strength: 1 } }), /half/);
});

test('stationaryStates1D on 10⁴ points: fast, residual small, a tunnelling doublet stays orthogonal and has the right parity', async () => {

  const n = 10000, L = 12, dx = L / n, xs = Array.from({ length: n }, (_, i) => -L / 2 + (i + 0.5) * dx);
  const V = Float64Array.from(xs, x => 6 * (x * x - 1) ** 2);               // symmetric double well, barrier 6
  const t0 = performance.now();
  const { energies, states } = stationaryStates1D({ V, dx, count: 4 });
  assert.ok(performance.now() - t0 < 500, `fast: ${performance.now() - t0} ms`);
  assert.ok(energies[1] - energies[0] > 0 && energies[1] < 6, `doublet below the barrier, splitting ${energies[1] - energies[0]}`);
  let ov = 0; for (let i = 0; i < n; i++) ov += states[0][i] * states[1][i] * dx;
  assert.ok(Math.abs(ov) < 1e-8, `overlap ${ov}`);
  // ground state even, first excited odd
  assert.ok(Math.abs(states[0][100] - states[0][n - 101]) < 1e-6 && Math.abs(states[1][100] + states[1][n - 101]) < 1e-6);
  // residual ‖Hψ − Eψ‖ for the ground state (finite-difference H)
  const tk = 1 / (2 * dx * dx); let r = 0;
  for (let i = 1; i < n - 1; i++) { const h = -tk * (states[0][i - 1] - 2 * states[0][i] + states[0][i + 1]) + V[i] * states[0][i]; r = Math.max(r, Math.abs(h - energies[0] * states[0][i])); }
  assert.ok(r < 1e-6, `residual ${r}`);
});
