import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  chebyshevNodes, chebyshevDiffMatrix, chebyshevInterpolant, chebyshevCoefficients, clenshawCurtisWeights, clenshawCurtis,
  chebyshevBVP, chebyshevEigen, eigenNonsymmetric, fourierDerivative, etdrk4Periodic, kdv, kuramotoSivashinsky, viscousBurgers,
  nls, sineGordon, nlsSoliton, nlsPeregrine, kdvSoliton, kdvTwoSoliton, kdvPhaseShifts, burgersColeHopf, sineGordonKink, sineGordonBreather,
} from '../core/spectral.mjs';
import { besselI } from '../core/numerics.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const grid = (n, x0, dx) => Float64Array.from({ length: n }, (_, i) => x0 + i * dx);
const maxErr = (xs, f) => { let e = 0; xs.forEach((x, i) => { e = Math.max(e, Math.abs(f(x, i))); }); return e; };

test('Chebyshev nodes, differentiation matrix and interpolant: exact for polynomials of degree ≤ n', () => {
  const n = 12, iv = [0, 3], { x, D, D2 } = chebyshevDiffMatrix(n, iv);
  assert.deepEqual(Array.from(x), Array.from(chebyshevNodes(n, iv)));
  close(x[0], 0, 0); close(x[n], 3, 1e-15); close(x[n / 2], 1.5, 1e-15, 'symmetric');
  const p = t => t ** 12 - 2 * t ** 7 + t, dp = t => 12 * t ** 11 - 14 * t ** 6 + 1, d2p = t => 132 * t ** 10 - 84 * t ** 5;
  const u = Float64Array.from(x, p);
  x.forEach((xi, i) => {
    let d1 = 0, d2 = 0; for (let j = 0; j <= n; j++) { d1 += D[i][j] * u[j]; d2 += D2[i][j] * u[j]; }
    close(d1, dp(xi), 1e-12 * 12 * 3 ** 11, 'first derivative');
    close(d2, d2p(xi), 1e-11 * 132 * 3 ** 10, 'second derivative');
  });
  const I = chebyshevInterpolant(u, iv);
  for (const t of [0.1, 1.234, 2.9]) close(I(t), p(t), 1e-13 * 3 ** 12);
  // T_3 has coefficient vector e_3
  const a = chebyshevCoefficients(Float64Array.from(chebyshevNodes(8), t => 4 * t ** 3 - 3 * t));
  a.forEach((v, k) => close(v, k === 3 ? 1 : 0, 1e-14));
});

test('Chebyshev spectral accuracy: derivative of e^x and the Runge function converge geometrically', () => {
  const errs = [8, 16, 24].map(n => {
    const { x, D } = chebyshevDiffMatrix(n), u = Float64Array.from(x, Math.exp);
    return maxErr(x, (xi, i) => D[i].reduce((s, d, j) => s + d * u[j], 0) - Math.exp(xi));
  });
  assert.ok(errs[0] < 1e-4 && errs[1] < 1e-12 && errs[2] < 1e-12, `${errs}`);
  // Runge 1/(1 + 25x²): interpolation error ∝ ρ^−n, ρ = (1 + √26)/5 (Trefethen ATAP ch. 8)
  const runge = t => 1 / (1 + 25 * t * t), rho = (1 + Math.sqrt(26)) / 5;
  const e = [40, 80].map(n => { const I = chebyshevInterpolant(Float64Array.from(chebyshevNodes(n), runge)); return maxErr(grid(301, -1, 2 / 300), t => I(t) - runge(t)); });
  close(Math.log(e[0] / e[1]) / 40, Math.log(rho), 0.03, 'measured convergence rate');
});

test('Clenshaw–Curtis: exact to degree n, positive weights, geometric convergence', () => {
  const { w } = clenshawCurtisWeights(10, [-1, 2]);
  assert.ok(w.every(v => v > 0));
  close(w.reduce((a, b) => a + b), 3, 1e-14);
  close(clenshawCurtis(t => t ** 10, -1, 2, 10), (2 ** 11 + 1) / 11, 1e-12);
  close(clenshawCurtis(Math.exp, -1, 1, 16), 2 * Math.sinh(1), 1e-15);
  const exact = Math.atan(4) / 2, e = [16, 32, 64].map(n => Math.abs(clenshawCurtis(t => 1 / (1 + 16 * t * t), -1, 1, n) - exact));
  assert.ok(e[0] > 1e4 * e[2] && e[2] < 1e-12, `${e}`);
});

test('Chebyshev BVP: u″ = e^{4x} (Trefethen Program 13) converges spectrally; Robin boundary', () => {
  const exact = t => (Math.exp(4 * t) - t * Math.sinh(4) - Math.cosh(4)) / 16;
  const errs = [8, 12, 16, 24].map(n => { const { x, u } = chebyshevBVP({ n, f: t => Math.exp(4 * t) }); return maxErr(x, (xi, i) => u[i] - exact(xi)); });
  for (let i = 1; i < errs.length; i++) assert.ok(errs[i] < 1e-2 * errs[i - 1] || errs[i] < 1e-13, `${errs}`);
  assert.ok(errs[3] < 1e-13);
  const sol = chebyshevBVP({ n: 24, f: t => Math.exp(4 * t) });
  assert.ok(sol.tail < 1e-14, 'coefficient tail flags resolution');
  close(sol.interpolant(0.3), exact(0.3), 1e-13);
  // u″ − u = 0, u(−1) = 1, u′(1) = 0  →  u = cosh(x − 1)/cosh 2
  const r = chebyshevBVP({ n: 24, c: -1, left: { value: 1 }, right: { alpha: 0, beta: 1, value: 0 } });
  assert.ok(maxErr(r.x, (xi, i) => r.u[i] - Math.cosh(xi - 1) / Math.cosh(2)) < 1e-13);
  assert.throws(() => chebyshevBVP({ n: 16, c: (Math.PI / 2) ** 2 }), /singular/, 'u″ + (π/2)²u = 0 has a nullspace');
});

test('Chebyshev eigenproblems: −u″ → (kπ/2)², harmonic oscillator → 2k + 1, spectral convergence', () => {
  const lam = k => (k * Math.PI / 2) ** 2;
  const err10 = [16, 24, 32].map(n => Math.abs(chebyshevEigen({ n, a: -1 }).values[9].re - lam(10)) / lam(10));
  assert.ok(err10[0] > 1e-3 && err10[1] < 1e-4 && err10[2] < 1e-9, `${err10}`);
  const { values, vectors, x } = chebyshevEigen({ n: 48, a: -1, vectors: 2 });
  for (let k = 1; k <= 16; k++) { close(values[k - 1].re, lam(k), 1e-10 * lam(k), `λ_${k}`); assert.equal(values[k - 1].im, 0); }
  // eigenfunctions: cos(πx/2), sin(πx)
  assert.ok(maxErr(x, (xi, i) => vectors[0].u[i] - Math.cos(Math.PI * xi / 2)) < 1e-10);
  const top = Math.max(...x.map(xi => Math.abs(Math.sin(Math.PI * xi))));            // max-norm on the nodes
  assert.ok(maxErr(x, (xi, i) => vectors[1].u[i] + Math.sin(Math.PI * xi) / top) < 1e-10, 'first lobe positive');
  const ho = chebyshevEigen({ n: 64, interval: [-8, 8], a: -1, c: t => t * t });
  for (let k = 0; k < 8; k++) close(ho.values[k].re, 2 * k + 1, 1e-10);
});

test('eigenNonsymmetric: complex pairs and a defective-free companion matrix', () => {
  // rotation-scaling block (eigenvalues 1 ± 2i) ⊕ [3]
  const v = eigenNonsymmetric([[1, -2, 0], [2, 1, 0], [5, 7, 3]]);
  close(v[0].re, 1, 1e-14); close(Math.abs(v[0].im), 2, 1e-14); close(v[0].im + v[1].im, 0, 1e-14); close(v[2].re, 3, 1e-14);
  // companion matrix of (x − 1)(x − 2)(x − 3)(x − 4)(x − 5)
  const c = [-120, 274, -225, 85, -15];                         // x⁵ − 15x⁴ + 85x³ − 225x² + 274x − 120
  const A = Array.from({ length: 5 }, (_, i) => Float64Array.from({ length: 5 }, (__, j) => (j === 4 ? -c[i] : i === j + 1 ? 1 : 0)));
  eigenNonsymmetric(A).forEach((e, k) => { close(e.re, k + 1, 1e-10); close(e.im, 0, 1e-10); });
});

test('Fourier pseudo-spectral derivative: exponential accuracy for exp(sin x)', () => {
  const n = 64, L = 2 * Math.PI, x = grid(n, 0, L / n), u = Float64Array.from(x, t => Math.exp(Math.sin(t)));
  const d1 = fourierDerivative(u, L), d2 = fourierDerivative(u, L, 2), d3 = fourierDerivative(u, L, 3);
  assert.ok(maxErr(x, (t, i) => d1[i] - Math.cos(t) * Math.exp(Math.sin(t))) < 1e-13);
  assert.ok(maxErr(x, (t, i) => d2[i] - (Math.cos(t) ** 2 - Math.sin(t)) * Math.exp(Math.sin(t))) < 1e-12);
  assert.ok(maxErr(x, (t, i) => d3[i] - (Math.cos(t) ** 3 - 3 * Math.sin(t) * Math.cos(t) - Math.cos(t)) * Math.exp(Math.sin(t))) < 1e-11);
  assert.throws(() => fourierDerivative(new Float64Array(10), 1), /power of two/);
});

test('NLS bright soliton: exact profile, norm to round-off, exact H and P, measured order 2 and 4', () => {
  const n = 512, dx = 64 / n, x = grid(n, -32, dx), opt = { amplitude: 1, velocity: 0.5 };
  const re0 = Float64Array.from(x, v => nlsSoliton(v, 0, opt).re), im0 = Float64Array.from(x, v => nlsSoliton(v, 0, opt).im);
  const err = (order, dt) => {
    const s = nls({ re0, im0, dx, dt, order }), N0 = s.norm(), H0 = s.hamiltonian();
    close(N0, 2, 1e-12, 'N = 2A'); close(H0, 0.25 - 1 / 3, 1e-12, 'H = Av² − A³/3'); close(s.momentum(), 1, 1e-12, 'P = 2Av');
    s.advanceTo(5);
    close(s.norm(), N0, 1e-13, 'unitary steps');
    return maxErr(x, (v, i) => { const e = nlsSoliton(v, s.t, opt); return Math.hypot(s.re[i] - e.re, s.im[i] - e.im); });
  };
  const e2 = [err(2, 0.02), err(2, 0.01)], e4 = [err(4, 0.04), err(4, 0.02)];
  close(Math.log2(e2[0] / e2[1]), 2, 0.1, 'Strang order'); close(Math.log2(e4[0] / e4[1]), 4, 0.2, 'Yoshida order');
  assert.ok(e4[1] < 1e-5);
});

test('NLS Peregrine breather: modulational focusing to peak amplitude 3 (exact rational solution)', () => {
  const n = 1024, Lx = 120, dx = Lx / n, x = grid(n, -Lx / 2, dx);
  const re0 = Float64Array.from(x, v => nlsPeregrine(v, -3).re), im0 = Float64Array.from(x, v => nlsPeregrine(v, -3).im);
  const s = nls({ re0, im0, dx, dt: 0.004, order: 4, t0: -3 });
  s.advanceTo(0);
  let peak = 0; for (let i = 0; i < n; i++) peak = Math.max(peak, Math.hypot(s.re[i], s.im[i]));
  close(peak, 3, 1e-4, 'max |ψ| = 3');
  assert.ok(maxErr(x, (v, i) => (Math.abs(v) < 10 ? Math.hypot(s.re[i] - nlsPeregrine(v, s.t).re, s.im[i] - nlsPeregrine(v, s.t).im) : 0)) < 1e-4);
});

test('KdV two-soliton collision: exact Hirota solution, Hirota phase shifts, invariants', () => {
  const n = 256, dx = 80 / n, x = grid(n, -40, dx), p = { k1: 1.5, k2: 1, x1: -20, x2: -8 };
  const s = kdv({ u0: Float64Array.from(x, v => kdvTwoSoliton(v, 0, p)), dx, dt: 0.004 });
  const M0 = s.mass(), P0 = s.momentum(), E0 = s.energy();
  // separated solitons: momentum ∫u² = (2/3)(k1³ + k2³); energy ∫(½u_x² − u³) = −(k1⁵ + k2⁵)/5
  close(P0, (2 / 3) * (1.5 ** 3 + 1), 1e-8); close(E0, -(1.5 ** 5 + 1) / 5, 1e-8); close(M0, 2 * (1.5 + 1), 1e-8);
  // peak positions by Newton on the trigonometric interpolant's derivative
  const peak = (u, lo, hi) => {
    let guess = lo, top = -Infinity;
    x.forEach((v, i) => { if (v >= lo && v <= hi && u[i] > top) { top = u[i]; guess = v; } });
    const L = n * dx, k = j => 2 * Math.PI * (j < n / 2 ? j : j - n) / L;
    const re = Float64Array.from(u), im = new Float64Array(n);
    for (let j = 0; j < n; j++) { let a = 0, b = 0; for (let m = 0; m < n; m++) { a += u[m] * Math.cos(k(j) * x[m]); b -= u[m] * Math.sin(k(j) * x[m]); } re[j] = a / n; im[j] = b / n; }
    let z = guess;
    for (let it = 0; it < 20; it++) {
      let d1 = 0, d2 = 0;
      for (let j = 0; j < n; j++) { if (j === n / 2) continue; const q = k(j), c = Math.cos(q * z), sn = Math.sin(q * z); d1 += q * (-re[j] * sn - im[j] * c); d2 += -q * q * (re[j] * c - im[j] * sn); }
      z -= d1 / d2;
    }
    return z;
  };
  const before = [peak(s.u, -30, -14), peak(s.u, -12, 0)];
  s.advanceTo(20);
  assert.ok(maxErr(x, (v, i) => s.u[i] - kdvTwoSoliton(v, s.t, p)) < 1e-6, 'exact two-soliton solution');
  const after = [peak(s.u, 20, 39), peak(s.u, 0, 19)], shift = kdvPhaseShifts(1.5, 1);
  close(after[0] - before[0] - 1.5 ** 2 * s.t, shift.fast, 2e-4, 'fast soliton jumps forward');
  close(after[1] - before[1] - s.t, shift.slow, 2e-4, 'slow soliton jumps back');
  close(s.mass(), M0, 1e-12); close(s.momentum(), P0, 1e-6); close(s.energy(), E0, 1e-6);
});

test('ETDRK4: fourth order in dt on a single KdV soliton; Courant guard', () => {
  const n = 128, dx = 40 / n, x = grid(n, -20, dx), u0 = Float64Array.from(x, v => kdvSoliton(v, 0, { k: 1.2, x0: -5 }));
  const run = dt => { const s = kdv({ u0, dx, dt }); s.advanceTo(4); return s.u; };
  const a = run(0.04), b = run(0.02), c = run(0.01);
  const d = (p, q) => maxErr(p, (_, i) => p[i] - q[i]);
  close(Math.log2(d(a, b) / d(b, c)), 4, 0.3, 'ETDRK4 order');
  assert.ok(maxErr(x, (v, i) => c[i] - kdvSoliton(v, 4, { k: 1.2, x0: -5 })) < 1e-6);
  assert.throws(() => kdv({ u0, dx, dt: 0.5 }).step(), /Courant/);
  assert.throws(() => etdrk4Periodic({ re0: new Float64Array(12), dx: 1, dt: 1, linear: () => [0, 0], nonlinear: () => {} }), /power of two/);
});

test('Kuramoto–Sivashinsky: exact linear growth k² − k⁴, energy budget, Kassam–Trefethen run', () => {
  const n = 128, L = 32 * Math.PI, dx = L / n, x = grid(n, 0, dx);
  // tiny single mode m = 11: k = 11/16, σ = k² − k⁴ (nonlinearity ~1e−10 relative)
  const kk = 11 / 16, small = kuramotoSivashinsky({ u0: Float64Array.from(x, v => 1e-10 * Math.cos(kk * v)), dx, dt: 0.25 });
  small.advanceTo(20);
  close(small.u[0] / 1e-10, Math.exp((kk ** 2 - kk ** 4) * 20), 1e-6, 'linear growth');
  // dE/dt = ∫(u_x² − u_xx²): Simpson's rule on the rate; the residual must shrink as dt⁴ (needs dealias:
  // the aliased product breaks ∫u·(u²)_x = 0 at the 1e−4 level on this grid)
  const u0 = Float64Array.from(x, v => Math.cos(v / 16) * (1 + Math.sin(v / 16)));
  const budget = dt => {
    const s = kuramotoSivashinsky({ u0, dx, dt, dealias: true }), m0 = s.mean(), E0 = s.energy(), steps = Math.round(20 / dt), rates = [s.dissipation()];
    for (let i = 0; i < steps; i++) { s.step(); rates.push(s.dissipation()); }
    let integral = 0;
    for (let i = 0; i < steps; i += 2) integral += dt / 3 * (rates[i] + 4 * rates[i + 1] + rates[i + 2]);
    close(s.mean(), m0, 1e-15, 'mean conserved');
    return Math.abs(s.energy() - E0 - integral);
  };
  const r1 = budget(0.1), r2 = budget(0.05);
  assert.ok(r2 < r1 / 10 && r2 < 1e-5, `energy-budget residual ${r1} → ${r2}`);
  // Kassam & Trefethen (2005) §6 run: h = 1/4 to t = 150
  const ref = kuramotoSivashinsky({ u0, dx, dt: 0.25 });
  ref.advanceTo(150);
  const amp = Math.max(...ref.u.map(Math.abs));
  assert.ok(amp > 1.5 && amp < 4, `saturated chaotic amplitude ${amp} (KT Fig. 6.1: |u| ≲ 3)`);
});

test('viscous Burgers against the exact Cole–Hopf solution (steep front, ν = 0.05)', () => {
  // Cole–Hopf quadrature against the independent Bessel-series form for u0 = sin x
  const nu0 = 0.5, a = 1 / (2 * nu0), series = (xx, t) => {
    let num = 0, den = besselI(0, a);
    for (let m = 1; m < 40; m++) { const e = besselI(m, a) * Math.exp(-nu0 * m * m * t); num += 4 * nu0 * m * e * Math.sin(m * xx); den += 2 * e * Math.cos(m * xx); }
    return num / den;
  };
  for (const xx of [0.5, 2, 3.5]) close(burgersColeHopf({ primitive: y => 1 - Math.cos(y), nu: nu0, umax: 1 }, xx, 0.7), series(xx, 0.7), 1e-12);
  const nu = 0.05, n = 512, L = 2 * Math.PI, dx = L / n, x = grid(n, 0, dx);
  const s = viscousBurgers({ u0: Float64Array.from(x, Math.sin), nu, dx, dt: 0.005 });
  const E0 = s.energy();
  s.advanceTo(1.5);
  assert.ok(maxErr(x.filter((_, i) => i % 4 === 0), (v, j) => s.u[4 * j] - burgersColeHopf({ primitive: y => 1 - Math.cos(y), nu, umax: 1 }, v, s.t)) < 1e-7);
  close(s.mean(), 0, 1e-15); assert.ok(s.energy() < E0);
});

test('sine-Gordon: moving kink and breather follow the exact solutions; energy 8γ and 16√(1 − ω²)', () => {
  const n = 512, dx = 80 / n, x = grid(n, -40, dx), p = { velocity: 0.5, x0: -10 };
  const kink = (dt, order) => {
    const s = sineGordon({ u0: Float64Array.from(x, v => sineGordonKink(v, 0, p).u), v0: Float64Array.from(x, v => sineGordonKink(v, 0, p).ut), dx, dt, winding: 1, order });
    const E0 = s.energy();
    close(E0, 8 / Math.sqrt(0.75), 1e-12, 'kink rest energy × γ');
    s.advanceTo(10);
    return { e: maxErr(x, (v, i) => s.u[i] - sineGordonKink(v, s.t, p).u), dE: Math.abs(s.energy() - E0) };
  };
  const a = kink(0.04, 2), b = kink(0.02, 2), c = kink(0.04, 4), d = kink(0.02, 4);
  close(Math.log2(a.e / b.e), 2, 0.1, 'Strang'); close(Math.log2(c.e / d.e), 4, 0.25, 'Yoshida');
  assert.ok(d.e < 1e-6 && d.dE < 1e-10);
  const q = { omega: 0.5 }, m = 256, h = 60 / m, y = grid(m, -30, h);
  const s = sineGordon({ u0: Float64Array.from(y, v => sineGordonBreather(v, 0, q).u), v0: Float64Array.from(y, v => sineGordonBreather(v, 0, q).ut), dx: h, dt: 0.02, order: 4 });
  close(s.energy(), 16 * Math.sqrt(0.75), 1e-12);
  s.advanceTo(4 * Math.PI);                                    // one breather period 2π/ω
  assert.ok(maxErr(y, (v, i) => s.u[i] - sineGordonBreather(v, s.t, q).u) < 1e-6);
  assert.throws(() => sineGordon({ u0: new Float64Array(8), dx: 1, dt: 0.1, winding: 0.5 }), /integer/);
});
