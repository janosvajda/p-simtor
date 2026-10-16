// Independent regressions from the scientific review. Numerical reference values were
// evaluated with mpmath at 65 decimal digits for the exact input doubles. Exact identities
// and analytic limits are used wherever possible; no production dependencies are added.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import * as ad from '../core/autodiff.mjs';
import { smallEigenvalues, classifyFixedPoint } from '../core/chaos.mjs';
import { eigenNonsymmetric, etdrk4Periodic, kdv, nls, sineGordon, kuramotoSivashinsky } from '../core/spectral.mjs';
import { stiffStepper, stiffnessRatio } from '../core/stiff.mjs';
import { sahaHydrogen, sahaRatio, halfIonisationTemperature } from '../core/stars.mjs';
import { qmcIntegrate } from '../core/qmc.mjs';
import { grating, fresnelCS, fresnelCoefficients, thinFilm } from '../core/optics.mjs';
import { levenbergMarquardt } from '../core/optimize.mjs';
import { airy, besselK, beta, oscillatorWavefunction, hyp2f1, hyp1f1, zeta, studentTPValue, polylog, expintEi } from '../core/special.mjs';
import { cosmology, soundHorizon, comovingDistance, lookbackTime, hubbleDistance, hubbleTime } from '../core/cosmology.mjs';
import { STA } from '../core/geometric-algebra.mjs';
import { torqueFreeExact } from '../core/rigid.mjs';

const close = (actual, expected, tolerance) => assert.ok(Number.isFinite(actual) && Math.abs(actual - expected) <= tolerance, `${actual} versus ${expected}, tolerance ${tolerance}`);
const relative = (actual, expected, tolerance = 1e-13) => close(actual, expected, Math.abs(expected) * tolerance);

test('Saha fractions remain finite when the Boltzmann factor underflows, retaining representable square roots', () => {
  assert.deepEqual(sahaHydrogen(100, 1e20), { x: 0, ne: 0 });
  assert.ok(sahaHydrogen(200, 1e20).x > 1e-170);
  for (const T of [200, 300, 1000, 5000, 1e5]) {
    const { x, ne } = sahaHydrogen(T, 1e20);
    assert.ok(x >= 0 && x <= 1 && Number.isFinite(ne));
    // In the representable range, satisfy the independently rearranged Saha equation.
    if (T >= 1000 && T <= 5000) relative(x * x / (1 - x), sahaRatio(T, 1e20), 5e-14);
  }
  assert.throws(() => sahaHydrogen(Infinity, 1e20), /finite/);
  const coldHalf = halfIonisationTemperature(1e-300);
  assert.ok(coldHalf < 300); close(sahaHydrogen(coldHalf, 1e-300).x, 0.5, 1e-12);
});

test('QMC refuses identical deterministic replicates instead of assigning a false zero error bar', () => {
  for (const sequence of ['sobol', 'halton']) assert.throws(() => qmcIntegrate(x => x[0] ** 2, 1, { n: 8, replicates: 4, sequence, scramble: 'none' }), /independent randomisations/);
  assert.throws(() => qmcIntegrate(() => 1, 0, { sequence: 'random' }), /dimension/);
  assert.throws(() => qmcIntegrate(() => 1, 1, { n: 2 ** 32 + 1 }), /2\^m/);
});

test('automatic differentiation retains exponentially small tanh derivatives and their second derivatives', () => {
  for (const x of [-20, 20, 350]) {
    const expected = 4 * Math.exp(-2 * Math.abs(x)) / (1 + Math.exp(-2 * Math.abs(x))) ** 2;
    relative(ad.derivative(ad.tanh, x).derivative, expected, 1e-14);
    relative(ad.gradient(([v]) => ad.tanh(v), [x]).gradient[0], expected, 1e-14);
    relative(ad.secondDerivative(ad.tanh, x).second, -2 * Math.tanh(x) * expected, 1e-14);
  }
});

test('atan2 derivatives are invariant to very large and very small coordinate units', () => {
  for (const scale of [1e-200, 1, 1e200]) {
    const g = ad.gradient(([x, y]) => ad.atan2(y, x), [scale, scale]).gradient;
    relative(g[0], -0.5 / scale); relative(g[1], 0.5 / scale);
    relative(ad.derivative(x => ad.atan2(scale, x), scale).derivative, -0.5 / scale);
  }
  const H = ad.hessian(([x, y]) => ad.atan2(y, x), [2, 1]).hessian;
  close(H[0][0], 4 / 25, 1e-15); close(H[0][1], -3 / 25, 1e-15); close(H[1][1], -4 / 25, 1e-15);
});

test('eigenvalue and stability diagnostics preserve repeated spectra instead of roots of a rounded polynomial', () => {
  for (const n of [3, 4, 8]) {
    const A = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_col, j) => (i === j ? 1 : 0)));
    for (const e of smallEigenvalues(A)) { assert.equal(e.re, 1); assert.equal(e.im, 0); }
    assert.equal(classifyFixedPoint(A, { kind: 'map' }).hyperbolic, false);
  }
  for (const scale of [1e-200, 1e200]) {
    const A = [[-2 * scale, -3 * scale], [3 * scale, -2 * scale]];
    for (const e of eigenNonsymmetric(A)) { relative(e.re, -2 * scale); relative(Math.abs(e.im), 3 * scale); }
    relative(stiffnessRatio(A).fastest, 2 * scale);
  }
});

test('all stiff methods safely accept callbacks that reuse derivative and Jacobian work arrays', () => {
  for (const method of ['radau5', 'rosenbrock', 'bdf']) {
    const run = reuse => {
      const out = new Float64Array(2), J = [new Float64Array(2), new Float64Array(2)];
      const s = stiffStepper({ method, y0: [1, 2], h: 0.01, rtol: 1e-9, atol: 1e-12,
        f: (_t, y) => { out[0] = -y[0]; out[1] = -2 * y[1]; return reuse ? out : Array.from(out); },
        jac: () => { J[0][0] = -1; J[1][1] = -2; return reuse ? J : J.map(r => Array.from(r)); } });
      assert.equal(s.advanceTo(1, 10000, { exact: true }), true);
      return s.y;
    };
    const a = run(false), b = run(true);
    assert.deepEqual(a, b);
    relative(b[0], Math.exp(-1), 2e-8); relative(b[1], 2 * Math.exp(-2), 2e-8);
  }
});

test('finite-difference stiff Jacobians remain meaningful for large dimensional state values', () => {
  const r = stiffnessRatio({ f: (_t, y) => [-3 * y[0]], y: [1e20] });
  relative(r.fastest, 3, 2e-8);
});

test('a fixed-step stiff integration still respects hmax and rejects invalid progress limits', () => {
  for (const method of ['radau5', 'rosenbrock', 'bdf']) {
    const s = stiffStepper({ f: (_t, y) => [-y[0]], y0: [1], method, h: 0.2, hmax: 0.01, adaptive: false });
    s.step(); close(s.t, 0.01, 1e-16);
    assert.throws(() => s.step(-1), /positive/);
  }
});

test('large gratings keep nonnegative intensities and the correct resolved peak shape', () => {
  const N = 1e10, base = { a: 0.01, d: 1, lambda: 1, N };
  relative(grating(Math.asin(0.5 / N), base), 4 / Math.PI ** 2, 2e-15);
  const zero = grating(Math.asin(2 / N), base);
  assert.ok(zero >= 0 && zero < 1e-28);
});

test('Fresnel phases retain their fractional square at large x, and huge finite x reaches the correct limit', () => {
  for (const [x, C, S] of [[1000, 0.4999999998986788, 0.4996816901138163], [123456.789, 0.4999992398522264, 0.5000024637082502], [1e10, 0.5, 0.499999999968169]]) {
    const f = fresnelCS(x); close(f.C, C, 2e-16); close(f.S, S, 2e-16);
  }
  assert.deepEqual(fresnelCS(1e200), { C: 0.5, S: 0.5 });
});

test('matched optical media transmit at grazing incidence and opaque films reach the semi-infinite limit', () => {
  const f = fresnelCoefficients(Math.PI / 2, 1.5, 1.5);
  assert.equal(f.Rs, 0); assert.equal(f.Rp, 0); assert.equal(f.Ts, 1); assert.equal(f.Tp, 1);
  for (const polarization of ['s', 'p']) {
    const film = thinFilm({ lambda: 5e-7, layers: [{ n: { re: 1.5, im: 1 }, d: 1e-3 }], ns: 1.5, polarization });
    relative(film.R, (0.5 ** 2 + 1) / (2.5 ** 2 + 1));
    assert.equal(film.T, 0); close(film.R + film.A, 1, 1e-15);
  }
  // An opaque lossless Bragg stack must also remain finite (ordinary transfer matrices overflow).
  const layers = Array.from({ length: 2400 }, (_, i) => ({ n: i % 2 ? 1.5 : 3, d: 1e-6 / (4 * (i % 2 ? 1.5 : 3)) }));
  const mirror = thinFilm({ lambda: 1e-6, layers, ns: 1 });
  close(mirror.R, 1, 1e-14); close(mirror.T, 0, 1e-28);
});

test('soft-L1 fitting retains residuals below sqrt(eps), and fixed/narrow bounds give a valid fit', () => {
  const fit = levenbergMarquardt(p => [p[0] - 1e-8], [0], { loss: 'soft_l1', jacobian: () => [[1]], gtol: 1e-20, xtol: 1e-16 });
  assert.equal(fit.converged, true); relative(fit.p[0], 1e-8, 1e-12);
  const bounded = levenbergMarquardt(p => [p[0] - 5e-9, p[1] - 3], [1e-9, 2], { lower: [0, 2], upper: [1e-8, 2], gtol: 1e-20, xtol: 1e-16, ftol: 0 });
  assert.equal(bounded.converged, true); relative(bounded.p[0], 5e-9, 1e-12);
  assert.equal(bounded.p[1], 2); assert.equal(bounded.active[1], true);
});

test('terminating hypergeometrics survive extreme polynomial cancellation with correct signs and magnitudes', () => {
  // Independent 65-digit evaluations, not the same floating-point term recurrence.
  for (const [a, b, c, x, expected] of [[-20, 10, 1.5, 0.9, -3.8475428234226724e-7], [-50, 20, 2.5, 0.5, 2.4533991657697238e-8], [-50, 20, 2.5, 0.99, 3.6286193032027904e-18]]) relative(hyp2f1(a, b, c, x), expected, 5e-16);
  relative(hyp1f1(-50, 3.5, 20), -0.6147655167329401, 5e-16);
  relative(hyp1f1(-50, 3.5, 100), 1.915357500425757e16, 5e-16);
  relative(hyp2f1(100, 50, 300, 1), 7620420027.851184, 5e-13);
  relative(beta(1e-200, 1e-200), 2e200, 2e-13);
});

test('special functions retain representable tails through intermediate exponential underflow/overflow', () => {
  relative(oscillatorWavefunction(100, 40), 1.0447593187323133e-252);
  relative(oscillatorWavefunction(2000, 80), 8.344467565167805e-234, 2e-12);
  close(oscillatorWavefunction(20, 40), 1.9e-322, Number.MIN_VALUE);
  // K_n references from the independent integral ∫ exp(−x cosh t) cosh(nt) dt at 65 digits.
  relative(besselK(300, 750), Number('4.45979832637210552946509706393105488e-302'), 2e-13);
  relative(besselK(1000, 800), 2.187306658024086e-103, 2e-13);
  relative(besselK(1000, 1000), 1.298180251466701e-233, 2e-13);
  relative(expintEi(716), Number('1.2605029106040893555309663378820853e308'));
  relative(airy(104.3).bi, 4.472500738060502e307, 2e-13);
  close(airy(107).aip, Number('-3.173718939782755292290276e-321'), Number.MIN_VALUE);
});

test('zeta is continuous at zero, Student tails do not disappear on squaring large statistics, and unreliable polylog expansions reject', () => {
  relative(zeta(-1e-10), -0.4999999999081062, 3e-16);
  relative(zeta(-1e-16), -0.4999999999999999, 3e-16);
  relative(studentTPValue(1e160, 1), 2 / Math.PI / 1e160);
  relative(studentTPValue(1e160, 0.5), 6.414019508284458e-81);
  close(studentTPValue(1e160, 2), 1e-320, Number.MIN_VALUE);
  assert.throws(() => polylog(2.000000001, 0.9999999999), /too close to an integer/);
});

test('dealiasing projects the initial state before nonlinear products can alias discarded modes', () => {
  const n = 16, dx = 2 * Math.PI / n;
  const s = etdrk4Periodic({ re0: Float64Array.from({ length: n }, (_, j) => Math.cos(6 * j * dx)), dx, dt: 0.001, linear: () => [0, 0], nonlinear: (r, _i, out) => r.forEach((v, j) => { out[j] = v * v; }), dealias: true });
  assert.ok(Math.max(...s.re.map(Math.abs)) < 1e-14);
  s.step(); assert.ok(Math.max(...s.re.map(Math.abs)) < 1e-14);
});

test('real KdV Nyquist data stay real; NLS kinetic energy and momentum use every propagated Fourier mode', () => {
  const n = 16, L = 2 * Math.PI, dx = L / n, u = Float64Array.from({ length: n }, (_, j) => (j % 2 ? -1 : 1));
  const k = kdv({ u0: u, dx, dt: 0.001 }); k.step();
  assert.ok(Math.max(...k.im.map(Math.abs)) < 1e-14);
  assert.deepEqual(k.u, u);
  const q = nls({ re0: u, dx, dt: 0.001, kappa: 0 });
  relative(q.hamiltonian(), 0.5 * 8 ** 2 * L); relative(q.momentum(), -8 * L);
  q.advanceTo(0.01); relative(q.hamiltonian(), 0.5 * 8 ** 2 * L);
});

test('real-field spectral energies and KS balance include the Nyquist gradient energy', () => {
  const n = 16, L = 2 * Math.PI, dx = L / n, u = Float64Array.from({ length: n }, (_, j) => (j % 2 ? -0.1 : 0.1));
  const s = sineGordon({ u0: u, dx, dt: 0.001 });
  relative(s.energy(), L * (0.5 * 64 * 0.01 + 1 - Math.cos(0.1)));
  const k = kuramotoSivashinsky({ u0: u, dx, dt: 0.0001 });
  relative(k.dissipation(), L * 0.01 * (64 - 4096));
});

test('cosmological inputs enforce the baryon budget and a sound horizon requires photons', () => {
  assert.throws(() => cosmology({ H0: 70, Om: 0.3, flat: true, Obh2: 1 }), /baryon density/);
  const noPhotons = cosmology({ H0: 70, Om: 0.3, flat: true, Tcmb: 0, Obh2: 0.02 });
  assert.throws(() => soundHorizon(1000, noPhotons), /photon density/);
});

test('cosmological distances and lookback times retain tiny nonzero redshifts', () => {
  const p = cosmology({ H0: 70, Om: 1, Or: 0, OL: 0 });
  for (const z of [1e-20, 1e-12, -1e-12, 0.01]) {
    relative(comovingDistance(z, p), -2 * Math.expm1(-0.5 * Math.log1p(z)) * hubbleDistance(p), 5e-15);
    relative(lookbackTime(z, p), -2 / 3 * Math.expm1(-1.5 * Math.log1p(z)) * hubbleTime(p), 5e-15);
  }
});

test('a nearly null multivector inverse does not amplify a neglected nonscalar remainder', () => {
  const a = STA.add(STA.vector([1, 1 - 1e-6, 0, 0]), STA.scalar(4e-15));
  const residual = STA.sub(STA.mul(a, STA.inverse(a)), STA.scalar(1));
  assert.ok(Math.max(...residual.map(Math.abs)) < 1e-12);
  for (const scale of [1e-200, 1e-155, 1e155, 1e200]) {
    const v = STA.vector([scale, 0, 0, 0]), inverse = STA.inverse(v);
    relative(inverse[1], 1 / scale);
    assert.ok(Math.max(...STA.sub(STA.mul(v, inverse), STA.scalar(1)).map(Math.abs)) < 1e-14);
  }
});

test('symmetric torque-free tops include stationary transverse states and arbitrarily slow precession', () => {
  for (const [inertia, L0] of [[[1, 1, 2], [0.4, 0.3, 0]], [[1, 2, 2], [0, 0.3, 0.4]]]) {
    const s = torqueFreeExact({ inertia, L0 });
    assert.equal(s.period, Infinity); assert.deepEqual(Array.from(s.at(10)), L0);
  }
  const s = torqueFreeExact({ inertia: [1, 1, 2], L0: [1, 0, 1e-12] });
  relative(s.period, 4 * Math.PI * 1e12);
  const L = s.at(s.period / 4); close(L[0], 0, 1e-15); close(L[1], 1, 1e-15);
});
