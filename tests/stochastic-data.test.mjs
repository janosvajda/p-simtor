import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gillespie, eulerMaruyama, langevinBAOAB, metropolis, autocorrelationTime, effectiveSampleSize } from '../core/stochastic.mjs';
import { interpolate, histogram, ksTest } from '../core/data.mjs';
import { rng, erf } from '../core/numerics.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const normalCdf = x => 0.5 * (1 + erf(x / Math.SQRT2));

test('Gillespie: radioactive decay ensemble mean is N₀e^(−kt) within binomial error', () => {
  const k = 0.3, N0 = 200, t = 2, runs = 400, rand = rng(11);
  let sum = 0;
  for (let r = 0; r < runs; r++) {
    const s = gillespie({ x0: [N0], reactions: [{ rate: x => k * x[0], change: [-1] }], rand });
    s.advanceTo(t);
    sum += s.x[0];
  }
  const mean = sum / runs, expected = N0 * Math.exp(-k * t);
  const sd = Math.sqrt(N0 * Math.exp(-k * t) * (1 - Math.exp(-k * t)) / runs);   // binomial
  close(mean, expected, 4 * sd, 'ensemble mean');
});

test('Gillespie: birth–death process reaches the Poisson stationary mean λ/μ', () => {
  const lambda = 20, mu = 1, rand = rng(5);
  const s = gillespie({ x0: [0], reactions: [{ rate: () => lambda, change: [1] }, { rate: x => mu * x[0], change: [-1] }], rand });
  s.advanceTo(10);
  const xs = [];
  for (let i = 0; i < 4000; i++) { s.advanceTo(s.t + 0.25); xs.push(s.x[0]); }
  const { mean, error } = effectiveSampleSize(xs);
  close(mean, lambda / mu, 4 * error, 'stationary mean');
});

test('Euler–Maruyama: Ornstein–Uhlenbeck stationary variance σ²/(2θ)', () => {
  const theta = 2, sigma = 1.5, rand = rng(3);
  const s = eulerMaruyama({ f: x => [-theta * x[0]], g: () => [sigma], x0: [0], dt: 0.002, rand });
  s.advanceTo(5);
  const xs = [];
  for (let i = 0; i < 20000; i++) { s.advanceTo(s.t + 0.05); xs.push(s.x[0] ** 2); }
  const { mean, error } = effectiveSampleSize(xs);
  close(mean, sigma * sigma / (2 * theta), 4 * error + 0.01, '⟨x²⟩ (plus O(dt) bias)');
});

test('Langevin BAOAB: harmonic oscillator samples ⟨x²⟩ = kT/k and ⟨v²⟩ = kT/m', () => {
  const kspring = 4, m = 2, kT = 1.3, rand = rng(9);
  const s = langevinBAOAB({ x0: [1], force: (x, F) => { F[0] = -kspring * x[0]; }, mass: m, gamma: 1, kT, dt: 0.05, rand });
  s.advanceTo(20);
  const x2 = [], v2 = [];
  for (let i = 0; i < 40000; i++) { s.step(); x2.push(s.x[0] ** 2); v2.push(s.v[0] ** 2); }
  const X = effectiveSampleSize(x2), V = effectiveSampleSize(v2);
  close(X.mean, kT / kspring, 4 * X.error, '⟨x²⟩');
  close(V.mean, kT / m, 4 * V.error + 0.01, '⟨v²⟩ (O(dt²) bias)');
});

test('Metropolis: samples a 2D Gaussian; KS test accepts the marginal; autocorrelation is honest', () => {
  const rand = rng(21);
  const s = metropolis({ logP: x => -0.5 * (x[0] ** 2 + (x[1] / 2) ** 2), x0: [0, 0], rand, stepSize: 1.2 });
  for (let i = 0; i < 2000; i++) s.step();
  const xs = [], ys = [];
  for (let i = 0; i < 60000; i++) { s.step(); xs.push(s.x[0]); ys.push(s.x[1] ** 2); }
  const X = effectiveSampleSize(xs), Y = effectiveSampleSize(ys);
  close(X.mean, 0, 4 * X.error);
  close(Y.mean, 4, 4 * Y.error, '⟨y²⟩ = σ_y²');
  assert.ok(X.tau > 1.5, 'Metropolis samples are correlated, and the estimator says so');
  const thinned = xs.filter((_, i) => i % Math.ceil(5 * X.tau) === 0);
  assert.ok(ksTest(thinned, normalCdf).p > 0.001, 'marginal is N(0, 1)');
  assert.ok(s.acceptance() > 0.2 && s.acceptance() < 0.9);
});

test('autocorrelation time: 1 for independent draws, (1+ρ)/(1−ρ) for an AR(1) chain', () => {
  const rand = rng(1);
  close(autocorrelationTime(Array.from({ length: 20000 }, () => rand.normal())), 1, 0.15);
  const rho = 0.8, chain = [0];
  for (let i = 1; i < 100000; i++) chain.push(rho * chain[i - 1] + Math.sqrt(1 - rho * rho) * rand.normal());
  close(autocorrelationTime(chain), (1 + rho) / (1 - rho), 1.0);
});

test('interpolation: hits the data, exact for straight lines, monotone data stay monotone', () => {
  const xs = [0, 1, 2, 3, 4, 5], ys = [0, 0.1, 0.2, 5, 5.1, 5.2];
  for (const method of ['linear', 'monotone']) {
    const f = interpolate(xs, ys, { method });
    xs.forEach((x, i) => close(f(x), ys[i], 1e-15));
    let prev = -Infinity;
    for (let x = 0; x <= 5; x += 0.01) { assert.ok(f(x) >= prev - 1e-12, `${method} monotone at ${x}`); prev = f(x); }
  }
  const line = interpolate([0, 2, 5], [1, 5, 11]);
  close(line(3.7), 1 + 2 * 3.7, 1e-12);
  assert.ok(Number.isNaN(line(6)) && Math.abs(interpolate([0, 1], [0, 1], { extrapolate: true })(2) - 2) < 1e-15);
});

test('histogram: density integrates to 1; KS test rejects the wrong model', () => {
  const rand = rng(4), xs = Array.from({ length: 5000 }, () => rand.normal());
  const h = histogram(xs, { bins: 40, range: [-5, 5] });
  close(h.values.reduce((s, v) => s + v * 0.25, 0), h.inRange / xs.length, 1e-12);
  assert.ok(ksTest(xs, normalCdf).p > 0.01);
  assert.ok(ksTest(xs, x => normalCdf(x / 1.3)).p < 1e-6, 'wrong width is detected');
});

test('metropolisSites: 4×4 periodic Ising matches exact enumeration of all 2¹⁶ states', async () => {
  const { metropolisSites } = await import('../core/stochastic.mjs');
  const L = 4, n = L * L, beta = 0.5;                     // kT/J = 2
  const nb = i => { const x = i % L, y = (i / L) | 0; return [((x + 1) % L) + L * y, ((x + L - 1) % L) + L * y, x + L * ((y + 1) % L), x + L * ((y + L - 1) % L)]; };
  const energy = st => { let E = 0; for (let i = 0; i < n; i++) E -= st[i] * (st[nb(i)[0]] + st[nb(i)[2]]); return E; };
  let Z = 0, sumE = 0, sumM = 0;
  const s = new Int8Array(n);
  for (let c = 0; c < 1 << n; c++) {
    for (let i = 0; i < n; i++) s[i] = (c >> i) & 1 ? 1 : -1;
    const E = energy(s), w = Math.exp(-beta * E), m = Math.abs(s.reduce((a, b) => a + b, 0)) / n;
    Z += w; sumE += w * E; sumM += w * m;
  }
  const exactE = sumE / Z / n, exactM = sumM / Z;
  const state = new Int8Array(n).fill(1), rand = rng(17);
  const mc = metropolisSites({
    state, rand, propose: i => -state[i],
    deltaLogP: (i, value) => { const h = nb(i).reduce((a, j) => a + state[j], 0); return -beta * (-(value - state[i]) * h); },
  });
  for (let k = 0; k < 2000; k++) mc.step();
  const es = [], ms = [];
  for (let k = 0; k < 60000; k++) { mc.step(); es.push(energy(state) / n); ms.push(Math.abs(state.reduce((a, b) => a + b, 0)) / n); }
  const E = effectiveSampleSize(es), M = effectiveSampleSize(ms);
  close(E.mean, exactE, 4 * E.error, '⟨E⟩/N');
  close(M.mean, exactM, 4 * M.error, '⟨|m|⟩');
  assert.ok(E.error < 0.01);
});

test('langevinBAOAB: H(t) − H(0) equals the heat from the bath, with O(dt²) error', () => {
  const k = 3, m = 1.5;
  const imbalance = dt => {
    const s = langevinBAOAB({ x0: [1.2], v0: [0.3], force: (x, F) => { F[0] = -k * x[0]; }, mass: m, gamma: 0.7, kT: 0.8, dt, rand: rng(5) });
    const H = () => 0.5 * m * s.v[0] ** 2 + 0.5 * k * s.x[0] ** 2, H0 = H();
    let worst = 0;
    while (s.t < 20 - 1e-9) { s.step(); worst = Math.max(worst, Math.abs(H() - H0 - s.heat)); }
    return worst;
  };
  const a = imbalance(0.02), b = imbalance(0.01);
  assert.ok(b < 5e-3, `imbalance ${b}`);
  assert.ok(a / b > 2.5, `shrinks with dt (${a} → ${b})`);
});
