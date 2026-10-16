// Minimisers checked on exact minima of published test functions (Moré, Garbow & Hillstrom 1981,
// "Testing unconstrained optimization software", ACM TOMS 7, 17: Rosenbrock #1, Beale #5,
// helical valley #7, Powell singular #13, extended Rosenbrock #21), Himmelblau (1972) and
// Rastrigin, plus measured convergence rates and finite-difference error orders.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  goldenSection, brentMinimize, nelderMead, wolfeLineSearch, bfgs, lbfgs, levenbergMarquardt,
  simulatedAnnealing, differentialEvolution, finiteGradient, finiteHessian, checkGradient,
} from '../core/optimize.mjs';
import { curveFit } from '../core/data.mjs';
import { leastSquares } from '../core/linalg.mjs';
import { rng } from '../core/numerics.mjs';
import * as AD from '../core/autodiff.mjs';

const close = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol, `${a} vs ${b} (tolerance ${tol})`);
const vclose = (a, b, tol) => { assert.equal(a.length, b.length); a.forEach((v, i) => close(v, b[i], tol)); };

const rosen = x => 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2;
const rosenGrad = x => [-400 * x[0] * (x[1] - x[0] ** 2) - 2 * (1 - x[0]), 200 * (x[1] - x[0] ** 2)];
const rosenRes = p => [10 * (p[1] - p[0] ** 2), 1 - p[0]];
const bealeRes = p => [1.5, 2.25, 2.625].map((y, i) => y - p[0] * (1 - p[1] ** (i + 1)));
const sumSq = r => x => r(x).reduce((a, v) => a + v * v, 0);
const powellRes = p => [p[0] + 10 * p[1], Math.sqrt(5) * (p[2] - p[3]), (p[1] - 2 * p[2]) ** 2, Math.sqrt(10) * (p[0] - p[3]) ** 2];
const helicalRes = p => {
  const theta = Math.atan(p[1] / p[0]) / (2 * Math.PI) + (p[0] < 0 ? 0.5 : 0);
  return [10 * (p[2] - 10 * theta), 10 * (Math.hypot(p[0], p[1]) - 1), p[2]];
};
const himmelblau = x => (x[0] ** 2 + x[1] - 11) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2;
const HIMMELBLAU_MINIMA = [[3, 2], [-2.805118086952745, 3.131312518250573], [-3.779310253377747, -3.283185991286170], [3.584428340330492, -1.848126526964404]];
const rastrigin = x => 10 * x.length + x.reduce((a, v) => a + v * v - 10 * Math.cos(2 * Math.PI * v), 0);
const nearestHimmelblau = x => Math.min(...HIMMELBLAU_MINIMA.map(m => Math.hypot(x[0] - m[0], x[1] - m[1])));

test('golden section and Brent locate exact 1D minima; Brent needs far fewer evaluations', () => {
  for (const minimise of [goldenSection, brentMinimize]) {
    const q = minimise(x => (x - 2) ** 2 + 1, 0, 5);
    assert.ok(q.converged); close(q.x, 2, 1e-7); close(q.f, 1, 1e-14);
    const c = minimise(Math.cos, 3, 4);
    close(c.x, Math.PI, 1e-7); close(c.f, -1, 1e-15);
    const k = minimise(x => Math.abs(x - 1 / 3), 0, 1);     // non-smooth: still bracketed
    close(k.x, 1 / 3, 1e-8);
    const edge = minimise(x => x, 0, 1);                       // monotone: converges to the endpoint
    close(edge.x, 0, 1e-9);
  }
  const g = goldenSection(Math.cos, 3, 4), b = brentMinimize(Math.cos, 3, 4);
  assert.ok(b.evaluations * 2 < g.evaluations, `${b.evaluations} vs ${g.evaluations}`);
  // golden section shrinks the bracket by exactly 1/φ per evaluation
  const width = g.bracket[1] - g.bracket[0];
  close(Math.log(width) / (g.evaluations - 1), Math.log((Math.sqrt(5) - 1) / 2), 0.02);
  assert.throws(() => brentMinimize(Math.cos, 1, 1), /a ≠ b/);
});

test('Nelder–Mead reaches Rosenbrock, Himmelblau and Beale minima and respects bounds', () => {
  const r = nelderMead(rosen, [-1.2, 1]);
  assert.ok(r.converged, r.reason); vclose(r.x, [1, 1], 1e-7);
  const h = nelderMead(himmelblau, [1, 1]);
  assert.ok(nearestHimmelblau(h.x) < 1e-7);
  const be = nelderMead(sumSq(bealeRes), [1, 1]);
  vclose(be.x, [3, 0.5], 1e-6);
  const bounded = nelderMead(rosen, [-1.2, 1], { upper: [0.5, 10] });
  vclose(bounded.x, [0.5, 0.25], 1e-6); close(bounded.f, 0.25, 1e-10);
  const n6 = x => x.reduce((a, v, i) => a + (i + 1) * (v - i) ** 2, 0);   // adaptive coefficients work in 6D
  vclose(nelderMead(n6, new Array(6).fill(0)).x, [0, 1, 2, 3, 4, 5], 1e-6);
});

test('strong-Wolfe line search satisfies both Wolfe conditions from long and short trial steps', () => {
  const x = [-1.2, 1], g = rosenGrad(x), d = g.map(v => -v), f0 = rosen(x), slope = -(g[0] ** 2 + g[1] ** 2);
  for (const alpha0 of [1, 1e-6]) {
    const ls = wolfeLineSearch(rosen, rosenGrad, x, d, { alpha0, c1: 1e-4, c2: 0.9 });
    assert.ok(ls.ok);
    assert.ok(ls.f <= f0 + 1e-4 * ls.alpha * slope);
    assert.ok(Math.abs(ls.g[0] * d[0] + ls.g[1] * d[1]) <= 0.9 * Math.abs(slope));
  }
  assert.throws(() => wolfeLineSearch(rosen, rosenGrad, x, g), /descent/);
});

test('BFGS and L-BFGS solve the Moré–Garbow–Hillstrom problems with superlinear convergence', () => {
  const errors = [];
  const b = bfgs(rosen, [-1.2, 1], { gradient: rosenGrad, onIteration: ({ x }) => errors.push(Math.hypot(x[0] - 1, x[1] - 1)) });
  assert.ok(b.converged, b.reason); vclose(b.x, [1, 1], 1e-10);
  // superlinear: the error ratio e_{k+1}/e_k tends to 0 over the last steps
  const tail = errors.filter(e => e > 1e-12 && e < 1e-2);
  const ratios = tail.slice(1).map((e, i) => e / tail[i]);
  assert.ok(Math.min(...ratios.slice(-3)) < 0.05, ratios.join(' '));
  const l = lbfgs(rosen, [-1.2, 1], { gradient: rosenGrad });
  assert.ok(l.converged); vclose(l.x, [1, 1], 1e-10);
  // finite-difference gradient: same minimum to the accuracy the differences allow
  const fd = bfgs(rosen, [-1.2, 1], { gtol: 1e-6 });
  assert.ok(fd.converged); vclose(fd.x, [1, 1], 1e-6);
  vclose(bfgs(sumSq(bealeRes), [1, 1], { gtol: 1e-7 }).x, [3, 0.5], 1e-6);
  vclose(bfgs(sumSq(helicalRes), [-1, 0, 0], { gtol: 1e-7 }).x, [1, 0, 0], 1e-6);
  // Powell singular: Hessian singular at x* (only linear convergence), still f → 0
  assert.ok(lbfgs(sumSq(powellRes), [3, -1, 0, 1], { gtol: 1e-10 }).f < 1e-12);
  // extended Rosenbrock, n = 100
  const n = 100, ext = x => { let s = 0; for (let i = 0; i < n; i += 2) s += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2; return s; };
  const extGrad = x => { const g = new Float64Array(n); for (let i = 0; i < n; i += 2) { const t = x[i + 1] - x[i] ** 2; g[i] = -400 * x[i] * t - 2 * (1 - x[i]); g[i + 1] = 200 * t; } return g; };
  const big = lbfgs(ext, Array.from({ length: n }, (_, i) => (i % 2 ? 1 : -1.2)), { gradient: extGrad });
  assert.ok(big.converged); assert.ok(big.x.every(v => Math.abs(v - 1) < 1e-8));
  // every Himmelblau start goes to one of its four exact minima
  const rand = rng(3);
  for (let k = 0; k < 10; k++) assert.ok(nearestHimmelblau(lbfgs(himmelblau, [10 * rand() - 5, 10 * rand() - 5], { gtol: 1e-9 }).x) < 1e-8);
});

test('projected L-BFGS finds bound-constrained minima and reports active bounds', () => {
  const r = lbfgs(rosen, [-1.2, 1], { gradient: rosenGrad, upper: [0.5, Infinity], gtol: 1e-6 });
  vclose(r.x, [0.5, 0.25], 1e-7); assert.deepEqual(r.active, [true, false]);
  // separable quadratic with the unconstrained minimum outside the box: x* = clip(c)
  // written as Σ (i+1)[(x−c)² − (b−c)²] = Σ (i+1)(x−b)(x+b−2c), b = clip(c), so f* = 0 exactly and
  // f resolves x to full precision (a large constant would hide changes below √ε)
  const c = [-2, 0.3, 5, 1], bc = c.map(v => Math.min(2, Math.max(0, v)));
  const q = x => x.reduce((a, v, i) => a + (i + 1) * (v - bc[i]) * (v + bc[i] - 2 * c[i]), 0);
  const qg = x => Array.from(x, (v, i) => 2 * (i + 1) * (v - c[i]));
  const b = lbfgs(q, [0.5, 0.5, 0.5, 0.5], { gradient: qg, lower: 0, upper: 2, gtol: 1e-12 });
  assert.ok(b.converged); vclose(b.x, [0, 0.3, 2, 1], 1e-10);
  assert.deepEqual(b.active, [true, false, true, false]);
});

test('Levenberg–Marquardt: exact minima, bounds, linear-model covariance, robust losses', () => {
  const r = levenbergMarquardt(rosenRes, [-1.2, 1]);
  assert.ok(r.converged); vclose(r.p, [1, 1], 1e-10);
  vclose(levenbergMarquardt(bealeRes, [1, 1]).p, [3, 0.5], 1e-9);
  assert.ok(levenbergMarquardt(powellRes, [3, -1, 0, 1]).cost < 1e-15);
  vclose(levenbergMarquardt(helicalRes, [-1, 0, 0]).p, [1, 0, 0], 1e-9);
  const bounded = levenbergMarquardt(rosenRes, [-1.2, 1], { upper: [0.5, 10] });
  assert.ok(bounded.converged); vclose(bounded.p, [0.5, 0.25], 1e-8); assert.deepEqual(bounded.active, [true, false]);
  assert.ok(Number.isNaN(bounded.errors[0]));
  // linear model with known σ: estimate and covariance are exactly the weighted LS ones
  const rand = rng(11), sigma = 0.1, xs = Array.from({ length: 30 }, (_, i) => i / 10);
  const ys = xs.map(x => 1.5 - 0.7 * x + sigma * rand.normal());
  const lin = levenbergMarquardt(p => xs.map((x, i) => (ys[i] - p[0] - p[1] * x) / sigma), [0, 0], { absoluteSigma: true });
  const exact = leastSquares(xs.map(x => [1, x]), ys).x;
  vclose(lin.p, exact, 1e-10);
  let s0 = 0, s1 = 0, s2 = 0;
  for (const x of xs) { s0 += 1; s1 += x; s2 += x * x; }
  const det = s0 * s2 - s1 * s1, w = sigma * sigma;
  close(lin.cov[0][0], w * s2 / det, 1e-12); close(lin.cov[1][1], w * s0 / det, 1e-12); close(lin.cov[0][1], -w * s1 / det, 1e-12);
  // agrees with curveFit (scaled covariance) on an exponential decay
  const ts = Array.from({ length: 40 }, (_, i) => i * 0.1), data = ts.map(t => 2 * Math.exp(-1.3 * t) + 0.02 * rand.normal());
  const model = (t, p) => p[0] * Math.exp(-p[1] * t);
  const cf = curveFit(model, ts, data, [1, 1]);
  const lm = levenbergMarquardt(p => ts.map((t, i) => data[i] - model(t, p)), [1, 1]);
  vclose(lm.p, cf.p, 1e-8); vclose(lm.errors, cf.errors, 1e-6 * Math.max(...cf.errors) + 1e-9);
  // robust loss: 20 % gross outliers barely move a soft-L1 fit, but wreck the plain one
  const xo = Array.from({ length: 50 }, (_, i) => i / 10), yo = xo.map((x, i) => 2 + 3 * x + 0.05 * rand.normal() + (i % 5 === 0 ? 30 : 0));
  const res = p => xo.map((x, i) => yo[i] - p[0] - p[1] * x);
  const plain = levenbergMarquardt(res, [0, 0]), robust = levenbergMarquardt(res, [0, 0], { loss: 'soft_l1', fScale: 0.1 });
  assert.ok(Math.abs(plain.p[0] - 2) > 3);
  close(robust.p[0], 2, 0.1); close(robust.p[1], 3, 0.03);
  for (const loss of ['huber', 'cauchy']) close(levenbergMarquardt(res, [0, 0], { loss, fScale: 0.1 }).p[1], 3, 0.03);
  assert.throws(() => levenbergMarquardt(res, [0, 0], { loss: 'l3' }), /unknown loss/);
});

test('Differential evolution and simulated annealing find global minima of multimodal functions (seeded)', () => {
  const lo = new Array(5).fill(-5.12), hi = new Array(5).fill(5.12);
  const de = differentialEvolution(rastrigin, lo, hi, { maxGenerations: 2000 });
  assert.ok(de.converged); assert.ok(de.f < 1e-9); assert.ok(de.x.every(v => Math.abs(v) < 1e-5));
  const again = differentialEvolution(rastrigin, lo, hi, { maxGenerations: 2000 });
  assert.deepEqual(again.x, de.x);                      // deterministic for a given seed
  const h = differentialEvolution(himmelblau, [-5, -5], [5, 5], { strategy: 'best1bin', polish: true });
  assert.ok(nearestHimmelblau(h.x) < 1e-7);
  // minimum on the boundary of the box
  vclose(differentialEvolution(x => x[0] + x[1] + x[2], [1, 1, 1], [2, 2, 2]).x, [1, 1, 1], 1e-8);
  for (const seed of [1, 2, 3]) {
    const sa = simulatedAnnealing(rastrigin, [3.3, -2.1], { lower: -5.12, upper: 5.12, seed, polish: true });
    assert.ok(sa.f < 1e-12, `seed ${seed}: ${sa.f}`); vclose(sa.x, [0, 0], 1e-6);
  }
});

test('finite differences: error orders 1, 2, 4 measured; Hessian matches Rosenbrock exactly', () => {
  // measured order: halving h divides the error by 2^order
  const f = x => Math.exp(x[0]), exact = Math.E;
  for (const [order, h, expect] of [[1, 1e-3, 2], [2, 1e-2, 4], [4, 1e-1, 16]]) {
    const e1 = Math.abs(finiteGradient(f, [1], { order, h })[0] - exact), e2 = Math.abs(finiteGradient(f, [1], { order, h: h / 2 })[0] - exact);
    close(e1 / e2, expect, 0.1 * expect);
  }
  // default steps reach the advertised accuracy
  const x = [0.3, 0.7], g = rosenGrad(x), scale = Math.hypot(...g);
  for (const [order, tol] of [[1, 1e-6], [2, 1e-9], [4, 1e-11]]) {
    const fd = finiteGradient(rosen, x, { order });
    assert.ok(Math.hypot(fd[0] - g[0], fd[1] - g[1]) / scale < tol, `order ${order}`);
  }
  const H = [[1200 * 0.09 - 400 * 0.7 + 2, -400 * 0.3], [-400 * 0.3, 200]];
  const Hg = finiteHessian(null, x, { gradient: rosenGrad }), Hf = finiteHessian(rosen, x);
  for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) { close(Hg[i][j], H[i][j], 1e-6); close(Hf[i][j], H[i][j], 1e-3); }
  assert.throws(() => finiteGradient(x0 => Math.sqrt(x0[0]), [0]), /not finite/);
});

test('checkGradient accepts correct gradients and flags a wrong component', () => {
  const rand = rng(5);
  for (let k = 0; k < 5; k++) {
    const x = [2 * rand() - 1, 2 * rand() - 1];
    assert.ok(checkGradient(rosen, rosenGrad, x, { rtol: 1e-9 }).ok);
  }
  const wrong = x => { const g = rosenGrad(x); g[1] *= 1.001; return g; };
  const c = checkGradient(rosen, wrong, [0.3, 0.7]);
  assert.equal(c.ok, false); assert.equal(c.worst, 1);
});

test('L-BFGS with autodiff gradients relaxes an LJ trimer to the exact minimum (equilateral, r = 2^(1/6), E = −3ε)', () => {
  // reaching gtol 1e-8 needs the approximate-Wolfe branch: f is at its rounding floor long before ∇f is
  const { add, sub, mul, pow, square } = AD;
  const energy = x => {
    let E = 0;
    for (let i = 0; i < x.length; i += 3) {
      for (let j = i + 3; j < x.length; j += 3) {
        const r2 = add(add(square(sub(x[i], x[j])), square(sub(x[i + 1], x[j + 1]))), square(sub(x[i + 2], x[j + 2])));
        const s6 = pow(r2, -3);
        E = add(E, mul(4, sub(square(s6), s6)));
      }
    }
    return E;
  };
  for (const minimise of [lbfgs, bfgs]) {
    const r = minimise(x => AD.gradient(energy, x).value, [0, 0, 0, 1.2, 0, 0, 0.5, 1, 0.1], { gradient: x => AD.gradient(energy, x).gradient, gtol: 1e-8 });
    assert.ok(r.converged, r.reason); close(r.f, -3, 1e-14);
    const d = (i, j) => Math.hypot(r.x[3 * i] - r.x[3 * j], r.x[3 * i + 1] - r.x[3 * j + 1], r.x[3 * i + 2] - r.x[3 * j + 2]);
    for (const [i, j] of [[0, 1], [0, 2], [1, 2]]) close(d(i, j), 2 ** (1 / 6), 1e-9);
  }
});
