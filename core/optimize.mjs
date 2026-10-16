// Minimisation: one-dimensional, derivative-free, quasi-Newton, nonlinear least squares with
// bounds and robust losses, global (seeded) search, and finite-difference derivatives.
// Plain ES module (browser + Node), no dependencies beyond core/. Deterministic for fixed seeds.
//
//   goldenSection / brentMinimize   1D minimum inside a bracket [a, b]
//   nelderMead                      derivative-free simplex (adaptive Gao–Han coefficients)
//   wolfeLineSearch                 strong-Wolfe line search (Nocedal & Wright Alg. 3.5/3.6)
//   bfgs / lbfgs                    quasi-Newton; lbfgs also takes box bounds (projected)
//   levenbergMarquardt              residual-vector least squares: bounds, robust loss, covariance
//   simulatedAnnealing              seeded Metropolis search with adaptive steps
//   differentialEvolution           seeded DE/rand/1/bin or DE/best/1/bin inside a box
//   finiteGradient / finiteHessian  central (or 4th-order) differences with scaled steps
//   checkGradient                   compare an analytic gradient against finite differences
//
// Every minimiser returns { x, f, iterations, evaluations, converged, reason, ... }. A
// minimiser certifies (at best) a local stationary point; it never proves a global minimum.
// Data: curveFit (core/data.mjs) remains the y = model(x, p) front end; levenbergMarquardt is
// the general residual-vector solver with bounds and robust losses.

import { rng } from './numerics.mjs';
import { solveLinear } from './linalg.mjs';

const EPS = Number.EPSILON;
const GOLDEN = (3 - Math.sqrt(5)) / 2;          // 0.381966…

const fn = (f, who, name = 'objective') => {
  if (typeof f !== 'function') throw new Error(`${who}: ${name} must be a function`);
};
const vec = (x, who, name = 'x0') => {
  if (!x || !(x.length > 0)) throw new Error(`${who}: ${name} must be a non-empty vector`);
  const v = Float64Array.from(x);
  if (!v.every(Number.isFinite)) throw new Error(`${who}: ${name} must be finite`);
  return v;
};
const posInt = (v, who, name) => {
  if (!(Number.isSafeInteger(v) && v > 0)) throw new Error(`${who}: ${name} must be a positive integer`);
};
const nonNeg = (v, who, name) => {
  if (!(Number.isFinite(v) && v >= 0)) throw new Error(`${who}: ${name} must be finite and non-negative`);
};
const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
const normInf = a => { let m = 0; for (const v of a) m = Math.max(m, Math.abs(v)); return m; };
const norm2 = a => Math.sqrt(dot(a, a));
const finiteValue = (v, who) => {
  if (typeof v !== 'number' || !Number.isFinite(v)) throw new Error(`${who}: objective is not finite at the start point`);
  return v;
};
const boxOf = (lower, upper, n, who) => {
  if (lower == null && upper == null) return null;
  const expand = (b, fill) => (b == null ? new Float64Array(n).fill(fill) : typeof b === 'number' ? new Float64Array(n).fill(b) : Float64Array.from(b));
  const lo = expand(lower, -Infinity), hi = expand(upper, Infinity);
  if (lo.length !== n || hi.length !== n) throw new Error(`${who}: bounds must have one entry per variable`);
  for (let i = 0; i < n; i++) if (Number.isNaN(lo[i]) || Number.isNaN(hi[i]) || !(lo[i] <= hi[i])) throw new Error(`${who}: need lower ≤ upper for every variable`);
  return { lo, hi };
};
const project = (x, box) => {
  if (box) for (let i = 0; i < x.length; i++) x[i] = Math.min(box.hi[i], Math.max(box.lo[i], x[i]));
  return x;
};

// ---------------------------------------------------------------- one dimension

const bracketOf = (f, a, b, who) => {
  fn(f, who);
  if (!(Number.isFinite(a) && Number.isFinite(b) && a !== b)) throw new Error(`${who}: need a finite interval with a ≠ b`);
  return a < b ? [a, b] : [b, a];
};

/**
 * Golden-section search for a minimum of f on [a, b]. Conditions: f unimodal on the interval
 * (otherwise a local minimum or an endpoint is returned). Robust and derivative-free; the
 * bracket shrinks by 0.618 per evaluation (linear convergence). Stops when the bracket is
 * ≤ 2(xtol·|x| + atol); xtol below √ε ≈ 1.5e-8 is futile because f is flat to rounding within
 * √ε·|x| of a smooth minimum. Default atol = 1e-10·(b − a).
 * Returns { x, f, iterations, evaluations, converged, reason, bracket }.
 */
export function goldenSection(f, a, b, { xtol = Math.sqrt(EPS), atol, maxIterations = 500 } = {}) {
  let [lo, hi] = bracketOf(f, a, b, 'goldenSection');
  const abs = atol ?? 1e-10 * (hi - lo);
  nonNeg(xtol, 'goldenSection', 'xtol'); nonNeg(abs, 'goldenSection', 'atol'); posInt(maxIterations, 'goldenSection', 'maxIterations');
  const r = 1 - GOLDEN;
  let c = hi - r * (hi - lo), d = lo + r * (hi - lo), fc = f(c), fd = f(d), evaluations = 2;
  for (let it = 0; it <= maxIterations; it++) {
    const x = fc <= fd ? c : d, fx = Math.min(fc, fd);
    if (hi - lo <= 2 * (xtol * Math.abs(x) + abs) || it === maxIterations) {
      const converged = hi - lo <= 2 * (xtol * Math.abs(x) + abs);
      return { x, f: fx, iterations: it, evaluations, converged, reason: converged ? 'bracket tolerance' : 'iteration budget exhausted', bracket: [lo, hi] };
    }
    if (fc <= fd || Number.isNaN(fd)) { hi = d; d = c; fd = fc; c = hi - r * (hi - lo); fc = f(c); }
    else { lo = c; c = d; fc = fd; d = lo + r * (hi - lo); fd = f(d); }
    evaluations++;
  }
  throw new Error('goldenSection: unreachable');
}

/**
 * Brent's minimiser (Brent 1973, "localmin"): parabolic interpolation safeguarded by golden
 * sections on [a, b]. Same conditions, tolerances and result as goldenSection; superlinear on
 * smooth minima (typically 2–4× fewer evaluations), never slower than golden section by more
 * than a constant factor. Accuracy limit: √ε relative in x (see goldenSection).
 */
export function brentMinimize(f, a, b, { xtol = Math.sqrt(EPS), atol, maxIterations = 500 } = {}) {
  let [lo, hi] = bracketOf(f, a, b, 'brentMinimize');
  const abs = atol ?? 1e-10 * (hi - lo);
  nonNeg(xtol, 'brentMinimize', 'xtol'); nonNeg(abs, 'brentMinimize', 'atol'); posInt(maxIterations, 'brentMinimize', 'maxIterations');
  let x = lo + GOLDEN * (hi - lo), w = x, v = x, fx = f(x), fw = fx, fv = fx, d = 0, e = 0, evaluations = 1;
  for (let it = 0; it < maxIterations; it++) {
    const m = (lo + hi) / 2, tol1 = xtol * Math.abs(x) + abs, tol2 = 2 * tol1;
    if (Math.abs(x - m) <= tol2 - (hi - lo) / 2) {
      return { x, f: fx, iterations: it, evaluations, converged: true, reason: 'bracket tolerance', bracket: [lo, hi] };
    }
    let golden = true;
    if (Math.abs(e) > tol1) {
      const r = (x - w) * (fx - fv);
      let q = (x - v) * (fx - fw), p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) p = -p; else q = -q;
      const eOld = e; e = d;
      if (Math.abs(p) < Math.abs(0.5 * q * eOld) && p > q * (lo - x) && p < q * (hi - x)) {
        d = p / q; golden = false;
        const u = x + d;
        if (u - lo < tol2 || hi - u < tol2) d = x < m ? tol1 : -tol1;
      }
    }
    if (golden) { e = (x < m ? hi : lo) - x; d = GOLDEN * e; }
    const u = x + (Math.abs(d) >= tol1 ? d : d > 0 ? tol1 : -tol1), fu = f(u);
    evaluations++;
    if (fu <= fx) {
      if (u < x) hi = x; else lo = x;
      v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
    } else {
      if (u < x) lo = u; else hi = u;
      if (fu <= fw || w === x) { v = w; fv = fw; w = u; fw = fu; }
      else if (fu <= fv || v === x || v === w) { v = u; fv = fu; }
    }
  }
  return { x, f: fx, iterations: maxIterations, evaluations, converged: false, reason: 'iteration budget exhausted', bracket: [lo, hi] };
}

// ---------------------------------------------------------------- Nelder–Mead

/**
 * Nelder–Mead simplex minimisation with the adaptive coefficients of Gao & Han (2012), which
 * keep it effective beyond ~5 dimensions. Derivative-free; for noisy, non-smooth or cheap
 * low-dimensional objectives (n ≲ 10). No convergence theory in general (it can stall on a
 * non-stationary point): restart from the result to confirm. Initial simplex x0 + step_i·e_i
 * (default step 0.1·max(|x0_i|, 1)). Optional box `lower`/`upper` (trial points are clipped).
 * Converged when the simplex ∞-diameter ≤ xtol AND the f spread ≤ ftol (both absolute).
 * Returns { x, f, iterations, evaluations, converged, reason, simplex }.
 */
export function nelderMead(f, x0, { step, xtol = 1e-8, ftol = 1e-12, maxIterations, maxEvaluations, adaptive = true, lower, upper } = {}) {
  fn(f, 'nelderMead');
  const start = vec(x0, 'nelderMead'), n = start.length, box = boxOf(lower, upper, n, 'nelderMead');
  const maxIt = maxIterations ?? 2000 * n, maxEv = maxEvaluations ?? 4000 * n;
  nonNeg(xtol, 'nelderMead', 'xtol'); nonNeg(ftol, 'nelderMead', 'ftol'); posInt(maxIt, 'nelderMead', 'maxIterations'); posInt(maxEv, 'nelderMead', 'maxEvaluations');
  const [alpha, beta, gamma, delta] = adaptive ? [1, 1 + 2 / n, 0.75 - 1 / (2 * n), 1 - 1 / n] : [1, 2, 0.5, 0.5];
  let evaluations = 0;
  const F = x => { evaluations++; const v = f(x); return Number.isNaN(v) ? Infinity : v; };
  project(start, box);
  const simplex = [start];
  for (let i = 0; i < n; i++) {
    const p = start.slice(), h = step == null ? 0.1 * Math.max(Math.abs(start[i]), 1) : typeof step === 'number' ? step : step[i];
    if (!(Number.isFinite(h) && h !== 0)) throw new Error('nelderMead: step must be finite and non-zero');
    p[i] += h;
    if (box && p[i] > box.hi[i]) p[i] = start[i] - h;
    simplex.push(project(p, box));
  }
  let values = simplex.map(F);
  finiteValue(values[0], 'nelderMead');
  const order = () => {
    const idx = values.map((_, i) => i).sort((i, j) => values[i] - values[j]);
    const s = idx.map(i => simplex[i]), v = idx.map(i => values[i]);
    for (let i = 0; i <= n; i++) { simplex[i] = s[i]; values[i] = v[i]; }
  };
  const along = (c, x, t) => project(Float64Array.from(c, (ci, k) => ci + t * (x[k] - ci)), box);
  let it = 0, reason = 'iteration budget exhausted', converged = false;
  for (; it < maxIt; it++) {
    order();
    let size = 0;
    for (let i = 1; i <= n; i++) for (let k = 0; k < n; k++) size = Math.max(size, Math.abs(simplex[i][k] - simplex[0][k]));
    if (size <= xtol && values[n] - values[0] <= ftol) { converged = true; reason = 'simplex tolerance'; break; }
    if (evaluations >= maxEv) { reason = 'evaluation budget exhausted'; break; }
    const c = new Float64Array(n);
    for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) c[k] += simplex[i][k] / n;
    const xr = along(c, simplex[n], -alpha), fr = F(xr);
    if (fr < values[0]) {
      const xe = along(c, xr, beta), fe = F(xe);
      if (fe < fr) { simplex[n] = xe; values[n] = fe; } else { simplex[n] = xr; values[n] = fr; }
      continue;
    }
    if (fr < values[n - 1]) { simplex[n] = xr; values[n] = fr; continue; }
    let shrink = true;
    if (fr < values[n]) {
      const xc = along(c, xr, gamma), fcv = F(xc);
      if (fcv <= fr) { simplex[n] = xc; values[n] = fcv; shrink = false; }
    } else {
      const xc = along(c, simplex[n], gamma), fcv = F(xc);
      if (fcv < values[n]) { simplex[n] = xc; values[n] = fcv; shrink = false; }
    }
    if (shrink) {
      for (let i = 1; i <= n; i++) { simplex[i] = along(simplex[0], simplex[i], delta); values[i] = F(simplex[i]); }
    }
  }
  order();
  values = values.slice();
  return { x: simplex[0], f: values[0], iterations: it, evaluations, converged, reason, simplex: simplex.map(p => p.slice()) };
}

// ---------------------------------------------------------------- finite differences

const stepFor = (x, i, h, scale, c) => {
  if (h != null) return typeof h === 'number' ? h : h[i];
  const s = scale == null ? 1 : typeof scale === 'number' ? scale : scale[i];
  return c * Math.max(Math.abs(x[i]), s);
};

/**
 * Finite-difference gradient. order 1 (forward, h = √ε·max(|x_i|, scale_i), error ~1e-8
 * relative), 2 (central, h = ε^{1/3}·…, error ~1e-10; default) or 4 (five-point, h = ε^{1/5}·…,
 * error ~1e-12). `scale`: typical magnitude of each variable (default 1) so steps are sensible
 * for x ≈ 0; `h` overrides the step (number or per-variable array). Steps are rounded so x ± h
 * is exactly representable. Throws if f is not finite at a stencil point (e.g. outside its
 * domain: use a smaller h or order 1 away from the boundary). Returns a Float64Array.
 */
export function finiteGradient(f, x, { order = 2, h, scale, f0 } = {}) {
  fn(f, 'finiteGradient');
  const x0 = vec(x, 'finiteGradient', 'x'), n = x0.length, g = new Float64Array(n);
  if (![1, 2, 4].includes(order)) throw new Error('finiteGradient: order must be 1, 2 or 4');
  const c = order === 1 ? Math.sqrt(EPS) : order === 2 ? Math.cbrt(EPS) : EPS ** 0.2;
  const at = (i, t) => { const p = x0.slice(); p[i] += t; const v = f(p); if (!Number.isFinite(v)) throw new Error('finiteGradient: objective not finite at a stencil point'); return v; };
  const base = order === 1 ? (f0 ?? f(x0)) : 0;
  if (order === 1 && !Number.isFinite(base)) throw new Error('finiteGradient: objective not finite at x');
  for (let i = 0; i < n; i++) {
    let hi = stepFor(x0, i, h, scale, c);
    hi = (x0[i] + hi) - x0[i];
    if (!(Number.isFinite(hi) && hi !== 0)) throw new Error('finiteGradient: step cannot be resolved at this x');
    if (order === 1) g[i] = (at(i, hi) - base) / hi;
    else if (order === 2) g[i] = (at(i, hi) - at(i, -hi)) / (2 * hi);
    else g[i] = (8 * (at(i, hi) - at(i, -hi)) - (at(i, 2 * hi) - at(i, -2 * hi))) / (12 * hi);
  }
  return g;
}

/**
 * Finite-difference Hessian (symmetric, rows as Float64Array). With `gradient` (analytic g(x)):
 * central differences of g, h = ε^{1/3}·max(|x_j|, scale_j), error ~1e-10 relative. Without:
 * second differences of f, h = ε^{1/4}·…, error ~1e-7 relative (cancellation-limited).
 * Cost 2n gradient or 2n² + 1 function evaluations. Same `scale` and failure rules as
 * finiteGradient.
 */
export function finiteHessian(f, x, { gradient, scale, h } = {}) {
  if (gradient === undefined) fn(f, 'finiteHessian');
  const x0 = vec(x, 'finiteHessian', 'x'), n = x0.length, H = Array.from({ length: n }, () => new Float64Array(n));
  const steps = Float64Array.from(x0, (_, i) => {
    const s = stepFor(x0, i, h, scale, gradient ? Math.cbrt(EPS) : EPS ** 0.25);
    const r = (x0[i] + s) - x0[i];
    if (!(Number.isFinite(r) && r !== 0)) throw new Error('finiteHessian: step cannot be resolved at this x');
    return r;
  });
  const shifted = (pairs) => { const p = x0.slice(); for (const [i, t] of pairs) p[i] += t; return p; };
  if (gradient) {
    fn(gradient, 'finiteHessian', 'gradient');
    for (let j = 0; j < n; j++) {
      const gp = gradient(shifted([[j, steps[j]]])), gm = gradient(shifted([[j, -steps[j]]]));
      for (let i = 0; i < n; i++) H[i][j] = (gp[i] - gm[i]) / (2 * steps[j]);
    }
    for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) { const m = (H[i][j] + H[j][i]) / 2; H[i][j] = m; H[j][i] = m; }
  } else {
    const F = p => { const v = f(p); if (!Number.isFinite(v)) throw new Error('finiteHessian: objective not finite at a stencil point'); return v; };
    const f0 = F(x0);
    for (let i = 0; i < n; i++) {
      const hi = steps[i];
      H[i][i] = (F(shifted([[i, hi]])) - 2 * f0 + F(shifted([[i, -hi]]))) / (hi * hi);
      for (let j = 0; j < i; j++) {
        const hj = steps[j];
        const v = (F(shifted([[i, hi], [j, hj]])) - F(shifted([[i, hi], [j, -hj]])) - F(shifted([[i, -hi], [j, hj]])) + F(shifted([[i, -hi], [j, -hj]]))) / (4 * hi * hj);
        H[i][j] = v; H[j][i] = v;
      }
    }
  }
  return H;
}

/**
 * Check an analytic gradient against a 4th-order finite-difference gradient at x. The error is
 * normwise: |g_analytic,i − g_numeric,i| / max(‖g_numeric‖∞, ‖g_analytic‖∞, tiny), so a component
 * that is wrong by more than rtol of the gradient's size fails. Returns
 * { ok, maxError, worst, analytic, numeric, errors }. Use at several random points; a pass at
 * one point (especially a symmetric one) proves little.
 */
export function checkGradient(f, gradient, x, { rtol = 1e-6, scale, h } = {}) {
  fn(f, 'checkGradient'); fn(gradient, 'checkGradient', 'gradient');
  const x0 = vec(x, 'checkGradient', 'x');
  const analytic = Float64Array.from(gradient(x0.slice()));
  if (analytic.length !== x0.length || !analytic.every(Number.isFinite)) throw new Error('checkGradient: gradient must return one finite value per variable');
  const numeric = finiteGradient(f, x0, { order: 4, scale, h });
  const size = Math.max(normInf(numeric), normInf(analytic), Number.MIN_VALUE);
  const errors = Float64Array.from(analytic, (a, i) => Math.abs(a - numeric[i]) / size);
  let worst = 0;
  for (let i = 1; i < errors.length; i++) if (errors[i] > errors[worst]) worst = i;
  return { ok: errors[worst] <= rtol, maxError: errors[worst], worst, analytic, numeric, errors };
}

// ---------------------------------------------------------------- line search

/** Minimiser of the cubic matching values and slopes at a and b; NaN if it does not exist. */
const cubicMin = (a, fa, da, b, fb, db) => {
  const d1 = da + db - 3 * (fa - fb) / (a - b), disc = d1 * d1 - da * db;
  if (!(disc >= 0)) return NaN;
  const d2 = Math.sign(b - a) * Math.sqrt(disc);
  return b - (b - a) * (db + d2 - d1) / (db - da + 2 * d2);
};

const wolfeCore = (phi, f0, d0, { alpha0, c1, c2, alphaMax, maxEvaluations }) => {
  let evaluations = 0;
  const ev = a => { evaluations++; const t = phi(a); if (!Number.isFinite(t.f) || !Number.isFinite(t.d)) t.f = Infinity; return t; };
  // Armijo, or near the rounding floor of f the approximate Wolfe condition of Hager & Zhang
  // (2005): f within 1e-12|f0| of f0 and φ′(α) ≤ (1 − 2c1)|φ′(0)|. Without it a search stalls
  // once the decrease of f is below its rounding error although ∇f can still be reduced.
  const fNoise = 1e-12 * Math.abs(f0);
  const sufficient = (a, t) => t.f <= f0 + c1 * a * d0 || (t.f <= f0 + fNoise && t.d <= (2 * c1 - 1) * d0);
  const zoom = (lo0, tlo0, hi0, thi0) => {
    let lo = lo0, tlo = tlo0, hi = hi0, thi = thi0;
    for (;;) {
      if (evaluations >= maxEvaluations) break;
      const span = hi - lo;
      let a = Number.isFinite(thi.f) ? cubicMin(lo, tlo.f, tlo.d, hi, thi.f, thi.d) : NaN;
      const left = Math.min(lo, hi), right = Math.max(lo, hi), margin = 0.1 * Math.abs(span);
      if (!(a > left + margin && a < right - margin)) a = lo + span / 2;
      if (a === lo || a === hi) break;
      const t = ev(a);
      if (!sufficient(a, t) || t.f > tlo.f + fNoise) { hi = a; thi = t; }
      else {
        if (Math.abs(t.d) <= -c2 * d0) return { alpha: a, t, evaluations, ok: true };
        if (t.d * (hi - lo) >= 0) { hi = lo; thi = tlo; }
        lo = a; tlo = t;
      }
    }
    return { alpha: lo, t: tlo, evaluations, ok: false };
  };
  let aPrev = 0, tPrev = { f: f0, d: d0 }, a = alpha0;
  for (let i = 0; ; i++) {
    if (evaluations >= maxEvaluations) break;
    const t = ev(a);
    if (!sufficient(a, t) || (i > 0 && t.f > tPrev.f + fNoise)) return zoom(aPrev, tPrev, a, t);
    if (Math.abs(t.d) <= -c2 * d0) return { alpha: a, t, evaluations, ok: true };
    if (t.d >= 0) return zoom(a, t, aPrev, tPrev);
    if (a >= alphaMax) return { alpha: a, t, evaluations, ok: false };
    aPrev = a; tPrev = t; a = Math.min(2 * a, alphaMax);
  }
  return { alpha: aPrev, t: tPrev, evaluations, ok: false };
};

/**
 * Strong-Wolfe line search along a descent direction d from x (Nocedal & Wright 2006,
 * Alg. 3.5 with cubic-interpolation zoom 3.6): finds α with
 *   f(x + αd) ≤ f(x) + c1·α·g·d   and   |g(x + αd)·d| ≤ c2·|g(x)·d|   (0 < c1 < c2 < 1).
 * gradient(x) → ∇f. Throws if d is not a descent direction. Non-finite trial values count as
 * "too long". Returns { alpha, x, f, g, evaluations, ok }; ok false means no Wolfe point was found
 * within maxEvaluations (the best sufficient-decrease point so far is returned, alpha may be 0).
 */
export function wolfeLineSearch(f, gradient, x, d, { f0, g0, alpha0 = 1, c1 = 1e-4, c2 = 0.9, alphaMax = 1e10, maxEvaluations = 40 } = {}) {
  fn(f, 'wolfeLineSearch'); fn(gradient, 'wolfeLineSearch', 'gradient');
  const x0 = vec(x, 'wolfeLineSearch', 'x'), dir = vec(d, 'wolfeLineSearch', 'd');
  if (dir.length !== x0.length) throw new Error('wolfeLineSearch: x and d must have equal length');
  if (!(c1 > 0 && c1 < c2 && c2 < 1)) throw new Error('wolfeLineSearch: need 0 < c1 < c2 < 1');
  const fx = f0 ?? f(x0), gx = g0 ? Float64Array.from(g0) : Float64Array.from(gradient(x0));
  const slope = dot(gx, dir);
  if (!(slope < 0)) throw new Error('wolfeLineSearch: d is not a descent direction (g·d ≥ 0)');
  const cache = new Map();
  const phi = a => {
    const p = Float64Array.from(x0, (v, i) => v + a * dir[i]), fv = f(p), gv = Float64Array.from(gradient(p));
    cache.set(a, { p, fv, gv });
    return { f: fv, d: dot(gv, dir) };
  };
  const r = wolfeCore(phi, fx, slope, { alpha0, c1, c2, alphaMax, maxEvaluations });
  const hit = cache.get(r.alpha);
  return { alpha: r.alpha, x: hit ? hit.p : x0, f: hit ? hit.fv : fx, g: hit ? hit.gv : gx, evaluations: r.evaluations, ok: r.ok };
}

// ---------------------------------------------------------------- quasi-Newton

const objectiveWithGradient = (f, gradient, who, gradOptions) => {
  fn(f, who);
  if (gradient !== undefined) fn(gradient, who, 'gradient');
  const counts = { f: 0, g: 0 };
  const fg = x => {
    counts.f++;
    const v = f(x);
    if (!Number.isFinite(v)) return { f: Infinity, g: null };
    counts.g++;
    let g;
    if (gradient) g = Float64Array.from(gradient(x));
    else {
      g = finiteGradient(y => { counts.f++; return f(y); }, x, gradOptions);
    }
    if (g.length !== x.length) throw new Error(`${who}: gradient must have one entry per variable`);
    return { f: v, g };
  };
  return { fg, counts };
};

/**
 * BFGS quasi-Newton minimisation with a strong-Wolfe line search (Nocedal & Wright 2006,
 * Alg. 6.1). Conditions: f smooth (C²) near the minimum; gradient(x) → ∇f optional (central
 * differences otherwise, which limit attainable gtol to ~1e-8·|f|/|x|). Superlinear local
 * convergence. Dense n×n inverse Hessian: use lbfgs beyond a few hundred variables.
 * Converged when ‖∇f‖∞ ≤ gtol (absolute; rescale f or loosen gtol for large |f|).
 * onIteration({ iteration, x, f, gradientNorm }) is called after every accepted step.
 * Returns { x, f, gradient, iterations, evaluations, gradientEvaluations, converged, reason,
 * inverseHessian } (inverseHessian approximates H⁻¹ at x; it is not an exact covariance).
 */
export function bfgs(f, x0, { gradient, gtol = 1e-8, maxIterations = 1000, c1 = 1e-4, c2 = 0.9, onIteration, gradientOptions } = {}) {
  const { fg, counts } = objectiveWithGradient(f, gradient, 'bfgs', gradientOptions);
  let x = vec(x0, 'bfgs'); const n = x.length;
  nonNeg(gtol, 'bfgs', 'gtol'); posInt(maxIterations, 'bfgs', 'maxIterations');
  let cur = fg(x);
  finiteValue(cur.f, 'bfgs');
  let H = null;
  const identity = () => Array.from({ length: n }, (_, i) => { const r = new Float64Array(n); r[i] = 1; return r; });
  H = identity();
  let it = 0, reason = 'iteration budget exhausted', fresh = true;
  for (; it < maxIterations; it++) {
    if (normInf(cur.g) <= gtol) { reason = 'gradient tolerance'; break; }
    let d = Float64Array.from(H, row => -dot(row, cur.g));
    if (!(dot(d, cur.g) < 0)) { H = identity(); fresh = true; d = Float64Array.from(cur.g, v => -v); }
    const alpha0 = fresh && it === 0 ? Math.min(1, 1 / normInf(cur.g)) : 1;
    const cache = new Map();
    const phi = a => { const p = Float64Array.from(x, (v, i) => v + a * d[i]), t = fg(p); cache.set(a, { p, t }); return { f: t.f, d: t.g ? dot(t.g, d) : NaN }; };
    const ls = wolfeCore(phi, cur.f, dot(cur.g, d), { alpha0, c1, c2, alphaMax: 1e10, maxEvaluations: 40 });
    const hit = cache.get(ls.alpha);
    if (!hit || !(hit.t.f < cur.f || (ls.ok && hit.t.f <= cur.f + 1e-12 * Math.abs(cur.f)))) {
      if (!fresh) { H = identity(); fresh = true; it--; continue; }
      reason = 'line search failed'; break;
    }
    const s = Float64Array.from(hit.p, (v, i) => v - x[i]), y = Float64Array.from(hit.t.g, (v, i) => v - cur.g[i]);
    const sy = dot(s, y);
    if (sy > EPS * norm2(s) * norm2(y)) {
      if (fresh) { const scale = sy / dot(y, y); for (let i = 0; i < n; i++) H[i][i] = scale; }
      const Hy = Float64Array.from(H, row => dot(row, y)), yHy = dot(y, Hy), a = (sy + yHy) / (sy * sy);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) H[i][j] += a * s[i] * s[j] - (Hy[i] * s[j] + s[i] * Hy[j]) / sy;
      fresh = false;
    }
    const stepSize = normInf(s);
    x = hit.p; cur = hit.t;
    onIteration?.({ iteration: it + 1, x: x.slice(), f: cur.f, gradientNorm: normInf(cur.g) });
    if (stepSize <= EPS * Math.max(1, normInf(x))) { it++; reason = normInf(cur.g) <= gtol ? 'gradient tolerance' : 'step too small'; break; }
  }
  return { x, f: cur.f, gradient: cur.g, iterations: it, evaluations: counts.f, gradientEvaluations: counts.g, converged: normInf(cur.g) <= gtol, reason, inverseHessian: H };
}

/**
 * Limited-memory BFGS (Nocedal 1980; two-loop recursion, `memory` pairs, default 10) for large
 * smooth problems: O(n·memory) work and storage per step. Unconstrained: strong-Wolfe line
 * search. With box bounds `lower`/`upper` (numbers or arrays, ±Infinity allowed): projected
 * L-BFGS — variables at a bound with the gradient pushing outwards are frozen, the quasi-Newton
 * step acts on the free variables, and a projected backtracking (Armijo) search keeps iterates
 * feasible. Converged when the projected gradient satisfies ‖P∇f‖∞ ≤ gtol. Same gradient,
 * onIteration and result conventions as bfgs (plus `active`: which variables sit on a bound).
 */
export function lbfgs(f, x0, { gradient, memory = 10, gtol = 1e-8, maxIterations = 5000, lower, upper, c1 = 1e-4, c2 = 0.9, onIteration, gradientOptions } = {}) {
  const { fg, counts } = objectiveWithGradient(f, gradient, 'lbfgs', gradientOptions);
  let x = vec(x0, 'lbfgs'); const n = x.length, box = boxOf(lower, upper, n, 'lbfgs');
  posInt(memory, 'lbfgs', 'memory'); nonNeg(gtol, 'lbfgs', 'gtol'); posInt(maxIterations, 'lbfgs', 'maxIterations');
  project(x, box);
  let cur = fg(x);
  finiteValue(cur.f, 'lbfgs');
  const pairs = [];
  const activeOf = (p, g) => Array.from(p, (v, i) => box != null && ((v <= box.lo[i] && g[i] > 0) || (v >= box.hi[i] && g[i] < 0)));
  const projGrad = (p, g) => { const act = activeOf(p, g); return Float64Array.from(g, (v, i) => (act[i] ? 0 : v)); };
  let it = 0, reason = 'iteration budget exhausted';
  for (; it < maxIterations; it++) {
    const pg = projGrad(x, cur.g);
    if (normInf(pg) <= gtol) { reason = 'gradient tolerance'; break; }
    const act = activeOf(x, cur.g);
    // two-loop recursion on the free subspace
    const q = pg.slice(), alphas = [];
    for (let k = pairs.length - 1; k >= 0; k--) { const { s, y, rho } = pairs[k], a = rho * dot(s, q); alphas[k] = a; for (let i = 0; i < n; i++) q[i] -= a * y[i]; }
    const last = pairs[pairs.length - 1], gamma = last ? dot(last.s, last.y) / dot(last.y, last.y) : 1;
    for (let i = 0; i < n; i++) q[i] *= gamma;
    for (let k = 0; k < pairs.length; k++) { const { s, y, rho } = pairs[k], b = rho * dot(y, q); for (let i = 0; i < n; i++) q[i] += (alphas[k] - b) * s[i]; }
    let d = Float64Array.from(q, (v, i) => (act[i] ? 0 : -v));
    if (!(dot(d, cur.g) < 0)) { pairs.length = 0; d = Float64Array.from(pg, v => -v); }
    const alpha0 = pairs.length === 0 ? Math.min(1, 1 / normInf(d)) : 1;
    let next = null;
    if (!box) {
      const cache = new Map();
      const phi = a => { const p = Float64Array.from(x, (v, i) => v + a * d[i]), t = fg(p); cache.set(a, { p, t }); return { f: t.f, d: t.g ? dot(t.g, d) : NaN }; };
      const ls = wolfeCore(phi, cur.f, dot(cur.g, d), { alpha0, c1, c2, alphaMax: 1e10, maxEvaluations: 40 });
      const hit = cache.get(ls.alpha);
      if (hit && (hit.t.f < cur.f || (ls.ok && hit.t.f <= cur.f + 1e-12 * Math.abs(cur.f)))) next = hit;
    } else {
      let a = alpha0;
      for (let k = 0; k < 50; k++) {
        const p = project(Float64Array.from(x, (v, i) => v + a * d[i]), box), t = fg(p);
        const decrease = dot(cur.g, Float64Array.from(p, (v, i) => v - x[i]));
        if (t.f <= cur.f + c1 * decrease && t.f < cur.f) { next = { p, t }; break; }
        a *= 0.5;
      }
    }
    if (!next) {
      if (pairs.length) { pairs.length = 0; it--; continue; }
      reason = 'line search failed'; break;
    }
    const s = Float64Array.from(next.p, (v, i) => v - x[i]), y = Float64Array.from(next.t.g, (v, i) => v - cur.g[i]), sy = dot(s, y);
    if (sy > EPS * dot(y, y)) { pairs.push({ s, y, rho: 1 / sy }); if (pairs.length > memory) pairs.shift(); }
    const stepSize = normInf(s);
    x = next.p; cur = next.t;
    onIteration?.({ iteration: it + 1, x: x.slice(), f: cur.f, gradientNorm: normInf(projGrad(x, cur.g)) });
    if (stepSize <= EPS * Math.max(1, normInf(x))) { it++; reason = 'step too small'; break; }
  }
  const pgNorm = normInf(projGrad(x, cur.g));
  if (pgNorm <= gtol) reason = 'gradient tolerance';
  return { x, f: cur.f, gradient: cur.g, iterations: it, evaluations: counts.f, gradientEvaluations: counts.g, converged: pgNorm <= gtol, reason, active: activeOf(x, cur.g) };
}

// ---------------------------------------------------------------- nonlinear least squares

const LOSSES = {
  linear: { rho: z => z, w: () => 1 },
  soft_l1: { rho: z => 2 * (z / (Math.sqrt(1 + z) + 1)), w: z => 1 / Math.sqrt(1 + z) },
  huber: { rho: z => (z <= 1 ? z : 2 * Math.sqrt(z) - 1), w: z => (z <= 1 ? 1 : 1 / Math.sqrt(z)) },
  cauchy: { rho: z => Math.log1p(z), w: z => 1 / (1 + z) },
};

/**
 * Levenberg–Marquardt for min ½ Σ ρ(r_i(p)²): residuals(p) → array r (length m ≥ number of
 * free parameters). Use when the problem is a residual vector (fits with several data sets,
 * implicit models, bounds, outliers); for y = model(x, p) data curveFit (core/data.mjs) is the
 * simpler front end. Options:
 *   jacobian(p) → m×n rows ∂r_i/∂p_j (default: central differences, one-sided at a bound);
 *   lower/upper: box bounds, enforced by projection; parameters pinned at a bound are reported
 *     in `active` and excluded from the step and the covariance (their cov rows are NaN);
 *   loss: 'linear' (default, Gauss–Newton/LM), 'soft_l1', 'huber' or 'cauchy' with scale fScale
 *     (ρ_c(z) = c²ρ(z/c²); residuals much larger than fScale are down-weighted, IRLS weights);
 *   absoluteSigma: true when residuals are already (y − model)/σ with trusted σ → cov = (JᵀWJ)⁻¹;
 *     false (default) scales cov by 2·cost/dof (errors estimated from the scatter).
 * Damping: Marquardt diagonal scaling with Nielsen's λ update. Stops on relative step ≤ xtol,
 * relative cost decrease ≤ ftol (after an accepted step) or ‖∇cost‖∞ ≤ gtol on free variables.
 * Robust-loss covariances are the usual sandwich-free approximation (JᵀWJ)⁻¹·s²: indicative only.
 * Returns { p, cost, residuals, chi2 (Σ r²), dof, jacobian, cov (null if singular), errors,
 * active, iterations, evaluations, jacobianEvaluations, converged, reason }.
 */
export function levenbergMarquardt(residuals, p0, { jacobian, lower, upper, loss = 'linear', fScale = 1, absoluteSigma = false, maxIterations = 200, xtol = 1e-10, ftol = 1e-14, gtol = 1e-12, lambda0 = 1e-3 } = {}) {
  fn(residuals, 'levenbergMarquardt', 'residuals');
  if (jacobian !== undefined) fn(jacobian, 'levenbergMarquardt', 'jacobian');
  let p = vec(p0, 'levenbergMarquardt', 'p0'); const n = p.length, box = boxOf(lower, upper, n, 'levenbergMarquardt');
  const L = LOSSES[loss];
  if (!L) throw new Error(`levenbergMarquardt: unknown loss "${loss}" (linear, soft_l1, huber, cauchy)`);
  if (!(fScale > 0 && Number.isFinite(fScale))) throw new Error('levenbergMarquardt: fScale must be positive');
  posInt(maxIterations, 'levenbergMarquardt', 'maxIterations');
  project(p, box);
  let evaluations = 0, jacobianEvaluations = 0;
  const R = q => { evaluations++; return Float64Array.from(residuals(q)); };
  const c2 = fScale * fScale;
  const costOf = r => { let s = 0; for (const v of r) s += c2 * L.rho(v * v / c2); return s / 2; };
  // Compare component changes before summing. Subtracting two total costs can erase a
  // resolvable improvement when another, constant residual dominates the objective.
  const costReduction = (before, after) => {
    let reduction = 0;
    for (let i = 0; i < before.length; i++) {
      const a = before[i], b = after[i], difference = (a - b) * (a + b);
      if (loss === 'linear') reduction += difference / 2;
      else if (loss === 'soft_l1') reduction += difference / (Math.sqrt(1 + (a / fScale) ** 2) + Math.sqrt(1 + (b / fScale) ** 2));
      else if (loss === 'cauchy') reduction += c2 / 2 * Math.log1p(difference / (c2 + b * b));
      else reduction += c2 / 2 * (L.rho(a * a / c2) - L.rho(b * b / c2));
    }
    return reduction;
  };
  let r = R(p);
  const m = r.length;
  if (!r.every(Number.isFinite)) throw new Error('levenbergMarquardt: residuals are not finite at p0');
  if (m < 1) throw new Error('levenbergMarquardt: need at least one residual');
  let cost = costOf(r);
  const jac = q => {
    jacobianEvaluations++;
    if (jacobian) {
      const J = jacobian(q);
      if (J.length !== m || J.some(row => row.length !== n)) throw new Error('levenbergMarquardt: jacobian must be m×n');
      return J.map(row => Float64Array.from(row));
    }
    const J = Array.from({ length: m }, () => new Float64Array(n)), rq = R(q);
    for (let j = 0; j < n; j++) {
      if (box && box.lo[j] === box.hi[j]) continue;  // fixed parameter: no feasible perturbation
      const h = Math.cbrt(EPS) * Math.max(Math.abs(q[j]), 1);
      const up = q.slice(), dn = q.slice();
      const hu = box ? Math.min(h, box.hi[j] - q[j]) : h;
      const hd = box ? Math.min(h, q[j] - box.lo[j]) : h;
      if (hu === 0 && hd === 0) throw new Error('levenbergMarquardt: bounds too narrow for finite-difference Jacobian');
      up[j] += hu; dn[j] -= hd;
      const ru = hu ? R(up) : rq, rd = hd ? R(dn) : rq, span = up[j] - dn[j];
      for (let i = 0; i < m; i++) J[i][j] = (ru[i] - rd[i]) / span;
    }
    return J;
  };
  const normal = (J, rv) => {
    const w = Float64Array.from(rv, v => L.w(v * v / c2));
    const A = Array.from({ length: n }, () => new Float64Array(n)), g = new Float64Array(n);
    for (let i = 0; i < m; i++) {
      const Ji = J[i], wi = w[i];
      for (let a = 0; a < n; a++) { const t = wi * Ji[a]; g[a] += t * rv[i]; for (let b = 0; b <= a; b++) A[a][b] += t * Ji[b]; }
    }
    for (let a = 0; a < n; a++) for (let b = 0; b < a; b++) A[b][a] = A[a][b];
    return { A, g };
  };
  let lambda = lambda0, nu = 2, it = 0, reason = 'iteration budget exhausted', converged = false;
  let J = jac(p), sys = normal(J, r), active = new Array(n).fill(false);
  const activeOf = g => Array.from(p, (v, i) => box != null && (box.lo[i] === box.hi[i] || (v <= box.lo[i] && g[i] > 0) || (v >= box.hi[i] && g[i] < 0)));
  for (; it < maxIterations; it++) {
    active = activeOf(sys.g);
    const free = [];
    for (let i = 0; i < n; i++) if (!active[i]) free.push(i);
    let gmax = 0;
    for (const i of free) gmax = Math.max(gmax, Math.abs(sys.g[i]));
    if (free.length === 0 || gmax <= gtol) { converged = true; reason = 'gradient tolerance'; break; }
    const diag = free.map(i => sys.A[i][i]), dmax = Math.max(...diag);
    let accepted = false;
    while (!accepted) {
      const M = free.map((i, a) => Float64Array.from(free, (j, b) => sys.A[i][j] + (a === b ? lambda * Math.max(diag[a], 1e-12 * dmax, 1e-300) : 0)));
      let delta;
      try { delta = solveLinear(M, free.map(i => -sys.g[i])); } catch { delta = null; }
      if (delta) {
        const trial = p.slice();
        free.forEach((i, a) => { trial[i] += delta[a]; });
        project(trial, box);
        const step = Float64Array.from(trial, (v, i) => v - p[i]);
        const rt = R(trial);
        const ct = rt.every(Number.isFinite) ? costOf(rt) : Infinity;
        const reduction = Number.isFinite(ct) ? costReduction(r, rt) : -Infinity;
        if (reduction > 0) {
          const Ad = dot(step, Float64Array.from(sys.A, row => dot(row, step)));
          const predicted = -dot(sys.g, step) - Ad / 2, gain = predicted > 0 ? reduction / predicted : 0;
          lambda *= Math.max(1 / 3, 1 - (2 * gain - 1) ** 3); nu = 2;
          const smallStep = norm2(step) <= xtol * (norm2(p) + xtol), smallDecrease = reduction <= ftol * cost;
          p = trial; r = rt; cost = ct; accepted = true;
          J = jac(p); sys = normal(J, r);
          if (smallStep) { converged = true; reason = 'step tolerance'; }
          else if (smallDecrease) { converged = true; reason = 'cost tolerance'; }
          continue;
        }
        // No decrease from a step already below the resolution of p: stationary to rounding.
        if (norm2(step) <= xtol * (norm2(p) + xtol)) { converged = true; reason = 'step tolerance (no further decrease)'; break; }
      }
      lambda *= nu; nu *= 2;
      if (!(lambda < 1e20)) break;
    }
    if (converged && !accepted) break;
    if (!accepted) { reason = 'cost cannot be decreased (λ overflow); check the result'; break; }
    if (converged) { it++; break; }
  }
  active = activeOf(sys.g);
  const free = [];
  for (let i = 0; i < n; i++) if (!active[i]) free.push(i);
  let chi2 = 0;
  for (const v of r) chi2 += v * v;
  const dof = m - free.length;
  let cov = null, errors = null;
  try {
    const Af = free.map(i => Float64Array.from(free, j => sys.A[i][j]));
    const cols = free.map((_, k) => solveLinear(Af, free.map((__, a) => (a === k ? 1 : 0))));
    const s2 = absoluteSigma ? 1 : dof > 0 ? 2 * cost / dof : NaN;
    cov = Array.from({ length: n }, () => new Float64Array(n).fill(NaN));
    free.forEach((i, a) => free.forEach((j, b) => { cov[i][j] = cols[b][a] * s2; }));
    errors = Float64Array.from({ length: n }, (_, i) => Math.sqrt(cov[i][i]));
  } catch { cov = null; errors = null; }
  return { p, cost, residuals: r, chi2, dof, jacobian: J, cov, errors, active, iterations: it, evaluations, jacobianEvaluations, converged, reason };
}

// ---------------------------------------------------------------- global (stochastic) search

const reflect = (v, lo, hi) => {
  if (!(Number.isFinite(lo) && Number.isFinite(hi))) return Math.min(hi, Math.max(lo, v));
  const w = hi - lo;
  if (w === 0) return lo;
  let t = (v - lo) % (2 * w);
  if (t < 0) t += 2 * w;
  return lo + (t <= w ? t : 2 * w - t);
};

/**
 * Simulated annealing (Kirkpatrick et al. 1983) for multimodal objectives: Metropolis moves
 * with Gaussian proposals whose per-coordinate step adapts towards ~40 % acceptance, under a
 * geometric cooling T_k = T0·(Tend/T0)^{k/N} over the evaluation budget N. Bounds reflect
 * proposals into the box. T0 defaults to the spread of f over 20 trial moves; Tend = 1e-6·T0.
 * Seeded and deterministic. Heuristic: no optimality certificate. Returns the best point seen;
 * polish: true refines it with Nelder–Mead (then converged/reason refer to that polish).
 * Verify by repeating with other seeds. Returns { x, f, evaluations, acceptance, converged,
 * reason, polished? }.
 */
export function simulatedAnnealing(f, x0, { lower, upper, seed = 1, maxEvaluations = 20000, T0, Tend, step, polish = false } = {}) {
  fn(f, 'simulatedAnnealing');
  const start = vec(x0, 'simulatedAnnealing'), n = start.length, box = boxOf(lower, upper, n, 'simulatedAnnealing');
  posInt(maxEvaluations, 'simulatedAnnealing', 'maxEvaluations');
  const rand = rng(seed);
  project(start, box);
  const steps = Float64Array.from(start, (v, i) => {
    if (step != null) return typeof step === 'number' ? step : step[i];
    if (box && Number.isFinite(box.hi[i] - box.lo[i])) return 0.1 * (box.hi[i] - box.lo[i]);
    return 0.1 * Math.max(Math.abs(v), 1);
  });
  const maxStep = Float64Array.from(steps, (s, i) => (box && Number.isFinite(box.hi[i] - box.lo[i]) ? box.hi[i] - box.lo[i] : 1e3 * s));
  let evaluations = 0;
  const F = x => { evaluations++; const v = f(x); return Number.isNaN(v) ? Infinity : v; };
  const propose = (x, scale = 1) => Float64Array.from(x, (v, i) => {
    const t = v + scale * steps[i] * rand.normal();
    return box ? reflect(t, box.lo[i], box.hi[i]) : t;
  });
  let x = start, fx = finiteValue(F(x), 'simulatedAnnealing');
  let T = T0;
  if (T == null) {
    const probes = Array.from({ length: 20 }, () => F(propose(x))).filter(Number.isFinite);
    const mean = probes.reduce((a, v) => a + v, 0) / Math.max(probes.length, 1);
    T = Math.sqrt(probes.reduce((a, v) => a + (v - mean) ** 2, 0) / Math.max(probes.length, 1)) || Math.max(Math.abs(fx), 1);
  }
  if (!(T > 0 && Number.isFinite(T))) throw new Error('simulatedAnnealing: T0 must be positive');
  const T1 = Tend ?? 1e-6 * T, budget = maxEvaluations - evaluations, cool = (T1 / T) ** (1 / Math.max(budget, 1));
  let best = x.slice(), fbest = fx, accepted = 0, window = 0, windowAccepted = 0, total = 0;
  for (;;) {
    if (evaluations >= maxEvaluations) break;
    const y = propose(x), fy = F(y);
    total++; window++;
    if (fy <= fx || rand() < Math.exp(-(fy - fx) / T)) {
      x = y; fx = fy; accepted++; windowAccepted++;
      if (fx < fbest) { best = x.slice(); fbest = fx; }
    }
    T *= cool;
    if (window === 50) {
      const rate = windowAccepted / window, factor = rate > 0.5 ? 1.25 : rate < 0.3 ? 0.8 : 1;
      for (let i = 0; i < n; i++) steps[i] = Math.min(maxStep[i], Math.max(steps[i] * factor, 1e-12 * Math.max(1, Math.abs(x[i]))));
      window = 0; windowAccepted = 0;
    }
  }
  const out = { x: best, f: fbest, evaluations, acceptance: accepted / Math.max(total, 1), converged: false, reason: 'evaluation budget exhausted (heuristic: no optimality certificate)' };
  if (polish) {
    const nm = nelderMead(f, best, { lower: box?.lo, upper: box?.hi, step: Float64Array.from(steps, s => Math.max(s, 1e-6)) });
    out.evaluations += nm.evaluations;
    if (nm.f <= fbest) { out.x = nm.x; out.f = nm.f; }
    out.converged = nm.converged; out.reason = `polished by Nelder–Mead: ${nm.reason}`; out.polished = nm;
  }
  return out;
}

/**
 * Differential evolution (Storn & Price 1997) inside the box [lower, upper] (all finite):
 * strategy 'rand1bin' (default, robust) or 'best1bin' (faster, greedier), mutation factor F
 * (number, or [Fmin, Fmax] dithered per generation; default [0.5, 1]), crossover CR (0.9),
 * population popSize (default max(15·n, 20)). Seeded and deterministic. Converged when the
 * population's f spread ≤ ftol + rtol·|f_best| (default 1e-10, 1e-10): the population has
 * collapsed onto one basin, usually (not provably) the global one. polish: true refines with
 * projected L-BFGS (finite-difference gradient). Returns { x, f, generations, evaluations,
 * converged, reason, population, values }.
 */
export function differentialEvolution(f, lower, upper, { seed = 1, popSize, F = [0.5, 1], CR = 0.9, strategy = 'rand1bin', maxGenerations = 1000, ftol = 1e-10, rtol = 1e-10, polish = false } = {}) {
  fn(f, 'differentialEvolution');
  const lo = vec(lower, 'differentialEvolution', 'lower'), hi = vec(upper, 'differentialEvolution', 'upper'), n = lo.length;
  if (hi.length !== n || lo.some((v, i) => !(v <= hi[i]))) throw new Error('differentialEvolution: need finite lower ≤ upper of equal length');
  if (!['rand1bin', 'best1bin'].includes(strategy)) throw new Error('differentialEvolution: strategy must be rand1bin or best1bin');
  const NP = popSize ?? Math.max(15 * n, 20);
  posInt(NP, 'differentialEvolution', 'popSize'); posInt(maxGenerations, 'differentialEvolution', 'maxGenerations');
  if (NP < 5) throw new Error('differentialEvolution: popSize must be ≥ 5');
  if (!(CR >= 0 && CR <= 1)) throw new Error('differentialEvolution: CR must be in [0, 1]');
  const [Fmin, Fmax] = typeof F === 'number' ? [F, F] : F;
  if (!(Fmin > 0 && Fmax >= Fmin && Fmax <= 2)) throw new Error('differentialEvolution: F must be in (0, 2]');
  const rand = rng(seed);
  let evaluations = 0;
  const E = x => { evaluations++; const v = f(x); return Number.isNaN(v) ? Infinity : v; };
  // Latin-hypercube start: each coordinate's NP strata used exactly once
  const pop = Array.from({ length: NP }, () => new Float64Array(n));
  for (let k = 0; k < n; k++) {
    const perm = Array.from({ length: NP }, (_, i) => i);
    for (let i = NP - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); [perm[i], perm[j]] = [perm[j], perm[i]]; }
    for (let i = 0; i < NP; i++) pop[i][k] = lo[k] + (hi[k] - lo[k]) * (perm[i] + rand()) / NP;
  }
  const vals = pop.map(E);
  let bestIdx = 0;
  for (let i = 1; i < NP; i++) if (vals[i] < vals[bestIdx]) bestIdx = i;
  const pick = (exclude) => { let r; do r = Math.floor(rand() * NP); while (exclude.includes(r)); return r; };
  let gen = 0, converged = false, reason = 'generation budget exhausted';
  for (; gen < maxGenerations; gen++) {
    let fmin = Infinity, fmax = -Infinity;
    for (const v of vals) { fmin = Math.min(fmin, v); fmax = Math.max(fmax, v); }
    if (fmax - fmin <= ftol + rtol * Math.abs(fmin)) { converged = true; reason = 'population collapsed (f spread tolerance)'; break; }
    const Fg = Fmin + (Fmax - Fmin) * rand();
    for (let i = 0; i < NP; i++) {
      const r0 = strategy === 'best1bin' ? bestIdx : pick([i]);
      const r1 = pick([i, r0]), r2 = pick([i, r0, r1]);
      const jr = Math.floor(rand() * n), trial = pop[i].slice();
      for (let k = 0; k < n; k++) {
        if (k === jr || rand() < CR) {
          let v = pop[r0][k] + Fg * (pop[r1][k] - pop[r2][k]);
          if (v < lo[k]) v = lo[k] + rand() * (pop[i][k] - lo[k]);
          else if (v > hi[k]) v = hi[k] - rand() * (hi[k] - pop[i][k]);
          trial[k] = v;
        }
      }
      const ft = E(trial);
      if (ft <= vals[i]) { pop[i] = trial; vals[i] = ft; if (ft < vals[bestIdx]) bestIdx = i; }
    }
  }
  const out = { x: pop[bestIdx].slice(), f: vals[bestIdx], generations: gen, evaluations, converged, reason, population: pop.map(p => p.slice()), values: vals.slice() };
  if (polish) {
    const lb = lbfgs(f, out.x, { lower: lo, upper: hi, gtol: 1e-10, maxIterations: 1000 });
    out.evaluations += lb.evaluations;
    if (lb.f <= out.f) { out.x = lb.x; out.f = lb.f; }
    out.polished = lb;
  }
  return out;
}
