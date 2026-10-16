// Implicit time integration for stiff systems  M y' = f(t, y)  (M = I unless a mass matrix is given).
//
//   stiffStepper    adaptive Radau IIA (order 5), Rosenbrock Rodas3 (order 3) or variable-order BDF (1–5)
//   stiffnessRatio  eigenvalues of the Jacobian and the stiffness ratio / index
//   odeAuto         decide stiff vs non-stiff for a problem and hand back a ready stepper
//
// The stepper mirrors odeStepper (core/ode.mjs): step(hLimit), advanceTo(t, maxSteps, { exact }),
// sampleAt(t), events with the same semantics, steps / rejected counters, plus work statistics.
// "Stiff" means the fastest decaying modes of the Jacobian are much faster than the solution you
// want to follow: explicit methods are then limited by stability (h·|λ_max| ≲ 3), not accuracy.

import { brent } from './numerics.mjs';
import { odeStepper } from './ode.mjs';
import { eigenNonsymmetric } from './spectral.mjs';

const EPS = Number.EPSILON;

// ---------------------------------------------------------------- dense linear algebra (private)

/** In-place LU with partial pivoting of a real n×n matrix stored row-major in a Float64Array. */
function luReal(A, n) {
  const piv = new Int32Array(n);
  for (let k = 0; k < n; k++) {
    let p = k, best = Math.abs(A[k * n + k]);
    for (let i = k + 1; i < n; i++) { const v = Math.abs(A[i * n + k]); if (v > best) { best = v; p = i; } }
    if (!(best > 0) || !Number.isFinite(best)) return null;
    piv[k] = p;
    if (p !== k) for (let j = 0; j < n; j++) { const t = A[k * n + j]; A[k * n + j] = A[p * n + j]; A[p * n + j] = t; }
    for (let i = k + 1; i < n; i++) {
      const f = A[i * n + k] / A[k * n + k];
      A[i * n + k] = f;
      if (f !== 0) for (let j = k + 1; j < n; j++) A[i * n + j] -= f * A[k * n + j];
    }
  }
  return { A, piv, n };
}
function solveReal({ A, piv, n }, b) {
  for (let k = 0; k < n; k++) { const p = piv[k]; if (p !== k) { const t = b[k]; b[k] = b[p]; b[p] = t; } }
  for (let i = 1; i < n; i++) { let s = b[i]; for (let j = 0; j < i; j++) s -= A[i * n + j] * b[j]; b[i] = s; }
  for (let i = n - 1; i >= 0; i--) { let s = b[i]; for (let j = i + 1; j < n; j++) s -= A[i * n + j] * b[j]; b[i] = s / A[i * n + i]; }
  return b;
}
/** Complex LU (separate real / imaginary parts), partial pivoting on |·|. */
function luComplex(R, I, n) {
  const piv = new Int32Array(n);
  for (let k = 0; k < n; k++) {
    let p = k, best = Math.hypot(R[k * n + k], I[k * n + k]);
    for (let i = k + 1; i < n; i++) { const v = Math.hypot(R[i * n + k], I[i * n + k]); if (v > best) { best = v; p = i; } }
    if (!(best > 0) || !Number.isFinite(best)) return null;
    piv[k] = p;
    if (p !== k) {
      for (let j = 0; j < n; j++) {
        let t = R[k * n + j]; R[k * n + j] = R[p * n + j]; R[p * n + j] = t;
        t = I[k * n + j]; I[k * n + j] = I[p * n + j]; I[p * n + j] = t;
      }
    }
    const scale = Math.max(Math.abs(R[k * n + k]), Math.abs(I[k * n + k]));
    const dr = R[k * n + k] / scale, di = I[k * n + k] / scale, dd = dr * dr + di * di;
    for (let i = k + 1; i < n; i++) {
      const ar = R[i * n + k] / scale, ai = I[i * n + k] / scale, fr = (ar * dr + ai * di) / dd, fi = (ai * dr - ar * di) / dd;
      R[i * n + k] = fr; I[i * n + k] = fi;
      if (fr !== 0 || fi !== 0) {
        for (let j = k + 1; j < n; j++) {
          const br = R[k * n + j], bi = I[k * n + j];
          R[i * n + j] -= fr * br - fi * bi; I[i * n + j] -= fr * bi + fi * br;
        }
      }
    }
  }
  return { R, I, piv, n };
}
function solveComplex({ R, I, piv, n }, br, bi) {
  for (let k = 0; k < n; k++) {
    const p = piv[k];
    if (p !== k) { let t = br[k]; br[k] = br[p]; br[p] = t; t = bi[k]; bi[k] = bi[p]; bi[p] = t; }
  }
  for (let i = 1; i < n; i++) {
    let sr = br[i], si = bi[i];
    for (let j = 0; j < i; j++) { const lr = R[i * n + j], li = I[i * n + j]; sr -= lr * br[j] - li * bi[j]; si -= lr * bi[j] + li * br[j]; }
    br[i] = sr; bi[i] = si;
  }
  for (let i = n - 1; i >= 0; i--) {
    let sr = br[i], si = bi[i];
    for (let j = i + 1; j < n; j++) { const ur = R[i * n + j], ui = I[i * n + j]; sr -= ur * br[j] - ui * bi[j]; si -= ur * bi[j] + ui * br[j]; }
    const scale = Math.max(Math.abs(R[i * n + i]), Math.abs(I[i * n + i]));
    const dr = R[i * n + i] / scale, di = I[i * n + i] / scale, dd = dr * dr + di * di;
    sr /= scale; si /= scale;
    br[i] = (sr * dr + si * di) / dd; bi[i] = (si * dr - sr * di) / dd;
  }
}

// ---------------------------------------------------------------- Jacobians

const toRows = (J, n, who) => {
  if (!J || J.length !== n) throw new Error(`${who}: jac(t, y) must return ${n} rows of ${n} numbers`);
  for (const row of J) if (!row || row.length !== n) throw new Error(`${who}: jac(t, y) must return ${n} rows of ${n} numbers`);
  return J;
};

/** Forward-difference Jacobian ∂f/∂y (column j perturbed by √eps·max(1e-5, |y_j|)). */
function fdJacobian(f, t, y, f0, n) {
  const baseline = Float64Array.from(f0);      // callbacks may reuse their output buffer
  const J = Array.from({ length: n }, () => new Float64Array(n)), yy = Float64Array.from(y);
  for (let j = 0; j < n; j++) {
    yy[j] = y[j] + Math.sqrt(EPS) * Math.max(1e-5, Math.abs(y[j]));
    const d = yy[j] - y[j];
    const fj = f(t, yy);
    for (let i = 0; i < n; i++) J[i][j] = (fj[i] - baseline[i]) / d;
    yy[j] = y[j];
  }
  return J;
}

/**
 * Eigenvalues of the Jacobian J = ∂f/∂y and the derived stiffness measures. Pass either a matrix
 * (array of rows) or { f, t, y, jac? } to evaluate J (finite differences when jac is absent).
 * Returns { eigenvalues: [{ re, im }], spectralRadius, fastest, slowest, ratio }:
 *   fastest = max(−Re λ) over decaying modes, slowest = min(−Re λ) over decaying modes whose
 *   −Re λ > zeroTolerance·spectralRadius (neutral modes such as conservation laws are skipped),
 *   ratio = fastest/slowest (the classical stiffness ratio; 1 if fewer than one decaying mode).
 * The ratio is a local diagnostic at one state: stiffness can appear later (see odeAuto).
 */
export function stiffnessRatio(input, { zeroTolerance = 1e-9 } = {}) {
  let J = input;
  if (!Array.isArray(input) && input && typeof input.f === 'function') {
    const { f, t = 0, y, jac } = input, n = y?.length;
    if (!n) throw new Error('stiffnessRatio: y must be a non-empty array');
    J = jac ? toRows(jac(t, Float64Array.from(y)), n, 'stiffnessRatio') : fdJacobian(f, t, y, Float64Array.from(f(t, Float64Array.from(y))), n);
  }
  const n = J?.length;
  if (!n || !Array.from(J).every(row => row?.length === n && Array.from(row).every(Number.isFinite))) throw new Error('stiffnessRatio: need a finite square matrix or { f, y }');
  const eigenvalues = eigenNonsymmetric(J);
  const spectralRadius = Math.max(...eigenvalues.map(e => Math.hypot(e.re, e.im)));
  const decay = eigenvalues.map(e => -e.re).filter(d => d > zeroTolerance * spectralRadius);
  const fastest = decay.length ? Math.max(...decay) : 0, slowest = decay.length ? Math.min(...decay) : 0;
  return { eigenvalues, spectralRadius, fastest, slowest, ratio: decay.length ? fastest / slowest : 1 };
}

// ---------------------------------------------------------------- Radau IIA constants (Hairer & Wanner RADAU5)

const SQ6 = Math.sqrt(6);
const C1 = (4 - SQ6) / 10, C2 = (4 + SQ6) / 10, C1M1 = C1 - 1, C2M1 = C2 - 1, C1MC2 = C1 - C2;
const DD1 = -(13 + 7 * SQ6) / 3, DD2 = (-13 + 7 * SQ6) / 3, DD3 = -1 / 3;
const RADAU = (() => {
  const c81 = Math.cbrt(81), c9 = Math.cbrt(9);
  const u1 = (6 + c81 - c9) / 30, alph0 = (12 - c81 + c9) / 60, beta0 = (c81 + c9) * Math.sqrt(3) / 60, cno = alph0 * alph0 + beta0 * beta0;
  return { u1: 1 / u1, alph: alph0 / cno, beta: beta0 / cno };
})();
const T11 = 0.09123239487089295, T12 = -0.1412552950209542, T13 = -0.030029194105147424;
const T21 = 0.241717932707107, T22 = 0.20412935229379994, T23 = 0.3829421127572619;
const T31 = 0.966048182615093;
const TI11 = 4.325579890063155, TI12 = 0.33919925181580984, TI13 = 0.5417705399358749;
const TI21 = -4.178718591551905, TI22 = -0.32768282076106237, TI23 = 0.47662355450055044;
const TI31 = -0.5028726349457868, TI32 = 2.571926949855605, TI33 = -0.5960392048282249;

// Rodas3 (Sandu et al. 1997): stiffly accurate, L-stable, order 3 with embedded order 2.
const ROS = { gamma: 0.5, A: [[], [0], [2, 0], [2, 0, 1]], C: [[], [4], [1, -1], [1, -1, -8 / 3]], M: [2, 0, 1, 1], E: [0, 0, 0, 1], alpha: [0, 0, 1, 1], gammas: [0.5, 1.5, 0, 0], newF: [true, false, true, true] };

const BDF_MAX = 5;
const BDF_GAMMA = [0, 1, 1.5, 11 / 6, 25 / 12, 137 / 60];
const BDF_ERR = [1, 1 / 2, 1 / 3, 1 / 4, 1 / 5, 1 / 6, 1 / 7];

/** Matrix R (order+1)² of the backward-difference step-change transformation (Shampine & Reichelt 1997). */
function bdfR(order, factor) {
  const M = Array.from({ length: order + 1 }, () => new Float64Array(order + 1));
  for (let j = 0; j <= order; j++) M[0][j] = 1;
  for (let i = 1; i <= order; i++) for (let j = 1; j <= order; j++) M[i][j] = (i - 1 - factor * j) / i;
  for (let i = 1; i <= order; i++) for (let j = 0; j <= order; j++) M[i][j] *= M[i - 1][j];
  return M;
}
function bdfChangeD(D, order, factor, n) {
  const R = bdfR(order, factor), U = bdfR(order, 1);
  const RU = R.map(row => Float64Array.from({ length: order + 1 }, (_, j) => row.reduce((s, r, k) => s + r * U[k][j], 0)));
  const out = Array.from({ length: order + 1 }, () => new Float64Array(n));
  for (let j = 0; j <= order; j++) for (let k = 0; k <= order; k++) { const c = RU[k][j]; if (c) for (let i = 0; i < n; i++) out[j][i] += c * D[k][i]; }
  for (let j = 0; j <= order; j++) D[j].set(out[j]);
}

/**
 * Adaptive implicit stepper for stiff systems M y′ = f(t, y).
 *   method 'radau5'     Radau IIA, 3 stages, order 5, L-stable, simplified Newton with the real/complex
 *                       splitting of Hairer & Wanner's RADAU5; supports a constant mass matrix (also
 *                       singular: index-1 DAEs, with consistent initial values). Best at tight tolerances.
 *   method 'rosenbrock' Rodas3 (Sandu et al. 1997), order 3(2), L-stable, stiffly accurate, linearly
 *                       implicit (one LU per step, no Newton); needs an accurate Jacobian every step.
 *                       Good at loose tolerances (rtol ≳ 1e-6).
 *   method 'bdf'        variable-order (1…maxOrder ≤ 5), variable-step BDF in backward-difference form
 *                       (Shampine & Reichelt 1997, as in MATLAB ode15s / SciPy BDF). Cheapest per step for
 *                       large smooth problems; BDF3–5 are not A-stable (avoid for eigenvalues near the
 *                       imaginary axis: use radau5).
 * f(t, y) returns y′ (array-like); jac(t, y) returns ∂f/∂y as n rows (else forward differences, n extra
 * f evaluations per Jacobian, ~1e-8 relative accuracy — fine for radau5/bdf, may limit rosenbrock at
 * tight tolerances). rtol, atol (number or per-component array) set the local error per step;
 * h is the first trial step (default: estimated). hmax caps steps (mutable s.hmax, as in odeStepper).
 * mass: constant n×n matrix (radau5 only). adaptive: false gives fixed steps h (order studies; the
 * error estimate is ignored; Newton failure throws). bdf option maxOrder (default 5).
 * events: [{ g(t, y), direction: 0 | +1 | −1, terminal }] — sign changes located with Brent on the dense
 * output exactly as in odeStepper; terminal events stop the integration (stopped = true).
 * Dense output: radau5 its collocation polynomial (order 3 uniformly, 5 at step ends); bdf its
 * interpolating polynomial; rosenbrock cubic Hermite from f at both ends (not available with a
 * singular mass matrix — hence mass is radau5-only).
 * Failure: throws on non-finite states, step-size underflow, or 60 consecutive rejections.
 * Returns { t, y, hmax, method, step(hLimit?), advanceTo(t, maxSteps, { exact }), sampleAt(t, out),
 * steps, rejected, events, stopped, stats: { fevals, jacobians, decompositions, solves, newton, order } }.
 */
export function stiffStepper({ f, jac, y0, t0 = 0, method = 'radau5', h, rtol = 1e-6, atol = 1e-9, hmax = Infinity, events: watch = [], mass, adaptive = true, maxOrder = BDF_MAX }) {
  const who = 'stiffStepper';
  if (typeof f !== 'function') throw new Error(`${who}: f(t, y) is required`);
  if (!y0?.length || !Array.from(y0).every(Number.isFinite)) throw new Error(`${who}: y0 must be non-empty and finite`);
  if (!['radau5', 'rosenbrock', 'bdf'].includes(method)) throw new Error(`${who}: unknown method "${method}" (radau5 | rosenbrock | bdf)`);
  if (!(rtol > 0 && Number.isFinite(rtol))) throw new Error(`${who}: rtol must be positive`);
  if (!Number.isFinite(t0)) throw new Error(`${who}: t0 must be finite`);
  if (!(hmax > 0 && (Number.isFinite(hmax) || hmax === Infinity))) throw new Error(`${who}: hmax must be positive`);
  const n = y0.length;
  const atolV = typeof atol === 'number' ? new Float64Array(n).fill(atol) : Float64Array.from(atol);
  if (atolV.length !== n || !atolV.every(v => v >= 0 && Number.isFinite(v))) throw new Error(`${who}: atol must be ≥ 0 (number or one per component)`);
  if (!adaptive && !(h > 0 && Number.isFinite(h))) throw new Error(`${who}: adaptive: false needs a positive step h`);
  if (h !== undefined && !(h > 0 && Number.isFinite(h))) throw new Error(`${who}: h must be positive`);
  if (!(Number.isInteger(maxOrder) && maxOrder >= 1 && maxOrder <= BDF_MAX)) throw new Error(`${who}: maxOrder must be 1…5`);
  let M = null;
  if (mass) {
    if (method !== 'radau5') throw new Error(`${who}: a mass matrix is supported by method 'radau5' only`);
    if (mass.length !== n || !Array.from(mass).every(r => r?.length === n && Array.from(r).every(Number.isFinite))) throw new Error(`${who}: mass must be a finite ${n}×${n} matrix`);
    M = Array.from(mass, r => Float64Array.from(r));
  }

  const stats = { fevals: 0, jacobians: 0, decompositions: 0, solves: 0, newton: 0, order: method === 'radau5' ? 5 : method === 'rosenbrock' ? 3 : 1 };
  const F = (t, y) => {
    stats.fevals++;
    const v = f(t, y);
    if (!v || v.length !== n) throw new Error(`${who}: f(t, y) must return ${n} values`);
    // Stage derivatives must remain distinct even when the callback reuses a work array.
    return Float64Array.from(v);
  };
  const jacobian = (t, y, fy) => {
    stats.jacobians++;
    if (jac) return Array.from(toRows(jac(t, y), n, who), row => Float64Array.from(row));
    const J = fdJacobian(f, t, y, fy, n);
    stats.fevals += n;
    return J;
  };
  const checkFinite = (v, t) => { for (let i = 0; i < n; i++) if (!Number.isFinite(v[i])) throw new Error(`${who}: non-finite state or derivative at t = ${t}`); };
  const Mv = (v, out) => { if (!M) { out.set(v); return out; } for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += M[i][j] * v[j]; out[i] = s; } return out; };
  const rms = (e, sc) => { let s = 0; for (let i = 0; i < n; i++) { const r = e[i] / sc[i]; s += r * r; } return Math.sqrt(s / n); };

  const y = Float64Array.from(y0);
  const f0 = Float64Array.from(F(t0, y));
  checkFinite(f0, t0);
  const prev = { t: t0, y: y.slice(), dy: f0.slice() };

  // initial step: h ≈ 0.01·‖y‖/‖f‖ in the weighted norm (Hairer, Nørsett & Wanner I, II.4)
  let hProp = h;
  if (hProp === undefined) {
    const sc = Float64Array.from(y, (v, i) => atolV[i] + rtol * Math.abs(v));
    const d0 = rms(y, sc), d1 = rms(f0, sc);
    hProp = d0 < 1e-5 || d1 < 1e-5 ? 1e-6 : 0.01 * d0 / d1;
    hProp = Math.min(hProp, hmax);
  }

  const s = {
    t: t0, y, steps: 0, rejected: 0, method, events: [], stopped: false, hmax, stats,
    step(hLimit = Infinity) {
      if (s.stopped) return;
      if (!(hLimit > 0) || !(s.hmax > 0)) throw new Error(`${who}: hLimit and hmax must be positive`);
      prev.t = s.t; prev.y.set(y); prev.dy.set(f0);
      impl.step(hLimit);
      checkFinite(y, s.t);
      s.steps++;
      if (watch.length) detectEvents();
    },
    advanceTo(target, maxSteps = Infinity, { exact = false } = {}) {
      let i = 0;
      while (s.t < target && !s.stopped) {
        if (i++ >= maxSteps) return false;
        s.step(exact ? target - s.t : Infinity);
        if (exact && Math.abs(s.t - target) <= 4 * EPS * Math.max(1, Math.abs(target))) s.t = target;
      }
      return true;
    },
    /** State at time t within the last step (method's dense output). */
    sampleAt(t, out = new Float64Array(n)) { return impl.sample(t, out); },
  };

  const newStep = (hh, who2) => { if (s.t + hh === s.t) throw new Error(`${who}: step size underflow at t = ${s.t} (${who2})`); };

  // ================================================================ Radau IIA (RADAU5)
  function makeRadau() {
    const quot = rtol > 0 ? atolV.map(a => a / rtol) : atolV;
    const rtolT = 0.1 * Math.pow(rtol, 2 / 3), atolT = quot.map(q => rtolT * q);
    const fnewt = Math.max(10 * EPS / rtolT, Math.min(0.03, Math.sqrt(rtolT)));
    const nit = 7, safe = 0.9, thet = 0.001, facl = 5, facr = 1 / 8;
    const scal = new Float64Array(n);
    const setScal = () => { for (let i = 0; i < n; i++) scal[i] = atolT[i] + rtolT * Math.abs(y[i]); };
    setScal();
    let J = null, jacCurrent = false, needJac = true, luH = NaN, E1 = null, E2 = null;
    let first = true, reject = false, theta = thet, faccon = 1, hold = hProp, hacc = 0, erracc = 0, naccpt = 0;
    const z1 = new Float64Array(n), z2 = new Float64Array(n), z3 = new Float64Array(n);
    const g1 = new Float64Array(n), g2 = new Float64Array(n), g3 = new Float64Array(n);   // transformed variables W = T⁻¹Z
    const tmp = new Float64Array(n), tmp2 = new Float64Array(n), tmp3 = new Float64Array(n), est = new Float64Array(n);
    const cont = [y.slice(), new Float64Array(n), new Float64Array(n), new Float64Array(n)];
    let contT = s.t, contH = hProp;

    function factor(hh) {
      const fac1 = RADAU.u1 / hh, alphn = RADAU.alph / hh, betan = RADAU.beta / hh;
      const A = new Float64Array(n * n), R = new Float64Array(n * n), I = new Float64Array(n * n);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          const m = M ? M[i][j] : (i === j ? 1 : 0);
          A[i * n + j] = fac1 * m - J[i][j];
          R[i * n + j] = alphn * m - J[i][j];
          I[i * n + j] = betan * m;
        }
      }
      E1 = luReal(A, n); E2 = luComplex(R, I, n);
      stats.decompositions += 2;
      if (!E1 || !E2) throw new Error(`${who}: singular iteration matrix at t = ${s.t}`);
      luH = hh;
    }

    function step(hLimit) {
      for (let attempt = 0; ; attempt++) {
        if (attempt > 60) throw new Error(`${who}: cannot meet tolerance at t = ${s.t}`);
        const hh = Math.min(hProp, s.hmax, hLimit);
        const clipped = adaptive && hLimit < Math.min(hProp, s.hmax), keep = hProp;
        newStep(hh, 'radau5');
        if (needJac) { J = jacobian(s.t, y, f0); jacCurrent = true; needJac = false; luH = NaN; }
        if (luH !== hh) factor(hh);
        const fac1 = RADAU.u1 / hh, alphn = RADAU.alph / hh, betan = RADAU.beta / hh;
        // starting values: extrapolate the previous collocation polynomial
        if (first) { z1.fill(0); z2.fill(0); z3.fill(0); g1.fill(0); g2.fill(0); g3.fill(0); } else {
          const c3q = hh / hold, c1q = C1 * c3q, c2q = C2 * c3q;
          for (let i = 0; i < n; i++) {
            const ak1 = cont[1][i], ak2 = cont[2][i], ak3 = cont[3][i];
            const a1 = c1q * (ak1 + (c1q - C2M1) * (ak2 + (c1q - C1M1) * ak3));
            const a2 = c2q * (ak1 + (c2q - C2M1) * (ak2 + (c2q - C1M1) * ak3));
            const a3 = c3q * (ak1 + (c3q - C2M1) * (ak2 + (c3q - C1M1) * ak3));
            z1[i] = a1; z2[i] = a2; z3[i] = a3;
            g1[i] = TI11 * a1 + TI12 * a2 + TI13 * a3;
            g2[i] = TI21 * a1 + TI22 * a2 + TI23 * a3;
            g3[i] = TI31 * a1 + TI32 * a2 + TI33 * a3;
          }
        }
        // simplified Newton iteration
        faccon = Math.pow(Math.max(faccon, EPS), 0.8);
        theta = Math.abs(thet);
        let newt = 0, dynold = 0, thqold = 0, failure = null;
        for (;;) {
          if (newt >= nit) { failure = 'newton'; break; }
          for (let i = 0; i < n; i++) tmp[i] = y[i] + z1[i];
          const fz1 = F(s.t + C1 * hh, tmp);
          for (let i = 0; i < n; i++) tmp[i] = y[i] + z2[i];
          const fz2 = F(s.t + C2 * hh, tmp);
          for (let i = 0; i < n; i++) tmp[i] = y[i] + z3[i];
          const fz3 = F(s.t + hh, tmp);
          for (let i = 0; i < n; i++) {
            const a1 = fz1[i], a2 = fz2[i], a3 = fz3[i];
            z1[i] = TI11 * a1 + TI12 * a2 + TI13 * a3;
            z2[i] = TI21 * a1 + TI22 * a2 + TI23 * a3;
            z3[i] = TI31 * a1 + TI32 * a2 + TI33 * a3;
          }
          Mv(g1, tmp); Mv(g2, tmp2); Mv(g3, tmp3);
          for (let i = 0; i < n; i++) {
            z1[i] -= fac1 * tmp[i];
            z2[i] -= alphn * tmp2[i] - betan * tmp3[i];
            z3[i] -= alphn * tmp3[i] + betan * tmp2[i];
          }
          solveReal(E1, z1); solveComplex(E2, z2, z3);
          stats.solves += 2; stats.newton++;
          newt++;
          let dyno = 0;
          for (let i = 0; i < n; i++) dyno += (z1[i] / scal[i]) ** 2 + (z2[i] / scal[i]) ** 2 + (z3[i] / scal[i]) ** 2;
          dyno = Math.sqrt(dyno / (3 * n));
          if (!Number.isFinite(dyno)) { failure = 'newton'; break; }
          if (newt > 1 && newt < nit) {
            const thq = dyno / dynold;
            theta = newt === 2 ? thq : Math.sqrt(thq * thqold);
            thqold = thq;
            if (theta < 0.99) {
              faccon = theta / (1 - theta);
              const dyth = faccon * dyno * Math.pow(theta, nit - 1 - newt) / fnewt;
              if (dyth >= 1) {
                const qnewt = Math.max(1e-4, Math.min(20, dyth));
                failure = { hhfac: 0.8 * Math.pow(qnewt, -1 / (4 + nit - 1 - newt)) };
                break;
              }
            } else { failure = 'newton'; break; }
          }
          dynold = Math.max(dyno, EPS);
          for (let i = 0; i < n; i++) {
            const w1 = g1[i] + z1[i], w2 = g2[i] + z2[i], w3 = g3[i] + z3[i];
            g1[i] = w1; g2[i] = w2; g3[i] = w3;
            z1[i] = T11 * w1 + T12 * w2 + T13 * w3;
            z2[i] = T21 * w1 + T22 * w2 + T23 * w3;
            z3[i] = T31 * w1 + w2;
          }
          if (faccon * dyno <= fnewt) break;
        }
        if (failure) {
          if (!adaptive) throw new Error(`${who}: Newton iteration failed with fixed step h = ${hh} at t = ${s.t}`);
          hProp = hh * (failure === 'newton' ? 0.5 : failure.hhfac);
          reject = true; s.rejected++;
          if (!jacCurrent) needJac = true;
          continue;
        }
        // error estimate (Hairer & Wanner IV.8, with the "first/reject" refinement)
        const hee1 = DD1 / hh, hee2 = DD2 / hh, hee3 = DD3 / hh;
        for (let i = 0; i < n; i++) tmp[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
        Mv(tmp, tmp2);
        for (let i = 0; i < n; i++) est[i] = tmp2[i] + f0[i];
        solveReal(E1, est); stats.solves++;
        let err = Math.max(rms(est, scal), 1e-10);
        if (err >= 1 && (first || reject)) {
          for (let i = 0; i < n; i++) tmp[i] = y[i] + est[i];
          const fe = F(s.t, tmp);
          for (let i = 0; i < n; i++) est[i] = fe[i] + tmp2[i];
          solveReal(E1, est); stats.solves++;
          err = Math.max(rms(est, scal), 1e-10);
        }
        if (!Number.isFinite(err)) err = 1e10;
        if (!adaptive) err = Math.min(err, 0.5);
        const fac = Math.min(safe, safe * (1 + 2 * nit) / (newt + 2 * nit));
        let q = Math.max(facr, Math.min(facl, Math.pow(err, 0.25) / fac));
        let hnew = hh / q;
        if (err < 1) {
          first = false; naccpt++;
          if (naccpt > 1) {
            let facgus = (hacc / hh) * Math.pow(err * err / erracc, 0.25) / safe;
            facgus = Math.max(facr, Math.min(facl, facgus));
            q = Math.max(q, facgus);
            hnew = hh / q;
          }
          hacc = hh; erracc = Math.max(1e-2, err);
          hold = hh;
          for (let i = 0; i < n; i++) {
            y[i] += z3[i];
            const c2 = (z2[i] - z3[i]) / C2M1, ak = (z1[i] - z2[i]) / C1MC2;
            let acont3 = z1[i] / C1;
            acont3 = (ak - acont3) / C2;
            const c3 = (ak - c2) / C1M1;
            cont[1][i] = c2; cont[2][i] = c3; cont[3][i] = c3 - acont3;
          }
          s.t += hh;
          cont[0].set(y); contT = s.t; contH = hh;
          f0.set(F(s.t, y));
          setScal();
          if (!adaptive) { hProp = h; jacCurrent = false; needJac = theta > thet; return; }
          hnew = Math.min(hnew, s.hmax);
          if (reject) hnew = Math.min(hnew, hh);
          reject = false;
          jacCurrent = false;
          if (clipped) { hProp = Math.max(hnew, keep); needJac = theta > thet; return; }
          const qt = hnew / hh;
          if (theta <= thet && qt >= 1 && qt <= 1.2) { hProp = hh; return; }   // keep h, J and LU
          hProp = hnew;
          needJac = theta > thet;
          return;
        }
        if (!adaptive) throw new Error(`${who}: unreachable`);
        reject = true; s.rejected++;
        hProp = first ? hh * 0.1 : hnew;
        if (!jacCurrent) needJac = true;
      }
    }
    function sample(t, out) {
      const sv = (t - contT) / contH;
      for (let i = 0; i < n; i++) out[i] = cont[0][i] + sv * (cont[1][i] + (sv - C2M1) * (cont[2][i] + (sv - C1M1) * cont[3][i]));
      return out;
    }
    return { step, sample };
  }

  // ================================================================ Rosenbrock (Rodas3)
  function makeRosenbrock() {
    const K = Array.from({ length: 4 }, () => new Float64Array(n));
    const ft = new Float64Array(n), yt = new Float64Array(n), ynew = new Float64Array(n), err = new Float64Array(n), sc = new Float64Array(n);
    let J = null, jacT = NaN, reject = false;
    function step(hLimit) {
      for (let attempt = 0; ; attempt++) {
        if (attempt > 60) throw new Error(`${who}: cannot meet tolerance at t = ${s.t}`);
        const hh = Math.min(hProp, s.hmax, hLimit);
        const clipped = adaptive && hLimit < Math.min(hProp, s.hmax), keep = hProp;
        newStep(hh, 'rosenbrock');
        if (jacT !== s.t) {
          J = jacobian(s.t, y, f0);
          const dt = Math.sqrt(EPS) * Math.max(Math.abs(s.t), hh);
          const f1 = F(s.t + dt, y);
          for (let i = 0; i < n; i++) ft[i] = (f1[i] - f0[i]) / dt;
          jacT = s.t;
        }
        const A = new Float64Array(n * n), g = 1 / (hh * ROS.gamma);
        for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) A[i * n + j] = (i === j ? g : 0) - J[i][j];
        const LU = luReal(A, n);
        stats.decompositions++;
        if (!LU) throw new Error(`${who}: singular iteration matrix at t = ${s.t}`);
        let fs = f0;
        for (let st = 0; st < 4; st++) {
          if (st > 0 && ROS.newF[st]) {
            for (let i = 0; i < n; i++) { let v = y[i]; for (let j = 0; j < st; j++) v += ROS.A[st][j] * K[j][i]; yt[i] = v; }
            fs = F(s.t + ROS.alpha[st] * hh, yt);
          }
          const k = K[st];
          for (let i = 0; i < n; i++) {
            let v = fs[i] + hh * ROS.gammas[st] * ft[i];
            for (let j = 0; j < st; j++) v += ROS.C[st][j] / hh * K[j][i];
            k[i] = v;
          }
          solveReal(LU, k); stats.solves++;
        }
        for (let i = 0; i < n; i++) {
          let v = y[i], e = 0;
          for (let j = 0; j < 4; j++) { v += ROS.M[j] * K[j][i]; e += ROS.E[j] * K[j][i]; }
          ynew[i] = v; err[i] = e;
          sc[i] = atolV[i] + rtol * Math.max(Math.abs(y[i]), Math.abs(v));
        }
        let e = rms(err, sc);
        if (!Number.isFinite(e)) e = 1e10;
        if (!adaptive) e = 0;
        let fac = e === 0 ? 6 : Math.min(6, Math.max(0.2, 0.9 * Math.pow(e, -1 / 3)));
        if (e <= 1) {
          if (reject) fac = Math.min(1, fac);
          reject = false;
          y.set(ynew); s.t += hh;
          f0.set(F(s.t, y));
          hProp = adaptive ? hh * fac : h;
          if (clipped) hProp = Math.max(hProp, keep);
          return;
        }
        reject = true; s.rejected++;
        hProp = hh * Math.max(0.1, fac);
      }
    }
    function sample(t, out) {
      const t0_ = prev.t, t1 = s.t, hs = t1 - t0_;
      if (!(hs !== 0)) { out.set(y); return out; }
      const u = Math.min(1, Math.max(0, (t - t0_) / hs)), u2 = u * u, u3 = u2 * u;
      const a = 2 * u3 - 3 * u2 + 1, b = u3 - 2 * u2 + u, c = -2 * u3 + 3 * u2, d = u3 - u2;
      for (let i = 0; i < n; i++) out[i] = a * prev.y[i] + b * hs * prev.dy[i] + c * y[i] + d * hs * f0[i];
      return out;
    }
    return { step, sample };
  }

  // ================================================================ BDF (variable order, backward differences)
  function makeBdf() {
    const newtonTol = Math.max(10 * EPS / rtol, Math.min(0.03, Math.sqrt(rtol)));
    const NEWTON_MAXITER = adaptive ? 4 : 30;
    const D = Array.from({ length: BDF_MAX + 3 }, () => new Float64Array(n));
    D[0].set(y);
    for (let i = 0; i < n; i++) D[1][i] = f0[i] * hProp;
    let hAbs = hProp, order = 1, nEqual = 0, LU = null, J = jacobian(s.t, y, f0), jacCurrent = true;
    const scale = new Float64Array(n), ypred = new Float64Array(n), psi = new Float64Array(n), d = new Float64Array(n), ynew = new Float64Array(n), dy = new Float64Array(n), rhs = new Float64Array(n);
    let denseH = hAbs, denseOrder = order;

    function newton(tNew, c) {
      d.fill(0); ynew.set(ypred);
      let normOld = null, converged = false, k = 0;
      for (; k < NEWTON_MAXITER; k++) {
        const fv = F(tNew, ynew);
        let ok = true;
        for (let i = 0; i < n; i++) if (!Number.isFinite(fv[i])) ok = false;
        if (!ok) break;
        for (let i = 0; i < n; i++) rhs[i] = c * fv[i] - psi[i] - d[i];
        dy.set(rhs); solveReal(LU, dy); stats.solves++; stats.newton++;
        const nrm = rms(dy, scale);
        const rate = normOld === null ? null : nrm / normOld;
        if (rate !== null && (rate >= 1 || Math.pow(rate, NEWTON_MAXITER - k) / (1 - rate) * nrm > newtonTol)) break;
        for (let i = 0; i < n; i++) { ynew[i] += dy[i]; d[i] += dy[i]; }
        if (nrm === 0 || (rate !== null && rate / (1 - rate) * nrm < newtonTol)) { converged = true; break; }
        normOld = nrm;
      }
      return { converged, iterations: k + 1 };
    }

    function step(hLimit) {
      if (hAbs > s.hmax) { bdfChangeD(D, order, s.hmax / hAbs, n); hAbs = s.hmax; nEqual = 0; LU = null; }
      let unclipped = 0;
      for (let attempt = 0; ; attempt++) {
        if (attempt > 60) throw new Error(`${who}: cannot meet tolerance at t = ${s.t}`);
        let hh = hAbs;
        if (hh > hLimit) { unclipped = Math.max(unclipped, hh); bdfChangeD(D, order, hLimit / hh, n); hh = hAbs = hLimit; nEqual = 0; LU = null; }
        newStep(hh, 'bdf');
        const tNew = s.t + hh;
        for (let i = 0; i < n; i++) { let v = 0; for (let j = 0; j <= order; j++) v += D[j][i]; ypred[i] = v; scale[i] = atolV[i] + rtol * Math.abs(v); }
        for (let i = 0; i < n; i++) { let v = 0; for (let j = 1; j <= order; j++) v += BDF_GAMMA[j] * D[j][i]; psi[i] = v / BDF_GAMMA[order]; }
        const c = hh / BDF_GAMMA[order];
        let res;
        for (;;) {
          if (LU === null || LU.c !== c) {
            const A = new Float64Array(n * n);
            for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) A[i * n + j] = (i === j ? 1 : 0) - c * J[i][j];
            LU = luReal(A, n); stats.decompositions++;
            if (!LU) throw new Error(`${who}: singular iteration matrix at t = ${s.t}`);
            LU.c = c;
          }
          res = newton(tNew, c);
          if (res.converged || jacCurrent) break;
          J = jacobian(tNew, ypred, F(tNew, ypred)); jacCurrent = true; LU = null;
        }
        if (!res.converged) {
          if (!adaptive) throw new Error(`${who}: Newton iteration failed with fixed step h = ${hh} at t = ${s.t}`);
          bdfChangeD(D, order, 0.5, n); hAbs = hh * 0.5; nEqual = 0; LU = null; s.rejected++;
          continue;
        }
        const safety = 0.9 * (2 * NEWTON_MAXITER + 1) / (2 * NEWTON_MAXITER + res.iterations);
        for (let i = 0; i < n; i++) scale[i] = atolV[i] + rtol * Math.abs(ynew[i]);
        for (let i = 0; i < n; i++) rhs[i] = BDF_ERR[order] * d[i];
        const errNorm = adaptive ? rms(rhs, scale) : 0;
        if (!(errNorm <= 1)) {
          const factor = Number.isFinite(errNorm) ? Math.max(0.2, safety * Math.pow(errNorm, -1 / (order + 1))) : 0.2;
          bdfChangeD(D, order, factor, n); hAbs = hh * factor; nEqual = 0; s.rejected++;
          continue;
        }
        // accepted
        nEqual++;
        s.t = tNew; y.set(ynew); jacCurrent = false;
        for (let i = 0; i < n; i++) { D[order + 2][i] = d[i] - D[order + 1][i]; D[order + 1][i] = d[i]; }
        for (let j = order; j >= 0; j--) for (let i = 0; i < n; i++) D[j][i] += D[j + 1][i];
        denseH = hh; denseOrder = order;
        f0.set(F(s.t, y));
        const restore = () => {                                   // undo an artificial clip (exact advance)
          if (adaptive && unclipped > hAbs) { bdfChangeD(D, order, unclipped / hAbs, n); hAbs = unclipped; nEqual = 0; LU = null; denseH = hAbs; denseOrder = order; }
        };
        if (nEqual < order + 1) { restore(); return; }
        if (!adaptive) {
          if (order < maxOrder) { order++; nEqual = 0; LU = null; }
          stats.order = order; denseH = hAbs; denseOrder = Math.min(denseOrder, order);
          return;
        }
        const errM = order > 1 ? (() => { for (let i = 0; i < n; i++) rhs[i] = BDF_ERR[order - 1] * D[order][i]; return rms(rhs, scale); })() : Infinity;
        const errP = order < maxOrder ? (() => { for (let i = 0; i < n; i++) rhs[i] = BDF_ERR[order + 1] * D[order + 2][i]; return rms(rhs, scale); })() : Infinity;
        const norms = [errM, errNorm, errP], factors = norms.map((e, k) => (e === 0 ? Infinity : Math.pow(e, -1 / (order + k))));
        let best = 0;
        for (let k = 1; k < 3; k++) if (factors[k] > factors[best]) best = k;
        order += best - 1;
        stats.order = order;
        const factor = Math.min(10, safety * factors[best]);
        bdfChangeD(D, order, factor, n);
        hAbs = hh * factor; nEqual = 0; LU = null;
        denseH = hAbs; denseOrder = order;
        restore();
        return;
      }
    }
    function sample(t, out) {
      out.set(D[0]);
      let p = 1;
      for (let j = 0; j < denseOrder; j++) {
        p *= (t - (s.t - denseH * j)) / (denseH * (1 + j));
        for (let i = 0; i < n; i++) out[i] += D[j + 1][i] * p;
      }
      return out;
    }
    return { step, sample };
  }

  const impl = method === 'radau5' ? makeRadau() : method === 'rosenbrock' ? makeRosenbrock() : makeBdf();

  const probe = new Float64Array(n);
  function detectEvents() {
    let firstHit = null;
    watch.forEach((e, index) => {
      const g0 = e.g(prev.t, prev.y), g1 = e.g(s.t, y), dir = e.direction ?? 0;
      const crosses = (g0 < 0 && g1 >= 0 && dir >= 0) || (g0 > 0 && g1 <= 0 && dir <= 0);
      if (!crosses) return;
      const te = g1 === 0 ? s.t : brent(tt => e.g(tt, s.sampleAt(tt, probe)), prev.t, s.t, 1e-14 * Math.max(1, Math.abs(s.t)));
      const hit = { index, t: te, y: s.sampleAt(te) };
      if (e.terminal && (!firstHit || te < firstHit.t)) firstHit = hit;
      s.events.push(hit);
    });
    if (firstHit) {
      y.set(firstHit.y); s.t = firstHit.t; s.stopped = true;
      s.events = s.events.filter(ev => ev.t <= firstHit.t);
    }
  }
  return s;
}

/**
 * Decide whether y′ = f(t, y) on [t0, tEnd] is stiff and return a suitable stepper.
 * Samples the Jacobian spectrum at t0, at optional extra `probes` [{ t, y }], and at the end of a
 * short explicit Dormand–Prince trial run (`trialSteps`, default 200) — stiffness often develops only
 * after a transient (Robertson kinetics is not stiff at t = 0). The stiffness index
 * S = max(−Re λ)·(tEnd − t0) estimates how many steps an explicit method would need for stability alone;
 * the problem is called stiff when S > threshold (default 1e3). This is a heuristic: a problem whose
 * fast modes are excited throughout (not just decaying) needs small steps with any method.
 * Returns { stiff, method: 'radau5' | 'dopri5', stiffnessIndex, fastest, ratio, reason, trialSpan, stepper }; the
 * stepper is a fresh stiffStepper (radau5) or odeStepper (dopri5) from (t0, y0) with rtol, atol.
 */
export function odeAuto({ f, jac, y0, t0 = 0, tEnd, rtol = 1e-6, atol = 1e-9, probes = [], trialSteps = 200, threshold = 1e3, events } = {}) {
  if (typeof f !== 'function') throw new Error('odeAuto: f(t, y) is required');
  if (!y0?.length || !Array.from(y0).every(Number.isFinite)) throw new Error('odeAuto: y0 must be non-empty and finite');
  if (!(Number.isFinite(tEnd) && tEnd > t0)) throw new Error('odeAuto: tEnd > t0 is required');
  const span = tEnd - t0, points = [{ t: t0, y: Array.from(y0) }, ...probes];
  let trialSpan = 0;
  try {
    const trial = odeStepper({ f, y0, t0, h: Math.min(span, 1e-6 * Math.max(1, span)), rtol, atol });
    for (let i = 0; i < trialSteps && trial.t < tEnd; i++) trial.step(tEnd - trial.t);
    trialSpan = trial.t - t0;
    points.push({ t: trial.t, y: Array.from(trial.y) });
  } catch { /* explicit trial failed (e.g. overflow from instability): judge by the spectrum alone */ }
  let fastest = 0, ratio = 1;
  for (const p of points) {
    const r = stiffnessRatio({ f, jac, t: p.t, y: p.y });
    if (r.fastest > fastest) { fastest = r.fastest; ratio = r.ratio; }
  }
  const stiffnessIndex = fastest * span;
  const stiff = stiffnessIndex > threshold;
  const reason = stiff
    ? `fastest decay rate ${fastest.toPrecision(3)} over a span ${span.toPrecision(3)} (index ${stiffnessIndex.toPrecision(3)} > ${threshold}): explicit steps would be stability-limited`
    : `stiffness index ${stiffnessIndex.toPrecision(3)} ≤ ${threshold}: an explicit method is efficient`;
  const stepper = stiff
    ? stiffStepper({ f, jac, y0, t0, method: 'radau5', rtol, atol, events })
    : odeStepper({ f, y0, t0, h: Math.min(span / 100, 1e-3 * Math.max(1, span)), rtol, atol, events });
  return { stiff, method: stiff ? 'radau5' : 'dopri5', stiffnessIndex, fastest, ratio, reason, trialSpan, stepper };
}
