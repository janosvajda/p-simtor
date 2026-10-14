// General time integration for finite-dimensional evolution, independent of subject.
//
//   odeStepper       y' = f(t, y)            any first-order system (RK4 or adaptive Dormand–Prince)
//   mechanicsStepper q'' = a(t, q)           position-only accelerations (symplectic leapfrog / Yoshida-4)
//
// Both keep the previous step so the state can be evaluated at any time between steps
// (cubic Hermite), which makes playback independent of frame rate: the step sequence
// depends only on the problem and the step parameters, never on when frames are drawn.

import { DORMAND_PRINCE as DP } from './numerics.mjs';

/** Cubic Hermite interpolation of values y with derivatives d between (t0, t1). */
function hermite(t0, y0, d0, t1, y1, d1, t, out) {
  const h = t1 - t0;
  if (!(h !== 0)) { out.set(y1); return out; }
  const u = Math.min(1, Math.max(0, (t - t0) / h)), u2 = u * u, u3 = u2 * u;
  const a = 2 * u3 - 3 * u2 + 1, b = u3 - 2 * u2 + u, c = -2 * u3 + 3 * u2, d = u3 - u2;
  for (let i = 0; i < y0.length; i++) out[i] = a * y0[i] + b * h * d0[i] + c * y1[i] + d * h * d1[i];
  return out;
}

const finiteArray = (a, who) => {
  if (!a.length || !Array.from(a).every(Number.isFinite)) throw new Error(`${who}: state must be non-empty and finite`);
};

/**
 * Stepper for y' = f(t, y).
 *   method 'rk4'    fixed step h (classical 4th order)
 *   method 'dopri5' adaptive Dormand–Prince 5(4) with error control (rtol, atol); h is the first trial step
 * f(t, y) returns the derivative (array-like of the same length).
 * Returns { t, y, step(), advanceTo(t, maxSteps), sampleAt(t, out), steps, rejected }.
 */
export function odeStepper({ f, y0, t0 = 0, method = 'dopri5', h, rtol = 1e-9, atol = 1e-12, hmax = Infinity }) {
  if (typeof f !== 'function') throw new Error('odeStepper: f(t, y) is required');
  finiteArray(y0, 'odeStepper');
  if (!(Number.isFinite(h) && h > 0)) throw new Error('odeStepper: a positive finite step h is required');
  if (method !== 'rk4' && method !== 'dopri5') throw new Error(`odeStepper: unknown method "${method}"`);
  const n = y0.length;
  const y = Float64Array.from(y0), dy = Float64Array.from(f(t0, y));
  const prev = { t: t0, y: y.slice(), dy: dy.slice() };
  const k = Array.from({ length: 7 }, () => new Float64Array(n));
  const tmp = new Float64Array(n), yn = new Float64Array(n);
  let hNext = h;
  const s = {
    t: t0, y, dy, steps: 0, rejected: 0, method,
    step() {
      prev.t = s.t; prev.y.set(y); prev.dy.set(dy);
      if (method === 'rk4') {
        const k1 = dy, k2 = new Float64Array(n), k3 = new Float64Array(n), k4 = new Float64Array(n);
        for (let i = 0; i < n; i++) tmp[i] = y[i] + h / 2 * k1[i];
        k2.set(f(s.t + h / 2, tmp));
        for (let i = 0; i < n; i++) tmp[i] = y[i] + h / 2 * k2[i];
        k3.set(f(s.t + h / 2, tmp));
        for (let i = 0; i < n; i++) tmp[i] = y[i] + h * k3[i];
        k4.set(f(s.t + h, tmp));
        for (let i = 0; i < n; i++) y[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        s.t += h;
      } else {
        for (let attempt = 0; ; attempt++) {
          if (attempt > 60) throw new Error(`odeStepper: cannot meet tolerance at t = ${s.t}`);
          const hh = Math.min(hNext, hmax);
          if (s.t + hh === s.t) throw new Error(`odeStepper: step underflow at t = ${s.t}`);
          k[0].set(dy);
          for (let st = 1; st < 7; st++) {
            for (let i = 0; i < n; i++) {
              let acc = y[i];
              for (let j = 0; j < st; j++) acc += hh * DP.a[st][j] * k[j][i];
              tmp[i] = acc;
            }
            k[st].set(f(s.t + DP.c[st] * hh, tmp));
          }
          let err = 0;
          for (let i = 0; i < n; i++) {
            let acc = y[i], e = 0;
            for (let j = 0; j < 7; j++) { acc += hh * DP.b[j] * k[j][i]; e += hh * DP.e[j] * k[j][i]; }
            yn[i] = acc;
            err = Math.max(err, Math.abs(e) / (atol + rtol * Math.max(Math.abs(y[i]), Math.abs(acc))));
          }
          if (!Number.isFinite(err)) throw new Error(`odeStepper: non-finite derivative or state at t = ${s.t}`);
          const fac = err === 0 ? 5 : Math.min(5, Math.max(0.2, 0.9 * Math.pow(err, -0.2)));
          hNext = hh * fac;
          if (err <= 1) { y.set(yn); s.t += hh; break; }
          s.rejected++;
        }
      }
      dy.set(f(s.t, y));
      if (!Array.from(y).every(Number.isFinite)) throw new Error(`odeStepper: state became non-finite at t = ${s.t}`);
      s.steps++;
    },
    /** Advance until t ≥ target. Returns false if maxSteps were used first. */
    advanceTo(target, maxSteps = Infinity) {
      let i = 0;
      while (s.t < target) { if (i++ >= maxSteps) return false; s.step(); }
      return true;
    },
    /** State at time t within the last step (cubic Hermite with the derivatives at both ends). */
    sampleAt(t, out = new Float64Array(n)) { return hermite(prev.t, prev.y, prev.dy, s.t, y, dy, t, out); },
  };
  return s;
}

const CBRT2 = Math.cbrt(2);
const YOSHIDA4 = [1 / (2 - CBRT2), -CBRT2 / (2 - CBRT2), 1 / (2 - CBRT2)];

/**
 * Symplectic stepper for q'' = a(t, q) with position-only accelerations (conservative,
 * separable mechanics). accel(t, q, out) writes the acceleration into out.
 *   method 'leapfrog' (kick–drift–kick, 2nd order) or 'yoshida4' (4th order composition)
 * Bounded energy error over long runs; time-reversible. Fixed step h.
 * Returns { t, q, v, step(), advanceTo(t, maxSteps), sampleAt(t, outQ, outV), steps }.
 */
export function mechanicsStepper({ accel, q0, v0, t0 = 0, h, method = 'yoshida4' }) {
  if (typeof accel !== 'function') throw new Error('mechanicsStepper: accel(t, q, out) is required');
  finiteArray(q0, 'mechanicsStepper'); finiteArray(v0, 'mechanicsStepper');
  if (q0.length !== v0.length) throw new Error('mechanicsStepper: q0 and v0 must have the same length');
  if (!(Number.isFinite(h) && h !== 0)) throw new Error('mechanicsStepper: a non-zero finite step h is required');
  const weights = method === 'leapfrog' ? [1] : method === 'yoshida4' ? YOSHIDA4 : null;
  if (!weights) throw new Error(`mechanicsStepper: unknown method "${method}"`);
  const n = q0.length;
  const q = Float64Array.from(q0), v = Float64Array.from(v0), a = new Float64Array(n);
  const prev = { t: t0, q: q.slice(), v: v.slice() };
  accel(t0, q, a);
  const s = {
    t: t0, q, v, steps: 0, h, method,
    step() {
      prev.t = s.t; prev.q.set(q); prev.v.set(v);
      for (const w of weights) {
        const dt = w * s.h;
        for (let i = 0; i < n; i++) v[i] += dt / 2 * a[i];
        for (let i = 0; i < n; i++) q[i] += dt * v[i];
        s.t += dt;
        accel(s.t, q, a);
        for (let i = 0; i < n; i++) v[i] += dt / 2 * a[i];
      }
      if (!q.every(Number.isFinite) || !v.every(Number.isFinite)) throw new Error(`mechanicsStepper: state became non-finite at t = ${s.t}`);
      s.steps++;
    },
    advanceTo(target, maxSteps = Infinity) {
      let i = 0;
      while ((s.h > 0 ? s.t < target : s.t > target)) { if (i++ >= maxSteps) return false; s.step(); }
      return true;
    },
    /** Positions (and velocities) at time t within the last step: cubic Hermite, O(h⁴). */
    sampleAt(t, outQ = new Float64Array(n), outV) {
      hermite(prev.t, prev.q, prev.v, s.t, q, v, t, outQ);
      if (outV) {
        // velocity: linear blend of end velocities (adequate for display; use step states for analysis)
        const span = s.t - prev.t, u = span ? Math.min(1, Math.max(0, (t - prev.t) / span)) : 1;
        for (let i = 0; i < n; i++) outV[i] = (1 - u) * prev.v[i] + u * v[i];
      }
      return outQ;
    },
  };
  return s;
}

/**
 * Observed order of accuracy of any solver: run(h) returns a result vector at a fixed end
 * time; the order is log2(|r(h) − r(h/2)| / |r(h/2) − r(h/4)|) using the max-norm.
 * Works for ODEs, PDEs and quadrature alike. Returns { order, errors: [e1, e2], results }.
 */
export function refinementOrder(run, h) {
  const r = [h, h / 2, h / 4].map(x => Array.from(run(x)));
  const diff = (a, b) => a.reduce((m, ai, i) => Math.max(m, Math.abs(ai - b[i])), 0);
  const e1 = diff(r[0], r[1]), e2 = diff(r[1], r[2]);
  return { order: e2 > 0 ? Math.log2(e1 / e2) : Infinity, errors: [e1, e2], results: r };
}

/**
 * Largest Lyapunov exponent of dy/dt = f(t, y): the average exponential growth rate of a small
 * perturbation (> 0: sensitive dependence on initial conditions; ≈ 0: regular motion).
 * With jvp(t, y, w) (the Jacobian–vector product J·w) it integrates the linearised equations;
 * otherwise a twin trajectory d0·max(1, |y|) away (Benettin's method). Both share the same adaptive steps and
 * are renormalised about every `interval`. Discard a `transient` before averaging.
 * Returns { lambda, history: [{ t, lambda }] } — a finite-time estimate: report how it settles.
 */
export function lyapunovExponent({ f, jvp, y0, t0 = 0, duration, interval = 1, transient = 0, d0 = 1e-8, h = 1e-3, rtol = 1e-10, atol = 1e-14 }) {
  if (typeof f !== 'function') throw new Error('lyapunovExponent: f(t, y) is required');
  finiteArray(y0, 'lyapunovExponent');
  if (!(duration > 0 && interval > 0 && transient >= 0 && d0 > 0)) throw new Error('lyapunovExponent: need duration, interval, d0 > 0 and transient ≥ 0');
  const n = y0.length, out = new Float64Array(2 * n), dir = new Float64Array(n).fill(1 / Math.sqrt(n));
  const aug = jvp
    ? (t, z) => { const y = z.subarray(0, n); out.set(f(t, y)); out.set(jvp(t, y, z.subarray(n)), n); return out; }
    : (t, z) => { out.set(f(t, z.subarray(0, n))); out.set(f(t, z.subarray(n)), n); return out; };
  const z = new Float64Array(2 * n);
  let eps = d0;                                     // twin offset, scaled with |y| so rounding cannot swamp it
  const restart = (y, u) => {
    z.set(y);
    eps = d0 * Math.max(1, Math.hypot(...y));
    for (let i = 0; i < n; i++) z[n + i] = jvp ? u[i] : y[i] + eps * u[i];
  };
  restart(y0, dir);
  const end = t0 + transient + duration, history = [];
  let t = t0, sum = 0, step = h;
  while (t < end - 1e-12 * interval) {
    const s = odeStepper({ f: aug, y0: z, t0: t, method: 'dopri5', h: step, rtol, atol });
    const stop = Math.min(t + interval, end);
    while (s.t < stop) s.step();
    step = Math.max((s.t - t) / Math.max(s.steps, 1), 1e-12);
    const y = s.y.subarray(0, n), sep = new Float64Array(n);
    for (let i = 0; i < n; i++) sep[i] = jvp ? s.y[n + i] : s.y[n + i] - y[i];
    const len = Math.hypot(...sep), growth = len / (jvp ? 1 : eps);
    if (!(growth > 0 && Number.isFinite(growth))) throw new Error('lyapunovExponent: perturbation collapsed or overflowed; shorten interval');
    if (s.t > t0 + transient + 1e-12 * interval) {
      sum += Math.log(growth) * (s.t - Math.max(t, t0 + transient)) / (s.t - t);
      history.push({ t: s.t, lambda: sum / (s.t - t0 - transient) });
    }
    restart(y, sep.map(v => v / len));
    t = s.t;
  }
  return { lambda: history.at(-1)?.lambda ?? NaN, history };
}

/**
 * Exact solution of first-order linear chains with branching (radioactive decay series, sequential
 * first-order reactions, compartment chains): species i is lost at rate rates[i] (0 = stable), and
 * links [[from, to, fraction]] send that fraction of from's losses to `to` (from < to).
 * Bateman's solution, generalised to branching: eigen-decomposition of the triangular rate matrix.
 * Rates must be distinct (Bateman's formula divides by their differences; for equal or nearly
 * equal rates use expm from core/linalg.mjs). Returns { at(t, x0) → amounts, matrix }.
 */
export function linearChain({ rates, links = [] }) {
  const n = rates.length;
  if (!n || !rates.every(r => Number.isFinite(r) && r >= 0)) throw new Error('linearChain: rates must be finite and ≥ 0');
  const A = Array.from({ length: n }, (_, i) => { const row = new Float64Array(n); row[i] = -rates[i]; return row; });
  const out = new Float64Array(n);
  for (const [from, to, fraction] of links) {
    if (!(Number.isInteger(from) && Number.isInteger(to) && from < to && to < n && fraction >= 0)) throw new Error('linearChain: links are [from, to, fraction] with from < to');
    A[to][from] += fraction * rates[from];
    out[from] += fraction;
  }
  if (out.some(f => f > 1 + 1e-12)) throw new Error('linearChain: branching fractions from one species exceed 1');
  for (let i = 0; i < n; i++) {for (let j = 0; j < i; j++) {
    if (Math.abs(rates[i] - rates[j]) <= 1e-9 * Math.max(rates[i], rates[j])) throw new Error('linearChain: equal rates; use expm from core/linalg.mjs');
  }}
  // eigenvector k (eigenvalue −λ_k) of the lower-triangular matrix, unit at position k
  const V = Array.from({ length: n }, (_, k) => {
    const v = new Float64Array(n); v[k] = 1;
    for (let i = k + 1; i < n; i++) { let s = 0; for (let j = k; j < i; j++) s += A[i][j] * v[j]; v[i] = s / (rates[i] - rates[k]); }
    return v;
  });
  return {
    matrix: A,
    at(t, x0) {
      if (x0.length !== n) throw new Error('linearChain: x0 must have one amount per species');
      const c = new Float64Array(n);                       // V c = x0, V unit lower triangular (columns = eigenvectors)
      for (let i = 0; i < n; i++) { let s = x0[i]; for (let k = 0; k < i; k++) s -= V[k][i] * c[k]; c[i] = s; }
      const x = new Float64Array(n);
      for (let k = 0; k < n; k++) { const e = c[k] * Math.exp(-rates[k] * t); for (let i = k; i < n; i++) x[i] += e * V[k][i]; }
      return x;
    },
  };
}
