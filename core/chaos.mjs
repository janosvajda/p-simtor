// Chaos and nonlinear dynamics, independent of subject: classic maps and flows (with Jacobians),
// bifurcation diagrams and Feigenbaum's δ, Lyapunov exponents and spectra (QR / Benettin),
// Poincaré sections, fractal dimensions, fixed points, periodic orbits and their stability.
//
// A MAP is { name, dim, params, map(x) → Float64Array, jacobian(x) → rows }.
// A FLOW is { name, dim, params, f(t, y) → Float64Array, jacobian(t, y) → rows, energy?(y) }:
// f plugs straight into odeStepper (core/ode.mjs). Matrices are arrays of rows (core/linalg.mjs).
// Everything is deterministic; results that are finite-time or finite-sample estimates say so,
// and should be shown with how they settle (history, fit range), not as a single number.

import { odeStepper } from './ode.mjs';
import { brent, linearFit } from './numerics.mjs';
import { newtonSolve } from './linalg.mjs';
import { eigenNonsymmetric } from './spectral.mjs';

const finite = (v, n, who) => {
  if (!v || v.length !== n || !Array.from(v).every(Number.isFinite)) throw new Error(`${who}: expected ${n} finite values`);
};
const num = (x, name, who) => { if (!Number.isFinite(x)) throw new Error(`${who}: ${name} must be finite`); return x; };
const isMap = (m, who) => { if (!m || typeof m.map !== 'function' || !(m.dim > 0)) throw new Error(`${who}: expected a map { dim, map, jacobian }`); };
const isFlow = (s, who) => { if (!s || typeof s.f !== 'function' || !(s.dim > 0)) throw new Error(`${who}: expected a flow { dim, f, jacobian }`); };
const TWO_PI = 2 * Math.PI;

// ---------------------------------------------------------------- maps

/** Logistic map x → r x (1 − x) on [0, 1] (0 ≤ r ≤ 4). r = 4: λ = ln 2; period doubling accumulates at r∞ = 3.5699456… */
export function logisticMap(r) {
  num(r, 'r', 'logisticMap');
  return { name: 'logistic', dim: 1, params: { r }, map: x => Float64Array.of(r * x[0] * (1 - x[0])), jacobian: x => [Float64Array.of(r * (1 - 2 * x[0]))] };
}

/**
 * Tent map x → μ min(x, 1 − x) (0 < μ ≤ 2). |f′| = μ everywhere, so λ = ln μ exactly. At μ = 2 every
 * binary floating-point orbit collapses to 0 within ~53 steps (a finite-precision artefact): use μ < 2.
 */
export function tentMap(mu) {
  num(mu, 'mu', 'tentMap');
  return { name: 'tent', dim: 1, params: { mu }, map: x => Float64Array.of(mu * Math.min(x[0], 1 - x[0])), jacobian: x => [Float64Array.of(x[0] < 0.5 ? mu : -mu)] };
}

/** Hénon map (x, y) → (1 − a x² + y, b x); a = 1.4, b = 0.3 gives the classic strange attractor. det J = −b. */
export function henonMap({ a = 1.4, b = 0.3 } = {}) {
  num(a, 'a', 'henonMap'); num(b, 'b', 'henonMap');
  return {
    name: 'henon', dim: 2, params: { a, b },
    map: ([x, y]) => Float64Array.of(1 - a * x * x + y, b * x),
    jacobian: ([x]) => [Float64Array.of(-2 * a * x, 1), Float64Array.of(b, 0)],
  };
}

/**
 * Chirikov standard map, state [θ, p]: p′ = p + K sin θ, θ′ = θ + p′. Area-preserving (det J = 1).
 * wrap: true reduces θ to [0, 2π) and p to [−π, π) (the torus, for plots); keep wrap: false (the lift)
 * for periodic-orbit searches and diffusion. Last KAM torus breaks at K ≈ 0.9716 (Greene 1979).
 */
export function standardMap(K, { wrap = true } = {}) {
  num(K, 'K', 'standardMap');
  const mod = (v, lo) => v - TWO_PI * Math.floor((v - lo) / TWO_PI);
  return {
    name: 'standard', dim: 2, params: { K, wrap },
    map: ([th, p]) => { const pn = p + K * Math.sin(th), tn = th + pn; return wrap ? Float64Array.of(mod(tn, 0), mod(pn, -Math.PI)) : Float64Array.of(tn, pn); },
    jacobian: ([th]) => { const k = K * Math.cos(th); return [Float64Array.of(1 + k, 1), Float64Array.of(k, 1)]; },
  };
}

/**
 * Sine circle map on the lift: θ → θ + Ω − (K/2π) sin 2πθ (take θ mod 1 to draw it). K < 1 is
 * invertible with Arnold tongues (mode locking); rotationNumber gives the devil's staircase.
 */
export function circleMap(Omega, K) {
  num(Omega, 'Omega', 'circleMap'); num(K, 'K', 'circleMap');
  return {
    name: 'circle', dim: 1, params: { Omega, K },
    map: x => Float64Array.of(x[0] + Omega - K / TWO_PI * Math.sin(TWO_PI * x[0])),
    jacobian: x => [Float64Array.of(1 - K * Math.cos(TWO_PI * x[0]))],
  };
}

/** n iterates after a discarded transient: Float64Array of n·dim values (point k at [k·dim, (k+1)·dim)). Throws if the orbit diverges. */
export function iterateMap(m, x0, n, { transient = 0 } = {}) {
  isMap(m, 'iterateMap'); finite(x0, m.dim, 'iterateMap');
  if (!(Number.isInteger(n) && n >= 0 && Number.isInteger(transient) && transient >= 0)) throw new Error('iterateMap: n and transient must be non-negative integers');
  let x = Float64Array.from(x0);
  for (let i = 0; i < transient; i++) x = m.map(x);
  const out = new Float64Array(n * m.dim);
  for (let k = 0; k < n; k++) {
    x = m.map(x);
    if (!x.every(Number.isFinite)) throw new Error(`iterateMap: orbit diverged at iterate ${k}`);
    out.set(x, k * m.dim);
  }
  return out;
}

/**
 * Bifurcation diagram: for each parameter value, `samples` iterates of one coordinate after a
 * `transient`, starting from x0 (or the previous parameter's last state with { follow: true },
 * which tracks one attractor through hysteresis). family(param) returns a map.
 * Returns [{ param, values: Float64Array }]. Diverging orbits give an empty values array.
 */
export function bifurcationDiagram({ family, params, x0, transient = 500, samples = 100, coordinate = 0, follow = false }) {
  if (typeof family !== 'function' || !params?.length) throw new Error('bifurcationDiagram: need family(param) → map and a list of params');
  let start = x0;
  return Array.from(params, param => {
    const m = family(param);
    let x = Float64Array.from(start), ok = true;
    const values = new Float64Array(samples);
    for (let i = 0; i < transient + samples && ok; i++) {
      x = m.map(x);
      if (!x.every(Number.isFinite)) ok = false;
      else if (i >= transient) values[i - transient] = x[coordinate];
    }
    if (follow && ok) start = x;
    return { param, values: ok ? values : new Float64Array(0) };
  });
}

/**
 * Superstable parameters R_0 < R_1 < … of the period-doubling cascade of a unimodal family
 * x → g(r, x) with critical point xc: R_n has xc on a cycle of period 2ⁿ, i.e. g^{2ⁿ}(R_n, xc) = xc.
 * r0, r1 are the first two (logistic: 2 and 1 + √5); later ones are bracketed from the geometric
 * extrapolation R_n + (R_n − R_{n−1})/δ_est and solved by Brent. Precision limits count to ~14.
 * They converge geometrically to r∞ with ratio Feigenbaum's δ.
 */
export function superstableParameters({ g, xc, r0, r1, count = 10 }) {
  if (typeof g !== 'function') throw new Error('superstableParameters: g(r, x) is required');
  num(xc, 'xc', 'superstableParameters'); num(r0, 'r0', 'superstableParameters'); num(r1, 'r1', 'superstableParameters');
  if (!(Number.isInteger(count) && count >= 2 && count <= 16)) throw new Error('superstableParameters: count must be an integer in [2, 16]');
  const R = [r0, r1];
  let delta = 4.669;
  for (let n = 2; n < count; n++) {
    const period = 2 ** n, d = (R[n - 1] - R[n - 2]) / delta;
    const h = r => { let x = xc; for (let i = 0; i < period; i++) x = g(r, x); return x - xc; };
    const lo = R[n - 1] + 0.4 * d, hi = R[n - 1] + 1.6 * d, steps = 48;
    let a = lo, fa = h(a), root = null;
    for (let k = 1; k <= steps && root === null; k++) {
      const b = lo + (hi - lo) * k / steps, fb = h(b);
      if (fa === 0) root = a;
      else if (fa * fb < 0) root = brent(h, a, b, 1e-16);
      a = b; fa = fb;
    }
    if (root === null) throw new Error(`superstableParameters: could not bracket R_${n}; check g, xc, r0, r1`);
    R.push(root);
    delta = (R[n - 1] - R[n - 2]) / (R[n] - R[n - 1]);
  }
  return Float64Array.from(R);
}

/**
 * Feigenbaum ratios δ_n = (r_n − r_{n−1}) / (r_{n+1} − r_n) from successive period-doubling (or
 * superstable) parameters; they converge to δ = 4.669201609102990… (Feigenbaum, J. Stat. Phys. 19,
 * 25 (1978)) for any smooth unimodal map with a quadratic maximum. Returns Float64Array(length − 2).
 */
export function feigenbaumDelta(points) {
  if (!points || points.length < 3) throw new Error('feigenbaumDelta: need at least three parameter values');
  return Float64Array.from({ length: points.length - 2 }, (_, i) => (points[i + 1] - points[i]) / (points[i + 2] - points[i + 1]));
}

/** Rotation number lim (θ_n − θ_0)/n of a circle-map lift, estimated from n iterates after a transient. */
export function rotationNumber(m, x0, n = 10000, { transient = 1000 } = {}) {
  isMap(m, 'rotationNumber');
  if (m.dim !== 1) throw new Error('rotationNumber: needs a one-dimensional lift');
  let x = Float64Array.from(x0);
  for (let i = 0; i < transient; i++) x = m.map(x);
  const start = x[0];
  for (let i = 0; i < n; i++) x = m.map(x);
  return (x[0] - start) / n;
}

/** Modified Gram–Schmidt on the columns of Q (rows-of-columns form: Q[j] is column j); returns diag(R). */
function orthonormalise(cols) {
  const k = cols.length, r = new Float64Array(k);
  for (let j = 0; j < k; j++) {
    const v = cols[j];
    for (let i = 0; i < j; i++) {
      let d = 0;
      for (let c = 0; c < v.length; c++) d += cols[i][c] * v[c];
      for (let c = 0; c < v.length; c++) v[c] -= d * cols[i][c];
    }
    let n = 0;
    for (let c = 0; c < v.length; c++) n = Math.hypot(n, v[c]);
    if (!(n > 0 && Number.isFinite(n))) throw new Error('Lyapunov spectrum: tangent vectors collapsed or overflowed; shorten the interval');
    for (let c = 0; c < v.length; c++) v[c] /= n;
    r[j] = n;
  }
  return r;
}

/**
 * Lyapunov spectrum of a map (all dim exponents, per iterate) by tangent-map QR (Benettin et al.,
 * Meccanica 15, 9 (1980)): Q ← J(x) Q, re-orthonormalised every step; λ_i = ⟨ln R_ii⟩.
 * 1-D maps reduce to ⟨ln |f′(x)|⟩. The sum equals ⟨ln |det J|⟩ (Hénon: ln b exactly; area-preserving: 0).
 * A finite-orbit estimate: report how history settles. Returns { exponents, sum, history: [{ n, exponents }] }.
 */
export function mapLyapunov(m, x0, { n = 100000, transient = 1000, record = 10 } = {}) {
  isMap(m, 'mapLyapunov'); finite(x0, m.dim, 'mapLyapunov');
  if (typeof m.jacobian !== 'function') throw new Error('mapLyapunov: the map needs a jacobian');
  const d = m.dim;
  let x = Float64Array.from(x0);
  for (let i = 0; i < transient; i++) x = m.map(x);
  let Q = Array.from({ length: d }, (_, j) => Float64Array.from({ length: d }, (_c, i) => (i === j ? 1 : 0)));
  const sums = new Float64Array(d), history = [], every = Math.max(1, Math.floor(n / record));
  for (let k = 1; k <= n; k++) {
    const J = m.jacobian(x);
    Q = Q.map(col => Float64Array.from(J, row => { let s = 0; for (let c = 0; c < d; c++) s += row[c] * col[c]; return s; }));
    const r = orthonormalise(Q);
    for (let i = 0; i < d; i++) sums[i] += Math.log(r[i]);
    x = m.map(x);
    if (!x.every(Number.isFinite)) throw new Error(`mapLyapunov: orbit diverged at iterate ${k}`);
    if (k % every === 0) history.push({ n: k, exponents: sums.map(s => s / k) });
  }
  const exponents = sums.map(s => s / n);
  return { exponents, sum: exponents.reduce((a, b) => a + b, 0), history };
}

// ---------------------------------------------------------------- flows

/** Lorenz system (σ, ρ, β) = (10, 28, 8/3) by default. Divergence −(σ + 1 + β) everywhere, so Σλ = −(σ+1+β) exactly. */
export function lorenz({ sigma = 10, rho = 28, beta = 8 / 3 } = {}) {
  [sigma, rho, beta].forEach(v => num(v, 'parameter', 'lorenz'));
  return {
    name: 'lorenz', dim: 3, params: { sigma, rho, beta },
    f: (_t, [x, y, z]) => Float64Array.of(sigma * (y - x), x * (rho - z) - y, x * y - beta * z),
    jacobian: (_t, [x, y, z]) => [Float64Array.of(-sigma, sigma, 0), Float64Array.of(rho - z, -1, -x), Float64Array.of(y, x, -beta)],
  };
}

/** Rössler system (a, b, c) = (0.2, 0.2, 5.7) by default: x′ = −y − z, y′ = x + a y, z′ = b + z (x − c). */
export function rossler({ a = 0.2, b = 0.2, c = 5.7 } = {}) {
  [a, b, c].forEach(v => num(v, 'parameter', 'rossler'));
  return {
    name: 'rossler', dim: 3, params: { a, b, c },
    f: (_t, [x, y, z]) => Float64Array.of(-y - z, x + a * y, b + z * (x - c)),
    jacobian: (_t, [x, , z]) => [Float64Array.of(0, -1, -1), Float64Array.of(1, a, 0), Float64Array.of(z, 0, x - c)],
  };
}

/**
 * Forced Duffing oscillator, state [x, v]: x″ + δ x′ + α x + β x³ = γ cos ωt. Defaults (δ, α, β, γ, ω) =
 * (0.3, −1, 1, 0.5, 1.2): the chaotic two-well oscillator (Guckenheimer & Holmes 1983, §2.2).
 * Non-autonomous: strobe it with a Poincaré section at t mod 2π/ω (section: (t) => sin(ω t / 2)).
 */
export function duffing({ delta = 0.3, alpha = -1, beta = 1, gamma = 0.5, omega = 1.2 } = {}) {
  [delta, alpha, beta, gamma, omega].forEach(v => num(v, 'parameter', 'duffing'));
  return {
    name: 'duffing', dim: 2, params: { delta, alpha, beta, gamma, omega },
    f: (t, [x, v]) => Float64Array.of(v, -delta * v - alpha * x - beta * x * x * x + gamma * Math.cos(omega * t)),
    jacobian: (_t, [x]) => [Float64Array.of(0, 1), Float64Array.of(-alpha - 3 * beta * x * x, -delta)],
    energy: ([x, v]) => 0.5 * v * v + 0.5 * alpha * x * x + 0.25 * beta * x ** 4,
  };
}

/**
 * Hénon–Heiles Hamiltonian, state [x, y, px, py]: H = ½(px² + py² + x² + y²) + x²y − y³/3.
 * Bounded motion for H < 1/6; regular at H = 1/12, largely chaotic at H = 1/6 (Hénon & Heiles, AJ 69, 73 (1964)).
 */
export function henonHeiles() {
  return {
    name: 'henon-heiles', dim: 4, params: {},
    f: (_t, [x, y, px, py]) => Float64Array.of(px, py, -x - 2 * x * y, -y - x * x + y * y),
    jacobian: (_t, [x, y]) => [Float64Array.of(0, 0, 1, 0), Float64Array.of(0, 0, 0, 1), Float64Array.of(-1 - 2 * y, -2 * x, 0, 0), Float64Array.of(-2 * x, -1 + 2 * y, 0, 0)],
    energy: ([x, y, px, py]) => 0.5 * (px * px + py * py + x * x + y * y) + x * x * y - y * y * y / 3,
  };
}

/**
 * Planar double pendulum (point masses m1, m2 on massless rods l1, l2, gravity g), state
 * [θ1, θ2, ω1, ω2], angles from the downward vertical. Exact Lagrangian equations; energy(y) is the
 * conserved total. The Jacobian is by central differences (relative step 1e-6, ~1e-9 accurate).
 */
export function doublePendulum({ m1 = 1, m2 = 1, l1 = 1, l2 = 1, g = 9.81 } = {}) {
  for (const v of [m1, m2, l1, l2]) if (!(v > 0)) throw new Error('doublePendulum: masses and lengths must be positive');
  num(g, 'g', 'doublePendulum');
  const f = (_t, [a, b, w1, w2]) => {
    const d = a - b, den = 2 * m1 + m2 - m2 * Math.cos(2 * d);
    return Float64Array.of(w1, w2,
      (-g * (2 * m1 + m2) * Math.sin(a) - m2 * g * Math.sin(a - 2 * b) - 2 * Math.sin(d) * m2 * (w2 * w2 * l2 + w1 * w1 * l1 * Math.cos(d))) / (l1 * den),
      (2 * Math.sin(d) * (w1 * w1 * l1 * (m1 + m2) + g * (m1 + m2) * Math.cos(a) + w2 * w2 * l2 * m2 * Math.cos(d))) / (l2 * den));
  };
  return {
    name: 'double-pendulum', dim: 4, params: { m1, m2, l1, l2, g }, f,
    jacobian: (t, y) => {
      const cols = [0, 1, 2, 3].map(j => {
        const e = 1e-6 * Math.max(1, Math.abs(y[j])), p = Float64Array.from(y), q = Float64Array.from(y);
        p[j] += e; q[j] -= e;
        const fp = f(t, p), fq = f(t, q);
        return fp.map((v, i) => (v - fq[i]) / (2 * e));
      });
      return [0, 1, 2, 3].map(i => Float64Array.from(cols, c => c[i]));
    },
    energy: ([a, b, w1, w2]) => 0.5 * (m1 + m2) * l1 * l1 * w1 * w1 + 0.5 * m2 * l2 * l2 * w2 * w2 + m2 * l1 * l2 * w1 * w2 * Math.cos(a - b)
      - (m1 + m2) * g * l1 * Math.cos(a) - m2 * g * l2 * Math.cos(b),
  };
}

/**
 * Full Lyapunov spectrum of a flow by Benettin's QR method: the state and dim tangent vectors
 * (Φ′ = J Φ) are integrated together (adaptive Dormand–Prince), re-orthonormalised every `interval`,
 * and λ_i = Σ ln R_ii / T after a `transient`. Requires system.jacobian. Exponents are sorted as
 * produced (descending for a typical run). Σλ equals the time-average divergence (Lorenz: −(σ+1+β)
 * exactly), a built-in check. Finite-time estimate: report the history.
 * Benchmark: Lorenz (10, 28, 8/3) λ = 0.9056, 0, −14.5723 (Sprott, Chaos and Time-Series Analysis, 2003).
 * Returns { exponents, sum, history: [{ t, exponents }] }.
 */
export function lyapunovSpectrum({ system, y0, t0 = 0, duration, interval = 1, transient = 0, h = 1e-3, rtol = 1e-10, atol = 1e-12 }) {
  isFlow(system, 'lyapunovSpectrum'); finite(y0, system.dim, 'lyapunovSpectrum');
  if (typeof system.jacobian !== 'function') throw new Error('lyapunovSpectrum: the flow needs a jacobian');
  if (!(duration > 0 && interval > 0 && transient >= 0)) throw new Error('lyapunovSpectrum: need duration, interval > 0 and transient ≥ 0');
  const n = system.dim, out = new Float64Array(n + n * n);
  const aug = (t, z) => {
    const y = z.subarray(0, n), J = system.jacobian(t, y);
    out.set(system.f(t, y));
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) { let s = 0; for (let c = 0; c < n; c++) s += J[i][c] * z[n + j * n + c]; out[n + j * n + i] = s; }
    }
    return out;
  };
  const z = new Float64Array(n + n * n);
  z.set(y0);
  for (let j = 0; j < n; j++) z[n + j * n + j] = 1;
  const sums = new Float64Array(n), history = [], start = t0 + transient, end = start + duration;
  let t = t0, step = h;
  while (t < end - 1e-9 * interval) {
    const stop = t < start ? Math.min(t + interval, start) : Math.min(t + interval, end);
    const s = odeStepper({ f: aug, y0: z, t0: t, method: 'dopri5', h: step, rtol, atol });
    s.advanceTo(stop, Infinity, { exact: true });
    step = Math.max((s.t - t) / Math.max(s.steps, 1), 1e-12);
    z.set(s.y);
    const cols = Array.from({ length: n }, (_, j) => z.subarray(n + j * n, n + (j + 1) * n));
    const r = orthonormalise(cols);
    if (t >= start - 1e-12 * interval) {
      for (let i = 0; i < n; i++) sums[i] += Math.log(r[i]);
      history.push({ t: stop, exponents: sums.map(v => v / (stop - start)) });
    }
    t = stop;
  }
  const exponents = sums.map(v => v / duration);
  return { exponents, sum: exponents.reduce((a, b) => a + b, 0), history };
}

/**
 * Poincaré section of a flow: the states where section(t, y) crosses zero in `direction` (+1 rising,
 * −1 falling, 0 both), located on odeStepper's dense output with Brent's method (event accuracy
 * O(h⁴) of the local step: keep rtol tight). Stops after `count` crossings or at tmax.
 * Example: Hénon–Heiles, section (t, y) => y[0] with direction +1 gives the (y, py) plane at x = 0, px > 0.
 * Returns { points: Float64Array[] (full states), times: Float64Array, stepper }.
 */
export function poincareSection({ system, y0, t0 = 0, section, direction = 1, count = 200, tmax = Infinity, h = 1e-2, rtol = 1e-11, atol = 1e-13 }) {
  isFlow(system, 'poincareSection'); finite(y0, system.dim, 'poincareSection');
  if (typeof section !== 'function') throw new Error('poincareSection: section(t, y) is required');
  if (!(Number.isInteger(count) && count > 0) || !(tmax > t0)) throw new Error('poincareSection: need a positive integer count and tmax > t0');
  const s = odeStepper({ f: system.f, y0, t0, h, rtol, atol, events: [{ g: section, direction }] });
  while (s.events.length < count && s.t < tmax) s.step();
  const hits = s.events.filter(e => e.t <= tmax).slice(0, count);
  return { points: hits.map(e => e.y), times: Float64Array.from(hits, e => e.t), stepper: s };
}

// ---------------------------------------------------------------- fractal dimensions

const asPoints = (points, who) => {
  if (!Array.isArray(points) || points.length < 2) throw new Error(`${who}: need an array of at least two points`);
  const d = points[0].length;
  for (const p of points) finite(p, d, who);
  return d;
};

/**
 * Correlation dimension D2 (Grassberger & Procaccia, PRL 50, 346 (1983)): C(r) = fraction of pairs
 * closer than r (Euclidean), D2 = slope of ln C against ln r over the scaling range [rMin, rMax]
 * (your judgement: plot `curve` and choose where it is straight). theiler excludes pairs closer
 * than that many samples in time (flows sampled finely). O(N²) pairs: N up to a few thousand.
 * Benchmarks: Hénon 1.21 ± 0.01, Lorenz 2.05 ± 0.01 (same paper). Finite-N estimates are biased
 * low at small r; quote the fit range. Returns { dimension, r2, curve: [{ r, C }] }.
 */
export function correlationDimension(points, { rMin, rMax, bins = 16, theiler = 0 } = {}) {
  asPoints(points, 'correlationDimension');
  if (!(rMin > 0 && rMax > rMin)) throw new Error('correlationDimension: need 0 < rMin < rMax');
  const N = points.length, lo = Math.log(rMin), span = Math.log(rMax) - lo;
  const hist = new Float64Array(bins + 1);                   // hist[k]: pairs with r < r_k, r_k = rMin (rMax/rMin)^(k/bins)
  let pairs = 0;
  for (let i = 0; i < N; i++) {
    const a = points[i];
    for (let j = i + 1 + theiler; j < N; j++) {
      const b = points[j];
      let d2 = 0;
      for (let c = 0; c < a.length; c++) { const e = a[c] - b[c]; d2 += e * e; }
      pairs++;
      if (d2 >= rMax * rMax) continue;
      const k = d2 < rMin * rMin ? 0 : Math.floor((0.5 * Math.log(d2) - lo) / span * bins) + 1;
      if (k <= bins) hist[k]++;
    }
  }
  if (!pairs) throw new Error('correlationDimension: no pairs left (theiler too large)');
  const curve = [];
  let cum = 0;
  for (let k = 0; k <= bins; k++) { cum += hist[k]; curve.push({ r: Math.exp(lo + span * k / bins), C: cum / pairs }); }
  const used = curve.filter(p => p.C > 0);
  if (used.length < 3) throw new Error('correlationDimension: too few pairs inside [rMin, rMax]; widen the range or add points');
  const fit = linearFit(used.map(p => Math.log(p.r)), used.map(p => Math.log(p.C)));
  return { dimension: fit.b, r2: fit.r2, curve };
}

/**
 * Box-counting (capacity) dimension: N(ε) = number of occupied grid cells of side ε, D0 = slope of
 * ln N against ln(1/ε) over the given sizes. Grid anchored at the origin; finite samples undercount
 * small boxes, so keep ε well above the mean point spacing. Returns { dimension, r2, counts }.
 * Benchmarks: middle-thirds Cantor set ln 2/ln 3; Hénon attractor 1.261 ± 0.003 (Russell, Hanson & Ott, PRL 45, 1175 (1980)).
 */
export function boxCountingDimension(points, { sizes }) {
  asPoints(points, 'boxCountingDimension');
  if (!sizes || sizes.length < 2 || !Array.from(sizes).every(e => e > 0)) throw new Error('boxCountingDimension: need at least two positive box sizes');
  const counts = Array.from(sizes, eps => {
    const cells = new Set();
    for (const p of points) cells.add(Array.from(p, v => Math.floor(v / eps)).join(','));
    return cells.size;
  });
  const fit = linearFit(Array.from(sizes, e => -Math.log(e)), counts.map(c => Math.log(c)));
  return { dimension: fit.b, r2: fit.r2, counts: Float64Array.from(counts) };
}

// ---------------------------------------------------------------- fixed points, periodic orbits, stability

/**
 * Eigenvalues of a small real square matrix (n ≤ 8) as [{ re, im }], sorted by descending real part.
 * Uses the balanced Hessenberg/Francis QR solver, without forming a characteristic polynomial
 * (which can destroy the accuracy of clustered and repeated roots). Backward error is at floating-
 * point precision for well-scaled matrices; forward accuracy depends on eigenvalue conditioning.
 * For symmetric matrices use eigenSymmetric (core/linalg.mjs).
 */
export function smallEigenvalues(J) {
  const n = J?.length;
  if (!(n >= 1 && n <= 8)) throw new Error('smallEigenvalues: matrix size must be 1…8');
  for (const row of J) finite(row, n, 'smallEigenvalues');
  return eigenNonsymmetric(J).sort((a, b) => b.re - a.re || b.im - a.im);
}

/**
 * Linear stability of a fixed point from its Jacobian. kind 'flow' uses Re λ (stable < 0), kind 'map'
 * uses |λ| (stable < 1); |Re λ| or ||λ| − 1| ≤ tol counts as marginal (non-hyperbolic: linearisation
 * inconclusive). Types: 'stable node' | 'stable focus' | 'unstable node' | 'unstable focus' | 'saddle' |
 * 'saddle-focus' | 'center' (flow, purely imaginary) | 'elliptic' (map, on the unit circle) | 'non-hyperbolic'.
 * Returns { eigenvalues, type, stable, hyperbolic, unstableDimension }.
 */
export function classifyFixedPoint(J, { kind = 'flow', tol = 1e-9 } = {}) {
  if (kind !== 'flow' && kind !== 'map') throw new Error("classifyFixedPoint: kind must be 'flow' or 'map'");
  const eigenvalues = smallEigenvalues(J);
  const growth = e => (kind === 'flow' ? e.re : Math.hypot(e.re, e.im) - 1);
  const g = eigenvalues.map(growth), complex = eigenvalues.some(e => e.im !== 0);
  const up = g.filter(v => v > tol).length, down = g.filter(v => v < -tol).length, marginal = g.length - up - down;
  let type;
  if (marginal) type = marginal === g.length && complex ? (kind === 'flow' ? 'center' : 'elliptic') : 'non-hyperbolic';
  else if (!up) type = complex ? 'stable focus' : 'stable node';
  else if (!down) type = complex ? 'unstable focus' : 'unstable node';
  else type = complex ? 'saddle-focus' : 'saddle';
  return { eigenvalues, type, stable: up === 0 && marginal === 0, hyperbolic: marginal === 0, unstableDimension: up };
}

/**
 * Fixed point of a flow (f(t, y) = 0 at time t) by damped Newton (newtonSolve, analytic Jacobian),
 * with its classification. Returns { x, converged, residualNorm, ...classifyFixedPoint }.
 * Newton finds the fixed point nearest the guess's basin: seed it from the formula when one exists.
 */
export function findFixedPoint(system, guess, { t = 0, tol = 1e-9 } = {}) {
  isFlow(system, 'findFixedPoint'); finite(guess, system.dim, 'findFixedPoint');
  const sol = newtonSolve(y => system.f(t, y), guess, { jacobian: system.jacobian ? y => system.jacobian(t, y) : undefined });
  return { x: sol.x, converged: sol.converged, residualNorm: sol.residualNorm, ...classifyFixedPoint(system.jacobian(t, sol.x), { kind: 'flow', tol }) };
}

/**
 * Periodic orbit of a map: x with F^p(x) = x by damped Newton on F^p(x) − x, Jacobian
 * DF^p = J(x_{p−1})…J(x_0) − 1. Use the unwrapped (lift) form of angle maps. Returns { x, orbit
 * (p points), converged, multipliers (eigenvalues of DF^p), ...classifyFixedPoint(DF^p, 'map') }.
 */
export function findPeriodicOrbit(m, x0, period, { tol = 1e-9 } = {}) {
  isMap(m, 'findPeriodicOrbit'); finite(x0, m.dim, 'findPeriodicOrbit');
  if (!(Number.isInteger(period) && period >= 1)) throw new Error('findPeriodicOrbit: period must be a positive integer');
  const d = m.dim;
  const run = x => {
    let y = Float64Array.from(x), D = Array.from({ length: d }, (_, i) => Float64Array.from({ length: d }, (_c, j) => (i === j ? 1 : 0)));
    const orbit = [];
    for (let k = 0; k < period; k++) {
      orbit.push(y);
      const J = m.jacobian(y);
      D = J.map(row => Float64Array.from({ length: d }, (_, j) => { let s = 0; for (let c = 0; c < d; c++) s += row[c] * D[c][j]; return s; }));
      y = m.map(y);
    }
    return { y, D, orbit };
  };
  const sol = newtonSolve(x => { const { y } = run(x); return y.map((v, i) => v - x[i]); }, x0, {
    jacobian: x => run(x).D.map((row, i) => row.map((v, j) => v - (i === j ? 1 : 0))),
  });
  const { D, orbit } = run(sol.x), cls = classifyFixedPoint(D, { kind: 'map', tol });
  return { x: sol.x, orbit, converged: sol.converged, residualNorm: sol.residualNorm, multipliers: cls.eigenvalues, ...cls };
}
