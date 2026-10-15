// General numerical toolkit: random numbers, FFT, ODE solvers, quadrature, root
// finding, special functions, linear algebra. Pure functions, browser + Node.
// Every routine is tested against known values in tests/numerics.test.mjs.
// Import these instead of writing them again in a page.

// ---------------------------------------------------------------- random numbers

/** Seeded PRNG (mulberry32), uniform in [0, 1). Deterministic for a given seed. */
export function rng(seed = 1) {
  let a = seed >>> 0;
  const next = () => {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  let spare = null;
  next.normal = () => {               // standard normal (Box–Muller, cached pair)
    if (spare != null) { const s = spare; spare = null; return s; }
    let u = 0; while (u === 0) u = next();
    const v = next(), r = Math.sqrt(-2 * Math.log(u));
    spare = r * Math.sin(2 * Math.PI * v);
    return r * Math.cos(2 * Math.PI * v);
  };
  return next;
}

/**
 * An independent seeded stream for item i (atom, particle, event): rngFor(seed, i) always gives
 * the same numbers for that item, whatever else was drawn, so one item can be re-traced alone.
 * The seed and index are mixed by a 32-bit hash (splitmix-style finaliser) before mulberry32.
 */
export function rngFor(seed, i) {
  let h = (Math.imul((seed >>> 0) ^ 0x9E3779B9, 0x85EBCA6B) + Math.imul(i >>> 0, 0xC2B2AE35)) >>> 0;
  h = Math.imul(h ^ (h >>> 16), 0x7FEB352D) >>> 0;
  h = Math.imul(h ^ (h >>> 15), 0x846CA68B) >>> 0;
  return rng((h ^ (h >>> 16)) >>> 0);
}

// ---------------------------------------------------------------- FFT

const fftTables = new Map();
/** Twiddle factors e^{−2πik/n} (k < n/2) and the bit-reversal permutation, cached per length. */
function fftTable(n) {
  let t = fftTables.get(n);
  if (!t) {
    const cos = new Float64Array(n >> 1), sin = new Float64Array(n >> 1), rev = new Uint32Array(n);
    for (let k = 0; k < n >> 1; k++) { cos[k] = Math.cos(2 * Math.PI * k / n); sin[k] = Math.sin(2 * Math.PI * k / n); }
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit; rev[i] = j;
    }
    t = { cos, sin, rev };
    if (fftTables.size > 32) fftTables.clear();
    fftTables.set(n, t);
  }
  return t;
}

/** In-place radix-2 complex FFT. re/im: Float64Array, length a power of two.
 *  Forward: X_k = Σ x_n e^{−2πikn/N}. Inverse (inverse = true) includes the 1/N factor. */
export function fft(re, im, inverse = false) {
  const n = re.length;
  if (n & (n - 1)) throw new Error('fft: length must be a power of two');
  if (n < 2) return;
  const { cos, sin, rev } = fftTable(n);
  for (let i = 1; i < n; i++) {
    const j = rev[i];
    if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
  }
  const sign = inverse ? 1 : -1;
  for (let len = 2; len <= n; len <<= 1) {
    const half = len >> 1, step = n / len;
    for (let i = 0; i < n; i += len) {
      for (let k = 0; k < half; k++) {
        const cr = cos[k * step], ci = sign * sin[k * step];
        const a = i + k, b = a + half;
        const xr = re[b] * cr - im[b] * ci, xi = re[b] * ci + im[b] * cr;
        re[b] = re[a] - xr; im[b] = im[a] - xi;
        re[a] += xr; im[a] += xi;
      }
    }
  }
  if (inverse) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
}

/** Angular wavenumbers k for an FFT grid of n points with spacing dx (FFT ordering). */
export function fftFreq(n, dx) {
  const k = new Float64Array(n);
  for (let i = 0; i < n; i++) k[i] = 2 * Math.PI * (i < n / 2 ? i : i - n) / (n * dx);
  return k;
}

// ---------------------------------------------------------------- ODEs

/** Classical RK4 step for y' = f(t, y). y: array; returns new array. */
export function rk4Step(f, t, y, h) {
  const k1 = f(t, y);
  const k2 = f(t + h / 2, y.map((v, i) => v + h / 2 * k1[i]));
  const k3 = f(t + h / 2, y.map((v, i) => v + h / 2 * k2[i]));
  const k4 = f(t + h, y.map((v, i) => v + h * k3[i]));
  return y.map((v, i) => v + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
}

// Dormand–Prince 5(4) coefficients (Dormand & Prince 1980; Hairer, Nørsett & Wanner, Table 5.2).
const DP = {
  c: [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1],
  a: [[], [1 / 5], [3 / 40, 9 / 40], [44 / 45, -56 / 15, 32 / 9],
    [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
    [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]],
  b: [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0],
  e: [71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40],
};

/** The Dormand–Prince 5(4) tableau (for steppers built on top of dopri5). */
export const DORMAND_PRINCE = DP;

/**
 * Adaptive Dormand–Prince 5(4) integration of y' = f(t, y) from t0 to t1.
 * Error control per component: |err_i| <= atol + rtol·max(|y_i|, |y_new_i|).
 * onStep(t, y) is called after every accepted step (return false to stop early).
 * Returns { t, y, steps, rejected, stopped }. Throws if tolerance cannot be met or
 * the step budget is exhausted; stopped is reserved for an intentional onStep stop.
 */
export function dopri5(f, t0, y0, t1, { rtol = 1e-9, atol = 1e-12, h0, hmax = Infinity, maxSteps = 1e6, onStep } = {}) {
  if (!Number.isFinite(t0) || !Number.isFinite(t1) || !y0.length || !Array.from(y0).every(Number.isFinite)) throw new Error('dopri5: finite times and a non-empty finite state are required');
  if (!(Number.isFinite(rtol) && rtol > 0 && Number.isFinite(atol) && atol > 0)) throw new Error('dopri5: tolerances must be finite and positive');
  if (!(hmax > 0) || (hmax !== Infinity && !Number.isFinite(hmax)) || (h0 != null && !(Number.isFinite(h0) && h0 !== 0))) throw new Error('dopri5: step sizes must be non-zero and finite (hmax may be Infinity)');
  if (!(Number.isSafeInteger(maxSteps) && maxSteps > 0)) throw new Error('dopri5: maxSteps must be a positive integer');
  const n = y0.length, dir = Math.sign(t1 - t0) || 1;
  let t = t0;
  const y = Float64Array.from(y0);
  let h = dir * Math.min(Math.abs(h0 ?? (t1 - t0) / 100), hmax);
  const k = Array.from({ length: 7 }, () => new Float64Array(n));
  const yt = new Float64Array(n), yn = new Float64Array(n);
  let steps = 0, rejected = 0, stopped = false;
  const maxAttempts = 4 * maxSteps;
  k[0].set(f(t, y));
  while (dir * (t1 - t) > 0 && steps < maxSteps && steps + rejected < maxAttempts) {
    if (dir * (t + h - t1) > 0) h = t1 - t;
    // This checks floating-point progress, not a dimensional time threshold. A
    // femtosecond step can be perfectly resolvable and still need error control.
    if (h === 0 || t + h === t) throw new Error(`dopri5: step underflow at t = ${t}; cannot meet tolerance`);
    for (let s = 1; s < 7; s++) {
      for (let i = 0; i < n; i++) {
        let acc = y[i];
        for (let j = 0; j < s; j++) acc += h * DP.a[s][j] * k[j][i];
        yt[i] = acc;
      }
      k[s].set(f(t + DP.c[s] * h, yt));
    }
    let err = 0;
    for (let i = 0; i < n; i++) {
      let acc = y[i], e = 0;
      for (let j = 0; j < 7; j++) { acc += h * DP.b[j] * k[j][i]; e += h * DP.e[j] * k[j][i]; }
      yn[i] = acc;
      const sc = atol + rtol * Math.max(Math.abs(y[i]), Math.abs(acc));
      err = Math.max(err, Math.abs(e) / sc);
    }
    if (!Number.isFinite(err)) throw new Error(`dopri5: non-finite derivative or state at t = ${t}`);
    if (err <= 1) {
      t += h; y.set(yn); k[0].set(k[6]); steps++;      // FSAL: k7 of this step = k1 of the next
      if (onStep && onStep(t, y) === false) { stopped = true; break; }
    } else rejected++;
    const fac = err === 0 ? 5 : Math.min(5, Math.max(0.2, 0.9 * Math.pow(err, -0.2)));
    h = dir * Math.min(Math.abs(h * fac), hmax);
  }
  if (!stopped && dir * (t1 - t) > 0) throw new Error(`dopri5: step budget exhausted at t = ${t}, before t = ${t1}`);
  return { t, y, steps, rejected, stopped };
}

// ---------------------------------------------------------------- quadrature

const GL_CACHE = new Map();
/** Gauss–Legendre nodes and weights on [−1, 1] (Newton iteration on P_n). */
export function gaussLegendreNodes(n) {
  if (GL_CACHE.has(n)) return GL_CACHE.get(n);
  const x = new Float64Array(n), w = new Float64Array(n);
  for (let i = 0; i < Math.ceil(n / 2); i++) {
    let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5)), dp = 0;
    for (let it = 0; it < 100; it++) {
      let p0 = 1, p1 = z;
      for (let k = 2; k <= n; k++) { const p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
      const pn = n === 1 ? z : p1, pm = n === 1 ? 1 : p0;
      dp = n * (z * pn - pm) / (z * z - 1);
      const dz = pn / dp;
      z -= dz;
      if (Math.abs(dz) < 1e-16) break;
    }
    x[i] = -z; x[n - 1 - i] = z;
    w[i] = w[n - 1 - i] = 2 / ((1 - z * z) * dp * dp);
  }
  const out = { x, w };
  GL_CACHE.set(n, out);
  return out;
}

/** ∫_a^b f(x) dx with n-point Gauss–Legendre (exact for polynomials of degree 2n−1). */
export function gaussLegendre(f, a, b, n = 32) {
  const { x, w } = gaussLegendreNodes(n);
  const m = (b - a) / 2, c = (a + b) / 2;
  let s = 0;
  for (let i = 0; i < n; i++) s += w[i] * f(m * x[i] + c);
  return s * m;
}

/** Adaptive Simpson quadrature to absolute tolerance tol. */
export function adaptiveSimpson(f, a, b, tol = 1e-10, maxDepth = 50) {
  const simpson = (fl, fmid, fr, lo, hi) => (hi - lo) / 6 * (fl + 4 * fmid + fr);
  const recurse = (lo, hi, fl, fmid, fr, whole, eps, depth) => {
    const mid = (lo + hi) / 2, lm = (lo + mid) / 2, rm = (mid + hi) / 2;
    const flm = f(lm), frm = f(rm);
    const left = simpson(fl, flm, fmid, lo, mid), right = simpson(fmid, frm, fr, mid, hi);
    if (depth <= 0 || Math.abs(left + right - whole) <= 15 * eps) return left + right + (left + right - whole) / 15;
    return recurse(lo, mid, fl, flm, fmid, left, eps / 2, depth - 1) + recurse(mid, hi, fmid, frm, fr, right, eps / 2, depth - 1);
  };
  const fa = f(a), fb = f(b), fm = f((a + b) / 2);
  return recurse(a, b, fa, fm, fb, simpson(fa, fm, fb, a, b), tol, maxDepth);
}

/** ∫_0^∞ f(x) dx via the substitution x = t/(1−t) and Gauss–Legendre panels. */
export function integrateToInfinity(f, { panels = 16, n = 32 } = {}) {
  const g = t => { const x = t / (1 - t); return f(x) / ((1 - t) * (1 - t)); };
  let s = 0;
  for (let p = 0; p < panels; p++) s += gaussLegendre(g, p / panels, (p + 1) / panels, n);   // nodes never hit t = 1
  return s;
}

// ---------------------------------------------------------------- roots

/** Brent's method: root of f in [a, b] where f(a)·f(b) < 0. */
export function brent(f, lo, hi, tol = 1e-14, maxIter = 200) {
  let a = lo, b = hi;                      // b: best estimate; a: previous; c: bracket partner
  let fa = f(a), fb = f(b);
  if (fa * fb > 0) throw new Error('brent: root not bracketed');
  let c = a, fc = fa, d = b - a, e = d;
  for (let i = 0; i < maxIter; i++) {
    if (fb * fc > 0) { c = a; fc = fa; d = e = b - a; }
    if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
    const tol1 = 2 * Number.EPSILON * Math.abs(b) + tol / 2, m = (c - b) / 2;
    if (Math.abs(m) <= tol1 || fb === 0) return b;
    if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
      let p, q, r;
      const s = fb / fa;
      if (a === c) { p = 2 * m * s; q = 1 - s; }
      else { q = fa / fc; r = fb / fc; p = s * (2 * m * q * (q - r) - (b - a) * (r - 1)); q = (q - 1) * (r - 1) * (s - 1); }
      if (p > 0) q = -q; else p = -p;
      if (2 * p < Math.min(3 * m * q - Math.abs(tol1 * q), Math.abs(e * q))) { e = d; d = p / q; }
      else { d = m; e = d; }
    } else { d = m; e = d; }
    a = b; fa = fb;
    b += Math.abs(d) > tol1 ? d : (m > 0 ? tol1 : -tol1);
    fb = f(b);
  }
  return b;
}

// ---------------------------------------------------------------- special functions

const LANCZOS = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
  -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];

/** ln Γ(x) for x > 0 (Lanczos, g = 7, n = 9; relative error ~1e-15). */
export function lgamma(z) {
  if (z <= 0 && Number.isInteger(z)) return Infinity;   // poles of Γ
  if (z < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * z))) - lgamma(1 - z);
  const x = z - 1;
  let a = LANCZOS[0];
  const t = x + 7.5;
  for (let i = 1; i < 9; i++) a += LANCZOS[i] / (x + i);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

export function gamma(x) {
  if (x <= 0 && Number.isInteger(x)) return NaN;        // poles of Γ
  if (x < 0.5) return Math.PI / (Math.sin(Math.PI * x) * gamma(1 - x));
  return Math.exp(lgamma(x));
}

export function factorial(n) {
  let f = 1;
  for (let i = 2; i <= n; i++) f *= i;
  return f;
}

/** Error function (Abramowitz–Stegun 7.1.26 is too coarse; this uses a series / continued fraction, |err| < 1e-14). */
export function erf(x) {
  const ax = Math.abs(x);
  if (ax > 6) return Math.sign(x);                      // erfc(6) ≈ 2e-17 < double precision near 1
  if (ax < 2.5) {
    let sum = ax, term = ax;
    for (let n = 1; n < 100; n++) {
      term *= -ax * ax / n;
      const add = term / (2 * n + 1);
      sum += add;
      if (Math.abs(add) < 1e-17 * Math.abs(sum)) break;
    }
    return Math.sign(x) * 2 / Math.sqrt(Math.PI) * sum;
  }
  // erfc continued fraction (Lentz)
  let f = ax, C = ax, D = 0;
  for (let n = 1; n < 200; n++) {
    const an = n / 2;
    D = ax + an * D; D = 1 / (D === 0 ? 1e-300 : D);
    C = ax + an / C; if (C === 0) C = 1e-300;
    const delta = C * D;
    f *= delta;
    if (Math.abs(delta - 1) < 1e-16) break;
  }
  const erfc = Math.exp(-ax * ax) / Math.sqrt(Math.PI) / f;
  return Math.sign(x) * (1 - erfc);
}

/**
 * Regularised lower incomplete gamma P(a, x) = γ(a, x)/Γ(a), a > 0, x ≥ 0. Series for x < a + 1,
 * Lentz continued fraction for Q otherwise (Numerical Recipes §6.2), |err| ≈ 1e-14.
 */
export function gammaP(a, x) {
  if (!(a > 0) || !(x >= 0)) throw new Error(`gammaP: need a > 0 and x ≥ 0 (got a = ${a}, x = ${x})`);
  return x < a + 1 ? gammaSeries(a, x) : 1 - gammaCF(a, x);
}

/** Regularised upper incomplete gamma Q(a, x) = 1 − P(a, x); the χ² upper-tail p-value is Q(dof/2, χ²/2). */
export function gammaQ(a, x) {
  if (!(a > 0) || !(x >= 0)) throw new Error(`gammaQ: need a > 0 and x ≥ 0 (got a = ${a}, x = ${x})`);
  return x < a + 1 ? 1 - gammaSeries(a, x) : gammaCF(a, x);
}

function gammaSeries(a, x) {
  if (x === 0) return 0;
  let ap = a, del = 1 / a, sum = del;
  for (let n = 0; n < 1000; n++) {
    ap++; del *= x / ap; sum += del;
    if (Math.abs(del) < Math.abs(sum) * 1e-16) break;
  }
  return sum * Math.exp(-x + a * Math.log(x) - lgamma(a));
}

function gammaCF(a, x) {
  if (x === Infinity) return 0;
  const tiny = 1e-300;
  let b = x + 1 - a, c = 1 / tiny, d = 1 / b, h = d;
  for (let i = 1; i < 1000; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b; if (Math.abs(d) < tiny) d = tiny;
    c = b + an / c; if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-16) break;
  }
  return Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
}

/** Modified Bessel function of the first kind I_n(x), integer n ≥ 0 (power series; fine for x ≲ 50). */
export function besselI(n, x) {
  let term = Math.pow(x / 2, n) / factorial(n), sum = term;
  for (let k = 1; k < 500; k++) {
    term *= (x * x / 4) / (k * (k + n));
    sum += term;
    if (Math.abs(term) < 1e-17 * Math.abs(sum)) break;
  }
  return sum;
}

/**
 * Bessel function of the first kind J_n(x), integer n ≥ 0, any real x. Miller's backward
 * recurrence normalised by J₀ + 2ΣJ₂ₖ = 1: relative accuracy also where J_n is tiny (x ≪ n),
 * absolute accuracy ~1e-15·√(x) for large x. J_n(−x) = (−1)ⁿJ_n(x).
 */
export function besselJ(n, x) {
  if (!(Number.isInteger(n) && n >= 0)) throw new Error('besselJ: order must be an integer ≥ 0');
  if (!Number.isFinite(x)) return Number.isNaN(x) ? NaN : 0;     // J_n(±∞) = 0
  if (x < 0) return (n % 2 ? -1 : 1) * besselJ(n, -x);
  if (x === 0) return n === 0 ? 1 : 0;
  const top = Math.max(n, x);
  const m = 2 * Math.ceil((top + 20 + Math.sqrt(60 * top)) / 2);    // even start index well above n and x
  let jp = 0, j = 1e-300, norm = 0, result = 0;
  for (let k = m; k > 0; k--) {
    const jm = (2 * k / x) * j - jp;                                  // J_{k−1} from J_k, J_{k+1}
    jp = j; j = jm;
    if (Math.abs(j) > 1e250) { j *= 1e-250; jp *= 1e-250; norm *= 1e-250; result *= 1e-250; }
    if (k - 1 === n) result = j;
    if ((k - 1) % 2 === 0 && k - 1 > 0) norm += 2 * j;
  }
  return result / (norm + j);
}

/**
 * First `count` positive zeros j_{n,1} < j_{n,2} < … of J_n (integer n ≥ 0): bracketed by a scan
 * finer than the zero spacing (always > π), refined with Brent. Mode frequencies of circular
 * membranes, waveguides and diffraction minima: ω_{n,s} = c·j_{n,s}/R.
 */
export function besselJZeros(n, count) {
  if (!(Number.isInteger(n) && n >= 0 && Number.isInteger(count) && count > 0)) throw new Error('besselJZeros: need integers n ≥ 0, count > 0');
  const zeros = [], step = 0.5;
  let a = n === 0 ? 1 : n, fa = besselJ(n, a);                        // j_{n,1} > n for n > 0, j_{0,1} > 2
  while (zeros.length < count) {
    const b = a + step, fb = besselJ(n, b);
    if (fa === 0) zeros.push(a);
    else if (fa * fb < 0) zeros.push(brent(x => besselJ(n, x), a, b, 1e-15));
    a = b; fa = fb;
  }
  return zeros;
}

/**
 * Complete elliptic integral of the first kind K(m) = ∫₀^{π/2} dθ / √(1 − m sin²θ), parameter
 * m = k² < 1, by the arithmetic–geometric mean: K = π / (2 AGM(1, √(1 − m))).
 */
export function ellipticK(m) {
  if (!(m < 1)) return m === 1 ? Infinity : NaN;
  let a = 1, b = Math.sqrt(1 - m);
  for (let i = 0; i < 40 && Math.abs(a - b) > 1e-16 * a; i++) [a, b] = [(a + b) / 2, Math.sqrt(a * b)];
  return Math.PI / (2 * a);
}

/** Generalised Laguerre polynomial L_n^(α)(x) by three-term recurrence. */
export function laguerre(n, alpha, x) {
  if (n === 0) return 1;
  let l0 = 1, l1 = 1 + alpha - x;
  for (let k = 1; k < n; k++) {
    const l2 = ((2 * k + 1 + alpha - x) * l1 - (k + alpha) * l0) / (k + 1);
    l0 = l1; l1 = l2;
  }
  return l1;
}

/** Associated Legendre P_l^m(x), m ≥ 0, with the Condon–Shortley phase (−1)^m. */
export function legendreP(l, m, x) {
  let pmm = 1;
  if (m > 0) {
    const s = Math.sqrt((1 - x) * (1 + x));
    let fact = 1;
    for (let i = 1; i <= m; i++) { pmm *= -fact * s; fact += 2; }
  }
  if (l === m) return pmm;
  let pmmp1 = x * (2 * m + 1) * pmm;
  if (l === m + 1) return pmmp1;
  let pll = 0;
  for (let ll = m + 2; ll <= l; ll++) {
    pll = ((2 * ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
    pmm = pmmp1; pmmp1 = pll;
  }
  return pll;
}

/**
 * Fully normalised associated Legendre function P̄_l^m(x) = N_lm P_l^m(x) with
 * N_lm = √((2l+1)/(4π) · (l−m)!/(l+m)!), Condon–Shortley phase included. Uses the standard
 * stable recurrences (no factorials), so it is accurate for large l and m.
 */
export function normalizedLegendre(l, m, x) {
  const s = Math.sqrt(Math.max(0, (1 - x) * (1 + x)));
  let pmm = 1 / Math.sqrt(4 * Math.PI);
  for (let i = 1; i <= m; i++) pmm *= -Math.sqrt((2 * i + 1) / (2 * i)) * s;
  if (l === m) return pmm;
  let pm1 = x * Math.sqrt(2 * m + 3) * pmm;
  if (l === m + 1) return pm1;
  let pl = 0;
  for (let ll = m + 2; ll <= l; ll++) {
    const a = Math.sqrt((4 * ll * ll - 1) / (ll * ll - m * m));
    const b = Math.sqrt(((ll - 1) ** 2 - m * m) / (4 * (ll - 1) ** 2 - 1));
    pl = a * (x * pm1 - b * pmm);
    pmm = pm1; pm1 = pl;
  }
  return pl;
}

/**
 * Real (tesseral) spherical harmonic, orthonormal on the sphere.
 * m > 0: √2 P̄_l^m cos(mφ); m < 0: √2 P̄_l^|m| sin(|m|φ); m = 0: P̄_l^0.
 * The Condon–Shortley phase is removed, so p_x ∝ +x, p_y ∝ +y.
 */
export function realYlm(l, m, theta, phi) {
  const am = Math.abs(m);
  const P = normalizedLegendre(l, am, Math.cos(theta)) * (am % 2 ? -1 : 1);
  if (m === 0) return P;
  return Math.SQRT2 * P * (m > 0 ? Math.cos(am * phi) : Math.sin(am * phi));
}

/**
 * Complex spherical harmonic Y_l^m(θ, φ) with the Condon–Shortley phase (the quantum-mechanics
 * convention): Y_l^m = P̄_l^|m|(cos θ) e^{imφ} for m ≥ 0 and Y_l^{−m} = (−1)^m conj(Y_l^m).
 * Orthonormal on the sphere; eigenfunctions of L² and L_z. Returns { re, im }.
 */
export function complexYlm(l, m, theta, phi) {
  if (!(Number.isInteger(l) && Number.isInteger(m) && l >= 0 && Math.abs(m) <= l)) throw new Error('complexYlm: need integers l ≥ 0, |m| ≤ l');
  const am = Math.abs(m), P = normalizedLegendre(l, am, Math.cos(theta));
  const sign = m < 0 && am % 2 ? -1 : 1;
  return { re: sign * P * Math.cos(m * phi), im: sign * P * Math.sin(m * phi) };
}

/**
 * Radial function R_nl(r) of a hydrogen-like atom (one electron, nuclear charge Z, point nucleus,
 * non-relativistic), normalised so ∫R²r²dr = 1:
 *   R_nl = √((2Z/na₀)³ (n−l−1)!/(2n (n+l)!)) e^{−ρ/2} ρ^l L_{n−l−1}^{(2l+1)}(ρ),  ρ = 2Zr/(na₀).
 * Lengths in units of a₀ unless a0 is given (use the reduced-mass a₀ for exact spectroscopy).
 */
export function hydrogenRadial(n, l, r, { Z = 1, a0 = 1 } = {}) {
  if (!(Number.isInteger(n) && Number.isInteger(l) && n >= 1 && l >= 0 && l < n)) throw new Error('hydrogenRadial: need integers n ≥ 1, 0 ≤ l < n');
  const rho = 2 * Z * r / (n * a0);
  const lnNorm = 1.5 * Math.log(2 * Z / (n * a0)) + 0.5 * (lgamma(n - l) - lgamma(n + l + 1) - Math.log(2 * n));
  return Math.exp(lnNorm - rho / 2) * rho ** l * laguerre(n - l - 1, 2 * l + 1, rho);
}

// ---------------------------------------------------------------- linear algebra

/** Solve a tridiagonal system (Thomas algorithm). a: sub, b: diag, c: super, d: rhs. Works for complex via callers. */
export function solveTridiagonal(a, b, c, d) {
  const n = b.length, cp = new Float64Array(n), dp = new Float64Array(n), x = new Float64Array(n);
  cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0];
  for (let i = 1; i < n; i++) {
    const m = b[i] - a[i] * cp[i - 1];
    cp[i] = c[i] / m;
    dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
  }
  x[n - 1] = dp[n - 1];
  for (let i = n - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];
  return x;
}

/**
 * All eigenvalues and eigenvectors of a real symmetric tridiagonal matrix (implicit QL with
 * Wilkinson shifts; Numerical Recipes "tqli"). diag: n diagonal entries; off: n−1 entries
 * coupling i and i+1. Returns { values, vectors } with vectors[i] the unit eigenvector of
 * values[i] (unsorted). Cost O(n³) with eigenvectors: fine to n ≈ 500 in a browser.
 * { vectors: false } returns eigenvalues only, in O(n²) (n = 10⁴ in ~3 s); for the lowest few use
 * lowestTridiagonal (bisection, O(n) per level).
 */
export function eigenTridiagonal(diag, off, { vectors: wantVectors = true } = {}) {
  const n = diag.length;
  if (off.length !== n - 1) throw new Error('eigenTridiagonal: off must have n − 1 entries');
  const d = Float64Array.from(diag), e = new Float64Array(n);
  e.set(off);
  const z = wantVectors ? new Float64Array(n * n) : null;
  if (z) for (let i = 0; i < n; i++) z[i * n + i] = 1;
  for (let l = 0; l < n; l++) {
    let iter = 0, m;
    do {
      for (m = l; m < n - 1; m++) {
        const dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
        if (Math.abs(e[m]) <= Number.EPSILON * dd) break;
      }
      if (m !== l) {
        if (iter++ === 60) throw new Error('eigenTridiagonal: no convergence');
        let g = (d[l + 1] - d[l]) / (2 * e[l]);
        let r = Math.hypot(g, 1);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? r : -r));
        let s = 1, c = 1, p = 0, i;
        for (i = m - 1; i >= l; i--) {
          let f = s * e[i];
          const b = c * e[i];
          r = Math.hypot(f, g);
          e[i + 1] = r;
          if (r === 0) { d[i + 1] -= p; e[m] = 0; break; }
          s = f / r; c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          if (z) {for (let k = 0; k < n; k++) {
            f = z[k * n + i + 1];
            z[k * n + i + 1] = s * z[k * n + i] + c * f;
            z[k * n + i] = c * z[k * n + i] - s * f;
          }}
        }
        if (r === 0 && i >= l) continue;
        d[l] -= p; e[l] = g; e[m] = 0;
      }
    } while (m !== l);
  }
  if (!z) return { values: Array.from(d) };
  const vectors = Array.from({ length: n }, (_v, i) => Float64Array.from({ length: n }, (_w, k) => z[k * n + i]));
  return { values: Array.from(d), vectors };
}

/**
 * The k lowest eigenvalues of a symmetric tridiagonal matrix by Sturm-sequence bisection:
 * O(n · 64) per eigenvalue, accurate to round-off, exact ordering (no missed or merged levels,
 * also for near-degenerate tunnelling doublets). Returns ascending values.
 */
export function lowestTridiagonal(diag, off, k) {
  const n = diag.length;
  if (off.length !== n - 1) throw new Error('lowestTridiagonal: off must have n − 1 entries');
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < n; i++) {                                 // Gershgorin bounds
    const r = (i > 0 ? Math.abs(off[i - 1]) : 0) + (i < n - 1 ? Math.abs(off[i]) : 0);
    lo = Math.min(lo, diag[i] - r); hi = Math.max(hi, diag[i] + r);
  }
  const e2 = Float64Array.from(off, v => v * v), tiny = Number.MIN_VALUE * 1e10 + 1e-300;
  const below = x => {                                           // number of eigenvalues < x
    let c = 0, q = diag[0] - x;
    if (q < 0) c++;
    for (let i = 1; i < n; i++) {
      q = diag[i] - x - e2[i - 1] / (q !== 0 ? q : tiny);
      if (q < 0) c++;
    }
    return c;
  };
  const out = [];
  for (let j = 0; j < Math.min(k, n); j++) {
    let a = j ? out[j - 1] : lo, b = hi;
    for (let it = 0; it < 200 && b - a > 2 * Number.EPSILON * Math.max(Math.abs(a), Math.abs(b), 1e-300); it++) {
      const m = 0.5 * (a + b);
      if (m === a || m === b) break;
      if (below(m) > j) b = m; else a = m;
    }
    out.push(0.5 * (a + b));
  }
  return out;
}

/** Least-squares fit y ≈ a + b x. Returns { a, b, r2 }. */
export function linearFit(xs, ys) {
  // two-pass (centred) sums: no cancellation when x has a large offset (e.g. timestamps)
  const n = xs.length;
  let mx = 0, my = 0;
  for (let i = 0; i < n; i++) { mx += xs[i]; my += ys[i]; }
  mx /= n; my /= n;
  let sxx = 0, sxy = 0, syy = 0;
  for (let i = 0; i < n; i++) { const dx = xs[i] - mx, dy = ys[i] - my; sxx += dx * dx; sxy += dx * dy; syy += dy * dy; }
  const b = sxy / sxx, a = my - b * mx;
  const r2 = syy === 0 ? 1 : (sxy * sxy) / (sxx * syy);
  return { a, b, r2 };
}

/**
 * Mean and standard error of a correlated series by binning (blocking) analysis. Returns the
 * largest error over block sizes with ≥ minBins blocks: conservative (≈ 10 % high for
 * uncorrelated data). Series too short for blocking fall back to the naive standard error.
 */
export function meanWithError(xs, minBins = 16) {
  const n = xs.length;
  const mean = xs.reduce((s, v) => s + v, 0) / n;
  if (n < minBins) {
    const v = n > 1 ? xs.reduce((s, x) => s + (x - mean) ** 2, 0) / (n - 1) : NaN;
    return { mean, error: Math.sqrt(v / n), naive: true };
  }
  let best = 0;
  for (let size = 1; n / size >= minBins; size *= 2) {
    const nb = Math.floor(n / size), means = [];
    for (let b = 0; b < nb; b++) { let s = 0; for (let i = 0; i < size; i++) s += xs[b * size + i]; means.push(s / size); }
    const m = means.reduce((s, v) => s + v, 0) / nb;
    const v = means.reduce((s, x) => s + (x - m) ** 2, 0) / (nb - 1);
    best = Math.max(best, Math.sqrt(v / nb));
  }
  return { mean, error: best };
}
