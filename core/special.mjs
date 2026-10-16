// Special functions of mathematical physics beyond core/numerics.mjs (which already has Γ, ln Γ,
// erf, the incomplete Γ, J_n, I_n, K(m), generalised Laguerre L_n^(α), associated Legendre and
// spherical harmonics). Plain real arithmetic, deterministic, no dependencies.
//
//   airy, airyAiZeros                      Ai, Ai′, Bi, Bi′ on the whole real line; zeros of Ai
//   sphericalBesselJ / sphericalBesselY    j_n, y_n (stable recurrences, small-x series)
//   besselY, besselK                        Y_n, K_n for integer n
//   hermiteH, oscillatorWavefunction       H_n and normalised ψ_n(ξ) for any n (no overflow)
//   hyp2f1, hyp1f1                         Gauss and Kummer hypergeometric functions
//   ellipticE, jacobiElliptic              complete E(m); sn, cn, dn, am
//   digamma, beta, logBeta, incompleteBeta, studentTCdf, studentTPValue
//   zeta, polylog                          Riemann ζ(s), Li_s(z)
//   expintE1, expintEn, expintEi           exponential integrals
//   wigner3j, wigner6j, clebschGordan      exact (BigInt) angular-momentum coupling
//   dawson                                 Dawson's integral F(x)
//
// Accuracy statements are relative unless they say "absolute"; oscillatory functions are only
// accurate in the absolute sense near their zeros (the zero itself is a rounding-sized target).
// Invalid input (non-finite arguments, poles, unsupported parameter regions) throws.

import { lgamma, gamma, ellipticK } from './numerics.mjs';

const EPS = Number.EPSILON;
const EULER = 0.5772156649015329;
const SQRT_PI = Math.sqrt(Math.PI);
const SQRT3 = Math.sqrt(3);

const need = (cond, message) => { if (!cond) throw new Error(message); };
const finite = (who, values) => {
  for (const [name, v] of Object.entries(values)) need(Number.isFinite(v), `${who}: ${name} must be a finite number (got ${v})`);
};
const isNonPosInt = v => Number.isInteger(v) && v <= 0;
const intArg = (who, name, v, min = 0) => need(Number.isInteger(v) && v >= min, `${who}: ${name} must be an integer ≥ ${min} (got ${v})`);

/** sin(πx) with exact argument reduction (exact zeros at integers). */
function sinpi(x) {
  if (Number.isInteger(x)) return 0;
  let r = x % 2;
  if (r > 1) r -= 2; else if (r < -1) r += 2;
  if (r > 0.5) r = 1 - r; else if (r < -0.5) r = -1 - r;
  return Math.sin(Math.PI * r);
}
/** cot(πx) with exact argument reduction. */
function cotpi(x) {
  const r = x - Math.round(x);
  return Math.cos(Math.PI * r) / Math.sin(Math.PI * r);
}
/** 1/Γ(z), zero at the poles. */
const rgamma = z => (isNonPosInt(z) ? 0 : 1 / gamma(z));
/** [sign, ln|Γ(z)|]; sign 0 at the poles (1/Γ = 0). */
function signLogGamma(z) {
  if (isNonPosInt(z)) return [0, Infinity];
  return [z > 0 || Math.ceil(-z) % 2 === 0 ? 1 : -1, lgamma(z)];
}

function gammaRatio(numerators, denominators) {
  let sign = 1, logValue = 0;
  for (const z of numerators) { const [s, l] = signLogGamma(z); sign *= s; logValue += l; }
  for (const z of denominators) { const [s, l] = signLogGamma(z); if (s === 0) return 0; sign *= s; logValue -= l; }
  return sign * Math.exp(logValue);
}

// ================================================================= Airy functions

const AI0 = 0.3550280538878172, AIP0 = -0.2588194037928068;
const AIRY_U = [1];
for (let k = 1; k < 60; k++) AIRY_U.push(AIRY_U[k - 1] * (6 * k - 5) * (6 * k - 3) * (6 * k - 1) / ((2 * k - 1) * 216 * k));
const AIRY_V = AIRY_U.map((u, k) => (k === 0 ? 1 : -(6 * k + 1) / (6 * k - 1) * u));
const AIRY_ASYMPTOTIC = 10;            // |x| ≥ 10: ζ ≥ 21, the smallest asymptotic term is ~e^{−2ζ} < 1e-18

/** Σ s_k c_k ζ^{−k} (s_k = (−1)^k if alternate), stopped at the smallest term. */
function asymptoticSum(coef, z, alternate, parity = -1) {
  // parity −1: all k; 0 / 1: only even / odd k, with sign (−1)^{(k−parity)/2}
  let sum = 0, last = Infinity, power = 1;
  for (let k = 0; k < coef.length; k++, power /= z) {
    if (parity >= 0 && k % 2 !== parity) continue;
    const sign = parity >= 0 ? (((k - parity) / 2) % 2 ? -1 : 1) : (alternate && k % 2 ? -1 : 1);
    const term = sign * coef[k] * power;
    if (Math.abs(term) > last) break;
    sum += term;
    if (Math.abs(term) <= 1e-17 * Math.abs(sum)) break;
    last = Math.abs(term);
  }
  return sum;
}

/** Maclaurin series (|x| ≲ 1 for Ai; any x ≥ 0 for Bi, whose terms are then all positive). */
function airySeries(x) {
  const x2 = x * x, x3 = x2 * x;
  let tf = 1, tg = x, dtf = 0, dtg = 1, f = 1, g = x, fp = 0, gp = 1;
  for (let k = 1; k < 400; k++) {
    dtf = tf * x2 / (3 * k - 1);
    dtg = tg * x2 / (3 * k);
    tf *= x3 / ((3 * k - 1) * 3 * k);
    tg *= x3 / (3 * k * (3 * k + 1));
    f += tf; g += tg; fp += dtf; gp += dtg;
    const scale = Math.abs(f) + Math.abs(g) + Math.abs(fp) + Math.abs(gp);
    if (Math.abs(tf) + Math.abs(tg) + Math.abs(dtf) + Math.abs(dtg) <= 1e-18 * scale) break;
  }
  const c1 = AI0, c2 = -AIP0;
  return { ai: c1 * f - c2 * g, aip: c1 * fp - c2 * gp, bi: SQRT3 * (c1 * f + c2 * g), bip: SQRT3 * (c1 * fp + c2 * gp) };
}

/** One Taylor step of y″ = x y from x0 by h: returns [y, y′]. */
function airyTaylorStep(x0, y, yp, h) {
  const h2 = h * h, h3 = h2 * h, scale = Math.abs(y) + Math.abs(yp * h);
  let tkm1 = 0, tk = y, tk1 = yp * h, sum = y + tk1, dsum = tk1;
  for (let k = 0; k < 300; k++) {
    const t2 = (x0 * h2 * tk + h3 * tkm1) / ((k + 2) * (k + 1));
    sum += t2; dsum += (k + 2) * t2;
    tkm1 = tk; tk = tk1; tk1 = t2;
    if (k > 4 && Math.abs(tk) + Math.abs(tk1) <= 1e-19 * Math.max(scale, Math.abs(sum))) break;
  }
  return [sum, dsum / h];
}

/** Integrate one solution pair (y, y′) of the Airy equation from x0 to x1 by Taylor steps ≤ 0.5. */
function airyTaylor(x0, x1, y, yp) {
  const n = Math.max(1, Math.ceil(Math.abs(x1 - x0) / 0.5)), h = (x1 - x0) / n;
  let v = [y, yp];
  for (let i = 0; i < n; i++) v = airyTaylorStep(x0 + i * h, v[0], v[1], h);
  return v;
}

function airyAsymptoticPositive(x) {
  const r = Math.sqrt(x), q = Math.sqrt(r), zt = 2 / 3 * x * r;
  const su = asymptoticSum(AIRY_U, zt, true), sv = asymptoticSum(AIRY_V, zt, true);
  const bu = asymptoticSum(AIRY_U, zt, false), bv = asymptoticSum(AIRY_V, zt, false);
  const em = Math.exp(-zt / 2), ep = Math.exp(zt / 2);
  return {
    ai: (em / (2 * SQRT_PI * q) * su) * em, aip: (-q / (2 * SQRT_PI) * em * sv) * em,
    bi: (ep / (SQRT_PI * q) * bu) * ep, bip: (q / SQRT_PI * ep * bv) * ep,
  };
}

function airyAsymptoticNegative(z) {
  const r = Math.sqrt(z), q = Math.sqrt(r), zt = 2 / 3 * z * r;
  const P = asymptoticSum(AIRY_U, zt, false, 0), Q = asymptoticSum(AIRY_U, zt, false, 1);
  const Pv = asymptoticSum(AIRY_V, zt, false, 0), Qv = asymptoticSum(AIRY_V, zt, false, 1);
  const cz = Math.cos(zt), sz = Math.sin(zt);
  const c = (cz + sz) * Math.SQRT1_2, s = (sz - cz) * Math.SQRT1_2;   // cos, sin of ζ − π/4
  return {
    ai: (c * P + s * Q) / (SQRT_PI * q), bi: (-s * P + c * Q) / (SQRT_PI * q),
    aip: q / SQRT_PI * (s * Pv - c * Qv), bip: q / SQRT_PI * (c * Pv + s * Qv),
  };
}

/**
 * Airy functions Ai(x), Ai′(x), Bi(x), Bi′(x) for any finite real x. Returns { ai, aip, bi, bip }.
 * Method: Maclaurin series for |x| ≤ 1; asymptotic expansions (DLMF 9.7.5–9.7.11) for |x| ≥ 10;
 * in between, Taylor-series integration of y″ = xy in steps ≤ 0.5 — forward from x = 0 for x < 0
 * (both solutions oscillate: stable), backward from x = 10 for Ai at x > 0 (the recessive solution
 * is stable backward), and the all-positive Maclaurin series for Bi at 0 < x < 10.
 * Accuracy: ~1e-14 relative for x ≥ 0 (Ai, Bi and derivatives); for x < 0 ~1e-14 absolute relative
 * to the envelope |x|^{−1/4}/√π (|x|^{1/4}/√π for derivatives), and for x < −10 the phase
 * ζ = (2/3)|x|^{3/2} carries ~eps·ζ absolute error. Ai underflows to 0 beyond x ≈ 104, Bi
 * overflows to Infinity there. Wronskian Ai·Bi′ − Ai′·Bi = 1/π is a convenient check.
 */
export function airy(x) {
  finite('airy', { x });
  if (Math.abs(x) <= 1) return airySeries(x);
  if (x > 0) {
    const bi = x < AIRY_ASYMPTOTIC ? airySeries(x) : airyAsymptoticPositive(x);
    let ai, aip;
    if (x >= AIRY_ASYMPTOTIC) ({ ai, aip } = bi);
    else {
      const start = airyAsymptoticPositive(AIRY_ASYMPTOTIC);
      [ai, aip] = airyTaylor(AIRY_ASYMPTOTIC, x, start.ai, start.aip);
    }
    return { ai, aip, bi: bi.bi, bip: bi.bip };
  }
  if (x <= -AIRY_ASYMPTOTIC) return airyAsymptoticNegative(-x);
  const [ai, aip] = airyTaylor(0, x, AI0, AIP0);
  const [bi, bip] = airyTaylor(0, x, SQRT3 * AI0, -SQRT3 * AIP0);
  return { ai, aip, bi, bip };
}

/**
 * The first `count` zeros a_1 > a_2 > … of Ai (all negative): asymptotic estimate
 * a_s ≈ −T(3π(4s−1)/8) (DLMF 9.9.6, 9.9.18) refined by Newton's method with Ai″ = x·Ai.
 * Accuracy ~1e-15 relative. Throws for count < 1 or non-integer.
 */
export function airyAiZeros(count) {
  intArg('airyAiZeros', 'count', count, 1);
  const zeros = [];
  for (let s = 1; s <= count; s++) {
    const t = 3 * Math.PI / 8 * (4 * s - 1), t2 = 1 / (t * t);
    let x = -Math.pow(t, 2 / 3) * (1 + t2 * (5 / 48 + t2 * (-5 / 36 + t2 * (77125 / 82944 - t2 * 108056875 / 6967296))));
    for (let i = 0; i < 20; i++) {
      const { ai, aip } = airy(x), dx = ai / aip;
      x -= dx;
      if (Math.abs(dx) <= 2 * EPS * Math.abs(x)) break;
    }
    zeros.push(x);
  }
  return zeros;
}

// ================================================================= Bessel family

/**
 * Spherical Bessel function of the first kind j_n(x) = √(π/2x) J_{n+½}(x), integer n ≥ 0, real x.
 * Power series for x² < 2n+3 (relative accuracy even where j_n ~ xⁿ is tiny), upward recurrence
 * for x > n, otherwise Miller's downward recurrence normalised to j₀ or j₁ (whichever is larger).
 * Accuracy ~1e-15 relative away from zeros (absolute ~1e-16/x near zeros). j_n(−x) = (−1)ⁿ j_n(x).
 */
export function sphericalBesselJ(n, x) {
  intArg('sphericalBesselJ', 'n', n); finite('sphericalBesselJ', { x });
  if (x < 0) return (n % 2 ? -1 : 1) * sphericalBesselJ(n, -x);
  if (x === 0) return n === 0 ? 1 : 0;
  if (x * x < 2 * n + 3) {
    let pre = 1;
    for (let k = 1; k <= n; k++) pre *= x / (2 * k + 1);
    if (pre === 0) return 0;
    let term = 1, sum = 1;
    const q = -x * x / 2;
    for (let k = 1; k < 200; k++) {
      term *= q / (k * (2 * n + 2 * k + 1));
      sum += term;
      if (Math.abs(term) <= 1e-17 * Math.abs(sum)) break;
    }
    return pre * sum;
  }
  const s = Math.sin(x), c = Math.cos(x), j0 = s / x;
  if (n === 0) return j0;
  const j1 = (j0 - c) / x;
  if (n === 1) return j1;
  if (x > n) {
    let a = j0, b = j1;
    for (let k = 1; k < n; k++) [a, b] = [b, (2 * k + 1) / x * b - a];
    return b;
  }
  const top = Math.max(n, x), m = Math.ceil(top + 20 + Math.sqrt(60 * top));
  let next = 0, cur = 1e-300, at = 0, m1 = 0;
  for (let k = m; k > 0; k--) {
    const prev = (2 * k + 1) / x * cur - next;          // j_{k−1}
    next = cur; cur = prev;
    if (Math.abs(cur) > 1e250) { cur *= 1e-250; next *= 1e-250; at *= 1e-250; m1 *= 1e-250; }
    if (k - 1 === n) at = cur;
    if (k - 1 === 1) m1 = cur;
  }
  const m0 = cur;
  return Math.abs(j0) >= Math.abs(j1) ? at * (j0 / m0) : at * (j1 / m1);
}

/**
 * Spherical Bessel function of the second kind y_n(x), integer n ≥ 0, x ≠ 0 (y_n(−x) = (−1)^{n+1} y_n(x)).
 * Closed forms for y₀, y₁ and upward recurrence (stable: y_n is the dominant solution).
 * Accuracy ~1e-15 relative away from zeros; overflows to −Infinity for x ≪ n. x = 0 gives −Infinity.
 */
export function sphericalBesselY(n, x) {
  intArg('sphericalBesselY', 'n', n); finite('sphericalBesselY', { x });
  if (x === 0) return -Infinity;
  if (x < 0) return (n % 2 ? 1 : -1) * sphericalBesselY(n, -x);
  const s = Math.sin(x), c = Math.cos(x);
  let a = -c / x;
  if (n === 0) return a;
  let b = -(c / x + s) / x;
  for (let k = 1; k < n && Number.isFinite(b); k++) [a, b] = [b, (2 * k + 1) / x * b - a];
  return b;
}

/** J_0 … J_m (m chosen internally ≥ max(M, x) + margin) for x > 0 by Miller's algorithm. */
function besselJArray(M, x) {
  const top = Math.max(M, x), m = 2 * Math.ceil((top + 20 + Math.sqrt(60 * top)) / 2);
  const J = new Float64Array(m + 2);
  J[m] = 1e-300;
  let norm = 0;
  for (let k = m; k > 0; k--) {
    J[k - 1] = (2 * k / x) * J[k] - J[k + 1];
    if (Math.abs(J[k - 1]) > 1e250) { for (let i = k - 1; i <= m; i++) J[i] *= 1e-250; norm *= 1e-250; }
    if ((k - 1) % 2 === 0 && k > 1) norm += 2 * J[k - 1];
  }
  norm += J[0];
  for (let i = 0; i < J.length; i++) J[i] /= norm;
  return J;
}

/** Hankel asymptotic P, Q for order ν (DLMF 10.17.3), stopped at the smallest term. */
function hankelPQ(nu, x) {
  const mu = 4 * nu * nu;
  let P = 0, Q = 0, term = 1, last = Infinity;
  for (let k = 0; k < 200; k++) {
    if (k > 0) term *= (mu - (2 * k - 1) ** 2) / (k * 8 * x);
    if (Math.abs(term) > last) break;
    const sign = Math.floor(k / 2) % 2 ? -1 : 1;
    if (k % 2 === 0) P += sign * term; else Q += sign * term;
    if (Math.abs(term) <= 1e-17 * Math.abs(P)) break;
    last = Math.abs(term);
  }
  return { P, Q };
}

/** [Y₀(x), Y₁(x)] for x > 0. */
function besselY01(x) {
  if (x < 1e-8) {
    const L = Math.log(x / 2) + EULER;
    return [2 / Math.PI * L, -2 / (Math.PI * x)];
  }
  if (x >= 25) {
    const c = Math.cos(x), s = Math.sin(x), amp = Math.sqrt(2 / (Math.PI * x));
    const out = [];
    for (const nu of [0, 1]) {
      const phi = (nu / 2 + 0.25) * Math.PI, cp = Math.cos(phi), sp = Math.sin(phi);
      const cosChi = c * cp + s * sp, sinChi = s * cp - c * sp;          // χ = x − φ
      const { P, Q } = hankelPQ(nu, x);
      out.push(amp * (P * sinChi + Q * cosChi));
    }
    return out;
  }
  const J = besselJArray(1, x), L = Math.log(x / 2) + EULER;
  let s0 = 0, s1 = 0;
  for (let k = 1; 2 * k + 1 < J.length; k++) {
    const sign = k % 2 ? -1 : 1;
    s0 += sign * J[2 * k] / k;
    s1 += sign * (J[2 * k - 1] - J[2 * k + 1]) / k;
  }
  return [2 / Math.PI * (L * J[0] - 2 * s0), 2 / Math.PI * (-J[0] / x + L * J[1] + s1)];
}

/**
 * Bessel function of the second kind Y_n(x), integer n (Y_{−n} = (−1)ⁿY_n), x > 0 (x = 0: −Infinity).
 * Y₀, Y₁ from Neumann series over Miller-normalised J_k (x < 25) or Hankel's asymptotic expansion
 * (x ≥ 25, error ~e^{−2x}), then upward recurrence (stable). Accuracy ~1e-15 absolute relative to the
 * envelope √(2/πx) in the oscillatory region (relative away from zeros), ~1e-15 relative for x < n.
 */
export function besselY(n, x) {
  need(Number.isInteger(n), `besselY: order must be an integer (got ${n})`);
  finite('besselY', { x });
  need(x >= 0, `besselY: x must be ≥ 0 (got ${x}); Y_n is complex for x < 0`);
  if (n < 0) return (n % 2 ? -1 : 1) * besselY(-n, x);
  if (x === 0) return -Infinity;
  let [a, b] = besselY01(x);
  if (n === 0) return a;
  for (let k = 1; k < n && Number.isFinite(b); k++) [a, b] = [b, 2 * k / x * b - a];
  return b;
}

/** [K₀(x), K₁(x)] for x > 0. */
function besselK01(x, scaled = false) {
  if (x <= 2) {
    const q = x * x / 4, L = Math.log(x / 2);
    let t = 1, I0 = 1, H = 0, S0 = 0;
    for (let k = 1; k < 100; k++) {
      t *= q / (k * k); H += 1 / k; I0 += t; S0 += H * t;
      if (t <= 1e-17 * I0) break;
    }
    let u = 1, I1 = 1, p1 = -EULER, p2 = 1 - EULER, S1 = p1 + p2;
    for (let k = 1; k < 100; k++) {
      u *= q / (k * (k + 1)); p1 += 1 / k; p2 += 1 / (k + 1);
      I1 += u; S1 += (p1 + p2) * u;
      if (u <= 1e-17 * I1) break;
    }
    return [-(L + EULER) * I0 + S0, 1 / x + L * (x / 2) * I1 - x / 4 * S1];
  }
  // Steed's continued fraction CF2 with μ = 0 (Numerical Recipes §6.7, Temme 1975)
  let b = 2 * (1 + x), d = 1 / b, h = d, delh = d, q1 = 0, q2 = 1;
  const a1 = 0.25;
  let q = a1, c = a1, a = -a1, s = 1 + q * delh, i = 2;
  for (; i < 10000; i++) {
    a -= 2 * (i - 1);
    c = -a * c / i;
    const qnew = (q1 - b * q2) / a;
    q1 = q2; q2 = qnew; q += c * qnew;
    b += 2; d = 1 / (b + a * d);
    delh = (b * d - 1) * delh; h += delh;
    const dels = q * delh;
    s += dels;
    if (Math.abs(dels / s) < 1e-17) break;
  }
  need(i < 10000, 'besselK: continued fraction did not converge');
  h *= a1;
  const k0 = Math.sqrt((Math.PI / 2) / x) * (scaled ? 1 : Math.exp(-x)) / s;
  return [k0, k0 * (1 + (0.5 - h) / x)];
}

/**
 * Modified Bessel function of the second kind K_n(x), integer n (K_{−n} = K_n), x > 0 (x = 0: Infinity).
 * K₀, K₁ by their logarithmic series (x ≤ 2) or Steed's continued fraction (x > 2), then upward
 * recurrence (stable, exponentially scaled for x > 2). Scaling preserves high-order values
 * even when K₀ and K₁ themselves underflow. Accuracy ~n·eps relative away from under/overflow.
 */
export function besselK(n, x) {
  need(Number.isInteger(n), `besselK: order must be an integer (got ${n})`);
  finite('besselK', { x });
  need(x >= 0, `besselK: x must be ≥ 0 (got ${x})`);
  if (n < 0) return besselK(-n, x);
  if (x === 0) return Infinity;
  let [a, b] = besselK01(x, x > 2), logScale = x > 2 ? -x : 0;
  for (let k = 1; k < n && Number.isFinite(b); k++) {
    if (Math.abs(b) > 1e200) { a *= 1e-200; b *= 1e-200; logScale += 200 * Math.LN10; }
    [a, b] = [b, a + 2 * k / x * b];
  }
  const v = n === 0 ? a : b;
  if (logScale === 0) return v;
  const half = Math.exp(logScale / 2);
  if (half > 0 && Number.isFinite(half)) return (v * half) * half;
  return Math.exp(Math.log(v) + logScale);
}

// ================================================================= Hermite

/**
 * Physicists' Hermite polynomial H_n(x) (H₀ = 1, H₁ = 2x, H_{k+1} = 2xH_k − 2kH_{k−1}).
 * Exact for integer x and small n; relative accuracy ~n·eps otherwise. Large n·|x| overflows
 * (H_n ~ (2x)ⁿ): use oscillatorWavefunction for normalised quantities.
 */
export function hermiteH(n, x) {
  intArg('hermiteH', 'n', n); finite('hermiteH', { x });
  if (n === 0) return 1;
  let a = 1, b = 2 * x;
  for (let k = 1; k < n; k++) [a, b] = [b, 2 * x * b - 2 * k * a];
  return b;
}

/**
 * Normalised harmonic-oscillator eigenfunction ψ_n(ξ) = (2ⁿ n! √π)^{−½} H_n(ξ) e^{−ξ²/2}
 * (∫ψ_m ψ_n dξ = δ_mn) in the dimensionless coordinate ξ = x √(mω/ħ); in physical units
 * ψ(x) = (mω/ħ)^{1/4} ψ_n(ξ). Normalised three-term recurrence with a running exponent, so it
 * neither overflows nor underflows prematurely for any n (tested to n = 2000) and |ξ| up to ~35 + √(2n).
 * Accuracy ~√n·eps relative to max|ψ_n| (absolute near nodes).
 */
export function oscillatorWavefunction(n, xi) {
  intArg('oscillatorWavefunction', 'n', n); finite('oscillatorWavefunction', { xi });
  let logScale = -xi * xi / 2, a = 0, b = Math.pow(Math.PI, -0.25);
  for (let k = 0; k < n; k++) {
    const c = Math.sqrt(2 / (k + 1)) * xi * b - Math.sqrt(k / (k + 1)) * a;
    a = b; b = c;
    const m = Math.abs(b);
    if (m > 1e150) { a /= m; b /= m; logScale += Math.log(m); }
  }
  if (b === 0) return b;
  // exp(logScale) may underflow before multiplication by the still-large recurrence value.
  // Split the exponential when possible, preserving subnormal as well as normal results.
  const half = Math.exp(logScale / 2);
  if (half > 0 && Number.isFinite(half)) return (b * half) * half;
  return Math.sign(b) * Math.exp(Math.log(Math.abs(b)) + logScale);
}

// ================================================================= hypergeometric functions

const CANCEL = 1e4;     // accept a sum only if its largest term is ≤ 1e4 × |sum| (≤ 4 digits lost)

// A terminating hypergeometric series can cancel by hundreds of digits, even at degree 50.
// Sum its rational terms exactly for the supplied binary64 parameters, then round only once.
function binaryFraction(x) {
  if (x === 0) return [0n, 1n];
  const buf = new DataView(new ArrayBuffer(8));
  buf.setFloat64(0, x);
  const bits = buf.getBigUint64(0), exp = Number((bits >> 52n) & 2047n);
  let num = (bits & ((1n << 52n) - 1n)) | (exp ? 1n << 52n : 0n);
  if (x < 0) num = -num;
  const power = exp ? exp - 1075 : -1074;
  return power >= 0 ? [num << BigInt(power), 1n] : [num, 1n << BigInt(-power)];
}
function hyperPolynomial(numerators, denominator, x, degree) {
  need(Number.isSafeInteger(degree) && degree <= 2048, 'hypergeometric polynomial: degree exceeds the exact-summation budget (2048)');
  const pars = numerators.map(binaryFraction), [cn, cd] = binaryFraction(denominator), [xn, xd] = binaryFraction(x);
  let tn = 1n, td = 1n, sn = 1n, sd = 1n;
  for (let k = 0; k < degree; k++) {
    const kk = BigInt(k);
    for (const [pn, pd] of pars) { tn *= pn + kk * pd; td *= pd; }
    tn *= cd * xn; td *= (cn + kk * cd) * BigInt(k + 1) * xd;
    if (tn === 0n) break;
    if (td < 0n) { tn = -tn; td = -td; }
    let g = gcdBig(tn, td); tn /= g; td /= g;
    g = gcdBig(sd, td);
    sn = sn * (td / g) + tn * (sd / g); sd *= td / g;
    g = gcdBig(sn, sd); sn /= g; sd /= g;
  }
  return (sn < 0n ? -1 : 1) * ratioToNumber(absBig(sn), sd);
}

/** Σ t_k with t_{k+1} = t_k·ratio(k); returns { sum, big } or null if not converged. */
function hyperSeries(ratio, x, maxTerms = 200000) {
  let term = 1, sum = 1, big = 1, small = 0;
  const tail = Math.abs(x) < 1 ? 1 / (1 - Math.abs(x)) : 1;
  for (let k = 0; k < maxTerms; k++) {
    term *= ratio(k);
    if (term === 0) return { sum, big };
    sum += term;
    big = Math.max(big, Math.abs(term));
    if (!Number.isFinite(sum)) return null;
    small = Math.abs(term) * tail <= 1e-17 * Math.abs(sum) ? small + 1 : 0;
    if (small >= 2) return { sum, big };
  }
  return null;
}

function f21Series(a, b, c, x) {
  return hyperSeries(k => (a + k) * (b + k) / ((c + k) * (k + 1)) * x, x);
}

/** 2F1 for 0 ≤ x < 1 (non-terminating a, b). */
function f21Unit(a, b, c, x, who) {
  if (isNonPosInt(a) || isNonPosInt(b)) return hyp2f1(a, b, c, x);
  if (x <= 0.5) {
    const r = f21Series(a, b, c, x);
    need(r && r.big <= CANCEL * Math.abs(r.sum), `${who}: series cancellation too large; unsupported parameters`);
    return r.sum;
  }
  const d = c - a - b;
  if (!Number.isInteger(d)) {                                   // DLMF 15.8.4 (1 − x transformation)
    need(Math.max(Math.abs(a), Math.abs(b), Math.abs(c), Math.abs(d)) < 160, `${who}: parameters too large for the 1 − x transformation`);
    const y = 1 - x;
    const r1 = f21Series(a, b, 1 - d, y), r2 = f21Series(c - a, c - b, d + 1, y);
    if (r1 && r2) {
      const t1 = gamma(c) * gamma(d) * rgamma(c - a) * rgamma(c - b) * r1.sum;
      const t2 = Math.pow(y, d) * gamma(c) * gamma(-d) * rgamma(a) * rgamma(b) * r2.sum;
      const value = t1 + t2, big = Math.max(Math.abs(t1), Math.abs(t2), Math.abs(t1) * r1.big / Math.abs(r1.sum), Math.abs(t2) * r2.big / Math.abs(r2.sum));
      if (Number.isFinite(value) && big <= CANCEL * Math.abs(value)) return value;
    }
  }
  const r = f21Series(a, b, c, x);                                // slow but valid for x < 1
  need(r && r.big <= CANCEL * Math.abs(r.sum), `${who}: x = ${x} too close to 1 for these parameters (c − a − b integer); unsupported`);
  return r.sum;
}

/**
 * Gauss hypergeometric function ₂F₁(a, b; c; x) for real parameters and real x ≤ 1.
 * Supported: terminating cases (a or b a non-positive integer: polynomial, any real x);
 * |x| ≤ ½ by the series; ½ < x < 1 by the 1 − x connection formula (DLMF 15.8.4) when c − a − b is
 * not an integer, otherwise the (slowly converging) series; x < −½ by Pfaff's transformation
 * (DLMF 15.8.1) onto (⅓, 1); x = 1 by Gauss's sum when c − a − b > 0.
 * Accuracy ~1e-14 relative typically, ~1e-12 for parameters of size ~50; every path checks the
 * cancellation it incurs and throws instead of returning fewer than ~12 correct digits.
 * Throws for x > 1 (branch cut, complex), c a non-positive integer (pole), divergence at x = 1,
 * or regions it cannot evaluate accurately. Terminating polynomials of degree ≤ 2048 use exact
 * rational arithmetic for the input doubles (one final rounding); larger degrees throw.
 */
export function hyp2f1(a, b, c, x) {
  const who = 'hyp2f1';
  finite(who, { a, b, c, x });
  const poly = isNonPosInt(a) || isNonPosInt(b);
  const degree = poly ? Math.max(isNonPosInt(a) ? a : -Infinity, isNonPosInt(b) ? b : -Infinity) : 0;
  need(!isNonPosInt(c) || (poly && degree > c), `${who}: c = ${c} is a non-positive integer (pole)`);
  if (x === 0) return 1;
  if (poly) {
    return hyperPolynomial([a, b], c, x, -degree);
  }
  need(x <= 1, `${who}: x = ${x} > 1 lies on the branch cut (complex result)`);
  if (x === 1) {
    need(c - a - b > 0, `${who}: diverges at x = 1 unless c − a − b > 0`);
    return gammaRatio([c, c - a - b], [c - a, c - b]);
  }
  if (x >= -0.5) {
    if (x >= 0) return f21Unit(a, b, c, x, who);
    const r = f21Series(a, b, c, x);
    need(r && r.big <= CANCEL * Math.abs(r.sum), `${who}: series cancellation too large; unsupported parameters`);
    return r.sum;
  }
  const z = x / (x - 1);                                          // Pfaff: z ∈ (1/3, 1)
  if (isNonPosInt(c - a) && !isNonPosInt(c - b)) return Math.pow(1 - x, -b) * f21Unit(c - a, b, c, z, who);
  return Math.pow(1 - x, -a) * f21Unit(a, c - b, c, z, who);
}

/** Asymptotic Σ (p)_k (q)_k / (k! w^k) stopped at the smallest term: { sum, exact } or null if not converged. */
function asymptoticPQ(p, q, w) {
  let term = 1, sum = 1;
  for (let k = 0; k < 500; k++) {
    const next = term * (p + k) * (q + k) / ((k + 1) * w);
    if (next === 0) return { sum, exact: true };
    if (Math.abs(next) > Math.abs(term)) return null;
    term = next; sum += term;
    if (Math.abs(term) <= 1e-17 * Math.abs(sum)) return { sum, exact: false };
  }
  return null;
}

/**
 * Large-|x| asymptotics of M(a, b, x) (DLMF 13.7.2 with Kummer's transformation), including the
 * exponentially recessive term with the Stokes-line average cos(πa) on the real axis. Accepted only
 * when that term is negligible (< 1e-16 relative) or the dominant series terminates (exact); else null.
 */
function m1f1Asymptotic(a, b, x) {
  const y = Math.abs(x), ly = Math.log(y);
  const [sb, lb] = signLogGamma(b), [sa, la] = signLogGamma(a), [sba, lba] = signLogGamma(b - a);
  // x > 0: dominant e^y y^{a−b}/Γ(a) Σ(1−a)_k(b−a)_k/(k!y^k), recessive cos(πa) y^{−a}/Γ(b−a) Σ(a)_k(a−b+1)_k/(k!(−y)^k)
  // x < 0: the same with a → b − a (Kummer), the e^{−y} factor moving to the recessive part
  const [p1, q1, p2, q2, sd, ld, sr, lr, cosr, ed, er] = x > 0
    ? [1 - a, b - a, a, a - b + 1, sa, la, sba, lba, sinpi(a + 0.5), y + (a - b) * ly, -a * ly]
    : [a, a - b + 1, b - a, 1 - a, sba, lba, sa, la, sinpi(b - a + 0.5), -a * ly, -y + (a - b) * ly];
  if (sd === 0) return null;
  const S1 = asymptoticPQ(p1, q1, y);
  if (S1 === null) return null;
  const dom = sb * sd * Math.exp(lb - ld + ed) * S1.sum;
  if (sr === 0 || cosr === 0) return dom;
  const relRec = Math.exp(ld - lr + er - ed);                      // |recessive/dominant| prefactor ratio
  if (relRec < 1e-17) return dom;
  if (!S1.exact) return null;
  const S2 = asymptoticPQ(p2, q2, -y);
  if (S2 === null) return null;
  return dom + cosr * sb * sr * Math.exp(lb - lr + er) * S2.sum;
}

function m1f1Series(a, b, x) {
  return hyperSeries(k => (a + k) / ((b + k) * (k + 1)) * x, 0, 100000);
}

/**
 * Kummer's confluent hypergeometric function M(a, b, x) = ₁F₁(a; b; x), real a, b, x.
 * Paths: polynomial when a is a non-positive integer (and e^x × polynomial when b − a is, x < 0);
 * the power series (x ≥ 0) or Kummer's transformation M = e^x M(b − a, b, −x) (x < 0) for |x| ≤ 700;
 * the large-|x| asymptotic expansion (DLMF 13.7.2) for |x| ≥ 30 whenever it is accurate to full precision
 * (its exponentially small second term is either negligible or the leading series terminates).
 * Accuracy ~1e-14 relative typically; each path checks its cancellation and the function throws
 * rather than return fewer than ~12 digits. Throws for b a non-positive integer (pole).
 * Overflow for x ≳ 700 returns Infinity.
 */
export function hyp1f1(a, b, x) {
  const who = 'hyp1f1';
  finite(who, { a, b, x });
  need(!isNonPosInt(b) || (isNonPosInt(a) && a > b), `${who}: b = ${b} is a non-positive integer (pole)`);
  if (x === 0) return 1;
  const poly = (aa, xx) => hyperPolynomial([aa], b, xx, -aa);
  if (isNonPosInt(a)) return poly(a, x);
  if (x < 0 && isNonPosInt(b - a)) return Math.exp(x) * poly(b - a, -x);
  if (Math.abs(x) >= 30) {
    const v = m1f1Asymptotic(a, b, x);
    if (v !== null) return v;
  }
  if (Math.abs(x) <= 700) {
    const tries = x >= 0 ? [[a, x, 1]] : [[b - a, -x, Math.exp(x)], [a, x, 1]];
    for (const [aa, xx, pre] of tries) {
      const r = m1f1Series(aa, b, xx);
      if (r && r.big <= CANCEL * Math.abs(r.sum)) return pre * r.sum;
    }
  }
  throw new Error(`${who}: cannot evaluate M(${a}, ${b}, ${x}) to ~12 digits (cancellation); unsupported region`);
}

// ================================================================= elliptic functions

/**
 * Complete elliptic integral of the second kind E(m) = ∫₀^{π/2} √(1 − m sin²θ) dθ, parameter m = k²
 * (same convention as ellipticK in core/numerics.mjs), m ≤ 1. Arithmetic–geometric mean with
 * E = K(1 − Σ 2^{n−1} c_n²) (A&S 17.6). Accuracy ~1e-15. E(1) = 1. Throws for m > 1.
 */
export function ellipticE(m) {
  finite('ellipticE', { m });
  need(m <= 1, `ellipticE: parameter m = k² must be ≤ 1 (got ${m})`);
  if (m === 1) return 1;
  let a = 1, b = Math.sqrt(1 - m), c = 0, sum = m / 2, pow = 0.5;
  for (let i = 0; i < 60; i++) {
    const an = (a + b) / 2;
    c = (i === 0 ? m : c * c) / (4 * an);                       // c_{n+1} = c_n²/(4a_{n+1}): no cancellation
    b = Math.sqrt(a * b); a = an;
    pow *= 2;
    sum += pow * c * c;
    if (Math.abs(c) <= EPS * a) break;
  }
  return Math.PI / (2 * a) * (1 - sum);
}

/** sn, cn, dn, am for 0 < m < 1 by the descending Landen / AGM scheme (A&S 16.4). */
function jacobiCore(u, m) {
  const a = [1], c = [Math.sqrt(m)];
  let b = Math.sqrt(1 - m), n = 0;
  while (c[n] > EPS * a[n] && n < 60) {
    const an = (a[n] + b) / 2;
    c.push(c[n] * c[n] / (4 * an));
    b = Math.sqrt(a[n] * b);
    a.push(an);
    n++;
  }
  let phi = Math.pow(2, n) * a[n] * u;
  for (let k = n; k > 0; k--) phi = (phi + Math.asin(c[k] * Math.sin(phi) / a[k])) / 2;
  const sn = Math.sin(phi), cn = Math.cos(phi);
  return { sn, cn, dn: Math.sqrt(cn * cn + (1 - m) * sn * sn), am: phi };   // dn² = cn² + m′sn²: no cancellation
}

/**
 * Jacobi elliptic functions sn(u|m), cn(u|m), dn(u|m) and the amplitude am(u|m) (sn = sin am,
 * cn = cos am), parameter m = k² (A&S ch. 16 convention), any real u and m.
 * 0 < m < 1: descending Landen transformation (A&S 16.4); m = 0, 1: circular / hyperbolic limits;
 * m > 1 and m < 0: reciprocal- and negative-parameter transformations (A&S 16.11, 16.10).
 * am is continuous and increasing for m < 1 (am(u + 2K) = am(u) + π); for m > 1 it is asin(sn).
 * Accuracy ~1e-15 absolute for sn, cn, dn with |u| ≲ 1e3·K (phase error grows ∝ |u|·eps).
 * Returns { sn, cn, dn, am }.
 */
export function jacobiElliptic(u, m) {
  finite('jacobiElliptic', { u, m });
  if (m === 0) return { sn: Math.sin(u), cn: Math.cos(u), dn: 1, am: u };
  if (m === 1) {
    const t = Math.tanh(u), sech = 1 / Math.cosh(u);
    return { sn: t, cn: sech, dn: sech, am: 2 * Math.atan(Math.tanh(u / 2)) };
  }
  if (m > 1) {
    const mu = Math.sqrt(m), r = jacobiCore(u * mu, 1 / m);
    const sn = r.sn / mu;
    return { sn, cn: r.dn, dn: r.cn, am: Math.asin(sn) };
  }
  if (m < 0) {
    const M = -m, mu = M / (1 + M), v = u * Math.sqrt(1 + M), r = jacobiCore(v, mu);
    const sn = r.sn / (Math.sqrt(1 + M) * r.dn), cn = r.cn / r.dn, dn = 1 / r.dn;
    const base = Math.atan2(sn, cn), K = ellipticK(m);
    const am = base + 2 * Math.PI * Math.round((Math.PI * u / (2 * K) - base) / (2 * Math.PI));
    return { sn, cn, dn, am };
  }
  return jacobiCore(u, m);
}

// ================================================================= Γ-related

/**
 * Digamma ψ(x) = Γ′(x)/Γ(x), real x ≠ 0, −1, −2, … (throws at the poles). Reflection
 * ψ(1−x) − ψ(x) = π cot πx for x < 0, upward recurrence to x ≥ 12, then the asymptotic series
 * (DLMF 5.11.2). Accuracy ~1e-15 absolute·(1 + |ln x|); relative near the positive root 1.4616…
 * is limited to the absolute error.
 */
export function digamma(x) {
  finite('digamma', { x });
  need(!isNonPosInt(x), `digamma: pole at x = ${x}`);
  if (x < 0) return digamma(1 - x) - Math.PI * cotpi(x);
  let shift = 0, y = x;
  while (y < 12) { shift -= 1 / y; y += 1; }
  const r = 1 / (y * y);
  const series = r * (1 / 12 - r * (1 / 120 - r * (1 / 252 - r * (1 / 240 - r * (1 / 132 - r * (691 / 32760 - r / 12))))));
  return shift + Math.log(y) - 0.5 / y - series;
}

const STIRLING = [1 / 12, -1 / 360, 1 / 1260, -1 / 1680, 1 / 1188, -691 / 360360, 1 / 156];
/** ln Γ(a + b) − ln Γ(a) for a ≥ 15, b > 0, without cancelling two large ln Γ values (Stirling). */
function lgammaDiff(a, b) {
  const ab = a + b;
  let corr = 0, pa = 1 / a, pab = 1 / ab;
  const ra = pa * pa, rab = pab * pab;
  for (const c of STIRLING) { corr += c * (pab - pa); pa *= ra; pab *= rab; }
  return (a - 0.5) * Math.log1p(b / a) + b * Math.log(ab) - b + corr;
}

/**
 * ln B(a, b) = ln Γ(a) + ln Γ(b) − ln Γ(a + b), a, b > 0. When max(a, b) ≥ 15 the difference
 * ln Γ(a+b) − ln Γ(max) is formed directly from Stirling's series (no cancellation), so a ≫ b
 * (e.g. Student t with ν ~ 1e9) stays accurate. Absolute accuracy ~eps·(1 + |ln Γ(min(a, b))|).
 */
export function logBeta(a, b) {
  finite('logBeta', { a, b });
  need(a > 0 && b > 0, `logBeta: need a, b > 0 (got ${a}, ${b})`);
  const hi = Math.max(a, b), lo = Math.min(a, b);
  if (hi < 15) return lgamma(a) + lgamma(b) - lgamma(a + b);
  return lgamma(lo) - lgammaDiff(hi, lo);
}

/**
 * Beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b), a, b > 0. Direct Γ ratio when a + b < 170 (~1e-14
 * relative), otherwise exp(logBeta) (relative error ~eps·|ln B|: ~1e-13 for |ln B| ~ 1e3).
 */
export function beta(a, b) {
  finite('beta', { a, b });
  need(a > 0 && b > 0, `beta: need a, b > 0 (got ${a}, ${b})`);
  if (a + b < 170) return (gamma(a) / gamma(a + b)) * gamma(b);
  return Math.exp(logBeta(a, b));
}

function betaCF(a, b, x, xc) {
  const tiny = 1e-300, qab = a + b, qap = a + 1, qam = a - 1;
  let c = 1, d = x > 0.5 ? (1 - b + qab * xc) / qap : 1 - qab * x / qap;   // 1 − (a+b)x/(a+1) without cancellation
  if (Math.abs(d) < tiny) d = tiny;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= 100000; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c; if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d; h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d; if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c; if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-16) return h;
  }
  throw new Error('incompleteBeta: continued fraction did not converge');
}

/** I_x(a, b) with the complement xc = 1 − x supplied exactly. */
function incBeta(x, xc, a, b) {
  if (x <= 0) return 0;
  if (xc <= 0) return 1;
  const lx = x > 0.5 ? Math.log1p(-xc) : Math.log(x), lxc = xc > 0.5 ? Math.log1p(-x) : Math.log(xc);
  const front = Math.exp(-logBeta(a, b) + a * lx + b * lxc);
  if (x < (a + 1) / (a + b + 2)) return front * betaCF(a, b, x, xc) / a;
  return 1 - front * betaCF(b, a, xc, x) / b;
}

/**
 * Regularised incomplete beta function I_x(a, b) = B(x; a, b)/B(a, b), 0 ≤ x ≤ 1, a, b > 0.
 * Lentz continued fraction with the symmetry I_x(a,b) = 1 − I_{1−x}(b,a) (Numerical Recipes §6.4).
 * Accuracy ~1e-14 for a, b ≲ 100; when both a and b are large (≳ 1e4) the prefactor x^a(1−x)^b/B
 * cancels and the error grows ~eps·min(a, b) (≈1e-12 at 1e4). a ≫ b (large-ν Student t) stays accurate.
 * Binomial tails: P(X ≥ k) = I_p(k, n−k+1); F-distribution and Student-t p-values (see studentTPValue).
 */
export function incompleteBeta(x, a, b) {
  finite('incompleteBeta', { x, a, b });
  need(a > 0 && b > 0, `incompleteBeta: need a, b > 0 (got ${a}, ${b})`);
  need(x >= 0 && x <= 1, `incompleteBeta: need 0 ≤ x ≤ 1 (got ${x})`);
  return incBeta(x, 1 - x, a, b);
}

/** Two-sided tail I_{ν/(ν+t²)}(ν/2, ½) = P(|T| ≥ |t|). */
function studentTail(t, nu) {
  const logRatio = 2 * Math.log(Math.abs(t)) - Math.log(nu);
  if (logRatio > Math.log(Number.MAX_VALUE)) {
    // x = ν/(ν+t²) underflows but x^(ν/2) need not (e.g. the Cauchy 1/t tail).
    const a = nu / 2, logX = -logRatio;
    return Math.exp(a * logX - Math.log(a) - logBeta(a, 0.5));
  }
  const ratio = Math.exp(logRatio);
  return incBeta(1 / (1 + ratio), ratio / (1 + ratio), nu / 2, 0.5);
}

/**
 * Cumulative distribution of Student's t with ν > 0 degrees of freedom (ν need not be an integer):
 * P(T ≤ t). Accuracy ~1e-14 absolute (the small tail itself is accurate relatively: use studentTPValue).
 */
export function studentTCdf(t, nu) {
  need(!Number.isNaN(t), 'studentTCdf: t must be a number');
  finite('studentTCdf', { nu });
  need(nu > 0, `studentTCdf: ν must be > 0 (got ${nu})`);
  if (!Number.isFinite(t)) return t > 0 ? 1 : 0;
  const tail = studentTail(t, nu) / 2;
  return t > 0 ? 1 - tail : tail;
}

/**
 * p-value of a Student-t statistic with ν > 0 degrees of freedom. tails = 2: P(|T| ≥ |t|);
 * tails = 1: upper tail P(T ≥ t). Small p-values keep full relative accuracy (~1e-13) for ν ≲ 1e6;
 * the absolute error grows beyond (~1e-10 at ν = 1e8): use the normal limit there.
 */
export function studentTPValue(t, nu, tails = 2) {
  need(!Number.isNaN(t), 'studentTPValue: t must be a number');
  finite('studentTPValue', { nu });
  need(nu > 0, `studentTPValue: ν must be > 0 (got ${nu})`);
  need(tails === 1 || tails === 2, 'studentTPValue: tails must be 1 or 2');
  if (!Number.isFinite(t)) return tails === 2 ? 0 : (t > 0 ? 0 : 1);
  const two = studentTail(t, nu);
  if (tails === 2) return two;
  return t >= 0 ? two / 2 : 1 - two / 2;
}

// ================================================================= ζ and polylogarithms

const BORWEIN = (() => {                                           // Borwein (2000), algorithm 2, n = 30
  const n = 30, d = [1];
  let T = 1;
  for (let i = 1; i <= n; i++) { T *= 4 * (n + i - 1) * (n - i + 1) / ((2 * i) * (2 * i - 1)); d.push(d[i - 1] + T); }
  return d.slice(0, n).map(dk => (d[n] - dk) / d[n]);
})();

/** Dirichlet η(s) = Σ (−1)^{k−1} k^{−s} for s > 0 (Borwein acceleration, error < 1e-20). */
function etaPositive(s) {
  let sum = 0;
  for (let k = 0; k < BORWEIN.length; k++) sum += (k % 2 ? -1 : 1) * BORWEIN[k] * Math.pow(k + 1, -s);
  return sum;
}

/**
 * Riemann zeta function ζ(s) for real s ≠ 1 (throws at the pole). s > 0: Borwein's accelerated
 * alternating series for η(s) and ζ = η/(1 − 2^{1−s}); s < 0: the functional equation
 * (log-scaled, so large |s| neither overflows nor underflows early); exact zeros at negative even
 * integers. Accuracy ~1e-15 relative for s > 0, ~1e-14·(1 + ln|s|) for s < 0.
 * ζ(2) = π²/6, ζ(3) = 1.2020569…, ζ(3/2) = 2.6123753… (Bose–Einstein condensation).
 */
export function zeta(s) {
  finite('zeta', { s });
  need(s !== 1, 'zeta: pole at s = 1');
  if (s === 0) return -0.5;
  // The reflection formula multiplies sin(πs/2) by the pole ζ(1−s), which loses
  // all relative accuracy when 1−s rounds to 1. Use the analytic expansion at zero.
  if (Math.abs(s) < 1e-7) return -0.5 - 0.5 * Math.log(2 * Math.PI) * s - 1.0031782279542925 * s * s;
  if (s > 0) return etaPositive(s) / -Math.expm1((1 - s) * Math.LN2);
  if (Number.isInteger(s) && s % 2 === 0) return 0;
  const sn = sinpi(s / 2), z1 = zeta(1 - s);
  if (1 - s < 170) return Math.pow(2, s) * Math.pow(Math.PI, s - 1) * sn * gamma(1 - s) * z1;
  return Math.sign(sn) * Math.exp(s * Math.LN2 + (s - 1) * Math.log(Math.PI) + Math.log(Math.abs(sn)) + lgamma(1 - s) + Math.log(z1));
}

/** η(s) = (1 − 2^{1−s}) ζ(s), all real s (η(1) = ln 2). */
function eta(s) {
  if (s === 1) return Math.LN2;
  if (s > 0) return etaPositive(s);
  return -Math.expm1((1 - s) * Math.LN2) * zeta(s);
}

/** Σ_{k≥1} z^k k^{−s}; null if not converged or cancelled. */
function polylogSeries(s, z, maxTerms) {
  let sum = 0, big = 0, p = 1, small = 0;
  for (let k = 1; k <= maxTerms; k++) {
    p *= z;
    const term = p * Math.pow(k, -s);
    sum += term; big = Math.max(big, Math.abs(term));
    small = Math.abs(term) <= 1e-17 * Math.abs(sum) && k > -s ? small + 1 : 0;
    if (small >= 2) return big <= CANCEL * Math.abs(sum) ? sum : null;
  }
  return null;
}

/** Li_s(e^μ) for −0.7 < μ < 0 (DLMF 25.12.12 and its integer-s limit). */
function polylogLog(s, mu) {
  if (!Number.isInteger(s)) {
    // Γ(1−s)(−μ)^(s−1) and ζ(s−k) contain cancelling poles near integer s.
    // Reject this expansion before parameter-rounding error contaminates the result.
    const distance = Math.abs(s - Math.round(s));
    need(distance >= 1e-4, `polylog: order s = ${s} too close to an integer for the logarithmic expansion; unsupported`);
  }
  let sum, big, start;
  if (Number.isInteger(s) && s >= 2) {
    let H = 0, fact = 1;
    for (let k = 1; k < s; k++) { H += 1 / k; fact *= k; }
    sum = Math.pow(mu, s - 1) / fact * (H - Math.log(-mu));
    start = 0;
  } else {
    sum = gamma(1 - s) * Math.pow(-mu, s - 1);
    start = 0;
  }
  big = Math.abs(sum);
  let p = 1, small = 0;
  for (let k = start; k < 200; k++) {
    if (k > 0) p *= mu / k;
    if (Number.isInteger(s) && s >= 2 && k === s - 1) continue;
    const term = zeta(s - k) * p;
    sum += term; big = Math.max(big, Math.abs(term));
    small = Math.abs(term) <= 1e-17 * Math.abs(sum) && k > s + 2 ? small + 1 : 0;
    if (small >= 3) break;
  }
  need(big <= CANCEL * Math.abs(sum), `polylog: cancellation near integer order s = ${s}; unsupported`);
  return sum;
}

/**
 * Polylogarithm Li_s(z) = Σ_{k≥1} z^k/k^s for real order s and real z ≤ 1 (Li_s(1) = ζ(s), s > 1).
 * Paths: power series (|z| ≤ ½); the expansion in μ = ln z about z = 1 (DLMF 25.12.12) for ½ < z < 1;
 * duplication Li_s(z) + Li_s(−z) = 2^{1−s} Li_s(z²) for −1 < z < −½ with s > 1, the alternating series
 * there for s ≤ 1; Li_s(−1) = −η(s); inversion formulae for z < −1 with integer s (DLMF 25.12.4-type).
 * Li₁(z) = −ln(1 − z). Bose–Einstein integrals: g_s(z) = Li_s(z); Fermi–Dirac: f_s(z) = −Li_s(−z).
 * Accuracy ~1e-14 typically (each path checks its cancellation); orders within ~1e-4 of an integer
 * in the μ-expansion lose digits and throw. Throws for z > 1 (complex), z = 1 with s ≤ 1, and z < −1
 * with non-integer s (unsupported).
 */
export function polylog(s, z) {
  const who = 'polylog';
  finite(who, { s, z });
  if (z === 0) return 0;
  need(z <= 1, `${who}: z = ${z} > 1 gives a complex value`);
  if (s === 1) return -Math.log1p(-z);
  if (z === 1) { need(s > 1, `${who}: Li_s(1) = ζ(s) diverges for s ≤ 1`); return zeta(s); }
  if (z === -1) return -eta(s);
  if (z < -1) {
    need(Number.isInteger(s), `${who}: z < −1 is supported only for integer order s (got s = ${s})`);
    if (s === 0) return z / (1 - z);
    if (s < 0) return -((-s) % 2 ? -1 : 1) * polylog(s, 1 / z);
    const L = Math.log(-z);
    let rest = 0, fact = 1;
    for (let k = 1; k <= s; k++) fact *= k;
    rest -= Math.pow(L, s) / fact;
    for (let k = 1; 2 * k <= s; k++) {
      let f = 1;
      for (let i = 2; i <= s - 2 * k; i++) f *= i;
      rest -= 2 * Math.pow(L, s - 2 * k) / f * eta(2 * k);
    }
    return -(s % 2 ? -1 : 1) * polylog(s, 1 / z) + rest;
  }
  if (Math.abs(z) <= 0.5) {
    const v = polylogSeries(s, z, 10000);
    need(v !== null, `${who}: series cancellation for s = ${s}, z = ${z}; unsupported`);
    return v;
  }
  if (z > 0) return polylogLog(s, Math.log(z));
  if (s > 1) return Math.pow(2, 1 - s) * polylog(s, z * z) - polylog(s, -z);
  const v = polylogSeries(s, z, 400000);
  need(v !== null, `${who}: z = ${z} too close to −1 for order s = ${s} ≤ 1; unsupported`);
  return v;
}

// ================================================================= exponential integrals

/**
 * Generalised exponential integral E_n(x) = ∫₁^∞ e^{−xt} t^{−n} dt, integer n ≥ 0, x ≥ 0
 * (x = 0 needs n ≥ 2: E_n(0) = 1/(n−1)). Series (x ≤ 1) or Lentz continued fraction (x > 1),
 * Numerical Recipes §6.3. Accuracy ~1e-15 relative; underflows to 0 for x ≳ 700.
 */
export function expintEn(n, x) {
  intArg('expintEn', 'n', n); finite('expintEn', { x });
  need(x >= 0, `expintEn: x must be ≥ 0 (got ${x}); E_n is complex for x < 0`);
  if (n === 0) return x === 0 ? Infinity : Math.exp(-x) / x;
  if (x === 0) { need(n >= 2, 'expintEn: E_1(0) diverges'); return 1 / (n - 1); }
  if (x > 1) {
    let b = x + n, c = 1e300, d = 1 / b, h = d;
    for (let i = 1; i < 10000; i++) {
      const an = -i * (n - 1 + i);
      b += 2;
      d = 1 / (an * d + b);
      c = b + an / c;
      const del = c * d;
      h *= del;
      if (Math.abs(del - 1) < 1e-16) return h * Math.exp(-x);
    }
    throw new Error('expintEn: continued fraction did not converge');
  }
  const nm1 = n - 1;
  let ans = nm1 !== 0 ? 1 / nm1 : -Math.log(x) - EULER, fact = 1;
  for (let i = 1; i < 1000; i++) {
    fact *= -x / i;
    let del;
    if (i !== nm1) del = -fact / (i - nm1);
    else {
      let psi = -EULER;
      for (let k = 1; k <= nm1; k++) psi += 1 / k;
      del = fact * (-Math.log(x) + psi);
    }
    ans += del;
    if (Math.abs(del) < Math.abs(ans) * 1e-17) return ans;
  }
  throw new Error('expintEn: series did not converge');
}

/** Exponential integral E₁(x) = ∫ₓ^∞ e^{−t}/t dt, x > 0 (E₁(0) = Infinity). See expintEn. */
export function expintE1(x) {
  finite('expintE1', { x });
  need(x >= 0, `expintE1: x must be ≥ 0 (got ${x}); for x < 0 use −expintEi(−x)`);
  return x === 0 ? Infinity : expintEn(1, x);
}

/**
 * Exponential integral Ei(x) = −PV∫_{−x}^∞ e^{−t}/t dt, real x ≠ 0 (Ei(0) = −Infinity; Ei(x) = −E₁(−x)
 * for x < 0). Positive x: power series γ + ln x + Σ x^k/(k·k!) for x ≤ 40, asymptotic e^x/x Σ k!/x^k
 * beyond. Accuracy ~1e-15 relative (absolute ~1e-16 near the root x₀ = 0.37250741…); overflows to
 * Infinity for x ≳ 710.
 */
export function expintEi(x) {
  finite('expintEi', { x });
  if (x === 0) return -Infinity;
  if (x < 0) return -expintE1(-x);
  if (x <= 40) {
    let term = 1, sum = 0;
    for (let k = 1; k < 500; k++) {
      term *= x / k;
      const add = term / k;
      sum += add;
      if (add < 1e-17 * sum) break;
    }
    return EULER + Math.log(x) + sum;
  }
  let term = 1, sum = 1;
  for (let k = 1; k < 100; k++) {
    const next = term * k / x;
    if (next > term) break;
    term = next; sum += term;
    if (term < 1e-17 * sum) break;
  }
  const half = Math.exp(x / 2);
  return ((half / x) * sum) * half;
}

// ================================================================= Dawson's integral

const DAWSON_H = 0.2, DAWSON_C = Array.from({ length: 18 }, (_, i) => Math.exp(-(((2 * i + 1) * DAWSON_H) ** 2)));

/**
 * Dawson's integral F(x) = e^{−x²} ∫₀ˣ e^{t²} dt (= (√π/2) e^{−x²} erfi x), any real x (odd).
 * Power series for |x| < 1, Rybicki's exponentially convergent sampling formula (h = 0.2,
 * truncation error ~e^{−(π/2h)²} ≈ 1e-27) for 1 ≤ |x| ≤ 50, asymptotic series 1/(2x) Σ (2k−1)!!/(2x²)^k
 * beyond. Accuracy ~1e-15 relative. Maximum F = 0.5410442246… at x = 0.9241388730….
 * Appears in plasma dispersion (Re Z), Voigt profiles and ion-acoustic responses.
 */
export function dawson(x) {
  finite('dawson', { x });
  const ax = Math.abs(x), sign = x < 0 ? -1 : 1;
  if (ax < 1) {
    let term = ax, sum = ax;
    const q = 2 * ax * ax;
    for (let k = 0; k < 100; k++) {
      term *= -q / (2 * k + 3);
      sum += term;
      if (Math.abs(term) <= 1e-17 * Math.abs(sum)) break;
    }
    return sign * sum;
  }
  if (ax > 50) {
    const r = 1 / (2 * ax * ax);
    let term = 1, sum = 1;
    for (let k = 1; k < 40; k++) { term *= (2 * k - 1) * r; sum += term; if (term < 1e-17 * sum) break; }
    return sign * sum / (2 * ax);
  }
  const n0 = 2 * Math.round(ax / (2 * DAWSON_H)), xp = ax - n0 * DAWSON_H;
  let e1 = Math.exp(2 * xp * DAWSON_H);
  const e2 = e1 * e1;
  let d1 = n0 + 1, d2 = d1 - 2, sum = 0;
  for (const c of DAWSON_C) {
    sum += c * (e1 / d1 + 1 / (d2 * e1));
    d1 += 2; d2 -= 2; e1 *= e2;
  }
  return sign * Math.exp(-xp * xp) * sum / SQRT_PI;
}

// ================================================================= Wigner symbols (exact)

const FACT = [1n];
function factBig(n) {
  while (FACT.length <= n) FACT.push(FACT[FACT.length - 1] * BigInt(FACT.length));
  return FACT[n];
}
const absBig = v => (v < 0n ? -v : v);
function gcdBig(u, v) {
  let a = absBig(u), b = absBig(v);
  while (b) [a, b] = [b, a % b];
  return a;
}
const bitLength = v => (v === 0n ? 0 : v.toString(2).length);
/** A/B (positive BigInts), rounded to binary64 with ties to even, including subnormals. */
function ratioToNumber(A, B) {
  if (A === 0n) return 0;
  let exponent = bitLength(A) - bitLength(B);
  if (exponent >= 0 ? A < (B << BigInt(exponent)) : (A << BigInt(-exponent)) < B) exponent--;
  if (exponent > 1023) return Infinity;
  if (exponent < -1075) return 0;
  const unit = Math.max(exponent - 52, -1074);
  const numerator = unit < 0 ? A << BigInt(-unit) : A;
  const denominator = unit > 0 ? B << BigInt(unit) : B;
  let q = numerator / denominator;
  const remainder = numerator % denominator, twiceRemainder = 2n * remainder;
  if (twiceRemainder > denominator || (twiceRemainder === denominator && (q & 1n))) q++;
  return Number(q) * 2 ** unit;
}
/** Doubled integer 2j for an integer or half-integer j. */
function twice(who, name, j) {
  need(Number.isFinite(j) && Number.isInteger(2 * j), `${who}: ${name} must be an integer or half-integer (got ${j})`);
  return 2 * j;
}
/** Triangle condition on doubled values: |a−b| ≤ c ≤ a+b and a+b+c even. */
const triangle = (a, b, c) => c >= Math.abs(a - b) && c <= a + b && (a + b + c) % 2 === 0;
/** Δ(abc) on doubled values as [numerator, denominator]. */
const delta = (a, b, c) => [factBig((a + b - c) / 2) * factBig((a - b + c) / 2) * factBig((-a + b + c) / 2), factBig((a + b + c) / 2 + 1)];

/** Exact square of a symbol: value = sign·√(num/den), reduced. */
function exactResult(sign, num, den) {
  if (num === 0n || sign === 0) return { sign: 0, num: 0n, den: 1n };
  const g = gcdBig(num, den);
  return { sign, num: num / g, den: den / g };
}
const toValue = r => (r.sign === 0 ? 0 : r.sign * Math.sqrt(ratioToNumber(r.num, r.den)));

function threeJExact(J1, J2, J3, M1, M2, M3) {
  const ok = M1 + M2 + M3 === 0 && [[J1, M1], [J2, M2], [J3, M3]].every(([j, m]) => Math.abs(m) <= j && (j + m) % 2 === 0) && triangle(J1, J2, J3);
  if (!ok) return { sign: 0, num: 0n, den: 1n };
  const j1 = J1 / 2, j2 = J2 / 2, j3 = J3 / 2, m1 = M1 / 2, m2 = M2 / 2, m3 = M3 / 2;
  const kmin = Math.max(0, j2 - j3 - m1, j1 - j3 + m2), kmax = Math.min(j1 + j2 - j3, j1 - m1, j2 + m2);
  let sn = 0n, sd = 1n;                                          // S = sn/sd
  for (let k = kmin; k <= kmax; k++) {
    const d = factBig(k) * factBig(j1 + j2 - j3 - k) * factBig(j1 - m1 - k) * factBig(j2 + m2 - k) * factBig(j3 - j2 + m1 + k) * factBig(j3 - j1 - m2 + k);
    sn = k % 2 ? sn * d - sd : sn * d + sd;
    sd *= d;
    const g = gcdBig(sn, sd);
    if (g > 1n) { sn /= g; sd /= g; }
  }
  if (sn === 0n) return { sign: 0, num: 0n, den: 1n };
  const [dn, dd] = delta(J1, J2, J3);
  const F = factBig(j1 + m1) * factBig(j1 - m1) * factBig(j2 + m2) * factBig(j2 - m2) * factBig(j3 + m3) * factBig(j3 - m3);
  const phase = (((j1 - j2 - m3) % 2) + 2) % 2 ? -1 : 1;
  return exactResult(phase * (sn < 0n ? -1 : 1), dn * F * sn * sn, dd * sd * sd);
}

/**
 * Wigner 3j symbol (j1 j2 j3; m1 m2 m3) for integer or half-integer arguments, by Racah's formula
 * in exact rational (BigInt) arithmetic: the only rounding is the final square root (~1 ulp).
 * Selection rules (m1+m2+m3 = 0, triangle, |m| ≤ j, j − m integer) give exactly 0. Practical for
 * j up to a few hundred. With { exact: true } returns { sign, num, den } (BigInt): value = sign·√(num/den).
 * Throws for arguments that are not integers or half-integers.
 */
export function wigner3j(j1, j2, j3, m1, m2, m3, { exact = false } = {}) {
  const w = 'wigner3j';
  const r = threeJExact(twice(w, 'j1', j1), twice(w, 'j2', j2), twice(w, 'j3', j3), twice(w, 'm1', m1), twice(w, 'm2', m2), twice(w, 'm3', m3));
  return exact ? r : toValue(r);
}

/**
 * Clebsch–Gordan coefficient ⟨j1 m1; j2 m2 | J M⟩ (Condon–Shortley phase) =
 * (−1)^{j1−j2+M} √(2J+1) (j1 j2 J; m1 m2 −M). Exact arithmetic as wigner3j; { exact: true } likewise.
 */
export function clebschGordan(j1, m1, j2, m2, J, M, { exact = false } = {}) {
  const w = 'clebschGordan';
  const J1 = twice(w, 'j1', j1), J2 = twice(w, 'j2', j2), JJ = twice(w, 'J', J);
  const r = threeJExact(J1, J2, JJ, twice(w, 'm1', m1), twice(w, 'm2', m2), -twice(w, 'M', M));
  if (r.sign === 0) return exact ? r : 0;
  const phase = ((((J1 - J2 + twice(w, 'M', M)) / 2) % 2) + 2) % 2 ? -1 : 1;
  const out = exactResult(phase * r.sign, r.num * BigInt(JJ + 1), r.den);
  return exact ? out : toValue(out);
}

/**
 * Wigner 6j symbol {j1 j2 j3; j4 j5 j6} by Racah's formula in exact rational arithmetic
 * (recoupling of three angular momenta). Zero unless the triads (j1 j2 j3), (j1 j5 j6), (j4 j2 j6),
 * (j4 j5 j3) satisfy the triangle rule with integer sums. { exact: true } as for wigner3j.
 */
export function wigner6j(j1, j2, j3, j4, j5, j6, { exact = false } = {}) {
  const w = 'wigner6j';
  const [a, b, c, d, e, f] = [j1, j2, j3, j4, j5, j6].map((j, i) => twice(w, `j${i + 1}`, j));
  const zero = { sign: 0, num: 0n, den: 1n };
  if (![[a, b, c], [a, e, f], [d, b, f], [d, e, c]].every(([x, y, z]) => x >= 0 && y >= 0 && z >= 0 && triangle(x, y, z))) return exact ? zero : 0;
  const A = [a + b + c, a + e + f, d + b + f, d + e + c].map(v => v / 2);
  const B = [a + b + d + e, b + c + e + f, c + a + f + d].map(v => v / 2);
  let sn = 0n, sd = 1n;
  for (let t = Math.max(...A); t <= Math.min(...B); t++) {
    const den = A.reduce((p, x) => p * factBig(t - x), 1n) * B.reduce((p, x) => p * factBig(x - t), 1n);
    const num = factBig(t + 1);
    sn = t % 2 ? sn * den - num * sd : sn * den + num * sd;
    sd *= den;
    const g = gcdBig(sn, sd);
    if (g > 1n) { sn /= g; sd /= g; }
  }
  if (sn === 0n) return exact ? zero : 0;
  let pn = 1n, pd = 1n;
  for (const [x, y, z] of [[a, b, c], [a, e, f], [d, b, f], [d, e, c]]) { const [n, m] = delta(x, y, z); pn *= n; pd *= m; }
  const r = exactResult(sn < 0n ? -1 : 1, pn * sn * sn, pd * sd * sd);
  return exact ? r : toValue(r);
}
