// Spectral methods and nonlinear waves (dimensionless model equations; supply your own units):
//   Chebyshev collocation on [a, b]   chebyshevNodes, chebyshevDiffMatrix, chebyshevInterpolant,
//                                    chebyshevCoefficients, clenshawCurtisWeights, clenshawCurtis,
//                                    chebyshevBVP, chebyshevEigen (+ eigenNonsymmetric: dense real QR)
//   Fourier pseudo-spectral (periodic, n a power of two)   fourierDerivative
//   ETDRK4 (Kassam & Trefethen, SIAM J. Sci. Comput. 26, 1214 (2005))   etdrk4Periodic, kdv,
//                                    kuramotoSivashinsky, viscousBurgers
//   Split-step Fourier (Strang, optional Yoshida 4th order)   nls, sineGordon
//   Exact solutions for verification   nlsSoliton, nlsPeregrine, kdvSoliton, kdvTwoSoliton,
//                                    kdvPhaseShifts, burgersColeHopf, sineGordonKink, sineGordonBreather
// Chebyshev nodes are Chebyshev–Lobatto points in ASCENDING order: x_j = −cos(jπ/n), j = 0…n.
// Fourier grids: x_j = x₀ + j·dx, j = 0…n−1, period n·dx; FFT conventions as numerics.mjs (fft, fftFreq).
// Dense matrices are arrays of Float64Array rows (linalg.mjs convention).

import { fft, fftFreq, gaussLegendreNodes } from './numerics.mjs';
import { solveLinear } from './linalg.mjs';

const EPS = Number.EPSILON;
const fail = (who, msg) => { throw new Error(`${who}: ${msg}`); };
const isPow2 = n => Number.isInteger(n) && n >= 2 && (n & (n - 1)) === 0;
const finiteArray = (v, who, name = 'values') => {
  if (!v || !(v.length > 0)) fail(who, `${name} must be a non-empty array`);
  for (let i = 0; i < v.length; i++) if (!Number.isFinite(v[i])) fail(who, `${name} must be finite`);
};
const interval = (iv, who) => {
  const [a, b] = iv ?? [-1, 1];
  if (!(Number.isFinite(a) && Number.isFinite(b) && b > a)) fail(who, 'interval must be [a, b] with finite a < b');
  return [a, b];
};
const order = (n, who, min = 1) => { if (!(Number.isInteger(n) && n >= min)) fail(who, `n must be an integer ≥ ${min}`); };
const coef = (v, x) => (typeof v === 'function' ? v(x) : v ?? 0);

// ------------------------------------------------------------------ Chebyshev collocation

/** n + 1 Chebyshev–Lobatto points on [a, b], ascending: x_j = a + (b − a)(1 − cos(jπ/n))/2 (symmetric to round-off). */
export function chebyshevNodes(n, iv = [-1, 1]) {
  order(n, 'chebyshevNodes');
  const [a, b] = interval(iv, 'chebyshevNodes');
  return Float64Array.from({ length: n + 1 }, (_, j) => a + (b - a) * (1 + Math.sin(Math.PI * (2 * j - n) / (2 * n))) / 2);
}

/**
 * Chebyshev differentiation matrix (Trefethen, Spectral Methods in MATLAB, cheb.m, with the
 * negative-sum trick for the diagonal) on n + 1 ascending nodes of [a, b]. (D·u)_i is the exact
 * derivative at x_i of the degree-n interpolant of u. D2 = D·D. Entries grow like n²/(b−a) and n⁴:
 * keep n ≲ 200 (round-off in D2 grows ~n⁴·ε). Returns { x, D, D2 }.
 */
export function chebyshevDiffMatrix(n, iv = [-1, 1]) {
  order(n, 'chebyshevDiffMatrix');
  const [a, b] = interval(iv, 'chebyshevDiffMatrix'), scale = 2 / (b - a);
  const s = Float64Array.from({ length: n + 1 }, (_, j) => Math.sin(Math.PI * (2 * j - n) / (2 * n)));   // on [−1, 1]
  const c = Float64Array.from({ length: n + 1 }, (_, j) => (j === 0 || j === n ? 2 : 1) * (j % 2 ? -1 : 1));
  const D = Array.from({ length: n + 1 }, () => new Float64Array(n + 1));
  for (let i = 0; i <= n; i++) {
    let sum = 0;
    for (let j = 0; j <= n; j++) {
      if (i === j) continue;
      const v = (c[i] / c[j]) / (s[i] - s[j]) * scale;
      D[i][j] = v; sum += v;
    }
    D[i][i] = -sum;
  }
  const D2 = Array.from({ length: n + 1 }, () => new Float64Array(n + 1));
  for (let i = 0; i <= n; i++) for (let k = 0; k <= n; k++) { const d = D[i][k]; if (d) for (let j = 0; j <= n; j++) D2[i][j] += d * D[k][j]; }
  return { x: chebyshevNodes(n, [a, b]), D, D2 };
}

/**
 * Barycentric interpolant through values at the n + 1 ascending Chebyshev nodes of [a, b]
 * (weights (−1)^j δ_j, δ = ½ at the ends): numerically stable for any n. Returns x => p(x);
 * exact at the nodes; outside [a, b] it extrapolates (unreliable).
 */
export function chebyshevInterpolant(values, iv = [-1, 1]) {
  finiteArray(values, 'chebyshevInterpolant');
  const n = values.length - 1;
  order(n, 'chebyshevInterpolant');
  const x = chebyshevNodes(n, iv), f = Float64Array.from(values);
  const w = Float64Array.from({ length: n + 1 }, (_, j) => (j % 2 ? -1 : 1) * (j === 0 || j === n ? 0.5 : 1));
  return xq => {
    let num = 0, den = 0;
    for (let j = 0; j <= n; j++) {
      const d = xq - x[j];
      if (d === 0) return f[j];
      const t = w[j] / d; num += t * f[j]; den += t;
    }
    return num / den;
  };
}

/**
 * Chebyshev coefficients a_k (k = 0…n) of the interpolant through values at the ascending nodes:
 * p(x) = Σ a_k T_k(x̂), x̂ = (2x − a − b)/(b − a). O(n²) direct sum (exact DCT-I). Their decay is
 * the resolution diagnostic: a smooth function is resolved when |a_n| ~ ε·max|a_k|.
 */
export function chebyshevCoefficients(values) {
  finiteArray(values, 'chebyshevCoefficients');
  const n = values.length - 1;
  order(n, 'chebyshevCoefficients');
  return Float64Array.from({ length: n + 1 }, (_, k) => {
    let sum = 0;
    for (let j = 0; j <= n; j++) sum += (j === 0 || j === n ? 0.5 : 1) * values[j] * Math.cos(Math.PI * k * (n - j) / n);
    return (k === 0 || k === n ? 1 : 2) * sum / n;
  });
}

/**
 * Clenshaw–Curtis quadrature (Trefethen clencurt.m): nodes = chebyshevNodes(n, [a, b]) and positive
 * weights w with Σ w_j f(x_j) exact for polynomials of degree ≤ n, and geometrically convergent
 * for analytic f. Returns { x, w }.
 */
export function clenshawCurtisWeights(n, iv = [-1, 1]) {
  order(n, 'clenshawCurtisWeights');
  const [a, b] = interval(iv, 'clenshawCurtisWeights'), half = (b - a) / 2;
  const w = new Float64Array(n + 1);
  const end = n % 2 === 0 ? 1 / (n * n - 1) : 1 / (n * n);
  w[0] = w[n] = end * half;
  for (let j = 1; j < n; j++) {
    const th = Math.PI * j / n;
    let v = 1;
    if (n % 2 === 0) {
      for (let k = 1; k < n / 2; k++) v -= 2 * Math.cos(2 * k * th) / (4 * k * k - 1);
      v -= Math.cos(n * th) / (n * n - 1);
    } else {
      for (let k = 1; k <= (n - 1) / 2; k++) v -= 2 * Math.cos(2 * k * th) / (4 * k * k - 1);
    }
    w[j] = 2 * v / n * half;
  }
  return { x: chebyshevNodes(n, [a, b]), w };
}

/** ∫_a^b f(x) dx by (n + 1)-point Clenshaw–Curtis. Compare n and 2n to estimate the error. */
export function clenshawCurtis(f, a, b, n = 64) {
  if (typeof f !== 'function') fail('clenshawCurtis', 'f must be a function');
  const { x, w } = clenshawCurtisWeights(n, [a, b]);
  let sum = 0;
  for (let j = 0; j <= n; j++) sum += w[j] * f(x[j]);
  if (!Number.isFinite(sum)) fail('clenshawCurtis', 'integrand is not finite at a node');
  return sum;
}

const boundaryRow = (spec, D, row, n, who) => {
  if (spec == null || !Number.isFinite(spec.value ?? 0)) fail(who, 'left and right must be { value } or { alpha, beta, value }');
  const alpha = spec.alpha ?? (spec.beta == null ? 1 : 0), beta = spec.beta ?? 0;
  if (!(Number.isFinite(alpha) && Number.isFinite(beta)) || (alpha === 0 && beta === 0)) fail(who, 'boundary needs alpha or beta non-zero');
  return Float64Array.from({ length: n + 1 }, (_, j) => (j === row ? alpha : 0) + beta * D[row][j]);
};

/**
 * Linear two-point boundary-value problem a(x)u'' + b(x)u' + c(x)u = f(x) on [a, b] by Chebyshev
 * collocation (n + 1 nodes, boundary rows replaced). a, b, c, f: numbers or functions of x.
 * left/right: { value } (u = value) or { alpha, beta, value } (αu + βu' = value, Robin/Neumann).
 * Converges geometrically for analytic data; check it by the coefficient tail (`tail` = |a_{n−1}| +
 * |a_n| relative to max|a_k|, should be ≲ 1e−12) or by repeating with larger n.
 * Throws if the collocation matrix is singular (e.g. an eigenvalue of the operator sits at 0).
 * Returns { x, u, interpolant, tail }.
 */
export function chebyshevBVP({ n = 32, interval: iv = [-1, 1], a = 1, b = 0, c = 0, f = 0, left = { value: 0 }, right = { value: 0 } } = {}) {
  order(n, 'chebyshevBVP', 2);
  const { x, D, D2 } = chebyshevDiffMatrix(n, iv);
  const A = [], rhs = new Float64Array(n + 1);
  for (let i = 0; i <= n; i++) {
    const ai = coef(a, x[i]), bi = coef(b, x[i]), ci = coef(c, x[i]);
    if (![ai, bi, ci].every(Number.isFinite)) fail('chebyshevBVP', 'coefficients must be finite on the interval');
    A.push(Float64Array.from({ length: n + 1 }, (_, j) => ai * D2[i][j] + bi * D[i][j] + (i === j ? ci : 0)));
    rhs[i] = coef(f, x[i]);
  }
  A[0] = boundaryRow(left, D, 0, n, 'chebyshevBVP'); rhs[0] = left.value ?? 0;
  A[n] = boundaryRow(right, D, n, n, 'chebyshevBVP'); rhs[n] = right.value ?? 0;
  finiteArray(rhs, 'chebyshevBVP', 'f');
  const u = solveLinear(A, rhs);
  const coeffs = chebyshevCoefficients(u);
  let peak = 0; for (const v of coeffs) peak = Math.max(peak, Math.abs(v));
  return { x, u, interpolant: chebyshevInterpolant(u, iv), tail: peak ? (Math.abs(coeffs[n - 1]) + Math.abs(coeffs[n])) / peak : 0 };
}

/**
 * Eigenvalues of a dense real (non-symmetric) matrix: balancing, Householder-free Gaussian reduction
 * to Hessenberg form and the Francis double-shift QR iteration (Numerical Recipes balanc/elmhes/hqr).
 * O(n³); n ≲ 300 in a browser. Returns [{ re, im }] sorted by re, then im. Throws if an eigenvalue
 * does not converge in 60 iterations.
 */
export function eigenNonsymmetric(A) {
  const n = A?.length;
  if (!(n > 0) || !A.every(r => r.length === n)) fail('eigenNonsymmetric', 'matrix must be square and non-empty');
  let scaleA = 0;
  for (const row of A) {
    for (const v of row) {
      if (!Number.isFinite(v)) fail('eigenNonsymmetric', 'entries must be finite');
      scaleA = Math.max(scaleA, Math.abs(v));
    }
  }
  if (scaleA === 0) return Array.from({ length: n }, () => ({ re: 0, im: 0 }));
  // QR forms products of entries. Normalise first so the same matrix in different units has
  // the same spectrum accuracy, without squaring huge/tiny dimensional values.
  const a = Array.from({ length: n + 1 }, (_, i) => Float64Array.from({ length: n + 1 }, (__, j) => (i && j ? A[i - 1][j - 1] / scaleA : 0)));
  // balance (powers of two: exact)
  for (let done = false; !done;) {
    done = true;
    for (let i = 1; i <= n; i++) {
      let r = 0, c = 0;
      for (let j = 1; j <= n; j++) if (j !== i) { c += Math.abs(a[j][i]); r += Math.abs(a[i][j]); }
      if (!c || !r) continue;
      let g = r / 2, f = 1;
      const s = c + r;
      while (c < g) { f *= 2; c *= 4; }
      g = r * 2;
      while (c > g) { f /= 2; c /= 4; }
      if ((c + r) / f < 0.95 * s) {
        done = false;
        for (let j = 1; j <= n; j++) { a[i][j] /= f; a[j][i] *= f; }
      }
    }
  }
  // reduction to upper Hessenberg form by stabilised elimination
  for (let m = 2; m < n; m++) {
    let x = 0, i = m;
    for (let j = m; j <= n; j++) if (Math.abs(a[j][m - 1]) > Math.abs(x)) { x = a[j][m - 1]; i = j; }
    if (i !== m) {
      for (let j = m - 1; j <= n; j++) { const t = a[i][j]; a[i][j] = a[m][j]; a[m][j] = t; }
      for (let j = 1; j <= n; j++) { const t = a[j][i]; a[j][i] = a[j][m]; a[j][m] = t; }
    }
    if (x) {
      for (let r = m + 1; r <= n; r++) {
        let y = a[r][m - 1];
        if (y) {
          y /= x; a[r][m - 1] = y;
          for (let j = m; j <= n; j++) a[r][j] -= y * a[m][j];
          for (let j = 1; j <= n; j++) a[j][m] += y * a[j][r];
        }
      }
    }
  }
  for (let i = 1; i <= n; i++) for (let j = 1; j < i - 1; j++) a[i][j] = 0;
  // Francis double-shift QR
  const wr = new Float64Array(n + 1), wi = new Float64Array(n + 1);
  let anorm = 0;
  for (let i = 1; i <= n; i++) for (let j = Math.max(i - 1, 1); j <= n; j++) anorm += Math.abs(a[i][j]);
  let nn = n, t = 0;
  let p = 0, q = 0, r = 0, s = 0, w = 0, x = 0, y = 0, z = 0;
  while (nn >= 1) {
    let its = 0, l;
    do {
      for (l = nn; l >= 2; l--) {
        s = Math.abs(a[l - 1][l - 1]) + Math.abs(a[l][l]);
        if (s === 0) s = anorm;
        if (Math.abs(a[l][l - 1]) <= EPS * s) { a[l][l - 1] = 0; break; }
      }
      x = a[nn][nn];
      if (l === nn) { wr[nn] = x + t; wi[nn--] = 0; } else {
        y = a[nn - 1][nn - 1]; w = a[nn][nn - 1] * a[nn - 1][nn];
        if (l === nn - 1) {
          p = 0.5 * (y - x); q = p * p + w; z = Math.sqrt(Math.abs(q)); x += t;
          if (q >= 0) {
            z = p + (p >= 0 ? Math.abs(z) : -Math.abs(z));
            wr[nn - 1] = wr[nn] = x + z;
            if (z) wr[nn] = x - w / z;
            wi[nn - 1] = wi[nn] = 0;
          } else { wr[nn - 1] = wr[nn] = x + p; wi[nn - 1] = -(wi[nn] = z); }
          nn -= 2;
        } else {
          if (its === 60) fail('eigenNonsymmetric', 'QR iteration did not converge');
          if (its === 10 || its === 20 || its === 40) {
            t += x;
            for (let i = 1; i <= nn; i++) a[i][i] -= x;
            s = Math.abs(a[nn][nn - 1]) + Math.abs(a[nn - 1][nn - 2]);
            y = x = 0.75 * s; w = -0.4375 * s * s;
          }
          ++its;
          let m;
          for (m = nn - 2; m >= l; m--) {
            z = a[m][m]; r = x - z; s = y - z;
            p = (r * s - w) / a[m + 1][m] + a[m][m + 1]; q = a[m + 1][m + 1] - z - r - s; r = a[m + 2][m + 1];
            s = Math.abs(p) + Math.abs(q) + Math.abs(r); p /= s; q /= s; r /= s;
            if (m === l) break;
            const u = Math.abs(a[m][m - 1]) * (Math.abs(q) + Math.abs(r));
            const v = Math.abs(p) * (Math.abs(a[m - 1][m - 1]) + Math.abs(z) + Math.abs(a[m + 1][m + 1]));
            if (u <= EPS * v) break;
          }
          for (let i = m + 2; i <= nn; i++) { a[i][i - 2] = 0; if (i !== m + 2) a[i][i - 3] = 0; }
          for (let k = m; k <= nn - 1; k++) {
            if (k !== m) {
              p = a[k][k - 1]; q = a[k + 1][k - 1]; r = k !== nn - 1 ? a[k + 2][k - 1] : 0;
              if ((x = Math.abs(p) + Math.abs(q) + Math.abs(r)) !== 0) { p /= x; q /= x; r /= x; }
            }
            const sq = Math.sqrt(p * p + q * q + r * r);
            if ((s = p >= 0 ? sq : -sq) !== 0) {
              if (k === m) { if (l !== m) a[k][k - 1] = -a[k][k - 1]; } else a[k][k - 1] = -s * x;
              p += s; x = p / s; y = q / s; z = r / s; q /= p; r /= p;
              for (let j = k; j <= nn; j++) {
                p = a[k][j] + q * a[k + 1][j];
                if (k !== nn - 1) { p += r * a[k + 2][j]; a[k + 2][j] -= p * z; }
                a[k + 1][j] -= p * y; a[k][j] -= p * x;
              }
              const mmin = nn < k + 3 ? nn : k + 3;
              for (let i = l; i <= mmin; i++) {
                p = x * a[i][k] + y * a[i][k + 1];
                if (k !== nn - 1) { p += z * a[i][k + 2]; a[i][k + 2] -= p * r; }
                a[i][k + 1] -= p * q; a[i][k] -= p;
              }
            }
          }
        }
      }
    } while (l < nn - 1);
  }
  const out = [];
  for (let i = 1; i <= n; i++) out.push({ re: wr[i] * scaleA, im: wi[i] * scaleA });
  return out.sort((u, v) => u.re - v.re || u.im - v.im);
}

/**
 * Eigenvalue problem a(x)u'' + b(x)u' + c(x)u = λu on [a, b] with u = 0 at both ends, by Chebyshev
 * collocation (interior (n−1)×(n−1) block, Trefethen Program 15) and eigenNonsymmetric.
 * Only roughly the lowest 2n/π modes are resolved (spurious large ones are expected): verify each
 * reported λ by repeating with a larger n. vectors: for the `vectors` lowest-|λ| real eigenvalues,
 * the eigenfunction on all n + 1 nodes (inverse iteration; max-norm 1, first lobe positive).
 * Returns { x, values: [{ re, im }] sorted by |λ|, vectors: [{ value, u }] }.
 */
export function chebyshevEigen({ n = 48, interval: iv = [-1, 1], a = 1, b = 0, c = 0, vectors = 0 } = {}) {
  order(n, 'chebyshevEigen', 3);
  const { x, D, D2 } = chebyshevDiffMatrix(n, iv);
  const m = n - 1, A = [];
  for (let i = 1; i < n; i++) {
    const ai = coef(a, x[i]), bi = coef(b, x[i]), ci = coef(c, x[i]);
    if (![ai, bi, ci].every(Number.isFinite)) fail('chebyshevEigen', 'coefficients must be finite on the interval');
    A.push(Float64Array.from({ length: m }, (_, j) => ai * D2[i][j + 1] + bi * D[i][j + 1] + (i === j + 1 ? ci : 0)));
  }
  const values = eigenNonsymmetric(A).sort((u, v) => Math.hypot(u.re, u.im) - Math.hypot(v.re, v.im) || u.im - v.im);
  const vecs = [];
  let scale = 0; for (const row of A) for (const v of row) scale = Math.max(scale, Math.abs(v));
  for (const { re, im } of values) {
    if (vecs.length >= vectors) break;
    if (im !== 0) continue;
    const shift = re + 1e-10 * scale;
    const B = A.map((row, i) => Float64Array.from(row, (v, j) => v - (i === j ? shift : 0)));
    let y = Float64Array.from({ length: m }, (_, i) => 1 + 0.1 * Math.sin(1.7 * i));
    for (let it = 0; it < 6; it++) {
      y = solveLinear(B, y);
      let peak = 0; for (const v of y) peak = Math.max(peak, Math.abs(v));
      for (let i = 0; i < m; i++) y[i] /= peak;
    }
    let peak = 0; for (const v of y) peak = Math.max(peak, Math.abs(v));
    const first = y.find(v => Math.abs(v) > 1e-3 * peak), sign = first > 0 ? 1 : -1;
    const u = new Float64Array(n + 1);
    for (let i = 0; i < m; i++) u[i + 1] = sign * y[i];
    vecs.push({ value: re, u });
  }
  return { x, values, vectors: vecs };
}

// ------------------------------------------------------------------ Fourier pseudo-spectral

/**
 * order-th derivative of periodic samples u (length n, a power of two, period `length`) by FFT:
 * multiply by (ik)^order; for odd orders the Nyquist mode is set to zero (keeps real data real).
 * Spectrally accurate for smooth periodic u; Gibbs oscillations for discontinuous u.
 */
export function fourierDerivative(u, length, ord = 1, out = new Float64Array(u?.length ?? 0)) {
  finiteArray(u, 'fourierDerivative', 'u');
  const n = u.length;
  if (!isPow2(n)) fail('fourierDerivative', 'length of u must be a power of two ≥ 2');
  if (!(Number.isFinite(length) && length > 0)) fail('fourierDerivative', 'period length must be positive');
  if (!(Number.isInteger(ord) && ord >= 0)) fail('fourierDerivative', 'order must be a non-negative integer');
  const re = Float64Array.from(u), im = new Float64Array(n), k = fftFreq(n, length / n);
  fft(re, im);
  for (let i = 0; i < n; i++) {
    if (ord % 2 === 1 && i === n / 2) { re[i] = im[i] = 0; continue; }
    // (ik)^ord = k^ord · i^ord
    const mag = k[i] ** ord, ph = ord % 4, r = re[i], m = im[i];
    if (ph === 0) { re[i] = mag * r; im[i] = mag * m; } else if (ph === 1) { re[i] = -mag * m; im[i] = mag * r; } else if (ph === 2) { re[i] = -mag * r; im[i] = -mag * m; } else { re[i] = mag * m; im[i] = -mag * r; }
  }
  fft(re, im, true);
  out.set(re);
  return out;
}

// ------------------------------------------------------------------ ETDRK4

/** ETDRK4 coefficients for one complex z = h·L by the Kassam–Trefethen contour mean (M points, radius 1). */
function etdCoefficients(zr, zi, h, M) {
  const acc = new Float64Array(8);   // Q, f1, f2, f3 (re, im)
  for (let j = 0; j < M; j++) {
    const th = 2 * Math.PI * (j + 0.5) / M, rr = zr + Math.cos(th), ri = zi + Math.sin(th);
    const ex = Math.exp(rr), er = ex * Math.cos(ri), ei = ex * Math.sin(ri);
    const eh = Math.exp(rr / 2), e2r = eh * Math.cos(ri / 2), e2i = eh * Math.sin(ri / 2);
    const d = rr * rr + ri * ri, invr = rr / d, invi = -ri / d;               // 1/r
    const r2r = rr * rr - ri * ri, r2i = 2 * rr * ri;                          // r²
    const inv3r = invr * (invr * invr - invi * invi) - invi * (2 * invr * invi), inv3i = invr * (2 * invr * invi) + invi * (invr * invr - invi * invi);
    const mul = (ar, ai, br, bi) => [ar * br - ai * bi, ar * bi + ai * br];
    let [qr, qi] = mul(e2r - 1, e2i, invr, invi);
    acc[0] += qr; acc[1] += qi;
    // f1 = (−4 − r + e^r(4 − 3r + r²))/r³
    let [tr, ti] = mul(er, ei, 4 - 3 * rr + r2r, -3 * ri + r2i);
    [qr, qi] = mul(-4 - rr + tr, -ri + ti, inv3r, inv3i); acc[2] += qr; acc[3] += qi;
    // f2 = (2 + r + e^r(r − 2))/r³
    [tr, ti] = mul(er, ei, rr - 2, ri);
    [qr, qi] = mul(2 + rr + tr, ri + ti, inv3r, inv3i); acc[4] += qr; acc[5] += qi;
    // f3 = (−4 − 3r − r² + e^r(4 − r))/r³
    [tr, ti] = mul(er, ei, 4 - rr, -ri);
    [qr, qi] = mul(-4 - 3 * rr - r2r + tr, -3 * ri - r2i + ti, inv3r, inv3i); acc[6] += qr; acc[7] += qi;
  }
  for (let i = 0; i < 8; i++) acc[i] *= h / M;
  return acc;
}

/**
 * Semilinear periodic evolution ∂u/∂t = L u + N(u) by ETDRK4 (Cox & Matthews 2002, with the
 * Kassam & Trefethen 2005 contour-integral coefficients: no cancellation for small or stiff hL).
 * L is diagonal in Fourier space: linear(k) → [Re, Im] of its symbol at angular wavenumber k.
 * N̂ = multiplier(k) · FFT[nonlinear(u)]: nonlinear(re, im, outRe, outIm, t) writes the physical-space
 * term (outputs pre-zeroed); multiplier(k) → [Re, Im] (default 1), e.g. −ik/2 for −(u²/2)_x.
 * dealias: project both the initial field and N̂ onto |mode| < n/3 (the 2/3 rule, exact
 * dealiasing for quadratic nonlinearities; higher-degree products need stricter truncation).
 * Fourth order in dt; stiff L allowed (no linear CFL). The explicit
 * nonlinear part has an RK4-like limit: pass speed(re, im) (the largest advection speed) and
 * steps throw when speed·dt/dx > maxCourant (default 0.85: spectral advection |λ|dt = π·speed·dt/dx ≤ 2.7;
 * a strongly damping L at high k, as in Kuramoto–Sivashinsky, tolerates more).
 * n a power of two. State: s.re/s.im (physical, refreshed each step). Throws on non-finite values.
 * Returns { t, re, im, steps, step(), advanceTo(t), spectrum() }.
 */
export function etdrk4Periodic({ re0, im0, dx, dt, linear, nonlinear, multiplier, speed, maxCourant = 0.85, dealias = false, t0 = 0, contourPoints = 32 }) {
  const who = 'etdrk4Periodic';
  finiteArray(re0, who, 're0');
  const n = re0.length;
  if (!isPow2(n)) fail(who, 'grid size must be a power of two ≥ 2');
  if (im0) { finiteArray(im0, who, 'im0'); if (im0.length !== n) fail(who, 'im0 must match re0'); }
  if (!(dx > 0 && Number.isFinite(dx)) || !(dt > 0 && Number.isFinite(dt))) fail(who, 'dx and dt must be positive and finite');
  if (!Number.isFinite(t0) || !(Number.isFinite(maxCourant) && maxCourant > 0)) fail(who, 't0 must be finite and maxCourant positive');
  if (!(Number.isInteger(contourPoints) && contourPoints >= 16)) fail(who, 'contourPoints must be an integer ≥ 16');
  if (typeof linear !== 'function' || typeof nonlinear !== 'function') fail(who, 'linear(k) and nonlinear(re, im, outRe, outIm, t) are required');
  const k = fftFreq(n, dx);
  const E = new Float64Array(2 * n), E2 = new Float64Array(2 * n), Q = new Float64Array(2 * n), F1 = new Float64Array(2 * n), F2 = new Float64Array(2 * n), F3 = new Float64Array(2 * n), G = new Float64Array(2 * n);
  for (let i = 0; i < n; i++) {
    const [lr, li] = linear(k[i]);
    if (!(Number.isFinite(lr) && Number.isFinite(li))) fail(who, 'linear symbol must be finite');
    const zr = dt * lr, zi = dt * li;
    E[2 * i] = Math.exp(zr) * Math.cos(zi); E[2 * i + 1] = Math.exp(zr) * Math.sin(zi);
    E2[2 * i] = Math.exp(zr / 2) * Math.cos(zi / 2); E2[2 * i + 1] = Math.exp(zr / 2) * Math.sin(zi / 2);
    const c = etdCoefficients(zr, zi, dt, contourPoints);
    Q[2 * i] = c[0]; Q[2 * i + 1] = c[1]; F1[2 * i] = c[2]; F1[2 * i + 1] = c[3];
    F2[2 * i] = c[4]; F2[2 * i + 1] = c[5]; F3[2 * i] = c[6]; F3[2 * i + 1] = c[7];
    const [gr, gi] = multiplier ? multiplier(k[i]) : [1, 0];
    const mode = i < n / 2 ? i : i - n, keep = !dealias || 3 * Math.abs(mode) < n;
    G[2 * i] = keep ? gr : 0; G[2 * i + 1] = keep ? gi : 0;
  }
  const vr = Float64Array.from(re0), vi = im0 ? Float64Array.from(im0) : new Float64Array(n);
  fft(vr, vi);
  if (dealias) {
    for (let i = 0; i < n; i++) {
      const mode = i < n / 2 ? i : i - n;
      if (3 * Math.abs(mode) >= n) { vr[i] = 0; vi[i] = 0; }
    }
  }
  const re = new Float64Array(n), im = new Float64Array(n), outR = new Float64Array(n), outI = new Float64Array(n);
  const buf = () => [new Float64Array(n), new Float64Array(n)];
  const [nvR, nvI] = buf(), [naR, naI] = buf(), [nbR, nbI] = buf(), [ncR, ncI] = buf(), [aR, aI] = buf(), [bR, bI] = buf(), [cR, cI] = buf();
  const physical = (sr, si) => { re.set(sr); im.set(si); fft(re, im, true); };
  const N = (sr, si, dr, di, t) => {
    physical(sr, si); outR.fill(0); outI.fill(0);
    nonlinear(re, im, outR, outI, t);
    fft(outR, outI);
    for (let i = 0; i < n; i++) { const g = G[2 * i], h = G[2 * i + 1]; dr[i] = g * outR[i] - h * outI[i]; di[i] = g * outI[i] + h * outR[i]; }
  };
  // d = A·x + B·y (complex, per mode)
  const lin = (Ac, xr, xi, Bc, yr, yi, dr, di) => {
    for (let i = 0; i < n; i++) {
      const p = Ac[2 * i], q = Ac[2 * i + 1], r = Bc[2 * i], w = Bc[2 * i + 1];
      const or = p * xr[i] - q * xi[i] + r * yr[i] - w * yi[i], oi = p * xi[i] + q * xr[i] + r * yi[i] + w * yr[i];
      dr[i] = or; di[i] = oi;
    }
  };
  const s = {
    t: t0, re, im, steps: 0,
    courant: 0,
    step() {
      if (speed) {
        s.courant = speed(re, im) * dt / dx;
        if (!(s.courant <= maxCourant)) fail(who, `advective Courant number speed·dt/dx = ${s.courant} exceeds ${maxCourant}; reduce dt`);
      }
      const t = s.t;
      N(vr, vi, nvR, nvI, t);
      lin(E2, vr, vi, Q, nvR, nvI, aR, aI);
      N(aR, aI, naR, naI, t + dt / 2);
      lin(E2, vr, vi, Q, naR, naI, bR, bI);
      N(bR, bI, nbR, nbI, t + dt / 2);
      for (let i = 0; i < n; i++) { nbR[i] = 2 * nbR[i] - nvR[i]; nbI[i] = 2 * nbI[i] - nvI[i]; }   // 2Nb − Nv (Nb kept below via c)
      lin(E2, aR, aI, Q, nbR, nbI, cR, cI);
      for (let i = 0; i < n; i++) { nbR[i] = (nbR[i] + nvR[i]) / 2; nbI[i] = (nbI[i] + nvI[i]) / 2; }   // restore Nb
      N(cR, cI, ncR, ncI, t + dt);
      let check = 0;
      for (let i = 0; i < n; i++) {
        const j = 2 * i;
        const sr = naR[i] + nbR[i], si = naI[i] + nbI[i];
        const xr = E[j] * vr[i] - E[j + 1] * vi[i] + F1[j] * nvR[i] - F1[j + 1] * nvI[i] + 2 * (F2[j] * sr - F2[j + 1] * si) + F3[j] * ncR[i] - F3[j + 1] * ncI[i];
        const xi = E[j] * vi[i] + E[j + 1] * vr[i] + F1[j] * nvI[i] + F1[j + 1] * nvR[i] + 2 * (F2[j] * si + F2[j + 1] * sr) + F3[j] * ncI[i] + F3[j + 1] * ncR[i];
        vr[i] = xr; vi[i] = xi; check += xr + xi;
      }
      if (!Number.isFinite(check)) fail(who, `solution became non-finite at t = ${t}; reduce dt or improve resolution`);
      physical(vr, vi);
      s.t = t + dt; s.steps++;
    },
    advanceTo(target, maxSteps = 1e7) {
      if (!(Number.isFinite(target))) fail(who, 'target time must be finite');
      let count = 0;
      while (s.t < target - 1e-9 * dt) { if (count++ >= maxSteps) return false; s.step(); }
      return true;
    },
    /** Fourier coefficients (FFT ordering, unnormalised as numerics.fft): copies { re, im, k }. */
    spectrum() { return { re: Float64Array.from(vr), im: Float64Array.from(vi), k: Float64Array.from(k) }; },
  };
  physical(vr, vi);
  return s;
}

const realField = (u0, who) => { finiteArray(u0, who, 'u0'); if (!isPow2(u0.length)) fail(who, 'grid size must be a power of two ≥ 2'); };
const maxAbs = v => { let m = 0; for (let i = 0; i < v.length; i++) { const a = Math.abs(v[i]); if (a > m) m = a; } return m; };
const sum = v => { let s = 0; for (let i = 0; i < v.length; i++) s += v[i]; return s; };
const sumSq = v => { let s = 0; for (let i = 0; i < v.length; i++) s += v[i] * v[i]; return s; };

/**
 * Korteweg–de Vries u_t + 6uu_x + u_xxx = 0 (periodic, n = 2^m) by ETDRK4: L = ik³, N̂ = −3ik·FFT(u²).
 * Solitons: kdvSoliton, kdvTwoSoliton. Invariants (conserved by the PDE; the scheme keeps them to
 * its truncation error, mass() to round-off): mass ∫u, momentum ∫u², energy ∫(½u_x² − u³).
 * Advective Courant 6·max|u|·dt/dx ≤ 0.85 enforced. Returns the etdrk4Periodic stepper with u
 * (= re) and mass(), momentum(), energy().
 */
export function kdv({ u0, dx, dt, dealias = false, t0 = 0 }) {
  realField(u0, 'kdv');
  const nyquist = -Math.PI / dx;
  const s = etdrk4Periodic({
    re0: u0, dx, dt, t0, dealias,
    linear: k => [0, k === nyquist ? 0 : k * k * k],
    multiplier: k => [0, k === nyquist ? 0 : -3 * k],
    nonlinear: (re, _im, out) => { for (let i = 0; i < re.length; i++) out[i] = re[i] * re[i]; },
    speed: re => 6 * maxAbs(re),
  });
  const L = u0.length * dx, ux = new Float64Array(u0.length);
  s.u = s.re;
  s.mass = () => sum(s.re) * dx;
  s.momentum = () => sumSq(s.re) * dx;
  s.energy = () => {
    fourierDerivative(s.re, L, 2, ux);
    let e = 0; for (let i = 0; i < ux.length; i++) e += -0.5 * s.re[i] * ux[i] - s.re[i] ** 3;
    return e * dx;
  };
  return s;
}

/**
 * Kuramoto–Sivashinsky u_t + uu_x + u_xx + u_xxxx = 0 (periodic) by ETDRK4: L = k² − k⁴,
 * N̂ = −(ik/2)·FFT(u²) (Kassam & Trefethen 2005 §6: L = 32π, n = 128, dt = ¼). Chaotic for
 * L ≳ 2π·(several): compare statistics, not trajectories, beyond the Lyapunov time.
 * mean() is conserved to round-off; energy() = ½∫u² obeys dE/dt = ∫(u_x² − u_xx²) = dissipation().
 * Courant max|u|·dt/dx ≤ maxCourant (default 2: the −k⁴ damping stabilises the high modes; the
 * Kassam–Trefethen run reaches ≈ 1.1) enforced. Returns the stepper with u, mean(), energy(), dissipation().
 */
export function kuramotoSivashinsky({ u0, dx, dt, dealias = false, t0 = 0, maxCourant = 2 }) {
  realField(u0, 'kuramotoSivashinsky');
  const nyquist = -Math.PI / dx;
  const s = etdrk4Periodic({
    re0: u0, dx, dt, t0, dealias, maxCourant,
    linear: k => [k * k - k ** 4, 0],
    multiplier: k => [0, k === nyquist ? 0 : -k / 2],
    nonlinear: (re, _im, out) => { for (let i = 0; i < re.length; i++) out[i] = re[i] * re[i]; },
    speed: re => maxAbs(re),
  });
  const n = u0.length, L = n * dx, d2 = new Float64Array(n);
  s.u = s.re;
  s.mean = () => sum(s.re) / n;
  s.energy = () => 0.5 * sumSq(s.re) * dx;
  s.dissipation = () => {
    fourierDerivative(s.re, L, 2, d2);
    let value = 0;
    for (let i = 0; i < n; i++) value += -s.re[i] * d2[i] - d2[i] * d2[i];
    return value * dx;
  };
  return s;
}

/**
 * Viscous Burgers u_t + uu_x = νu_xx (periodic) by ETDRK4: L = −νk², N̂ = −(ik/2)FFT(u²).
 * Exact benchmark: burgersColeHopf. Resolve the front: dx ≲ ν/max|u| (else Gibbs ripples);
 * Courant max|u|·dt/dx ≤ 0.85 enforced. mean() conserved to round-off.
 * Returns the stepper with u, mean(), energy() (½∫u², decreasing).
 */
export function viscousBurgers({ u0, nu, dx, dt, dealias = false, t0 = 0 }) {
  realField(u0, 'viscousBurgers');
  if (!(nu > 0 && Number.isFinite(nu))) fail('viscousBurgers', 'viscosity nu must be positive');
  const nyquist = -Math.PI / dx;
  const s = etdrk4Periodic({
    re0: u0, dx, dt, t0, dealias,
    linear: k => [-nu * k * k, 0],
    multiplier: k => [0, k === nyquist ? 0 : -k / 2],
    nonlinear: (re, _im, out) => { for (let i = 0; i < re.length; i++) out[i] = re[i] * re[i]; },
    speed: re => maxAbs(re),
  });
  s.u = s.re;
  s.mean = () => sum(s.re) / u0.length;
  s.energy = () => 0.5 * sumSq(s.re) * dx;
  return s;
}

// ------------------------------------------------------------------ split-step Fourier

const YOSHIDA = (() => { const w1 = 1 / (2 - 2 ** (1 / 3)); return [w1, 1 - 2 * w1, w1]; })();

/**
 * Nonlinear Schrödinger iψ_t + ½ψ_xx + κ|ψ|²ψ − Vψ = 0 (κ = 1 focusing, −1 defocusing; periodic,
 * n = 2^m) by split-step Fourier: Strang (order 2) or Yoshida's triple-jump composition (order 4).
 * Each sub-step is unitary, so norm() = ∫|ψ|² is conserved to round-off; hamiltonian() =
 * ∫[½|ψ_x|² + V|ψ|² − (κ/2)|ψ|⁴] and momentum() = ∫Im(ψ*ψ_x) (V = 0) drift only at O(dt^order).
 * Norm preservation does not guarantee perturbation stability or accuracy. Resolve the nonlinear
 * phase and spatial spectrum, and refine dt and dx together (numerical sideband instabilities
 * are possible). Returns { t, re, im, step(), advanceTo(), norm(), hamiltonian(), momentum() }.
 */
export function nls({ re0, im0, dx, dt, kappa = 1, V, order: ord = 2, t0 = 0 }) {
  const who = 'nls';
  finiteArray(re0, who, 're0');
  const n = re0.length;
  if (!isPow2(n)) fail(who, 'grid size must be a power of two ≥ 2');
  if (im0) { finiteArray(im0, who, 'im0'); if (im0.length !== n) fail(who, 'im0 must match re0'); }
  if (V) { finiteArray(V, who, 'V'); if (V.length !== n) fail(who, 'V must match re0'); }
  if (!(dx > 0 && dt > 0 && Number.isFinite(dx) && Number.isFinite(dt))) fail(who, 'dx and dt must be positive and finite');
  if (!Number.isFinite(kappa)) fail(who, 'kappa must be finite');
  if (ord !== 2 && ord !== 4) fail(who, 'order must be 2 or 4');
  const k = fftFreq(n, dx);
  const weights = ord === 2 ? [1] : YOSHIDA;
  const phases = weights.map(w => { const h = w * dt; return [Float64Array.from(k, q => Math.cos(-q * q * h / 2)), Float64Array.from(k, q => Math.sin(-q * q * h / 2))]; });
  const re = Float64Array.from(re0), im = im0 ? Float64Array.from(im0) : new Float64Array(n);
  const wr = new Float64Array(n), wi = new Float64Array(n);
  const potentialKick = h => {
    for (let i = 0; i < n; i++) {
      const ph = (kappa * (re[i] * re[i] + im[i] * im[i]) - (V ? V[i] : 0)) * h, c = Math.cos(ph), sn = Math.sin(ph), a = re[i];
      re[i] = a * c - im[i] * sn; im[i] = a * sn + im[i] * c;
    }
  };
  const strang = (w, [C, S]) => {
    const h = w * dt;
    potentialKick(h / 2);
    fft(re, im);
    for (let i = 0; i < n; i++) { const a = re[i]; re[i] = a * C[i] - im[i] * S[i]; im[i] = a * S[i] + im[i] * C[i]; }
    fft(re, im, true);
    potentialKick(h / 2);
  };
  const s = {
    t: t0, re, im, steps: 0,
    step() {
      for (let j = 0; j < weights.length; j++) strang(weights[j], phases[j]);
      let check = 0; for (let i = 0; i < n; i++) check += re[i] + im[i];
      if (!Number.isFinite(check)) fail(who, `ψ became non-finite at t = ${s.t}`);
      s.t += dt; s.steps++;
    },
    advanceTo(target, maxSteps = 1e7) {
      let count = 0;
      while (s.t < target - 1e-9 * dt) { if (count++ >= maxSteps) return false; s.step(); }
      return true;
    },
    norm() { let a = 0; for (let i = 0; i < n; i++) a += re[i] * re[i] + im[i] * im[i]; return a * dx; },
    hamiltonian() {
      wr.set(re); wi.set(im); fft(wr, wi);
      let kinetic = 0, potential = 0;
      for (let i = 0; i < n; i++) {
        kinetic += 0.5 * k[i] * k[i] * (wr[i] * wr[i] + wi[i] * wi[i]);
        const p = re[i] * re[i] + im[i] * im[i];
        potential += (V ? V[i] : 0) * p - 0.5 * kappa * p * p;
      }
      return (kinetic / n + potential) * dx;
    },
    momentum() {
      wr.set(re); wi.set(im); fft(wr, wi);
      let p = 0; for (let i = 0; i < n; i++) p += k[i] * (wr[i] * wr[i] + wi[i] * wi[i]);
      return p * dx / n;
    },
  };
  return s;
}

/**
 * Sine-Gordon u_tt = u_xx − sin u on a periodic grid (n = 2^m, period n·dx) with topological
 * charge `winding` Q: u(x + L) = u(x) + 2πQ (Q = 1 for one kink). The periodic part w = u − 2πQj/n
 * is evolved by a symplectic splitting: exact Fourier rotation of the free wave u_tt = u_xx and
 * kicks v −= h·sin u (Strang, order 2; order 4 by Yoshida composition). The nonlinear kick
 * imposes a stability restriction: the uniform small-amplitude mode needs dt < 2 for Strang;
 * energy() = ∫[½v² + ½u_x² + 1 − cos u] stays bounded with O(dt^order) oscillation (exact kink
 * energy 8/√(1 − v²), breather 16√(1 − ω²)). u0 must already carry the 2πQ ramp (u0[n−1] ≈ u0[0] +
 * 2πQ − slope·dx). Returns { t, u, v, step(), advanceTo(), energy(), charge }.
 */
export function sineGordon({ u0, v0, dx, dt, winding = 0, order: ord = 2, t0 = 0 }) {
  const who = 'sineGordon';
  realField(u0, who);
  const n = u0.length;
  if (v0) { finiteArray(v0, who, 'v0'); if (v0.length !== n) fail(who, 'v0 must match u0'); }
  if (!(dx > 0 && dt > 0 && Number.isFinite(dx) && Number.isFinite(dt))) fail(who, 'dx and dt must be positive and finite');
  if (!Number.isInteger(winding)) fail(who, 'winding must be an integer');
  if (ord !== 2 && ord !== 4) fail(who, 'order must be 2 or 4');
  const L = n * dx, k = fftFreq(n, dx), ramp = Float64Array.from({ length: n }, (_, j) => 2 * Math.PI * winding * j / n);
  const weights = ord === 2 ? [1] : YOSHIDA;
  const rot = weights.map(w => { const h = w * dt; return [Float64Array.from(k, q => Math.cos(Math.abs(q) * h)), Float64Array.from(k, q => (q === 0 ? h : Math.sin(Math.abs(q) * h) / Math.abs(q))), Float64Array.from(k, q => -Math.abs(q) * Math.sin(Math.abs(q) * h))]; });
  const u = Float64Array.from(u0), v = v0 ? Float64Array.from(v0) : new Float64Array(n);
  const wr = new Float64Array(n), wi = new Float64Array(n), vr = new Float64Array(n), vi = new Float64Array(n);
  const kick = h => { for (let i = 0; i < n; i++) v[i] -= h * Math.sin(u[i]); };
  const drift = ([C, S, M]) => {
    for (let i = 0; i < n; i++) { wr[i] = u[i] - ramp[i]; wi[i] = 0; vr[i] = v[i]; vi[i] = 0; }
    fft(wr, wi); fft(vr, vi);
    for (let i = 0; i < n; i++) {
      const ar = wr[i], ai = wi[i];
      wr[i] = C[i] * ar + S[i] * vr[i]; wi[i] = C[i] * ai + S[i] * vi[i];
      vr[i] = M[i] * ar + C[i] * vr[i]; vi[i] = M[i] * ai + C[i] * vi[i];
    }
    fft(wr, wi, true); fft(vr, vi, true);
    for (let i = 0; i < n; i++) { u[i] = wr[i] + ramp[i]; v[i] = vr[i]; }
  };
  const s = {
    t: t0, u, v, steps: 0, charge: winding,
    step() {
      for (let j = 0; j < weights.length; j++) { const h = weights[j] * dt; kick(h / 2); drift(rot[j]); kick(h / 2); }
      let check = 0; for (let i = 0; i < n; i++) check += u[i] + v[i];
      if (!Number.isFinite(check)) fail(who, `field became non-finite at t = ${s.t}`);
      s.t += dt; s.steps++;
    },
    advanceTo(target, maxSteps = 1e7) {
      let count = 0;
      while (s.t < target - 1e-9 * dt) { if (count++ >= maxSteps) return false; s.step(); }
      return true;
    },
    energy() {
      for (let i = 0; i < n; i++) { wr[i] = u[i] - ramp[i]; wi[i] = 0; }
      fft(wr, wi);
      let gradient = 0;
      for (let i = 0; i < n; i++) gradient += k[i] * k[i] * (wr[i] * wr[i] + wi[i] * wi[i]);
      const slope = 2 * Math.PI * winding / L;
      let e = 0;
      for (let i = 0; i < n; i++) e += 0.5 * v[i] * v[i] + 1 - Math.cos(u[i]);
      return (e + 0.5 * gradient / n) * dx + 0.5 * slope * slope * L;
    },
  };
  return s;
}

// ------------------------------------------------------------------ exact solutions

/** NLS bright soliton (κ = 1): ψ = A sech(A(x − x0 − vt)) exp(i[vx + ½(A² − v²)t + φ0]). Returns { re, im }. N = 2A, P = 2Av, H = Av² − A³/3. */
export function nlsSoliton(x, t, { amplitude: A = 1, velocity: v = 0, x0 = 0, phase = 0 } = {}) {
  const r = A / Math.cosh(A * (x - x0 - v * t)), ph = v * (x - x0) + 0.5 * (A * A - v * v) * t + phase;
  return { re: r * Math.cos(ph), im: r * Math.sin(ph) };
}

/** Peregrine breather (κ = 1, unit background): ψ = [1 − 4(1 + 2it)/(1 + 4x² + 4t²)]e^{it}; peak |ψ| = 3 at x = t = 0. */
export function nlsPeregrine(x, t) {
  const d = 1 + 4 * x * x + 4 * t * t, ar = 1 - 4 / d, ai = -8 * t / d, c = Math.cos(t), s = Math.sin(t);
  return { re: ar * c - ai * s, im: ar * s + ai * c };
}

/** KdV (u_t + 6uu_x + u_xxx = 0) one-soliton u = (k²/2) sech²(k(x − x0 − k²t)/2): speed k², height k²/2. */
export function kdvSoliton(x, t, { k = 1, x0 = 0 } = {}) {
  return 0.5 * k * k / Math.cosh(0.5 * k * (x - x0 - k * k * t)) ** 2;
}

/**
 * KdV two-soliton (Hirota): u = 2 ∂²ₓ ln F, F = 1 + e^{η1} + e^{η2} + A e^{η1+η2}, ηᵢ = kᵢ(x − xᵢ) − kᵢ³t,
 * A = ((k1 − k2)/(k1 + k2))². xᵢ is where soliton i would be at t = 0 without the other; the phase
 * shifts are kdvPhaseShifts(k1, k2). Evaluated in a scaled form that cannot overflow.
 */
export function kdvTwoSoliton(x, t, { k1, k2, x1 = 0, x2 = 0 }) {
  if (!(k1 > 0 && k2 > 0) || k1 === k2) fail('kdvTwoSoliton', 'k1, k2 must be positive and distinct');
  const A = ((k1 - k2) / (k1 + k2)) ** 2, e1 = k1 * (x - x1) - k1 ** 3 * t, e2 = k2 * (x - x2) - k2 ** 3 * t, e12 = e1 + e2 + Math.log(A);
  // F = Σ e^{θ_m} over θ = 0, η1, η2, η1+η2+ln A with slopes 0, k1, k2, k1+k2: ln F'' = Var of slopes under weights e^θ
  const th = [0, e1, e2, e12], sl = [0, k1, k2, k1 + k2], top = Math.max(...th);
  let w = 0, m1 = 0, m2 = 0;
  for (let i = 0; i < 4; i++) { const e = Math.exp(th[i] - top); w += e; m1 += e * sl[i]; m2 += e * sl[i] * sl[i]; }
  return 2 * (m2 / w - (m1 / w) ** 2);
}

/** Asymptotic position shifts of KdV solitons after a collision (k1 > k2): the faster moves forward by Δ1 = (2/k1)ln((k1+k2)/(k1−k2)), the slower back by Δ2 = (2/k2)ln(…). */
export function kdvPhaseShifts(k1, k2) {
  if (!(k1 > k2 && k2 > 0)) fail('kdvPhaseShifts', 'need k1 > k2 > 0');
  const l = Math.log((k1 + k2) / (k1 - k2));
  return { fast: 2 * l / k1, slow: -2 * l / k2 };
}

/**
 * Exact viscous Burgers solution u_t + uu_x = νu_xx on the line by the Cole–Hopf transform:
 * u(x, t) = ∫((x − y)/t)K dy / ∫K dy, K = exp(−[F(y) + (x − y)²/(2t)]/(2ν)), F(y) = ∫₀^y u₀.
 * primitive = F (for periodic u₀ with zero mean F is periodic). umax bounds |u₀| (sets the
 * domain of dependence |x − y| ≤ umax·t + 40√(νt)). Gauss–Legendre panels of width √(νt)/2,
 * exponent shifted by its minimum (no overflow for small ν). Accurate to ~1e−12 relative.
 */
export function burgersColeHopf({ primitive, nu, umax }, x, t) {
  if (typeof primitive !== 'function' || !(nu > 0) || !(umax >= 0) || !(t > 0)) fail('burgersColeHopf', 'need primitive(y), nu > 0, umax ≥ 0 and t > 0');
  const width = umax * t + 40 * Math.sqrt(nu * t), panel = Math.sqrt(nu * t) / 2, panels = Math.ceil(2 * width / panel);
  const h = 2 * width / panels, { x: nodes, w: weights } = gaussLegendreNodes(16);
  const G = y => primitive(y) + (x - y) ** 2 / (2 * t);
  let gmin = Infinity;
  for (let p = 0; p < panels; p++) for (const z of nodes) gmin = Math.min(gmin, G(x - width + h * (p + (z + 1) / 2)));
  let num = 0, den = 0;
  for (let p = 0; p < panels; p++) {
    for (let q = 0; q < nodes.length; q++) {
      const y = x - width + h * (p + (nodes[q] + 1) / 2), K = weights[q] * Math.exp(-(G(y) - gmin) / (2 * nu));
      num += (x - y) / t * K; den += K;
    }
  }
  return num / den;
}

/** Sine-Gordon kink (charge +1; −1 for an antikink): u = 4 arctan(exp(±γ(x − x0 − vt))), γ = 1/√(1 − v²). Returns { u, ut }. Energy 8γ. */
export function sineGordonKink(x, t, { velocity: v = 0, x0 = 0, charge = 1 } = {}) {
  if (!(Math.abs(v) < 1)) fail('sineGordonKink', '|velocity| must be < 1');
  const g = 1 / Math.sqrt(1 - v * v), xi = charge * g * (x - x0 - v * t);
  return { u: 4 * Math.atan(Math.exp(xi)), ut: -2 * charge * g * v / Math.cosh(xi) };
}

/** Sine-Gordon breather at rest: u = 4 arctan(√(1 − ω²) sin(ωt) / (ω cosh(√(1 − ω²)(x − x0)))), 0 < ω < 1. Returns { u, ut }. Energy 16√(1 − ω²). */
export function sineGordonBreather(x, t, { omega: w = 0.5, x0 = 0 } = {}) {
  if (!(w > 0 && w < 1)) fail('sineGordonBreather', 'omega must be in (0, 1)');
  const q = Math.sqrt(1 - w * w), ch = Math.cosh(q * (x - x0)), s = q * Math.sin(w * t) / (w * ch);
  return { u: 4 * Math.atan(s), ut: 4 * q * Math.cos(w * t) / ch / (1 + s * s) };
}
