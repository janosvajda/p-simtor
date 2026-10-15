// Working with data, independent of subject: interpolation that does not invent wiggles,
// histograms with uncertainties, goodness-of-fit tests and error bars on derived quantities.
//
//   interpolate       linear / monotone cubic through tabulated points
//   histogram         one-shot histogram with Poisson errors
//   createHistogram   incremental (weighted) histogram with under/overflow and density errors
//   ksTest            Kolmogorov–Smirnov test against a model CDF
//   chiSquareTest     Pearson χ² of binned counts against expectations, with p-value
//   jackknife         blocked jackknife: value, error and bias of any function of means
//   bootstrap         blocked bootstrap: error and confidence quantiles of any function of means
//   meanSem / pull    mean ± standard error of independent samples; (measured − expected)/error
//   jackknifeRefit    blocked jackknife of any estimator that needs the samples themselves (a refit)
//   curveFit          nonlinear least squares (Levenberg–Marquardt) with parameter covariance

import { rng, gammaQ } from './numerics.mjs';
import { solveLinear } from './linalg.mjs';

/**
 * Interpolate tabulated data (xs strictly increasing). method 'linear' or 'monotone'
 * (Fritsch–Carlson monotone cubic: no overshoot between data points, so it never invents
 * maxima that the data do not have). Outside [x₀, x_last] returns NaN unless extrapolate: true
 * (then linear from the end segment). Returns a function x → y.
 */
export function interpolate(xs, ys, { method = 'monotone', extrapolate = false } = {}) {
  const n = xs.length;
  if (n < 2 || ys.length !== n) throw new Error('interpolate: need ≥ 2 points and equal-length xs, ys');
  for (let i = 1; i < n; i++) if (!(xs[i] > xs[i - 1])) throw new Error('interpolate: xs must be strictly increasing');
  const h = Array.from({ length: n - 1 }, (_, i) => xs[i + 1] - xs[i]);
  const delta = h.map((hi, i) => (ys[i + 1] - ys[i]) / hi);
  let m = null;
  if (method === 'monotone') {
    m = new Float64Array(n);
    m[0] = delta[0]; m[n - 1] = delta[n - 2];
    for (let i = 1; i < n - 1; i++) m[i] = delta[i - 1] * delta[i] <= 0 ? 0 : (delta[i - 1] + delta[i]) / 2;
    for (let i = 0; i < n - 1; i++) {
      if (delta[i] === 0) { m[i] = 0; m[i + 1] = 0; continue; }
      const a = m[i] / delta[i], b = m[i + 1] / delta[i], s = a * a + b * b;
      if (s > 9) { const t = 3 / Math.sqrt(s); m[i] = t * a * delta[i]; m[i + 1] = t * b * delta[i]; }
    }
  } else if (method !== 'linear') throw new Error(`interpolate: unknown method "${method}"`);
  return x => {
    if (x < xs[0] || x > xs[n - 1]) {
      if (!extrapolate) return NaN;
      const i = x < xs[0] ? 0 : n - 2;
      return ys[i] + delta[i] * (x - xs[i]);
    }
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; xs[mid] <= x ? (lo = mid) : (hi = mid); }
    const t = (x - xs[lo]) / h[lo];
    if (!m) return ys[lo] + t * (ys[hi] - ys[lo]);
    const t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * ys[lo] + (t3 - 2 * t2 + t) * h[lo] * m[lo] + (-2 * t3 + 3 * t2) * ys[hi] + (t3 - t2) * h[lo] * m[hi];
  };
}

/**
 * Histogram with Poisson uncertainties, normalised as a probability density if density: true.
 * Returns { edges, centers, counts, values, errors }.
 */
export function histogram(samples, { bins = 30, range, density = true } = {}) {
  const [lo, hi] = range ?? [Math.min(...samples), Math.max(...samples)];
  if (!(hi > lo)) throw new Error('histogram: empty or degenerate range');
  const w = (hi - lo) / bins, counts = new Float64Array(bins);
  let inRange = 0;
  for (const s of samples) {
    if (s < lo || s > hi) continue;
    counts[Math.min(bins - 1, Math.floor((s - lo) / w))]++; inRange++;
  }
  const norm = density ? 1 / (samples.length * w) : 1;
  return {
    edges: Array.from({ length: bins + 1 }, (_, i) => lo + i * w),
    centers: Array.from({ length: bins }, (_, i) => lo + (i + 0.5) * w),
    counts, inRange,
    values: Array.from(counts, c => c * norm),
    errors: Array.from(counts, c => Math.sqrt(c) * norm),
  };
}

/**
 * One-sample Kolmogorov–Smirnov test of samples against a model CDF. Returns { D, p }:
 * D is the largest CDF difference; p the asymptotic probability of a D at least this large if
 * the samples came from the model (assumes independent samples; thin correlated chains first).
 */
export function ksTest(samples, cdf) {
  const s = Float64Array.from(samples).sort(), n = s.length;
  if (n < 5) throw new Error('ksTest: need at least 5 samples');
  let D = 0;
  for (let i = 0; i < n; i++) {
    const F = cdf(s[i]);
    D = Math.max(D, Math.abs((i + 1) / n - F), Math.abs(F - i / n));
  }
  const l = (Math.sqrt(n) + 0.12 + 0.11 / Math.sqrt(n)) * D;      // Stephens' correction
  let p = 0;
  for (let k = 1; k <= 100; k++) p += 2 * (k % 2 ? 1 : -1) * Math.exp(-2 * k * k * l * l);
  return { D, p: Math.min(1, Math.max(0, p)) };
}

/**
 * Incremental histogram over range [lo, hi] with `bins` equal bins and optional weights.
 * add(x, w = 1); x outside the range goes to underflow / overflow, NaN to nan. total sums every
 * weight added (in range or not). errors() = √Σw² per bin (Poisson for unit weights);
 * density(total) → { values: counts/(total·width), errors } (default total: everything added).
 */
export function createHistogram({ bins = 30, range } = {}) {
  if (!(Number.isInteger(bins) && bins > 0)) throw new Error('createHistogram: bins must be a positive integer');
  const [lo, hi] = range ?? [];
  if (!(Number.isFinite(lo) && Number.isFinite(hi) && hi > lo)) throw new Error('createHistogram: range [lo, hi] with finite hi > lo is required');
  const width = (hi - lo) / bins;
  const h = {
    bins, range: [lo, hi], width,
    counts: new Float64Array(bins), sumW2: new Float64Array(bins),
    total: 0, underflow: 0, overflow: 0, nan: 0,
    edges: Array.from({ length: bins + 1 }, (_, i) => lo + i * width),
    centers: Array.from({ length: bins }, (_, i) => lo + (i + 0.5) * width),
    add(x, w = 1) {
      h.total += w;
      if (Number.isNaN(x)) h.nan += w;
      else if (x < lo) h.underflow += w;
      else if (x > hi) h.overflow += w;
      else { const i = Math.min(bins - 1, Math.floor((x - lo) / width)); h.counts[i] += w; h.sumW2[i] += w * w; }
      return h;
    },
    errors: () => Array.from(h.sumW2, Math.sqrt),
    density(total = h.total) {
      const k = 1 / (total * width);
      return { values: Array.from(h.counts, c => c * k), errors: Array.from(h.sumW2, s => Math.sqrt(s) * k) };
    },
  };
  return h;
}

/**
 * Pearson χ² test of observed counts against expected counts: Σ(o−e)²/e, or Σ((o−e)/σ)² when
 * errors σ are given. Bins with e ≤ 0 (or σ ≤ 0) are skipped and counted. dof = used bins − ddof
 * (ddof: constraints, e.g. 1 if the expectation was normalised to the total, plus fitted parameters).
 * Returns { chi2, dof, p (upper tail Q(dof/2, χ²/2)), z = (χ² − dof)/√(2 dof), skipped }.
 * For Pearson's form the χ² law needs e ≳ 5 per bin; merge sparse bins first.
 */
export function chiSquareTest(observed, expected, { ddof = 0, errors } = {}) {
  if (observed.length !== expected.length || (errors && errors.length !== observed.length)) throw new Error('chiSquareTest: observed, expected (and errors) must have equal length');
  let chi2 = 0, used = 0, skipped = 0;
  for (let i = 0; i < observed.length; i++) {
    const e = expected[i], s = errors ? errors[i] : Math.sqrt(e);
    if (!(e > 0) || !(s > 0)) { skipped++; continue; }
    chi2 += ((observed[i] - e) / s) ** 2; used++;
  }
  const dof = used - ddof;
  if (!(dof > 0)) throw new Error(`chiSquareTest: no degrees of freedom left (${used} usable bins, ddof = ${ddof})`);
  return { chi2, dof, p: Number.isFinite(chi2) ? gammaQ(dof / 2, chi2 / 2) : NaN, z: (chi2 - dof) / Math.sqrt(2 * dof), skipped };
}

/** Mean, standard error and sample standard deviation of INDEPENDENT samples (use meanWithError from numerics for correlated series). */
export function meanSem(xs) {
  const n = xs.length;
  if (!n) throw new Error('meanSem: no samples');
  let mean = 0;
  for (const x of xs) mean += x;
  mean /= n;
  let ss = 0;
  for (const x of xs) ss += (x - mean) ** 2;
  const sd = n > 1 ? Math.sqrt(ss / (n - 1)) : NaN;
  return { mean, sem: sd / Math.sqrt(n), sd, n };
}

/** Pull (measured − expected)/error: how many standard errors a measurement is from the expectation. */
export function pull(measured, error, expected) { return (measured - expected) / error; }

// Group samples (numbers or equal-length vectors) into B consecutive blocks of equal size,
// dropping the remainder at the START (it is the least thermalised part of a chain).
function blockSums(samples, bins, who) {
  const n = samples?.length ?? 0;
  if (n < 2) throw new Error(`${who}: need at least 2 samples`);
  const B = bins ?? Math.min(n, 50);
  if (!(Number.isInteger(B) && B >= 2 && B <= n)) throw new Error(`${who}: bins must be an integer in [2, ${n}] (got ${B})`);
  const scalar = typeof samples[0] === 'number', dim = scalar ? 1 : samples[0].length;
  if (!scalar && !(dim > 0)) throw new Error(`${who}: samples must be numbers or non-empty arrays`);
  const size = Math.floor(n / B), start = n - B * size;
  const sums = Array.from({ length: B }, () => new Float64Array(dim)), total = new Float64Array(dim);
  for (let b = 0; b < B; b++) {
    for (let i = start + b * size, end = i + size; i < end; i++) {
      const s = samples[i];
      if (scalar) sums[b][0] += s;
      else { if (s.length !== dim) throw new Error(`${who}: sample ${i} has length ${s.length}, expected ${dim}`); for (let k = 0; k < dim; k++) sums[b][k] += s[k]; }
    }
    for (let k = 0; k < dim; k++) total[k] += sums[b][k];
  }
  // mean of a set of blocks, in the samples' own shape (number or Float64Array)
  const meanOf = (vec, count) => scalar ? vec[0] / (count * size) : vec.map(v => v / (count * size));
  return { B, size, scalar, dim, sums, total, meanOf };
}

const asArray = v => (typeof v === 'number' ? [v] : Array.from(v));

/**
 * Blocked jackknife error of θ = estimator(mean). samples: numbers, or arrays of observables
 * measured together (e.g. all Wilson loops of one sweep). Consecutive samples form `bins` blocks
 * (default min(n, 50)); the estimator sees the full mean and each leave-one-block-out mean
 * (a number or a Float64Array, like the samples) and returns a number or an array. Blocks much
 * longer than the autocorrelation time τ_int also make the error honest for correlated chains.
 * Returns { value: θ(full mean), error: √((B−1)/B Σ(θᵢ−θ̄)²), bias: (B−1)(θ̄−θ), blocks: B }
 * (arrays when the estimator returns arrays). A NaN from any block gives a NaN error, never a
 * silently smaller one: show it as "unresolved".
 */
export function jackknife(samples, estimator, { bins } = {}) {
  const { B, sums, total, meanOf } = blockSums(samples, bins, 'jackknife');
  const raw = estimator(meanOf(total.slice(), B)), vector = typeof raw !== 'number', value = asArray(raw);
  const m = value.length, thetas = sums.map(s => {
    const t = asArray(estimator(meanOf(total.map((v, k) => v - s[k]), B - 1)));
    if (t.length !== m) throw new Error('jackknife: estimator returned a different length for a block');
    return t;
  });
  const error = new Array(m), bias = new Array(m);
  for (let k = 0; k < m; k++) {
    const mean = thetas.reduce((a, t) => a + t[k], 0) / B;
    error[k] = Math.sqrt((B - 1) / B * thetas.reduce((a, t) => a + (t[k] - mean) ** 2, 0));
    bias[k] = (B - 1) * (mean - value[k]);
  }
  return vector ? { value, error, bias, blocks: B } : { value: value[0], error: error[0], bias: bias[0], blocks: B };
}

/**
 * Blocked bootstrap of θ = estimator(mean), same sample conventions as jackknife: each of
 * `resamples` replicas draws B blocks with replacement. rand defaults to rng(1), so results are
 * reproducible. Returns { value: θ(full mean), error: standard deviation over replicas,
 * quantiles(p): replica quantile(s) for p in [0, 1] (e.g. 0.16 and 0.84 for a 68 % interval),
 * blocks, resamples }; arrays when the estimator returns arrays. NaN replicas give NaN.
 */
export function bootstrap(samples, estimator, { resamples = 1000, bins, rand = rng(1) } = {}) {
  if (typeof rand !== 'function') throw new Error('bootstrap: rand must be a function (e.g. rng(seed) from core/numerics.mjs)');
  if (!(Number.isInteger(resamples) && resamples >= 2)) throw new Error('bootstrap: resamples must be an integer ≥ 2');
  const { B, dim, sums, total, meanOf } = blockSums(samples, bins, 'bootstrap');
  const raw = estimator(meanOf(total.slice(), B)), vector = typeof raw !== 'number', value = asArray(raw);
  const m = value.length, reps = Array.from({ length: m }, () => new Float64Array(resamples));
  const acc = new Float64Array(dim);
  for (let r = 0; r < resamples; r++) {
    acc.fill(0);
    for (let b = 0; b < B; b++) { const s = sums[Math.floor(rand() * B)]; for (let k = 0; k < dim; k++) acc[k] += s[k]; }
    const t = asArray(estimator(meanOf(acc.slice(), B)));
    if (t.length !== m) throw new Error('bootstrap: estimator returned a different length for a replica');
    for (let k = 0; k < m; k++) reps[k][r] = t[k];
  }
  const error = reps.map(x => {
    const mu = x.reduce((a, v) => a + v, 0) / resamples;
    return Math.sqrt(x.reduce((a, v) => a + (v - mu) ** 2, 0) / (resamples - 1));
  });
  const sorted = reps.map(x => (x.some(Number.isNaN) ? null : x.slice().sort()));
  const quant = p => sorted.map(s => {
    if (!s) return NaN;
    const pos = Math.min(1, Math.max(0, p)) * (resamples - 1), i = Math.floor(pos), f = pos - i;
    return i + 1 < resamples ? s[i] * (1 - f) + s[i + 1] * f : s[i];
  });
  const quantiles = p => (Array.isArray(p) ? p.map(q => (vector ? quant(q) : quant(q)[0])) : (vector ? quant(p) : quant(p)[0]));
  return vector ? { value, error, quantiles, blocks: B, resamples } : { value: value[0], error: error[0], quantiles, blocks: B, resamples };
}

/**
 * Blocked jackknife for estimators that need the samples themselves, not just their means (a fit,
 * a median, a spectral estimate): fit(subset) is re-run on the data with each block left out.
 * samples: an array of anything (points, segments, events); consecutive samples form `bins`
 * blocks (default min(n, 20)); the remainder at the start is dropped. fit returns a number or an
 * array. Returns { value, error, bias, blocks } like jackknife. Costs B + 1 fits.
 */
export function jackknifeRefit(samples, fit, { bins } = {}) {
  const n = samples.length, B = Math.min(n, bins ?? 20);
  if (!(B >= 2)) throw new Error('jackknifeRefit: need at least 2 samples and 2 blocks');
  const size = Math.floor(n / B), used = samples.slice(n - B * size);
  const asArr = v => (typeof v === 'number' ? [v] : Array.from(v));
  const raw = fit(used), vector = typeof raw !== 'number', value = asArr(raw);
  const thetas = Array.from({ length: B }, (_, b) => asArr(fit([...used.slice(0, b * size), ...used.slice((b + 1) * size)])));
  const mean = value.map((_, k) => thetas.reduce((a, t) => a + t[k], 0) / B);
  const error = value.map((_, k) => Math.sqrt((B - 1) / B * thetas.reduce((a, t) => a + (t[k] - mean[k]) ** 2, 0)));
  const bias = value.map((v, k) => (B - 1) * (mean[k] - v));
  return vector ? { value, error, bias, blocks: B } : { value: value[0], error: error[0], bias: bias[0], blocks: B };
}

/**
 * Nonlinear least squares by Levenberg–Marquardt: minimise Σ ((y_i − model(x_i, p)) / σ_i)².
 * model(x, p) → number; jacobian(x, p) → ∂model/∂p (optional; central differences otherwise).
 * sigma: per-point errors (number or array). With sigma the covariance is (JᵀWJ)⁻¹; without it,
 * that matrix is scaled by χ²/dof (errors estimated from the scatter). Returns
 * { p, errors, cov, chi2, dof, p_value (only with sigma), iterations, converged }.
 * Start values matter: a damped oscillation needs a frequency guess within ~20 %.
 */
export function curveFit(model, xs, ys, p0, { sigma, jacobian, maxIterations = 200, tol = 1e-12 } = {}) {
  const n = xs.length, m = p0.length;
  if (ys.length !== n || !(n > m)) throw new Error('curveFit: need equal-length xs, ys and more points than parameters');
  const w = Array.from({ length: n }, (_, i) => { const s = sigma == null ? 1 : typeof sigma === 'number' ? sigma : sigma[i]; if (!(s > 0)) throw new Error('curveFit: sigma must be positive'); return 1 / (s * s); });
  const resid = p => xs.map((x, i) => ys[i] - model(x, p));
  const chi2Of = r => r.reduce((a, v, i) => a + w[i] * v * v, 0);
  const jac = p => xs.map(x => {
    if (jacobian) return Array.from(jacobian(x, p));
    return p.map((pk, k) => {
      const h = 1e-6 * Math.max(1, Math.abs(pk)), up = p.slice(), dn = p.slice();
      up[k] += h; dn[k] -= h;
      return (model(x, up) - model(x, dn)) / (2 * h);
    });
  });
  const normal = (J, r) => {
    const A = Array.from({ length: m }, () => new Array(m).fill(0)), g = new Array(m).fill(0);
    for (let i = 0; i < n; i++) for (let a = 0; a < m; a++) { g[a] += w[i] * J[i][a] * r[i]; for (let b = 0; b < m; b++) A[a][b] += w[i] * J[i][a] * J[i][b]; }
    return { A, g };
  };
  let p = Array.from(p0), r = resid(p), chi2 = chi2Of(r), lambda = 1e-3, it = 0, converged = false;
  if (!Number.isFinite(chi2)) throw new Error('curveFit: the model is not finite at the start values');
  for (; it < maxIterations; it++) {
    const J = jac(p), { A, g } = normal(J, r);
    let improved = false;
    for (let tries = 0; tries < 30 && !improved; tries++) {
      const M = A.map((row, a) => row.map((v, b) => (a === b ? v * (1 + lambda) + 1e-300 : v)));
      let dp;
      try { dp = solveLinear(M, g); } catch { lambda *= 10; continue; }
      const trial = p.map((v, k) => v + dp[k]), rt = resid(trial), ct = chi2Of(rt);
      if (Number.isFinite(ct) && ct <= chi2) {
        const rel = (chi2 - ct) / Math.max(chi2, 1e-300), step = Math.max(...dp.map((d, k) => Math.abs(d) / Math.max(Math.abs(trial[k]), 1e-300)));
        p = trial; r = rt; chi2 = ct; lambda = Math.max(lambda / 10, 1e-12); improved = true;
        if (rel < tol || step < tol) converged = true;
      } else lambda *= 10;
    }
    if (!improved) { converged = true; break; }   // no downhill step left: at the minimum (to round-off)
    if (converged) break;
  }
  const { A } = normal(jac(p), r), dof = n - m;
  const inv = Array.from({ length: m }, (_a, k) => solveLinear(A, Array.from({ length: m }, (_b, j) => (j === k ? 1 : 0))));
  const scale = sigma == null ? chi2 / dof : 1;
  const cov = inv.map(row => Array.from(row, v => v * scale));
  return { p, errors: cov.map((row, k) => Math.sqrt(row[k])), cov, chi2, dof, p_value: sigma == null ? undefined : gammaQ(dof / 2, chi2 / 2), iterations: it, converged };
}

