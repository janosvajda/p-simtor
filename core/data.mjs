// Working with data, independent of subject: interpolation that does not invent wiggles,
// histograms with uncertainties, and a goodness-of-fit test against a model distribution.

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
