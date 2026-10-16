// Quasi-Monte Carlo: low-discrepancy point sets for integration over the unit cube, their
// randomisation for honest error bars, and discrepancy measures for small sets.
// Plain ES module (browser + Node). Deterministic for fixed seeds.
//
//   sobol                 Sobol sequence, Joe–Kuo (2008) direction numbers, up to 21 dimensions;
//                         scramble: 'none' | 'shift' (random digital shift) | 'owen' (nested
//                         uniform scrambling, hash-based: Burley 2020 / Laine–Karras)
//   halton                Halton sequence in the first d prime bases; scramble: 'none' |
//                         'permute' (random digit permutations) | 'shift' (Cranley–Patterson)
//   qmcIntegrate          randomised-QMC integration with replicate error bars (also plain MC)
//   starDiscrepancy       exact star discrepancy D* by grid enumeration (small sets only)
//   l2StarDiscrepancy     exact L2 star discrepancy (Warnock's formula), any size
//   SOBOL_DIRECTIONS      the [s, a, m₁…m_s] table used for dimensions 2…21
//
// Why: for smooth integrands QMC error falls like ~N⁻¹ (scrambled Sobol: RMSE ~N⁻³ᐟ² for very
// smooth ones) against N⁻¹ᐟ² for Monte Carlo. A single deterministic QMC estimate has no usable
// error estimate; average R independently scrambled replicates and quote their standard error.
//
// References: S. Joe & F. Y. Kuo, "Constructing Sobol sequences with better two-dimensional
// projections", SIAM J. Sci. Comput. 30 (2008) 2635 (file new-joe-kuo-6.21201); A. B. Owen,
// "Scrambling Sobol' and Niederreiter–Xing points", J. Complexity 14 (1998) 466; B. Burley,
// "Practical hash-based Owen scrambling", JCGT 9 (2020); T. T. Warnock (1972) for the L2 formula;
// A. Genz (1984) test integrands.

import { rng, rngFor } from './numerics.mjs';

/** Joe–Kuo new-joe-kuo-6.21201 rows for dimensions 2…21: [degree s, a (polynomial bits), m₁…m_s]. */
export const SOBOL_DIRECTIONS = Object.freeze([
  [1, 0, 1], [2, 1, 1, 3], [3, 1, 1, 3, 1], [3, 2, 1, 1, 1], [4, 1, 1, 1, 3, 3], [4, 4, 1, 3, 5, 13],
  [5, 2, 1, 1, 5, 5, 17], [5, 4, 1, 1, 5, 5, 5], [5, 7, 1, 1, 7, 11, 19], [5, 11, 1, 1, 5, 1, 1],
  [5, 13, 1, 1, 1, 3, 11], [5, 14, 1, 3, 5, 5, 31], [6, 1, 1, 3, 3, 9, 7, 49], [6, 13, 1, 1, 1, 15, 21, 21],
  [6, 16, 1, 3, 1, 13, 27, 49], [6, 19, 1, 1, 1, 15, 7, 5], [6, 22, 1, 3, 1, 15, 13, 25], [6, 25, 1, 1, 5, 5, 19, 61],
  [7, 1, 1, 3, 7, 11, 23, 15, 103], [7, 4, 1, 3, 7, 13, 13, 15, 69],
].map(row => Object.freeze(row)));
export const SOBOL_MAX_DIM = SOBOL_DIRECTIONS.length + 1;

const BITS = 32, TWO32 = 2 ** 32;

/** Direction integers V[k] (k = 0…31, scaled by 2³²) for Sobol dimension j (0-based). */
function directionIntegers(j) {
  const V = new Uint32Array(BITS);
  if (j === 0) { for (let k = 0; k < BITS; k++) V[k] = 2 ** (BITS - 1 - k); return V; }
  const [s, a, ...m] = SOBOL_DIRECTIONS[j - 1];
  for (let k = 0; k < Math.min(s, BITS); k++) V[k] = m[k] * 2 ** (BITS - 1 - k);
  for (let k = s; k < BITS; k++) {
    let v = V[k - s] ^ (V[k - s] >>> s);
    for (let i = 1; i < s; i++) if ((a >>> (s - 1 - i)) & 1) v ^= V[k - i];
    V[k] = v;
  }
  return V;
}

const reverseBits = x => {
  let v = x >>> 0;
  v = ((v >>> 1) & 0x55555555) | ((v & 0x55555555) << 1);
  v = ((v >>> 2) & 0x33333333) | ((v & 0x33333333) << 2);
  v = ((v >>> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
  v = ((v >>> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
  return ((v >>> 16) | (v << 16)) >>> 0;
};
// Laine–Karras hash: each output bit depends only on itself and lower bits, so applied to
// bit-reversed values it permutes every dyadic subinterval within its parent (nested scrambling).
const laineKarras = (x, seed) => {
  let v = (x + seed) >>> 0;
  v ^= Math.imul(v, 0x6c50b47c); v ^= Math.imul(v, 0xb82f1e52);
  v ^= Math.imul(v, 0xc7afe638); v ^= Math.imul(v, 0x8d22f6e6);
  return v >>> 0;
};
const owenScramble = (x, seed) => reverseBits(laineKarras(reverseBits(x), seed));
const uint32From = rand => Math.floor(rand() * TWO32) >>> 0;

const dimension = (dim, max, who) => {
  if (!(Number.isSafeInteger(dim) && dim >= 1 && dim <= max)) throw new Error(`${who}: dimension must be an integer in 1…${max}`);
};
const indexCheck = (i, who) => {
  if (!(Number.isSafeInteger(i) && i >= 0)) throw new Error(`${who}: index must be a non-negative integer`);
};

/**
 * Sobol low-discrepancy sequence in [0, 1)^dim, dim ≤ 21 (Joe–Kuo 2008 direction numbers,
 * Gray-code order, 32-bit resolution: at most 2³² points). Point 0 is the origin when
 * unscrambled. Properties: every block of 2^m consecutive points starting at a multiple of 2^m
 * is a (t, m, d)-net; in particular each coordinate takes each value j/2^m once. Use n = 2^m
 * points (balance is lost for other n).
 * scramble: 'none' (default), 'shift' (XOR with a random 32-bit word per coordinate) or 'owen'
 * (hash-based nested uniform scrambling; preserves the net property and makes each point
 * uniform on the cube, so independent seeds give unbiased replicates). Scrambled coordinates are
 * offset to the centre of their 2⁻³² cell (never exactly 0).
 * Returns { dim, point(i) → Float64Array, next() → Float64Array, points(n) → Float64Array[] }
 * (next/points continue from `skip`).
 */
export function sobol(dim, { scramble = 'none', seed = 1, skip = 0 } = {}) {
  dimension(dim, SOBOL_MAX_DIM, 'sobol'); indexCheck(skip, 'sobol');
  if (!['none', 'shift', 'owen'].includes(scramble)) throw new Error('sobol: scramble must be none, shift or owen');
  const V = Array.from({ length: dim }, (_, j) => directionIntegers(j));
  const rand = rng(seed), keys = Uint32Array.from({ length: dim }, () => uint32From(rand));
  const offset = scramble === 'none' ? 0 : 0.5;
  const finish = (x, j) => {
    const y = scramble === 'shift' ? (x ^ keys[j]) >>> 0 : scramble === 'owen' ? owenScramble(x, keys[j]) : x;
    return (y + offset) / TWO32;
  };
  const raw = i => {
    indexCheck(i, 'sobol');
    if (i >= TWO32) throw new Error('sobol: index beyond 2³² points');
    const out = new Uint32Array(dim);
    let g = (i ^ Math.floor(i / 2)) >>> 0;
    for (let k = 0; g; k++, g >>>= 1) if (g & 1) for (let j = 0; j < dim; j++) out[j] ^= V[j][k];
    return out;
  };
  const state = raw(skip);
  let index = skip;
  const point = i => { const x = raw(i); return Float64Array.from(x, (v, j) => finish(v, j)); };
  const next = () => {
    if (index >= TWO32) throw new Error('sobol: sequence exhausted (2³² points)');
    const out = Float64Array.from(state, (v, j) => finish(v, j));
    // Gray code: point index+1 differs from point index in the direction of index's lowest zero bit
    let c = 0, t = index;
    while (t & 1) { t = Math.floor(t / 2); c++; }
    for (let j = 0; j < dim; j++) state[j] ^= V[j][c];
    index++;
    return out;
  };
  const points = n => { indexCheck(n, 'sobol'); return Array.from({ length: n }, next); };
  return { dim, point, next, points };
}

const PRIMES = [2];
const primeAt = k => {
  for (let c = PRIMES[PRIMES.length - 1] + 1; PRIMES.length <= k; c++) {
    if (PRIMES.every(p => p * p > c || c % p !== 0)) PRIMES.push(c);
  }
  return PRIMES[k];
};

/**
 * Halton sequence in [0, 1)^dim: coordinate j is the radical inverse of the index in the j-th
 * prime base. Any dim, but beyond ~10 dimensions unscrambled Halton has strongly correlated
 * pairs of coordinates: scramble then (or use sobol). scramble: 'none' (default; point 0 is
 * the origin, consider skip ≥ 1), 'permute' (an independent random permutation of the digits
 * 0…b−1 at every digit position, digits taken to double precision: Owen 2017's random digit
 * permutation) or 'shift' (Cranley–Patterson random shift mod 1, which does not preserve
 * stratification). Returns the same interface as sobol.
 */
export function halton(dim, { scramble = 'none', seed = 1, skip = 0 } = {}) {
  dimension(dim, 1000, 'halton'); indexCheck(skip, 'halton');
  if (!['none', 'permute', 'shift'].includes(scramble)) throw new Error('halton: scramble must be none, permute or shift');
  const rand = rng(seed), bases = Array.from({ length: dim }, (_, j) => primeAt(j));
  const digits = bases.map(b => Math.ceil(53 * Math.LN2 / Math.log(b)));
  const perms = scramble === 'permute' ? bases.map((b, j) => Array.from({ length: digits[j] }, () => {
    const p = Array.from({ length: b }, (_, i) => i);
    for (let i = b - 1; i > 0; i--) { const r = Math.floor(rand() * (i + 1)); [p[i], p[r]] = [p[r], p[i]]; }
    return p;
  })) : null;
  const shifts = scramble === 'shift' ? bases.map(() => rand()) : null;
  const coord = (i, j) => {
    const b = bases[j];
    let x = 0, scale = 1 / b, n = i;
    if (perms) {
      for (let k = 0; k < digits[j]; k++, scale /= b) { const d = n % b; n = (n - d) / b; x += perms[j][k][d] * scale; }
    } else {
      while (n > 0) { const d = n % b; n = (n - d) / b; x += d * scale; scale /= b; }
    }
    if (shifts) { x += shifts[j]; if (x >= 1) x -= 1; }
    return Math.min(x, 1 - Number.EPSILON / 2);
  };
  const point = i => { indexCheck(i, 'halton'); return Float64Array.from({ length: dim }, (_, j) => coord(i, j)); };
  let index = skip;
  const next = () => point(index++);
  const points = n => { indexCheck(n, 'halton'); return Array.from({ length: n }, next); };
  return { dim, point, next, points };
}

/**
 * Randomised quasi-Monte Carlo integral of f over the box [lower, upper] (default unit cube):
 * R = `replicates` independent randomisations of an n-point set; value = mean of the R
 * estimates, error = their sample standard deviation / √R (an estimate of the standard error;
 * its square is unbiased for the estimator variance when replicates are independent).
 * Deterministic scramble: 'none' is rejected: identical replicates cannot estimate uncertainty.
 * Digital randomisations have finite resolution (Sobol: 32 bits); their error bar does not
 * measure the bias between the corresponding grid integral and the continuous integral.
 * sequence: 'sobol' (default; n must be a power of two), 'halton' or 'random' (plain Monte
 * Carlo with the same error bar, for comparison). f(x) receives a fresh Float64Array.
 * Conditions: f square-integrable; the QMC rate needs f smooth (bounded variation) — kinks or
 * discontinuities not aligned with the axes degrade it towards N⁻¹ᐟ²–N⁻³ᐟ⁴. Error bars from
 * R ≳ 10 replicates; cost n·R evaluations. Returns { value, error, n, replicates, estimates }.
 */
export function qmcIntegrate(f, dim, { n = 4096, replicates = 16, sequence = 'sobol', scramble, seed = 1, lower, upper } = {}) {
  if (typeof f !== 'function') throw new Error('qmcIntegrate: f must be a function');
  dimension(dim, sequence === 'sobol' ? SOBOL_MAX_DIM : 1000, 'qmcIntegrate');
  if (!(Number.isSafeInteger(n) && n > 0)) throw new Error('qmcIntegrate: n must be a positive integer');
  if (!(Number.isSafeInteger(replicates) && replicates >= 2)) throw new Error('qmcIntegrate: need at least 2 replicates for an error bar');
  if (sequence === 'sobol' && (n > TWO32 || !Number.isInteger(Math.log2(n)))) throw new Error('qmcIntegrate: Sobol needs n = 2^m ≤ 2^32 points');
  if (sequence !== 'random' && scramble === 'none') throw new Error('qmcIntegrate: independent randomisations are required for an error bar; use sobol/halton directly for deterministic points');
  const lo = lower ? Float64Array.from(lower) : new Float64Array(dim), hi = upper ? Float64Array.from(upper) : new Float64Array(dim).fill(1);
  if (lo.length !== dim || hi.length !== dim || lo.some((v, i) => !(Number.isFinite(v) && Number.isFinite(hi[i]) && hi[i] > v))) throw new Error('qmcIntegrate: need finite lower < upper per dimension');
  let volume = 1;
  for (let i = 0; i < dim; i++) volume *= hi[i] - lo[i];
  const estimates = new Float64Array(replicates);
  for (let r = 0; r < replicates; r++) {
    let source;
    if (sequence === 'sobol') source = sobol(dim, { scramble: scramble ?? 'owen', seed: seed * 7919 + r });
    else if (sequence === 'halton') source = halton(dim, { scramble: scramble ?? 'permute', seed: seed * 7919 + r });
    else if (sequence === 'random') { const rand = rngFor(seed, r); source = { next: () => Float64Array.from({ length: dim }, () => rand()) }; }
    else throw new Error('qmcIntegrate: sequence must be sobol, halton or random');
    let s = 0;
    for (let i = 0; i < n; i++) {
      const u = source.next();
      for (let k = 0; k < dim; k++) u[k] = lo[k] + (hi[k] - lo[k]) * u[k];
      const v = f(u);
      if (!Number.isFinite(v)) throw new Error('qmcIntegrate: integrand is not finite at a sample point');
      s += v;
    }
    estimates[r] = volume * s / n;
  }
  const mean = estimates.reduce((a, v) => a + v, 0) / replicates;
  const variance = estimates.reduce((a, v) => a + (v - mean) ** 2, 0) / (replicates - 1);
  return { value: mean, error: Math.sqrt(variance / replicates), n, replicates, estimates };
}

const pointSet = (points, who) => {
  const N = points?.length, d = points?.[0]?.length;
  if (!(N > 0 && d > 0)) throw new Error(`${who}: need a non-empty array of points`);
  for (const p of points) {
    if (p.length !== d) throw new Error(`${who}: all points need the same dimension`);
    for (const v of p) if (!(v >= 0 && v <= 1)) throw new Error(`${who}: coordinates must lie in [0, 1]`);
  }
  return [N, d];
};

/**
 * Exact star discrepancy D* = sup_y |#{x ∈ [0, y)}/N − vol[0, y)| by enumerating the grid of
 * critical boxes (each upper corner coordinate one of the points' coordinates or 1; Niederreiter
 * 1972). Cost ~ (N+1)^d·N·d: fine for N ≲ 1000 in 2D, ≲ 100 in 3D; throws above maxWork
 * (default 1e8) — use l2StarDiscrepancy for larger sets. Returns D*.
 */
export function starDiscrepancy(points, { maxWork = 1e8 } = {}) {
  const [N, d] = pointSet(points, 'starDiscrepancy');
  const grids = Array.from({ length: d }, (_, k) => {
    const vals = [...new Set(points.map(p => p[k]))];
    vals.push(1);
    return Float64Array.from(new Set(vals)).sort();
  });
  const cells = grids.reduce((a, g) => a * g.length, 1);
  if (cells * N * d > maxWork) throw new Error(`starDiscrepancy: ${cells} boxes × ${N} points exceeds maxWork; use l2StarDiscrepancy`);
  // Sort points along the first axis so the innermost loop can stop early.
  const pts = points.map(p => Float64Array.from(p)).sort((p, q) => p[0] - q[0]);
  const idx = new Int32Array(d);
  let D = 0;
  for (let c = 0; c < cells; c++) {
    let t = c;
    for (let k = d - 1; k >= 0; k--) { idx[k] = t % grids[k].length; t = (t - idx[k]) / grids[k].length; }
    let vol = 1;
    for (let k = 0; k < d; k++) vol *= grids[k][idx[k]];
    const y0 = grids[0][idx[0]];
    let open = 0, closed = 0;
    for (const p of pts) {
      if (p[0] > y0) break;
      let le = true, lt = true;
      for (let k = 0; k < d; k++) { const y = grids[k][idx[k]]; if (p[k] > y) { le = false; break; } if (p[k] === y) lt = false; }
      if (le) { closed++; if (lt) open++; }
    }
    D = Math.max(D, vol - open / N, closed / N - vol);
  }
  return D;
}

/**
 * Exact L2 star discrepancy T* = (∫ (#{x ∈ [0, y)}/N − vol[0, y))² dy)^{1/2} by Warnock's
 * formula, O(N²·d). A cheap, smooth quality measure for any set; random points have
 * E[T*²] = (2⁻ᵈ − 3⁻ᵈ)/N. Returns T*.
 */
export function l2StarDiscrepancy(points) {
  const [N, d] = pointSet(points, 'l2StarDiscrepancy');
  let a = 0, b = 0;
  for (let i = 0; i < N; i++) {
    let prod = 1;
    for (let k = 0; k < d; k++) prod *= 1 - points[i][k] ** 2;
    a += prod;
    for (let j = 0; j < N; j++) {
      let q = 1;
      for (let k = 0; k < d; k++) q *= 1 - Math.max(points[i][k], points[j][k]);
      b += q;
    }
  }
  const T2 = 3 ** -d - 2 ** (1 - d) * a / N + b / (N * N);
  return Math.sqrt(Math.max(0, T2));
}
