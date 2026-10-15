// Stable probabilities, information measures, sample covariances and finite-spectrum
// canonical equilibrium. No physical model is inferred from a list of numbers.

import { erf } from './numerics.mjs';

function sum(values) {
  let result = 0, correction = 0;
  for (const value of values) {
    const next = result + value;
    correction += Math.abs(result) >= Math.abs(value) ? (result - next) + value : (value - next) + result;
    result = next;
  }
  return result + correction;
}
function finite(value, who) {
  if (!Number.isFinite(value)) throw new Error(`${who}: numerical range exceeded; rescale the problem`);
  return value;
}
function logMaximum(values) {
  let maximum = -Infinity;
  for (const value of values) {
    if (!(Number.isFinite(value) || value === -Infinity)) throw new Error('log weights: values must be finite or minus infinity');
    maximum = Math.max(maximum, value);
  }
  return maximum;
}
function distribution(p, who) {
  if (!p?.length || !Array.from(p).every(value => Number.isFinite(value) && value >= 0 && value <= 1) || Math.abs(sum(p) - 1) > 1e-12) throw new Error(`${who}: probabilities must be non-negative and normalized`);
}

/** ln(sum(exp(logWeights))) by maximum shifting. Empty/all-impossible support returns -Infinity. */
export function logSumExp(logWeights) {
  const maximum = logMaximum(logWeights);
  if (maximum === -Infinity) return -Infinity;
  let skippedMaximum = false;
  const otherWeights = [];
  for (const value of logWeights) {
    if (!skippedMaximum && value === maximum) { skippedMaximum = true; continue; }
    otherWeights.push(Math.exp(value - maximum));
  }
  return finite(maximum + Math.log1p(sum(otherWeights)), 'logSumExp');
}

/** Normalize non-negative finite weights without overflowing their sum. */
export function normalizeWeights(weights) {
  let maximum = 0;
  if (!weights?.length) throw new Error('normalizeWeights: non-empty weights required');
  for (const value of weights) {
    if (!(Number.isFinite(value) && value >= 0)) throw new Error('normalizeWeights: weights must be finite and non-negative');
    maximum = Math.max(maximum, value);
  }
  if (maximum === 0) throw new Error('normalizeWeights: total weight must be positive');
  const scaled = Float64Array.from(weights, value => value / maximum), total = sum(scaled);
  return Float64Array.from(scaled, value => value / total);
}

/** Stable softmax: subtract maximum and divide by the shifted sum, not exp(x-logSumExp(x)). */
export function normalizeLogWeights(logWeights) {
  const maximum = logMaximum(logWeights);
  if (maximum === -Infinity) throw new Error('normalizeLogWeights: at least one possible state is required');
  return normalizeWeights(Float64Array.from(logWeights, value => Math.exp(value - maximum)));
}

/** Shannon entropy of a normalized discrete distribution, in nats (base e) or bits (base 2). */
export function shannonEntropy(p, { base = Math.E } = {}) {
  distribution(p, 'shannonEntropy');
  if (!(Number.isFinite(base) && base > 1)) throw new Error('shannonEntropy: logarithm base must exceed one');
  return sum(Array.from(p, value => value === 0 ? 0 : -value * Math.log(value))) / Math.log(base);
}

/** KL divergence D(p||q), in nats. Positive p at zero q correctly returns +Infinity. */
export function relativeEntropy(p, q) {
  distribution(p, 'relativeEntropy'); distribution(q, 'relativeEntropy');
  if (p.length !== q.length) throw new Error('relativeEntropy: support lengths must match');
  const terms = [];
  for (let i = 0; i < p.length; i++) {
    if (p[i] === 0) { terms.push(q[i]); continue; }
    if (q[i] === 0) return Infinity;
    const relativeDifference = (p[i] - q[i]) / q[i];
    // For normalized distributions, subtracting the linear term sums to zero. Evaluating
    // the remaining convex term avoids losing a tiny KL signal to normalization round-off.
    if (Math.abs(relativeDifference) < 1e-4) {
      const r = relativeDifference;
      terms.push(q[i] * r * r * (0.5 + r * (-1 / 6 + r * (1 / 12 + r * (-1 / 20 + r / 30)))));
    } else terms.push(p[i] * (Math.log(p[i]) - Math.log(q[i])) - (p[i] - q[i]));
  }
  return sum(terms);
}

/**
 * Mean and covariance of joint row samples. ddof=1 gives sample covariance, ddof=0 population.
 * Centers about the first sample before summation, retaining variance at large mean offsets.
 * Returns { mean, covariance, count, ddof }; this is not an independent-sample error estimate.
 */
export function covarianceMatrix(samples, { ddof = 1 } = {}) {
  const count = samples?.length, dimension = samples?.[0]?.length;
  if (!(count > 0 && dimension > 0 && Number.isSafeInteger(ddof) && ddof >= 0 && ddof < count)) throw new Error('covarianceMatrix: non-empty samples and count > ddof >= 0 required');
  for (const row of samples) if (row.length !== dimension || !Array.from(row).every(Number.isFinite)) throw new Error('covarianceMatrix: equal-length finite row samples required');
  const origin = samples[0], offsets = new Float64Array(dimension), means = new Float64Array(dimension), scales = new Float64Array(dimension);
  for (let j = 0; j < dimension; j++) {
    offsets[j] = finite(sum(samples.map(row => finite(row[j] - origin[j], 'covarianceMatrix') / count)), 'covarianceMatrix');
    means[j] = finite(origin[j] + offsets[j], 'covarianceMatrix');
    for (const row of samples) scales[j] = Math.max(scales[j], Math.abs((row[j] - origin[j]) - offsets[j]));
  }
  const covariance = Array.from({ length: dimension }, () => new Float64Array(dimension));
  for (let j = 0; j < dimension; j++) { for (let k = j; k < dimension; k++) {
    const terms = samples.map(row => scales[j] === 0 || scales[k] === 0 ? 0 : (((row[j] - origin[j]) - offsets[j]) / scales[j]) * (((row[k] - origin[k]) - offsets[k]) / scales[k]) / (count - ddof));
    const scaled = sum(terms), root = Math.sqrt(Math.abs(scaled));
    covariance[j][k] = covariance[k][j] = finite(Math.sign(scaled) * (scales[j] * root) * (scales[k] * root), 'covarianceMatrix');
  } }
  return { mean: means, covariance, count, ddof };
}

/**
 * Canonical distribution over temperature-independent discrete energy levels with positive
 * integer degeneracies. kT is an ENERGY (caller supplies k_B*T). Entropy is S/k_B in nats,
 * heatCapacity is C/k_B. Includes degeneracies within levels, not just entropy of level labels.
 * Returns a shifted partition function with its reference energy so extreme offsets never
 * require exp(-E/kT). It does not compute translational/other missing degrees of freedom.
 */
export function canonicalEnsemble(energies, { kT, degeneracies } = {}) {
  if (!(energies?.length > 0 && Array.from(energies).every(Number.isFinite) && Number.isFinite(kT) && kT > 0)) throw new Error('canonicalEnsemble: finite energy levels and positive finite kT required');
  const n = energies.length, g = degeneracies ? Array.from(degeneracies) : new Array(n).fill(1);
  if (g.length !== n || !g.every(value => Number.isSafeInteger(value) && value > 0)) throw new Error('canonicalEnsemble: positive integer degeneracy for each level required');
  let referenceEnergy = Infinity;
  for (const value of energies) referenceEnergy = Math.min(referenceEnergy, value);
  const gaps = new Float64Array(n);
  const logWeights = Float64Array.from(energies, (value, i) => {
    const difference = value - referenceEnergy;
    const gap = Number.isFinite(difference) ? difference / kT : value / kT - referenceEnergy / kT;
    gaps[i] = gap; return Math.log(g[i]) - gap;
  });
  const probabilities = normalizeLogWeights(logWeights), logShiftedPartition = logSumExp(logWeights);
  let referenceIndex = 0;
  for (let i = 1; i < n; i++) if (probabilities[i] > probabilities[referenceIndex]) referenceIndex = i;
  const origin = energies[referenceIndex], differences = new Float64Array(n);
  let scale = 0;
  for (let i = 0; i < n; i++) { if (probabilities[i] > 0) {
    differences[i] = finite(energies[i] - origin, 'canonicalEnsemble');
    scale = Math.max(scale, Math.abs(differences[i]));
  } }
  const normalized = Float64Array.from(differences, value => scale === 0 ? 0 : value / scale);
  const meanOffset = sum(Array.from(normalized, (value, i) => probabilities[i] * value));
  const meanEnergy = finite(origin + scale * meanOffset, 'canonicalEnsemble');
  const scaledVariance = sum(Array.from(normalized, (value, i) => probabilities[i] * (value - meanOffset) ** 2));
  const deviation = finite(scale * Math.sqrt(scaledVariance), 'canonicalEnsemble');
  const energyVariance = finite(deviation * deviation, 'canonicalEnsemble');
  const heatCapacity = finite((deviation / kT) ** 2, 'canonicalEnsemble');
  const entropy = logShiftedPartition + sum(Array.from(probabilities, (value, i) => value === 0 ? 0 : value * gaps[i]));
  const freeEnergy = finite(referenceEnergy - kT * logShiftedPartition, 'canonicalEnsemble');
  return { probabilities, meanEnergy, energyVariance, heatCapacity, entropy, freeEnergy, referenceEnergy, logShiftedPartition };
}

/**
 * Speeds of gas molecules of mass m (kg) at temperature T (K), with σ = √(kT/m):
 *   maxwellSpeed(T, m)   speeds in the gas, f(v) ∝ v² e^{−v²/2σ²}
 *   effusiveSpeed(T, m)  speeds of molecules leaving a small hole (beam from an oven), ∝ v³ e^{−v²/2σ²}
 * Each returns { pdf(v), cdf(v), sample(rand), mean, mostProbable, rms, sigma }. Sampling is exact:
 * σ√χ²₃ (three normals) and σ√χ²₄ = σ√(−2 ln(u₁u₂)).
 */
const KB = 1.380649e-23;
export function maxwellSpeed(T, m) {
  const sigma = Math.sqrt(KB * T / m);
  if (!(sigma > 0)) throw new Error('maxwellSpeed: T and m must be positive');
  const c = Math.sqrt(2 / Math.PI) / sigma ** 3;
  return {
    sigma, mean: sigma * Math.sqrt(8 / Math.PI), mostProbable: sigma * Math.SQRT2, rms: sigma * Math.sqrt(3),
    pdf: v => (v < 0 ? 0 : c * v * v * Math.exp(-v * v / (2 * sigma * sigma))),
    cdf: v => { if (v <= 0) return 0; const x = v / sigma; return erf(x / Math.SQRT2) - Math.sqrt(2 / Math.PI) * x * Math.exp(-x * x / 2); },
    sample: rand => { let s2 = 0; for (let k = 0; k < 3; k++) { const z = normalFrom(rand); s2 += z * z; } return sigma * Math.sqrt(s2); },
  };
}
export function effusiveSpeed(T, m) {
  const sigma = Math.sqrt(KB * T / m);
  if (!(sigma > 0)) throw new Error('effusiveSpeed: T and m must be positive');
  return {
    sigma, mean: sigma * 3 * Math.sqrt(Math.PI / 8), mostProbable: sigma * Math.sqrt(3), rms: sigma * 2,
    pdf: v => (v < 0 ? 0 : v ** 3 / (2 * sigma ** 4) * Math.exp(-v * v / (2 * sigma * sigma))),
    cdf: v => { if (v <= 0) return 0; const x2 = v * v / (sigma * sigma); return 1 - (1 + x2 / 2) * Math.exp(-x2 / 2); },
    sample: rand => { let u1 = 0, u2 = 0; while (u1 === 0) u1 = rand(); while (u2 === 0) u2 = rand(); return sigma * Math.sqrt(-2 * Math.log(u1 * u2)); },
  };
}
function normalFrom(rand) {
  if (rand.normal) return rand.normal();
  let u = 0; while (u === 0) u = rand();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rand());
}
