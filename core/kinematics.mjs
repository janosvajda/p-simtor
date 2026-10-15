// Relativistic kinematics and Monte Carlo event sampling, independent of subject. Units are the
// caller's with c = 1; 4-vectors are [E, px, py, pz] with metric (+, −, −, −). All randomness
// comes from a seeded rand = rng(seed) (core/numerics.mjs), so every event is reproducible.
//
//   dot4 / mass / add4         Minkowski product, invariant mass, sum of 4-vectors
//   boost / restFrameBoost     active Lorentz boost by velocity β; β that brings P to rest
//   randomUnitVector           uniform direction on the unit sphere S²
//   randomRotation             uniform (Haar) random rotation in SO(3), via a unit quaternion
//   orthonormalFrame           right-handed frame [e1, e2, n] around a direction n
//   twoBodyDecay               isotropic M → m1 m2 in the rest frame (Källén momentum)
//   rambo                      n-body phase space, flat in Lorentz-invariant measure (Kleiss, Stirling, Ellis 1986)
//   acceptReject               unweighting by accept/reject, failing loudly if a weight exceeds wmax

import { brent, lgamma } from './numerics.mjs';

const needRand = (rand, who) => { if (typeof rand !== 'function') throw new Error(`${who}: pass rand = rng(seed) from core/numerics.mjs`); };
const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const unit = (v, who) => { const r = Math.hypot(v[0], v[1], v[2]); if (!(r > 0)) throw new Error(`${who}: zero-length vector`); return [v[0] / r, v[1] / r, v[2] / r]; };

/** Minkowski product p·q = E_p E_q − p⃗·q⃗. */
export function dot4(p, q) { return p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3]; }

/** Invariant mass √max(0, p·p) (round-off can make p·p slightly negative for massless p). */
export function mass(p) { return Math.sqrt(Math.max(0, dot4(p, p))); }

/** Sum of any number of 4-vectors. */
export function add4(...ps) {
  const s = [0, 0, 0, 0];
  for (const p of ps) for (let k = 0; k < 4; k++) s[k] += p[k];
  return s;
}

/**
 * Active Lorentz boost of p by velocity β = [βx, βy, βz] (|β| < 1): a particle at rest is given
 * velocity β, i.e. (m, 0⃗) → (γm, γmβ⃗). Equivalently, p as seen from a frame moving with −β.
 * boost(boost(p, β), −β) = p.
 */
export function boost(p, beta) {
  const b2 = dot3(beta, beta);
  if (!(b2 < 1)) throw new Error(`boost: need |β| < 1 (got |β|² = ${b2})`);
  if (b2 === 0) return [p[0], p[1], p[2], p[3]];
  const g = 1 / Math.sqrt(1 - b2), bp = beta[0] * p[1] + beta[1] * p[2] + beta[2] * p[3];
  const k = (g - 1) * bp / b2 + g * p[0];
  return [g * (p[0] + bp), p[1] + k * beta[0], p[2] + k * beta[1], p[3] + k * beta[2]];
}

/** Velocity β = −P⃗/E such that boost(P, β) is at rest (P must be timelike with E > 0). */
export function restFrameBoost(P) {
  if (!(P[0] > 0) || !(dot4(P, P) > 0)) throw new Error('restFrameBoost: P must be timelike with E > 0');
  return [-P[1] / P[0], -P[2] / P[0], -P[3] / P[0]];
}

/** Uniform random direction on S² (cos θ uniform in [−1, 1], φ uniform). */
export function randomUnitVector(rand) {
  needRand(rand, 'randomUnitVector');
  const c = 2 * rand() - 1, s = Math.sqrt(Math.max(0, 1 - c * c)), phi = 2 * Math.PI * rand();
  return [s * Math.cos(phi), s * Math.sin(phi), c];
}

/** Uniform random rotation (Haar measure on SO(3)): 3×3 row-major [[..],[..],[..]] from a uniform unit quaternion (Shoemake 1992). */
export function randomRotation(rand) {
  needRand(rand, 'randomRotation');
  const u1 = rand(), u2 = 2 * Math.PI * rand(), u3 = 2 * Math.PI * rand();
  const a = Math.sqrt(1 - u1), b = Math.sqrt(u1);
  const x = a * Math.sin(u2), y = a * Math.cos(u2), z = b * Math.sin(u3), w = b * Math.cos(u3);
  return [
    [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
    [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
    [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
  ];
}

/** Right-handed orthonormal frame [e1, e2, n̂] with n̂ = n/|n|, e1 ⟂ n̂ and e1 × e2 = n̂. */
export function orthonormalFrame(n) {
  const nz = unit(n, 'orthonormalFrame');
  const ax = Math.abs(nz[0]) <= Math.abs(nz[1]) && Math.abs(nz[0]) <= Math.abs(nz[2]) ? [1, 0, 0] : Math.abs(nz[1]) <= Math.abs(nz[2]) ? [0, 1, 0] : [0, 0, 1];
  const e1 = unit(cross(ax, nz), 'orthonormalFrame');
  return [e1, cross(nz, e1), nz];
}

/** Two-body momentum in the rest frame of M: √λ(M², m1², m2²)/(2M), λ the Källén function. */
function twoBodyMomentum(M, m1, m2) {
  return Math.sqrt(Math.max(0, (M * M - (m1 + m2) ** 2) * (M * M - (m1 - m2) ** 2))) / (2 * M);
}

/** Isotropic two-body decay M → m1 + m2 in the rest frame of M. Returns [p1, p2] (back to back). */
export function twoBodyDecay(M, m1, m2, rand) {
  needRand(rand, 'twoBodyDecay');
  if (!(m1 >= 0 && m2 >= 0 && M >= m1 + m2)) throw new Error(`twoBodyDecay: kinematically forbidden (M = ${M} < m1 + m2 = ${m1 + m2})`);
  const k = twoBodyMomentum(M, m1, m2), n = randomUnitVector(rand);
  return [[Math.sqrt(m1 * m1 + k * k), k * n[0], k * n[1], k * n[2]], [Math.sqrt(m2 * m2 + k * k), -k * n[0], -k * n[1], -k * n[2]]];
}

/**
 * RAMBO (Kleiss, Stirling, Ellis, Comput. Phys. Commun. 40 (1986) 359): n ≥ 2 momenta with total
 * (√s, 0, 0, 0), distributed uniformly in Lorentz-invariant phase space for massless particles.
 * weight is the phase-space volume element with the usual conventions,
 * ∫ Πᵢ d³pᵢ/((2π)³ 2Eᵢ) (2π)⁴ δ⁴(P − Σpᵢ): constant (2π)^{4−3n} (π/2)^{n−1} s^{n−2}/((n−1)!(n−2)!)
 * for massless; with masses the momenta are rescaled by ξ (Σ√(mᵢ² + ξ²|pᵢ|²) = √s) and the
 * weight varies event by event (the average weight is the massive phase-space volume).
 * Returns { momenta: [[E, px, py, pz], ...], weight }.
 */
export function rambo(n, sqrtS, rand, masses) {
  needRand(rand, 'rambo');
  if (!(Number.isInteger(n) && n >= 2)) throw new Error('rambo: n must be an integer ≥ 2');
  if (!(sqrtS > 0)) throw new Error('rambo: √s must be positive');
  if (masses != null && (masses.length !== n || !Array.from(masses).every(m => m >= 0))) throw new Error('rambo: masses must be n non-negative numbers');
  const msum = masses ? masses.reduce((a, m) => a + m, 0) : 0;
  if (!(msum < sqrtS)) throw new Error(`rambo: Σm = ${msum} ≥ √s = ${sqrtS}`);
  // massless isotropic momenta q with energy density q0 e^{−q0}
  const q = [], Q = [0, 0, 0, 0];
  for (let i = 0; i < n; i++) {
    const c = 2 * rand() - 1, s = Math.sqrt(1 - c * c), f = 2 * Math.PI * rand();
    let u = 0; while (u === 0) u = rand() * rand();
    const e = -Math.log(u), qi = [e, e * s * Math.cos(f), e * s * Math.sin(f), e * c];
    q.push(qi); for (let k = 0; k < 4; k++) Q[k] += qi[k];
  }
  // conformal transformation to total momentum (√s, 0⃗)
  const M = Math.sqrt(dot4(Q, Q)), b = [-Q[1] / M, -Q[2] / M, -Q[3] / M], x = sqrtS / M, g = Q[0] / M, a = 1 / (1 + g);
  const p = q.map(qi => {
    const bq = b[0] * qi[1] + b[1] * qi[2] + b[2] * qi[3], h = a * bq + qi[0];
    return [x * (g * qi[0] + bq), x * (qi[1] + h * b[0]), x * (qi[2] + h * b[1]), x * (qi[3] + h * b[2])];
  });
  const s = sqrtS * sqrtS;
  let logW = (4 - 3 * n) * Math.log(2 * Math.PI) + (n - 1) * Math.log(Math.PI / 2) + (n - 2) * Math.log(s) - lgamma(n) - lgamma(n - 1);
  if (!masses || msum === 0) return { momenta: p, weight: Math.exp(logW) };
  // massive: rescale 3-momenta by ξ, recompute energies, multiply by the Jacobian
  const xi = brent(z => p.reduce((acc, pi, i) => acc + Math.sqrt(masses[i] ** 2 + z * z * pi[0] * pi[0]), 0) - sqrtS, 0, 1);
  let w2 = 0, w3 = 0;
  const k = p.map((pi, i) => {
    const E = Math.sqrt(masses[i] ** 2 + xi * xi * pi[0] * pi[0]), kk = xi * pi[0];
    w2 += Math.log(kk / E); w3 += kk * kk / E;
    return [E, xi * pi[1], xi * pi[2], xi * pi[3]];
  });
  logW += (2 * n - 3) * Math.log(xi) + w2 + Math.log(sqrtS / w3);
  return { momenta: k, weight: Math.exp(logW) };
}

/**
 * Unweighted sampling by accept/reject: draws x = sample() and accepts it with probability
 * weight(x)/wmax. Throws if any weight exceeds wmax (or is negative / NaN), since a wrong bound
 * would otherwise bias the result silently. Returns { value, tries }.
 */
export function acceptReject(sample, weight, wmax, rand, { maxTries = 1e6 } = {}) {
  needRand(rand, 'acceptReject');
  if (!(wmax > 0 && Number.isFinite(wmax))) throw new Error('acceptReject: wmax must be finite and positive');
  for (let tries = 1; tries <= maxTries; tries++) {
    const x = sample(), w = weight(x);
    if (!(w >= 0)) throw new Error(`acceptReject: weight ${w} is negative or NaN`);
    if (w > wmax) throw new Error(`acceptReject: weight ${w} exceeds wmax = ${wmax} (raise the bound; results so far would be biased)`);
    if (rand() * wmax < w) return { value: x, tries };
  }
  throw new Error(`acceptReject: no acceptance in ${maxTries} tries`);
}
