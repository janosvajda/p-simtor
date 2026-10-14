import { test } from 'node:test';
import assert from 'node:assert/strict';
import { logSumExp, normalizeWeights, normalizeLogWeights, shannonEntropy, relativeEntropy, covarianceMatrix, canonicalEnsemble } from '../core/probability.mjs';

const close = (a, b, tol = 1e-12) => assert.ok(Math.abs(a - b) <= tol, `${a} vs ${b} (tolerance ${tol})`);
const relativeClose = (a, b, tol = 1e-12) => close(a / b, 1, tol);

test('stable weights work where direct exponentials/sums overflow or underflow', () => {
  const p = normalizeLogWeights([-1000, -1001, -Infinity]);
  close(p[0], 1 / (1 + Math.exp(-1))); close(p[1], 1 - p[0]); close(p[2], 0, 0);
  close(logSumExp([-1000, -1001]), -1000 + Math.log1p(Math.exp(-1)));
  close(logSumExp([1000, 1001]), 1001 + Math.log1p(Math.exp(-1)));
  assert.deepEqual([...normalizeLogWeights([1e16, 1e16])], [0.5, 0.5]);
  assert.deepEqual([...normalizeWeights([1e308, 1e308])], [0.5, 0.5]);
  relativeClose(logSumExp([0, -100]), Math.exp(-100), 1e-15);
  assert.equal(logSumExp([]), -Infinity); assert.equal(logSumExp([-Infinity]), -Infinity);
});

test('normalization rejects impossible support and invalid weights', () => {
  assert.throws(() => normalizeWeights([0, 0]), /positive/);
  assert.throws(() => normalizeWeights([1, -1]), /non-negative/);
  assert.throws(() => normalizeLogWeights([-Infinity, -Infinity]), /possible state/);
  assert.throws(() => logSumExp([NaN]), /finite/);
  assert.throws(() => normalizeLogWeights([Infinity, 1]), /finite/);
});

test('Shannon entropy and KL divergence reproduce exact distributions, zeros and near equality', () => {
  close(shannonEntropy([1, 0]), 0, 0); close(shannonEntropy([0.25, 0.25, 0.25, 0.25]), Math.log(4));
  close(shannonEntropy([0.25, 0.25, 0.25, 0.25], { base: 2 }), 2);
  close(relativeEntropy([0.25, 0.75], [0.5, 0.5]), 0.25 * Math.log(0.5) + 0.75 * Math.log(1.5));
  close(relativeEntropy([0.4, 0.6], [0.4, 0.6]), 0, 0);
  assert.equal(relativeEntropy([1, 0], [0, 1]), Infinity);
  const delta = 1e-8;
  relativeClose(relativeEntropy([0.5 + delta, 0.5 - delta], [0.5, 0.5]), 2 * delta * delta, 2e-8);
  assert.throws(() => shannonEntropy([0.2, 0.2]), /normalized/);
  assert.throws(() => relativeEntropy([1], [0.5, 0.5]), /lengths/);
  assert.throws(() => shannonEntropy([1], { base: 1 }), /base/);
});

test('joint covariance retains correlations and variance at large offsets', () => {
  const samples = [-2, -1, 1, 2].map(x => [1e12 + x, -1e12 + 2 * x + 3]);
  const result = covarianceMatrix(samples);
  close(result.mean[0], 1e12, 0); close(result.mean[1], -1e12 + 3, 0);
  close(result.covariance[0][0], 10 / 3); close(result.covariance[0][1], 20 / 3); close(result.covariance[1][1], 40 / 3);
  const population = covarianceMatrix(samples, { ddof: 0 });
  close(population.covariance[0][0] + 2 * population.covariance[0][1] + population.covariance[1][1], 22.5);
  // The mean rounds to a representable neighbour; centering about it would give variance 4.
  close(covarianceMatrix([[1e16], [1e16 + 2]]).covariance[0][0], 2, 1e-15);
  const rare = Array.from({ length: 100 }, (_, i) => [i === 0 ? 1e155 : 0]);
  relativeClose(covarianceMatrix(rare).covariance[0][0], 1e308, 5e-15);
  close(covarianceMatrix([[3]], { ddof: 0 }).covariance[0][0], 0, 0);
  assert.throws(() => covarianceMatrix([[3]]), /ddof/);
  assert.throws(() => covarianceMatrix([[1, 2], [1]]), /row samples/);
});

test('finite-spectrum canonical calculation matches exact binomial-level thermodynamics', () => {
  const epsilon = 2, kT = 0.8, g = [1, 2, 1], r = 1 / (1 + Math.exp(epsilon / kT));
  const s = canonicalEnsemble([0, epsilon, 2 * epsilon], { kT, degeneracies: g });
  [((1 - r) ** 2), 2 * r * (1 - r), r * r].forEach((value, i) => close(s.probabilities[i], value));
  close(s.meanEnergy, 2 * epsilon * r);
  close(s.energyVariance, 2 * epsilon * epsilon * r * (1 - r));
  close(s.heatCapacity, s.energyVariance / (kT * kT));
  close(s.entropy, -2 * (r * Math.log(r) + (1 - r) * Math.log1p(-r)));
  close(s.freeEnergy, -2 * kT * Math.log1p(Math.exp(-epsilon / kT)));
  close(s.entropy, (s.meanEnergy - s.freeEnergy) / kT);
  const probabilityFormula = -s.probabilities.reduce((sum, p, i) => sum + p * Math.log(p / g[i]), 0);
  close(s.entropy, probabilityFormula);
  const shifted = canonicalEnsemble([1e12, 1e12 + epsilon, 1e12 + 2 * epsilon], { kT, degeneracies: g });
  assert.deepEqual(shifted.probabilities, s.probabilities);
  close(shifted.entropy, s.entropy, 0); close(shifted.energyVariance, s.energyVariance, 0);
  close(shifted.freeEnergy - 1e12, s.freeEnergy, 1e-4);
});

test('canonical heat capacity equals the independently differentiated mean energy', () => {
  const energies = [-3, 0, 4, 8], g = [2, 1, 3, 1], kT = 1.7, h = 1e-5;
  const at = value => canonicalEnsemble(energies, { kT: value, degeneracies: g });
  const derivative = (at(kT + h).meanEnergy - at(kT - h).meanEnergy) / (2 * h);
  close(at(kT).heatCapacity, derivative, 1e-10);
  const freeDerivative = -(at(kT + h).freeEnergy - at(kT - h).freeEnergy) / (2 * h);
  close(at(kT).entropy, freeDerivative, 1e-10);
});

test('canonical low/high temperature limits retain degeneracy and tiny entropy', () => {
  const cold = canonicalEnsemble([0, 1], { kT: 0.01 });
  relativeClose(cold.entropy, 101 * Math.exp(-100), 2e-15);
  relativeClose(cold.freeEnergy, -0.01 * Math.exp(-100), 2e-15);
  const ground = canonicalEnsemble([3, 4], { kT: 1e-300, degeneracies: [2, 1] });
  close(ground.meanEnergy, 3, 0); close(ground.entropy, Math.log(2), 0); close(ground.energyVariance, 0, 0);
  const hot = canonicalEnsemble([0, 1, 2], { kT: 1e200, degeneracies: [1, 2, 1] });
  close(hot.entropy, Math.log(4)); close(hot.meanEnergy, 1); close(hot.energyVariance, 0.5);
  assert.throws(() => canonicalEnsemble([0, 1], { kT: 0 }), /positive/);
  assert.throws(() => canonicalEnsemble([0, 1], { kT: 1, degeneracies: [1, 0.5] }), /integer/);
  assert.throws(() => canonicalEnsemble([-1e308, 1e308], { kT: 1e308 }), /range/);
});
