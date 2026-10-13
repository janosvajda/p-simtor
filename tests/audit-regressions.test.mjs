// Regression tests for issues found in the independent math audit (October 2026).
// Each test reproduces the original failure and must keep passing.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { solveKepler, keplerRelative, keplerBodies } from '../core/kepler.mjs';
import { dopri5, besselI, linearFit, meanWithError, erf, lgamma, gamma, realYlm, integrateToInfinity, gaussLegendre, rng } from '../core/numerics.mjs';
import { createSystem, createStepper, drift, invariants } from '../core/nbody.mjs';
import { lagrangeTriangle } from '../core/presets.mjs';
import { planetState, PLANETS } from '../core/solar-system.mjs';
import { G_AU_MSUN_YR } from '../core/constants.mjs';
import { computeAccelerations } from '../core/nbody.mjs';

test('Kepler solver converges for negative mean anomaly at high eccentricity', () => {
  for (const e of [0.8, 0.99, 0.999999]) {
    for (let M = -20; M <= 20; M += 0.37) {
      const E = solveKepler(M, e);
      const Mr = ((M % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
      assert.ok(Math.abs(E - e * Math.sin(E) - Mr) < 1e-12, `e=${e} M=${M}: E=${E}`);
    }
  }
  // and the orbit at negative time is the mirror image of positive time
  const p = keplerRelative(1, 0.95, 1, 0.3), q = keplerRelative(1, 0.95, 1, -0.3);
  assert.ok(Math.abs(p.x - q.x) < 1e-12 && Math.abs(p.y + q.y) < 1e-12);
  assert.throws(() => solveKepler(1, 1));
  assert.throws(() => keplerBodies(1, 1, 1, 1));
});

test('dopri5 stops with an error on NaN instead of hanging', () => {
  assert.throws(() => dopri5((t, y) => [t > 0.5 ? NaN : y[0]], 0, [1], 1), /non-finite/);
});

test('besselI is correct for negative arguments: I_n(−x) = (−1)ⁿ I_n(x)', () => {
  assert.ok(Math.abs(besselI(1, -2) + 1.590636854637329) < 1e-14);
  assert.ok(Math.abs(besselI(2, -2) - 0.688948447698738) < 1e-14);
});

test('linearFit is exact with a large x offset; constant y is handled', () => {
  const xs = Array.from({ length: 100 }, (_, i) => 1e8 + i), ys = xs.map((_, i) => 3 + 2 * i);
  const f = linearFit(xs, ys);
  assert.ok(Math.abs(f.b - 2) < 1e-9 && Math.abs(f.r2 - 1) < 1e-12, JSON.stringify(f));
  assert.ok(Math.abs(f.a + 2 * 1e8 - 3) / 2e8 < 1e-12);
  const flat = linearFit([1, 2, 3], [5, 5, 5]);
  assert.equal(flat.b, 0); assert.equal(flat.r2, 1);
});

test('meanWithError never reports zero error for a short series', () => {
  const r = meanWithError([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
  assert.ok(r.naive && Math.abs(r.error - Math.sqrt(55 / 6 / 10)) < 1e-12);
  // long uncorrelated series: conservative but within ~25 % of the true standard error
  const g = rng(2), xs = Array.from({ length: 4096 }, () => g.normal());
  const m = meanWithError(xs);
  assert.ok(m.error > 0.95 / 64 && m.error < 1.25 / 64, `${m.error}`);
});

test('erf at ±∞ and far tails; Γ at its poles', () => {
  assert.equal(erf(Infinity), 1); assert.equal(erf(-Infinity), -1); assert.equal(erf(7), 1);
  assert.equal(lgamma(-1), Infinity); assert.equal(lgamma(0), Infinity);
  assert.ok(Number.isNaN(gamma(-2)));
  assert.ok(Math.abs(gamma(-0.5) + 2 * Math.sqrt(Math.PI)) < 1e-12, 'non-integer negatives still work');
});

test('real spherical harmonics stay finite and normalised at high l', () => {
  const v = realYlm(90, 90, 1, 0.3);
  assert.ok(Number.isFinite(v) && v !== 0, `${v}`);
  // φ part oscillates as cos²(30φ): needs enough Gauss points (256) to resolve it
  const norm = gaussLegendre(th => Math.sin(th) * gaussLegendre(ph => realYlm(40, 30, th, ph) ** 2, 0, 2 * Math.PI, 256), 0, Math.PI, 96);
  assert.ok(Math.abs(norm - 1) < 1e-9, `${norm}`);
});

test('integrateToInfinity keeps the far tail', () => {
  assert.ok(Math.abs(integrateToInfinity(x => 1 / (1 + x * x)) - Math.PI / 2) < 1e-13);
});

test('N-body: coincident bodies are rejected; marginally bound systems give a finite energy drift', () => {
  assert.throws(() => createStepper(createSystem([{ m: 1, x: [0, 0, 0], v: [0, 0, 0] }, { m: 1, x: [0, 0, 0], v: [1, 0, 0] }])), /not usable/);
  // parabolic two-body: E₀ = 0 exactly
  const s = createSystem([{ m: 1, x: [-0.5, 0, 0], v: [0, -1, 0] }, { m: 1, x: [0.5, 0, 0], v: [0, 1, 0] }]);
  const ref = invariants(s);
  assert.ok(Math.abs(ref.E) < 1e-15);
  assert.ok(Number.isFinite(drift(s, ref).dE));
});

test('lagrangeTriangle honours G', () => {
  const s = createSystem(lagrangeTriangle([1, 1, 1], 2, 4), { G: 4 });
  computeAccelerations(s);
  const omega2 = 4 * 3 / 8;
  for (let i = 0; i < 3; i++) assert.ok(Math.abs(s.a[3 * i] + omega2 * s.x[3 * i]) < 1e-12);
});

test('Solar System units: one Julian year of mean motion for Earth is consistent with G in AU–M☉–yr', () => {
  // Kepler with the osculating a = 1.00000261 AU: P = 1.0000174 yr (sidereal) × (1.00000261)^1.5
  // = 1.0000213 yr. The former G = 4π² gave 1.0000024 yr, 1.5e-5 short.
  const p = PLANETS.find(x => x.name === 'Earth');
  const P = 2 * Math.PI * Math.sqrt(p.a ** 3 / (G_AU_MSUN_YR * (1 + 1 / p.sunOverM)));
  assert.ok(Math.abs(P - 1.0000174 * p.a ** 1.5) < 1e-6, `${P}`);
  assert.ok(Number.isFinite(planetState(p).x[0]));
});

test('round-off eccentricity from a circular orbit is accepted (regression: e = −1e−16 crashed)', () => {
  // e from the eccentricity vector of an exactly circular orbit, as a page would compute it
  const mu = 1.3, r = [0.7, 0.4, 0], v0 = Math.sqrt(mu / Math.hypot(...r));
  const v = [-v0 * r[1] / Math.hypot(...r), v0 * r[0] / Math.hypot(...r), 0];
  const rn = Math.hypot(...r), v2 = v[0] ** 2 + v[1] ** 2, rv = r[0] * v[0] + r[1] * v[1];
  const evec = [0, 1].map(k => ((v2 - mu / rn) * r[k] - rv * v[k]) / mu);
  const e = Math.hypot(...evec);                                       // round-off only: ~1e−16
  assert.ok(e < 1e-12, `${e}`);
  for (const eps of [e, -e, -1e-16, 5e-13]) assert.ok(Number.isFinite(solveKepler(1.2, eps)), `${eps}`);
  assert.ok(Number.isFinite(keplerRelative(1, -3e-16, 1, 0.4).x));
  assert.equal(keplerBodies(1, 1, 1, -1e-16).length, 2);
  assert.throws(() => solveKepler(1, -0.01), /0 ≤ e < 1/);          // real negatives still rejected
  assert.throws(() => solveKepler(1, 1));
});
