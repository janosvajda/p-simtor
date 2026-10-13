import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  createLattice, sweep, plaquette, wilsonLoop, creutzRatio, gaugeTransform, action, qmul, qdag,
} from '../core/lattice-gauge.mjs';
import { besselI, meanWithError } from '../core/numerics.mjs';

function measure(lat, therm, n, f) {
  for (let i = 0; i < therm; i++) sweep(lat);
  const xs = [];
  for (let i = 0; i < n; i++) { sweep(lat); xs.push(f(lat)); }
  return meanWithError(xs, 8);
}

test('SU(2) quaternion algebra: unitarity and associativity', () => {
  const a = [0.5, 0.5, 0.5, 0.5], b = [0.6, 0, 0.8, 0], c = [0, 0, 0, 1];
  const ab_c = qmul(qmul(a, b), c), a_bc = qmul(a, qmul(b, c));
  for (let i = 0; i < 4; i++) assert.ok(Math.abs(ab_c[i] - a_bc[i]) < 1e-15);
  const one = qmul(a, qdag(a));
  assert.ok(Math.abs(one[0] - 1) < 1e-15 && Math.abs(one[1]) < 1e-15);
});

test('2D SU(2) is exactly solvable: ⟨P⟩ = I₂(β)/I₁(β) and W(R,T) = ⟨P⟩^{RT}', () => {
  for (const beta of [1.5, 4]) {
    const lat = createLattice({ L: 16, d: 2, beta, seed: 11 });
    const u = besselI(2, beta) / besselI(1, beta);
    const P = measure(lat, 40, 160, plaquette);
    assert.ok(Math.abs(P.mean - u) < 4 * P.error + 1e-3, `β=${beta}: ${P.mean} ± ${P.error} vs ${u}`);
    const W = measure(lat, 0, 160, l => wilsonLoop(l, 2, 1));
    assert.ok(Math.abs(W.mean - u ** 2) < 4 * W.error + 2e-3, `W(2,1) ${W.mean} ± ${W.error} vs ${u ** 2}`);
  }
});

test('4D strong coupling: ⟨P⟩ = u + 4u⁵ + … resolves the 4u⁵ term at β = 1', () => {
  // u = I₂(β)/I₁(β); at β = 1, 4u⁵ = 3.2e-3, several standard errors: plain u is excluded.
  const beta = 1, u = besselI(2, beta) / besselI(1, beta);
  const lat = createLattice({ L: 4, d: 4, beta, seed: 5, start: 'hot' });
  const P = measure(lat, 30, 800, plaquette);
  assert.ok(Math.abs(P.mean - (u + 4 * u ** 5)) < 4 * P.error, `${P.mean} ± ${P.error} vs ${u + 4 * u ** 5}`);
  assert.ok(Math.abs(P.mean - u) > 3 * P.error, 'leading order alone is ruled out');
});

test('observables are gauge invariant', () => {
  const lat = createLattice({ L: 4, d: 4, beta: 2.3, seed: 2, start: 'hot' });
  for (let i = 0; i < 3; i++) sweep(lat);
  const P = plaquette(lat), W = wilsonLoop(lat, 2, 2);
  gaugeTransform(lat);
  assert.ok(Math.abs(plaquette(lat) - P) < 1e-12);
  assert.ok(Math.abs(wilsonLoop(lat, 2, 2) - W) < 1e-12);
});

test('over-relaxation is microcanonical: it leaves the action unchanged', () => {
  const lat = createLattice({ L: 4, d: 3, beta: 2.0, seed: 9, start: 'hot' });
  sweep(lat);
  const S = action(lat);
  sweep(lat, { heatbath: false, overrelax: 3 });
  assert.ok(Math.abs(action(lat) - S) / S < 1e-10, `${action(lat)} vs ${S}`);
});

test('confinement: Wilson loops obey an area law (positive string tension) in 3D', () => {
  const lat = createLattice({ L: 8, d: 3, beta: 3.0, seed: 4 });
  const W = {};
  for (let i = 0; i < 30; i++) sweep(lat);
  const acc = {};
  for (let n = 0; n < 40; n++) {
    sweep(lat);
    for (const [R, T] of [[1, 1], [2, 1], [1, 2], [2, 2]]) acc[`${R},${T}`] = (acc[`${R},${T}`] ?? 0) + wilsonLoop(lat, R, T) / 40;
  }
  Object.assign(W, acc);
  const chi = creutzRatio((R, T) => W[`${R},${T}`], 2, 2);
  assert.ok(chi > 0.02 && chi < 1, `χ(2,2) = ${chi}`);
  // perimeter law would give χ → 0; area law gives χ ≈ σa² > 0
});
