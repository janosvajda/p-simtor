import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  createLattice, sweep, plaquette, wilsonLoop, wilsonLoops, polyakovLoop, plaquetteField, apeSmear,
  strongCouplingPlaquette, gaugeTransform, siteCoords, siteIndex, qmul, qmulInto, qdag,
} from '../core/lattice-gauge.mjs';
import { besselI, meanWithError } from '../core/numerics.mjs';

const rough = (L = 4, d = 4, seed = 5) => { const lat = createLattice({ L, d, beta: 2.2, start: 'hot', seed }); for (let i = 0; i < 3; i++) sweep(lat); return lat; };

test('qmulInto matches qmul, with daggers', () => {
  const a = [0.5, 0.5, 0.5, 0.5], b = [0.6, 0, 0.8, 0];
  const out = new Float64Array(4);
  for (const [da, db] of [[false, false], [true, false], [false, true], [true, true]]) {
    const ref = qmul(da ? qdag(a) : a, db ? qdag(b) : b);
    qmulInto(out, 0, a, 0, b, 0, da, db);
    for (let k = 0; k < 4; k++) assert.ok(Math.abs(out[k] - ref[k]) < 1e-16);
  }
});

test('wilsonLoops equals wilsonLoop for every R×T, also on an anisotropic lattice', () => {
  const lat = createLattice({ L: [5, 4, 6, 5], beta: 2.3, start: 'hot', seed: 3 });
  for (let i = 0; i < 2; i++) sweep(lat);
  const W = wilsonLoops(lat, { Rmax: 3, Tmax: 3 });
  for (let R = 1; R <= 3; R++) for (let T = 1; T <= 3; T++) assert.ok(Math.abs(W(R, T) - wilsonLoop(lat, R, T)) < 1e-14, `${R}×${T}`);
  assert.ok(Math.abs(W(1, 1) - plaquette(lat)) < 1e-14);
  assert.equal(W.values.length, 9);
  assert.throws(() => W(4, 1), /outside/);
  assert.throws(() => wilsonLoops(lat, { Rmax: 4 }), /wrap/);
});

test('wilsonLoops on chosen planes: spatial R × time T loops are gauge invariant', () => {
  const lat = rough();
  const planes = [[0, 3], [1, 3], [2, 3]];
  const before = wilsonLoops(lat, { Rmax: 2, Tmax: 3, planes }).values.slice();
  gaugeTransform(lat);
  const after = wilsonLoops(lat, { Rmax: 2, Tmax: 3, planes }).values;
  for (let i = 0; i < before.length; i++) assert.ok(Math.abs(before[i] - after[i]) < 1e-13);
});

test('polyakovLoop: gauge invariant per site; 1 on a cold lattice', () => {
  const cold = createLattice({ L: 4, d: 4, beta: 2 });
  assert.equal(polyakovLoop(cold).mean, 1);
  const lat = rough();
  const p = polyakovLoop(lat);
  assert.equal(p.values.length, 4 ** 3);
  gaugeTransform(lat);
  const q = polyakovLoop(lat);
  for (let i = 0; i < p.values.length; i++) assert.ok(Math.abs(p.values[i] - q.values[i]) < 1e-13);
});

test('plaquetteField: mean is 1 − ⟨P⟩, gauge invariant at every site and plane', () => {
  const lat = rough(4, 3);
  const f = plaquetteField(lat);
  assert.equal(f.planes.length, 3);
  const mean = f.values.reduce((a, b) => a + b, 0) / f.values.length;
  assert.ok(Math.abs(mean - (1 - plaquette(lat))) < 1e-14);
  gaugeTransform(lat);
  const g = plaquetteField(lat);
  for (let i = 0; i < f.values.length; i++) assert.ok(Math.abs(f.values[i] - g.values[i]) < 1e-13);
});

test('siteCoords and siteIndex are inverse (axis 0 fastest)', () => {
  const lat = createLattice({ L: [3, 4, 5], beta: 1 });
  for (let s = 0; s < lat.V; s++) assert.equal(siteIndex(lat, siteCoords(lat, s)), s);
  assert.deepEqual(siteCoords(lat, 1), [1, 0, 0]);
  assert.deepEqual(siteCoords(lat, 3), [0, 1, 0]);
});

test('apeSmear: α = 0 is the identity, links stay in SU(2), lat untouched, time links copied', () => {
  const lat = rough();
  const copy = Float64Array.from(lat.U);
  const same = apeSmear(lat, { alpha: 0, steps: 2 });
  for (let i = 0; i < copy.length; i++) assert.ok(Math.abs(same.U[i] - copy[i]) < 1e-15);
  const sm = apeSmear(lat, { alpha: 0.5, steps: 3 });
  assert.deepEqual(lat.U, copy);
  for (let i = 0; i < lat.V * lat.D; i++) {
    assert.ok(Math.abs(Math.hypot(sm.U[4 * i], sm.U[4 * i + 1], sm.U[4 * i + 2], sm.U[4 * i + 3]) - 1) < 1e-14);
    if (i % lat.D === lat.D - 1) for (let k = 0; k < 4; k++) assert.equal(sm.U[4 * i + k], lat.U[4 * i + k]);
  }
});

test('apeSmear: gauge covariant (smeared loops gauge invariant) and smoothing (spatial plaquette rises)', () => {
  const lat = rough();
  const spatial = l => wilsonLoops(l, { Rmax: 1, planes: [[0, 1], [0, 2], [1, 2]] })(1, 1);
  const p0 = spatial(lat);
  const sm = apeSmear(lat, { alpha: 0.5, steps: 4 });
  const p1 = spatial(sm);
  assert.ok(p1 > p0 + 0.1, `${p0} → ${p1}`);
  const W = wilsonLoops(sm, { Rmax: 2, Tmax: 2, planes: [[0, 1], [0, 2], [1, 2]] }).values.slice();
  gaugeTransform(lat);
  const W2 = wilsonLoops(apeSmear(lat, { alpha: 0.5, steps: 4 }), { Rmax: 2, Tmax: 2, planes: [[0, 1], [0, 2], [1, 2]] }).values;
  for (let i = 0; i < W.length; i++) assert.ok(Math.abs(W[i] - W2[i]) < 1e-12);
});

test('strongCouplingPlaquette: exact in 2D, matches 3D Monte Carlo at small β', () => {
  assert.ok(Math.abs(strongCouplingPlaquette(1.3, 2) - besselI(2, 1.3) / besselI(1, 1.3)) < 1e-15);
  const lat = createLattice({ L: 6, d: 3, beta: 0.5, seed: 9 });
  for (let i = 0; i < 20; i++) sweep(lat);
  const xs = []; for (let i = 0; i < 160; i++) { sweep(lat); xs.push(plaquette(lat)); }
  const P = meanWithError(xs, 8);
  assert.ok(Math.abs(P.mean - strongCouplingPlaquette(0.5, 3)) < 4 * P.error + 1e-3, `${P.mean} ± ${P.error}`);
});

import { haarSU2, qexp, qlog, centreTransform, landauGauge } from '../core/lattice-gauge.mjs';
import { rng } from '../core/numerics.mjs';
import { ksTest } from '../core/data.mjs';

test('qexp/qlog are inverse; haarSU2 has the Haar law of ½Tr: density (2/π)√(1−a²)', () => {
  const v = [0.3, -1.1, 0.7], back = qlog(qexp(v));
  for (let k = 0; k < 3; k++) assert.ok(Math.abs(back[k] - v[k]) < 1e-14);
  const r = rng(4), a0 = Array.from({ length: 4000 }, () => haarSU2(r)[0]);
  const cdf = a => 0.5 + (a * Math.sqrt(1 - a * a) + Math.asin(a)) / Math.PI;
  assert.ok(ksTest(a0, cdf).p > 0.01);
});

test('gaugeTransform with a given g returns it and keeps loops; unfixed links average to zero (Elitzur)', () => {
  const lat = rough();
  const W0 = wilsonLoops(lat, { Rmax: 2 }).values.slice();
  const r = rng(7), g = Array.from({ length: lat.V }, () => haarSU2(r));
  assert.equal(gaugeTransform(lat, g), g);
  const W1 = wilsonLoops(lat, { Rmax: 2 }).values;
  for (let i = 0; i < W0.length; i++) assert.ok(Math.abs(W0[i] - W1[i]) < 1e-13);
  let mean = 0; for (let i = 0; i < lat.V * lat.D; i++) mean += lat.U[4 * i];
  assert.ok(Math.abs(mean / (lat.V * lat.D)) < 0.05, `⟨½Tr U⟩ = ${mean / (lat.V * lat.D)}`);
});

test('centreTransform: action unchanged, Polyakov loops flip sign', () => {
  const lat = rough();
  const p0 = plaquette(lat), L0 = polyakovLoop(lat).values.slice();
  centreTransform(lat, { slice: 2 });
  assert.ok(Math.abs(plaquette(lat) - p0) < 1e-14);
  const L1 = polyakovLoop(lat).values;
  for (let i = 0; i < L0.length; i++) assert.ok(Math.abs(L1[i] + L0[i]) < 1e-14);
});

test('landauGauge: converges, raises the link functional, leaves every plaquette unchanged, g reproduces it', () => {
  const lat = createLattice({ L: 4, d: 4, beta: 2.4, start: 'hot', seed: 2 });
  for (let i = 0; i < 10; i++) sweep(lat);
  const original = Float64Array.from(lat.U), field = plaquetteField(lat).values.slice();
  const res = landauGauge(lat, { tol: 1e-12 });
  assert.ok(res.converged && res.theta < 1e-12, `θ = ${res.theta} after ${res.iterations}`);
  assert.ok(res.F > 0.5, `F = ${res.F}`);
  const after = plaquetteField(lat).values;
  for (let i = 0; i < field.length; i++) assert.ok(Math.abs(after[i] - field[i]) < 1e-12);
  const copy = { ...lat, U: original };
  gaugeTransform(copy, res.g);
  for (let i = 0; i < original.length; i++) assert.ok(Math.abs(copy.U[i] - lat.U[i]) < 1e-10);
});

import { topology, bpstInstanton } from '../core/lattice-gauge.mjs';

test('topology on an exact BPST instanton: Q = +1 within 2 % (improved clover beats plain), anti-instanton gives −Q', () => {
  const inst = bpstInstanton({ L: 12, rho: 3 }), anti = bpstInstanton({ L: 12, rho: 3, anti: true });
  const imp = topology(inst), plain = topology(inst, { improved: false }), a = topology(anti);
  assert.ok(Math.abs(imp.Q - 1) < 0.02, `Q = ${imp.Q}`);
  assert.ok(Math.abs(imp.Q - 1) < Math.abs(plain.Q - 1), `improved ${imp.Q} vs plain ${plain.Q}`);
  assert.ok(Math.abs(a.Q + imp.Q) < 1e-10, `anti ${a.Q}`);
  assert.ok(imp.S >= Math.abs(imp.Q) - 1e-9, 'S ≥ |Q|');
});

test('topology: zero on a cold lattice; densities gauge invariant site by site', () => {
  const cold = topology(createLattice({ L: 4, d: 4, beta: 2 }));
  assert.ok(Math.abs(cold.Q) < 1e-15 && Math.abs(cold.S) < 1e-15);
  const lat = rough(4, 4);
  const t0 = topology(lat);
  gaugeTransform(lat, rng(5));
  const t1 = topology(lat);
  for (let i = 0; i < lat.V; i++) assert.ok(Math.abs(t1.q[i] - t0.q[i]) < 1e-13 && Math.abs(t1.s[i] - t0.s[i]) < 1e-13);
});

test('cooling: action never rises, the instanton survives and becomes self-dual (S → |Q|)', () => {
  const lat = bpstInstanton({ L: 12, rho: 3 });
  let last = topology(lat).S;
  for (let k = 0; k < 6; k++) {
    sweep(lat, { cool: true });
    const t = topology(lat);
    assert.ok(t.S <= last + 1e-12, `S rose at cooling sweep ${k}`);
    last = t.S;
  }
  const t = topology(lat);
  assert.ok(Math.abs(t.Q - 1) < 0.02 && t.S / t.Q < 1.1, `Q = ${t.Q}, S = ${t.S}`);
});

test('sweep over site ranges: the sweep counter advances once; cooling in pieces equals one full cooling sweep', () => {
  for (const opts of [{}, { cool: true }]) {
    const a = createLattice({ L: 4, beta: 2.3, seed: 3, start: 'hot' }), b = createLattice({ L: 4, beta: 2.3, seed: 3, start: 'hot' });
    sweep(a, opts);
    sweep(b, { ...opts, sites: [0, 100] }); assert.equal(b.sweeps, 0);
    sweep(b, { ...opts, sites: [100, 256] }); assert.equal(b.sweeps, 1);
    if (!opts.cool) continue;   // heat bath + over-relaxation per range is a different (equally valid) update order
    assert.ok(a.U.every((v, i) => v === b.U[i]));
  }
});
