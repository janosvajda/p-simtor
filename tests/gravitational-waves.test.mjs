import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  chirpMass, symmetricMassRatio, quadrupoleLuminosity, circularInspiral, chirpRateSI, timeToCoalescenceSI,
  gwFrequencySI, iscoFrequencySI, eccentricityEnhancement, periodDecay,
  quadrupoleTidalAmplitude, binaryTidalField, radiatedPowerFromField,
} from '../core/gravitational-waves.mjs';
import { GMsun, c } from '../core/constants.mjs';

const rel = (a, b) => Math.abs(a / b - 1);
const binaryQ = (mu, a) => { const A = mu * a * a / 2; return { re: [A, 0, 0, 0, -A, 0, 0, 0, 0], im: [0, A, 0, A, 0, 0, 0, 0, 0] }; };

test('chirp mass of GW150914 (36 + 29 M☉) is inside the published 28 (+1.8 −1.5) M☉', () => {
  const Mc = chirpMass(36, 29);
  assert.ok(Mc > 26.5 && Mc < 29.8, `${Mc}`);
  assert.equal(symmetricMassRatio(1, 1), 0.25);
});

test('Hulse–Taylor pulsar: GR period decay −2.4026e-12 (Weisberg & Huang 2016)', () => {
  const Pdot = periodDecay({ Pb: 0.322997448918 * 86400, e: 0.6171340, mp: 1.4398, mc: 1.3886 });
  assert.ok(rel(Pdot, -2.40263e-12) < 1e-3, `${Pdot}`);
  assert.equal(eccentricityEnhancement(0), 1);
});

test('circular inspiral solves da/dt = −(64/5) m₁m₂M/a³ and dφ/dt = Ω (centred differences, order 2)', () => {
  const ins = circularInspiral({ m1: 0.6, m2: 0.4, a0: 15 });
  const t = 0.4 * ins.tEnd, errs = [];
  for (const h of [4, 2, 1]) {
    const dadt = (ins.a(t + h) - ins.a(t - h)) / (2 * h), dphi = (ins.phi(t + h) - ins.phi(t - h)) / (2 * h);
    const a = ins.a(t);
    errs.push(Math.max(rel(dadt, -64 / 5 * 0.24 / a ** 3), rel(dphi, ins.omega(t))));
  }
  assert.ok(errs[2] < 1e-6, `${errs}`);
  const p = Math.log2(errs[1] / errs[2]);
  assert.ok(p > 1.7, `order ${p}`);
  assert.ok(Math.abs(ins.a(ins.tEnd) - 6) < 1e-9);
  // energy balance: d/dt(−m₁m₂/2a) = −L
  const dE = (0.24 / (2 * ins.a(t - 1)) - 0.24 / (2 * ins.a(t + 1))) / 2;
  assert.ok(rel(dE, -ins.luminosity(t)) < 1e-5, `${dE} vs ${-ins.luminosity(t)}`);
});

test('equal masses from 14 M to 6 M: 2900 M and 12.84 orbits', () => {
  const ins = circularInspiral({ a0: 14 });
  assert.ok(Math.abs(ins.tEnd - (14 ** 4 - 6 ** 4) / (4 * 64 / 5 * 0.25)) < 1e-9);
  assert.ok(Math.abs(ins.orbitsLeft(0) - 12.8355) < 1e-3, `${ins.orbitsLeft(0)}`);
});

test('SI chirp: time to coalescence agrees with the geometric-unit inspiral; df/dt = 3f/(8τ)', () => {
  const Msol = 65, m1 = 36 / 65, m2 = 29 / 65, Mc = chirpMass(36, 29);
  const f = 35, a = (GMsun * Msol / (Math.PI * f) ** 2) ** (1 / 3);
  assert.ok(rel(gwFrequencySI(a, Msol), f) < 1e-12);
  const aM = a / (GMsun * Msol / c ** 2), tsec = circularInspiral({ m1, m2, a0: aM }).tCoalesce * GMsun * Msol / c ** 3;
  assert.ok(rel(timeToCoalescenceSI(f, Mc), tsec) < 1e-9, `${timeToCoalescenceSI(f, Mc)} vs ${tsec}`);
  assert.ok(rel(chirpRateSI(f, Mc), 3 * f / (8 * timeToCoalescenceSI(f, Mc))) < 1e-12);
  // the ISCO frequency of 65 M☉ is about 68 Hz (the real signal peaked near 250 Hz, through the merger)
  assert.ok(Math.abs(iscoFrequencySI(65) - 67.6) < 0.2, `${iscoFrequencySI(65)}`);
});

test('tidal field: static limit equals ∂ᵢ∂ⱼ of the Newtonian quadrupole potential', () => {
  const Phi = (x, y, z) => -1.5 * (x * x - y * y) / Math.hypot(x, y, z) ** 5;
  const p = [1.3, 0.7, 0.4], h = 1e-3;
  const A = quadrupoleTidalAmplitude(...p, 1e-7, { re: [1, 0, 0, 0, -1, 0, 0, 0, 0], im: new Array(9).fill(0) });
  for (let i = 0; i < 3; i++) {for (let j = 0; j < 3; j++) {
    const at = (si, sj) => { const q = [...p]; q[i] += si * h; q[j] += sj * h; return Phi(...q); };
    const d2 = (at(1, 1) - at(1, -1) - at(-1, 1) + at(-1, -1)) / (4 * h * h);
    assert.ok(Math.abs(d2 - A.re[3 * i + j]) < 1e-5, `E${i}${j}: ${d2} vs ${A.re[3 * i + j]}`);
  }}
});

test('tidal field is traceless, symmetric and divergence-free in vacuum (residual falls at order 2)', () => {
  const Q = binaryQ(0.25, 10), w = 2 * 10 ** -1.5, p = [17, -9, 6];
  const E = (x, y, z) => quadrupoleTidalAmplitude(x, y, z, w, Q, { re: new Float64Array(9), im: new Float64Array(9) });
  const e = E(...p);
  const scale = Math.max(...e.re.map(Math.abs), ...e.im.map(Math.abs));
  assert.ok(Math.abs(e.re[0] + e.re[4] + e.re[8]) < 1e-12 * scale && Math.abs(e.im[0] + e.im[4] + e.im[8]) < 1e-12 * scale);
  for (const [i, j] of [[0, 1], [0, 2], [1, 2]]) assert.ok(Math.abs(e.re[3 * i + j] - e.re[3 * j + i]) < 1e-12 * scale);
  const residual = h => {
    let m = 0;
    for (let i = 0; i < 3; i++) {
      let s = 0;
      // the full field Re[Ẽ e^{iωr}] (the amplitude has the outgoing phase e^{iωr} factored out)
      const full = (q, k) => { const r = Math.hypot(...q), e2 = E(...q); return e2.re[k] * Math.cos(w * r) - e2.im[k] * Math.sin(w * r); };
      for (let j = 0; j < 3; j++) { const a = [...p], b = [...p]; a[j] += h; b[j] -= h; s += (full(a, 3 * i + j) - full(b, 3 * i + j)) / (2 * h); }
      m = Math.max(m, Math.abs(s));
    }
    return m / (scale / Math.hypot(...p));
  };
  const r1 = residual(0.4), r2 = residual(0.2);
  assert.ok(Math.log2(r1 / r2) > 1.7, `order ${Math.log2(r1 / r2)} (${r1}, ${r2})`);
});

test('energy carried by the computed field equals the orbit energy loss (quadrupole formula)', () => {
  for (const [m1, m2, a] of [[0.5, 0.5, 14], [0.8, 0.2, 9]]) {
    const mu = m1 * m2, w = 2 * a ** -1.5;
    const P = radiatedPowerFromField(w, binaryQ(mu, a));
    assert.ok(rel(P, quadrupoleLuminosity(a, m1, m2)) < 1e-7, `${P}`);
  }
});

test('binary tidal field far above the orbit equals −½ ḧ₊ of the strain at the retarded time', () => {
  const ins = circularInspiral({ a0: 12 });
  const t = 200, z = 4000, E1 = binaryTidalField(ins, t, 0, 0, z);
  // face-on: E_xx = −½ ḧ₊ = 2ω² · ½h₊ (h ∝ cos 2φ), compare with the strain at the same retarded time
  const h = ins.strain(t, z), w = 2 * ins.omega(t - z);
  assert.ok(Math.abs(E1[0] - 0.5 * w * w * h.plus) < 2e-3 * Math.abs(0.5 * w * w * 4 * ins.mu / (ins.a(t - z) * z)), `${E1[0]} vs ${0.5 * w * w * h.plus}`);
});
