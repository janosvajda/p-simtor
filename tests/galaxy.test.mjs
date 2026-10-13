import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  G, TIME_UNIT_MYR, R0, V0, MW, normalise, potential, acceleration, circularSpeed, haloMass, haloDensity,
  NFW_SCALE, eilers2019, createDisk, stepDisk, energies, angularMomenta,
} from '../core/galaxy.mjs';
import { G as Gsi, Msun } from '../core/constants.mjs';
import { gaussLegendre, rng, brent } from '../core/numerics.mjs';

const KPC = 3.085677581491367e19;            // m (IAU 2012: 1 pc = 648000/π AU)
const MYR = 1e6 * 365.25 * 86400;
const all = { bulge: true, disk: true, halo: true };
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('G in kpc (km/s)² / M☉ and the time unit kpc/(km/s) follow from SI constants', () => {
  assert.ok(rel(G, Gsi * Msun / KPC / 1e6) < 2e-4, `${G} vs ${Gsi * Msun / KPC / 1e6}`);
  assert.ok(rel(TIME_UNIT_MYR, KPC / 1e3 / MYR) < 1e-8);
});

test('MWPotential2014 normalisation: v_c(8 kpc) = 220 km/s with 5/60/35 % force fractions', () => {
  const v = circularSpeed(R0);
  assert.ok(Math.abs(v.total - V0) < 1e-10);
  assert.ok(Math.abs(v.bulge ** 2 / V0 ** 2 - 0.05) < 1e-12);
  assert.ok(Math.abs(v.disk ** 2 / V0 ** 2 - 0.60) < 1e-12);
  assert.ok(Math.abs(v.halo ** 2 / V0 ** 2 - 0.35) < 1e-12);
  const P = normalise({ fBulge: 0, fDisk: 1, fHalo: 0 });
  assert.equal(P.Mb, 0);
  // published masses (Bovy 2015, Table 1): disk 6.8e10 M☉, bulge 0.5e10 M☉
  assert.ok(rel(MW.Md, 6.8e10) < 0.02, `${MW.Md}`);
  assert.ok(rel(MW.Mb, 0.5e10) < 0.1, `${MW.Mb}`);
});

test('analytic accelerations equal −∇Φ everywhere, including off the plane', () => {
  const a = new Float64Array(3), h = 1e-5;
  const P = (x, y, z) => potential(Math.hypot(x, y), z, all);
  for (const [x, y, z] of [[8, 0, 0], [3, -4, 1.5], [0.3, 0.2, -0.05], [20, 15, 8], [1e-3, 0, 2]]) {
    acceleration(x, y, z, all, MW, a);
    const g = [-(P(x + h, y, z) - P(x - h, y, z)) / (2 * h), -(P(x, y + h, z) - P(x, y - h, z)) / (2 * h), -(P(x, y, z + h) - P(x, y, z - h)) / (2 * h)];
    const n = Math.hypot(...g);
    for (let k = 0; k < 3; k++) assert.ok(Math.abs(a[k] - g[k]) / n < 1e-6, `(${x},${y},${z}) component ${k}: ${a[k]} vs ${g[k]}`);
  }
});

test('the centre is regular: finite potential and zero force at r = 0', () => {
  const phi0 = potential(0, 0, all), phiEps = potential(1e-9, 0, all);
  assert.ok(Number.isFinite(phi0) && Math.abs(phi0 - phiEps) / Math.abs(phi0) < 1e-6, `${phi0} vs ${phiEps}`);
  const a = acceleration(0, 0, 0, all, MW, new Float64Array(3));
  assert.deepEqual(Array.from(a), [0, 0, 0]);
});

test('circular speed equals √(R·|a_R|) from the acceleration (midplane)', () => {
  const a = new Float64Array(3);
  for (const R of [0.5, 3, 8, 15, 30]) {
    acceleration(R, 0, 0, all, MW, a);
    assert.ok(rel(circularSpeed(R).total, Math.sqrt(-R * a[0])) < 1e-12, `R=${R}`);
  }
});

test('Hernquist bulge: v² = GMr/(r + a)² (closed form, a = 0.6 kpc)', () => {
  for (const r of [0.2, 1, 5]) assert.ok(rel(circularSpeed(r).bulge ** 2, G * MW.Mb * r / (r + 0.6) ** 2) < 1e-12);
});

test('NFW: density integrates to the enclosed-mass formula, and v_c peaks at 2.163 r_s', () => {
  for (const r of [1, 16, 100]) {
    const M = gaussLegendre(s => 4 * Math.PI * s * s * haloDensity(s), 1e-12, r, 96);
    assert.ok(rel(M, haloMass(r)) < 1e-7, `r=${r}`);
  }
  const vh = r => circularSpeed(r).halo;
  const rPeak = brent(r => (vh(r * 1.0001) - vh(r * 0.9999)), 10, 80, 1e-10);
  assert.ok(Math.abs(rPeak / NFW_SCALE - 2.16258) < 1e-3, `${rPeak / NFW_SCALE}`);   // Navarro et al. 1997
});

test('without the halo the outer rotation curve falls toward Keplerian, with it stays flat', () => {
  const noDM = { bulge: true, disk: true, halo: false };
  const ratioNo = circularSpeed(25, noDM).total / circularSpeed(10, noDM).total;
  const ratioDM = circularSpeed(25).total / circularSpeed(10).total;
  assert.ok(ratioNo < 0.75 && ratioDM > 0.85, `${ratioNo} ${ratioDM}`);
  // far out, visible matter is a point mass: v → √(G M_vis / R)
  const Rfar = 500, vis = MW.Mb + MW.Md;
  assert.ok(rel(circularSpeed(Rfar, noDM).total, Math.sqrt(G * vis / Rfar)) < 1e-3);
});

test('Eilers et al. (2019) fit: 229.0 km/s at R = 8.122 kpc, slope −1.7 km/s/kpc', () => {
  assert.equal(eilers2019(8.122), 229);
  assert.ok(Math.abs(eilers2019(18.122) - 212) < 1e-9);
});

test('star disk: exponential sampling, circular start, exact L_z and 2nd-order energy', () => {
  const d = createDisk({ n: 4000, comps: all, rand: rng(5), dispersion: 0 });
  // truncated Gamma(2, 2.6) on [1, 25]: mean ≈ 5.4 kpc (analytic: ∫R·R e^{−R/Rd} / ∫R e^{−R/Rd})
  const num = gaussLegendre(R => R * R * Math.exp(-R / 2.6), 1, 25, 64), den = gaussLegendre(R => R * Math.exp(-R / 2.6), 1, 25, 64);
  let mean = 0;
  for (let i = 0; i < d.n; i++) mean += d.R0s[i] / d.n;
  assert.ok(Math.abs(mean - num / den) < 0.1, `${mean} vs ${num / den}`);
  // with zero dispersion every star starts on the circular speed
  const i = 7, R = d.R0s[i], v = Math.hypot(d.v[3 * i], d.v[3 * i + 1]);
  assert.ok(rel(v, circularSpeed(R).total) < 1e-12);
  const err = dt => {
    const s = createDisk({ n: 40, comps: all, rand: rng(9) });
    for (let k = 0; k < Math.round(0.2 / dt); k++) stepDisk(s, dt);
    const E = energies(s), L = angularMomenta(s);
    let e = 0, l = 0;
    for (let j = 0; j < s.n; j++) { e = Math.max(e, Math.abs((E[j] - s.E0[j]) / s.E0[j])); l = Math.max(l, Math.abs(L[j] - s.L0[j]) / Math.abs(s.L0[j])); }
    return { e, l };
  };
  const a = err(1e-3), b = err(5e-4);
  assert.ok(Math.abs(Math.log2(a.e / b.e) - 2) < 0.4, `order ${Math.log2(a.e / b.e)}`);
  assert.ok(a.l < 1e-12 && b.l < 1e-12);
});
