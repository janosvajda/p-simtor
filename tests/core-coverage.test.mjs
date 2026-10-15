// Independent checks for every remaining public function in core/: each test uses a
// different route to the same answer (identities, symmetry, limits, reversibility, data).

import { test } from 'node:test';
import assert from 'node:assert/strict';
import * as nb from '../core/nbody.mjs';
import { keplerBodies, keplerRelative } from '../core/kepler.mjs';
import { lagrangeTriangle, THREE_BODY_PRESETS } from '../core/presets.mjs';
import { verifyScenario } from '../core/nbody-verify.mjs';
import { PLANETS, planetState, keplerPeriod, solarSystemScenarios, MOONS } from '../core/solar-system.mjs';
import { gaussLegendreNodes, factorial } from '../core/numerics.mjs';
import {
  cieXYZbar, spectrumToXYZ, chromaticity, xyzToLinearSRGB, srgbGamma, wavelengthColor, peakWavelength,
  WIEN_B, WIEN_X, ZETA3, c1L, c2, wienApprox, planckLambda,
} from '../core/radiation.mjs';
import {
  turningPoint, orbitalOmega, orbitalUt, orbitalEnergy, orbitalL, R_S, R_PHOTON, B_CRIT, M,
  deflectionExact, tracePhoton, shadowAngle,
} from '../core/schwarzschild.mjs';
import { createLattice, sweep, staple, plaquette, action, halfTr, qmul, wilsonLoop } from '../core/lattice-gauge.mjs';
import * as K from '../core/constants.mjs';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

// ---------------------------------------------------------------- constants (CODATA 2018 relations)

test('CODATA constants are mutually consistent', () => {
  assert.ok(rel(K.hbar, K.h / (2 * Math.PI)) < 1e-9);
  assert.ok(rel(K.eps0 * K.mu0 * K.c * K.c, 1) < 1e-9);
  assert.ok(rel(K.alpha, K.e * K.e / (4 * Math.PI * K.eps0 * K.hbar * K.c)) < 1e-9);
  assert.ok(rel(K.a0, K.hbar / (K.me * K.c * K.alpha)) < 1e-9);
  assert.ok(rel(K.Eh, K.alpha ** 2 * K.me * K.c ** 2) < 1e-9);
  assert.ok(rel(K.Ry * K.eV, K.Eh / 2) < 1e-9);
  assert.ok(rel(K.mpc2MeV, K.mp * K.c ** 2 / K.MeV) < 1e-8);
  assert.ok(rel(K.mec2MeV, K.me * K.c ** 2 / K.MeV) < 1e-8);
  assert.ok(rel(K.hbarcMeVfm, K.hbar * K.c / K.MeV * 1e15) < 1e-8);
  assert.ok(rel(K.NA * K.u, 1e-3) < 1e-9, 'molar mass constant 1 g/mol');
  assert.ok(K.mn > K.mp && rel(K.mn - K.mp, 2.30557435e-30) < 1e-6, 'neutron–proton mass difference');
  // Msun derives from the same G (no CODATA-edition mix); G in AU–M☉–yr is consistent with SI
  assert.ok(rel(K.G * K.Msun, K.GMsun) < 1e-15);
  assert.ok(rel(K.G_AU_MSUN_YR, K.GMsun * K.yearJulian ** 2 / K.AU ** 3) < 1e-15);
  assert.ok(rel(K.G_AU_MSUN_YR, 4 * Math.PI ** 2) < 1e-4, '≈ 4π² (Gaussian year)');
  assert.equal(K.hbar, K.h / (2 * Math.PI));
});

// ---------------------------------------------------------------- N-body internals

const three = THREE_BODY_PRESETS.pythagorean.bodies;

test('forces obey Newton’s third law and equal −∇(−U)', () => {
  const s = nb.createSystem(three);
  const U = nb.computeAccelerations(s);
  for (let k = 0; k < 3; k++) assert.ok(Math.abs(s.m[0] * s.a[k] + s.m[1] * s.a[3 + k] + s.m[2] * s.a[6 + k]) < 1e-13);
  assert.ok(rel(U, nb.forceFunction(s)) < 1e-15);
  // −∂V/∂x₀ with V = −U, by central differences
  const h = 1e-6, x0 = s.x[0];
  s.x[0] = x0 + h; const Up = nb.forceFunction(s);
  s.x[0] = x0 - h; const Um = nb.forceFunction(s);
  s.x[0] = x0;
  assert.ok(rel(s.m[0] * s.a[0], (Up - Um) / (2 * h)) < 1e-8);
});

test('kinetic energy, momentum, angular momentum and centre of mass by hand', () => {
  const s = nb.createSystem([{ m: 2, x: [1, 0, 0], v: [0, 3, 0] }, { m: 1, x: [0, 2, 0], v: [1, 0, 0] }]);
  assert.equal(nb.kineticEnergy(s), 0.5 * 2 * 9 + 0.5 * 1 * 1);
  assert.deepEqual(nb.totalMomentum(s), [1, 6, 0]);
  assert.deepEqual(nb.angularMomentum(s), [0, 0, 2 * (1 * 3) + 1 * (0 * 0 - 2 * 1)]);
  const c = nb.centerOfMass(s);
  assert.ok(rel(c.x[0], 2 / 3) < 1e-15 && rel(c.x[1], 2 / 3) < 1e-15 && c.M === 3);
  const inv = nb.invariants(s);
  assert.ok(rel(inv.E, inv.T - inv.U) < 1e-15);
  const d = nb.drift(s, inv);
  assert.ok(d.dE === 0 && d.dP === 0 && d.dL === 0);
});

test('cloneSystem is a deep copy', () => {
  const s = nb.createSystem(three), t = nb.cloneSystem(s);
  t.x[0] = 99;
  assert.notEqual(s.x[0], 99);
});

test('every integrator is time-reversible: h then −h returns to the start', () => {
  for (const [name, spec] of Object.entries(nb.INTEGRATORS)) {
    const s = nb.toBarycentric(nb.createSystem(THREE_BODY_PRESETS['figure-eight'].bodies));
    const ctx = { B: -nb.totalEnergy(s) };
    const x0 = s.x.slice(), h = name.startsWith('logh') ? 0.02 : 0.01;
    for (let i = 0; i < 50; i++) spec.step(s, h, ctx);
    for (let i = 0; i < 50; i++) spec.step(s, -h, ctx);
    assert.ok(nb.stateDistance(s.x, x0) < 1e-11, `${name}: ${nb.stateDistance(s.x, x0)}`);
    assert.ok(Math.abs(s.t) < 1e-12);
  }
});

test('characteristic times: crossing time, pair time, and the default step', () => {
  const k = nb.createSystem(keplerBodies(1, 1e-9, 1, 0));
  assert.ok(rel(2 * Math.PI * nb.pairTime(k), keplerRelative(1, 0, 1 + 1e-9, 0).period) < 1e-12, 'pair time = P/2π');
  // extreme mass ratio: the crossing time blows up but the default step must not
  assert.ok(nb.crossingTime(k) > 1e3 * nb.pairTime(k));
  const dt = nb.defaultStep(k, 'leapfrog', 1000);
  assert.ok(rel(dt, 2 * Math.PI * nb.pairTime(k) / 1000) < 1e-12);
  const hLog = nb.defaultStep(k, 'logh', 1000);
  assert.ok(rel(hLog, dt * nb.forceFunction(k)) < 1e-12, 'LogH: h = dt·U');
});

test('keplerBodies: barycentre at rest, pericentre distance a(1−e), vis-viva speed', () => {
  const [b1, b2] = keplerBodies(3, 1, 2, 0.4);
  for (let k = 0; k < 3; k++) {
    assert.ok(Math.abs(b1.m * b1.x[k] + b2.m * b2.x[k]) < 1e-14);
    assert.ok(Math.abs(b1.m * b1.v[k] + b2.m * b2.v[k]) < 1e-14);
  }
  assert.ok(rel(b2.x[0] - b1.x[0], 2 * 0.6) < 1e-15);
  assert.ok(rel(b2.v[1] - b1.v[1], Math.sqrt(4 * (2 / 1.2 - 1 / 2))) < 1e-14);
});

test('Lagrange triangle: equilateral, and every body accelerates centripetally at Ω²r', () => {
  const m = [1, 2, 3], bodies = lagrangeTriangle(m, 1.5);
  const s = nb.createSystem(bodies);
  nb.computeAccelerations(s);
  const d = (i, j) => Math.hypot(...[0, 1, 2].map(k => s.x[3 * i + k] - s.x[3 * j + k]));
  assert.ok(rel(d(0, 1), 1.5) < 1e-14 && rel(d(1, 2), 1.5) < 1e-14 && rel(d(0, 2), 1.5) < 1e-14);
  const omega2 = 6 / 1.5 ** 3;
  for (let i = 0; i < 3; i++) for (let k = 0; k < 2; k++) assert.ok(Math.abs(s.a[3 * i + k] + omega2 * s.x[3 * i + k]) < 1e-12);
});

test('verifyScenario produces a complete report', () => {
  const r = verifyScenario({ bodies: THREE_BODY_PRESETS['figure-eight'].bodies, duration: 2 });
  assert.ok(r.conservation.dE < 1e-10 && r.determinism.identical && r.frameIndependence.identical);
  assert.equal(r.benchmarks.length, 3);
  assert.ok(Math.abs(r.convergence.observedOrder - 4) < 0.3);
});

// ---------------------------------------------------------------- Solar System elements

test('planet state vectors: inclination and node recovered from r × v; Kepler periods (NASA)', () => {
  const nasa = { Mercury: 0.2408467, Venus: 0.61519726, Earth: 1.0000174, Mars: 1.8808476, Jupiter: 11.862615, Saturn: 29.447498, Uranus: 84.016846, Neptune: 164.79132 };
  for (const p of PLANETS) {
    const { x, v } = planetState(p);
    const L = [x[1] * v[2] - x[2] * v[1], x[2] * v[0] - x[0] * v[2], x[0] * v[1] - x[1] * v[0]];
    // atan2 keeps full precision for tiny inclinations (acos near 1 does not)
    const inc = Math.atan2(Math.hypot(L[0], L[1]), L[2]) * 180 / Math.PI;
    assert.ok(Math.abs(inc - Math.abs(p.I)) < 1e-9, `${p.name} inclination ${inc}`);
    if (Math.abs(p.I) > 0.01) {
      const node = (Math.atan2(L[0], -L[1]) * 180 / Math.PI + 360) % 360;    // ascending node from L
      assert.ok(Math.abs(node - p.O) < 1e-8, `${p.name} node ${node} vs ${p.O}`);
    }
    // osculating vs mean periods differ slightly for the outer giants
    assert.ok(rel(keplerPeriod(p), nasa[p.name]) < 0.003, `${p.name} ${keplerPeriod(p)}`);
  }
  const sc = solarSystemScenarios();
  assert.deepEqual(sc.map(s => s.key), ['inner', 'outer', 'all']);
  assert.equal(sc[2].bodies.length, 9);
});

// ---------------------------------------------------------------- numerics leftovers

test('Gauss–Legendre nodes: symmetric, inside (−1, 1), weights sum to 2; factorial', () => {
  for (const n of [1, 2, 5, 32]) {
    const { x, w } = gaussLegendreNodes(n);
    assert.ok(Math.abs(w.reduce((a, b) => a + b, 0) - 2) < 1e-13, `n=${n}`);
    for (let i = 0; i < n; i++) { assert.ok(Math.abs(x[i] + x[n - 1 - i]) < 1e-14); assert.ok(Math.abs(x[i]) < 1); }
  }
  assert.equal(factorial(0), 1);
  assert.equal(factorial(10), 3628800);
});

// ---------------------------------------------------------------- colour pipeline

test('CIE 1931 fit: equal-energy white lands at (1/3, 1/3); D65 white maps to sRGB (1, 1, 1)', () => {
  const [x, y] = chromaticity(spectrumToXYZ(() => 1));
  assert.ok(Math.abs(x - 1 / 3) < 0.004 && Math.abs(y - 1 / 3) < 0.004, `${x}, ${y}`);
  const rgb = xyzToLinearSRGB([0.95047, 1, 1.08883]);
  for (const v of rgb) assert.ok(Math.abs(v - 1) < 1e-3, `${rgb}`);
  // ȳ (luminous efficiency) peaks near 555 nm
  let best = 0, at = 0;
  for (let nm = 400; nm <= 700; nm++) if (cieXYZbar(nm)[1] > best) { best = cieXYZbar(nm)[1]; at = nm; }
  assert.ok(Math.abs(at - 556) < 6, `${at}`);
});

test('sRGB transfer curve: endpoints, continuity at the knee, monotonic', () => {
  assert.equal(srgbGamma(0), 0);
  assert.ok(Math.abs(srgbGamma(1) - 1) < 1e-12);
  const k = 0.0031308;
  assert.ok(Math.abs(srgbGamma(k - 1e-12) - srgbGamma(k + 1e-12)) < 1e-6);
  for (let v = 0; v < 1; v += 0.01) assert.ok(srgbGamma(v + 0.01) > srgbGamma(v));
});

test('spectral colours: violet/blue short, green middle, red long', () => {
  const [r1, , b1] = wavelengthColor(450), [r2, g2, b2] = wavelengthColor(530), [r3, , b3] = wavelengthColor(650);
  assert.ok(b1 > r1 && g2 > r2 && g2 > b2 && r3 > b3);
});

test('radiation helper constants', () => {
  assert.ok(rel(peakWavelength(5772), WIEN_B / 5772) < 1e-15);
  assert.ok(Math.abs(WIEN_X - 4.965114231744276) < 1e-13);
  assert.ok(Math.abs(ZETA3 - 1.2020569031595942) < 1e-16);
  assert.ok(rel(c1L, 2 * K.h * K.c ** 2) < 1e-15 && rel(c2, 1.438776877e-2) < 1e-9, 'c₂ (CODATA 2018)');
  assert.ok(rel(wienApprox(5e-7, 3000), planckLambda(5e-7, 3000)) < 1e-4, 'Wien limit at hc/λkT ≈ 9.6');
});

// ---------------------------------------------------------------- Schwarzschild leftovers

test('Schwarzschild edge cases: captured rays and observers inside the photon sphere', () => {
  assert.ok(Number.isNaN(turningPoint(5)) && Number.isNaN(turningPoint(B_CRIT)));
  assert.equal(deflectionExact(5), Infinity);
  // just above the critical value the quadrature still agrees with the ray tracer
  assert.ok(Math.abs(deflectionExact(5.2, 256) - tracePhoton(5.2, { dphi: 2e-4, maxPhi: 80 }).deflection) < 1e-6);
  // shadow angle: 90° at the photon sphere, > 90° inside it, undefined inside the horizon
  assert.ok(Math.abs(shadowAngle(3) - Math.PI / 2) < 1e-7);
  assert.ok(shadowAngle(2.5) > Math.PI / 2);
  assert.ok(Number.isNaN(shadowAngle(1.5)));
  // continuity across r₀ = 3M
  assert.ok(Math.abs(shadowAngle(3 - 1e-9) - shadowAngle(3 + 1e-9)) < 1e-3);
});

test('Schwarzschild: horizon, photon sphere, turning point, circular orbits', () => {
  assert.equal(R_S, 2 * M); assert.equal(R_PHOTON, 3 * M);
  assert.equal(-1 / R_PHOTON + 3 * M / R_PHOTON ** 2, 0, 'u = 1/3M is a fixed point of u″ = −u + 3Mu²');
  for (const b of [6, 10, 100]) {
    const u = turningPoint(b);
    assert.ok(Math.abs(1 / (b * b) - u * u + 2 * M * u ** 3) < 1e-15 && u < 1 / (3 * M));
  }
  assert.ok(rel(orbitalOmega(10), Math.sqrt(1 / 1000)) < 1e-15, 'Kepler’s law holds exactly in t');
  assert.ok(rel(orbitalUt(6), Math.sqrt(2)) < 1e-15);
  // E and L satisfy the effective-potential circular-orbit conditions: E² = (1−2M/r)(1+L²/r²)
  for (const r of [6, 10, 50]) assert.ok(rel(orbitalEnergy(r) ** 2, (1 - 2 / r) * (1 + orbitalL(r) ** 2 / r ** 2)) < 1e-14);
  assert.ok(B_CRIT > R_PHOTON);
});

// ---------------------------------------------------------------- lattice leftovers

test('staple: changing one link changes the action by −½Tr(ΔU·A)', () => {
  const lat = createLattice({ L: 3, d: 4, beta: 2, seed: 3, start: 'hot' });
  sweep(lat);
  const s = 5, mu = 2, o = 4 * (s * lat.D + mu);
  const A = staple(lat, s, mu);
  const before = action(lat), old = lat.U.slice(o, o + 4);
  const neu = [0.6, 0, 0.8, 0];
  lat.U.set(neu, o);
  const dS = action(lat) - before;
  const expected = -(halfTr(qmul(neu, A)) - halfTr(qmul(old, A)));
  assert.ok(Math.abs(dS - expected) < 1e-12, `${dS} vs ${expected}`);
  assert.ok(Math.abs(plaquette(lat)) <= 1);
  assert.throws(() => wilsonLoop(lat, 3, 1), /lattice extent/);
});

test('moon counts: every planet listed, Mercury and Venus have none', () => {
  for (const p of PLANETS) assert.ok(Number.isInteger(MOONS.count[p.name]), p.name);
  assert.equal(MOONS.count.Mercury, 0); assert.equal(MOONS.count.Venus, 0); assert.equal(MOONS.count.Earth, 1);
  assert.match(MOONS.asOf, /^\d{4}-\d{2}$/);
});

test('everyday constants are consistent with each other and with published values', async () => {
  const C = await import('../core/constants.mjs');
  assert.ok(Math.abs(C.R - C.kB * C.NA) < 1e-8, 'R = kB·NA');
  assert.ok(Math.abs(C.p0 * C.Mair / (C.R * C.T0isa) - C.rho0isa) < 5e-4, 'ISA density from the ideal gas law');
  assert.ok(Math.abs(Math.sqrt(C.gammaAir * C.R * C.T0isa / C.Mair) - 340.294) < 0.01, 'ISA speed of sound 340.294 m/s');
  assert.ok(Math.abs(C.GMmoon / C.Rmoon ** 2 - 1.624) < 0.002, 'Moon surface gravity 1.62 m/s²');
  assert.ok(Math.abs(C.GMearth / C.RearthMean ** 2 - C.gn) < 0.03, 'Earth surface gravity ≈ g_n (rotation and shape aside)');
});
