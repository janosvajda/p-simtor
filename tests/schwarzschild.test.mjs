import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  tracePhoton, deflectionExact, deflectionWeak, B_CRIT, R_ISCO, shadowAngle, impactFromAngle,
  orbitalEnergy, orbitalL, effectivePotential, diskRedshift, diskTemperature,
} from '../core/schwarzschild.mjs';

test('ray tracer agrees with exact quadrature deflection', () => {
  for (const b of [6, 8, 10, 30]) {
    const a = tracePhoton(b).deflection, e = deflectionExact(b);
    assert.ok(Math.abs(a - e) < 1e-8, `b=${b}: ${a} vs ${e}`);
  }
});

test('weak-field deflection → 4M/b + (15π/4)(M/b)² + (128/3)(M/b)³', () => {
  const b = 1000, a = tracePhoton(b).deflection;
  assert.ok(Math.abs(a - deflectionWeak(b)) / a < 1e-6);
  assert.ok(Math.abs(a - 4 / b) / a > 1e-3, 'higher orders are resolved');
});

test('capture threshold is b = 3√3 M (shadow edge)', () => {
  let lo = 5, hi = 5.5;
  for (let i = 0; i < 50; i++) { const m = (lo + hi) / 2; tracePhoton(m, { dphi: 2e-4, maxPhi: 80 }).captured ? (lo = m) : (hi = m); }
  assert.ok(Math.abs(lo - B_CRIT) < 1e-9, `${lo}`);
});

test('tracer converges at 4th order in dφ', () => {
  const ex = deflectionExact(8, 256);
  const e1 = Math.abs(tracePhoton(8, { dphi: 0.02 }).deflection - ex);
  const e2 = Math.abs(tracePhoton(8, { dphi: 0.01 }).deflection - ex);
  // the exit-angle interpolation is 2nd order, so the observed order sits between 2 and 4
  assert.ok(Math.log2(e1 / e2) > 1.8, `order ${Math.log2(e1 / e2)}`);
});

test('shadow angle formula is consistent with the impact-parameter map', () => {
  const ro = 30;
  assert.ok(Math.abs(impactFromAngle(ro, shadowAngle(ro)) - B_CRIT) < 1e-12);
});

test('ISCO at 6M: the effective-potential minimum disappears there', () => {
  // circular orbits: dV/dr = 0; stability requires d²V/dr² > 0. At r = 6M they coincide.
  const r = R_ISCO, L = orbitalL(r), h = 1e-3;
  const d2 = (effectivePotential(r + h, L) - 2 * effectivePotential(r, L) + effectivePotential(r - h, L)) / (h * h);
  assert.ok(Math.abs(d2) < 1e-6, `${d2}`);
  assert.ok(Math.abs(orbitalEnergy(6) - Math.sqrt(8 / 9)) < 1e-15, 'E_isco = √(8/9): 5.7% binding efficiency');
});

test('disk redshift: approaching side blueshifted, receding redshifted, gravitational redshift face-on', () => {
  const r = 10, L = 7;
  assert.ok(diskRedshift(r, +L) > diskRedshift(r, -L), 'co-rotating photons are blueshifted relative to counter-rotating');
  assert.ok(Math.abs(diskRedshift(r, 0) - Math.sqrt(1 - 3 / r)) < 1e-15);
  assert.equal(diskTemperature(5, 1e7), 0, 'no emission inside the ISCO');
  assert.ok(Math.abs(diskTemperature(6 * 49 / 36, 1e7) - 1e7) < 1e-6, 'peak temperature at r = (49/36) r_isco');
});

test('relativistic orbits: exact precession = orbit integration (strong field) = weak-field series (far away)', async () => {
  const { precessionPerOrbit, precessionWeak, relativisticOrbit, orbitConstants, precessionSI } = await import('../core/schwarzschild.mjs');
  for (const [p, e] of [[10, 0.5], [20, 0.2], [100, 0.7]]) {
    const o = relativisticOrbit(p, e, { orbits: 2 });
    assert.ok(Math.abs(o.periapsides[1] - 2 * Math.PI - precessionPerOrbit(p, e)) < 1e-7, `p=${p} e=${e}`);
  }
  for (const e of [0, 0.3, 0.9]) assert.ok(Math.abs(precessionPerOrbit(1e4, e) / precessionWeak(1e4, e) - 1) < 1e-6);
  // orbit constants put both turning points on the effective potential: E² = (1 − 2u)(1 + L²u²)
  const { L2, E2 } = orbitConstants(12, 0.4);
  for (const u of [1.4 / 12, 0.6 / 12]) assert.ok(Math.abs(E2 - (1 - 2 * u) * (1 + L2 * u * u)) < 1e-14);
  // separatrix p = 6 + 2e: the orbit whirls forever
  assert.equal(precessionPerOrbit(6.8, 0.4), Infinity);
  // Mercury: 42.98″ per century (JPL elements, IAU GM☉)
  const perOrbit = precessionSI({ GM: 1.3271244e20, a: 5.7909e10, e: 0.205630 });
  assert.ok(Math.abs(perOrbit * 206264.806 * 36525 / 87.9691 - 42.98) < 0.01);
});

test('relativistic orbit carries time: Newtonian period is Kepler’s, relativistic period → Kepler far away', async () => {
  const { relativisticOrbit } = await import('../core/schwarzschild.mjs');
  for (const [p, e] of [[50, 0.3], [400, 0.6]]) {
    const a = p / (1 - e * e), kepler = 2 * Math.PI * a ** 1.5;          // M = 1 units
    const n = relativisticOrbit(p, e, { orbits: 2, newtonian: true });
    assert.ok(Math.abs(n.periapsisTimes[1] / kepler - 1) < 1e-9, `Newton p=${p}: ${n.periapsisTimes[1]} vs ${kepler}`);
    const g = relativisticOrbit(p, e, { orbits: 2 });
    const rel = g.periapsisTimes[1] / kepler - 1;
    assert.ok(rel > 0 && rel < 20 / p, `GR p=${p}: ${rel}`);              // first-order correction O(M/p)
    // at(t): returns the periapsis distance at periapsis times, and positions on the orbit in between
    const q = g.at(g.periapsisTimes[1]);
    assert.ok(Math.abs(q.r - p / (1 + e)) / q.r < 1e-5);   // linear interpolation between samples: ~1e-5, far below a pixel
    const mid = g.at(g.periapsisTimes[1] / 2);
    assert.ok(mid.r > p / (1 + e) && mid.r <= p / (1 - e) * (1 + 1e-9));
  }
});
