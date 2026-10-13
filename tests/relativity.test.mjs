import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  GM_EARTH, W0, L_G, GPS_A, GPS_F0, R_EARTH, OMEGA_EARTH, groundRate, orbitRateTerms, fractionalOffset,
  circularOffset, F_ICD, eccentricityTerm, MICRO_PER_DAY,
} from '../core/relativity.mjs';
import { c } from '../core/constants.mjs';
import { solveKepler } from '../core/kepler.mjs';
import { brent } from '../core/numerics.mjs';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const DAY = 86400;

test('ground clock: W₀/c² is the IAU 2000 defining constant L_G', () => {
  assert.ok(rel(W0 / (c * c), L_G) < 1e-9);
  assert.equal(groundRate(), 1 - W0 / (c * c));
  // W₀ is close to GM/R plus the centrifugal term at the equator (crude check of magnitude)
  const approx = GM_EARTH / R_EARTH + 0.5 * (OMEGA_EARTH * R_EARTH) ** 2;
  assert.ok(rel(W0, approx) < 2e-3, `${W0} vs ${approx}`);
});

test('orbit rate terms add up and match the fractional offset', () => {
  const r = 2e7, v2 = 1.5e7;
  const t = orbitRateTerms(r, v2);
  assert.equal(t.rate, 1 - t.grav - t.kin);
  // (rate − groundRate) subtracts two numbers ≈ 1: round-off ~1e-16 is the honest tolerance
  assert.ok(Math.abs(fractionalOffset(r, v2) - (t.rate - groundRate())) < 1e-15);
  assert.ok(Math.abs(fractionalOffset(GPS_A, GM_EARTH / GPS_A) - circularOffset(GPS_A)) < 1e-24);
});

test('GPS: +45.8 (gravity) − 7.2 (speed) = 38.6 µs/day; factory frequency 10.22999999543 MHz (IS-GPS-200)', () => {
  const grav = (W0 - GM_EARTH / GPS_A) / (c * c) * MICRO_PER_DAY;
  const kin = -GM_EARTH / GPS_A / 2 / (c * c) * MICRO_PER_DAY;
  assert.ok(Math.abs(grav - 45.79) < 0.01 && Math.abs(kin + 7.21) < 0.01, `${grav} ${kin}`);
  assert.ok(Math.abs(circularOffset(GPS_A) * MICRO_PER_DAY - 38.575) < 0.001);
  assert.ok(Math.abs(GPS_F0 * (1 - circularOffset(GPS_A)) - 10.22999999543e6) < 1e-4);
});

test('F constant: −2√μ/c² equals the IS-GPS-200 value with the ICD’s μ = 3.986005e14', () => {
  assert.ok(rel(-2 * Math.sqrt(3.986005e14) / (c * c), -4.442807633e-10) < 1e-9);
  assert.ok(rel(F_ICD, -2 * Math.sqrt(GM_EARTH) / (c * c)) < 1e-15);
});

test('periodic eccentricity term: integrating the clock rate along an exact Kepler orbit gives F e √a sin E', () => {
  const a = GPS_A, e = 0.02, n = Math.sqrt(GM_EARTH / a ** 3), T = 2 * Math.PI / n;
  // state on the Kepler orbit at time t (perigee at t = 0)
  const state = t => {
    const E = solveKepler(n * t, e), r = a * (1 - e * Math.cos(E));
    const v2 = GM_EARTH * (2 / r - 1 / a);                 // vis-viva
    return { E, r, v2 };
  };
  // Simpson integration of the periodic rate (offset − mean drift)
  const N = 4000, h = T / N;
  let acc = 0, worst = 0;
  const f = t => { const s = state(t); return fractionalOffset(s.r, s.v2) - circularOffset(a); };
  for (let k = 0; k < N; k += 2) {
    acc += h / 3 * (f(k * h) + 4 * f((k + 1) * h) + f((k + 2) * h));
    const E = state((k + 2) * h).E;
    worst = Math.max(worst, Math.abs(acc - eccentricityTerm(e, a, E)));
  }
  const amp = -F_ICD * e * Math.sqrt(a);
  assert.ok(Math.abs(amp * 1e9 - 45.79) < 0.05, `${amp * 1e9} ns`);
  assert.ok(worst < 1e-4 * amp, `${worst} vs ${amp}`);
});

test('ISS clocks run slow: about −25 µs/day, ~4.5 ms over a 6-month stay', () => {
  const a = R_EARTH + 400e3;
  const perDay = circularOffset(a) * MICRO_PER_DAY;
  assert.ok(perDay < -20 && perDay > -30, `${perDay}`);
  const halfYearMs = perDay * 182.6 / 1000;
  assert.ok(halfYearMs < -4 && halfYearMs > -5.5, `${halfYearMs} ms`);
});

test('gravity and speed effects cancel at a = 3GM/(2W₀) ≈ 9 545 km (altitude ≈ 3 170 km)', () => {
  const a = brent(circularOffset, R_EARTH + 1e3, GPS_A);
  assert.ok(Math.abs(a - 1.5 * GM_EARTH / W0) < 1e-3);
  assert.ok(Math.abs((a - R_EARTH) / 1e3 - 3166) < 10, `${(a - R_EARTH) / 1e3} km`);
});

test('GPS orbit: period is half a sidereal day', () => {
  const T = 2 * Math.PI * Math.sqrt(GPS_A ** 3 / GM_EARTH);
  assert.ok(Math.abs(T - 86164.0905 / 2) < 1, `${T}`);
  assert.ok(Math.abs(DAY - 86400) === 0);
});
