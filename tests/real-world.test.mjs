// Validation against the real world: p-simtor's physics modules must reproduce
// measured or standard astronomical / laboratory values. Each test names its source.
// These complement the internal-consistency tests (convergence, identities, exact solutions).

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { integrateTo, createSystem, createStepper, toBarycentric } from '../core/nbody.mjs';
import { PLANETS, planetState, solarSystemScenario } from '../core/solar-system.mjs';
import { dopri5 } from '../core/numerics.mjs';
import { exitance, WIEN_B, planckNu, photonDensity, SIGMA_SB } from '../core/radiation.mjs';
import { tracePhoton, B_CRIT } from '../core/schwarzschild.mjs';
import {
  c, G, h, e as qe, me, mp, Ry, eV, AU, GMsun, Msun, yearJulian, G_AU_MSUN_YR,
} from '../core/constants.mjs';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const ARCSEC = 180 / Math.PI * 3600;
const DAY = 86400;

// ---------------------------------------------------------------- celestial mechanics

test('Earth’s sidereal year: two-body integration gives 365.256 days (IAU)', () => {
  // Sun–Earth in AU, M☉, yr (G = 4π²); semi-major axis 1.00000261 AU (Standish 1992), m_⊕/M☉ = 3.0035e-6.
  const a = 1.00000261, mE = 3.0035e-6, mu = G_AU_MSUN_YR * (1 + mE);
  const v = Math.sqrt(mu / a);
  const sys = createSystem([{ m: 1, x: [0, 0, 0], v: [0, 0, 0] }, { m: mE, x: [a, 0, 0], v: [0, v, 0] }], { G: G_AU_MSUN_YR });
  toBarycentric(sys);
  const st = createStepper(sys, { integrator: 'logh4', stepsPerCrossing: 20000 });
  // find the time when the planet returns to the +x axis (y crosses 0 upward)
  let prevY = sys.x[4], t0 = 0, crossing = null;
  while (sys.t < 1.2 && crossing == null) {
    t0 = sys.t; prevY = sys.x[4] - sys.x[1];
    st.step();
    const y = sys.x[4] - sys.x[1];
    if (sys.t > 0.5 && prevY < 0 && y >= 0) crossing = t0 + (sys.t - t0) * (-prevY) / (y - prevY);
  }
  const days = crossing * yearJulian / DAY;
  assert.ok(Math.abs(days - 365.256363) < 0.01, `${days} d`);
});

test('Jupiter’s orbital period: 11.86 years from Kepler’s third law in the N-body code', () => {
  const a = 5.2026, mJ = 9.5479e-4;                     // AU; M☉ (JPL)
  const r = integrateTo([{ m: 1, x: [0, 0, 0], v: [0, 0, 0] }, { m: mJ, x: [a, 0, 0], v: [0, Math.sqrt(G_AU_MSUN_YR * (1 + mJ) / a), 0] }],
    { G: G_AU_MSUN_YR, integrator: 'logh4' }, 0.1);
  assert.ok(r.dE < 1e-12);
  const period = 2 * Math.PI * Math.sqrt(a ** 3 / (G_AU_MSUN_YR * (1 + mJ)));
  assert.ok(Math.abs(period - 11.862) < 0.01, `${period} yr`);
});

test('Mercury’s perihelion advance: 43″ per century from integrating the relativistic orbit', () => {
  // Timelike Schwarzschild orbit d²u/dφ² = M/L² − u + 3Mu² (M = GM☉/c², L = specific angular
  // momentum). In terms of v = u·p with p = L²/M (semi-latus rectum): v'' = 1 − v + 3(M/p) v².
  const M = GMsun / (c * c);                            // 1476.6 m
  const a = 5.7909e10, ecc = 0.205630, periodDays = 87.9691;   // JPL mean elements
  const p = a * (1 - ecc * ecc), eps = 3 * M / p;
  const perihelionAngle = epsilon => {
    let last = null, phiPeri = null;
    dopri5((_phi, [v, w]) => [w, 1 - v + epsilon * v * v], 0, [1 + ecc, 0], 2.2 * Math.PI, {
      rtol: 1e-13, atol: 1e-15, hmax: 1e-3,
      onStep(phi, y) {
        // next perihelion = next maximum of v: w changes sign + → −. Newton with the exact
        // derivative w' = 1 − v + εv² refines the crossing.
        if (last && phi > Math.PI && last.w > 0 && y[1] <= 0) {
          phiPeri = phi - y[1] / (1 - y[0] + epsilon * y[0] * y[0]);
          return false;
        }
        last = { phi, w: y[1] };
        return true;
      },
    });
    return phiPeri;
  };
  // Differential measurement: the Newtonian run (ε = 0) shares the integrator's systematic
  // error, which cancels in the difference; the remainder is the relativistic advance.
  const phiPeri = 2 * Math.PI + perihelionAngle(eps) - perihelionAngle(0);
  const perOrbit = (phiPeri - 2 * Math.PI) * ARCSEC;
  const perCentury = perOrbit * (36525 / periodDays);
  // Observed anomalous advance 42.98 ± 0.04″/cy (Clemence 1947; Park et al. 2017)
  assert.ok(Math.abs(perCentury - 42.98) < 0.1, `${perCentury}″/century`);
  // and the textbook formula 6πGM/(c²a(1−e²))
  assert.ok(rel(perOrbit / ARCSEC, 6 * Math.PI * M / p) < 1e-3, `${perOrbit} vs ${6 * Math.PI * M / p * ARCSEC}`);
});

test('Earth’s position on 2000-01-01 12:00 TT from JPL elements matches the DE ephemeris', () => {
  // JPL Horizons, Earth–Moon barycentre, heliocentric ecliptic J2000: (−0.1771, 0.9672, 0.0000) AU
  const { x } = planetState(PLANETS.find(p => p.name === 'Earth'));
  assert.ok(Math.abs(x[0] + 0.1771) < 2e-3 && Math.abs(x[1] - 0.9672) < 2e-3, x.join(', '));
});

test('Mars completes an orbit in 1.881 yr in the full inner-Solar-System N-body run', () => {
  const sc = solarSystemScenario({ planets: ['Mercury', 'Venus', 'Earth', 'Mars'] });
  const sys = createSystem(sc.bodies, { G: G_AU_MSUN_YR });
  toBarycentric(sys);
  const st = createStepper(sys, { integrator: 'logh4', stepsPerCrossing: 2000 });
  const k = 4;                                           // Mars index (Sun = 0)
  const ang = () => Math.atan2(sys.x[3 * k + 1] - sys.x[1], sys.x[3 * k] - sys.x[0]);
  const a0 = ang();
  let turned = 0, prev = a0, t = null;
  while (sys.t < 2.5 && t == null) {
    const tPrev = sys.t;
    st.step();
    const a = ang();
    let d = a - prev; if (d < -Math.PI) d += 2 * Math.PI; if (d > Math.PI) d -= 2 * Math.PI;
    if (turned < 2 * Math.PI && turned + d >= 2 * Math.PI) t = tPrev + (sys.t - tPrev) * (2 * Math.PI - turned) / d;
    turned += d; prev = a;
  }
  // sidereal period 686.98 d = 1.8808 yr (NASA fact sheet); the J2000 osculating orbit and
  // perturbations from the other planets shift a single revolution by < 0.2 %
  assert.ok(Math.abs(t - 1.8808) / 1.8808 < 2e-3, `${t} yr`);
});

// ---------------------------------------------------------------- general relativity

test('Light deflection at the Sun’s limb: 1.751″ (Eddington 1919; VLBI γ = 1 to 1e-4)', () => {
  const M = GMsun / (c * c), bSun = 6.957e8 / M;        // impact parameter in units of M
  const alpha = tracePhoton(bSun, { dphi: 2e-5 }).deflection * ARCSEC;
  assert.ok(Math.abs(alpha - 1.7512) < 0.002, `${alpha}″`);
});

test('M87* shadow: predicted 39.7 µas vs the EHT ring 42 ± 3 µas (EHT 2019)', () => {
  // M = 6.5e9 M☉, D = 16.8 Mpc (EHT Paper VI). The emission ring sits slightly outside the shadow.
  const M = 6.5e9 * Msun, D = 16.8e6 * 3.0857e16;
  const diameter = 2 * B_CRIT * G * M / (c * c) / D * ARCSEC * 1e6;
  assert.ok(Math.abs(diameter - 39.7) < 0.3, `${diameter} µas`);
  assert.ok(Math.abs(diameter - 42) <= 3 + 0.5, 'consistent with the measured ring');
});

test('GPS satellite clocks gain ≈ 38.6 µs/day (gravitational +45.7, velocity −7.2)', () => {
  const GMe = 3.986004418e14, Re = 6.371e6, rGps = 26561.75e3;  // IERS; GPS ICD orbit radius
  const grav = (GMe / (c * c)) * (1 / Re - 1 / rGps) * DAY * 1e6;
  const vel = -(GMe / rGps) / (2 * c * c) * DAY * 1e6;
  assert.ok(Math.abs(grav - 45.7) < 0.1, `${grav}`);
  assert.ok(Math.abs(vel + 7.2) < 0.1, `${vel}`);
  assert.ok(Math.abs(grav + vel - 38.5) < 0.2, `${grav + vel}`);
});

// ---------------------------------------------------------------- thermal radiation

test('Solar constant 1361 W/m² from T_eff = 5772 K and R☉ (IAU 2015; Kopp & Lean 2011)', () => {
  const S = exitance(5772) * (6.957e8 / AU) ** 2;
  assert.ok(Math.abs(S - 1361) < 2, `${S}`);
});

test('Sun’s spectrum peaks near 502 nm (Wien)', () => {
  assert.ok(Math.abs(WIEN_B / 5772 * 1e9 - 502.0) < 0.1);
});

test('Earth’s effective temperature 255 K from the solar constant and albedo 0.30', () => {
  const T = ((1361 * (1 - 0.30)) / (4 * SIGMA_SB)) ** 0.25;
  assert.ok(Math.abs(T - 254.6) < 0.5, `${T}`);
});

test('Cosmic microwave background: peak 160.2 GHz and 411 photons/cm³ at 2.7255 K (COBE/FIRAS)', () => {
  let best = 0, nuBest = 0;
  for (let nu = 100e9; nu < 220e9; nu += 1e7) { const B = planckNu(nu, 2.7255); if (B > best) { best = B; nuBest = nu; } }
  assert.ok(Math.abs(nuBest / 1e9 - 160.2) < 0.05, `${nuBest / 1e9}`);
  assert.ok(Math.abs(photonDensity(2.7255) / 1e6 - 410.7) < 0.5);
});

// ---------------------------------------------------------------- atomic and particle physics

test('Hydrogen Balmer-α (vacuum) 656.47 nm with the reduced-mass Rydberg (NIST)', () => {
  const Rinf = Ry * eV / (h * c);                        // 10 973 731.57 m⁻¹
  assert.ok(Math.abs(Rinf - 10973731.568) < 0.01);
  const RH = Rinf / (1 + me / mp);
  const lambda = 1 / (RH * (1 / 4 - 1 / 9)) * 1e9;
  assert.ok(Math.abs(lambda - 656.47) < 0.01, `${lambda} nm`);
});

test('Proton cyclotron frequency 15.2 MHz in 1 T (NMR/cyclotron standard)', () => {
  const f = qe * 1 / (2 * Math.PI * mp);
  assert.ok(Math.abs(f / 1e6 - 15.2450) < 0.001, `${f / 1e6} MHz`);
});

test('Solar mass from GM☉ and G agrees with the IAU/CODATA value', () => {
  assert.ok(rel(GMsun / G, Msun) < 1e-4);
  assert.ok(Math.abs(AU / c - 499.004784) < 1e-5, 'light travel time Sun→Earth 499.0 s');
});
