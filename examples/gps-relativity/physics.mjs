// GPS and relativity: a satellite clock and a ground clock, ticking at rates set by
// gravity (general relativity) and speed (special relativity).
// The orbit is integrated numerically (Yoshida 4th order, Earth as a point mass); the clock
// offset is integrated along it: d(τ_sat − τ_ground)/dt = (W₀ − GM/r − v²/2)/c².

import { createSystem, createStepper } from 'psim:core/nbody.mjs';
import {
  GM_EARTH, W0, L_G, GPS_A, GPS_F0, R_EARTH, fractionalOffset, circularOffset, orbitRateTerms,
  eccentricityTerm, MICRO_PER_DAY,
} from 'psim:core/relativity.mjs';
import { solveKepler } from 'psim:core/kepler.mjs';
import { c } from 'psim:core/constants.mjs';

export const HOUR = 3600;

export const PARAMS = [
  { key: 'a', label: 'Orbit radius (semi-major axis)', min: 6800, max: 60000, value: 26561.75, log: true, unit: 'km', digits: 6,
    hint: 'GPS: 26 561.75 km (two orbits per sidereal day). Low orbits run slow; high orbits run fast.' },
  { key: 'e', label: 'Eccentricity', min: 0, max: 0.05, value: 0.01, unit: '', digits: 3,
    hint: 'GPS satellites have e ≲ 0.02; eccentricity adds a periodic clock wobble of tens of nanoseconds.' },
  { key: 'dt', label: 'Integrator step', min: 1, max: 120, value: 10, unit: 's', digits: 3 },
];

export function createSim(p) {
  const a = p.a * 1e3, e = p.e;
  // start at perigee: r = a(1 − e), v = √(GM(1 + e)/(a(1 − e)))
  const rp = a * (1 - e), vp = Math.sqrt(GM_EARTH * (1 + e) / rp);
  const sys = createSystem([
    { m: 1, x: [0, 0, 0], v: [0, 0, 0] },
    { m: 1e-30, x: [rp, 0, 0], v: [0, vp, 0] },
  ], { G: GM_EARTH });
  const st = createStepper(sys, { integrator: 'yoshida4', h: p.dt });
  let offset = 0, grav = 0, kin = 0;               // accumulated seconds: total, gravity, velocity
  const rateNow = () => {
    const x = sys.x[3] - sys.x[0], y = sys.x[4] - sys.x[1], z = sys.x[5] - sys.x[2];
    const vx = sys.v[3] - sys.v[0], vy = sys.v[4] - sys.v[1], vz = sys.v[5] - sys.v[2];
    const r = Math.hypot(x, y, z), v2 = vx * vx + vy * vy + vz * vz;
    return { r, v2, total: fractionalOffset(r, v2), ...orbitRateTerms(r, v2) };
  };
  let last = rateNow();
  const history = [[0, 0, 0, 0]];                  // t (h), offset (µs), gravity part, velocity part
  const E0 = last.v2 / 2 - GM_EARTH / last.r;
  const period = 2 * Math.PI * Math.sqrt(a ** 3 / GM_EARTH);
  const n = 2 * Math.PI / period;
  const sim = {
    params: p, sys, a, e, period,
    get t() { return sys.t / HOUR; },
    step() {
      const t0 = sys.t;
      st.step();
      const now = rateNow(), h = sys.t - t0;
      // trapezoidal rule on the smooth rate (2nd order; error ≪ 1 ps per day at 10 s steps)
      offset += h * (last.total + now.total) / 2;
      grav += h * ((W0 / (c * c) - last.grav) + (W0 / (c * c) - now.grav)) / 2;
      kin -= h * (last.kin + now.kin) / 2;
      last = now;
      if (sys.steps % Math.max(1, Math.round(300 / p.dt)) === 0) history.push([sys.t / HOUR, offset * 1e6, grav * 1e6, kin * 1e6]);
    },
    offset: () => offset,
    parts: () => ({ grav, kin }),
    rate: () => last,
    history,
    /** Linear (secular) drift: what remains after removing the periodic eccentricity term. */
    secularRate: () => circularOffset(a),
    /**
     * Analytic periodic term at the current time. With vis-viva, the periodic rate is
     * −2GM(1/r − 1/a)/c² = −(2/c²) d(r·v)/dt, so the offset is −2 r·v/c² = F e √a sin E.
     */
    periodicAnalytic() {
      const E = solveKepler(n * sys.t, e);
      return eccentricityTerm(e, a, E);
    },
    position() { return [sys.x[3] - sys.x[0], sys.x[4] - sys.x[1]]; },
    diagnostics() {
      const r = rateNow(), En = r.v2 / 2 - GM_EARTH / r.r;
      return [
        { key: 'dE', label: 'Orbit energy |ΔE/E₀|', value: Math.abs((En - E0) / E0), kind: 'error', good: 1e-9, warn: 1e-6 },
        { key: 'off', label: 'Satellite clock ahead of ground clock', value: offset * 1e6, kind: 'value', unit: 'µs' },
        { key: 'err', label: 'Ranging error if not corrected (c × offset)', value: c * offset / 1000, kind: 'value', unit: 'km' },
        { key: 'rate', label: 'Average drift', value: circularOffset(a) * MICRO_PER_DAY, kind: 'value', unit: 'µs/day' },
        { key: 'T', label: 'Orbital period', value: period / HOUR, kind: 'value', unit: 'h' },
        { key: 'f', label: 'Factory clock setting for 10.23 MHz in orbit', value: 0, kind: 'value', text: `${(GPS_F0 * (1 - circularOffset(a)) / 1e6).toFixed(11)} MHz` },
      ];
    },
    state: () => ({ t: sys.t, offset, x: Array.from(sys.x.subarray(3, 6)) }),
    record: () => ({ header: ['t_hours', 'offset_us', 'gravity_part_us', 'velocity_part_us'], rows: history }),
  };
  return sim;
}

export function verify() {
  const checks = [];
  checks.push({ label: 'Ground clock rate W₀/c² vs IAU defining constant L_G', value: `${(W0 / (c * c)).toExponential(9)} vs ${L_G.toExponential(9)}`, pass: Math.abs(W0 / (c * c) / L_G - 1) < 1e-9 });
  const off = circularOffset(GPS_A);
  const f = GPS_F0 * (1 - off);
  checks.push({ label: 'Factory frequency offset vs GPS ICD 10.22999999543 MHz', value: `${(f / 1e6).toFixed(11)} MHz`, pass: Math.abs(f - 10.22999999543e6) < 1e-4 });
  checks.push({ label: 'Net drift of a GPS clock (textbook ≈ 38.6 µs/day)', value: `${(off * MICRO_PER_DAY).toFixed(3)} µs/day`, pass: Math.abs(off * MICRO_PER_DAY - 38.6) < 0.1 });
  // ICD constant uses the older μ = 3.986005e14 m³/s²
  const Ficd = -2 * Math.sqrt(3.986005e14) / (c * c);
  checks.push({ label: 'Relativistic constant F vs IS-GPS-200 (−4.442807633×10⁻¹⁰ s/√m, μ = 3.986005×10¹⁴)', value: Ficd.toExponential(9), pass: Math.abs(Ficd / -4.442807633e-10 - 1) < 1e-9 });
  // numerical integration over one day vs the analytic drift (circular orbit)
  const day = createSim({ a: GPS_A / 1e3, e: 0, dt: 10 });
  while (day.sys.t < 86400 - 1e-9) day.step();
  const rel = Math.abs(day.offset() / (off * day.sys.t) - 1);
  checks.push({ label: 'Integrated clock offset over 1 day vs (W₀ − 3GM/2a)/c² × t', value: `${(day.offset() * 1e6).toFixed(6)} µs (rel. diff ${rel.toExponential(1)})`, pass: rel < 1e-6 });
  // eccentric orbit: periodic part vs F e √a sin E
  const ecc = createSim({ a: GPS_A / 1e3, e: 0.02, dt: 5 });
  let worst = 0, amp = 0;
  while (ecc.sys.t < ecc.period) {
    ecc.step();
    const periodic = ecc.offset() - circularOffset(ecc.a) * ecc.sys.t;
    const analytic = ecc.periodicAnalytic();
    worst = Math.max(worst, Math.abs(periodic - analytic)); amp = Math.max(amp, Math.abs(analytic));
  }
  checks.push({ label: 'Eccentric orbit (e = 0.02): periodic clock term vs F e √a sin E (IS-GPS-200)', value: `amplitude ${(amp * 1e9).toFixed(2)} ns, max diff ${(worst * 1e12).toFixed(1)} ps`, pass: worst < 0.01 * amp });
  const T = day.period;
  checks.push({ label: 'GPS orbital period vs half a sidereal day (43 082.05 s)', value: `${T.toFixed(2)} s`, pass: Math.abs(T - 43082.045) < 1 });
  return {
    summary: 'Clock rates, orbit and GPS constants checked against IAU, IERS and the GPS interface specification (IS-GPS-200):',
    checks,
    notes: ['Earth is a point mass here; its oblateness (J₂) changes GPS clock rates by ~10⁻¹³, far below the 4.5×10⁻¹⁰ effect shown.'],
  };
}

export { R_EARTH, GM_EARTH, W0, circularOffset, MICRO_PER_DAY };
