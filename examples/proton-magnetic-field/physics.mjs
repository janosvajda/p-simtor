// Proton in uniform static fields B = B ẑ, E = E ŷ (SI units).
// Relativistic Boris pusher (Boris 1970; Birdsall & Langdon, "Plasma Physics via
// Computer Simulation", §4-3, 4-7). Time exposed to the page in nanoseconds.
//
// Exact references used for verification:
//   gyroradius   r = γ m v⊥ / (|q| B)          (E = 0)
//   gyroperiod   T = 2π γ m / (|q| B)          (E = 0)
//   E×B drift    v_d = E / B                    (E ≪ cB, guiding-centre average)
// Boris rotates the velocity by θ = 2 atan(Ω Δt / 2) per step instead of Ω Δt, so its
// discrete orbit is a circle of radius r·√(1 + (ΩΔt/2)²): a known, testable error.

import { c, e as qe, mp, MeV } from 'psim:core/constants.mjs';

export const PARAMS = [
  { key: 'B', label: 'Magnetic field B', min: 0.01, max: 10, value: 1, log: true, unit: 'T' },
  { key: 'KE', label: 'Kinetic energy', min: 0.001, max: 10000, value: 1, log: true, unit: 'MeV',
    hint: 'Above ~100 MeV the proton is relativistic: γ grows and the gyroperiod lengthens.' },
  { key: 'pitch', label: 'Pitch angle to B', min: 1, max: 90, value: 70, unit: '°', digits: 2 },
  { key: 'E', label: 'Electric field E (along y)', min: 0, max: 1e6, value: 0, unit: 'V/m', digits: 3,
    hint: 'Non-zero E adds an E×B drift along x at v = E/B.' },
  { key: 'stepsPerGyration', label: 'Steps per gyration', min: 8, max: 4096, value: 256, log: true, digits: 2 },
];

const NS = 1e9;

export function analytic({ B, KE, pitch }) {
  const gamma = 1 + KE * MeV / (mp * c * c);
  const beta = Math.sqrt(1 - 1 / (gamma * gamma));
  const v = beta * c;
  const vperp = v * Math.sin(pitch * Math.PI / 180);
  const vpar = v * Math.cos(pitch * Math.PI / 180);
  const omega = qe * B / (gamma * mp);          // relativistic gyrofrequency, rad/s
  return { gamma, beta, v, vperp, vpar, omega, period: 2 * Math.PI / omega, radius: gamma * mp * vperp / (qe * B) };
}

/** Create the simulation. State in SI; sim.t exposed in ns. */
export function createSim(p) {
  const a = analytic(p);
  const dt = a.period / p.stepsPerGyration;
  const qm = qe / mp;
  // q > 0 in B = +ẑ gyrates clockwise seen from +z: at (0, −r) it moves in −x, so the
  // Lorentz force q v × B points to +y, towards the orbit centre at the origin.
  // Boris is a leapfrog: u lives at half steps, rotating by θ = 2·atan(ΩΔt/2) per step on a
  // discrete circle of radius r·√(1 + (ΩΔt/2)²). Starting x on that circle and u half a step
  // back (rotated by +θ/2) centres the E = 0 orbit exactly on the origin.
  const half = a.omega * dt / 2, theta = 2 * Math.atan(half);
  const x = [0, -a.radius * Math.sqrt(1 + half * half), 0];
  const ux0 = -a.gamma * a.vperp;
  const u = [ux0 * Math.cos(theta / 2), ux0 * Math.sin(theta / 2), a.gamma * a.vpar];  // u = γv at t = −Δt/2
  const E = [0, p.E, 0], Bv = [0, 0, p.B];
  let tS = 0, steps = 0;
  const gamma0 = a.gamma;
  const trail = [[...x]];
  const rows = [];

  function gammaOf(w) { return Math.sqrt(1 + (w[0] ** 2 + w[1] ** 2 + w[2] ** 2) / (c * c)); }

  function step() {
    const h = qm * dt / 2;
    const um = [u[0] + h * E[0], u[1] + h * E[1], u[2] + h * E[2]];
    const g = gammaOf(um);
    const t = [h * Bv[0] / g, h * Bv[1] / g, h * Bv[2] / g];
    const t2 = t[0] ** 2 + t[1] ** 2 + t[2] ** 2;
    const up = [um[0] + um[1] * t[2] - um[2] * t[1], um[1] + um[2] * t[0] - um[0] * t[2], um[2] + um[0] * t[1] - um[1] * t[0]];
    const s = t.map(q => 2 * q / (1 + t2));
    const uplus = [um[0] + up[1] * s[2] - up[2] * s[1], um[1] + up[2] * s[0] - up[0] * s[2], um[2] + up[0] * s[1] - up[1] * s[0]];
    for (let k = 0; k < 3; k++) u[k] = uplus[k] + h * E[k];
    const gn = gammaOf(u);
    for (let k = 0; k < 3; k++) x[k] += u[k] / gn * dt;
    tS += dt; steps++;
    if (steps % Math.max(1, Math.round(p.stepsPerGyration / 128)) === 0) {
      trail.push([...x]);
      if (trail.length > 6000) trail.splice(0, trail.length - 6000);
    }
    if (rows.length < 100000 && steps % Math.max(1, Math.round(p.stepsPerGyration / 64)) === 0) {
      rows.push([tS * NS, ...x, ...u.map(q => q / gn)]);
    }
  }

  return {
    get t() { return tS * NS; },
    dt, analytic: a, x, u, trail, params: p,
    step,
    gamma: () => gammaOf(u),
    diagnostics() {
      const g = gammaOf(u);
      const list = [];
      if (!p.E) list.push({ key: 'dg', label: '|Δγ/γ₀| (energy)', value: Math.abs(g - gamma0) / gamma0, kind: 'error', good: 1e-12, warn: 1e-9 });
      list.push(
        { key: 'r', label: 'Gyroradius r = γmv⊥/(qB)', value: a.radius, kind: 'value', unit: 'm' },
        { key: 'T', label: 'Gyroperiod T = 2πγm/(qB)', value: a.period * NS, kind: 'value', unit: 'ns' },
        { key: 'g', label: 'Lorentz factor γ', value: g, kind: 'value' },
        { key: 'beta', label: 'Speed v/c', value: Math.sqrt(1 - 1 / (g * g)), kind: 'value' },
      );
      if (p.E) list.push({ key: 'vd', label: 'E×B drift E/B', value: p.E / p.B, kind: 'value', unit: 'm/s' });
      list.push({ key: 'steps', label: 'Steps', value: steps, kind: 'value', text: steps.toLocaleString() });
      return list;
    },
    state: () => ({ x: [...x], v: u.map(q => q / gammaOf(u)), steps }),
    record: () => ({ header: ['t_ns', 'x_m', 'y_m', 'z_m', 'vx_m_s', 'vy_m_s', 'vz_m_s'], rows }),
  };
}

// ------------------------------------------------------------------ verification

function run(p, gyrations) {
  const sim = createSim(p);
  const n = Math.round(gyrations * p.stepsPerGyration);
  let rMax = 0, rMin = Infinity;
  for (let i = 0; i < n; i++) {
    sim.step();
    const r = Math.hypot(sim.x[0], sim.x[1]);
    rMax = Math.max(rMax, r); rMin = Math.min(rMin, r);
  }
  return { sim, rMax, rMin };
}

/** Phase error after one analytic period, as an angle (E = 0). */
function phaseError(p) {
  const { sim } = run(p, 1);
  const ang = Math.atan2(sim.x[1], sim.x[0]);
  return Math.abs(Math.atan2(Math.sin(ang + Math.PI / 2), Math.cos(ang + Math.PI / 2)));
}

export function verify() {
  const checks = [];
  const base = { B: 1, KE: 1, pitch: 90, E: 0, stepsPerGyration: 256 };

  // 1. Energy: a pure magnetic field does no work; Boris conserves |u| to round-off.
  const e1 = run(base, 200);
  const dg = Math.abs(e1.sim.gamma() - e1.sim.analytic.gamma) / e1.sim.analytic.gamma;
  checks.push({ label: 'Energy conserved in pure B (200 gyrations)', value: `|Δγ/γ| = ${dg.toExponential(2)}`, pass: dg < 1e-13 });

  // 2. Orbit radius matches Boris' known discrete radius r·√(1+(ΩΔt/2)²), and → r as Δt → 0.
  const a = e1.sim.analytic, w = a.omega * e1.sim.dt;
  const rPred = a.radius * Math.sqrt(1 + (w / 2) ** 2);
  const rMeas = (e1.rMax + e1.rMin) / 2;
  const relR = Math.abs(rMeas - rPred) / rPred;
  checks.push({ label: 'Orbit radius vs exact discrete-Boris radius', value: `rel. diff ${relR.toExponential(2)}`, pass: relR < 1e-9 });
  const rErr = Math.abs(rMeas - a.radius) / a.radius;
  checks.push({ label: 'Orbit radius vs analytic γmv⊥/(qB), 256 steps/gyration', value: `rel. err ${rErr.toExponential(2)}`, pass: rErr < 1e-4 });

  // 2b. Orbit centred on the origin (half-step start): extremes of the discrete circle.
  {
    const cs = createSim({ ...base, stepsPerGyration: 32 });
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (let i = 0; i < 32 * 100; i++) {
      cs.step();
      x0 = Math.min(x0, cs.x[0]); x1 = Math.max(x1, cs.x[0]); y0 = Math.min(y0, cs.x[1]); y1 = Math.max(y1, cs.x[1]);
    }
    const off = Math.hypot((x0 + x1) / 2, (y0 + y1) / 2) / cs.analytic.radius;
    checks.push({ label: 'Orbit centred on the guiding centre (32 steps/gyration)', value: `offset ${off.toExponential(1)} r`, pass: off < 1e-4 });
  }

  // 3. Convergence: phase error after one period scales as Δt² (2nd-order method).
  const pe = [128, 256, 512].map(s => phaseError({ ...base, stepsPerGyration: s }));
  const order = Math.log2(pe[1] / pe[2]);
  checks.push({ label: 'Phase-error convergence order (expected 2)', value: order.toFixed(2), pass: Math.abs(order - 2) < 0.1 });

  // 4. Relativistic gyroperiod at 1 GeV: measure the period from orbit crossings.
  const rel = { ...base, KE: 1000, stepsPerGyration: 2048 };
  const sim = createSim(rel);
  const crossings = [];
  let prevX = sim.x[0], tPrev = 0;
  while (crossings.length < 6) {
    tPrev = sim.t; prevX = sim.x[0];
    sim.step();
    if (prevX < 0 && sim.x[0] >= 0) crossings.push(tPrev + (sim.t - tPrev) * (-prevX) / (sim.x[0] - prevX));
  }
  const Tmeas = (crossings[5] - crossings[0]) / 5;
  const Tex = sim.analytic.period * 1e9;
  const relT = Math.abs(Tmeas - Tex) / Tex;
  checks.push({ label: `Relativistic gyroperiod at 1 GeV (γ = ${sim.analytic.gamma.toFixed(4)})`, value: `${Tmeas.toFixed(4)} ns vs ${Tex.toFixed(4)} ns`, pass: relT < 1e-5 });

  // 5. E×B drift: guiding-centre velocity averaged over whole gyrations equals E/B.
  const drift = { B: 1, KE: 0.01, pitch: 90, E: 1e4, stepsPerGyration: 512 };
  // Guiding centre R = x + (γm/(qB)) (w × B̂), with w = v − v_E the velocity in the drift
  // frame: removes the gyration exactly, so ΔR/Δt is the drift velocity.
  const ds = createSim(drift);
  const vE = drift.E / drift.B;
  const gc = () => {
    const g = ds.gamma(), vx = ds.u[0] / g, vy = ds.u[1] / g;
    const k = g * mp / (qe * drift.B);
    // (w × ẑ) = (w_y, −w_x, 0), with w = (v_x − v_E, v_y)
    return [ds.x[0] + k * vy, ds.x[1] - k * (vx - vE)];
  };
  const R0 = gc(), t0 = ds.t;
  for (let i = 0; i < 50 * drift.stepsPerGyration; i++) ds.step();
  const R1 = gc();
  const vd = (R1[0] - R0[0]) / ((ds.t - t0) * 1e-9);
  const relD = Math.abs(vd - drift.E / drift.B) / (drift.E / drift.B);
  checks.push({ label: 'E×B drift speed vs E/B (non-relativistic)', value: `${vd.toFixed(1)} vs ${(drift.E / drift.B).toFixed(1)} m/s`, pass: relD < 1e-3 });

  return {
    summary: 'Relativistic Boris pusher checked against exact single-particle results:',
    checks,
    notes: ['Radiation reaction (synchrotron loss) is neglected; for protons below TeV energies in ≤10 T it is many orders of magnitude below the integration error.'],
  };
}
