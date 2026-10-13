// Disk galaxy with and without a dark-matter halo: thousands of stars integrated in a
// Milky-Way mass model (core/galaxy.mjs), plus the rotation curve that reveals dark matter.

import {
  createDisk, stepDisk, energies, angularMomenta, circularSpeed, potential, acceleration,
  haloMass, haloDensity, eilers2019, MW, G, R0, V0, TIME_UNIT_MYR, NFW_SCALE,
} from 'psim:core/galaxy.mjs';
import { rng, gaussLegendre } from 'psim:core/numerics.mjs';

export const DT = 2.5e-4;                          // time units (0.244 Myr): ~900 steps per solar orbit

export const PARAMS = [
  { key: 'halo', label: 'Dark-matter halo', value: 'on', options: [
    { value: 'on', label: 'On: bulge + disk + dark halo' },
    { value: 'off', label: 'Off: only visible matter (bulge + disk)' },
  ] },
  { key: 'velocities', label: 'Initial star velocities', value: 'observed', options: [
    { value: 'observed', label: 'Full-model circular speeds' },
    { value: 'equilibrium', label: 'Circular for the chosen potential' },
  ], hint: 'With the halo off and full-model starting velocities, the disk expands and some outer stars escape; most stars remain bound.' },
  { key: 'n', label: 'Number of stars', min: 500, max: 8000, value: 3000, log: true, digits: 2, integer: true },
  { key: 'seed', label: 'Random seed', min: 1, max: 99, value: 7, digits: 2, integer: true },
];

const comps = p => ({ bulge: true, disk: true, halo: p.halo === 'on' });

export function createSim(p) {
  const n = Math.round(p.n);
  const c = comps(p);
  const disk = createDisk({ n, comps: c, velocityModel: p.velocities, rand: rng(Math.round(p.seed)) });
  const sim = {
    params: p, disk, comps: c,
    get t() { return disk.t * TIME_UNIT_MYR; },
    step: () => stepDisk(disk, DT),
    /** Mean tangential speed of the simulated stars in radial bins (the "measured" curve). */
    measuredCurve(bins = 24, Rmax = 24) {
      const sum = new Float64Array(bins), cnt = new Float64Array(bins);
      for (let i = 0; i < disk.n; i++) {
        const x = disk.x[3 * i], y = disk.x[3 * i + 1], R = Math.hypot(x, y);
        const k = Math.floor(R / Rmax * bins);
        if (k < 0 || k >= bins) continue;
        sum[k] += (x * disk.v[3 * i + 1] - y * disk.v[3 * i]) / R; cnt[k]++;
      }
      const pts = [];
      for (let k = 1; k < bins; k++) if (cnt[k] >= 8) pts.push([(k + 0.5) * Rmax / bins, sum[k] / cnt[k]]);
      return pts;
    },
    escaped() {
      let e = 0;
      for (let i = 0; i < disk.n; i++) if (Math.hypot(disk.x[3 * i], disk.x[3 * i + 1], disk.x[3 * i + 2]) > 40) e++;
      return e / disk.n;
    },
    diagnostics() {
      const E = energies(disk), L = angularMomenta(disk);
      let dE = 0, dL = 0;
      for (let i = 0; i < disk.n; i++) {
        // normalise by max(|E₀|, ½V₀²): unbound stars have E₀ ≈ 0, which would inflate a relative error
        dE = Math.max(dE, Math.abs(E[i] - disk.E0[i]) / Math.max(Math.abs(disk.E0[i]), 0.5 * V0 * V0));
        dL = Math.max(dL, Math.abs(L[i] - disk.L0[i]) / Math.abs(disk.L0[i]));
      }
      return [
        { key: 'dE', label: 'Worst star |ΔE| / max(|E₀|, ½V₀²)', value: dE, kind: 'error', good: 1e-4, warn: 1e-2 },
        { key: 'dL', label: 'Worst star |ΔL_z/L_z| (axisymmetry)', value: dL, kind: 'error', good: 1e-10, warn: 1e-6 },
        { key: 'esc', label: 'Stars beyond 40 kpc', value: 0, kind: 'value', text: `${(100 * sim.escaped()).toFixed(1)} %` },
        { key: 'unbound', label: 'Unbound stars (specific energy ≥ 0)', value: 100 * Array.from(E).filter(q => q >= 0).length / disk.n, kind: 'value', unit: '%' },
        { key: 'vc', label: 'Model v_c at R₀ = 8 kpc', value: circularSpeed(R0, c).total, kind: 'value', unit: 'km/s' },
        { key: 'orb', label: 'Orbital period at R₀', value: 2 * Math.PI * R0 / circularSpeed(R0, c).total * TIME_UNIT_MYR, kind: 'value', unit: 'Myr' },
      ];
    },
    state: () => {
      let r2 = 0;
      for (let i = 0; i < disk.n; i++) r2 += disk.x[3 * i] ** 2 + disk.x[3 * i + 1] ** 2 + disk.x[3 * i + 2] ** 2;
      return { t: disk.t, n: disk.n, rmsRadius: Math.sqrt(r2 / disk.n), escaped: sim.escaped(), first: Array.from(disk.x.subarray(0, 12)) };
    },
    record() {
      const rows = [];
      for (let R = 0.5; R <= 30; R += 0.5) {
        const v = circularSpeed(R), nd = circularSpeed(R, { bulge: true, disk: true, halo: false });
        rows.push([R, v.total, v.bulge, v.disk, v.halo, nd.total, eilers2019(R)]);
      }
      return { header: ['R_kpc', 'vc_total', 'vc_bulge', 'vc_disk', 'vc_halo', 'vc_no_dark_matter', 'vc_measured_Eilers2019'], rows };
    },
  };
  return sim;
}

export function verify() {
  const checks = [];
  const v8 = circularSpeed(R0);
  checks.push({ label: 'Normalisation: v_c(8 kpc) = 220 km/s, force fractions 5/60/35 % (MWPotential2014)', value: `${v8.total.toFixed(6)} km/s`, pass: Math.abs(v8.total - V0) < 1e-9 });
  checks.push({ label: 'Disk mass vs Bovy (2015) 6.8×10¹⁰ M☉', value: `${(MW.Md / 1e10).toFixed(2)}×10¹⁰ M☉`, pass: Math.abs(MW.Md / 6.8e10 - 1) < 0.02 });
  // analytic forces vs numerical gradient of the potential
  const a = new Float64Array(3), h = 1e-5, all = { bulge: true, disk: true, halo: true };
  let worst = 0;
  for (const [x, y, z] of [[5, 3, 0.4], [0.7, -1.2, 0.1], [15, 2, -2]]) {
    acceleration(x, y, z, all, MW, a);
    const P = (X, Y, Z) => potential(Math.hypot(X, Y), Z, all);
    const g = [-(P(x + h, y, z) - P(x - h, y, z)) / (2 * h), -(P(x, y + h, z) - P(x, y - h, z)) / (2 * h), -(P(x, y, z + h) - P(x, y, z - h)) / (2 * h)];
    for (let k = 0; k < 3; k++) worst = Math.max(worst, Math.abs(a[k] - g[k]) / Math.hypot(...g));
  }
  checks.push({ label: 'Analytic forces = −∇Φ (bulge + Miyamoto–Nagai + NFW)', value: `max rel. diff ${worst.toExponential(1)}`, pass: worst < 1e-7 });
  // NFW enclosed mass vs integral of the density
  const r = 30, Mnum = gaussLegendre(s => 4 * Math.PI * s * s * haloDensity(s), 1e-9, r, 64);
  checks.push({ label: 'NFW M(<30 kpc): formula vs ∫4πr²ρ dr', value: `rel. diff ${Math.abs(Mnum / haloMass(r) - 1).toExponential(1)}`, pass: Math.abs(Mnum / haloMass(r) - 1) < 1e-6 });
  // against the measured Milky Way rotation curve
  let dev = 0;
  for (let R = 5; R <= 25; R += 1) dev = Math.max(dev, Math.abs(circularSpeed(R).total / eilers2019(R) - 1));
  checks.push({ label: 'Model vs measured rotation curve 5–25 kpc (Eilers et al. 2019)', value: `max deviation ${(100 * dev).toFixed(1)} %`, pass: dev < 0.1, state: dev > 0.05 ? 'warn' : undefined });
  const noDM = circularSpeed(20, { bulge: true, disk: true, halo: false }).total;
  checks.push({ label: 'Without dark matter, v_c(20 kpc) falls far below the measurement', value: `${noDM.toFixed(0)} vs ${eilers2019(20).toFixed(0)} km/s measured`, pass: noDM < 0.7 * eilers2019(20) });
  // leapfrog: Lz exact, energy 2nd order
  const errAt = dt => {
    const d = createDisk({ n: 50, comps: all, rand: rng(3) });
    const steps = Math.round(0.25 / dt);
    for (let i = 0; i < steps; i++) stepDisk(d, dt);
    const E = energies(d); let m = 0;
    for (let i = 0; i < d.n; i++) m = Math.max(m, Math.abs((E[i] - d.E0[i]) / d.E0[i]));
    return { m, d };
  };
  const e1 = errAt(1e-3), e2 = errAt(5e-4);
  checks.push({ label: 'Leapfrog energy error order (expected 2)', value: Math.log2(e1.m / e2.m).toFixed(2), pass: Math.abs(Math.log2(e1.m / e2.m) - 2) < 0.4 });
  const L = angularMomenta(e2.d); let dl = 0;
  for (let i = 0; i < e2.d.n; i++) dl = Math.max(dl, Math.abs(L[i] - e2.d.L0[i]) / Math.abs(e2.d.L0[i]));
  checks.push({ label: 'L_z conserved (axisymmetric potential, 245 Myr)', value: dl.toExponential(1), pass: dl < 1e-12 });
  return {
    summary: 'Mass model, forces and integrator checked against exact identities, Bovy (2015) and the measured Milky Way:',
    checks,
    notes: ['The model is normalised to 220 km/s at 8 kpc (MWPotential2014); modern measurements give ≈229 km/s, so it runs 4–8 % low.'],
  };
}

export { circularSpeed, eilers2019, haloDensity, R0, G, NFW_SCALE, TIME_UNIT_MYR, DT as STEP };
