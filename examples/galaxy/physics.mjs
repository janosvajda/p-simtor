// Disk galaxy with and without a dark-matter halo: two copies of the same star disk integrated
// side by side in a Milky-Way mass model (core/galaxy.mjs), one with the halo and one with the
// visible matter only, plus their live rotation curves (the evidence for dark matter).
// Stars are test particles (no star–star gravity), so this model cannot form spiral arms.

import {
  createDisk, stepDisk, energies, angularMomenta, circularSpeed, potential, acceleration,
  haloMass, haloDensity, eilers2019, MW, G, R0, V0, TIME_UNIT_MYR, NFW_SCALE,
} from 'psim:core/galaxy.mjs';
import { rng, gaussLegendre } from 'psim:core/numerics.mjs';

export const DT = 2.5e-4;                          // time units (0.244 Myr): ~900 steps per solar orbit

export const PARAMS = [
  { key: 'n', label: 'Stars in each galaxy', min: 500, max: 8000, value: 3000, log: true, digits: 2, integer: true },
  { key: 'twin', label: 'Starting speeds of the galaxy without the halo', value: 'observed', options: [
    { value: 'observed', label: 'Measured speeds (same as the left galaxy)' },
    { value: 'equilibrium', label: 'Slow enough for visible matter to hold' },
  ], hint: 'With measured speeds and no halo, the disk expands and some outer stars escape. The slow option stays together but orbits far below the measured speeds.' },
  { key: 'seed', label: 'Random seed', min: 1, max: 99, value: 7, digits: 2, integer: true },
];

const comps = p => ({ bulge: true, disk: true, halo: p.halo === 'on' });

/** One galaxy, params { halo: 'on' | 'off', velocities: 'observed' | 'equilibrium', n, seed }. */
export function createRun(p) {
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

/**
 * The page's model: two copies of the same disk (same seed, so the same stars), stepped together.
 * `halo` feels bulge + disk + dark-matter halo with the full-model (measured-like) speeds; `none`
 * feels the visible matter only and starts with params.twin speeds.
 */
export function createSim(p) {
  const n = Math.round(p.n), seed = Math.round(p.seed);
  const A = createRun({ halo: 'on', velocities: 'observed', n, seed });
  const B = createRun({ halo: 'off', velocities: p.twin ?? 'observed', n, seed });
  return {
    params: p, runs: [A, B], halo: A, none: B,
    get t() { return A.t; },
    step() { A.step(); B.step(); },
    diagnostics() {
      const a = A.diagnostics(), b = B.diagnostics(), get = (d, k) => d.find(q => q.key === k);
      return [
        { key: 'dE', label: 'Worst star |ΔE| / max(|E₀|, ½V₀²), both galaxies', value: Math.max(get(a, 'dE').value, get(b, 'dE').value), kind: 'error', good: 1e-4, warn: 1e-2 },
        { key: 'dL', label: 'Worst star |ΔL_z/L_z|, both galaxies', value: Math.max(get(a, 'dL').value, get(b, 'dL').value), kind: 'error', good: 1e-10, warn: 1e-6 },
        { key: 'escA', label: 'With halo: stars beyond 40 kpc', value: 0, kind: 'value', text: get(a, 'esc').text },
        { key: 'escB', label: 'Without halo: stars beyond 40 kpc', value: 0, kind: 'value', text: get(b, 'esc').text },
        { key: 'unbB', label: 'Without halo: unbound stars', value: get(b, 'unbound').value, kind: 'value', unit: '%' },
        { key: 'vc', label: 'v_c at R₀ = 8 kpc with / without halo', value: 0, kind: 'value', text: `${get(a, 'vc').value.toFixed(1)} / ${get(b, 'vc').value.toFixed(1)} km/s` },
        { key: 'orb', label: 'Orbital period at R₀ (with halo)', value: get(a, 'orb').value, kind: 'value', unit: 'Myr' },
      ];
    },
    state: () => ({ t: A.disk.t, halo: A.state(), none: B.state() }),
    record: () => A.record(),
  };
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
  // the page's evidence: the stars' own rotation curves with and without the halo, after 500 Myr
  const tw = createSim({ n: 1500, twin: 'observed', seed: 5 });
  const steps500 = Math.round(500 / (DT * TIME_UNIT_MYR));
  for (let k = 0; k < steps500; k++) tw.step();
  const band = (run, ref) => { const c = run.measuredCurve(24, 24).filter(([R]) => R >= 8 && R <= 16); return Math.max(...c.map(([R, v]) => Math.abs(v / ref(R) - 1))); };
  const devHalo = band(tw.halo, R => circularSpeed(R).total);
  const lowNone = Math.max(...tw.none.measuredCurve(24, 24).filter(([R]) => R >= 12 && R <= 20).map(([R, v]) => v / eilers2019(R)));
  checks.push({ label: 'With halo after 500 Myr: stars’ mean v_φ at 8–16 kpc vs circular speed (1500 stars)', value: `max deviation ${(100 * devHalo).toFixed(1)} %`, pass: devHalo < 0.06 });
  checks.push({ label: 'Without halo after 500 Myr: stars’ mean v_φ at 12–20 kpc vs measured (Eilers et al. 2019)', value: `at most ${(100 * lowNone).toFixed(0)} % of the measured speed`, pass: lowNone < 0.8 });
  return {
    summary: 'Mass model, forces and integrator checked against exact identities, Bovy (2015) and the measured Milky Way:',
    checks,
    notes: ['The model is normalised to 220 km/s at 8 kpc (MWPotential2014); modern measurements give ≈229 km/s, so it runs 4–8 % low.'],
  };
}

export { circularSpeed, eilers2019, haloDensity, R0, G, NFW_SCALE, TIME_UNIT_MYR, DT as STEP };
