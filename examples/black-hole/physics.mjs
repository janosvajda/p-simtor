// Schwarzschild black hole seen by a distant static observer: exact null geodesics
// (ray traced on the GPU with the same RK4 scheme as core/schwarzschild.mjs), a thin
// accretion disk emitting black-body light, and exact gravitational + Doppler redshift.
// Geometric units G = c = 1, lengths in M.

import {
  tracePhoton, deflectionExact, deflectionWeak, B_CRIT, R_ISCO, R_PHOTON, R_S, shadowAngle,
  orbitalEnergy, orbitalL, effectivePotential, diskRedshift, diskTemperature, impactFromAngle,
} from 'psim:core/schwarzschild.mjs';
import { blackbodyColor, WIEN_B } from 'psim:core/radiation.mjs';
import { G as Gsi, c as csi, Msun } from 'psim:core/constants.mjs';

export const PARAMS = [
  { key: 'ro', label: 'Observer distance r₀', min: 8, max: 400, value: 60, log: true, unit: 'M', digits: 3 },
  { key: 'incl', label: 'Inclination (0° = face-on)', min: 0, max: 89, value: 82, unit: '°', digits: 3 },
  { key: 'rout', label: 'Disk outer radius', min: 8, max: 60, value: 20, unit: 'M', digits: 3 },
  { key: 'Tmax', label: 'Peak disk temperature', min: 2000, max: 1e7, value: 5200, log: true, unit: 'K', digits: 3,
    hint: 'Stellar-mass black-hole disks reach ~10⁷ K (X-rays) and look blue-white; lower values show the colour gradient.' },
  { key: 'effects', label: 'Light shift', value: 'full', options: [
    { value: 'full', label: 'Full: gravitational + Doppler shift and beaming' },
    { value: 'grav', label: 'Gravitational redshift only (no disk rotation)' },
    { value: 'none', label: 'None (geometry only, for comparison)' },
  ] },
  { key: 'massSun', label: 'Black-hole mass (for physical units)', min: 1, max: 1e10, value: 10, log: true, unit: 'M☉', digits: 3 },
  { key: 'rings', label: 'Radius rings on the disk (every 2M)', value: 'on', options: [{ value: 'on', label: 'On: show disk geometry' }, { value: 'off', label: 'Off' }] },
  { key: 'fov', label: 'Field of view', min: 10, max: 90, value: 30, unit: '°', digits: 3 },
  { key: 'exposure', label: 'Display exposure (not physics)', min: 0.05, max: 20, value: 2.5, log: true, digits: 2 },
];

/** Blackbody colour lookup table for the shader: log-spaced T, sRGB-encoded bytes (brightest channel = 255). */
export function colorLUT(n = 256, T0 = 500, T1 = 2e5) {
  const data = new Uint8Array(n * 4);
  for (let i = 0; i < n; i++) {
    const T = T0 * Math.pow(T1 / T0, i / (n - 1));
    const { rgb } = blackbodyColor(T);
    data.set([rgb[0], rgb[1], rgb[2], 255], i * 4);
  }
  return { data, n, T0, T1 };
}

export function createSim(p) {
  const shadow = shadowAngle(p.ro);
  const metres = Gsi * p.massSun * Msun / (csi * csi);         // GM/c² in metres
  return {
    t: 0, params: p, timeDependent: false,
    shadow, metres,
    diagnostics() {
      return [
        { key: 'sh', label: 'Shadow angular radius (exact)', value: shadow * 180 / Math.PI, kind: 'value', unit: '°' },
        { key: 'bc', label: 'Shadow edge impact parameter 3√3 M', value: B_CRIT, kind: 'value', unit: 'M' },
        { key: 'rs', label: 'Event horizon r_s = 2GM/c²', value: 2 * metres / 1000, kind: 'value', unit: 'km' },
        { key: 'rp', label: 'Photon sphere 3M', value: 3 * metres / 1000, kind: 'value', unit: 'km' },
        { key: 'ri', label: 'Inner disk edge (ISCO) 6M', value: 6 * metres / 1000, kind: 'value', unit: 'km' },
        { key: 'eff', label: 'Radiative efficiency 1 − E_isco', value: 1 - orbitalEnergy(6), kind: 'value', text: `${(100 * (1 - orbitalEnergy(6))).toFixed(2)} %` },
        { key: 'lp', label: 'Disk peak emission λ_max = b/T', value: WIEN_B / p.Tmax * 1e9, kind: 'value', unit: 'nm' },
      ];
    },
    state: () => ({ ro: p.ro, incl: p.incl, shadowDeg: shadow * 180 / Math.PI }),
  };
}

export function verify() {
  const checks = [];
  for (const b of [6, 10, 30]) {
    const a = tracePhoton(b, { dphi: 1e-3 }).deflection, e = deflectionExact(b);
    checks.push({ label: `Light bending at b = ${b}M: RK4 tracer vs exact quadrature`, value: `${a.toFixed(9)} vs ${e.toFixed(9)} rad`, pass: Math.abs(a - e) < 1e-8 });
  }
  const bw = 1000, aw = tracePhoton(bw).deflection;
  checks.push({ label: 'Weak field b = 1000M vs 4M/b + (15π/4)(M/b)² + (128/3)(M/b)³', value: `rel. diff ${(Math.abs(aw - deflectionWeak(bw)) / aw).toExponential(1)}`, pass: Math.abs(aw - deflectionWeak(bw)) / aw < 1e-6 });
  let lo = 5, hi = 5.5;
  for (let i = 0; i < 48; i++) { const m = (lo + hi) / 2; tracePhoton(m, { dphi: 2e-4, maxPhi: 80 }).captured ? (lo = m) : (hi = m); }
  checks.push({ label: 'Capture threshold vs 3√3 M (shadow edge)', value: `${lo.toFixed(10)} vs ${B_CRIT.toFixed(10)}`, pass: Math.abs(lo - B_CRIT) < 1e-8 });
  // Sun deflection: 1.75 arcsec (Eddington 1919), from the exact deflection integral in solar units
  const bSun = 6.957e8 / (Gsi * Msun / (csi * csi));
  const arc = deflectionExact(bSun, 256) * 206264.806;
  checks.push({ label: 'Starlight grazing the Sun (4GM/bc²)', value: `${arc.toFixed(3)}″ (Eddington 1919: 1.75″)`, pass: Math.abs(arc - 1.751) < 0.01 });
  const L6 = orbitalL(6), h = 1e-3;
  const d2 = (effectivePotential(6 + h, L6) - 2 * effectivePotential(6, L6) + effectivePotential(6 - h, L6)) / (h * h);
  checks.push({ label: 'ISCO at 6M (marginal stability, V″ = 0)', value: `V″ = ${d2.toExponential(1)}`, pass: Math.abs(d2) < 1e-6 });
  checks.push({ label: 'Binding efficiency 1 − √(8/9)', value: `${(100 * (1 - orbitalEnergy(6))).toFixed(3)} %`, pass: Math.abs(1 - orbitalEnergy(6) - (1 - Math.sqrt(8 / 9))) < 1e-15 });
  const g0 = diskRedshift(10, 0);
  checks.push({ label: 'Redshift of a face-on emitter at 10M: √(1 − 3M/r)', value: g0.toFixed(12), pass: Math.abs(g0 - Math.sqrt(0.7)) < 1e-15 });
  return {
    summary: 'Geodesics, orbits and redshift checked against exact general-relativistic results:',
    checks,
    notes: [
      'Disk emission uses the Shakura–Sunyaev temperature profile (Newtonian zero-torque form) as a modelling approximation; geodesics and redshifts are exact.',
      'The GPU image is checked in the browser against the exact shadow radius (psim check).',
    ],
  };
}

export { B_CRIT, R_ISCO, R_PHOTON, R_S, diskTemperature, impactFromAngle };
