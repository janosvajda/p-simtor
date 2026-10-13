// Black-body radiation: Planck's law and everything that follows from it.
// Exact formulas evaluated with CODATA constants (tier: analytical), plus numerical
// integration checks. Lengths in m internally; wavelengths shown in nm.

import {
  planckLambda, rayleighJeans, wienApprox, WIEN_B, WIEN_X, WIEN_X_NU, SIGMA_SB, exitance, bandRadiance,
  photonDensity, MEAN_PHOTON_ENERGY_KT, blackbodyColor, chromaticity, spectrumToXYZ, c2, c1L,
} from 'psim:core/radiation.mjs';
import { kB, h, c, eV } from 'psim:core/constants.mjs';

export const PRESETS = [
  { label: 'Cosmic microwave background', T: 2.7255, note: 'Fixsen 2009 (COBE/FIRAS)' },
  { label: 'Human body', T: 310 },
  { label: 'Incandescent lamp (CIE illuminant A)', T: 2856 },
  { label: 'Sun (effective temperature)', T: 5772, note: 'IAU 2015 nominal' },
  { label: 'Sirius A', T: 9940 },
  { label: 'Hot O-type star', T: 40000 },
];

export const PARAMS = [
  { key: 'T', label: 'Temperature T', min: 1, max: 100000, value: 5772, log: true, unit: 'K', digits: 4 },
  { key: 'scale', label: 'Spectrum axes', value: 'log', options: [
    { value: 'log', label: 'Log–log (all temperatures)' },
    { value: 'lin', label: 'Linear (shape near the peak)' },
  ] },
];

export const VISIBLE = [380e-9, 750e-9];

export function createSim(p) {
  const T = p.T;
  const lmax = WIEN_B / T;
  const color = blackbodyColor(T);
  const visibleFraction = bandRadiance(T, VISIBLE[0], VISIBLE[1]) / (exitance(T) / Math.PI);
  return {
    t: 0, params: p, T, lmax, color, visibleFraction,
    B: l => planckLambda(l, T),
    diagnostics() {
      return [
        { key: 'lmax', label: 'Peak wavelength λ_max = b/T', value: lmax * 1e9, kind: 'value', unit: 'nm' },
        { key: 'M', label: 'Radiant exitance σT⁴', value: exitance(T), kind: 'value', unit: 'W/m²' },
        { key: 'vis', label: 'Fraction emitted as visible light', value: visibleFraction, kind: 'value', text: `${(100 * visibleFraction).toPrecision(3)} %` },
        { key: 'n', label: 'Photon density in a cavity', value: photonDensity(T), kind: 'value', unit: 'm⁻³' },
        { key: 'Eph', label: 'Mean photon energy 2.701 kT', value: MEAN_PHOTON_ENERGY_KT * kB * T / eV, kind: 'value', unit: 'eV' },
        { key: 'xy', label: 'CIE 1931 chromaticity (x, y)', value: 0, kind: 'value', text: color.xy.map(v => v.toFixed(4)).join(', ') },
      ];
    },
    state: () => ({ T, lmax, visibleFraction, rgb: color.rgb }),
  };
}

export function verify() {
  const checks = [];
  const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
  checks.push({ label: 'Wien constant b = hc/(x k), x = 5(1 − e⁻ˣ) root, vs CODATA 2018', value: `${WIEN_B.toExponential(9)} m K`, pass: rel(WIEN_B, 2.897771955e-3) < 1e-9 });
  checks.push({ label: 'σ = 2π⁵k⁴/(15h³c²) vs CODATA 2018', value: `${SIGMA_SB.toExponential(9)} W m⁻² K⁻⁴`, pass: rel(SIGMA_SB, 5.670374419e-8) < 1e-9 });
  const T = 5772, I = bandRadiance(T, 1e-8, 1e-2);
  checks.push({ label: 'Numerical ∫B_λ dλ vs σT⁴/π (Gauss–Legendre, 10 nm–1 cm)', value: `rel. err ${rel(I, exitance(T) / Math.PI).toExponential(1)}`, pass: rel(I, exitance(T) / Math.PI) < 1e-9 });
  // Sun: luminosity from the IAU nominal T_eff and radius → solar constant at 1 AU
  const Rsun = 6.957e8, AU = 1.495978707e11, L = 4 * Math.PI * Rsun * Rsun * exitance(5772);
  const S = L / (4 * Math.PI * AU * AU);
  checks.push({ label: 'Solar constant from T_eff = 5772 K (IAU 2015) vs measured 1361 W/m² (Kopp & Lean 2011)', value: `${S.toFixed(1)} W/m²`, pass: Math.abs(S - 1361) < 2 });
  checks.push({ label: 'Solar luminosity vs IAU nominal 3.828×10²⁶ W', value: `${L.toExponential(4)} W`, pass: rel(L, 3.828e26) < 1e-3 });
  // CMB: peak of B_ν at 2.7255 K is 160.2 GHz (Fixsen 2009)
  const nuPeak = WIEN_X_NU * kB * 2.7255 / h;
  checks.push({ label: 'CMB spectral peak (B_ν) at 2.7255 K vs 160.2 GHz (COBE/FIRAS)', value: `${(nuPeak / 1e9).toFixed(1)} GHz`, pass: Math.abs(nuPeak / 1e9 - 160.2) < 0.1 });
  checks.push({ label: 'CMB photon density vs 411 cm⁻³', value: `${(photonDensity(2.7255) / 1e6).toFixed(1)} cm⁻³`, pass: Math.abs(photonDensity(2.7255) / 1e6 - 410.7) < 1 });
  // Colour: Planckian locus at 6500 K and illuminant A (limited by the CMF fit)
  const [x65, y65] = blackbodyColor(6500).xy;
  checks.push({ label: 'Chromaticity at 6500 K vs Planckian locus (0.3135, 0.3236)', value: `(${x65.toFixed(4)}, ${y65.toFixed(4)})`, pass: Math.abs(x65 - 0.3135) < 0.002 && Math.abs(y65 - 0.3236) < 0.002 });
  const [xA, yA] = chromaticity(spectrumToXYZ(nm => planckLambda(nm * 1e-9, 2856 * c2 / 1.435e-2)));
  checks.push({ label: 'CIE illuminant A vs (0.44757, 0.40745)', value: `(${xA.toFixed(4)}, ${yA.toFixed(4)}) ±0.002 CMF fit`, pass: Math.abs(xA - 0.44757) < 0.003 && Math.abs(yA - 0.40745) < 0.003, state: 'warn' });
  // Limits
  checks.push({ label: 'Rayleigh–Jeans limit at λ = 1 m, 300 K', value: `rel. diff ${rel(planckLambda(1, 300), rayleighJeans(1, 300)).toExponential(1)}`, pass: rel(planckLambda(1, 300), rayleighJeans(1, 300)) < 1e-4 });
  checks.push({ label: 'Wien limit at λ = 100 nm, 300 K', value: `rel. diff ${rel(planckLambda(1e-7, 300), wienApprox(1e-7, 300)).toExponential(1)}`, pass: rel(planckLambda(1e-7, 300), wienApprox(1e-7, 300)) < 1e-12 });
  return {
    summary: 'Planck’s law, its integrals and colour checked against CODATA and astronomical measurements:',
    checks,
    notes: ['Colour uses the Wyman–Sloan–Shirley (2013) analytic fit to the CIE 1931 2° observer (±0.002 in chromaticity).'],
  };
}

export { WIEN_B, WIEN_X, SIGMA_SB, c1L, c2, exitance, rayleighJeans, wienApprox, planckLambda, kB, h, c };
