// Charged particles passing through matter: mean energy loss and range (PDG Review of Particle
// Physics, "Passage of particles through matter", eq. 34.5, with the Sternheimer density effect).
// Units: energies in MeV, mass stopping power in MeV cm²/g, ranges in g/cm² (divide by ρ for cm).
//
//   MATERIALS                          ρ, Z/A, mean excitation I and Sternheimer parameters (PDG tables)
//   betheStoppingPower(βγ, { mass, z, material })   mean collision loss of a heavy particle (M ≫ mₑ)
//   densityEffect(βγ, material)        Sternheimer δ(βγ)
//   csdaRange(T, { mass, z, material, Tmin })       ∫ dT / S from Tmin to T (continuous-slowing-down range)
//
// Validity: 0.1 ≲ βγ ≲ 1000 for heavy particles (muons, pions, protons, alphas with z = 2).
// Below βγ ≈ 0.1 shell and Barkas corrections matter (the plain formula overestimates the range of
// a few-MeV muon in hydrogen by ~6 %); electrons and positrons need Møller/Bhabha and
// bremsstrahlung instead. This is the MEAN loss: real tracks fluctuate (Landau straggling).

import { adaptiveSimpson } from './numerics.mjs';

const ME = 0.51099895000;            // electron mass, MeV (CODATA 2018)
const K = 0.307075;                  // 4π N_A r_e² mₑc², MeV mol⁻¹ cm² (PDG)

/**
 * Material data from the PDG "Atomic and nuclear properties" tables: ρ (g/cm³), ZA = ⟨Z/A⟩ (mol/g),
 * I (MeV), radiation length X0 (g/cm²) and the Sternheimer density-effect parameters
 * (a, k, x0, x1, Cbar, delta0). Checked in tests against the PDG minimum ionisation.
 */
export const MATERIALS = Object.freeze({
  liquidHydrogen: Object.freeze({ name: 'liquid hydrogen', rho: 0.0708, ZA: 0.99212, I: 21.8e-6, X0: 63.04, a: 0.1348, k: 5.6249, x0: 0.4400, x1: 1.8856, Cbar: 3.0977, delta0: 0, minIon: 4.034 }),
  water: Object.freeze({ name: 'water (liquid)', rho: 1.0, ZA: 0.55509, I: 75.0e-6, X0: 36.08, a: 0.0912, k: 3.4773, x0: 0.2400, x1: 2.8004, Cbar: 3.5017, delta0: 0, minIon: 1.992 }),
  silicon: Object.freeze({ name: 'silicon', rho: 2.329, ZA: 0.49848, I: 173.0e-6, X0: 21.82, a: 0.1492, k: 3.2546, x0: 0.2015, x1: 2.8716, Cbar: 4.4355, delta0: 0.14, minIon: 1.664 }),
  lead: Object.freeze({ name: 'lead', rho: 11.35, ZA: 0.39575, I: 823.0e-6, X0: 6.37, a: 0.0936, k: 3.1608, x0: 0.3776, x1: 3.8073, Cbar: 6.2018, delta0: 0.14, minIon: 1.122 }),
  air: Object.freeze({ name: 'dry air (sea level)', rho: 1.205e-3, ZA: 0.49919, I: 85.7e-6, X0: 36.62, a: 0.1091, k: 3.3994, x0: 1.7418, x1: 4.2759, Cbar: 10.5961, delta0: 0, minIon: 1.815 }),
});

/** Sternheimer density-effect correction δ(βγ). */
export function densityEffect(bg, mat) {
  const x = Math.log10(bg), { a, k, x0, x1, Cbar, delta0 } = mat;
  if (x >= x1) return 2 * Math.LN10 * x - Cbar;
  if (x >= x0) return 2 * Math.LN10 * x - Cbar + a * (x1 - x) ** k;
  return delta0 * 10 ** (2 * (x - x0));
}

/** Mean mass stopping power −⟨dE/dx⟩/ρ (MeV cm²/g) of a particle of charge z·e and mass (MeV) at βγ. */
export function betheStoppingPower(bg, { mass, z = 1, material }) {
  if (!(bg > 0 && mass > 0)) throw new Error('betheStoppingPower: need βγ > 0 and mass > 0');
  const mat = typeof material === 'string' ? MATERIALS[material] : material;
  if (!mat) throw new Error(`betheStoppingPower: unknown material ${material}`);
  const g = Math.sqrt(1 + bg * bg), b2 = bg * bg / (g * g), r = ME / mass;
  const Wmax = 2 * ME * bg * bg / (1 + 2 * g * r + r * r);
  return K * z * z * mat.ZA / b2 * (0.5 * Math.log(2 * ME * bg * bg * Wmax / (mat.I * mat.I)) - b2 - densityEffect(bg, mat) / 2);
}

/** Continuous-slowing-down range (g/cm²) from kinetic energy T down to Tmin (MeV). */
export function csdaRange(T, { mass, z = 1, material, Tmin = 1e-3 * mass }) {
  const S = t => betheStoppingPower(Math.sqrt(t * (t + 2 * mass)) / mass, { mass, z, material });
  if (!(T > Tmin)) return 0;
  // integrate in ln T for accuracy over decades
  return adaptiveSimpson(u => { const t = Math.exp(u); return t / S(t); }, Math.log(Tmin), Math.log(T), 1e-10);
}
