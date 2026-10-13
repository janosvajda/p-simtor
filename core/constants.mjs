// Physical constants, SI units. CODATA 2018 recommended values
// (Tiesinga et al., Rev. Mod. Phys. 93, 025010, 2021). "exact" = defined by the 2019 SI.
// Import from here instead of typing constants by hand; cite the source in the page.

export const CODATA = 'CODATA 2018';

export const c = 299792458;                 // m s^-1, exact
export const h = 6.62607015e-34;            // J s, exact
export const hbar = 6.62607015e-34 / (2 * Math.PI);  // J s, exact (h/2π)
export const e = 1.602176634e-19;           // C, exact
export const kB = 1.380649e-23;             // J K^-1, exact
export const NA = 6.02214076e23;            // mol^-1, exact
export const eps0 = 8.8541878128e-12;       // F m^-1
export const mu0 = 1.25663706212e-6;        // N A^-2
export const G = 6.67430e-11;               // m^3 kg^-1 s^-2 (relative uncertainty 2.2e-5)
export const me = 9.1093837015e-31;         // kg
export const mp = 1.67262192369e-27;        // kg
export const mn = 1.67492749804e-27;        // kg
export const u = 1.66053906660e-27;         // kg (atomic mass constant)
export const a0 = 5.29177210903e-11;        // m (Bohr radius)
export const alpha = 7.2973525693e-3;       // fine-structure constant
export const Eh = 4.3597447222071e-18;      // J (Hartree energy)
export const Ry = 13.605693122994;          // eV (Rydberg energy)
export const eV = 1.602176634e-19;          // J, exact
export const MeV = 1.602176634e-13;         // J, exact
export const mpc2MeV = 938.27208816;        // MeV (proton rest energy)
export const mec2MeV = 0.51099895000;       // MeV (electron rest energy)
export const hbarcMeVfm = 197.3269804;      // MeV fm

// Astronomy (IAU 2012/2015 nominal values)
export const AU = 1.495978707e11;           // m, exact (IAU 2012)
export const GMsun = 1.3271244e20;          // m^3 s^-2 (IAU 2015 nominal)
export const Msun = 1.3271244e20 / 6.67430e-11;     // kg, GM☉/G (G-limited, ±2.2e-5)
export const yearJulian = 365.25 * 86400;   // s
/**
 * G in AU³ M☉⁻¹ yr⁻² with the Julian year and IAU GM☉: 39.47693 (≈ 4π², which is exact only
 * for the Gaussian year 2π/k). Using this value keeps AU–yr runs consistent with SI ones.
 */
export const G_AU_MSUN_YR = 1.3271244e20 * (365.25 * 86400) ** 2 / 1.495978707e11 ** 3;
