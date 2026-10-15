// Physical constants, SI units. CODATA 2018 recommended values
// (Tiesinga et al., Rev. Mod. Phys. 93, 025010, 2021). "exact" = defined by the 2019 SI.
// Import from here instead of typing constants by hand; cite the source in the page.
// PDG (end of file): particle-physics inputs as { value, error, unit, source }, from the PDG 2024
// Review of Particle Physics or CODATA 2022 as each entry's source says.

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
export const mAlphaMeV = 3727.3794066;      // MeV (alpha-particle rest energy)
export const muB = 1.602176634e-19 * (6.62607015e-34 / (2 * Math.PI)) / (2 * 9.1093837015e-31); // J T^-1, Bohr magneton eħ/2mₑ (CODATA 9.2740100783e-24)
export const ge = 2.00231930436256;         // |g| of the free electron (g_e = −2.00231930436256)

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

// Everyday mechanics and the Earth–Moon system
export const gn = 9.80665;                  // m s^-2, standard gravity (exact by definition, 3rd CGPM 1901)
export const GMearth = 3.986004e14;         // m^3 s^-2 (IAU 2015 nominal terrestrial mass parameter)
export const Rearth = 6.3781e6;             // m, equatorial (IAU 2015 nominal); mean radius 6.3710e6 m
export const RearthMean = 6.3710e6;         // m (IUGG mean radius)
export const GMmoon = 4.9028e12;            // m^3 s^-2 (JPL DE430: 4902.800066 km^3 s^-2)
export const Rmoon = 1.7374e6;              // m, mean radius (NASA lunar fact sheet)
export const earthMoonDistance = 3.844e8;   // m, mean (NASA lunar fact sheet)

// Gases and the atmosphere
export const R = 8.314462618;               // J mol^-1 K^-1, molar gas constant = kB·NA (exact)
export const Mair = 28.9644e-3;             // kg mol^-1, dry air (U.S. Standard Atmosphere 1976 / ISA)
export const gammaAir = 1.4;                // cp/cv of dry air near room temperature (measured 1.401 at 20 °C)
export const p0 = 101325;                   // Pa, standard atmosphere (exact)
export const T0isa = 288.15;                // K, ISA sea-level temperature
export const rho0isa = 1.2250;              // kg m^-3, ISA sea-level density (= p0·Mair/(R·T0isa))

// Particle-physics inputs: { value, error, unit, source }. PDG 2024 = Review of Particle Physics,
// Navas et al., Phys. Rev. D 110, 030001 (2024); CODATA 2022 where it is the primary source.
// Quark masses are MS-bar running masses except the top (direct-measurement mass).
const pdg = (value, error, unit, source, note) => Object.freeze(note ? { value, error, unit, source, note } : { value, error, unit, source });
const RPP = 'PDG 2024 (RPP)';
export const PDG = Object.freeze({
  me: pdg(0.51099895069, 0.00000000016, 'MeV', 'CODATA 2022'),
  mmu: pdg(105.6583755, 0.0000023, 'MeV', 'CODATA 2022'),
  tauMu: pdg(2.1969811e-6, 0.0000022e-6, 's', RPP),
  amu: pdg(116592059e-11, 22e-11, '', RPP, 'experimental world average (BNL + FNAL 2023), a_μ = (g−2)/2'),
  mp: pdg(938.27208943, 0.00000029, 'MeV', 'CODATA 2022'),
  mpiCharged: pdg(139.57039, 0.00018, 'MeV', RPP),
  mpi0: pdg(134.9768, 0.0005, 'MeV', RPP),
  tauPiCharged: pdg(2.6033e-8, 0.0005e-8, 's', RPP, 'π± mean life (cτ = 7.80 m)'),
  tauPi0: pdg(8.43e-17, 0.13e-17, 's', RPP, 'π⁰ mean life'),
  mKCharged: pdg(493.677, 0.015, 'MeV', RPP),
  tauKCharged: pdg(1.2380e-8, 0.0020e-8, 's', RPP, 'K± mean life'),
  GF: pdg(1.1663788e-5, 0.0000006e-5, 'GeV^-2', RPP, 'G_F/(ħc)³ from the muon lifetime'),
  mW: pdg(80.3692, 0.0133, 'GeV', RPP),
  mZ: pdg(91.1876, 0.0021, 'GeV', RPP),
  mH: pdg(125.20, 0.11, 'GeV', RPP),
  sin2thetaW: pdg(0.23129, 0.00004, '', RPP, 'MS-bar, at μ = M_Z'),
  alphaS_MZ: pdg(0.1180, 0.0009, '', RPP, 'MS-bar, nf = 5, world average'),
  alphaInv0: pdg(137.035999177, 0.000000021, '', 'CODATA 2022', 'α(0)⁻¹, Thomson limit'),
  hbarc: pdg(197.3269804, 0, 'MeV fm', 'CODATA 2022', 'exact in the 2019 SI (shown truncated)'),
  mc_mc: pdg(1.2730, 0.0046, 'GeV', RPP, 'MS-bar m_c(m_c)'),
  mb_mb: pdg(4.183, 0.007, 'GeV', RPP, 'MS-bar m_b(m_b)'),
  mt: pdg(172.57, 0.29, 'GeV', RPP, 'direct measurements (≈ pole mass)'),
  michelRho: pdg(0.74979, 0.00026, '', RPP, 'Michel parameter ρ (SM: 3/4)'),
  michelDelta: pdg(0.75047, 0.00034, '', RPP, 'Michel parameter δ (SM: 3/4)'),
  michelXi: pdg(1.0009, 0.0016, '', RPP, 'Michel parameter ξ (SM: 1); asymmetric error +0.0016 −0.0007, larger side given'),
  michelEta: pdg(0.057, 0.034, '', RPP, 'Michel parameter η (SM: 0)'),
});

// Standard atomic weights (IUPAC; abridged values where the natural variation is large), g/mol.
// Multiply by u·1000 (or divide by NA·1000) for the mass of one atom in kg.
export const ATOMIC_WEIGHT = Object.freeze({
  H: 1.008, He: 4.002602, C: 12.011, N: 14.007, O: 15.999, Na: 22.98976928, Al: 26.9815384, Si: 28.085,
  Cl: 35.45, Ar: 39.95, K: 39.0983, Fe: 55.845, Cu: 63.546, Rb: 85.4678, Ag: 107.8682, Cs: 132.90545196,
  W: 183.84, Au: 196.966570, Hg: 200.592, Pb: 207.2, U: 238.02891,
});
