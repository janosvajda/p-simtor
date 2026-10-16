// Homogeneous cosmology: the Friedmann–Lemaître–Robertson–Walker background.
//
// A model is { H0 (km s⁻¹ Mpc⁻¹), Om (matter), Or (radiation: photons + massless neutrinos),
// OL (cosmological constant), Ok = 1 − Om − Or − OL (curvature) } built by `cosmology()`.
// Units: times in Gyr (Julian years), distances in Mpc, H in km s⁻¹ Mpc⁻¹, temperatures in K.
//
//   E(z)² = Or(1+z)⁴ + Om(1+z)³ + Ok(1+z)² + OL         (H = H0·E)
//
// Integrals are done in s = √a, where every integrand is smooth (a ∝ t^{2/3} and a ∝ t^{1/2}
// near the Big Bang become polynomial), with Gauss–Legendre on geometric panels: relative
// accuracy ~1e-13 for the presets. Models without an expanding history back to a = 0 (bounce,
// pure de Sitter) throw from the functions that need the Big Bang.
//
// Constants: pc = 648000/π au (IAU 2015 Resolution B2), Julian year, CODATA 2018 via constants.mjs.
// Planck 2018 parameters: Planck Collaboration VI, A&A 641, A6 (2020), arXiv:1807.06209,
// Table 2, column "TT,TE,EE+lowE+lensing+BAO".

import { c, G, kB, hbar, AU, yearJulian } from './constants.mjs';
import { gaussLegendre, brent } from './numerics.mjs';
import { odeStepper } from './ode.mjs';

/** Megaparsec in m: 10⁶ × 648000/π au (IAU 2015 B2). */
export const MPC = AU * 648000 / Math.PI * 1e6;
/** Gigayear (Julian) in s. */
export const GYR = 1e9 * yearJulian;
/** Speed of light in km/s. */
export const C_KMS = c / 1e3;
/** 1/(1 km s⁻¹ Mpc⁻¹) in Gyr ≈ 977.79: Hubble time t_H = HUBBLE_TIME_UNIT / H0. */
export const HUBBLE_TIME_UNIT = MPC / 1e3 / GYR;
/** CMB temperature today, K (Fixsen, ApJ 707, 916 (2009): 2.72548 ± 0.00057). */
export const T_CMB = 2.7255;
/** Massless-neutrino energy per photon energy per effective species: (7/8)(4/11)^{4/3}. */
export const NEUTRINO_FACTOR = 7 / 8 * Math.pow(4 / 11, 4 / 3);

/**
 * Photon density parameter times h² from the CMB temperature: Ωγh² = 8πG·(π²/15)(kT)⁴/(ħc)³ / (3 (100 km/s/Mpc)² c²).
 * 2.4728e-5 at 2.7255 K (PDG Astrophysical constants: 2.473e-5).
 */
export function omegaGammaH2(Tcmb = T_CMB) {
  if (!(Tcmb >= 0 && Number.isFinite(Tcmb))) throw new Error('omegaGammaH2: Tcmb must be finite and ≥ 0');
  const u = Math.PI ** 2 / 15 * (kB * Tcmb) ** 4 / (hbar * c) ** 3;   // J m⁻³
  const H100 = 1e5 / MPC;                                            // s⁻¹
  return 8 * Math.PI * G * u / (3 * H100 * H100 * c * c);
}

/**
 * Build a model. Arguments: H0 (km/s/Mpc, > 0), Om (≥ 0), OL (finite) or flat: true (OL = 1 − Om − Or),
 * Or (≥ 0; default from Tcmb and Neff: Ωγ(1 + 0.2271 Neff)), Tcmb (K, default 2.7255), Neff (default 3.046),
 * optional Obh2 (baryons, for the sound horizon). Ok = 1 − Om − Or − OL. Returns a frozen object
 * { H0, h, Om, Or, OL, Ok, Tcmb, Neff, Obh2, Ogamma, name }.
 */
export function cosmology({ H0, Om, OL, Or, Tcmb = T_CMB, Neff = 3.046, flat = false, Obh2, name = 'custom' } = {}) {
  if (!(Number.isFinite(H0) && H0 > 0)) throw new Error('cosmology: H0 (km/s/Mpc) must be finite and > 0');
  if (!(Number.isFinite(Om) && Om >= 0)) throw new Error('cosmology: Om must be finite and ≥ 0');
  if (!(Number.isFinite(Neff) && Neff >= 0)) throw new Error('cosmology: Neff must be finite and ≥ 0');
  const h = H0 / 100;
  const Ogamma = omegaGammaH2(Tcmb) / (h * h);
  const Orad = Or ?? Ogamma * (1 + NEUTRINO_FACTOR * Neff);
  if (!(Number.isFinite(Orad) && Orad >= 0)) throw new Error('cosmology: Or must be finite and ≥ 0');
  const lam = flat ? 1 - Om - Orad : OL;
  if (!Number.isFinite(lam)) throw new Error('cosmology: give OL or flat: true');
  if (Obh2 != null && !(Number.isFinite(Obh2) && Obh2 > 0)) throw new Error('cosmology: Obh2 must be > 0');
  if (Obh2 != null && Obh2 > Om * h * h) throw new Error('cosmology: baryon density Obh2 cannot exceed the total matter density Om h²');
  const Ok = flat ? 0 : 1 - Om - Orad - lam;
  return Object.freeze({ name, H0, h, Om, Or: Orad, OL: lam, Ok: Math.abs(Ok) < 1e-15 ? 0 : Ok, Tcmb, Neff, Obh2, Ogamma });
}

/**
 * Planck 2018 (TT,TE,EE+lowE+lensing+BAO; arXiv:1807.06209 Table 2): H0 = 67.66, Ωm = 0.3111,
 * Ωbh² = 0.02242, flat. Radiation: photons at 2.7255 K + Neff = 3.046 massless neutrinos, so
 * ΩΛ = 1 − Ωm − Ωr = 0.68881 (Table 2 rounds to 0.6889). Planck's one 0.06 eV neutrino is counted in
 * Ωm here and also as radiation via Neff (≈0.2 % effect on H at z ~ 1000; disclose).
 */
export const PLANCK18 = cosmology({ H0: 67.66, Om: 0.3111, flat: true, Obh2: 0.02242, name: 'Planck 2018' });
/** Published Planck 2018 derived values (same column) with 1σ errors, for verification. */
export const PLANCK18_DERIVED = Object.freeze({
  age: { value: 13.787, error: 0.020, unit: 'Gyr' },
  zEq: { value: 3387, error: 21, unit: '' },
  zStar: { value: 1089.80, error: 0.21, unit: '' },
  rStar: { value: 144.57, error: 0.22, unit: 'Mpc' },
  zDrag: { value: 1060.01, error: 0.29, unit: '' },
  rDrag: { value: 147.21, error: 0.23, unit: 'Mpc' },
  Och2: { value: 0.11933, error: 0.00091, unit: '' },
  source: 'Planck 2018 VI, arXiv:1807.06209, Table 2 (TT,TE,EE+lowE+lensing+BAO)',
});

const checkModel = (p, who) => {
  if (!p || !Number.isFinite(p.H0) || !Number.isFinite(p.Om) || !Number.isFinite(p.Or) || !Number.isFinite(p.OL)) throw new Error(`${who}: pass a model from cosmology()`);
};
const checkZ = (z, who) => { if (!(Number.isFinite(z) && z > -1)) throw new Error(`${who}: z must be finite and > −1`); };
// Q(a) = a⁴E² = Or + Om a + Ok a² + OL a⁴
const Q = (p, a) => p.Or + a * (p.Om + a * (p.Ok + a * a * p.OL));

/** Dimensionless Hubble rate E(z) = H(z)/H0. Throws where E² ≤ 0 (no classical expanding solution). */
export function E(z, p) {
  checkModel(p, 'E'); checkZ(z, 'E');
  const x = 1 + z, e2 = x * x * (x * (x * p.Or + p.Om) + p.Ok) + p.OL;
  if (!(e2 > 0)) throw new Error(`E: E² = ${e2} ≤ 0 at z = ${z}: this model does not reach this redshift`);
  return Math.sqrt(e2);
}
/** H(z) in km s⁻¹ Mpc⁻¹. */
export const hubbleParameter = (z, p) => p.H0 * E(z, p);
/** Hubble time 1/H0 in Gyr. */
export const hubbleTime = p => HUBBLE_TIME_UNIT / p.H0;
/** Hubble distance c/H0 in Mpc. */
export const hubbleDistance = p => C_KMS / p.H0;

// ∫_lo^hi f(s) ds for smooth f, geometric panels towards lo (resolves the radiation→matter knee).
function quad(f, lo, hi) {
  let sum = 0, b = hi;
  while (b / 2 > lo && b > hi * 1e-12) { sum += gaussLegendre(f, b / 2, b, 24); b /= 2; }
  return sum + gaussLegendre(f, lo, b, 24);
}
// integrand 2 s^k / sqrt(Q(s²)); throws where Q ≤ 0 (expansion did not pass through this a)
function sIntegral(p, k, sLo, sHi, who) {
  return quad(s => {
    const q = Q(p, s * s);
    if (!(q > 0)) throw new Error(`${who}: E² ≤ 0 at a = ${s * s}: no expanding history through this scale factor`);
    return 2 * Math.pow(s, k) / Math.sqrt(q);
  }, sLo, sHi);
}
const hasBigBang = p => p.Or > 0 || p.Om > 0 || p.Ok > 0;

/**
 * Cosmic time since the Big Bang at redshift z (default 0: the age), Gyr: t = ∫₀^a da/(aH).
 * Throws for models without a Big Bang (Or = Om = 0 and Ok ≤ 0: de Sitter, bounces).
 * Checks: EdS 2/(3H0); Milne 1/H0; Planck 2018 13.787 ± 0.020 Gyr.
 */
export function age(p, z = 0) {
  checkModel(p, 'age'); checkZ(z, 'age');
  if (!hasBigBang(p)) throw new Error('age: this model has no Big Bang (Or = Om = 0, Ok ≤ 0): age is infinite');
  return sIntegral(p, 3, 0, Math.sqrt(1 / (1 + z)), 'age') * hubbleTime(p);
}

/** Lookback time to redshift z, Gyr: ∫_a^1 da/(aH). Finite for every model that reached z. */
export function lookbackTime(z, p) {
  checkModel(p, 'lookbackTime'); checkZ(z, 'lookbackTime');
  // Keep tiny redshifts in the integration variable: sqrt(1/(1+z)) can round to one.
  if (Math.abs(z) < 0.1) return gaussLegendre(zz => 1 / ((1 + zz) * E(zz, p)), 0, z, 24) * hubbleTime(p);
  const s = Math.sqrt(1 / (1 + z));
  return (z >= 0 ? 1 : -1) * sIntegral(p, 3, Math.min(s, 1), Math.max(s, 1), 'lookbackTime') * hubbleTime(p);
}

/** Line-of-sight comoving distance D_C = (c/H0)∫₀^z dz/E, Mpc (Hogg 1999 eq. 15). */
export function comovingDistance(z, p) {
  checkModel(p, 'comovingDistance'); checkZ(z, 'comovingDistance');
  if (Math.abs(z) < 0.1) return gaussLegendre(zz => 1 / E(zz, p), 0, z, 24) * hubbleDistance(p);
  const s = Math.sqrt(1 / (1 + z));
  return (z >= 0 ? 1 : -1) * sIntegral(p, 1, Math.min(s, 1), Math.max(s, 1), 'comovingDistance') * hubbleDistance(p);
}

/** Transverse comoving distance D_M (Hogg eq. 16): sinh/sin of D_C for open/closed curvature, Mpc. */
export function transverseComovingDistance(z, p) {
  const dc = comovingDistance(z, p), dh = hubbleDistance(p), k = Math.sqrt(Math.abs(p.Ok));
  if (p.Ok > 0) return dh / k * Math.sinh(k * dc / dh);
  if (p.Ok < 0) return dh / k * Math.sin(k * dc / dh);
  return dc;
}
/** Angular-diameter distance D_A = D_M/(1+z), Mpc (physical size / angle). */
export const angularDiameterDistance = (z, p) => transverseComovingDistance(z, p) / (1 + z);
/** Luminosity distance D_L = (1+z) D_M, Mpc (flux = L/(4π D_L²)). */
export const luminosityDistance = (z, p) => (1 + z) * transverseComovingDistance(z, p);
/** Distance modulus μ = 5 log₁₀(D_L / 10 pc), mag. z must be > 0. */
export function distanceModulus(z, p) {
  if (!(z > 0)) throw new Error('distanceModulus: z must be > 0');
  return 5 * Math.log10(luminosityDistance(z, p) * 1e5);
}

/**
 * Comoving particle horizon at redshift z (default 0), Mpc: the farthest comoving distance light
 * can have travelled since the Big Bang. Infinity when Or = Om = 0 (no horizon).
 */
export function particleHorizon(p, z = 0) {
  checkModel(p, 'particleHorizon'); checkZ(z, 'particleHorizon');
  if (!(p.Or > 0 || p.Om > 0)) return Infinity;
  return sIntegral(p, 1, 0, Math.sqrt(1 / (1 + z)), 'particleHorizon') * hubbleDistance(p);
}

/** Matter–radiation equality: 1 + z_eq = Om/Or. Throws if Or = 0. */
export function zMatterRadiation(p) {
  checkModel(p, 'zMatterRadiation');
  if (!(p.Or > 0)) throw new Error('zMatterRadiation: Or = 0, no radiation era');
  return p.Om / p.Or - 1;
}
/** Matter–Λ equality: 1 + z = (OL/Om)^{1/3}. Throws unless OL > 0 and Om > 0. */
export function zMatterLambda(p) {
  if (!(p.OL > 0 && p.Om > 0)) throw new Error('zMatterLambda: needs OL > 0 and Om > 0');
  return Math.cbrt(p.OL / p.Om) - 1;
}
/** Deceleration parameter q(z) = [Or(1+z)⁴ + Om(1+z)³/2 − OL] / E². */
export function decelerationParameter(z, p) {
  const x = 1 + z, e = E(z, p);
  return (p.Or * x ** 4 + p.Om * x ** 3 / 2 - p.OL) / (e * e);
}
/**
 * Deceleration→acceleration transition z where q = 0 (root by Brent on [−0.99, 1e4]).
 * Without radiation: 1 + z = (2 OL/Om)^{1/3}. Throws if OL ≤ 0 (never accelerates) or q < 0 throughout.
 */
export function zAcceleration(p) {
  checkModel(p, 'zAcceleration');
  if (!(p.OL > 0)) throw new Error('zAcceleration: OL ≤ 0, the expansion never accelerates');
  const g = z => { const x = 1 + z; return p.Or * x ** 4 + p.Om * x ** 3 / 2 - p.OL; };
  if (g(1e4) < 0) throw new Error('zAcceleration: accelerating at all z < 1e4');
  return brent(g, -0.99, 1e4, 1e-13);
}

/** CMB (any decoupled photon gas) temperature at redshift z: T = T0 (1+z), K. */
export function cmbTemperature(z, T0 = T_CMB) {
  checkZ(z, 'cmbTemperature');
  return T0 * (1 + z);
}

/**
 * Linear growth factor D(a) of matter perturbations (growing mode), normalised so D → a deep in
 * matter domination, by the Heath (1977) integral D = (5Om/2) E(a) ∫₀^a da'/(a'E)³.
 * Exact for matter + Λ + curvature; radiation is ignored here (error ~ a_eq/a at early times;
 * not a radiation-era growth solution). Requires Om > 0. EdS: D = a exactly.
 */
export function growthFactor(a, p) {
  checkModel(p, 'growthFactor');
  if (!(a > 0 && Number.isFinite(a))) throw new Error('growthFactor: a must be finite and > 0');
  if (!(p.Om > 0)) throw new Error('growthFactor: needs Om > 0');
  const pm = { ...p, Or: 0 };
  const q = s => { const v = Q(pm, s * s); if (!(v > 0)) throw new Error('growthFactor: E² ≤ 0 on the history'); return v; };
  const I = quad(s => 2 * s ** 7 / Math.pow(q(s), 1.5), 0, Math.sqrt(a));    // ∫ a³/Q^{3/2} da
  const Ea = Math.sqrt(Q(pm, a)) / (a * a);
  return 2.5 * p.Om * Ea * I;
}

/** Growth rate f = dln D/dln a from the Heath form (same conditions as growthFactor); EdS: f = 1. */
export function growthRate(a, p) {
  const D = growthFactor(a, p), pm = { ...p, Or: 0 };
  const e2 = Q(pm, a) / a ** 4;
  const dlnE = (-3 * p.Om / a ** 3 - 2 * p.Ok / a ** 2) / (2 * e2);
  // f = dlnE/dlna + (5Om/2)/(a² E² D)
  return dlnE + 2.5 * p.Om / (a * a * e2 * D);
}

/**
 * Carroll, Press & Turner (ARA&A 30, 499 (1992), eq. 29) fitting formula for the growth factor:
 * D ≈ a · (5/2)Ωm(a) / [Ωm(a)^{4/7} − ΩΛ(a) + (1 + Ωm(a)/2)(1 + ΩΛ(a)/70)]. APPROXIMATE (~1 %); same
 * normalisation as growthFactor. Radiation ignored.
 */
export function growthFactorCPT(a, p) {
  if (!(a > 0)) throw new Error('growthFactorCPT: a must be > 0');
  const e2 = p.Om / a ** 3 + p.Ok / a ** 2 + p.OL;
  const om = p.Om / a ** 3 / e2, ol = p.OL / e2;
  return a * 2.5 * om / (Math.pow(om, 4 / 7) - ol + (1 + om / 2) * (1 + ol / 70));
}

/**
 * Comoving sound horizon of the photon–baryon fluid at redshift z, Mpc:
 * r_s = ∫_z^∞ c_s dz/H, c_s = c/√(3(1+R)), R = 3ρ_b/(4ρ_γ) = (3Ωb/4Ωγ)a. Needs p.Obh2.
 * With Planck's z_drag = 1060.01 it reproduces r_drag = 147.21 Mpc to ~0.2 % (massive-ν approximation).
 */
export function soundHorizon(z, p) {
  checkModel(p, 'soundHorizon'); checkZ(z, 'soundHorizon');
  if (!(p.Obh2 > 0)) throw new Error('soundHorizon: the model needs Obh2 (baryon density)');
  if (!(p.Ogamma > 0 && p.Or >= p.Ogamma)) throw new Error('soundHorizon: needs a nonzero photon density included in the background radiation');
  const Rcoef = 3 * p.Obh2 / (4 * p.Ogamma * p.h * p.h);
  return quad(s => { const a = s * s; return 2 * s / Math.sqrt(Q(p, a)) / Math.sqrt(3 * (1 + Rcoef * a)); }, 0, Math.sqrt(1 / (1 + z))) * hubbleDistance(p);
}

/**
 * Sound horizon at the drag epoch, FITTING FORMULA of Aubourg et al., PRD 92, 123516 (2015), eq. 16
 * (0.021 % for Neff = 3.046 near Planck values): r_d = 55.154 exp[−72.3(ων+0.0006)²] / (ωcb^0.25351 ωb^0.12807) Mpc.
 * ωb = Ωb h², ωcb = (Ωc + Ωb) h², ων = Σmν/93.14 eV (default 0.06 eV → 0.000644).
 */
export function soundHorizonDragFit({ obh2, ocbh2, onuh2 = 0.06 / 93.14 }) {
  if (!(obh2 > 0 && ocbh2 > obh2 && onuh2 >= 0)) throw new Error('soundHorizonDragFit: need 0 < obh2 < ocbh2 and onuh2 ≥ 0');
  return 55.154 * Math.exp(-72.3 * (onuh2 + 0.0006) ** 2) / (Math.pow(ocbh2, 0.25351) * Math.pow(obh2, 0.12807));
}

/**
 * Solve the Friedmann acceleration equation ä = H0²[−Or/a³ − Om/(2a²) + OL a] for a(t) with dopri5
 * (constraint ȧ² = H0²Q(a)/a² sets the start). origin 'bang': start at aStart (default 1e-4) at
 * t = age(aStart), times since the Big Bang; origin 'today': a = 1 at t = 0 (for models without
 * a Big Bang, e.g. de Sitter). Runs to tMax (Gyr) or a Big Crunch (a falls back to aStart).
 * Returns { t: Float64Array (Gyr), a, H (km/s/Mpc), crunch: Gyr | null } with n uniform samples.
 * Checks: EdS a ∝ t^{2/3}, de Sitter e^{H0 t}, Milne H0 t, closed cycloid (tests).
 */
export function friedmannTrajectory(p, { tMax, n = 200, origin = 'bang', aStart = 1e-4, rtol = 1e-10 } = {}) {
  checkModel(p, 'friedmannTrajectory');
  if (!(Number.isFinite(tMax))) throw new Error('friedmannTrajectory: tMax (Gyr) is required');
  if (!(Number.isInteger(n) && n >= 2)) throw new Error('friedmannTrajectory: n must be an integer ≥ 2');
  if (origin !== 'bang' && origin !== 'today') throw new Error('friedmannTrajectory: origin must be "bang" or "today"');
  const tH = hubbleTime(p);
  let a0 = 1, t0 = 0;
  if (origin === 'bang') {
    if (!(aStart > 0 && aStart < 1)) throw new Error('friedmannTrajectory: 0 < aStart < 1');
    a0 = aStart; t0 = age(p, 1 / aStart - 1) / tH;
  }
  if (!(tMax / tH > t0)) throw new Error('friedmannTrajectory: tMax must be after the start time');
  const q0 = Q(p, a0);
  if (!(q0 > 0)) throw new Error('friedmannTrajectory: E² ≤ 0 at the start');
  const f = (_t, y) => [y[1], -p.Or / y[0] ** 3 - p.Om / (2 * y[0] * y[0]) + p.OL * y[0]];
  const floor = origin === 'bang' ? aStart : 1e-4;
  const st = odeStepper({ f, y0: [a0, Math.sqrt(q0) / a0], t0, h: 1e-3 * Math.max(t0, 1e-6), rtol, atol: 1e-14,
    events: [{ g: (_t, y) => y[0] - floor, direction: -1, terminal: true }] });
  const t = new Float64Array(n), a = new Float64Array(n), H = new Float64Array(n);
  const buf = new Float64Array(2);
  let m = 0;
  for (let i = 0; i < n; i++) {
    const ti = t0 + (tMax / tH - t0) * i / (n - 1);
    if (!st.advanceTo(ti, 1e6)) throw new Error('friedmannTrajectory: step budget exhausted');
    if (st.stopped && ti > st.t) break;
    st.sampleAt(ti, buf);
    t[m] = ti * tH; a[m] = buf[0]; H[m] = p.H0 * buf[1] / buf[0]; m++;
  }
  return { t: t.slice(0, m), a: a.slice(0, m), H: H.slice(0, m), crunch: st.stopped ? st.t * tH : null };
}
