// Stellar structure and nuclear astrophysics: polytropes, white dwarfs, nuclear binding,
// ionisation, energy generation and the scaling laws of stars. SI units unless stated.
//
//   laneEmden(n)            θ'' + (2/ξ)θ' + θⁿ = 0, θ(0) = 1, θ'(0) = 0  (dopri5 + event at θ = 0)
//   polytrope({ n, M, R })  star of given mass and radius with P = Kρ^{1+1/n}; ideal-gas temperature
//   whiteDwarf({ xc })      Chandrasekhar's T = 0 degenerate-electron star (exact EOS, integrated)
//   chandrasekharMass(μe)   limiting mass (ω₃√(3π)/2)(ħc/G)^{3/2}/(μe m)²
//   semfBinding(Z, A)       Bethe–Weizsäcker binding energy (Rohlf coefficients)
//   sahaHydrogen(T, nH)     ionised fraction of pure hydrogen in LTE
//   epsPP, epsCNO           power-law energy generation rates (APPROXIMATE, labelled)
//   main-sequence scalings, Eddington luminosity, Kelvin–Helmholtz and nuclear timescales
//
// Solar values: IAU 2015 Resolution B3 nominal R☉ = 6.957e8 m, L☉ = 3.828e26 W, T_eff = 5772 K,
// GM☉ = 1.3271244e20 m³ s⁻² (M☉ = GM☉/G from constants.mjs). Other constants CODATA 2018.

import { G, c, h, hbar, kB, me, mp, u, a0, alpha, eV, Msun, yearJulian, NUCLIDE_MASS } from './constants.mjs';
import { brent, gaussLegendre } from './numerics.mjs';
import { odeStepper } from './ode.mjs';
import { SIGMA_SB } from './radiation.mjs';

/** IAU 2015 B3 nominal solar radius, m. */
export const R_SUN = 6.957e8;
/** IAU 2015 B3 nominal solar luminosity, W. */
export const L_SUN = 3.828e26;
/** IAU 2015 B3 nominal solar effective temperature, K (rounded; derived from L☉, R☉ and σ). */
export const TEFF_SUN = 5772;
/** Neutron mass in u (CODATA 2018). */
export const NEUTRON_MASS_U = 1.00866491595;
/** Atomic mass unit energy u c² in MeV (CODATA 2018). */
export const U_MEV = 931.49410242;
/** Thomson cross-section (8π/3) r_e², r_e = α² a₀, m² (CODATA 2018: 6.6524587321e-29). */
export const SIGMA_T = 8 * Math.PI / 3 * (alpha * alpha * a0) ** 2;
/** Hydrogen ionisation energy, eV (NIST ASD: 13.598434599702). */
export const CHI_H = 13.598434599702;

/**
 * Standard Solar Model reference values (BP2000: Bahcall, Pinsonneault & Basu, ApJ 555, 990 (2001)),
 * rounded; for comparison with the APPROXIMATE polytrope model only.
 */
export const SSM_CENTRE = Object.freeze({ Tc: 1.57e7, rhoC: 1.52e5, Pc: 2.34e16, source: 'BP2000 Standard Solar Model, ApJ 555, 990 (2001), rounded' });

/**
 * AME2020 binding energy per nucleon, keV (Wang et al., Chin. Phys. C 45, 030003 (2021)).
 * Ni-62 has the largest B/A of all nuclides, then Fe-58 and Fe-56.
 */
export const AME2020_BA = Object.freeze({
  H2: { Z: 1, A: 2, BA: 1112.2831 }, He4: { Z: 2, A: 4, BA: 7073.9156 }, C12: { Z: 6, A: 12, BA: 7680.1446 },
  O16: { Z: 8, A: 16, BA: 7976.2072 }, Ca40: { Z: 20, A: 40, BA: 8551.3046 }, Fe56: { Z: 26, A: 56, BA: 8790.3563 },
  Fe58: { Z: 26, A: 58, BA: 8792.2510 }, Ni62: { Z: 28, A: 62, BA: 8794.5532 }, Pb208: { Z: 82, A: 208, BA: 7867.4530 },
  U238: { Z: 92, A: 238, BA: 7570.1262 },
});

// ---------------------------------------------------------------- Lane–Emden

// cubic Hermite on stored nodes (x, value, derivative)
function hermiteTable(xs, ys, ds, x) {
  let lo = 0, hi = xs.length - 1;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[hi]) return ys[hi];
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  const hh = xs[hi] - xs[lo], t = (x - xs[lo]) / hh, t2 = t * t, t3 = t2 * t;
  return (2 * t3 - 3 * t2 + 1) * ys[lo] + (t3 - 2 * t2 + t) * hh * ds[lo] + (-2 * t3 + 3 * t2) * ys[hi] + (t3 - t2) * hh * ds[hi];
}

const LE_CACHE = new Map();
/**
 * Lane–Emden solution for polytropic index n (0 ≤ n ≤ 10), integrated with dopri5 (rtol default
 * 1e-12) from a series start at ξ = 1e-3 to the first zero ξ₁ (n < 5) or to xiMax (n ≥ 5: no zero).
 * Returns { n, xi1, dtheta1 = θ'(ξ₁), omega = −ξ₁²θ'(ξ₁) (mass constant), densityRatio = ρc/ρ̄
 * = −ξ₁/(3θ'(ξ₁)), Wn = 1/(4π(n+1)θ'(ξ₁)²) (Pc = Wn GM²/R⁴), theta(ξ), dtheta(ξ) }.
 * For n ≥ 5: xi1 = Infinity, densityRatio = Infinity, omega/dtheta1 evaluated at xiMax.
 * Exact: n = 0 (1 − ξ²/6, ξ₁ = √6), n = 1 (sin ξ/ξ, ξ₁ = π), n = 5 ((1 + ξ²/3)^{−1/2}).
 * Throws if n < 5 and no zero is found before xiMax (default 1e4).
 */
export function laneEmden(n, { rtol = 1e-12, xiMax } = {}) {
  if (!(Number.isFinite(n) && n >= 0 && n <= 10)) throw new Error('laneEmden: n must be in [0, 10]');
  if (!(rtol > 0 && rtol < 1e-3)) throw new Error('laneEmden: rtol must be in (0, 1e-3)');
  const xiEnd = xiMax ?? (n < 5 ? 1e4 : 50);
  if (!(xiEnd > 1e-2)) throw new Error('laneEmden: xiMax must be > 0.01');
  const key = `${n}|${rtol}|${xiEnd}`;
  if (LE_CACHE.has(key)) return LE_CACHE.get(key);
  const pw = th => (n === 0 ? 1 : Math.pow(Math.max(th, 0), n));
  const f = (xi, y) => [y[1], -pw(y[0]) - 2 * y[1] / xi];
  const x0 = 1e-3, k = n * (8 * n - 5);
  const y0 = [1 - x0 * x0 / 6 + n * x0 ** 4 / 120 - k * x0 ** 6 / 15120, -x0 / 3 + n * x0 ** 3 / 30 - k * x0 ** 5 / 2520];
  const xs = [0, x0], th = [1, y0[0]], dth = [0, y0[1]], d2 = [-1 / 3, f(x0, y0)[1]];
  const st = odeStepper({ f, y0, t0: x0, h: 1e-3, rtol, atol: 1e-15,
    events: n < 5 ? [{ g: (_x, y) => y[0], direction: -1, terminal: true }] : [] });
  let steps = 0;
  while (!st.stopped && st.t < xiEnd) {
    if (++steps > 1e6) throw new Error('laneEmden: step budget exhausted');
    st.step(xiEnd - st.t);
    xs.push(st.t); th.push(st.stopped ? 0 : st.y[0]); dth.push(st.y[1]); d2.push(st.dy[1]);
  }
  if (n < 5 && !st.stopped) throw new Error(`laneEmden: no zero of θ before ξ = ${xiEnd} (n = ${n} close to 5?)`);
  const xi1 = n < 5 ? st.t : Infinity, xe = st.t, d1 = st.y[1];
  const out = Object.freeze({
    n, xi1, dtheta1: d1, omega: -xe * xe * d1,
    densityRatio: n < 5 ? -xi1 / (3 * d1) : Infinity,
    Wn: n < 5 ? 1 / (4 * Math.PI * (n + 1) * d1 * d1) : Infinity,
    xiEnd: xe,
    theta: xi => (xi > xe ? NaN : hermiteTable(xs, th, dth, xi)),
    dtheta: xi => (xi > xe ? NaN : hermiteTable(xs, dth, d2, xi)),
  });
  LE_CACHE.set(key, out);
  return out;
}

/**
 * Polytropic star of index n (0 ≤ n < 5), mass M (kg) and radius R (m); ideal-gas temperature with
 * mean molecular weight mu (radiation pressure neglected). Returns { n, M, R, rhoMean, rhoC, Pc, Tc,
 * K, alpha (length scale), profile(r) → { rho, P, T, m } } (r in [0, R]).
 * APPROXIMATE as a model of a real star: the n = 3 Sun (Eddington's standard model) gives
 * ρc ≈ 76 g cm⁻³ and Tc ≈ 1.2e7 K vs 152 g cm⁻³ and 1.57e7 K in the Standard Solar Model.
 */
export function polytrope({ n, M, R, mu = 0.62 }) {
  if (!(n >= 0 && n < 5)) throw new Error('polytrope: n must be in [0, 5) for a finite radius');
  if (!(M > 0 && R > 0 && mu > 0 && Number.isFinite(M) && Number.isFinite(R))) throw new Error('polytrope: M, R, mu must be finite and > 0');
  const le = laneEmden(n);
  const rhoMean = 3 * M / (4 * Math.PI * R ** 3), rhoC = le.densityRatio * rhoMean;
  const Pc = le.Wn * G * M * M / R ** 4, Tc = mu * u * Pc / (kB * rhoC);
  const scale = R / le.xi1;
  const K = n === 0 ? NaN : Pc / Math.pow(rhoC, 1 + 1 / n);
  const profile = r => {
    if (!(r >= 0 && r <= R * (1 + 1e-12))) throw new Error('polytrope.profile: r must be in [0, R]');
    const xi = Math.min(r / scale, le.xi1), t = Math.max(le.theta(xi), 0);
    return { rho: rhoC * (n === 0 ? 1 : t ** n), P: Pc * t ** (n + 1), T: Tc * t, m: Math.min(M, 4 * Math.PI * scale ** 3 * rhoC * (-xi * xi * le.dtheta(xi))) };
  };
  return { n, M, R, mu, rhoMean, rhoC, Pc, Tc, K, alpha: scale, profile };
}

// ---------------------------------------------------------------- white dwarfs

/**
 * Chandrasekhar limiting mass M_Ch = ω₃ (√(3π)/2) (ħc/G)^{3/2} / (μe m)², kg, with ω₃ = 2.01824 from
 * laneEmden(3). mUnit: mass per unit μe (default the atomic mass unit u, so μe = 2 for He, C, O:
 * 1.4563 M☉ = 5.825/μe² M☉). With m = m_p one gets 1.4354 M☉ (older texts' "1.44").
 */
export function chandrasekharMass(mue = 2, { mUnit = u } = {}) {
  if (!(mue > 0 && mUnit > 0)) throw new Error('chandrasekharMass: mue and mUnit must be > 0');
  return laneEmden(3).omega * Math.sqrt(3 * Math.PI) / 2 * Math.pow(hbar * c / G, 1.5) / (mue * mUnit) ** 2;
}

// Degenerate electron gas: P = A[x(2x²−3)√(1+x²) + 3 asinh x], ρ = B x³, x = p_F/(m_e c)
const wdA = () => Math.PI * me ** 4 * c ** 5 / (3 * h ** 3);
const wdB = mue => mue * u * 8 * Math.PI / 3 * (me * c / h) ** 3;

/**
 * Zero-temperature white dwarf with the exact (relativistic, ideal, fully degenerate) electron EOS
 * (Chandrasekhar 1935). Give the central relativity parameter xc = p_F/(m_e c) or the central density
 * rhoC (kg m⁻³). Integrates y = √(1+x²): dy/dr' = −μ/r'², dμ/dr' = r'²(y²−1)^{3/2} to y = 1 (in w = y − 1).
 * Limits (tested): xc ≪ 1 → n = 3/2 polytrope; xc → ∞ → M_Ch, R → 0. Coulomb, inverse-β and GR
 * corrections are not included (a few % near M_Ch). Returns { mass (kg), radius (m), xc, rhoC, mue }.
 */
export function whiteDwarf({ xc, rhoC, mue = 2, rtol = 1e-11 } = {}) {
  if (!(mue > 0)) throw new Error('whiteDwarf: mue must be > 0');
  const A = wdA(), B = wdB(mue);
  const x = xc ?? (rhoC > 0 ? Math.cbrt(rhoC / B) : NaN);
  if (!(Number.isFinite(x) && x > 0 && x < 1e6)) throw new Error('whiteDwarf: give xc in (0, 1e6) or rhoC > 0');
  const L = Math.sqrt(2 * A / (Math.PI * G)) / B, M0 = 4 * Math.PI * B * L ** 3;
  // state w = y − 1 = x²/(√(1+x²) + 1) keeps full relative precision in the non-relativistic limit
  const wc = x * x / (Math.sqrt(1 + x * x) + 1), Lc = Math.sqrt(6 * wc / x ** 3), r0 = 1e-4 * Lc;
  const f = (r, s) => [-s[1] / (r * r), r * r * Math.pow(Math.max(s[0] * (s[0] + 2), 0), 1.5)];
  const st = odeStepper({ f, y0: [wc - x ** 3 * r0 * r0 / 6, x ** 3 * r0 ** 3 / 3], t0: r0, h: r0, rtol, atol: 1e-300,
    events: [{ g: (_r, s) => s[0], direction: -1, terminal: true }] });
  if (!st.advanceTo(1e3 * Lc, 1e6) || !st.stopped) throw new Error('whiteDwarf: surface not reached');
  return { mass: M0 * st.y[1], radius: L * st.t, xc: x, rhoC: B * x ** 3, mue };
}

/** White dwarf of given mass M (kg, below M_Ch): root of whiteDwarf(xc).mass = M in log xc. */
export function whiteDwarfOfMass(M, { mue = 2 } = {}) {
  const Mch = chandrasekharMass(mue);
  if (!(M > 0 && M < 0.999 * Mch)) throw new Error('whiteDwarfOfMass: need 0 < M < 0.999 M_Ch');
  const lx = brent(l => whiteDwarf({ xc: Math.exp(l), mue }).mass - M, Math.log(1e-6), Math.log(1e5), 1e-10);
  return whiteDwarf({ xc: Math.exp(lx), mue });
}

// ---------------------------------------------------------------- nuclei

/** Rohlf (Modern Physics from α to Z⁰, Wiley 1994) SEMF coefficients, MeV: volume, surface, Coulomb, asymmetry, pairing. */
export const SEMF_ROHLF = Object.freeze({ aV: 15.75, aS: 17.8, aC: 0.711, aA: 23.7, aP: 11.18 });

const checkZA = (Z, A, who) => {
  if (!(Number.isInteger(A) && A >= 1 && Number.isInteger(Z) && Z >= 0 && Z <= A)) throw new Error(`${who}: need integers A ≥ 1 and 0 ≤ Z ≤ A`);
};

/**
 * Bethe–Weizsäcker binding energy B(Z, A), MeV: aV A − aS A^{2/3} − aC Z(Z−1)/A^{1/3} − aA(A−2Z)²/A + δ,
 * δ = +aP/√A (even–even), −aP/√A (odd–odd), 0 (odd A). ~1 % for A ≳ 30; poor for light nuclei (He-4 −20 %).
 */
export function semfBinding(Z, A, k = SEMF_ROHLF) {
  checkZA(Z, A, 'semfBinding');
  const N = A - Z, pair = A % 2 === 1 ? 0 : (Z % 2 === 0 ? 1 : -1) * k.aP / Math.sqrt(A);
  return k.aV * A - k.aS * A ** (2 / 3) - k.aC * Z * (Z - 1) / Math.cbrt(A) - k.aA * (N - Z) ** 2 / A + pair;
}

/** Continuous Z minimising the atomic mass at fixed A (valley of stability), including m_n − m_H. */
export function mostStableZ(A, k = SEMF_ROHLF) {
  if (!(A >= 1)) throw new Error('mostStableZ: A must be ≥ 1');
  const dm = (NEUTRON_MASS_U - NUCLIDE_MASS.H1) * U_MEV, c3 = Math.cbrt(A);
  return (4 * k.aA + k.aC / c3 + dm) / (2 * k.aC / c3 + 8 * k.aA / A);
}

/** Binding energy from atomic masses (u): B = Z m(¹H) + N m_n − m(atom), MeV. */
export function bindingFromMass(Z, A, massU) {
  checkZA(Z, A, 'bindingFromMass');
  if (!(massU > 0)) throw new Error('bindingFromMass: mass (u) must be > 0');
  return (Z * NUCLIDE_MASS.H1 + (A - Z) * NEUTRON_MASS_U - massU) * U_MEV;
}

// ---------------------------------------------------------------- ionisation

/**
 * Saha ratio n_{i+1}/n_i = (2 U_{i+1}/U_i) / (n_e λ_e³) · exp(−χ/kT), λ_e = h/√(2π m_e kT).
 * T in K, ne in m⁻³, chi in eV, uRatio = U_{i+1}/U_i (hydrogen ground states: 1/2). LTE, ideal gas,
 * no pressure ionisation or continuum lowering.
 */
function logSahaRatio(T, ne, chi, uRatio) {
  if (![T, ne, chi, uRatio].every(v => Number.isFinite(v) && v > 0)) throw new Error('sahaRatio: T, ne, chi, uRatio must be finite and > 0');
  return Math.LN2 + Math.log(uRatio) - Math.log(ne)
    + 1.5 * (Math.log(2 * Math.PI * me * kB) + Math.log(T) - 2 * Math.log(h))
    - (chi * eV / kB) / T;
}
export function sahaRatio(T, ne, { chi = CHI_H, uRatio = 0.5 } = {}) {
  return Math.exp(logSahaRatio(T, ne, chi, uRatio));
}

/** Pure hydrogen of total nucleus density nH (m⁻³) at T (K): ionised fraction x from x²/(1−x) = S/nH. Returns { x, ne }. */
export function sahaHydrogen(T, nH) {
  const logS = logSahaRatio(T, nH, CHI_H, 0.5);          // log(x²/(1−x))
  // Solve using sqrt(s) or 1/s: s itself can underflow while x ≈ sqrt(s) is representable.
  const rootS = logS < 0 ? Math.exp(logS / 2) : 0;
  const x = logS < 0 ? 2 * rootS / (rootS + Math.hypot(rootS, 2)) : 2 / (1 + Math.sqrt(1 + 4 * Math.exp(-logS)));
  return { x, ne: x * nH };
}

/** Temperature (K) at which pure hydrogen of density nH (m⁻³) is 50 % ionised (Saha). */
export function halfIonisationTemperature(nH) {
  if (!(nH > 0 && nH < 1e30)) throw new Error('halfIonisationTemperature: nH must be in (0, 1e30) m⁻³');
  // x = 1/2 means S/nH = 1/2. The logarithmic equation also resolves very dilute gas,
  // whose half-ionisation temperature is below the old arbitrary 300 K bracket.
  return brent(T => logSahaRatio(T, nH, CHI_H, 0.5) + Math.LN2, 1, 1e6, 1e-10);
}

// ---------------------------------------------------------------- energy generation

const t6 = T => { if (!(T > 0)) throw new Error('T must be > 0'); return T / 1e6; };
/**
 * pp-chain power per unit mass, W kg⁻¹: ε ≈ 1.08e-12 ρ X² T₆⁴ (ρ in kg m⁻³). POWER-LAW APPROXIMATION
 * valid near T₆ ≈ 15 (Carroll & Ostlie, An Introduction to Modern Astrophysics, 2nd ed., §10.3).
 */
export const epsPP = (rho, X, T) => 1.08e-12 * rho * X * X * t6(T) ** 4;
/** CNO-cycle power per unit mass, W kg⁻¹: ε ≈ 8.24e-31 ρ X X_CNO T₆^{19.9}. POWER-LAW APPROXIMATION near T₆ ≈ 15 (same source). */
export const epsCNO = (rho, X, XCNO, T) => 8.24e-31 * rho * X * XCNO * t6(T) ** 19.9;

// ---------------------------------------------------------------- global stellar relations

/** L = 4πR²σT⁴, W (R in m, T in K). */
export function stefanBoltzmannLuminosity(R, T) {
  if (!(R > 0 && T > 0)) throw new Error('stefanBoltzmannLuminosity: R and T must be > 0');
  return 4 * Math.PI * R * R * SIGMA_SB * T ** 4;
}
/** T_eff = (L/(4πR²σ))^{1/4}, K. */
export function effectiveTemperature(L, R) {
  if (!(L > 0 && R > 0)) throw new Error('effectiveTemperature: L and R must be > 0');
  return Math.pow(L / (4 * Math.PI * R * R * SIGMA_SB), 0.25);
}

/**
 * Eddington luminosity 4πGMc/κ, W, electron-scattering opacity κ = σ_T(1+X)/(2m_p) (fully ionised,
 * hydrogen mass fraction X, default 1 → 1.257e31 W per M☉). Pass kappa (m² kg⁻¹) to override.
 */
export function eddingtonLuminosity(M, { X = 1, kappa } = {}) {
  if (!(M > 0)) throw new Error('eddingtonLuminosity: M must be > 0');
  const k = kappa ?? SIGMA_T * (1 + X) / (2 * mp);
  if (!(k > 0)) throw new Error('eddingtonLuminosity: opacity must be > 0');
  return 4 * Math.PI * G * M * c / k;
}

/** Kelvin–Helmholtz time q·GM²/(RL), s (q = 1 convention; a polytrope has q = 3/(5−n)). Sun: 31 Myr. */
export function kelvinHelmholtzTime(M, R, L, { q = 1 } = {}) {
  if (!(M > 0 && R > 0 && L > 0 && q > 0)) throw new Error('kelvinHelmholtzTime: M, R, L, q must be > 0');
  return q * G * M * M / (R * L);
}

/** Nuclear time f·ε·Mc²/L, s (ε = 0.007 for H → He, f = 0.1 of the mass burnt). Sun: ≈ 10 Gyr. */
export function nuclearTime(M, L, { efficiency = 0.007, coreFraction = 0.1 } = {}) {
  if (!(M > 0 && L > 0 && efficiency > 0 && coreFraction > 0)) throw new Error('nuclearTime: arguments must be > 0');
  return coreFraction * efficiency * M * c * c / L;
}

/**
 * Main-sequence mass–luminosity relation, L/L☉ for M/M☉ (APPROXIMATE piecewise power law, ~factor 2
 * scatter; e.g. Duric, Advanced Astrophysics (CUP 2004) §1.5): 0.23M^2.3 (M < 0.43), M⁴ (< 2),
 * 1.4M^3.5 (< 55), 32000M (above).
 */
export function mainSequenceLuminosity(m) {
  if (!(m > 0.08 && m < 300)) throw new Error('mainSequenceLuminosity: M/M☉ must be in (0.08, 300)');
  if (m < 0.43) return 0.23 * m ** 2.3;
  if (m < 2) return m ** 4;
  if (m < 55) return 1.4 * m ** 3.5;
  return 32000 * m;
}
/** Main-sequence radius R/R☉ ≈ (M/M☉)^0.8 (APPROXIMATE, ~30 % scatter, for 0.1–20 M☉). */
export function mainSequenceRadius(m) {
  if (!(m > 0.08 && m < 300)) throw new Error('mainSequenceRadius: M/M☉ must be in (0.08, 300)');
  return m ** 0.8;
}
/** Main-sequence lifetime ≈ 10 Gyr (M/M☉)/(L/L☉), Gyr (scaled from the Sun; APPROXIMATE). */
export const mainSequenceLifetime = m => 10 * m / mainSequenceLuminosity(m);

/** Seconds per Julian year (for converting timescales). */
export const YEAR = yearJulian;
/** Solar mass, kg (GM☉/G; re-exported for convenience). */
export const M_SUN = Msun;

/** ∫₀^R 4πr²ρ dr of a polytrope by Gauss–Legendre (consistency check of the profile), kg. */
export function polytropeMass(star, n = 64) {
  return gaussLegendre(r => 4 * Math.PI * r * r * star.profile(r).rho, 0, star.R, n);
}
