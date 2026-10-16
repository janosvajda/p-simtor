// Gravitational waves from a compact binary, to leading (quadrupole) order in the weak field.
// Geometric units G = c = 1 for the binary functions (masses, lengths and times in the same unit,
// usually the total mass M: 1 M = GM/c² of length or GM/c³ of time); SI helpers say so in their name.
//
// What this covers, and what it does not:
//   - the circular inspiral driven by the radiated energy (Peters 1964): exact solution of the
//     leading-order equations, valid while v² = M/a ≪ 1 (stop at the ISCO, a = 6M, or earlier);
//   - the leading-order strain h₊, h× far away (Maggiore 2008, §4.1);
//   - the eccentric-orbit period decay (Peters & Mathews 1963), the binary-pulsar test;
//   - the linearised tidal field E_ij = R_0i0j of the binary's rotating mass quadrupole in the near,
//     transition and wave zones (the weak-field "tendex" field, Nichols et al. 2011, PRD 84, 124014).
//   NOT covered: post-Newtonian corrections (≈ v² ≈ 10–20 % near the ISCO), spins, merger and
//   ring-down (numerical relativity), the monopole (Newtonian) part of the field.

import { GMsun, c as C_LIGHT } from './constants.mjs';

/** Chirp mass (m₁m₂)^{3/5} / (m₁+m₂)^{1/5}: the mass combination the inspiral signal measures. */
export const chirpMass = (m1, m2) => (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2;
/** Symmetric mass ratio η = m₁m₂/(m₁+m₂)², ¼ for equal masses. */
export const symmetricMassRatio = (m1, m2) => m1 * m2 / (m1 + m2) ** 2;

/** Power radiated by a circular binary of separation a: (32/5) m₁²m₂²(m₁+m₂)/a⁵ (G = c = 1). */
export const quadrupoleLuminosity = (a, m1, m2) => 32 / 5 * (m1 * m2) ** 2 * (m1 + m2) / a ** 5;

/**
 * Circular inspiral, exact solution of da/dt = −(64/5) m₁m₂M/a³ with Kepler's Ω = √(M/a³):
 *   a(t) = (a₀⁴ − 4βt)^{1/4},  β = (64/5) m₁m₂M,  φ(t) = (2√M / 5β)(a₀^{5/2} − a^{5/2}).
 * Returns functions of t (clamped at the end, tEnd = time to reach aEnd, default the ISCO 6M) plus
 * tCoalesce (a → 0 formally; the model fails long before). Phase φ is the orbital phase; the wave
 * phase is 2φ and the wave frequency 2Ω (f_gw = Ω/π).
 */
export function circularInspiral({ m1 = 0.5, m2 = 0.5, a0, aEnd = 6 * (m1 + m2) } = {}) {
  const M = m1 + m2, beta = 64 / 5 * m1 * m2 * M;
  if (!(a0 > aEnd)) throw new Error('a0 must exceed aEnd');
  const tEnd = (a0 ** 4 - aEnd ** 4) / (4 * beta);
  const a = t => (a0 ** 4 - 4 * beta * Math.min(Math.max(t, 0), tEnd)) ** 0.25;
  const phi = t => 2 * Math.sqrt(M) / (5 * beta) * (a0 ** 2.5 - a(t) ** 2.5);
  const omega = t => Math.sqrt(M / a(t) ** 3);
  return {
    m1, m2, M, mu: m1 * m2 / M, eta: m1 * m2 / M ** 2, a0, aEnd, beta, tEnd,
    tCoalesce: a0 ** 4 / (4 * beta),
    a, phi, omega,
    tOf: aa => (a0 ** 4 - aa ** 4) / (4 * beta),
    luminosity: t => quadrupoleLuminosity(a(t), m1, m2),
    /** Orbits left until the end, from time t. */
    orbitsLeft: t => (phi(tEnd) - phi(t)) / (2 * Math.PI),
    /**
     * Strain at distance r and inclination ι (0 = face-on), at retarded time t − r, from
     * h_ij = (2/r) Ï_ij^TT with φ measured from the x axis (the line of nodes): h₊ = h_xx, h× = h_xy,
     *   h₊ = −(4μM/(a r)) (1+cos²ι)/2 cos 2φ,  h× = −(4μM/(a r)) cos ι sin 2φ.
     * (Many texts drop the sign into a phase offset; it matters only when comparing with E_ij.)
     */
    strain(t, r, iota = 0) {
      const tr = t - r, A = -4 * m1 * m2 / (a(tr) * r), ci = Math.cos(iota), p = 2 * phi(tr);
      return { plus: A * (1 + ci * ci) / 2 * Math.cos(p), cross: A * ci * Math.sin(p) };
    },
  };
}

// ---------------------------------------------------------------- SI helpers (real systems)

/** Chirp rate df/dt (Hz/s) of the wave frequency f (Hz) for chirp mass Mc (M☉): (96/5) π^{8/3} (GMc/c³)^{5/3} f^{11/3}. */
export function chirpRateSI(f, McSol) {
  const tc = GMsun * McSol / C_LIGHT ** 3;   // GMc/c³ in s
  return 96 / 5 * Math.PI ** (8 / 3) * tc ** (5 / 3) * f ** (11 / 3);
}
/** Time (s) from wave frequency f (Hz) to formal coalescence, chirp mass Mc (M☉): (5/256) (GMc/c³)^{−5/3} (πf)^{−8/3}. */
export function timeToCoalescenceSI(f, McSol) {
  const tc = GMsun * McSol / C_LIGHT ** 3;
  return 5 / 256 * tc ** (-5 / 3) * (Math.PI * f) ** (-8 / 3);
}
/** Wave frequency (Hz) of a circular binary of total mass `Msol` (M☉) at separation a (m): Ω/π. */
export const gwFrequencySI = (a, Msol) => Math.sqrt(GMsun * Msol / a ** 3) / Math.PI;
/** Wave frequency (Hz) at the ISCO (a = 6GM/c²) of total mass `Msol` solar masses: c³/(6^{3/2} π G M). */
export const iscoFrequencySI = Msol => C_LIGHT ** 3 / (6 ** 1.5 * Math.PI * GMsun * Msol);

/** Peters–Mathews eccentricity enhancement (1 + 73e²/24 + 37e⁴/96) / (1 − e²)^{7/2}. */
export const eccentricityEnhancement = e => (1 + 73 / 24 * e * e + 37 / 96 * e ** 4) / (1 - e * e) ** 3.5;

/**
 * Orbital period decay dP_b/dt (dimensionless, s/s) of an eccentric binary (Peters & Mathews 1963):
 *   −(192π/5) (2π T☉/P_b)^{5/3} f(e) m_p m_c (m_p+m_c)^{−1/3},  masses in M☉, P_b in s, T☉ = GM☉/c³.
 */
export function periodDecay({ Pb, e, mp, mc }) {
  const Tsun = GMsun / C_LIGHT ** 3;
  return -192 * Math.PI / 5 * (2 * Math.PI * Tsun / Pb) ** (5 / 3) * eccentricityEnhancement(e) * mp * mc * (mp + mc) ** (-1 / 3);
}

// ---------------------------------------------------------------- linearised tidal field (tendexes)

const cr = new Float64Array(14), ci = new Float64Array(14), nr = new Float64Array(14), ni = new Float64Array(14);
const PR = new Float64Array(5), PI = new Float64Array(5);
// F_n(r) = ((1/r) d/dr)^n [e^{iωr}/r] with e^{iωr} removed: Σ_m c_m r^{−m}, n = 0…4, into PR/PI
function radial(r, w) {
  cr.fill(0); ci.fill(0); cr[1] = 1;
  const rp = 1 / r;
  for (let n = 0; n <= 4; n++) {
    let sr = 0, si = 0, p = rp;
    for (let m = 1; m <= 1 + 2 * n; m++) { sr += cr[m] * p; si += ci[m] * p; p *= rp; }
    PR[n] = sr; PI[n] = si;
    if (n === 4) break;
    nr.fill(0); ni.fill(0);
    for (let m = 1; m <= 1 + 2 * n; m++) {   // (1/r) d/dr of c_m r^{−m} e^{iωr}
      nr[m + 1] += -w * ci[m]; ni[m + 1] += w * cr[m];
      nr[m + 2] += -m * cr[m]; ni[m + 2] += -m * ci[m];
    }
    cr.set(nr); ci.set(ni);
  }
}

/**
 * Complex amplitude Ẽ_ij of the tidal tensor E_ij = R_0i0j of an oscillating, traceless mass
 * quadrupole I_ij(t) = Re[Q_ij e^{−iωt}] at the origin, at point (x, y, z) (linearised gravity,
 * outgoing waves, exact in all zones): E_ij(t) = Re[Ẽ_ij e^{−iωt}]. Q is given as row-major
 * 3×3 { re, im } (symmetric, traceless). Returns { re, im } Float64Array(9), row-major.
 */
export function quadrupoleTidalAmplitude(x, y, z, w, Q, out = { re: new Float64Array(9), im: new Float64Array(9) }) {
  radial(Math.hypot(x, y, z), w);
  const X = [x, y, z], Qr = Q.re, Qi = Q.im;
  const qr = [0, 1, 2].map(i => Qr[3 * i] * x + Qr[3 * i + 1] * y + Qr[3 * i + 2] * z);
  const qi = [0, 1, 2].map(i => Qi[3 * i] * x + Qi[3 * i + 1] * y + Qi[3 * i + 2] * z);
  const sr = x * qr[0] + y * qr[1] + z * qr[2], si = x * qi[0] + y * qi[1] + z * qi[2];
  const w2 = w * w, w4 = w2 * w2;
  for (let i = 0; i < 3; i++) {for (let j = 0; j < 3; j++) {
    const k = 3 * i + j, d = i === j;
    const xqr = X[i] * qr[j] + X[j] * qr[i], xqi = X[i] * qi[j] + X[j] * qi[i];
    // Ẽ = −ω⁴Q F0 − 2ω²Q F1 + (−ω²(xq+qx) − Q + ½ω²δ s) F2 + (−(xq+qx) − ½δ s) F3 − ½ x x s F4
    const c = [
      [-w4 * Qr[k], -w4 * Qi[k]],
      [-2 * w2 * Qr[k], -2 * w2 * Qi[k]],
      [-w2 * xqr - Qr[k] + (d ? 0.5 * w2 * sr : 0), -w2 * xqi - Qi[k] + (d ? 0.5 * w2 * si : 0)],
      [-xqr - (d ? 0.5 * sr : 0), -xqi - (d ? 0.5 * si : 0)],
      [-0.5 * X[i] * X[j] * sr, -0.5 * X[i] * X[j] * si],
    ];
    let re = 0, im = 0;
    for (let n = 0; n < 5; n++) { re += c[n][0] * PR[n] - c[n][1] * PI[n]; im += c[n][0] * PI[n] + c[n][1] * PR[n]; }
    out.re[k] = re; out.im[k] = im;
  }}
  return out;
}

/**
 * The tidal tensor E_ij (row-major Float64Array(9)) of a circular inspiral `ins` (from
 * circularInspiral, binary in the z = 0 plane) at time t and point (x, y, z): the rotating
 * quadrupole I = (μa²/2)·Re[(1, i; i, −1) e^{−2iφ}], evaluated at the retarded time t − r with the
 * slowly changing a and ω taken adiabatically (valid while the orbit changes little per wave period).
 * Units G = c = 1. Excludes the static monopole field M/r of the total mass.
 */
export function binaryTidalField(ins, t, x, y, z, out = new Float64Array(9)) {
  const r = Math.hypot(x, y, z), tr = t - r, a = ins.a(tr), amp = ins.mu * a * a / 2;
  const Q = { re: [amp, 0, 0, 0, -amp, 0, 0, 0, 0], im: [0, amp, 0, amp, 0, 0, 0, 0, 0] };
  const A = quadrupoleTidalAmplitude(x, y, z, 2 * ins.omega(tr), Q);
  const P = 2 * ins.phi(tr), c = Math.cos(P), s = Math.sin(P);
  // Re[(Ẽ evaluated with e^{iωr} removed) · e^{−i·2φ(t−r)}]: the retardation sits in φ(t − r)
  for (let k = 0; k < 9; k++) out[k] = c * A.re[k] + s * A.im[k];
  return out;
}

/** Radiated power from the far field, ∮ r² Σ|Ẽ_ij|²/(16πω²) dΩ (time average), by Gauss–Legendre × trapezoid quadrature. */
export function radiatedPowerFromField(w, Q, { r = 2e4 / w, nTheta = 24, nPhi = 16 } = {}) {
  // local Gauss–Legendre nodes on [−1, 1] (Newton on P_n)
  const xs = [], ws = [];
  for (let i = 1; i <= nTheta; i++) {
    let x = Math.cos(Math.PI * (i - 0.25) / (nTheta + 0.5)), dp = 0;
    for (let it = 0; it < 100; it++) {
      let p0 = 1, p1 = x;
      for (let k = 2; k <= nTheta; k++) { const p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
      dp = nTheta * (x * p1 - p0) / (x * x - 1);
      const dx = p1 / dp; x -= dx; if (Math.abs(dx) < 1e-16) break;
    }
    xs.push(x); ws.push(2 / ((1 - x * x) * dp * dp));
  }
  const buf = { re: new Float64Array(9), im: new Float64Array(9) };
  let P = 0;
  for (let i = 0; i < nTheta; i++) {
    const ct = xs[i], st = Math.sqrt(1 - ct * ct);
    for (let j = 0; j < nPhi; j++) {
      const ph = 2 * Math.PI * j / nPhi;
      const A = quadrupoleTidalAmplitude(r * st * Math.cos(ph), r * st * Math.sin(ph), r * ct, w, Q, buf);
      let s2 = 0; for (let k = 0; k < 9; k++) s2 += A.re[k] ** 2 + A.im[k] ** 2;
      P += ws[i] * (2 * Math.PI / nPhi) * r * r * s2 / (16 * Math.PI * w * w);
    }
  }
  return P;
}
