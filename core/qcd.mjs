// Running strong coupling α_s(μ) in the MS-bar scheme, for any gauge group and flavour count.
//
//   const run = runAlphaS();                  // PDG 2024 α_s(M_Z) = 0.1180, 4 loops, c/b/t thresholds
//   run.at(10);                               // α_s at Q = 10 GeV;  run.nf(10) → 5
//
// Convention: a = α_s/(4π),  da/d ln μ² = −Σ_k β_k a^{k+2}  (k = 0..3), so for SU(3)
// β0 = 11 − 2nf/3, β1 = 102 − 38nf/3, β2 = 2857/2 − 5033nf/18 + 325nf²/54 and
// β3(nf = 5) ≈ 4826.16 (van Ritbergen, Vermaseren & Larin, Phys. Lett. B 400 (1997) 379).
//
// Flavour thresholds: nf switches at μ = m_q with α_s kept CONTINUOUS. The 2–4-loop MS-bar
// decoupling constants (Chetyrkin, Kniehl & Steinhauser), which make α_s jump at O(α_s³), are
// NOT applied; this simplification is not quantified here and pages must not claim its size.
//
// Perturbation theory fails below ~1 GeV (α_s ≳ 0.5): pages must shade that region, and
// runAlphaS throws below the Landau pole where the truncated series diverges.
//
// Exports: suN, betaCoefficients, runAlphaS, alphaS1Loop, lambdaQCD1Loop, exactTwoLoop.

import { dopri5, brent } from './numerics.mjs';
import { PDG } from './constants.mjs';

const ZETA3 = 1.2020569031595942;
const ALPHA_MAX = 100;                    // α_s above this is treated as the Landau pole

/**
 * Casimirs and quartic invariants of SU(N), fundamental quarks:
 * { CA, CF, TF, dAA, dFA, dFF } with d_XY = d_X^{abcd} d_Y^{abcd} / N_A.
 */
export function suN(N) {
  if (!(Number.isInteger(N) && N >= 2)) throw new Error('suN: N must be an integer >= 2');
  return {
    CA: N, CF: (N * N - 1) / (2 * N), TF: 1 / 2,
    dAA: N * N * (N * N + 36) / 24, dFA: N * (N * N + 6) / 48, dFF: (N ** 4 - 6 * N * N + 18) / (96 * N * N),
  };
}

const groupOf = g => {
  if (g === 'SU3' || g == null) return suN(3);
  if (typeof g === 'object' && [g.CA, g.CF, g.TF].every(Number.isFinite)) return g;
  throw new Error("betaCoefficients: group must be 'SU3' or { CA, CF, TF[, dAA, dFA, dFF] }");
};

/**
 * [β0, β1, β2, β3] in the convention above. β0–β2 from the Casimir formulas (β2: Tarasov,
 * Vladimirov & Zharkov 1980; Larin & Vermaseren 1993); β3 in Casimir form from van Ritbergen,
 * Vermaseren & Larin (1997), which needs dAA, dFA, dFF (given by suN); NaN if they are missing.
 */
export function betaCoefficients(nf, { group = 'SU3' } = {}) {
  if (!(Number.isFinite(nf) && nf >= 0)) throw new Error('betaCoefficients: nf must be a finite number >= 0');
  const { CA, CF, TF, dAA, dFA, dFF } = groupOf(group);
  const T = TF * nf, z = ZETA3;
  const b0 = 11 / 3 * CA - 4 / 3 * T;
  const b1 = 34 / 3 * CA * CA - 4 * CF * T - 20 / 3 * CA * T;
  const b2 = 2857 / 54 * CA ** 3 + 2 * CF * CF * T - 205 / 9 * CF * CA * T - 1415 / 27 * CA * CA * T
    + 44 / 9 * CF * T * T + 158 / 27 * CA * T * T;
  const b3 = ![dAA, dFA, dFF].every(Number.isFinite) ? NaN
    : CA ** 4 * (150653 / 486 - 44 / 9 * z) + dAA * (-80 / 9 + 704 / 3 * z)
    + CA ** 3 * T * (-39143 / 81 + 136 / 3 * z) + CA * CA * CF * T * (7073 / 243 - 656 / 9 * z)
    + CA * CF * CF * T * (-4204 / 27 + 352 / 9 * z) + nf * dFA * (512 / 9 - 1664 / 3 * z) + 46 * CF ** 3 * T
    + CA * CA * T * T * (7930 / 81 + 224 / 9 * z) + CF * CF * T * T * (1352 / 27 - 704 / 9 * z)
    + CA * CF * T * T * (17152 / 243 + 448 / 9 * z) + nf * nf * dFF * (-704 / 9 + 512 / 3 * z)
    + 424 / 243 * CA * T ** 3 + 1232 / 243 * CF * T ** 3;
  return [b0, b1, b2, b3];
}

/**
 * α_s(Q) by integrating the RG equation in t = ln μ² with dopri5 (rtol 1e-12), both up and
 * down from mu0. thresholds: quark masses (GeV) where nf steps by one above 3 light flavours
 * (continuous matching, see header); pass nf to run at fixed flavour number instead.
 * Accepted steps are cached; each query integrates only from the nearest cached node.
 * Returns { at(Q), nf(Q), loops, thresholds } where thresholds is the sorted mass list.
 */
export function runAlphaS({
  alpha0 = PDG.alphaS_MZ.value, mu0 = PDG.mZ.value, loops = 4,
  thresholds = { c: PDG.mc_mc.value, b: PDG.mb_mb.value, t: PDG.mt.value }, nf: nfFixed, group = 'SU3',
} = {}) {
  if (!(Number.isFinite(alpha0) && alpha0 > 0)) throw new Error('runAlphaS: alpha0 must be finite and > 0');
  if (!(Number.isFinite(mu0) && mu0 > 0)) throw new Error('runAlphaS: mu0 must be finite and > 0 (GeV)');
  if (![1, 2, 3, 4].includes(loops)) throw new Error('runAlphaS: loops must be 1, 2, 3 or 4');
  const masses = nfFixed != null ? [] : Object.values(thresholds ?? {}).map(Number).sort((x, y) => x - y);
  if (!masses.every(m => Number.isFinite(m) && m > 0)) throw new Error('runAlphaS: thresholds must be finite positive masses (GeV)');
  const nfLow = nfFixed ?? 3;
  const tb = masses.map(m => Math.log(m * m));
  const betas = tb.map((_, i) => nfLow + i).concat(nfLow + tb.length).map(n => betaCoefficients(n, { group }).slice(0, loops));
  if (betas.some(b => !b.every(Number.isFinite))) throw new Error('runAlphaS: β3 unknown for this group; use loops <= 3 or give dAA, dFA, dFF');
  const region = t => { let r = 0; while (r < tb.length && t >= tb[r]) r++; return r; };
  const t0 = Math.log(mu0 * mu0), r0 = region(t0), a0 = alpha0 / (4 * Math.PI);
  const nodes = [];                                    // per region: sorted [t, a] pairs

  const evalIn = (t, r) => {
    if (!nodes[r]) nodes[r] = [r === r0 ? [t0, a0] : r > r0 ? [tb[r - 1], evalIn(tb[r - 1], r - 1)] : [tb[r], evalIn(tb[r], r + 1)]];
    const list = nodes[r];
    let lo = 0, hi = list.length - 1;                  // nearest node by binary search
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (list[m][0] <= t) lo = m; else hi = m; }
    const [ts, as] = Math.abs(list[lo][0] - t) <= Math.abs(list[hi][0] - t) ? list[lo] : list[hi];
    if (ts === t) return as;
    const b = betas[r], added = [];
    const f = (_, y) => { const a = y[0]; let s = 0; for (let k = b.length - 1; k >= 0; k--) s = s * a + b[k]; return [-s * a * a]; };
    let last = ts, res;
    const pole = () => new Error(`runAlphaS: Q = ${Math.exp(t / 2).toPrecision(4)} GeV is beyond the Landau pole `
      + `(α_s diverges near μ ≈ ${Math.exp(last / 2).toPrecision(4)} GeV at ${loops} loops, nf = ${nfLow + r}); perturbative running is undefined there`);
    try {
      res = dopri5(f, ts, [as], t, {
        rtol: 1e-12, atol: 1e-30, h0: (t - ts) / 50,
        onStep: (tt, y) => { last = tt; if (!(y[0] > 0 && 4 * Math.PI * y[0] < ALPHA_MAX)) return false; added.push([tt, y[0]]); },
      });
    } catch { throw pole(); }
    if (res.stopped) throw pole();
    for (const p of added) list.push(p);
    list.sort((x, y) => x[0] - y[0]);
    return res.y[0];
  };

  const check = Q => { if (!(Number.isFinite(Q) && Q > 0)) throw new Error('runAlphaS: Q must be finite and > 0 (GeV)'); return Math.log(Q * Q); };
  return {
    at(Q) { const t = check(Q); return 4 * Math.PI * evalIn(t, region(t)); },
    nf(Q) { return nfLow + region(check(Q)); },
    loops, thresholds: masses,
  };
}

/** One-loop closed form α_s(Q) = 4π / (β0 ln(Q²/Λ²)); throws for Q <= Λ. */
export function alphaS1Loop(Q, Lambda, nf, { group = 'SU3' } = {}) {
  if (!(Lambda > 0 && Q > Lambda)) throw new Error('alphaS1Loop: need Q > Lambda > 0 (Landau pole at Q = Lambda)');
  return 4 * Math.PI / (betaCoefficients(nf, { group })[0] * Math.log(Q * Q / (Lambda * Lambda)));
}

/** One-loop Λ from α_s(μ): Λ = μ exp(−2π / (β0 α_s)). */
export function lambdaQCD1Loop(alpha, mu, nf, { group = 'SU3' } = {}) {
  if (!(alpha > 0 && mu > 0)) throw new Error('lambdaQCD1Loop: alpha and mu must be > 0');
  return mu * Math.exp(-2 * Math.PI / (betaCoefficients(nf, { group })[0] * alpha));
}

/**
 * Exact two-loop solution at fixed nf (for validating the integrator): with b = β1/β0,
 * 1/a + b ln(a/(1+ba)) = 1/a0 + b ln(a0/(1+b·a0)) + β0 ln(Q²/μ0²), solved for a by brent in ln a.
 */
export function exactTwoLoop(alpha0, mu0, Q, nf, { group = 'SU3' } = {}) {
  if (!(alpha0 > 0 && mu0 > 0 && Q > 0)) throw new Error('exactTwoLoop: alpha0, mu0 and Q must be > 0');
  const [b0, b1] = betaCoefficients(nf, { group }), b = b1 / b0;
  if (!(b0 > 0 && b > 0)) throw new Error('exactTwoLoop: needs β0 > 0 and β1 > 0');
  const a0 = alpha0 / (4 * Math.PI), F = a => 1 / a + b * Math.log(a / (1 + b * a));
  const C = F(a0) + b0 * Math.log(Q * Q / (mu0 * mu0));
  const g = x => { const a = Math.exp(x); return 1 / a - b * Math.log(1 / a + b) - C; };   // F(a) − C, overflow-safe
  let x;
  try { x = brent(g, Math.log(a0) - 60, Math.log(a0) + 60, 1e-15); } catch { throw new Error('exactTwoLoop: Q is beyond the Landau pole (no real solution)'); }
  return 4 * Math.PI * Math.exp(x);
}
