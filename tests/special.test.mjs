import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  airy, airyAiZeros, sphericalBesselJ, sphericalBesselY, besselY, besselK, hermiteH, oscillatorWavefunction,
  hyp2f1, hyp1f1, ellipticE, jacobiElliptic, digamma, beta, logBeta, incompleteBeta, studentTCdf, studentTPValue,
  zeta, polylog, expintE1, expintEn, expintEi, dawson, wigner3j, wigner6j, clebschGordan,
} from '../core/special.mjs';
import { gaussLegendre, legendreP, besselJ, besselI, ellipticK, laguerre, gammaP, erf, lgamma, gamma } from '../core/numerics.mjs';

const close = (a, b, tol, msg = '') => assert.ok(Math.abs(a - b) <= tol, `${msg} ${a} vs ${b} (tol ${tol})`);
const rel = (a, b, tol, msg = '') => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg} ${a} vs ${b} (rel ${tol})`);
/** Composite 32-point Gauss–Legendre on [a, b] with `panels` panels. */
const quad = (f, a, b, panels = 20) => { let s = 0; const w = (b - a) / panels; for (let i = 0; i < panels; i++) s += gaussLegendre(f, a + i * w, a + (i + 1) * w, 32); return s; };
const EULER = 0.5772156649015329;

// ------------------------------------------------------------------ Airy

test('Airy: A&S Table 10.11 values, Wronskian 1/π everywhere, continuity across method boundaries', () => {
  const table = [   // x, Ai, Ai′, Bi, Bi′  (A&S 10.4.4–5 and Table 10.11)
    [0, 0.3550280538878172, -0.2588194037928068, 0.6149266274460008, 0.4482883573538264],
    [1, 0.13529241631288141, -0.15914744129679328, 1.2074235949528713, 0.9324359333927756],
    [-1, 0.5355608832923521, -0.01016056711664521, 0.10399738949694461, 0.5923756264227924],
    [2, 0.03492413042327438, -0.053090384433653715, 3.2980949999782148, 4.1006820499328915],
  ];
  for (const [x, ai, aip, bi, bip] of table) {
    const r = airy(x);
    rel(r.ai, ai, 1e-14, `Ai(${x})`); rel(r.aip, aip, 1e-13, `Ai'(${x})`); rel(r.bi, bi, 1e-14, `Bi(${x})`); rel(r.bip, bip, 1e-14, `Bi'(${x})`);
  }
  rel(airy(10).ai, 1.1047532552898687e-10, 1e-13, 'Ai(10)');
  for (let x = -60; x <= 60; x += 0.37) {
    const r = airy(x), w = (r.ai * r.bip - r.aip * r.bi) * Math.PI;
    close(w, 1, x < 0 ? 1e-13 * Math.sqrt(Math.abs(x)) : 1e-13, `Wronskian at ${x}`);
  }
  for (const b of [-10, -1, 1, 10]) {
    const lo = airy(b - 1e-13), hi = airy(b + 1e-13), env = b < 0 ? Math.pow(-b, 0.25) : 1;
    for (const k of ['ai', 'aip', 'bi', 'bip']) close(lo[k], hi[k], 1e-12 * Math.max(env, Math.abs(hi[k])), `${k} at ${b}`);
  }
});

test('Airy: Ai″ = x Ai (derivative consistency) and Ai(x) = (1/π)√(x/3) K_{1/3}(ζ) by quadrature', () => {
  for (const x of [-12, -7.3, -3, -0.4, 0.6, 3.3, 8, 11]) {
    const h = 1e-4, d = (airy(x + h).aip - airy(x - h).aip) / (2 * h), r = airy(x);
    close(d, x * r.ai, 1e-7 * Math.max(1, Math.abs(x)), `Ai'' at ${x}`);
    const db = (airy(x + h).bip - airy(x - h).bip) / (2 * h);
    close(db, x * r.bi, 1e-7 * Math.max(1, Math.abs(x * r.bi)), `Bi'' at ${x}`);
  }
  for (const x of [0.5, 2, 5, 9, 15]) {
    const zt = 2 / 3 * x ** 1.5;
    const K13 = quad(t => Math.exp(-zt * Math.cosh(t)) * Math.cosh(t / 3), 0, 8, 40);
    rel(airy(x).ai, Math.sqrt(x / 3) * K13 / Math.PI, 1e-13, `Ai(${x}) vs integral`);
  }
});

test('Airy zeros: A&S Table 10.13 (a1…a5, Ai′(a1)), no zero skipped', () => {
  const ref = [-2.338107410459767, -4.087949444130971, -5.520559828095551, -6.786708090071759, -7.944133587120853];
  const z = airyAiZeros(30);
  ref.forEach((v, i) => rel(z[i], v, 1e-14, `a${i + 1}`));
  rel(airy(z[0]).aip, 0.7012108227206914, 1e-13, "Ai'(a1)");
  for (let i = 0; i < z.length; i++) close(airy(z[i]).ai, 0, 1e-14 * Math.pow(-z[i], -0.25) * 10, `Ai(a${i + 1})`);
  // count sign changes on a fine grid from 0 down to a30 − 0.1
  let changes = 0, prev = airy(0).ai;
  for (let x = -0.01; x > z[29] - 0.1; x -= 0.01) { const v = airy(x).ai; if (v * prev < 0) changes++; prev = v; }
  assert.equal(changes, 30);
  assert.throws(() => airyAiZeros(0), /count/);
  assert.throws(() => airy(NaN), /finite/);
});

// ------------------------------------------------------------------ Bessel family

test('spherical Bessel j_n, y_n: closed forms, small x, Rayleigh integral, Σ(2n+1)j_n² = 1, Wronskian', () => {
  for (const x of [0.01, 0.3, 1, 2.5, 7, 30]) {
    const s = Math.sin(x), c = Math.cos(x);
    rel(sphericalBesselJ(2, x), (3 / x ** 3 - 1 / x) * s - 3 * c / x ** 2, x < 0.1 ? 1e-6 : 1e-13, `j2(${x})`);   // closed form cancels at small x
    rel(sphericalBesselY(2, x), (-3 / x ** 3 + 1 / x) * c - 3 * s / x ** 2, 1e-13, `y2(${x})`);
  }
  // small x: j_n ≈ xⁿ/(2n+1)!! (1 − x²/(2(2n+3)))
  for (const n of [0, 3, 10, 100]) {
    const x = 1e-3;
    let pre = 1;
    for (let k = 1; k <= n; k++) pre *= x / (2 * k + 1);
    rel(sphericalBesselJ(n, x), pre * (1 - x * x / (2 * (2 * n + 3))), 1e-13, `j_${n}(1e-3)`);
  }
  // Rayleigh: j_n(x) = ½ ∫₋₁¹ cos(xt − nπ/2) P_n(t) dt
  for (const n of [0, 1, 5, 10, 30]) {
    for (const x of [0.7, 5, 20, 45]) {
      const I = 0.5 * quad(t => Math.cos(x * t - n * Math.PI / 2) * legendreP(n, 0, t), -1, 1, 20);
      close(sphericalBesselJ(n, x), I, 1e-14, `Rayleigh j_${n}(${x})`);
    }
  }
  for (const x of [0.5, 3, 17.3, 80]) {
    let sum = 0;
    for (let n = 0; n < 200; n++) sum += (2 * n + 1) * sphericalBesselJ(n, x) ** 2;
    close(sum, 1, 1e-13, `Σ(2n+1)j_n²(${x})`);
    for (const n of [0, 1, 7, 40]) {      // DLMF 10.50.3: j_{n+1} y_n − j_n y_{n+1} = x^{−2}
      const w = sphericalBesselJ(n + 1, x) * sphericalBesselY(n, x) - sphericalBesselJ(n, x) * sphericalBesselY(n + 1, x);
      rel(w, 1 / (x * x), 1e-12, `Wronskian n=${n} x=${x}`);
    }
  }
  close(sphericalBesselJ(3, -2), -sphericalBesselJ(3, 2), 1e-16);
  close(sphericalBesselY(2, -2), -sphericalBesselY(2, 2), 1e-16);
  assert.equal(sphericalBesselJ(0, 0), 1); assert.equal(sphericalBesselJ(4, 0), 0); assert.equal(sphericalBesselY(1, 0), -Infinity);
});

test('Bessel Y_n: A&S Table 9.1 values and the Wronskian with J_n (both sides of x = 25)', () => {
  rel(besselY(0, 1), 0.08825696421567696, 1e-13);
  rel(besselY(1, 1), -0.7812128213002887, 1e-14);
  rel(besselY(0, 10), 0.055671167283599395, 1e-12);
  rel(besselY(1, 10), 0.24901542420695388, 1e-13);
  rel(besselY(2, 1), -1.6506826068162543, 1e-14);
  for (const x of [1e-6, 0.01, 0.9, 5, 13.7, 24.9, 25.1, 60, 300]) {
    for (const n of [0, 1, 2, 7, 20]) {   // DLMF 10.5.2: J_{n+1}Y_n − J_n Y_{n+1} = 2/(πx)
      const w = besselJ(n + 1, x) * besselY(n, x) - besselJ(n, x) * besselY(n + 1, x);
      if (Number.isFinite(w)) rel(w, 2 / (Math.PI * x), 1e-12, `Y Wronskian n=${n} x=${x}`);
    }
  }
  close(besselY(0, 25 - 1e-14), besselY(0, 25 + 1e-14), 1e-14, 'Neumann/Hankel switch');
  close(besselY(1, 25 - 1e-14), besselY(1, 25 + 1e-14), 1e-14, 'Neumann/Hankel switch');
  rel(besselY(-3, 2), -besselY(3, 2), 1e-15);
  assert.equal(besselY(0, 0), -Infinity);
  assert.throws(() => besselY(0, -1), /complex/);
  assert.throws(() => besselY(0.5, 1), /integer/);
});

test('modified Bessel K_n: A&S Table 9.8 values, integral representation, Wronskian with I_n', () => {
  rel(besselK(0, 1), 0.42102443824070834, 1e-14);
  rel(besselK(1, 1), 0.6019072301972346, 1e-14);
  rel(besselK(0, 2), 0.11389387274953344, 1e-14);
  rel(besselK(1, 2), 0.13986588181652243, 1e-14);
  for (const x of [0.1, 0.5, 1.999, 2.001, 3, 10, 40]) {
    for (const n of [0, 1, 3]) {        // K_n(x) = ∫₀^∞ e^{−x cosh t} cosh(nt) dt
      const T = Math.acosh(Math.max(1, 800 / x));
      rel(besselK(n, x), quad(t => Math.exp(-x * Math.cosh(t)) * Math.cosh(n * t), 0, T, 60), 1e-13, `K_${n}(${x})`);
    }
    if (x <= 30) for (const n of [0, 1, 4]) rel(besselI(n, x) * besselK(n + 1, x) + besselI(n + 1, x) * besselK(n, x), 1 / x, 1e-13, `I/K Wronskian ${n} ${x}`);
  }
  rel(besselK(-2, 1.5), besselK(2, 1.5), 0);
  assert.equal(besselK(0, 0), Infinity);
  assert.equal(besselK(0, 800), 0);
});

// ------------------------------------------------------------------ Hermite

test('Hermite H_n and oscillator wavefunctions: closed forms, ψ_n(0), orthonormality to n = 100, n = 2000 normalised', () => {
  for (const x of [-1.3, 0, 0.5, 2]) {
    close(hermiteH(3, x), 8 * x ** 3 - 12 * x, 1e-12);
    close(hermiteH(4, x), 16 * x ** 4 - 48 * x * x + 12, 1e-12);
  }
  assert.equal(hermiteH(10, 1), 8224);                                 // A&S Table 22.12 (H_10(1))
  const psi = oscillatorWavefunction;
  for (const m of [0, 1, 5, 50]) {        // ψ_{2m}(0) = π^{−1/4} (−1)^m √((2m)!) / (2^m m!)
    const v = Math.pow(Math.PI, -0.25) * (m % 2 ? -1 : 1) * Math.exp(0.5 * lgamma(2 * m + 1) - m * Math.LN2 - lgamma(m + 1));
    rel(psi(2 * m, 0), v, 1e-13, `ψ_${2 * m}(0)`);
    close(psi(2 * m + 1, 0), 0, 1e-300);
  }
  for (const [m, n] of [[0, 0], [3, 3], [3, 5], [40, 40], [40, 41], [99, 100], [100, 100]]) {
    const I = quad(x => psi(m, x) * psi(n, x), -22, 22, 120);
    close(I, m === n ? 1 : 0, 1e-12, `<${m}|${n}>`);
  }
  const N = quad(x => psi(2000, x) ** 2, -72, 72, 500);
  close(N, 1, 1e-10, 'norm n=2000');
  const far = psi(1000, 60);                                         // deep in the forbidden region, no overflow
  assert.ok(far !== 0 && Math.abs(far) < 1e-100 && Number.isFinite(far), `${far}`);
  // physical check: ⟨x²⟩ = n + 1/2 for ψ_n
  close(quad(x => x * x * psi(30, x) ** 2, -22, 22, 120), 30.5, 1e-10);
});

// ------------------------------------------------------------------ hypergeometric

test('₂F₁: elementary closed forms in every region, Gauss sum, K(m), Legendre polynomials', () => {
  for (const x of [-5, -0.9, -0.3, 0.3, 0.7, 0.95]) rel(hyp2f1(1, 1, 2, x), -Math.log1p(-x) / x, 1e-13, `log at ${x}`);
  for (const x of [-30, -3, -0.6, 0.2, 0.6, 0.9, 0.999]) rel(hyp2f1(0.3, 1.7, 1.7, x), Math.pow(1 - x, -0.3), 1e-13, `binomial at ${x}`);
  for (const x of [0.2, 0.6, 0.9, 0.99]) rel(hyp2f1(0.5, 0.5, 1.5, x * x), Math.asin(x) / x, 1e-13, `asin at ${x}`);
  for (const x of [0.5, 1, 3, 30]) rel(hyp2f1(0.5, 1, 1.5, -x * x), Math.atan(x) / x, 1e-13, `atan at ${x}`);
  rel(hyp2f1(0.3, 0.7, 2.5, 1), gamma(2.5) * gamma(1.5) / (gamma(2.2) * gamma(1.8)), 1e-14, 'Gauss');
  for (const m of [0.1, 0.5, 0.9]) rel(Math.PI / 2 * hyp2f1(0.5, 0.5, 1, m), ellipticK(m), 1e-13, `K(${m})`);
  // non-integer c − a − b connection formula vs. direct series at x = 0.6 (both paths valid)
  rel(hyp2f1(0.37, 1.21, 2.9, 0.75), hyp2f1(0.37, 1.21, 2.9, 0.75 - 1e-12), 1e-11);
  // terminating series (accuracy relative to its largest term: here ~1e4 × eps at x = −0.8, n = 10)
  for (const n of [0, 3, 10]) for (const x of [-0.8, 0.1, 0.95]) close(hyp2f1(-n, n + 1, 1, (1 - x) / 2), legendreP(n, 0, x), 1e-10, `P_${n}(${x})`);
  assert.throws(() => hyp2f1(1, 1, 2, 1.5), /branch cut/);
  assert.throws(() => hyp2f1(1, 1, -2, 0.3), /pole/);
  assert.throws(() => hyp2f1(1, 1, 2, 1), /diverges/);
  assert.throws(() => hyp2f1(1, 1, 2, 1 - 1e-9), /too close to 1/);
});

test('₁F₁: closed forms, erf and incomplete Γ relations (series, Kummer and asymptotic paths), Laguerre', () => {
  for (const x of [-40, -5, 0.3, 5, 40, 200]) {
    rel(hyp1f1(0.7, 0.7, x), Math.exp(x), 1e-13, `exp ${x}`);
    rel(hyp1f1(1, 2, x), Math.expm1(x) / x, 1e-13, `(e^x−1)/x ${x}`);
  }
  for (const x of [0.5, 2, 4, 10, 60]) rel(hyp1f1(0.5, 1.5, -x * x), Math.sqrt(Math.PI) * erf(x) / (2 * x), 1e-13, `erf ${x}`);
  for (const [a, x] of [[2.5, 1], [2.5, 10], [7.3, 40], [0.4, 80]]) {   // P(a, x) = xᵃ e^{−x} M(1, a+1, x) / Γ(a+1)
    rel(Math.pow(x, a) * Math.exp(-x) * hyp1f1(1, a + 1, x) / gamma(a + 1), gammaP(a, x), 1e-12, `gammaP ${a} ${x}`);
  }
  for (const [n, al] of [[3, 0.5], [8, 2], [12, 0]]) {for (const x of [0.2, 3, 9]) {
    const binom = Math.exp(lgamma(n + al + 1) - lgamma(n + 1) - lgamma(al + 1));
    close(binom * hyp1f1(-n, al + 1, x), laguerre(n, al, x), 1e-11 * Math.max(1, Math.abs(laguerre(n, al, x))), `L_${n}^${al}(${x})`);
  }}
  rel(hyp1f1(0.3, 2.2, -7), Math.exp(-7) * hyp1f1(1.9, 2.2, 7), 1e-13, 'Kummer');
  assert.throws(() => hyp1f1(1, -2, 1), /pole/);
});

// ------------------------------------------------------------------ elliptic

test('complete elliptic E(m): A&S Table 17.1, Legendre relation, quadrature', () => {
  close(ellipticE(0), Math.PI / 2, 1e-16); assert.equal(ellipticE(1), 1);
  rel(ellipticE(0.5), 1.3506438810476755, 1e-15);
  for (const m of [0.01, 0.1, 0.5, 0.9, 0.999999]) {
    const K = ellipticK(m), Kp = ellipticK(1 - m), E = ellipticE(m), Ep = ellipticE(1 - m);
    close(E * Kp + Ep * K - K * Kp, Math.PI / 2, 1e-14, `Legendre relation m=${m}`);   // A&S 17.3.13
    rel(E, quad(t => Math.sqrt(1 - m * Math.sin(t) ** 2), 0, Math.PI / 2, 200), 1e-14, `E(${m})`);
  }
  rel(ellipticE(-2), quad(t => Math.sqrt(1 + 2 * Math.sin(t) ** 2), 0, Math.PI / 2, 4), 1e-14);
  assert.throws(() => ellipticE(1.5), /≤ 1/);
});

test('Jacobi elliptic sn, cn, dn, am: inversion of F(φ|m), identities, quarter period, derivative, limits', () => {
  const F = (phi, m) => quad(t => 1 / Math.sqrt(1 - m * Math.sin(t) ** 2), 0, phi, 60);
  for (const m of [-3, -0.5, 0.1, 0.5, 0.9, 0.999, 2.5]) {
    const phis = m > 1 ? [0.2, 0.5] : [0.3, 1.2, 2.9, 7];
    for (const phi of phis) {
      const u = F(phi, m), r = jacobiElliptic(u, m);
      close(r.sn, Math.sin(phi), 2e-14, `sn m=${m} φ=${phi}`);
      close(r.cn, Math.cos(phi), 2e-14, `cn m=${m} φ=${phi}`);
      close(r.dn, Math.sqrt(1 - m * Math.sin(phi) ** 2), 2e-14, `dn m=${m} φ=${phi}`);
      if (m <= 1) close(r.am, phi, 2e-14, `am m=${m} φ=${phi}`);
    }
    for (const u of [0.37, 1.9, -4.2]) {
      const r = jacobiElliptic(u, m), h = 1e-5;
      close(r.sn ** 2 + r.cn ** 2, 1, 1e-15); close(r.dn ** 2 + m * r.sn ** 2, 1, 2e-15);
      const d = (jacobiElliptic(u + h, m).sn - jacobiElliptic(u - h, m).sn) / (2 * h);
      close(d, r.cn * r.dn, 1e-9, `sn' m=${m}`);
    }
    if (m < 1) {
      const K = ellipticK(m), q = jacobiElliptic(K, m);
      close(q.sn, 1, 1e-15); close(q.cn, 0, 1e-15); close(q.dn, Math.sqrt(1 - m), 1e-15);
      close(jacobiElliptic(0.7 + 4 * K, m).sn, jacobiElliptic(0.7, m).sn, 1e-14, 'period 4K');
      close(jacobiElliptic(0.7 + 2 * K, m).am, jacobiElliptic(0.7, m).am + Math.PI, 1e-14, 'am(u+2K)');
    }
  }
  const c = jacobiElliptic(0.8, 0);
  close(c.sn, Math.sin(0.8), 1e-16); close(c.cn, Math.cos(0.8), 1e-16);
  const hyp = jacobiElliptic(0.8, 1);
  close(hyp.sn, Math.tanh(0.8), 1e-16); close(hyp.dn, 1 / Math.cosh(0.8), 1e-16);
  close(jacobiElliptic(0.8, 1 - 1e-12).sn, Math.tanh(0.8), 1e-11);
});

// ------------------------------------------------------------------ Γ-related

test('digamma: DLMF 5.4 special values, recurrence, reflection, derivative of ln Γ, positive root', () => {
  close(digamma(1), -EULER, 1e-15);
  close(digamma(0.5), -EULER - 2 * Math.LN2, 1e-15);
  close(digamma(0.25), -EULER - Math.PI / 2 - 3 * Math.LN2, 1e-15);
  let H = 0;
  for (let n = 1; n <= 30; n++) { H += 1 / n; close(digamma(n + 1), H - EULER, 1e-14); }
  for (const x of [-7.3, -2.5, -0.2, 0.01, 1.3, 7.7, 40, 1e6]) {
    close(digamma(x + 1), digamma(x) + 1 / x, 1e-13 * Math.max(1, Math.abs(digamma(x))), `recurrence ${x}`);
    if (x > 0.5 && x < 100) close((lgamma(x + 1e-5) - lgamma(x - 1e-5)) / 2e-5, digamma(x), 1e-8, `lnΓ' ${x}`);
  }
  close(digamma(1 - 0.3) - digamma(0.3), Math.PI / Math.tan(0.3 * Math.PI), 1e-14);
  close(digamma(1.4616321449683622), 0, 2e-16);                    // DLMF 5.4.vi
  assert.throws(() => digamma(-3), /pole/);
});

test('beta, incomplete beta and Student t: exact cases, binomial tails, t tables', () => {
  rel(beta(0.5, 0.5), Math.PI, 1e-14);
  rel(beta(4, 6), 6 * 120 / 362880, 1e-14);
  rel(beta(300, 1), 1 / 300, 1e-12);
  close(logBeta(2, 3), Math.log(1 / 12), 1e-14);
  for (const x of [0, 0.1, 0.5, 0.93, 1]) {
    close(incompleteBeta(x, 2.5, 1), x ** 2.5, 1e-15);
    close(incompleteBeta(x, 1, 3.5), 1 - (1 - x) ** 3.5, 1e-15);
    close(incompleteBeta(x, 2.2, 3.7), 1 - incompleteBeta(1 - x, 3.7, 2.2), 1e-15);
  }
  // binomial tail P(X ≥ k) = I_p(k, n − k + 1)
  const n = 20, p = 0.3;
  for (const k of [1, 5, 12, 20]) {
    let tail = 0;
    for (let j = k; j <= n; j++) tail += Math.exp(lgamma(n + 1) - lgamma(j + 1) - lgamma(n - j + 1)) * p ** j * (1 - p) ** (n - j);
    rel(incompleteBeta(p, k, n - k + 1), tail, 1e-12, `binomial k=${k}`);
  }
  for (const t of [-30, -1, 0, 0.5, 4]) {
    close(studentTCdf(t, 1), 0.5 + Math.atan(t) / Math.PI, 1e-15, `Cauchy ${t}`);
    close(studentTCdf(t, 2), 0.5 + t / (2 * Math.sqrt(2 + t * t)), 1e-15, `ν=2 ${t}`);
  }
  close(studentTPValue(2.228138851986274, 10), 0.05, 1e-13);       // t_{0.975,10} (A&S Table 26.10)
  close(studentTPValue(Math.tan(0.475 * Math.PI), 1), 0.05, 1e-14);   // t_{0.975,1} = tan(0.475π)
  // ν → ∞: normal quantile z_{0.975} plus the leading 1/ν correction (t³ + t)φ(t)/(2ν) of the two-sided tail
  const z = 1.959963984540054, corr = (z ** 3 + z) * Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI) / 2;
  for (const nu of [1e4, 1e6]) close(studentTPValue(z, nu), 0.05 + corr / nu, 3 / nu ** 2 + 1e-13, `ν = ${nu}`);
  close(studentTPValue(2.228138851986274, 10, 1), 0.025, 1e-13);
  rel(studentTPValue(40, 5), 2 * (1 - studentTCdf(40, 5)), 1e-6);
  assert.ok(studentTPValue(1e4, 3) > 0 && studentTPValue(1e4, 3) < 1e-11);
  assert.throws(() => incompleteBeta(1.2, 1, 1), /0 ≤ x ≤ 1/);
});

// ------------------------------------------------------------------ ζ and polylog

test('zeta: ζ(2) = π²/6, ζ(4) = π⁴/90, Apéry ζ(3), ζ(3/2), negative integers, pole behaviour', () => {
  rel(zeta(2), Math.PI ** 2 / 6, 1e-15);
  rel(zeta(4), Math.PI ** 4 / 90, 1e-15);
  rel(zeta(3), 1.2020569031595942, 1e-15);
  rel(zeta(1.5), 2.612375348685488, 1e-15);
  rel(zeta(0.5), -1.4603545088095868, 1e-15);
  rel(zeta(-0.5), -0.20788622497735457, 1e-14);
  assert.equal(zeta(0), -0.5);
  rel(zeta(-1), -1 / 12, 1e-14); rel(zeta(-3), 1 / 120, 1e-14); rel(zeta(-7), 1 / 240, 1e-14);
  assert.equal(zeta(-2), 0); assert.equal(zeta(-40), 0);
  const s1 = 1 + 1e-7;
  rel(zeta(s1) * (s1 - 1), 1 + EULER * (s1 - 1), 1e-14, 'Laurent ζ ≈ 1/(s−1) + γ');
  rel(zeta(60), 1 + 2 ** -60, 1e-16);
  assert.ok(Number.isFinite(zeta(-201.5)));
  assert.throws(() => zeta(1), /pole/);
});

test('polylog: Li₂ special values and reflection, Li₃(½), Bose ζ(3/2), duplication, inversion, Li₁, Li₀', () => {
  rel(polylog(2, 0.5), Math.PI ** 2 / 12 - Math.LN2 ** 2 / 2, 1e-15);
  rel(polylog(2, -1), -(Math.PI ** 2) / 12, 1e-15);
  rel(polylog(3, 0.5), 7 / 8 * zeta(3) - Math.PI ** 2 / 12 * Math.LN2 + Math.LN2 ** 3 / 6, 1e-15);
  rel(polylog(1.5, 1), 2.612375348685488, 1e-15);
  for (const z of [0.1, 0.3, 0.6, 0.9, 0.99]) rel(polylog(2, z) + polylog(2, 1 - z), Math.PI ** 2 / 6 - Math.log(z) * Math.log(1 - z), 1e-14, `Euler reflection ${z}`);
  for (const z of [-0.9, 0.2, 0.7]) rel(polylog(1, z), -Math.log1p(-z), 1e-15);
  for (const z of [-3, -0.7, 0.4, 0.8]) { rel(polylog(0, z), z / (1 - z), 1e-14); rel(polylog(-1, z), z / (1 - z) ** 2, 1e-13); }
  for (const s of [0.5, 2.5, -1.5]) {for (const z of [0.6, 0.8, 0.95]) {   // duplication with z < 0 evaluated by the alternating series (s ≤ 1)
    rel(polylog(s, z) + polylog(s, -z), 2 ** (1 - s) * polylog(s, z * z), 1e-13, `duplication s=${s} z=${z}`);
  }}
  for (const z of [-1.5, -3, -40]) {   // Li₂(z) = −∫₀^z ln(1−u)/u du
    rel(polylog(2, z), quad(u => -Math.log1p(-u) / u, z, 0, 40) * -1, 1e-13, `Li2(${z})`);
  }
  rel(polylog(3, -5) - polylog(3, -1 / 5), -(Math.PI ** 2) / 6 * Math.log(5) - Math.log(5) ** 3 / 6, 1e-14);
  rel(polylog(2.5, -1), -(1 - 2 ** -1.5) * zeta(2.5), 1e-15);
  close(polylog(1.5, -0.999999), polylog(1.5, -1), 1e-6);
  // series vs μ-expansion on both sides of z = ½
  rel(polylog(2.7, 0.5), polylog(2.7, 0.5 + 1e-13), 1e-12);
  assert.throws(() => polylog(2, 1.5), /complex/);
  assert.throws(() => polylog(2.5, -3), /integer/);
  assert.throws(() => polylog(0.5, 1), /diverges/);
});

// ------------------------------------------------------------------ exponential integrals

test('exponential integrals: A&S Table 5.1 values, quadrature, recurrence, Ei root and derivative', () => {
  rel(expintE1(1), 0.21938393439552029, 1e-15);
  rel(expintEi(1), 1.8951178163559368, 1e-15);
  rel(expintE1(10), 0.000004156968929685325, 1e-14);
  rel(expintEi(10), 2492.2289762418777, 1e-14);
  rel(expintE1(0.01), 4.03792957653811, 2e-15);
  for (const x of [0.05, 0.5, 1, 3, 20]) {
    for (const n of [1, 2, 5]) {       // E_n(x) = ∫₀¹ e^{−x/u} u^{n−2} du
      rel(expintEn(n, x), quad(u => (u === 0 ? 0 : Math.exp(-x / u) * u ** (n - 2)), 0, 1, 200), 1e-12, `E_${n}(${x}) quadrature`);
      rel(n * expintEn(n + 1, x), Math.exp(-x) - x * expintEn(n, x), 1e-13, `recurrence ${n} ${x}`);
    }
  }
  assert.equal(expintEn(3, 0), 0.5);
  close(expintEi(0.3725074107813666), 0, 1e-15);              // DLMF 6.13
  for (const x of [-3, 0.2, 2, 39.99, 40.01, 100]) {
    const h = 1e-6 * Math.max(1, Math.abs(x));
    rel((expintEi(x + h) - expintEi(x - h)) / (2 * h), Math.exp(x) / x, 1e-7, `Ei' ${x}`);
  }
  rel(expintEi(40 - 1e-14), expintEi(40 + 1e-14), 1e-13, 'series/asymptotic switch');
  close(expintEi(-2), -expintE1(2), 0);
  assert.throws(() => expintE1(-1), /complex|≥ 0/);
});

// ------------------------------------------------------------------ Dawson

test('Dawson F(x): erfi series, ODE F′ = 1 − 2xF, maximum, asymptotics, continuity', () => {
  const erfiForm = x => {           // e^{−x²} Σ x^{2k+1}/(k!(2k+1)) (positive terms)
    let term = x, sum = x;
    for (let k = 1; k < 300; k++) { term *= x * x / k; sum += term / (2 * k + 1); if (term < 1e-18 * sum) break; }
    return Math.exp(-x * x) * sum;
  };
  for (const x of [0.1, 0.5, 0.999, 1.5, 2, 3.3, 5]) rel(dawson(x), erfiForm(x), 1e-14, `F(${x})`);
  rel(dawson(1), 0.5380795069127684, 1e-15);
  for (const x of [0.3, 1.7, 6, 30]) {
    const h = 1e-5, d = (dawson(x + h) - dawson(x - h)) / (2 * h);
    close(d, 1 - 2 * x * dawson(x), 1e-9, `ODE at ${x}`);
  }
  const xm = 0.9241388730045917;
  close(2 * xm * dawson(xm), 1, 1e-15, 'max');
  rel(dawson(xm), 0.5410442246351816, 1e-15);
  const r = 1 / (2 * 100 ** 2);
  rel(dawson(100), (1 + r + 3 * r * r + 15 * r ** 3 + 105 * r ** 4) / 200, 1e-15);
  for (const b of [1, 50]) rel(dawson(b - 1e-12), dawson(b + 1e-12), 1e-12);
  close(dawson(-2), -dawson(2), 0);
});

// ------------------------------------------------------------------ Wigner symbols

test('Wigner 3j / Clebsch–Gordan: tabulated values (exact), orthogonality, symmetry, selection rules', () => {
  rel(wigner3j(1, 1, 0, 0, 0, 0), -1 / Math.sqrt(3), 1e-15);
  assert.deepEqual(wigner3j(1, 1, 0, 0, 0, 0, { exact: true }), { sign: -1, num: 1n, den: 3n });
  assert.deepEqual(wigner3j(2, 2, 2, 0, 0, 0, { exact: true }), { sign: -1, num: 2n, den: 35n });
  for (const j of [0.5, 1, 2.5, 7]) for (let m = -j; m <= j; m++) rel(wigner3j(j, j, 0, m, -m, 0), ((j - m) % 2 ? -1 : 1) / Math.sqrt(2 * j + 1), 1e-15);
  close(clebschGordan(0.5, 0.5, 0.5, -0.5, 1, 0), Math.SQRT1_2, 1e-16);
  close(clebschGordan(0.5, 0.5, 0.5, -0.5, 0, 0), Math.SQRT1_2, 1e-16);
  close(clebschGordan(0.5, -0.5, 0.5, 0.5, 0, 0), -Math.SQRT1_2, 1e-16);
  close(clebschGordan(1, 1, 1, 0, 2, 1), Math.SQRT1_2, 1e-16);
  rel(clebschGordan(1, 1, 1, -1, 0, 0), 1 / Math.sqrt(3), 1e-15);
  rel(clebschGordan(1, 0, 1, 0, 0, 0), -1 / Math.sqrt(3), 1e-15);
  rel(clebschGordan(1, 1, 0.5, -0.5, 1.5, 0.5), Math.sqrt(1 / 3), 1e-15);
  // CG unitarity for j1 = 2, j2 = 3/2
  const j1 = 2, j2 = 1.5;
  for (let J = 0.5; J <= 3.5; J++) {for (let Jp = 0.5; Jp <= 3.5; Jp++) {for (let M = -Math.min(J, Jp); M <= Math.min(J, Jp); M++) {
    let s = 0;
    for (let m1 = -j1; m1 <= j1; m1++) { const m2 = M - m1; if (Math.abs(m2) <= j2) s += clebschGordan(j1, m1, j2, m2, J, M) * clebschGordan(j1, m1, j2, m2, Jp, M); }
    close(s, J === Jp ? 1 : 0, 1e-14, `CG orthogonality ${J} ${Jp} ${M}`);
  }}}
  // 3j orthogonality at large j (exact arithmetic keeps it to rounding)
  let s = 0;
  for (let m1 = -40; m1 <= 40; m1++) s += 121 * wigner3j(40, 35, 60, m1, 5 - m1, -5) ** 2;      // (2j3 + 1) Σ (3j)² = 1
  close(s, 1, 1e-13, 'orthogonality j=40');
  // symmetries: even column permutation, odd with (−1)^{j1+j2+j3}, m → −m
  const v = wigner3j(3, 2.5, 1.5, -1, 0.5, 0.5);
  close(wigner3j(2.5, 1.5, 3, 0.5, 0.5, -1), v, 1e-16);
  close(wigner3j(2.5, 3, 1.5, 0.5, -1, 0.5), -v, 1e-16);
  close(wigner3j(3, 2.5, 1.5, 1, -0.5, -0.5), -v, 1e-16);
  assert.equal(wigner3j(1, 1, 3, 0, 0, 0), 0);      // triangle
  assert.equal(wigner3j(1, 1, 1, 0, 0, 0), 0);      // odd J with all m = 0
  assert.equal(wigner3j(1, 1, 1, 1, 1, 0), 0);      // m sum
  assert.throws(() => wigner3j(0.3, 1, 1, 0, 0, 0), /half-integer/);
});

test('Wigner 6j: closed forms, {1 1 1; 1 1 1} = 1/6, orthogonality', () => {
  for (const [a, b, c] of [[1, 1, 1], [0.5, 1, 1.5], [2, 3, 4], [2.5, 1.5, 3]]) {   // {a b c; 0 c b} = (−1)^{a+b+c}/√((2b+1)(2c+1))
    rel(wigner6j(a, b, c, 0, c, b), ((a + b + c) % 2 ? -1 : 1) / Math.sqrt((2 * b + 1) * (2 * c + 1)), 1e-15);
  }
  assert.deepEqual(wigner6j(1, 1, 1, 1, 1, 1, { exact: true }), { sign: 1, num: 1n, den: 36n });
  close(wigner6j(0.5, 0.5, 1, 0.5, 0.5, 1), 1 / 6, 1e-16);
  // Σ_x (2x+1)(2f+1) {a b x; c d f}{a b x; c d f′} = δ_{ff′}
  const [a, b, c, d] = [2, 1.5, 2.5, 1];
  for (const f of [1, 2, 3]) {for (const fp of [1, 2, 3]) {
    let s = 0;
    for (let x = 0.5; x <= 3.5; x++) s += (2 * x + 1) * (2 * f + 1) * wigner6j(a, b, x, c, d, f) * wigner6j(a, b, x, c, d, fp);
    close(s, f === fp ? 1 : 0, 1e-14, `6j orthogonality ${f} ${fp}`);
  }}
  assert.equal(wigner6j(1, 1, 3, 1, 1, 1), 0);
});
