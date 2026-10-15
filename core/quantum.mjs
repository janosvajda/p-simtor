// Finite-dimensional quantum systems: spins, qubits, excitons, few-level molecules. Complex dense
// matrices, tensor products, spin operators for any j, a complex Hermitian eigensolver, exact
// propagators and open-system (Lindblad) dynamics of the density matrix. Pure functions.
//
//   cmat(n), identity(n), fromReal(rows), cadd, cscale, cmul, dagger, kron, trace, expectation
//   spinOps(j)                        { x, y, z, plus, minus } for spin j (ħ = 1), dimension 2j + 1
//   embed(op, k, dims)                op acting on subsystem k of a tensor product with these dims
//   eigHermitian(H)                   complex Jacobi: ascending real values, unitary V (columns = eigenvectors)
//   propagator(eig, t, hbar)          U = V e^{−iEt/ħ} V† (exact for time-independent H)
//   evolve(rho, U)                    U ρ U†
//   lindblad({ H, jumps, rho0 })      stepper for dρ/dt = −i[H, ρ]/ħ + Σ γ(LρL† − ½{L†L, ρ}), adaptive dopri5
//   thermalState(H, kT)               e^{−H/kT}/Z (detailed balance reference)
//   checkDensity(rho)                 { trace, hermiticity, minEigenvalue } for verify()
//
// A complex matrix is { n, re, im } with row-major Float64Arrays of length n². Units: ħ = 1 unless
// hbar is passed; energies and rates then share one unit (e.g. cm⁻¹ with time in cm, or rad/µs).

import { odeStepper } from './ode.mjs';

export const cmat = n => ({ n, re: new Float64Array(n * n), im: new Float64Array(n * n) });
export function identity(n) { const A = cmat(n); for (let i = 0; i < n; i++) A.re[i * n + i] = 1; return A; }
/** From real rows (and optional imaginary rows). */
export function fromReal(rows, imRows) {
  const n = rows.length, A = cmat(n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { A.re[i * n + j] = rows[i][j]; if (imRows) A.im[i * n + j] = imRows[i][j]; }
  return A;
}
const same = (A, B, who) => { if (A.n !== B.n) throw new Error(`${who}: matrices must have the same size`); };
export function cadd(A, B, b = 1) { same(A, B, 'cadd'); const C = cmat(A.n); for (let i = 0; i < A.re.length; i++) { C.re[i] = A.re[i] + b * B.re[i]; C.im[i] = A.im[i] + b * B.im[i]; } return C; }
/** s·A for a real s or a complex [sRe, sIm]. */
export function cscale(A, s) {
  const [sr, si] = typeof s === 'number' ? [s, 0] : s, C = cmat(A.n);
  for (let i = 0; i < A.re.length; i++) { C.re[i] = sr * A.re[i] - si * A.im[i]; C.im[i] = sr * A.im[i] + si * A.re[i]; }
  return C;
}
export function cmul(A, B) {
  same(A, B, 'cmul');
  const n = A.n, C = cmat(n);
  for (let i = 0; i < n; i++) {for (let k = 0; k < n; k++) {
    const ar = A.re[i * n + k], ai = A.im[i * n + k];
    if (ar === 0 && ai === 0) continue;
    for (let j = 0; j < n; j++) {
      const br = B.re[k * n + j], bi = B.im[k * n + j];
      C.re[i * n + j] += ar * br - ai * bi; C.im[i * n + j] += ar * bi + ai * br;
    }
  }}
  return C;
}
export function dagger(A) { const n = A.n, C = cmat(n); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { C.re[j * n + i] = A.re[i * n + j]; C.im[j * n + i] = -A.im[i * n + j]; } return C; }
export function kron(A, B) {
  const n = A.n * B.n, C = cmat(n);
  for (let i = 0; i < A.n; i++) {for (let j = 0; j < A.n; j++) {
    const ar = A.re[i * A.n + j], ai = A.im[i * A.n + j];
    if (ar === 0 && ai === 0) continue;
    for (let k = 0; k < B.n; k++) {for (let l = 0; l < B.n; l++) {
      const br = B.re[k * B.n + l], bi = B.im[k * B.n + l], o = (i * B.n + k) * n + j * B.n + l;
      C.re[o] = ar * br - ai * bi; C.im[o] = ar * bi + ai * br;
    }}
  }}
  return C;
}
export function trace(A) { let r = 0, i = 0; for (let k = 0; k < A.n; k++) { r += A.re[k * A.n + k]; i += A.im[k * A.n + k]; } return { re: r, im: i }; }
/** Tr(ρ O): the expectation value of O in the state ρ (real for Hermitian O). */
export function expectation(rho, O) {
  same(rho, O, 'expectation');
  const n = rho.n; let r = 0;
  for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) r += rho.re[i * n + k] * O.re[k * n + i] - rho.im[i * n + k] * O.im[k * n + i];
  return r;
}

/** Spin operators S_x, S_y, S_z, S₊, S₋ for spin j (integer or half-integer), basis m = j, j−1, …, −j. */
export function spinOps(j) {
  const d = Math.round(2 * j + 1);
  if (!(d >= 2 && Math.abs(d - (2 * j + 1)) < 1e-12)) throw new Error('spinOps: j must be ½, 1, 3/2, …');
  const z = cmat(d), plus = cmat(d), minus = cmat(d);
  for (let a = 0; a < d; a++) {
    const m = j - a;
    z.re[a * d + a] = m;
    if (a > 0) plus.re[(a - 1) * d + a] = Math.sqrt(j * (j + 1) - m * (m + 1));     // S₊|m⟩ = √(j(j+1) − m(m+1)) |m+1⟩
  }
  for (let a = 0; a < d; a++) for (let b = 0; b < d; b++) minus.re[a * d + b] = plus.re[b * d + a];
  const x = cscale(cadd(plus, minus), 0.5), y = cscale(cadd(plus, minus, -1), [0, -0.5]);
  return { x, y, z, plus, minus, dim: d };
}

/** The operator op on subsystem k of a tensor product of subsystems with dimensions dims. */
export function embed(op, k, dims) {
  let M = null;
  dims.forEach((d, i) => { const f = i === k ? op : identity(d); M = M ? kron(M, f) : f; });
  return M;
}

/**
 * Eigen-decomposition of a complex Hermitian matrix by cyclic complex Jacobi rotations: each
 * off-diagonal element is first made real by a phase on its column, then removed by a real
 * rotation. Handles degenerate eigenvalues (vectors stay orthonormal). Returns { values (ascending),
 * vectors } where vectors is a cmat whose column i is the eigenvector of values[i].
 */
export function eigHermitian(H, { rtol = 1e-14, maxSweeps = 100 } = {}) {
  const n = H.n, re = Float64Array.from(H.re), im = Float64Array.from(H.im);
  for (let i = 0; i < n; i++) {for (let j = i; j < n; j++) {
    const dr = re[i * n + j] - re[j * n + i], di = im[i * n + j] + im[j * n + i];
    if (Math.hypot(dr, di) > 1e-10 * (1 + Math.abs(re[i * n + j]) + Math.abs(im[i * n + j]))) throw new Error('eigHermitian: matrix must be Hermitian');
  }}
  for (let i = 0; i < n; i++) im[i * n + i] = 0;
  const V = identity(n), vr = V.re, vi = V.im;
  let scale = 0; for (let i = 0; i < n * n; i++) scale = Math.hypot(scale, re[i], im[i]);
  const off = () => { let s = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) s = Math.hypot(s, re[i * n + j], im[i * n + j]); return s; };
  for (let sweep = 0; off() > rtol * (scale || 1); sweep++) {
    if (sweep === maxSweeps) throw new Error('eigHermitian: no convergence');
    for (let p = 0; p < n; p++) {for (let q = p + 1; q < n; q++) {
      const cr = re[p * n + q], ci = im[p * n + q], mag = Math.hypot(cr, ci);
      if (mag <= 1e-300 || mag < 1e-18 * (scale || 1)) continue;
      // 1. phase e^{−iφ} on column q (and e^{iφ} on row q) makes A_pq real = |c|
      const ur = cr / mag, ui = -ci / mag;                    // e^{−iφ}
      for (let k = 0; k < n; k++) {
        let a = re[k * n + q], b = im[k * n + q];               // column q × e^{−iφ}
        re[k * n + q] = a * ur - b * ui; im[k * n + q] = a * ui + b * ur;
        a = re[q * n + k]; b = im[q * n + k];                   // row q × e^{iφ}
        re[q * n + k] = a * ur + b * ui; im[q * n + k] = -a * ui + b * ur;
        a = vr[k * n + q]; b = vi[k * n + q];
        vr[k * n + q] = a * ur - b * ui; vi[k * n + q] = a * ui + b * ur;
      }
      im[q * n + q] = 0;
      // 2. real Jacobi rotation in the (p, q) plane
      const app = re[p * n + p], aqq = re[q * n + q], delta = (aqq - app) / 2;
      const t = delta === 0 ? 1 : mag / (delta + Math.sign(delta) * Math.hypot(delta, mag));
      const c = 1 / Math.hypot(1, t), s = t * c;
      re[p * n + p] = app - t * mag; re[q * n + q] = aqq + t * mag;
      re[p * n + q] = re[q * n + p] = im[p * n + q] = im[q * n + p] = 0;
      for (let k = 0; k < n; k++) {
        if (k !== p && k !== q) {
          const ur2 = re[k * n + p], ui2 = im[k * n + p], wr = re[k * n + q], wi = im[k * n + q];
          re[k * n + p] = c * ur2 - s * wr; im[k * n + p] = c * ui2 - s * wi;
          re[k * n + q] = s * ur2 + c * wr; im[k * n + q] = s * ui2 + c * wi;
          re[p * n + k] = re[k * n + p]; im[p * n + k] = -im[k * n + p];
          re[q * n + k] = re[k * n + q]; im[q * n + k] = -im[k * n + q];
        }
        const ar = vr[k * n + p], ai = vi[k * n + p], br = vr[k * n + q], bi = vi[k * n + q];
        vr[k * n + p] = c * ar - s * br; vi[k * n + p] = c * ai - s * bi;
        vr[k * n + q] = s * ar + c * br; vi[k * n + q] = s * ai + c * bi;
      }
    }}
  }
  const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => re[a * n + a] - re[b * n + b]);
  const vectors = cmat(n);
  order.forEach((col, i) => { for (let k = 0; k < n; k++) { vectors.re[k * n + i] = vr[k * n + col]; vectors.im[k * n + i] = vi[k * n + col]; } });
  return { values: Float64Array.from(order, i => re[i * n + i]), vectors };
}

/** U(t) = V diag(e^{−iE_k t/ħ}) V† from an eigHermitian result. */
export function propagator({ values, vectors }, t, hbar = 1) {
  const n = values.length, D = cmat(n);
  for (let k = 0; k < n; k++) { const ph = -values[k] * t / hbar; D.re[k * n + k] = Math.cos(ph); D.im[k * n + k] = Math.sin(ph); }
  return cmul(cmul(vectors, D), dagger(vectors));
}
export const evolve = (rho, U) => cmul(cmul(U, rho), dagger(U));

/** Thermal (Gibbs) state e^{−H/kT}/Z. */
export function thermalState(H, kT) {
  const { values, vectors } = eigHermitian(H), n = H.n, D = cmat(n), e0 = values[0];
  let Z = 0; for (const e of values) Z += Math.exp(-(e - e0) / kT);
  values.forEach((e, k) => { D.re[k * n + k] = Math.exp(-(e - e0) / kT) / Z; });
  return cmul(cmul(vectors, D), dagger(vectors));
}

/**
 * Open-system dynamics, Lindblad form: dρ/dt = −i[H, ρ]/ħ + Σ_k γ_k (L_k ρ L_k† − ½{L_k†L_k, ρ}).
 * jumps: [{ L (cmat), rate γ }] (dephasing: L = |i⟩⟨i|; decay i → j: L = |j⟩⟨i|; a trap/sink is a
 * decay into an extra level). Trace-preserving and completely positive. Adaptive Dormand–Prince on
 * the 2n² real entries. Returns { t, rho (cmat, current), step(), advanceTo(t), population(i) }.
 */
export function lindblad({ H, jumps = [], rho0, hbar = 1, t0 = 0, rtol = 1e-10, atol = 1e-12, h = 1e-3 }) {
  const n = H.n, N = n * n;
  const Ls = jumps.map(({ L, rate }) => { if (!(rate >= 0)) throw new Error('lindblad: rates must be ≥ 0'); return { L, Ld: dagger(L), LdL: cmul(dagger(L), L), rate }; });
  const unpack = y => ({ n, re: y.subarray(0, N), im: y.subarray(N) });
  const f = (_t, y) => {
    const rho = unpack(y), out = new Float64Array(2 * N);
    const comm = cadd(cmul(H, rho), cmul(rho, H), -1);          // [H, ρ]
    for (let i = 0; i < N; i++) { out[i] = comm.im[i] / hbar; out[N + i] = -comm.re[i] / hbar; }   // −i[H, ρ]/ħ
    for (const { L, Ld, LdL, rate } of Ls) {
      if (!rate) continue;
      const a = cmul(cmul(L, rho), Ld), b = cmul(LdL, rho), c = cmul(rho, LdL);
      for (let i = 0; i < N; i++) { out[i] += rate * (a.re[i] - 0.5 * (b.re[i] + c.re[i])); out[N + i] += rate * (a.im[i] - 0.5 * (b.im[i] + c.im[i])); }
    }
    return out;
  };
  const y0 = new Float64Array(2 * N); y0.set(rho0.re); y0.set(rho0.im, N);
  const s = odeStepper({ f, y0, t0, method: 'dopri5', h, rtol, atol });
  const api = {
    get t() { return s.t; },
    get rho() { return { n, re: Float64Array.from(s.y.subarray(0, N)), im: Float64Array.from(s.y.subarray(N)) }; },
    population: i => s.y[i * n + i],
    step: () => s.step(),
    advanceTo: (t, opts) => s.advanceTo(t, Infinity, opts),
    stepper: s,
  };
  return api;
}

/** Trace, Hermiticity error and smallest eigenvalue of a density matrix (all should be 1, ~0, ≥ −round-off). */
export function checkDensity(rho) {
  const n = rho.n; let herm = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) herm = Math.max(herm, Math.hypot(rho.re[i * n + j] - rho.re[j * n + i], rho.im[i * n + j] + rho.im[j * n + i]));
  const sym = cadd(rho, dagger(rho)); const half = cscale(sym, 0.5);
  return { trace: trace(rho).re, hermiticity: herm, minEigenvalue: eigHermitian(half).values[0] };
}
