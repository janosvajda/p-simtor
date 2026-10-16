// Geometric (Clifford) algebra Cl(p, q): multivectors, products, rotors. Plain ES module.
//
//   cliffordAlgebra(p, q)   an algebra object: products, involutions, grades, exp, inverse, dual
//   CL3, STA                Cl(3,0) (Euclidean 3D, basis e1 e2 e3) and the spacetime algebra
//                           Cl(1,3) (basis e0 e1 e2 e3, signature + − − −, c = 1)
//   rotor3, rotor3Between, rotateVector3, rotorToMatrix3, cross3
//                           3D rotations as rotors R = exp(−I n̂ θ/2), v ↦ R v R̃
//   boostRotor, fourVector, fourVelocity, spatialVelocity
//                           Lorentz boosts as rotors R = exp(−e0∧n̂ φ/2) in STA
//   fieldBivector, fieldComponents, fieldInvariants
//                           Maxwell field F = E + I B (σ_k = e_k e0, I = e0e1e2e3); F² = E²−B² + 2I E·B
//
// Representation: a multivector is a Float64Array of 2ⁿ coefficients indexed by blade bitmask
// (bit i set ⇔ basis vector i is a factor, factors in increasing order: index 0b101 = e1e3
// in CL3). Products are exact sums of coefficient products (rounding only). Conventions follow
// Doran & Lasenby, "Geometric Algebra for Physicists" (CUP 2003): inner = left contraction,
// reverse ~, rotors normalised to R R̃ = 1, active transformations x ↦ R x R̃.

import { solveLinear } from './linalg.mjs';

const popcount = x => { let c = 0, v = x; while (v) { c += v & 1; v >>>= 1; } return c; };
// sign from reordering blade a·b into canonical order (count transpositions)
const reorderSign = (a, b) => {
  let s = 0, v = a >>> 1;
  while (v) { s += popcount(v & b); v >>>= 1; }
  return s & 1 ? -1 : 1;
};

/**
 * The Clifford algebra Cl(p, q) on n = p + q ≤ 8 basis vectors: the first p square to +1, the
 * next q to −1. firstIndex labels basis vectors in names (default 1: e1…en). Returns an object
 * whose methods take and return Float64Array multivectors (inputs are never modified):
 *   zero(), scalar(s), vector(components), blade(name | mask, c), basis(i), pseudoscalar()
 *   add, sub, scale(a, s), neg, mul (geometric), wedge (outer), inner (left contraction a⌋b),
 *   scalarProduct (⟨ab⟩₀), commutator (½(ab − ba))
 *   grade(a, k), grades(a) (non-zero grades), reverse, involute, conjugate, dual (a I⁻¹)
 *   norm2 (⟨a ã⟩₀; may be negative in mixed signature), inverse, exp, sandwich(R, x) = R x R̃
 *   equals(a, b, tol), toString(a)
 * plus n, p, q, size, names, gradeOf (per blade index).
 */
export function cliffordAlgebra(p, q = 0, { firstIndex = 1 } = {}) {
  if (!(Number.isSafeInteger(p) && Number.isSafeInteger(q) && p >= 0 && q >= 0 && p + q >= 1 && p + q <= 8)) {
    throw new Error('cliffordAlgebra: need integers p, q ≥ 0 with 1 ≤ p + q ≤ 8');
  }
  const n = p + q, size = 1 << n;
  const square = Array.from({ length: n }, (_, i) => (i < p ? 1 : -1));
  const gradeOf = Uint8Array.from({ length: size }, (_, i) => popcount(i));
  const sign = new Int8Array(size * size);
  for (let a = 0; a < size; a++) {
    for (let b = 0; b < size; b++) {
      let s = reorderSign(a, b);
      for (let i = 0, common = a & b; common; i++, common >>>= 1) if (common & 1) s *= square[i];
      sign[a * size + b] = s;
    }
  }
  const names = Array.from({ length: size }, (_, m) => {
    if (m === 0) return '1';
    let s = 'e';
    for (let i = 0; i < n; i++) if (m & (1 << i)) s += String(i + firstIndex);
    return s;
  });
  const indexOfName = new Map(names.map((nm, i) => [nm, i]));
  const mv = (a, who) => {
    if (!(a instanceof Float64Array) || a.length !== size) throw new Error(`${who}: expected a multivector of this algebra (Float64Array(${size}))`);
    return a;
  };
  const zero = () => new Float64Array(size);
  const product = (a, b, keep, who) => {
    mv(a, who); mv(b, who);
    const out = zero();
    for (let i = 0; i < size; i++) {
      const ai = a[i];
      if (ai === 0) continue;
      for (let j = 0; j < size; j++) {
        const bj = b[j];
        if (bj === 0 || !keep(i, j)) continue;
        out[i ^ j] += sign[i * size + j] * ai * bj;
      }
    }
    return out;
  };
  const involution = (a, f, who) => { mv(a, who); return Float64Array.from(a, (v, i) => (f(gradeOf[i]) ? -v : v)); };
  const A = {
    n, p, q, size, names, gradeOf,
    zero,
    scalar: s => { const out = zero(); out[0] = s; return out; },
    vector: comps => {
      if (!comps || comps.length !== n) throw new Error(`vector: need ${n} components`);
      const out = zero();
      for (let i = 0; i < n; i++) out[1 << i] = comps[i];
      return out;
    },
    blade: (id, c = 1) => {
      const m = typeof id === 'string' ? indexOfName.get(id) : id;
      if (!(Number.isSafeInteger(m) && m >= 0 && m < size)) throw new Error(`blade: unknown blade ${id} (names: ${names.join(' ')})`);
      const out = zero(); out[m] = c; return out;
    },
    basis: i => A.blade(1 << (i - firstIndex)),
    pseudoscalar: () => A.blade(size - 1),
    add: (a, b) => { mv(a, 'add'); mv(b, 'add'); return Float64Array.from(a, (v, i) => v + b[i]); },
    sub: (a, b) => { mv(a, 'sub'); mv(b, 'sub'); return Float64Array.from(a, (v, i) => v - b[i]); },
    scale: (a, s) => { mv(a, 'scale'); return Float64Array.from(a, v => v * s); },
    neg: a => A.scale(a, -1),
    /** Geometric product ab. */
    mul: (a, b) => product(a, b, () => true, 'mul'),
    /** Outer product a∧b (grade-raising part: blades with no common factor). */
    wedge: (a, b) => product(a, b, (i, j) => (i & j) === 0, 'wedge'),
    /** Left contraction a⌋b: the grade (s − r) part of (r-blade)(s-blade); zero when r > s. */
    inner: (a, b) => product(a, b, (i, j) => (i & j) === i, 'inner'),
    scalarProduct: (a, b) => {
      mv(a, 'scalarProduct'); mv(b, 'scalarProduct');
      let s = 0;
      for (let i = 0; i < size; i++) if (a[i] !== 0 && b[i] !== 0) s += sign[i * size + i] * a[i] * b[i];
      return s;
    },
    commutator: (a, b) => A.scale(A.sub(A.mul(a, b), A.mul(b, a)), 0.5),
    grade: (a, k) => { mv(a, 'grade'); return Float64Array.from(a, (v, i) => (gradeOf[i] === k ? v : 0)); },
    grades: (a, tol = 0) => { mv(a, 'grades'); return [...new Set(Array.from(a).flatMap((v, i) => (Math.abs(v) > tol ? [gradeOf[i]] : [])))].sort((x, y) => x - y); },
    /** Reverse ã: reverses the order of vector factors; sign (−1)^{k(k−1)/2} on grade k. */
    reverse: a => involution(a, k => (k * (k - 1) / 2) % 2 === 1, 'reverse'),
    /** Grade involution â: sign (−1)^k. */
    involute: a => involution(a, k => k % 2 === 1, 'involute'),
    /** Clifford conjugate ā = reverse of involution: sign (−1)^{k(k+1)/2}. */
    conjugate: a => involution(a, k => (k * (k + 1) / 2) % 2 === 1, 'conjugate'),
    dual: a => A.mul(a, A.inverse(A.pseudoscalar())),
    norm2: a => A.scalarProduct(a, A.reverse(a)),
    /**
     * Inverse: ã/⟨aã⟩₀ when a ã is a scalar (vectors, blades, versors: exact); otherwise by
     * solving the 2ⁿ×2ⁿ linear system of left multiplication. Throws if a is not invertible.
     */
    inverse: a => {
      mv(a, 'inverse');
      let scale = 0;
      for (const v of a) scale = Math.max(scale, Math.abs(v));
      if (!Number.isFinite(scale)) throw new Error('inverse: coefficients must be finite');
      if (scale !== 0 && scale !== 1) return Float64Array.from(A.inverse(Float64Array.from(a, v => v / scale)), v => v / scale);
      const ar = A.reverse(a), aa = A.mul(a, ar), s = aa[0];
      let rest = 0;
      for (let i = 1; i < size; i++) rest = Math.max(rest, Math.abs(aa[i]));
      // The neglected part is divided by s as well. Compare against |s|, not |a|²:
      // near a null multivector those scales can differ by many orders of magnitude.
      if (s !== 0 && rest <= 1e-14 * Math.abs(s)) return Float64Array.from(ar, v => v / s);
      const L = Array.from({ length: size }, (_r, row) => Float64Array.from({ length: size }, (_c, col) => {
        // (a · e_col) coefficient on blade row
        const j = row ^ col;
        return sign[j * size + col] * a[j];
      }));
      let x;
      try { x = solveLinear(L, Float64Array.from({ length: size }, (_, i) => (i === 0 ? 1 : 0))); }
      catch { throw new Error('inverse: multivector is not invertible'); }
      return Float64Array.from(x);
    },
    /**
     * Exponential e^a. If a² is a scalar s (any simple bivector, vector or null element):
     * closed form cos/cosh/1 + a — exact to rounding. Otherwise (e.g. a non-simple bivector
     * such as a general F) Taylor series with scaling and squaring, ~1e-15 relative.
     */
    exp: a => {
      mv(a, 'exp');
      if (a[0] !== 0) { const rest = a.slice(); rest[0] = 0; return A.scale(A.exp(rest), Math.exp(a[0])); }   // scalars commute
      const a2 = A.mul(a, a);
      let rest = 0, mag = 0;
      for (let i = 1; i < size; i++) rest = Math.max(rest, Math.abs(a2[i]));
      for (const v of a2) mag = Math.max(mag, Math.abs(v));
      if (rest <= 1e-15 * Math.max(mag, 1e-300)) {
        const s = a2[0], out = A.scale(a, 1);
        if (s < 0) { const t = Math.sqrt(-s); for (let i = 0; i < size; i++) out[i] *= Math.sin(t) / t; out[0] += Math.cos(t); }
        else if (s > 0) { const t = Math.sqrt(s); for (let i = 0; i < size; i++) out[i] *= Math.sinh(t) / t; out[0] += Math.cosh(t); }
        else out[0] += 1;
        return out;
      }
      let normA = 0;
      for (const v of a) normA += Math.abs(v);
      const k = Math.max(0, Math.ceil(Math.log2(normA)) + 1), x = A.scale(a, 2 ** -k);
      let term = A.scalar(1), sum = A.scalar(1);
      for (let m = 1; m < 40; m++) {
        term = A.scale(A.mul(term, x), 1 / m);
        sum = A.add(sum, term);
        if (term.every(v => Math.abs(v) <= 1e-17 * Math.max(1, Math.abs(sum[0])))) break;
      }
      for (let i = 0; i < k; i++) sum = A.mul(sum, sum);
      return sum;
    },
    /** R x R̃ (for a rotor R R̃ = 1 this is the rotation/boost of x; grades are preserved). */
    sandwich: (R, x) => A.mul(A.mul(R, x), A.reverse(R)),
    equals: (a, b, tol = 1e-12) => { mv(a, 'equals'); mv(b, 'equals'); return a.every((v, i) => Math.abs(v - b[i]) <= tol); },
    toString: a => {
      mv(a, 'toString');
      const terms = [];
      for (let i = 0; i < size; i++) if (a[i] !== 0) terms.push(i === 0 ? `${a[i]}` : `${a[i]}${names[i]}`);
      return terms.length ? terms.join(' + ').replace(/\+ -/g, '- ') : '0';
    },
  };
  return Object.freeze(A);
}

/** Euclidean 3D algebra Cl(3,0): e1, e2, e3 square to +1; pseudoscalar I = e123, I² = −1. */
export const CL3 = cliffordAlgebra(3, 0);
/** Spacetime algebra Cl(1,3): e0² = +1 (time), e1² = e2² = e3² = −1; I = e0123, I² = −1. */
export const STA = cliffordAlgebra(1, 3, { firstIndex: 0 });

const vec3 = (v, who) => {
  if (!v || v.length !== 3 || !Array.from(v).every(Number.isFinite)) throw new Error(`${who}: need a finite 3-vector`);
  return Array.from(v);
};
const unit3 = (v, who) => {
  const u = vec3(v, who), r = Math.hypot(...u);
  if (!(r > 0)) throw new Error(`${who}: direction must be non-zero`);
  return u.map(c => c / r);
};

/**
 * 3D rotor for a right-handed rotation by angle (radians) about axis (any non-zero length):
 * R = exp(−I n̂ θ/2) = cos(θ/2) − I n̂ sin(θ/2) in CL3. Rotors compose by R₂R₁ (first R₁);
 * R and −R give the same rotation (double cover). Equivalent to a unit quaternion.
 */
export function rotor3(axis, angle) {
  if (!Number.isFinite(angle)) throw new Error('rotor3: angle must be finite');
  const nhat = unit3(axis, 'rotor3'), B = CL3.mul(CL3.pseudoscalar(), CL3.vector(nhat));
  return CL3.exp(CL3.scale(B, -angle / 2));
}

/** Rotor taking direction a to direction b by the smallest rotation (throws if antiparallel). */
export function rotor3Between(a, b) {
  const u = unit3(a, 'rotor3Between'), v = unit3(b, 'rotor3Between');
  const R = CL3.add(CL3.scalar(1), CL3.mul(CL3.vector(v), CL3.vector(u)));
  const n2 = CL3.norm2(R);
  if (!(n2 > 1e-24)) throw new Error('rotor3Between: antiparallel directions have no unique rotation');
  return CL3.scale(R, 1 / Math.sqrt(n2));
}

/** Rotate a 3-vector by a CL3 rotor: R v R̃, returned as [x, y, z]. */
export function rotateVector3(R, v) {
  const out = CL3.sandwich(R, CL3.vector(vec3(v, 'rotateVector3')));
  return [out[1], out[2], out[4]];
}

/** The 3×3 rotation matrix (rows) of a CL3 rotor: column j is R e_j R̃. */
export function rotorToMatrix3(R) {
  const cols = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(e => rotateVector3(R, e));
  return [0, 1, 2].map(i => cols.map(c => c[i]));
}

/** Cross product from the outer product: a × b = −I (a ∧ b). */
export function cross3(a, b) {
  const w = CL3.wedge(CL3.vector(vec3(a, 'cross3')), CL3.vector(vec3(b, 'cross3')));
  const c = CL3.neg(CL3.mul(CL3.pseudoscalar(), w));
  return [c[1], c[2], c[4]];
}

/** STA vector t e0 + x e1 + y e2 + z e3 from [t, x, y, z]. Its square is t² − x² − y² − z². */
export function fourVector(x) {
  if (!x || x.length !== 4) throw new Error('fourVector: need [t, x, y, z]');
  return STA.vector(x);
}

/** Four-velocity γ(e0 + v) of 3-velocity v (|v| < 1, units of c). */
export function fourVelocity(v) {
  const u = vec3(v, 'fourVelocity'), s2 = u[0] ** 2 + u[1] ** 2 + u[2] ** 2;
  if (!(s2 < 1)) throw new Error('fourVelocity: need |v| < 1 (units of c)');
  const g = 1 / Math.sqrt(1 - s2);
  return fourVector([g, g * u[0], g * u[1], g * u[2]]);
}

/** 3-velocity [x/t, y/t, z/t] of a timelike STA vector. */
export function spatialVelocity(u) {
  if (!(u[1] !== 0)) throw new Error('spatialVelocity: need a vector with a non-zero time component');
  return [u[2] / u[1], u[4] / u[1], u[8] / u[1]];
}

/**
 * Boost rotor R = exp(−(e0 ∧ n̂) φ/2), φ = atanh|v|: R e0 R̃ = γ(e0 + v), i.e. it actively
 * boosts a particle at rest to 3-velocity v (|v| < 1, units of c). Boosts along one line
 * compose by adding rapidities; non-collinear boosts compose to a boost times a (Thomas–Wigner)
 * rotation.
 */
export function boostRotor(v) {
  const u = vec3(v, 'boostRotor'), speed = Math.hypot(...u);
  if (!(speed < 1)) throw new Error('boostRotor: need |v| < 1 (units of c)');
  if (speed === 0) return STA.scalar(1);
  const phi = Math.atanh(speed), K = STA.wedge(STA.basis(0), STA.vector([0, u[0] / speed, u[1] / speed, u[2] / speed]));
  return STA.exp(STA.scale(K, -phi / 2));
}

const SIGMA = [1, 2, 3].map(k => STA.mul(STA.basis(k), STA.basis(0)));      // σ_k = e_k e0 (relative vectors)
const ISIGMA = SIGMA.map(s => STA.mul(STA.pseudoscalar(), s));

/** Electromagnetic field bivector F = E + I B with E = Σ E_k σ_k, B = Σ B_k σ_k (Gaussian/natural units, c = 1). */
export function fieldBivector(E, B) {
  const e = vec3(E, 'fieldBivector'), b = vec3(B, 'fieldBivector');
  let F = STA.zero();
  for (let k = 0; k < 3; k++) F = STA.add(F, STA.add(STA.scale(SIGMA[k], e[k]), STA.scale(ISIGMA[k], b[k])));
  return F;
}

/** { E, B } relative to the e0 observer of a field bivector (inverse of fieldBivector). */
export function fieldComponents(F) {
  const coef = X => STA.scalarProduct(F, STA.reverse(X)) / STA.scalarProduct(X, STA.reverse(X));
  return { E: SIGMA.map(coef), B: ISIGMA.map(coef) };
}

/**
 * Lorentz invariants from F² = (E² − B²) + 2 I (E·B): returns { EsqMinusBsq, EdotB } read off the
 * scalar and pseudoscalar parts of F². Unchanged by any boost or rotation F ↦ R F R̃.
 */
export function fieldInvariants(F) {
  const F2 = STA.mul(F, F);
  return { EsqMinusBsq: F2[0], EdotB: F2[STA.size - 1] / 2 };
}
