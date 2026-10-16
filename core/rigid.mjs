// Rigid bodies: inertia tensors, principal axes, torque-free motion (exact Jacobi-elliptic solution
// and a structure-preserving integrator) and the heavy symmetric (Lagrange) top.
//
// Units: any consistent set (SI: kg, m, s; moments in kg m², angular momentum in kg m²/s).
// Body frame = principal axes of the body. Orientation q is a unit quaternion (core/quaternion.mjs)
// mapping body to space coordinates; the space z axis points up and gravity is along −z.
// Matrices are arrays of three Float64Array rows, as in core/linalg.mjs.

import { eigenSymmetric } from './linalg.mjs';
import { ellipticK, brent } from './numerics.mjs';
import { jacobiElliptic } from './special.mjs';
import { qmul, qnormalize, qExp, qRotate, qToMatrix, qFromMatrix, qFromAxisAngle, qconj } from './quaternion.mjs';

const finite = (v, n, who) => {
  if (!v || v.length !== n || !Array.from(v).every(Number.isFinite)) throw new Error(`${who}: expected ${n} finite values`);
};
const positive = (x, name, who) => { if (!(Number.isFinite(x) && x > 0)) throw new Error(`${who}: ${name} must be positive and finite`); };
const diag = (a, b, c) => [Float64Array.of(a, 0, 0), Float64Array.of(0, b, 0), Float64Array.of(0, 0, c)];
const cross = (a, b) => Float64Array.of(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
const checkMatrix = (I, who) => { if (!I || I.length !== 3) throw new Error(`${who}: expected a 3×3 matrix`); for (const r of I) finite(r, 3, who); };

// ---------------------------------------------------------------- inertia tensors (about the centre of mass)

/** Solid cuboid of mass m and full side lengths [a, b, c] along x, y, z: I = m/12 diag(b²+c², a²+c², a²+b²). */
export function inertiaBox(m, [a, b, c]) {
  positive(m, 'mass', 'inertiaBox'); positive(a, 'side', 'inertiaBox'); positive(b, 'side', 'inertiaBox'); positive(c, 'side', 'inertiaBox');
  return diag(m * (b * b + c * c) / 12, m * (a * a + c * c) / 12, m * (a * a + b * b) / 12);
}

/** Sphere of mass m, radius r: solid 2/5 m r², thin shell ({ shell: true }) 2/3 m r². */
export function inertiaSphere(m, r, { shell = false } = {}) {
  positive(m, 'mass', 'inertiaSphere'); positive(r, 'radius', 'inertiaSphere');
  const I = (shell ? 2 / 3 : 2 / 5) * m * r * r;
  return diag(I, I, I);
}

/**
 * Cylinder (axis z) of mass m, outer radius r, length h, optional inner radius (thick tube, uniform
 * density): I_z = ½ m (r² + r_in²), I_x = I_y = m (3(r² + r_in²) + h²)/12.
 */
export function inertiaCylinder(m, r, h, { innerRadius = 0 } = {}) {
  positive(m, 'mass', 'inertiaCylinder'); positive(r, 'radius', 'inertiaCylinder'); positive(h, 'length', 'inertiaCylinder');
  if (!(innerRadius >= 0 && innerRadius < r)) throw new Error('inertiaCylinder: innerRadius must be in [0, r)');
  const s = r * r + innerRadius * innerRadius;
  return diag(m * (3 * s + h * h) / 12, m * (3 * s + h * h) / 12, m * s / 2);
}

/**
 * Solid right circular cone (axis z) of mass m, base radius r, height h, about its centre of mass,
 * which lies on the axis h/4 above the base (3h/4 below the apex):
 * I_z = 3/10 m r², I_x = I_y = 3/20 m r² + 3/80 m h². Returns { inertia, comFromBase: h/4 }.
 */
export function inertiaCone(m, r, h) {
  positive(m, 'mass', 'inertiaCone'); positive(r, 'radius', 'inertiaCone'); positive(h, 'height', 'inertiaCone');
  const t = 3 / 20 * m * r * r + 3 / 80 * m * h * h;
  return { inertia: diag(t, t, 3 / 10 * m * r * r), comFromBase: h / 4 };
}

/**
 * Parallel-axis (Huygens–Steiner) theorem: inertia about a point P from the inertia about the centre
 * of mass, I_P = I_cm + m (|d|² 1 − d dᵀ), d = position of the centre of mass relative to P.
 */
export function parallelAxis(Icm, m, d) {
  checkMatrix(Icm, 'parallelAxis'); finite(d, 3, 'parallelAxis');
  if (!Number.isFinite(m)) throw new Error('parallelAxis: mass must be finite');
  const d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  return [0, 1, 2].map(i => Float64Array.from([0, 1, 2], j => Icm[i][j] + m * ((i === j ? d2 : 0) - d[i] * d[j])));
}

/**
 * Composite body from parts [{ mass, inertia (about the part's own centre of mass, in the part frame),
 * position (part centre of mass in the body frame), orientation? (unit quaternion part → body) }].
 * A part may have negative mass to cut a cavity (its inertia then negative too); the total mass
 * must be positive. Returns { mass, com, inertia } with the inertia about the composite centre of mass.
 */
export function compositeBody(parts) {
  if (!Array.isArray(parts) || !parts.length) throw new Error('compositeBody: need a non-empty array of parts');
  let M = 0;
  const com = new Float64Array(3);
  for (const p of parts) {
    if (!Number.isFinite(p.mass)) throw new Error('compositeBody: each part needs a finite mass');
    finite(p.position, 3, 'compositeBody'); checkMatrix(p.inertia, 'compositeBody');
    M += p.mass;
    for (let k = 0; k < 3; k++) com[k] += p.mass * p.position[k];
  }
  if (!(M > 0)) throw new Error('compositeBody: total mass must be positive');
  for (let k = 0; k < 3; k++) com[k] /= M;
  const I = [new Float64Array(3), new Float64Array(3), new Float64Array(3)];
  for (const p of parts) {
    let Ip = p.inertia;
    if (p.orientation) {
      const R = qToMatrix(p.orientation);
      Ip = [0, 1, 2].map(i => Float64Array.from([0, 1, 2], j => { let s = 0; for (let k = 0; k < 3; k++) for (let l = 0; l < 3; l++) s += R[i][k] * p.inertia[k][l] * R[j][l]; return s; }));
    }
    const shifted = parallelAxis(Ip, p.mass, [0, 1, 2].map(k => p.position[k] - com[k]));
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) I[i][j] += shifted[i][j];
  }
  return { mass: M, com, inertia: I };
}

/**
 * Principal moments (ascending) and axes of a symmetric inertia tensor via eigenSymmetric.
 * axes[k] is the unit principal axis k (in the input frame), made right-handed (axes[2] = axes[0] × axes[1]);
 * rotation has the axes as columns (principal → input frame) and quaternion is the same rotation.
 * Throws if the tensor is not symmetric or not positive definite.
 */
export function principalAxes(I) {
  checkMatrix(I, 'principalAxes');
  const { values, vectors } = eigenSymmetric(I);
  if (!(values[0] > 0)) throw new Error('principalAxes: inertia tensor must be positive definite');
  const axes = [Float64Array.from(vectors[0]), Float64Array.from(vectors[1])];
  axes.push(cross(axes[0], axes[1]));
  const rotation = [0, 1, 2].map(i => Float64Array.from(axes, a => a[i]));
  return { moments: Float64Array.from(values), axes, rotation, quaternion: qFromMatrix(rotation) };
}

// ---------------------------------------------------------------- elliptic functions

/** Carlson's symmetric integral R_F(x, y, z) by duplication (Carlson, Numer. Algorithms 10, 13 (1995)). */
function carlsonRF(x0, y0, z0) {
  let x = x0, y = y0, z = z0;
  for (let i = 0; i < 200; i++) {
    const A = (x + y + z) / 3, r = Math.max(Math.abs(A - x), Math.abs(A - y), Math.abs(A - z)) / A;
    if (r < 1e-3) break;
    const sx = Math.sqrt(x), sy = Math.sqrt(y), sz = Math.sqrt(z), lam = sx * sy + sy * sz + sz * sx;
    x = (x + lam) / 4; y = (y + lam) / 4; z = (z + lam) / 4;
  }
  const A = (x + y + z) / 3, X = 1 - x / A, Y = 1 - y / A, Z = -X - Y;
  const E2 = X * Y - Z * Z, E3 = X * Y * Z;
  return (1 - E2 / 10 + E3 / 14 + E2 * E2 / 24 - 3 * E2 * E3 / 44) / Math.sqrt(A);
}

/**
 * Incomplete elliptic integral of the first kind F(φ | m) = ∫₀^φ dθ / √(1 − m sin²θ), parameter
 * m = k² ∈ [0, 1], any real φ (quasi-periodic: F(φ + jπ) = F(φ) + 2jK). Carlson R_F, ~1e-15 relative.
 * For m = 1 only |φ| < π/2 is finite. It inverts the amplitude: F(am(u | m) | m) = u.
 */
export function ellipticF(phi, m) {
  if (!Number.isFinite(phi) || !(m >= 0 && m <= 1)) throw new Error('ellipticF: need finite φ and 0 ≤ m ≤ 1');
  const j = Math.round(phi / Math.PI), r = phi - j * Math.PI, s = Math.sin(r), c = Math.cos(r);
  const base = s * carlsonRF(c * c, 1 - m * s * s, 1);
  if (j === 0) return base;
  if (m === 1) return Math.sign(phi) * Infinity;
  return base + 2 * j * ellipticK(m);
}

// Jacobi elliptic functions: the general implementation in special.mjs (re-exported here)
export { jacobiElliptic };

// ---------------------------------------------------------------- torque-free motion: exact solution

const sortedMoments = (inertia, who) => {
  finite(inertia, 3, who);
  const [I1, I2, I3] = inertia;
  if (!(I1 > 0 && I1 <= I2 && I2 <= I3)) throw new Error(`${who}: principal moments must satisfy 0 < I1 ≤ I2 ≤ I3`);
  if (I3 > I1 + I2 * (1 + 1e-12)) throw new Error(`${who}: moments violate the triangle inequality I3 ≤ I1 + I2`);
  return [I1, I2, I3];
};

/**
 * Exact torque-free motion of the body angular momentum L(t) for principal moments I1 ≤ I2 ≤ I3
 * (Jacobi elliptic functions; Landau & Lifshitz, Mechanics §37). With 2E = Σ L_i²/I_i and
 * L² = |L|²: for L² > 2E I2 the polhode circles axis 3 (Ω1 ∝ cn, Ω2 ∝ sn, Ω3 ∝ dn), for L² < 2E I2
 * it circles axis 1 (roles of 1 and 3 exchanged). Any initial L0 (body components) is matched by
 * a phase shift and the π-rotation symmetries of Euler's equations.
 * Returns { at(t) → body L, period (of L(t): 4K(m)/λ), m (= k²), lambda (rate of τ = λt), axis
 * (0 or 2: the axis the motion circles), energy, L2 }. Throws on the separatrix L² = 2E I2
 * (relative 1e-13), where the period is infinite, for triaxial bodies. Symmetric tops use
 * their elementary rotation, including stationary states; spherical bodies give L constant.
 */
export function torqueFreeExact({ inertia, L0 }) {
  const I = sortedMoments(inertia, 'torqueFreeExact');
  finite(L0, 3, 'torqueFreeExact');
  const L2 = L0[0] ** 2 + L0[1] ** 2 + L0[2] ** 2, E2 = L0[0] ** 2 / I[0] + L0[1] ** 2 / I[1] + L0[2] ** 2 / I[2];
  const constant = { at: () => Float64Array.from(L0), period: Infinity, m: 0, lambda: 0, axis: 2, energy: E2 / 2, L2 };
  if (I[0] === I[2] || L2 === 0) return constant;
  if (I[0] === I[1] || I[1] === I[2]) {
    const axis = I[0] === I[1] ? 2 : 0, j = (axis + 1) % 3, k = (axis + 2) % 3;
    const rate = (1 / I[j] - 1 / I[axis]) * L0[axis];
    const stationary = rate === 0 || (L0[j] === 0 && L0[k] === 0);
    return { ...constant, axis, lambda: Math.abs(rate), period: stationary ? Infinity : 2 * Math.PI / Math.abs(rate),
      at(t) {
        const out = Float64Array.from(L0), c = Math.cos(rate * t), s = Math.sin(rate * t);
        out[j] = L0[j] * c - L0[k] * s; out[k] = L0[j] * s + L0[k] * c;
        return out;
      } };
  }
  const D = L2 - E2 * I[1];
  if (Math.abs(D) <= 1e-13 * L2) throw new Error('torqueFreeExact: on the separatrix L² = 2E·I2 (infinite period)');
  const [a, b, c, eps] = D > 0 ? [0, 1, 2, 1] : [2, 1, 0, -1];
  const Ia = I[a], Ib = I[b], Ic = I[c];
  const A = Math.sqrt(Math.max(0, (L2 - E2 * Ic) / (Ia * (Ia - Ic))));
  const B = Ib === Ic ? 0 : Math.sqrt(Math.max(0, (L2 - E2 * Ic) / (Ib * (Ib - Ic))));
  const C = Math.sqrt(Math.max(0, (L2 - E2 * Ia) / (Ic * (Ic - Ia))));
  const lambda = Math.sqrt((Ic - Ib) * (Ic - Ia) / (Ia * Ib)) * C;
  const m = A > 0 && B > 0 ? Math.min(1, Math.max(0, eps * (Ib - Ia) * A * B / (Ic * C * lambda))) : 0;
  const w0 = [L0[0] / I[0], L0[1] / I[1], L0[2] / I[2]];
  const flip = w0[c] < 0 ? -1 : 1;                       // π rotation about axis b: Ω_a, Ω_c → −Ω_a, −Ω_c
  const phi0 = Math.atan2(B > 0 ? w0[b] / B : 0, A > 0 ? flip * w0[a] / A : 1);
  const tau0 = ellipticF(phi0, m);
  return {
    period: 4 * ellipticK(m) / lambda, m, lambda, axis: c, energy: E2 / 2, L2,
    at(t) {
      const { sn, cn, dn } = jacobiElliptic(tau0 + lambda * t, m);
      const w = new Float64Array(3);
      w[a] = flip * A * cn; w[b] = B * sn; w[c] = flip * C * dn;
      return Float64Array.from(w, (v, i) => v * I[i]);
    },
  };
}

// ---------------------------------------------------------------- structure-preserving integrator

const CBRT2 = Math.cbrt(2);
const YOSHIDA4 = [1 / (2 - CBRT2), -CBRT2 / (2 - CBRT2), 1 / (2 - CBRT2)];

/** Split H = Σ L_i²/(2I_i) = |L|²/(2I_r) + Σ_{i≠r} c_i L_i², each piece exactly solvable. */
function kineticPieces(I) {
  // reference axis: one of a repeated pair when possible, so a symmetric top has a single, commuting extra piece
  const r = I[0] === I[1] ? 0 : I[1] === I[2] ? 1 : I[0] === I[2] ? 0 : 1;
  const pieces = [{ axis: -1, c: 1 / (2 * I[r]) }];
  for (let i = 0; i < 3; i++) if (i !== r && I[i] !== I[r]) pieces.push({ axis: i, c: 1 / (2 * I[i]) - 1 / (2 * I[r]) });
  return pieces;
}

function rigidStepper(who, { inertia, L0, q0 = [1, 0, 0, 0], h, t0 = 0, method = 'strang', mgl = 0, chi = [0, 0, 1] }) {
  finite(inertia, 3, who); finite(L0, 3, who); finite(q0, 4, who); finite(chi, 3, who);
  if (!inertia.every(v => v > 0)) throw new Error(`${who}: principal moments must be positive`);
  if (!(Number.isFinite(h) && h !== 0)) throw new Error(`${who}: a non-zero finite step h is required`);
  if (!Number.isFinite(mgl)) throw new Error(`${who}: mgl must be finite`);
  const weights = method === 'strang' ? [1] : method === 'yoshida4' ? YOSHIDA4 : null;
  if (!weights) throw new Error(`${who}: unknown method "${method}" (use 'strang' or 'yoshida4')`);
  const I = Float64Array.from(inertia), pieces = kineticPieces(I), chiU = Float64Array.from(chi);
  const q = qnormalize(q0), L = Float64Array.from(L0), tq = new Float64Array(4), rot = new Float64Array(4);
  const up = Float64Array.of(0, 0, 1);
  const gammaInto = out => qRotate(qconj(q), up, out);       // space "up" in body coordinates
  const G = new Float64Array(3);

  function kineticPiece(p, dt) {
    if (p.axis < 0) {                                         // |L|²/(2I): rotate the body about L, L fixed
      const f = p.c * dt;
      q.set(qmul(q, qExp([0, f * L[0], f * L[1], f * L[2]]), tq));
      return;
    }
    const i = p.axis, j = (i + 1) % 3, k = (i + 2) % 3, th = 2 * p.c * L[i] * dt;
    const cs = Math.cos(th), sn = Math.sin(th), Lj = L[j], Lk = L[k];
    L[j] = cs * Lj + sn * Lk; L[k] = -sn * Lj + cs * Lk;       // body L rotates by −θ about e_i
    rot.fill(0); rot[0] = Math.cos(th / 2); rot[1 + i] = Math.sin(th / 2);
    q.set(qmul(q, rot, tq));                                  // body rotates by +θ about e_i
  }
  function kinetic(dt) {
    const n = pieces.length;
    for (let i = 0; i < n - 1; i++) kineticPiece(pieces[i], dt / 2);
    kineticPiece(pieces[n - 1], dt);
    for (let i = n - 2; i >= 0; i--) kineticPiece(pieces[i], dt / 2);
  }
  function potential(dt) {                                    // L += dt · mgl Γ × χ (orientation fixed)
    if (mgl === 0) return;
    gammaInto(G);
    const tq3 = cross(G, chiU);
    for (let k = 0; k < 3; k++) L[k] += dt * mgl * tq3[k];
  }

  const s = {
    t: t0, q, L, h, method, steps: 0,
    step() {
      for (const w of weights) { const dt = w * s.h; potential(dt / 2); kinetic(dt); potential(dt / 2); }
      q.set(qnormalize(q));
      s.t += s.h; s.steps++;
      if (!L.every(Number.isFinite)) throw new Error(`${who}: state became non-finite at t = ${s.t}`);
    },
    /**
     * Step until t reaches target: stops at the first step end at or past it, where a step end within
     * 1e-6 h of target counts as reaching it (accumulated round-off in t). Returns false if maxSteps ran out.
     */
    advanceTo(target, maxSteps = Infinity) {
      let i = 0;
      const slack = 1e-6 * Math.abs(s.h);
      while (s.h > 0 ? s.t < target - slack : s.t > target + slack) { if (i++ >= maxSteps) return false; s.step(); }
      return true;
    },
    /** Body angular velocity Ω = I⁻¹ L. */
    omega() { return Float64Array.from(L, (v, i) => v / I[i]); },
    /** Total energy Σ L_i²/(2I_i) + mgl Γ·χ. */
    energy() { gammaInto(G); return L[0] ** 2 / (2 * I[0]) + L[1] ** 2 / (2 * I[1]) + L[2] ** 2 / (2 * I[2]) + mgl * (G[0] * chiU[0] + G[1] * chiU[1] + G[2] * chiU[2]); },
    /** Angular momentum in space coordinates, R(q) L (conserved without gravity; its z component always). */
    spatialMomentum() { return qRotate(q, L); },
    /** Body → space rotation matrix. */
    rotation() { return qToMatrix(q); },
    /** Space "up" (+z) in body coordinates, Γ = Rᵀ e_z. */
    gamma() { return gammaInto(new Float64Array(3)); },
  };
  return s;
}

/**
 * Torque-free rigid body: Euler's equations L̇ = L × I⁻¹L with orientation q̇ = ½ q Ω, by
 * Lie–Poisson splitting (H = |L|²/2I_r + Σ c_i L_i², each piece an exact rotation; Strang
 * composition, McLachlan, Phys. Rev. Lett. 71, 3043 (1993); Dullweber, Leimkuhler & McLachlan,
 * J. Chem. Phys. 107, 5840 (1997)). Symplectic and time-reversible: |L| and the spatial angular
 * momentum are conserved to round-off, the energy error stays bounded, O(h²) ('strang') or O(h⁴)
 * ('yoshida4'). For a symmetric body (two equal moments) the pieces commute and the step is exact.
 * inertia = principal moments [I1, I2, I3] (any order), L0 = body angular momentum, q0 = orientation.
 * Returns { t, q, L, h, step(), advanceTo(t, maxSteps), omega(), energy(), spatialMomentum(), rotation(), gamma() }.
 * Set s.h = −h to run backwards (retraces to round-off).
 */
export function freeRigidBody({ inertia, L0, q0, h, t0, method }) {
  return rigidStepper('freeRigidBody', { inertia, L0, q0, h, t0, method });
}

/**
 * Heavy top with a fixed point: principal moments about the PIVOT, mgl = mass·g·(pivot → centre-of-mass
 * distance) (J), chi = unit body direction from pivot to centre of mass (default the body z axis).
 * Potential V = mgl (R χ)·e_z; body torque mgl Γ × χ, Γ = Rᵀe_z. Splitting V(h/2) T(h) V(h/2) with
 * T the free-body step above (Lie–Poisson integrator for the heavy top; Touma & Wisdom, AJ 107, 1189 (1994)).
 * Conserved to round-off: |Γ| = 1, the Casimir L·Γ (= vertical spatial angular momentum) and, for a
 * symmetric top with χ on the symmetry axis, L3. Energy error bounded (O(h²) or O(h⁴)).
 * Returns the same stepper as freeRigidBody (energy() includes the potential).
 */
export function heavyTop({ inertia, mgl, chi, L0, q0, h, t0, method }) {
  if (!(Number.isFinite(mgl))) throw new Error('heavyTop: mgl (mass·g·length) must be finite');
  if (chi && Math.abs(Math.hypot(...chi) - 1) > 1e-12) throw new Error('heavyTop: chi must be a unit vector');
  return rigidStepper('heavyTop', { inertia, mgl, chi, L0, q0, h, t0, method });
}

// ---------------------------------------------------------------- Lagrange (heavy symmetric) top

/**
 * Initial state of a symmetric top (I1 = I2 about the pivot, I3 about the symmetry axis = body z)
 * from z-x-z Euler angles: orientation q = Rz(φ) Rx(θ) Rz(ψ), θ = tilt of the axis from vertical.
 * Body angular velocity (Goldstein §4.9): ω = (φ̇ sinθ sinψ + θ̇ cosψ, φ̇ sinθ cosψ − θ̇ sinψ, φ̇ cosθ + ψ̇).
 * Returns { q, L (body), omega (body) } ready for heavyTop({ inertia: [I1, I1, I3], ... }).
 */
export function lagrangeTopInitial({ I1, I3, theta, phi = 0, psi = 0, thetaDot = 0, phiDot = 0, psiDot = 0 }) {
  positive(I1, 'I1', 'lagrangeTopInitial'); positive(I3, 'I3', 'lagrangeTopInitial');
  finite([theta, phi, psi, thetaDot, phiDot, psiDot], 6, 'lagrangeTopInitial');
  const q = qmul(qmul(qFromAxisAngle([0, 0, 1], phi), qFromAxisAngle([1, 0, 0], theta)), qFromAxisAngle([0, 0, 1], psi));
  const st = Math.sin(theta), omega = Float64Array.of(
    phiDot * st * Math.sin(psi) + thetaDot * Math.cos(psi),
    phiDot * st * Math.cos(psi) - thetaDot * Math.sin(psi),
    phiDot * Math.cos(theta) + psiDot);
  return { q, omega, L: Float64Array.of(I1 * omega[0], I1 * omega[1], I3 * omega[2]) };
}

/**
 * Exact nutation of the Lagrange top (Goldstein §5.7). With u = cos θ the motion obeys
 * u̇² = f(u) = (1 − u²)(α − βu) − (b − a u)², a = p_ψ/I1, b = p_φ/I1, α = 2(E − ½I3ω3²)/I1, β = 2mgl/I1.
 * f = β(u − u1)(u2 − u)(u3 − u) with u1 ≤ u2 ≤ 1 < u3, so θ oscillates between acos(u2) and acos(u1)
 * with period T = 4K(m)/√(β(u3 − u1)), m = (u2 − u1)/(u3 − u1).
 * Needs mgl > 0 and 0 < θ < π. Returns { thetaMin, thetaMax, period, roots: [u1, u2, u3], pPhi, pPsi, energy }.
 */
export function lagrangeTopNutation({ I1, I3, mgl, theta, thetaDot = 0, phiDot = 0, psiDot = 0 }) {
  positive(I1, 'I1', 'lagrangeTopNutation'); positive(I3, 'I3', 'lagrangeTopNutation'); positive(mgl, 'mgl', 'lagrangeTopNutation');
  finite([theta, thetaDot, phiDot, psiDot], 4, 'lagrangeTopNutation');
  if (!(theta > 0 && theta < Math.PI)) throw new Error('lagrangeTopNutation: need 0 < θ < π (a vertical top does not nutate)');
  const u0 = Math.cos(theta), s2 = 1 - u0 * u0, w3 = phiDot * u0 + psiDot;
  const pPsi = I3 * w3, pPhi = I1 * phiDot * s2 + pPsi * u0;
  const energy = 0.5 * I1 * (thetaDot ** 2 + phiDot ** 2 * s2) + 0.5 * I3 * w3 * w3 + mgl * u0;
  const a = pPsi / I1, b = pPhi / I1, alpha = 2 * (energy - 0.5 * I3 * w3 * w3) / I1, beta = 2 * mgl / I1;
  const C2 = -(alpha + a * a) / beta, C1 = (2 * a * b - beta) / beta, C0 = (alpha - b * b) / beta;
  const f = u => ((u + C2) * u + C1) * u + C0;               // f(u)/β, monic cubic
  let hi = 2;
  while (f(hi) <= 0) hi *= 2;
  const u3 = f(1) >= 0 ? 1 : brent(f, 1, hi, 1e-15);
  const p = C2 + u3, r = C1 + p * u3;                        // f/β = (u − u3)(u² + p u + r)
  const disc = Math.max(0, p * p - 4 * r), sq = Math.sqrt(disc);
  let u1 = (-p - sq) / 2, u2 = (-p + sq) / 2;
  const polish = x => { let y = x; for (let i = 0; i < 3; i++) { const d = (3 * y + 2 * C2) * y + C1; if (Math.abs(d) > 1e-8) y -= f(y) / d; } return y; };
  u1 = Math.max(-1, Math.min(u0, polish(u1))); u2 = Math.min(1, Math.max(u0, polish(u2)));
  const m = (u2 - u1) / (u3 - u1);
  return {
    thetaMin: Math.acos(u2), thetaMax: Math.acos(u1), period: 4 * ellipticK(m) / Math.sqrt(beta * (u3 - u1)),
    roots: [u1, u2, u3], pPhi, pPsi, energy,
  };
}

/**
 * Steady precession of a symmetric top at fixed tilt θ with spin ω3: the precession rates Ω solving
 * I1 cosθ Ω² − I3 ω3 Ω + mgl = 0 (θ̈ = 0 in Lagrange's equation). slow → mgl/(I3 ω3) for a fast top
 * (the gyroscope formula, returned as `approx`); fast ≈ I3 ω3/(I1 cosθ). Throws if
 * I3²ω3² < 4 I1 mgl cosθ (spin too slow: no steady precession at that tilt).
 */
export function steadyPrecession({ I1, I3, mgl, omega3, theta }) {
  positive(I1, 'I1', 'steadyPrecession'); positive(I3, 'I3', 'steadyPrecession');
  finite([mgl, omega3, theta], 3, 'steadyPrecession');
  if (omega3 === 0) throw new Error('steadyPrecession: omega3 must be non-zero');
  const c = Math.cos(theta), p = I3 * omega3, disc = p * p - 4 * I1 * mgl * c;
  if (disc < 0) throw new Error('steadyPrecession: spin too slow for steady precession at this tilt');
  const big = p + Math.sign(p) * Math.sqrt(disc);
  return { slow: 2 * mgl / big, fast: c === 0 ? Infinity : big / (2 * I1 * c), approx: mgl / p };
}
