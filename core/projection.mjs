// Seeing more than three dimensions: rotations in 4D, projections to 3D, and maps of the
// 3-sphere S³ (unit quaternions, e.g. SU(2) group elements). Pure maths, no rendering: feed the
// results to ui/gl.mjs. Every projection distorts something; name it beside the view.
//
//   rotation4({ xy, xz, xw, yz, yw, zw })   4×4 rotation (row-major) from plane angles, applied in that order
//   apply4(R, p)                            R·p for a 4-vector
//   project4(p, { mode, eye })              4D → 3D: 'perspective' from a 4D eye at w = eye (scale eye/(eye − w)),
//                                           or 'orthographic' (drop w). Connectivity is exact; lengths and angles are not
//   stereographic(q), inverseStereographic(x)   S³ ↔ ℝ³ from the pole −1 = [−1, 0, 0, 0] (conformal; the pole goes to ∞)
//   axisAngle(q)                            S³ → ball of radius π: θn̂ for q = [cos θ, sin θ n̂] (every point finite;
//                                           q and −q land at antipodal shell points when θ → π)
//   hopfMap(q), hopfFibre(p, n)             S³ → S² and the circle of S³ above a point p of S²; the fibre
//                                           through q is { exp(i t σ₃)·q }: left multiplication by a U(1) subgroup
//
// Coordinates: 4-vectors [x, y, z, w]; quaternions [a0, a1, a2, a3] ↔ a0 + i a·σ (as in lattice-gauge.mjs).

const PLANES = [['xy', 0, 1], ['xz', 0, 2], ['xw', 0, 3], ['yz', 1, 2], ['yw', 1, 3], ['zw', 2, 3]];

/** Product of plane rotations, in the order xy, xz, xw, yz, yw, zw (angles in radians). */
export function rotation4(angles = {}) {
  let R = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
  for (const [name, i, j] of PLANES) {
    const a = angles[name] ?? 0;
    if (!a) continue;
    const c = Math.cos(a), s = Math.sin(a), G = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
    G[4 * i + i] = c; G[4 * i + j] = -s; G[4 * j + i] = s; G[4 * j + j] = c;
    const out = new Array(16).fill(0);
    for (let r = 0; r < 4; r++) for (let k = 0; k < 4; k++) for (let m = 0; m < 4; m++) out[4 * r + k] += G[4 * r + m] * R[4 * m + k];
    R = out;
  }
  return R;
}

export function apply4(R, p) {
  return [0, 1, 2, 3].map(r => R[4 * r] * p[0] + R[4 * r + 1] * p[1] + R[4 * r + 2] * p[2] + R[4 * r + 3] * p[3]);
}

/** 4D → 3D. Perspective needs every point in front of the eye (w < eye); it throws otherwise. */
export function project4(p, { mode = 'perspective', eye = 3 } = {}) {
  if (mode === 'orthographic') return [p[0], p[1], p[2]];
  if (mode !== 'perspective') throw new Error("project4: mode must be 'perspective' or 'orthographic'");
  const d = eye - p[3];
  if (!(d > 0)) throw new Error(`project4: point at w = ${p[3]} is behind the 4D eye at w = ${eye}`);
  const f = eye / d;
  return [f * p[0], f * p[1], f * p[2]];
}

/** S³ → ℝ³ from the pole −1: x = a/(1 + a0). Points near −1 go far away (Infinity at the pole). */
export function stereographic(q) {
  const n = Math.hypot(q[0], q[1], q[2], q[3]), d = 1 + q[0] / n;
  return [q[1] / n / d, q[2] / n / d, q[3] / n / d];
}
export function inverseStereographic(x) {
  const r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2], d = 1 + r2;
  return [(1 - r2) / d, 2 * x[0] / d, 2 * x[1] / d, 2 * x[2] / d];
}

/** S³ → closed ball of radius π: θn̂ with q = [cos θ, sin θ n̂], θ = angle from the identity. */
export function axisAngle(q) {
  const s = Math.hypot(q[1], q[2], q[3]), th = Math.atan2(s, q[0]), f = s > 1e-300 ? th / s : 1;
  return [f * q[1], f * q[2], f * q[3]];
}

/** Hopf map S³ → S²: with z₁ = a0 + i a3, z₂ = a2 + i a1, p = (2 Re z₁z̄₂, 2 Im z₁z̄₂, |z₁|² − |z₂|²). */
export function hopfMap(q) {
  const [x1, y1, x2, y2] = [q[0], q[3], q[2], q[1]];
  return [2 * (x1 * x2 + y1 * y2), 2 * (y1 * x2 - x1 * y2), x1 * x1 + y1 * y1 - x2 * x2 - y2 * y2];
}

/** The great circle of S³ (n quaternions, closed) that hopfMap sends to the unit vector p. */
export function hopfFibre(p, n = 64) {
  const l = Math.hypot(...p), [p1, p2, p3] = p.map(v => v / l), out = [];
  for (let k = 0; k <= n; k++) {
    const t = 2 * Math.PI * k / n, c = Math.cos(t), s = Math.sin(t);
    let z1, z2;   // [re, im]
    if (p3 > -1 + 1e-12) {
      const a = Math.sqrt((1 + p3) / 2), b = 1 / Math.sqrt(2 * (1 + p3));
      z1 = [a * c, a * s];
      z2 = [b * (p1 * c + p2 * s), b * (p1 * s - p2 * c)];          // e^{it}(p1 − i p2)/√(2(1+p3))
    } else { z1 = [0, 0]; z2 = [c, s]; }
    out.push([z1[0], z2[1], z2[0], z1[1]]);                          // back to [a0, a1, a2, a3]
  }
  return out;
}
