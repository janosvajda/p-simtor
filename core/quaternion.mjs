// Quaternions for 3-D rotations, independent of subject.
//
// Convention: q = [w, x, y, z] (scalar first, as in core/projection.mjs). A unit quaternion
// q = [cos(θ/2), sin(θ/2) n̂] is the ACTIVE rotation by angle θ (right-handed) about the unit axis n̂:
// v' = q v q*. For a body with orientation q, qToMatrix(q) maps body coordinates to space
// coordinates, the body angular velocity ω obeys q̇ = ½ q ω and the space one q̇ = ½ ω q.
// q and −q are the same rotation; functions that return a canonical form choose w ≥ 0.
// All functions are pure, deterministic and return new Float64Array values (matrices are arrays
// of three Float64Array rows, as in core/linalg.mjs). Invalid input throws.

const finite = (v, n, who) => {
  if (!v || v.length !== n || !Array.from(v).every(Number.isFinite)) throw new Error(`${who}: expected ${n} finite values`);
};

/** Hamilton product a ⊗ b (rotation b first, then a). Returns out (new Float64Array(4) by default). */
export function qmul(a, b, out = new Float64Array(4)) {
  const [aw, ax, ay, az] = a, [bw, bx, by, bz] = b;
  out[0] = aw * bw - ax * bx - ay * by - az * bz;
  out[1] = aw * bx + ax * bw + ay * bz - az * by;
  out[2] = aw * by - ax * bz + ay * bw + az * bx;
  out[3] = aw * bz + ax * by - ay * bx + az * bw;
  return out;
}

/** Conjugate q* = [w, −x, −y, −z] (the inverse rotation for a unit quaternion). */
export function qconj(q) { return Float64Array.of(q[0], -q[1], -q[2], -q[3]); }

/** Euclidean norm |q|. */
export function qnorm(q) { return Math.hypot(q[0], q[1], q[2], q[3]); }

/** q / |q|. Throws for a zero or non-finite quaternion (no rotation is defined). */
export function qnormalize(q, out = new Float64Array(4)) {
  finite(q, 4, 'qnormalize');
  const n = qnorm(q);
  if (!(n > 0 && Number.isFinite(n))) throw new Error('qnormalize: zero quaternion has no direction');
  for (let i = 0; i < 4; i++) out[i] = q[i] / n;
  return out;
}

/** Multiplicative inverse q* / |q|² (any non-zero quaternion). */
export function qinverse(q) {
  finite(q, 4, 'qinverse');
  const n2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
  if (!(n2 > 0)) throw new Error('qinverse: zero quaternion has no inverse');
  return Float64Array.of(q[0] / n2, -q[1] / n2, -q[2] / n2, -q[3] / n2);
}

/**
 * Unit quaternion for a rotation by `angle` (radians, right-handed) about `axis` (any length > 0;
 * it is normalised). Throws for a zero or non-finite axis.
 */
export function qFromAxisAngle(axis, angle) {
  finite(axis, 3, 'qFromAxisAngle');
  if (!Number.isFinite(angle)) throw new Error('qFromAxisAngle: angle must be finite');
  const n = Math.hypot(axis[0], axis[1], axis[2]);
  if (!(n > 0)) throw new Error('qFromAxisAngle: axis must be non-zero');
  const s = Math.sin(angle / 2) / n;
  return Float64Array.of(Math.cos(angle / 2), s * axis[0], s * axis[1], s * axis[2]);
}

/**
 * Axis and angle of the rotation q (normalised first): { axis (unit), angle ∈ [0, π] }.
 * The angle uses atan2, accurate near 0 and π. For the identity the axis is [1, 0, 0] by convention.
 */
export function qToAxisAngle(q) {
  const [w, x, y, z] = qnormalize(q), sgn = w < 0 ? -1 : 1;
  const s = Math.hypot(x, y, z), angle = 2 * Math.atan2(s, Math.abs(w));
  if (s === 0) return { axis: Float64Array.of(1, 0, 0), angle: 0 };
  return { axis: Float64Array.of(sgn * x / s, sgn * y / s, sgn * z / s), angle };
}

/**
 * Rotation matrix R(q) (rows) with R v = q v q*: orthogonal, det +1. q is normalised first, so
 * a slightly drifted quaternion still gives an orthogonal matrix.
 */
export function qToMatrix(q) {
  const [w, x, y, z] = qnormalize(q);
  return [
    Float64Array.of(1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)),
    Float64Array.of(2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)),
    Float64Array.of(2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)),
  ];
}

/**
 * Unit quaternion (w ≥ 0) of a proper rotation matrix (Shepperd's method: divides by the largest
 * of 4w², 4x², 4y², 4z², so it is accurate for every angle including π). Throws if R is not
 * orthogonal with det +1 to `tolerance` (max-norm of RᵀR − I and |det − 1|).
 */
export function qFromMatrix(R, { tolerance = 1e-9 } = {}) {
  if (!R || R.length !== 3) throw new Error('qFromMatrix: expected a 3×3 matrix');
  for (const row of R) finite(row, 3, 'qFromMatrix');
  let off = 0;
  for (let i = 0; i < 3; i++) {for (let j = 0; j < 3; j++) {
    let s = 0;
    for (let k = 0; k < 3; k++) s += R[k][i] * R[k][j];
    off = Math.max(off, Math.abs(s - (i === j ? 1 : 0)));
  }}
  const det = R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);
  if (off > tolerance || Math.abs(det - 1) > tolerance) throw new Error('qFromMatrix: matrix is not a proper rotation (orthogonal, det +1)');
  const tr = R[0][0] + R[1][1] + R[2][2];
  const cand = [1 + tr, 1 + R[0][0] - R[1][1] - R[2][2], 1 - R[0][0] + R[1][1] - R[2][2], 1 - R[0][0] - R[1][1] + R[2][2]];
  const k = cand.indexOf(Math.max(...cand)), r = Math.sqrt(cand[k]), f = 0.5 / r;
  let q;
  if (k === 0) q = [r / 2, (R[2][1] - R[1][2]) * f, (R[0][2] - R[2][0]) * f, (R[1][0] - R[0][1]) * f];
  else if (k === 1) q = [(R[2][1] - R[1][2]) * f, r / 2, (R[0][1] + R[1][0]) * f, (R[0][2] + R[2][0]) * f];
  else if (k === 2) q = [(R[0][2] - R[2][0]) * f, (R[0][1] + R[1][0]) * f, r / 2, (R[1][2] + R[2][1]) * f];
  else q = [(R[1][0] - R[0][1]) * f, (R[0][2] + R[2][0]) * f, (R[1][2] + R[2][1]) * f, r / 2];
  if (q[0] < 0) q = q.map(v => -v);
  return qnormalize(q);
}

/**
 * Rotate the 3-vector v by the UNIT quaternion q: q v q* (v' = v + w t + u × t, t = 2 u × v).
 * q is not renormalised here (speed); a non-unit q scales the result by |q|².
 */
export function qRotate(q, v, out = new Float64Array(3)) {
  const [w, x, y, z] = q, [vx, vy, vz] = v;
  const tx = 2 * (y * vz - z * vy), ty = 2 * (z * vx - x * vz), tz = 2 * (x * vy - y * vx);
  out[0] = vx + w * tx + (y * tz - z * ty);
  out[1] = vy + w * ty + (z * tx - x * tz);
  out[2] = vz + w * tz + (x * ty - y * tx);
  return out;
}

/**
 * Spherical linear interpolation between unit quaternions a (t = 0) and b (t = 1) along the
 * shorter arc (b is negated if a·b < 0): constant angular velocity, unit norm for all t.
 * Falls back to normalised linear interpolation when the arc is below 1e-6 rad.
 */
export function qSlerp(a, b, t) {
  finite(a, 4, 'qSlerp'); finite(b, 4, 'qSlerp');
  if (!Number.isFinite(t)) throw new Error('qSlerp: t must be finite');
  let d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
  const sgn = d < 0 ? -1 : 1;
  d = Math.min(1, Math.abs(d));
  const theta = Math.acos(d);
  let wa, wb;
  if (theta < 1e-6) { wa = 1 - t; wb = t; } else { const s = Math.sin(theta); wa = Math.sin((1 - t) * theta) / s; wb = Math.sin(t * theta) / s; }
  return qnormalize([0, 1, 2, 3].map(i => wa * a[i] + sgn * wb * b[i]));
}

/** Quaternion exponential e^q = e^w [cos|v|, sin|v| v/|v|] (any quaternion; exp of a pure quaternion is unit). */
export function qExp(q) {
  finite(q, 4, 'qExp');
  const n = Math.hypot(q[1], q[2], q[3]), e = Math.exp(q[0]);
  const s = n > 1e-8 ? Math.sin(n) / n : 1 - n * n / 6;
  return Float64Array.of(e * Math.cos(n), e * s * q[1], e * s * q[2], e * s * q[3]);
}

/**
 * Principal quaternion logarithm [ln|q|, θ v/|v|], θ = atan2(|v|, w) ∈ [0, π]. exp(log q) = q.
 * Throws for q = 0. For q on the negative real axis the direction is undefined: [ln|q|, π, 0, 0].
 */
export function qLog(q) {
  finite(q, 4, 'qLog');
  const n = qnorm(q), v = Math.hypot(q[1], q[2], q[3]);
  if (!(n > 0)) throw new Error('qLog: logarithm of zero');
  const th = Math.atan2(v, q[0]);
  if (v === 0) return Float64Array.of(Math.log(n), q[0] < 0 ? Math.PI : 0, 0, 0);
  return Float64Array.of(Math.log(n), th * q[1] / v, th * q[2] / v, th * q[3] / v);
}

/**
 * Orientation after time dt under angular velocity ω (rad per time unit, 3-vector), exact for
 * constant ω: body frame q ⊗ exp(½ ω dt) (solves q̇ = ½ q ω), space frame exp(½ ω dt) ⊗ q
 * (solves q̇ = ½ ω q). For time-varying ω it is the first-order Lie–Euler step (use midpoint ω
 * for second order). The result is renormalised, so long runs do not drift off unit norm.
 */
export function qIntegrate(q, omega, dt, { frame = 'body' } = {}) {
  finite(q, 4, 'qIntegrate'); finite(omega, 3, 'qIntegrate');
  if (!Number.isFinite(dt)) throw new Error('qIntegrate: dt must be finite');
  if (frame !== 'body' && frame !== 'space') throw new Error(`qIntegrate: frame must be 'body' or 'space', got "${frame}"`);
  const e = qExp([0, omega[0] * dt / 2, omega[1] * dt / 2, omega[2] * dt / 2]);
  return qnormalize(frame === 'body' ? qmul(q, e) : qmul(e, q));
}
