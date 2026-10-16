import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  qmul, qconj, qnorm, qnormalize, qinverse, qFromAxisAngle, qToAxisAngle, qToMatrix, qFromMatrix,
  qRotate, qSlerp, qExp, qLog, qIntegrate,
} from '../core/quaternion.mjs';
import { rng } from '../core/numerics.mjs';
import { odeStepper } from '../core/ode.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const closeVec = (a, b, tol, msg) => a.forEach((v, i) => close(v, b[i], tol, `${msg ?? ''}[${i}]`));
const sameRotation = (a, b, tol) => { const s = Math.sign(a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]) || 1; closeVec(a, b.map(v => s * v), tol); };
const rand = rng(7);
const randomQ = () => qnormalize([rand.normal(), rand.normal(), rand.normal(), rand.normal()]);
const matMul = (A, B) => A.map(row => [0, 1, 2].map(j => row[0] * B[0][j] + row[1] * B[1][j] + row[2] * B[2][j]));
const matVec = (A, v) => A.map(row => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);

test('Hamilton product: i² = j² = k² = ijk = −1, norm is multiplicative, inverse and conjugate', () => {
  const i = [0, 1, 0, 0], j = [0, 0, 1, 0], k = [0, 0, 0, 1];
  for (const u of [i, j, k]) assert.deepEqual(Array.from(qmul(u, u)), [-1, 0, 0, 0]);
  assert.deepEqual(Array.from(qmul(qmul(i, j), k)), [-1, 0, 0, 0]);
  assert.deepEqual(Array.from(qmul(i, j)), [0, 0, 0, 1]);
  const a = [1, -2, 0.5, 3], b = [0.3, 1, -1, 2];
  close(qnorm(qmul(a, b)), qnorm(a) * qnorm(b), 1e-13);
  closeVec(qmul(a, qinverse(a)), [1, 0, 0, 0], 1e-15);
  closeVec(qconj(qmul(a, b)), qmul(qconj(b), qconj(a)), 1e-14);
  assert.throws(() => qnormalize([0, 0, 0, 0]), /zero/);
  assert.throws(() => qinverse([0, 0, 0, 0]), /zero/);
  assert.throws(() => qnormalize([1, NaN, 0, 0]), /finite/);
});

test('rotation matrices are orthogonal with det +1 and compose like quaternions', () => {
  for (let n = 0; n < 50; n++) {
    const a = randomQ(), b = randomQ(), A = qToMatrix(a), B = qToMatrix(b);
    const AtA = matMul(A.map((_, r) => A.map(row => row[r])), A);
    AtA.forEach((row, r) => closeVec(row, [0, 1, 2].map(c => (r === c ? 1 : 0)), 4e-15));
    const det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    close(det, 1, 1e-14);
    matMul(A, B).forEach((row, r) => closeVec(row, qToMatrix(qmul(a, b))[r], 1e-14, 'R(ab) = R(a)R(b)'));
    const v = [rand.normal(), rand.normal(), rand.normal()];
    closeVec(qRotate(a, v), matVec(A, v), 1e-14, 'qRotate = R v');
    const pv = qmul(qmul(a, [0, ...v]), qconj(a));
    closeVec(qRotate(a, v), pv.slice(1), 1e-14, 'q v q*');
    sameRotation(qFromMatrix(A), a, 1e-14);
  }
});

test('axis–angle: right-handed, round trip, and accurate near 0 and π', () => {
  const q = qFromAxisAngle([0, 0, 2], Math.PI / 2);
  closeVec(qRotate(q, [1, 0, 0]), [0, 1, 0], 1e-15, '+90° about z takes x to y');
  for (const angle of [1e-9, 0.3, 2, Math.PI - 1e-9, Math.PI]) {
    const axis = [1, -2, 0.5].map(v => v / Math.hypot(1, 2, 0.5));
    const r = qToAxisAngle(qFromAxisAngle(axis, angle));
    close(r.angle, angle, 1e-15 * Math.max(1, angle) + 1e-16, 'angle');
    if (angle < Math.PI) closeVec(r.axis, axis, 1e-7, 'axis');
    const R = qToMatrix(qFromAxisAngle(axis, angle));
    sameRotation(qFromMatrix(R), qFromAxisAngle(axis, angle), 1e-14);   // Shepperd branch near π
  }
  assert.equal(qToAxisAngle([1, 0, 0, 0]).angle, 0);
  assert.throws(() => qFromAxisAngle([0, 0, 0], 1), /non-zero/);
  assert.throws(() => qFromMatrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]]), /proper rotation/);
});

test('slerp: endpoints, unit norm, constant angular speed along the shorter arc', () => {
  const a = randomQ(), b = randomQ();
  closeVec(qSlerp(a, b, 0), a, 1e-15);
  sameRotation(qSlerp(a, b, 1), b, 1e-15);
  const total = qToAxisAngle(qmul(qconj(a), b)).angle;
  assert.ok(total <= Math.PI);
  for (const t of [0.1, 0.25, 0.5, 0.9]) {
    const s = qSlerp(a, b, t);
    close(qnorm(s), 1, 1e-15);
    close(qToAxisAngle(qmul(qconj(a), s)).angle, t * total, 1e-12, `angle at t = ${t}`);
  }
  sameRotation(qSlerp(a, b.map(v => -v), 0.3), qSlerp(a, b, 0.3), 1e-15);   // −b is the same rotation
  closeVec(qSlerp(a, a, 0.4), a, 1e-15);
});

test('exp and log are inverse; exp of a pure quaternion is the half-angle rotation', () => {
  for (let n = 0; n < 20; n++) {
    const q = [rand.normal(), rand.normal(), rand.normal(), rand.normal()];
    closeVec(qExp(qLog(q)), q, 1e-13);
  }
  const v = [0, 0.3, -0.4, 1.2];
  closeVec(qLog(qExp(v)), v, 1e-15);
  const th = 2 * Math.hypot(0.3, -0.4, 1.2);
  sameRotation(qExp(v), qFromAxisAngle([0.3, -0.4, 1.2], th), 1e-15);
  close(qnorm(qExp([0, 1e-12, 0, 0])), 1, 1e-16);
  assert.throws(() => qLog([0, 0, 0, 0]), /zero/);
});

test('qIntegrate: exact for constant ω (steps compose), body vs space frame, matches q̇ = ½qω', () => {
  const q0 = randomQ(), w = [0.7, -1.1, 0.4];
  let q = q0;
  for (let n = 0; n < 1000; n++) q = qIntegrate(q, w, 0.01);
  sameRotation(q, qIntegrate(q0, w, 10), 1e-12);
  sameRotation(q, qmul(q0, qFromAxisAngle(w, 10 * Math.hypot(...w))), 1e-12);
  // body-frame ω equals space-frame R(q)ω
  sameRotation(qIntegrate(q0, w, 0.8), qIntegrate(q0, qRotate(q0, w), 0.8, { frame: 'space' }), 1e-14);
  // time-dependent ω: Lie–Euler with midpoint ω converges to the RK solution of q̇ = ½ q ω(t)
  const omega = t => [Math.cos(t), 0.5 * t, -0.3];
  const f = (t, y) => { const d = qmul(y, [0, ...omega(t)]); return d.map(v => v / 2); };
  const ref = odeStepper({ f, y0: q0, h: 0.01, rtol: 1e-12, atol: 1e-14 });
  ref.advanceTo(2, Infinity, { exact: true });
  const run = h => { let p = q0; for (let n = 0; n < Math.round(2 / h); n++) p = qIntegrate(p, omega((n + 0.5) * h), h); return p; };
  const e1 = Math.max(...run(0.02).map((v, i) => Math.abs(v - ref.y[i]))), e2 = Math.max(...run(0.01).map((v, i) => Math.abs(v - ref.y[i])));
  close(Math.log2(e1 / e2), 2, 0.1, 'second order with midpoint ω');
  assert.throws(() => qIntegrate(q0, w, 0.1, { frame: 'lab' }), /frame/);
});
