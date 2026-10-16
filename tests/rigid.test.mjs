import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  inertiaBox, inertiaSphere, inertiaCylinder, inertiaCone, parallelAxis, compositeBody, principalAxes,
  ellipticF, jacobiElliptic, torqueFreeExact, freeRigidBody, heavyTop, lagrangeTopInitial,
  lagrangeTopNutation, steadyPrecession,
} from '../core/rigid.mjs';
import { qFromAxisAngle, qRotate, qmul, qToMatrix, qnormalize } from '../core/quaternion.mjs';
import { ellipticK, gaussLegendre } from '../core/numerics.mjs';
import { refinementOrder } from '../core/ode.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const closeVec = (a, b, tol, msg) => Array.from(a).forEach((v, i) => close(v, b[i], tol, `${msg ?? ''}[${i}]`));
const closeMat = (A, B, tol, msg) => A.forEach((row, i) => closeVec(row, B[i], tol, `${msg ?? ''}[${i}]`));

// ---------------------------------------------------------------- inertia

test('inertia tensors: textbook formulas, symmetry limits, independent slice integration', () => {
  closeMat(inertiaBox(6, [1, 1, 1]), [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1e-15, 'cube m a²/6');
  closeMat(inertiaSphere(5, 2), [[8, 0, 0], [0, 8, 0], [0, 0, 8]], 1e-14);
  close(inertiaSphere(3, 1, { shell: true })[0][0], 2, 1e-15);
  // cylinder with h² = 3r² is isotropic; a thin-walled tube tends to m r² about its axis
  const cyl = inertiaCylinder(2, 1, Math.sqrt(3));
  close(cyl[0][0], cyl[2][2], 1e-15); close(cyl[2][2], 1, 1e-15);
  close(inertiaCylinder(1, 1, 1, { innerRadius: 0.999999 })[2][2], 1, 1e-6);
  // cone: integrate thin discs (radius ρ(z) = r(1 − z/h), z from the base) about the centre of mass
  const m = 2, r = 0.7, h = 1.9, rho = 3 * m / (Math.PI * r * r * h);
  const disc = z => Math.PI * (r * (1 - z / h)) ** 2 * rho;   // dm/dz
  const zc = gaussLegendre(z => z * disc(z), 0, h) / m;
  const Iz = gaussLegendre(z => 0.5 * disc(z) * (r * (1 - z / h)) ** 2, 0, h);
  const Ix = gaussLegendre(z => disc(z) * ((r * (1 - z / h)) ** 2 / 4 + (z - zc) ** 2), 0, h);
  const cone = inertiaCone(m, r, h);
  close(cone.comFromBase, zc, 1e-14); close(cone.inertia[2][2], Iz, 1e-14); close(cone.inertia[0][0], Ix, 1e-14);
  assert.throws(() => inertiaBox(-1, [1, 1, 1]), /mass/);
  assert.throws(() => inertiaCylinder(1, 1, 1, { innerRadius: 1 }), /innerRadius/);
});

test('parallel axis and composites: two half-boxes make the box, a rod about its end is m L²/3', () => {
  const half = inertiaBox(1, [1, 2, 3]);                    // each half: 1 × 2 × 3, glued along x
  const body = compositeBody([{ mass: 1, inertia: half, position: [-0.5, 0, 0] }, { mass: 1, inertia: half, position: [0.5, 0, 0] }]);
  closeMat(body.inertia, inertiaBox(2, [2, 2, 3]), 1e-14);
  closeVec(body.com, [0, 0, 0], 1e-15);
  const rod = inertiaCylinder(3, 1e-9, 2);
  close(parallelAxis(rod, 3, [0, 0, 1])[0][0], 3 * 4 / 3, 1e-12);
  // a rotated part: a box turned 90° about z swaps its x and y moments
  const turned = compositeBody([{ mass: 1, inertia: inertiaBox(1, [1, 2, 3]), position: [0, 0, 0], orientation: qFromAxisAngle([0, 0, 1], Math.PI / 2) }]);
  closeMat(turned.inertia, inertiaBox(1, [2, 1, 3]), 1e-15);
  // a cavity: a sphere minus a concentric smaller sphere equals the thick-shell formula
  const hollow = compositeBody([{ mass: 8, inertia: inertiaSphere(8, 2), position: [0, 0, 0] }, { mass: -1, inertia: inertiaSphere(1, 1).map(row => row.map(v => -v)), position: [0, 0, 0] }]);
  close(hollow.inertia[0][0], 2 / 5 * 7 * (2 ** 5 - 1) / (2 ** 3 - 1), 1e-13);
  assert.throws(() => compositeBody([]), /non-empty/);
});

test('principal axes recover moments and a right-handed frame of a rotated tensor', () => {
  const q = qnormalize([0.3, -0.5, 0.7, 0.2]), R = qToMatrix(q), D = [1, 2.5, 3];
  const I = [0, 1, 2].map(i => [0, 1, 2].map(j => R[i][0] * D[0] * R[j][0] + R[i][1] * D[1] * R[j][1] + R[i][2] * D[2] * R[j][2]));
  const p = principalAxes(I);
  closeVec(p.moments, D, 1e-13);
  p.axes.forEach((a, k) => close(Math.abs(a[0] * R[0][k] + a[1] * R[1][k] + a[2] * R[2][k]), 1, 1e-12, `axis ${k}`));
  const Rp = p.rotation;
  const det = Rp[0][0] * (Rp[1][1] * Rp[2][2] - Rp[1][2] * Rp[2][1]) - Rp[0][1] * (Rp[1][0] * Rp[2][2] - Rp[1][2] * Rp[2][0]) + Rp[0][2] * (Rp[1][0] * Rp[2][1] - Rp[1][1] * Rp[2][0]);
  close(det, 1, 1e-13);
  closeVec(qRotate(p.quaternion, [0, 0, 1]), p.axes[2], 1e-13);
  assert.throws(() => principalAxes([[1, 0, 0], [0, -1, 0], [0, 0, 1]]), /positive definite/);
});

// ---------------------------------------------------------------- elliptic functions

test('Jacobi elliptic functions: A&S 16.5 special values, limits, identities, F inverts am', () => {
  for (const m of [0.1, 0.5, 0.81, 0.99, 0.999999]) {
    const K = ellipticK(m), kp = Math.sqrt(1 - m);
    const atK = jacobiElliptic(K, m);
    close(atK.sn, 1, 1e-14); close(atK.cn, 0, 1e-8); close(atK.dn, kp, 1e-8, 'dn(K) = k′');
    // A&S 16.5.2: sn(K/2) = 1/√(1 + k′), cn(K/2) = √(k′/(1 + k′)), dn(K/2) = √k′
    const half = jacobiElliptic(K / 2, m);
    close(half.sn, 1 / Math.sqrt(1 + kp), 1e-13); close(half.cn, Math.sqrt(kp / (1 + kp)), 1e-13); close(half.dn, Math.sqrt(kp), 1e-13);
    for (const u of [0.3, 1.7, -2.2, 11.4, 4 * K + 0.5]) {
      const { sn, cn, dn, am } = jacobiElliptic(u, m);
      close(sn * sn + cn * cn, 1, 1e-14); close(dn * dn + m * sn * sn, 1, 1e-14);
      close(ellipticF(am, m), u, 1e-12 * Math.max(1, Math.abs(u)), `F(am(u)) = u, m = ${m}`);
      const d = 1e-5, p = jacobiElliptic(u + d, m), q = jacobiElliptic(u - d, m);
      close((p.sn - q.sn) / (2 * d), cn * dn, 1e-8, 'sn′ = cn dn');
    }
  }
  const s = jacobiElliptic(0.8, 0); close(s.sn, Math.sin(0.8), 1e-16);
  const t = jacobiElliptic(0.8, 1); close(t.sn, Math.tanh(0.8), 1e-16); close(t.dn, 1 / Math.cosh(0.8), 1e-16);
  // F by quadrature, and F(π/2) = K
  close(ellipticF(1.1, 0.7), gaussLegendre(th => 1 / Math.sqrt(1 - 0.7 * Math.sin(th) ** 2), 0, 1.1, 40), 1e-14);
  close(ellipticF(Math.PI / 2, 0.7), ellipticK(0.7), 1e-14);
  // the general version (special.mjs) also covers m > 1 by the reciprocal-parameter transform (A&S 16.11)
  { const k = Math.sqrt(1.5), r = jacobiElliptic(k, 1 / 1.5); assert.ok(Math.abs(jacobiElliptic(1, 1.5).sn - r.sn / k) < 1e-14); }
});

// ---------------------------------------------------------------- torque-free motion

const I123 = [1, 2, 3];

test('exact torque-free solution satisfies Euler’s equations, conserves E and |L|, matches Landau’s k²', () => {
  for (const L0 of [[0.3, 0.2, 1], [1, 0.3, -0.2], [-0.4, -0.9, 0.6], [0.2, -0.5, -0.7]]) {
    const sol = torqueFreeExact({ inertia: I123, L0 });
    closeVec(sol.at(0), L0, 1e-13, 'initial state');
    for (const t of [0.7, 3.1, 10]) {
      const L = sol.at(t), d = 1e-5, Lp = sol.at(t + d), Lm = sol.at(t - d);
      const w = L.map((v, i) => v / I123[i]);
      const rhs = [L[1] * w[2] - L[2] * w[1], L[2] * w[0] - L[0] * w[2], L[0] * w[1] - L[1] * w[0]];   // L̇ = L × Ω
      closeVec(Lp.map((v, i) => (v - Lm[i]) / (2 * d)), rhs, 1e-8, 'Euler');
      close(L.reduce((s, v) => s + v * v, 0), sol.L2, 1e-13);
      close(L.reduce((s, v, i) => s + v * v / (2 * I123[i]), 0), sol.energy, 1e-13);
    }
    closeVec(sol.at(sol.period), L0, 1e-11, 'periodic');
    const [I1, I2, I3] = I123, E2 = 2 * sol.energy, M2 = sol.L2;
    const k2 = M2 > E2 * I2 ? (I2 - I1) * (E2 * I3 - M2) / ((I3 - I2) * (M2 - E2 * I1)) : (I3 - I2) * (E2 * I1 - M2) / ((I2 - I1) * (M2 - E2 * I3));
    close(sol.m, k2, 1e-13, 'Landau §37 k²');
  }
  assert.throws(() => torqueFreeExact({ inertia: I123, L0: [0, 1, 0] }), /separatrix/);
  assert.throws(() => torqueFreeExact({ inertia: [3, 2, 1], L0: [0, 1, 0] }), /I1 ≤ I2 ≤ I3/);
});

test('Lie–Poisson integrator: |L| and spatial L to round-off, bounded energy, orders 2 and 4, reversible', () => {
  const L0 = [0.3, 0.9, 0.4], sol = torqueFreeExact({ inertia: I123, L0 });
  const s = freeRigidBody({ inertia: I123, L0, h: 0.02 });
  const S0 = s.spatialMomentum(), E0 = s.energy(), n0 = Math.hypot(...L0);
  let eFirst = 0, eSecond = 0;
  for (let i = 1; i <= 50000; i++) {
    s.step();
    const dE = Math.abs(s.energy() - E0);
    if (i <= 25000) eFirst = Math.max(eFirst, dE); else eSecond = Math.max(eSecond, dE);
  }
  close(Math.hypot(...s.L), n0, 1e-13, '|L|');
  closeVec(s.spatialMomentum(), S0, 1e-12, 'spatial L');
  assert.ok(eFirst < 1e-4 && eSecond < 1.2 * eFirst, `bounded energy error, no drift: ${eFirst} ${eSecond}`);
  // against the exact Jacobi-elliptic solution, with the observed order
  for (const [method, order] of [['strang', 2], ['yoshida4', 4]]) {
    const run = h => { const r = freeRigidBody({ inertia: I123, L0, h, method }); r.advanceTo(4); return r.L; };
    close(refinementOrder(run, 0.04).order, order, 0.15, method);
    closeVec(run(0.005), sol.at(4), method === 'strang' ? 1e-4 : 1e-9, `${method} vs exact`);
  }
  // time reversal retraces the orientation and momentum
  const r = freeRigidBody({ inertia: I123, L0, q0: qnormalize([1, 0.2, -0.1, 0.3]), h: 0.05 });
  const q0 = Float64Array.from(r.q);
  for (let i = 0; i < 400; i++) r.step();
  r.h = -0.05;
  for (let i = 0; i < 400; i++) r.step();
  closeVec(r.L, L0, 1e-12); closeVec(r.q, q0, 1e-12);
  assert.throws(() => freeRigidBody({ inertia: I123, L0, h: 0.1, method: 'euler' }), /unknown method/);
});

test('symmetric body: splitting is exact; the symmetry axis cones about L at rate |L|/I1', () => {
  const I = [1, 1, 2.5], L0 = [0.4, -0.3, 1.2], q0 = qnormalize([0.9, 0.1, 0.3, -0.2]);
  const s = freeRigidBody({ inertia: I, L0, q0, h: 0.37 });   // large step on purpose
  const Ls = s.spatialMomentum(), n = Math.hypot(...Ls), axis0 = qRotate(q0, [0, 0, 1]);
  s.advanceTo(37);
  const expected = qRotate(qFromAxisAngle(Ls, n * s.t / I[0]), axis0);
  closeVec(qRotate(s.q, [0, 0, 1]), expected, 1e-12, 'body axis in space');
  // body-frame precession of L about e3 at (I3 − I1)/I1 · ω3 (exactly): L̇ = L × Ω = κ (L2, −L1, 0), κ = (1/I3 − 1/I1) L3
  const w3 = L0[2] / I[2], ang = (I[2] - I[0]) / I[0] * w3 * s.t;
  closeVec(s.L, [L0[0] * Math.cos(ang) - L0[1] * Math.sin(ang), L0[0] * Math.sin(ang) + L0[1] * Math.cos(ang), L0[2]], 1e-12);
});

test('Dzhanibekov effect: intermediate-axis flips recur at half the elliptic period 2K(m)/λ', () => {
  const L0 = [0.01, 1, 0.01], exact = torqueFreeExact({ inertia: I123, L0 });
  const s = freeRigidBody({ inertia: I123, L0, h: 0.002, method: 'yoshida4' });
  const flips = [];
  let prev = s.L[1];
  while (flips.length < 4) {
    s.step();
    if (Math.sign(s.L[1]) !== Math.sign(prev)) flips.push(s.t - s.h * s.L[1] / (s.L[1] - prev));
    prev = s.L[1];
  }
  for (let i = 1; i < flips.length; i++) close(flips[i] - flips[i - 1], exact.period / 2, 1e-6, 'flip interval');
  assert.ok(exact.m > 0.99, `near the separatrix: m = ${exact.m}`);
  // stable rotation about the largest axis does not flip
  const st = freeRigidBody({ inertia: I123, L0: [0.01, 0.01, 1], h: 0.01 });
  st.advanceTo(200);
  assert.ok(st.L[2] > 0.999);
});

// ---------------------------------------------------------------- heavy top

const top = { I1: 1, I3: 0.5, mgl: 1 };

test('heavy top: Casimirs |Γ|, L·Γ and L3 to round-off, energy bounded', () => {
  const init = lagrangeTopInitial({ ...top, theta: 0.6, psiDot: 4, phiDot: 0.3, thetaDot: 0.2 });
  const s = heavyTop({ inertia: [1, 1, 0.5], mgl: 1, L0: init.L, q0: init.q, h: 0.005 });
  const G0 = s.gamma(), C0 = s.L[0] * G0[0] + s.L[1] * G0[1] + s.L[2] * G0[2], E0 = s.energy(), L30 = s.L[2];
  let worst = 0;
  for (let i = 0; i < 20000; i++) { s.step(); worst = Math.max(worst, Math.abs(s.energy() - E0)); }
  const G = s.gamma();
  close(Math.hypot(...G), 1, 1e-14, '|Γ|');
  close(s.L[0] * G[0] + s.L[1] * G[1] + s.L[2] * G[2], C0, 1e-12, 'L·Γ');
  close(s.spatialMomentum()[2], C0, 1e-12, 'vertical spatial L');
  close(s.L[2], L30, 1e-12, 'L3');
  assert.ok(worst < 1e-4, `energy ${worst}`);
});

test('Lagrange top initial state: q̇ = ½ q ω for the Euler-angle rates', () => {
  const a = { ...top, theta: 0.7, phi: 0.4, psi: 1.1, thetaDot: 0.3, phiDot: -0.8, psiDot: 2.5 }, d = 1e-6;
  const at = dt => lagrangeTopInitial({ ...a, theta: a.theta + a.thetaDot * dt, phi: a.phi + a.phiDot * dt, psi: a.psi + a.psiDot * dt }).q;
  const qdot = at(d).map((v, i) => (v - at(-d)[i]) / (2 * d)), init = lagrangeTopInitial(a);
  closeVec(qdot, qmul(init.q, [0, ...init.omega]).map(v => v / 2), 1e-8);
  closeVec(qRotate(init.q, [0, 0, 1]), [Math.sin(0.7) * Math.sin(0.4), -Math.sin(0.7) * Math.cos(0.4), Math.cos(0.7)], 1e-15);
});

const axisAngles = s => { const e = qRotate(s.q, [0, 0, 1]); return { theta: Math.acos(e[2]), phi: Math.atan2(e[0], -e[1]) }; };

test('Lagrange top nutation: integrated θ range and period match the cubic-root / elliptic-K solution', () => {
  const a = { ...top, theta: 0.5, psiDot: 6 };               // released with no precession: cusped nutation
  const exact = lagrangeTopNutation(a), init = lagrangeTopInitial(a);
  close(exact.thetaMin, 0.5, 1e-12, 'starts at a turning point');
  const s = heavyTop({ inertia: [1, 1, 0.5], mgl: 1, L0: init.L, q0: init.q, h: 0.001, method: 'yoshida4' });
  let tMax = 0, thMax = 0, thMin = Infinity, prevTh = 0.5, rising = true;
  const minima = [];
  while (minima.length < 3) {
    s.step();
    const { theta } = axisAngles(s);
    thMax = Math.max(thMax, theta); thMin = Math.min(thMin, theta);
    if (rising && theta < prevTh) { rising = false; tMax = s.t; }
    if (!rising && theta > prevTh) { rising = true; minima.push(s.t); }
    prevTh = theta;
  }
  assert.ok(tMax > 0);
  close(thMax, exact.thetaMax, 1e-6, 'θ_max'); close(thMin, exact.thetaMin, 1e-6, 'θ_min');
  close((minima[2] - minima[0]) / 2, exact.period, 2e-3, 'nutation period');
  assert.throws(() => lagrangeTopNutation({ ...a, theta: 0 }), /vertical/);
});

test('steady precession: θ stays fixed at the slow root; a fast top precesses at ≈ mgl/(I3 ω3)', () => {
  const theta = 0.8, omega3 = 8, sp = steadyPrecession({ ...top, omega3, theta });
  close(top.I1 * Math.cos(theta) * sp.slow ** 2 - top.I3 * omega3 * sp.slow + top.mgl, 0, 1e-13);
  close(top.I1 * Math.cos(theta) * sp.fast ** 2 - top.I3 * omega3 * sp.fast + top.mgl, 0, 1e-12);
  const init = lagrangeTopInitial({ ...top, theta, phiDot: sp.slow, psiDot: omega3 - sp.slow * Math.cos(theta) });
  const s = heavyTop({ inertia: [1, 1, 0.5], mgl: 1, L0: init.L, q0: init.q, h: 0.002, method: 'yoshida4' });
  s.advanceTo(20);
  const { theta: th, phi } = axisAngles(s);
  close(th, theta, 1e-8, 'no nutation');
  close(phi, ((sp.slow * s.t + Math.PI) % (2 * Math.PI)) - Math.PI, 1e-7, 'precession angle');
  // dropped fast top (φ̇0 = 0): mean precession over whole nutation periods → mgl/(I3ω3)
  const fast = { ...top, theta: 1, psiDot: 60 }, nut = lagrangeTopNutation(fast), f0 = lagrangeTopInitial(fast);
  const f = heavyTop({ inertia: [1, 1, 0.5], mgl: 1, L0: f0.L, q0: f0.q, h: 2e-4, method: 'yoshida4' });
  f.advanceTo(40 * nut.period);
  let phiEnd = axisAngles(f).phi;
  if (phiEnd < 0) phiEnd += 2 * Math.PI;
  const rate = phiEnd / f.t, approx = steadyPrecession({ ...top, omega3: 60, theta: 1 }).approx;
  assert.ok(Math.abs(rate / approx - 1) < 0.01, `mean precession ${rate} vs ${approx}`);
  assert.throws(() => steadyPrecession({ ...top, omega3: 0.1, theta: 0.2 }), /too slow/);
});
