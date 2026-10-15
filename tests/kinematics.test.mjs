import { test } from 'node:test';
import assert from 'node:assert/strict';
import { dot4, mass, add4, boost, restFrameBoost, randomUnitVector, randomRotation, orthonormalFrame, twoBodyDecay, rambo, acceptReject } from '../core/kinematics.mjs';
import { ksTest, meanSem } from '../core/data.mjs';
import { rng, adaptiveSimpson } from '../core/numerics.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const rotate = (R, p) => [p[0], ...R.map(r => r[0] * p[1] + r[1] * p[2] + r[2] * p[3])];
const det3 = R => R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);

test('dot4, mass, add4', () => {
  assert.equal(dot4([5, 1, 2, 3], [2, 1, 0, -1]), 10 - 1 + 3);
  close(mass([5, 3, 0, 0]), 4, 1e-15);
  assert.equal(mass([1, 1 + 1e-16, 0, 0]), 0);
  assert.deepEqual(add4([1, 2, 3, 4], [1, 1, 1, 1], [0, -2, 0, 1]), [2, 1, 4, 6]);
});

test('boost: rest particle gets velocity β; inverse boost; rest-frame boost; invariance of dot4', () => {
  const m = 2, beta = [0.3, -0.4, 0.5], g = 1 / Math.sqrt(1 - 0.5);
  const moved = boost([m, 0, 0, 0], beta);
  moved.forEach((v, k) => close(v, k ? g * m * beta[k - 1] : g * m, 1e-14, `component ${k}`));
  const rand = rng(1);
  for (let t = 0; t < 50; t++) {
    const p = [0, rand() - 0.5, rand() - 0.5, rand() - 0.5]; p[0] = Math.hypot(p[1], p[2], p[3], rand());
    const q = [0, rand() - 0.5, rand() - 0.5, rand() - 0.5]; q[0] = Math.hypot(q[1], q[2], q[3], rand());
    const b = randomUnitVector(rand).map(v => v * 0.95 * rand());
    boost(boost(p, b), b.map(v => -v)).forEach((v, k) => close(v, p[k], 1e-12, 'boost ∘ inverse'));
    close(dot4(boost(p, b), boost(q, b)), dot4(p, q), 1e-12, 'p·q under boost');
    const R = randomRotation(rand);
    close(dot4(rotate(R, p), rotate(R, q)), dot4(p, q), 1e-12, 'p·q under rotation');
    const P = add4(p, q), rest = boost(P, restFrameBoost(P));
    for (let k = 1; k < 4; k++) close(rest[k], 0, 1e-12, 'rest-frame momentum');
    close(rest[0], mass(P), 1e-12, 'rest-frame energy = invariant mass');
  }
  assert.throws(() => boost([1, 0, 0, 0], [1, 0, 0]), /boost/);
});

test('randomRotation: orthogonal, det +1, uniform; randomUnitVector isotropic', () => {
  const rand = rng(2), N = 20000;
  let mz = [0, 0, 0], z2 = 0;
  for (let i = 0; i < N; i++) {
    const R = randomRotation(rand);
    if (i < 200) {
      for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) close(R[0][a] * R[0][b] + R[1][a] * R[1][b] + R[2][a] * R[2][b], a === b ? 1 : 0, 1e-14, 'RᵀR = 1');
      close(det3(R), 1, 1e-14, 'det R');
    }
    const z = [R[0][2], R[1][2], R[2][2]];                 // R·ẑ
    mz = mz.map((v, k) => v + z[k] / N); z2 += z[2] * z[2] / N;
  }
  mz.forEach(v => close(v, 0, 4 / Math.sqrt(3 * N), '⟨R·ẑ⟩'));
  close(z2, 1 / 3, 4 * Math.sqrt(4 / 45 / N), '⟨(R·ẑ)_z²⟩');
  const cz = Array.from({ length: 5000 }, () => randomUnitVector(rand)[2]);
  assert.ok(ksTest(cz, c => (c + 1) / 2).p > 0.01, 'cos θ uniform');
});

test('orthonormalFrame is right-handed with e1 ⟂ n', () => {
  for (const n of [[0, 0, 1], [0, 0, -3], [1, 2, 3], [-5, 0.1, 0], [0, 1e-9, 1]]) {
    const [e1, e2, e3] = orthonormalFrame(n), r = Math.hypot(...n);
    const d = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    e3.forEach((v, k) => close(v, n[k] / r, 1e-15));
    close(d(e1, e1), 1, 1e-14); close(d(e2, e2), 1, 1e-14); close(d(e1, e3), 0, 1e-14); close(d(e2, e3), 0, 1e-14); close(d(e1, e2), 0, 1e-14);
    const c = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
    c.forEach((v, k) => close(v, e3[k], 1e-14, 'e1 × e2 = n̂'));
  }
  assert.throws(() => orthonormalFrame([0, 0, 0]), /orthonormalFrame/);
});

test('twoBodyDecay: Källén momentum, masses, conservation, isotropy', () => {
  const M = 5.28, m1 = 0.494, m2 = 0.140, rand = rng(3), N = 20000;
  const k = Math.sqrt((M ** 2 - (m1 + m2) ** 2) * (M ** 2 - (m1 - m2) ** 2)) / (2 * M);
  let c1 = 0, c2 = 0;
  for (let i = 0; i < N; i++) {
    const [p1, p2] = twoBodyDecay(M, m1, m2, rand);
    if (i < 100) {
      close(Math.hypot(p1[1], p1[2], p1[3]), k, 1e-13, '|p|');
      close(mass(p1), m1, 1e-12); close(mass(p2), m2, 1e-12);
      add4(p1, p2).forEach((v, j) => close(v, j ? 0 : M, 1e-13, 'total'));
    }
    const c = p1[3] / k; c1 += c / N; c2 += c * c / N;
  }
  close(c1, 0, 4 / Math.sqrt(3 * N), '⟨cos θ⟩');
  close(c2, 1 / 3, 4 * Math.sqrt(4 / 45 / N), '⟨cos² θ⟩');
  assert.throws(() => twoBodyDecay(1, 0.6, 0.5, rand), /twoBodyDecay/);
});

test('rambo massless: conservation, masses, constant weight = phase-space volume', () => {
  const rand = rng(4), sqrtS = 91.2;
  for (const n of [2, 3, 5, 8]) {
    for (let t = 0; t < 50; t++) {
      const { momenta, weight } = rambo(n, sqrtS, rand);
      assert.equal(momenta.length, n);
      add4(...momenta).forEach((v, j) => close(v, j ? 0 : sqrtS, 1e-12 * sqrtS, `n=${n} total`));
      momenta.forEach(p => close(dot4(p, p) / (p[0] * p[0]), 0, 1e-12, 'massless'));
      if (n === 2) close(weight, 1 / (8 * Math.PI), 1e-15, 'Φ₂ = 1/(8π)');
      if (n === 3) close(weight / (sqrtS ** 2 / (256 * Math.PI ** 3)), 1, 1e-12, 'Φ₃ = s/(256π³)');
    }
  }
});

test('rambo n=3 massless: x = 2E/√s has density 2x on [0,1] (uniform Dalitz plot)', () => {
  // Uniform Dalitz plot in (x1, x2) on the triangle x1, x2 ≤ 1 ≤ x1 + x2  ⇒  marginal ∝ x1, CDF x².
  const rand = rng(5), xs = [], ys = [];
  for (let i = 0; i < 5000; i++) {
    const { momenta } = rambo(3, 10, rand);
    xs.push(2 * momenta[0][0] / 10); ys.push(2 * momenta[2][0] / 10);
  }
  assert.ok(ksTest(xs, x => Math.min(1, Math.max(0, x)) ** 2).p > 0.01, 'x₁');
  assert.ok(ksTest(ys, x => Math.min(1, Math.max(0, x)) ** 2).p > 0.01, 'x₃');
  assert.ok(ksTest(xs, x => Math.min(1, Math.max(0, x))).p < 1e-6, 'a flat marginal is rejected');
});

test('rambo massive: conservation, masses, mean weight = massive phase-space volume', () => {
  const rand = rng(6);
  // n = 2: weight is exactly |k|/(4π√s)
  const M = 3, m1 = 1, m2 = 0.5, k = Math.sqrt((M ** 2 - (m1 + m2) ** 2) * (M ** 2 - (m1 - m2) ** 2)) / (2 * M);
  for (let t = 0; t < 20; t++) {
    const { momenta, weight } = rambo(2, M, rand, [m1, m2]);
    close(weight, k / (4 * Math.PI * M), 1e-13, 'Φ₂ massive');
    close(mass(momenta[0]), m1, 1e-12); close(mass(momenta[1]), m2, 1e-12);
  }
  // n = 3: ⟨weight⟩ = Dalitz area / (128 π³ s)
  const S = 5, ms = [0.5, 1, 1.5], N = 40000, ws = [];
  for (let i = 0; i < N; i++) {
    const { momenta, weight } = rambo(3, S, rand, ms);
    if (i < 200) {
      add4(...momenta).forEach((v, j) => close(v, j ? 0 : S, 1e-12, 'total'));
      momenta.forEach((p, j) => close(mass(p), ms[j], 1e-10, 'mass'));
    }
    ws.push(weight);
  }
  const [a, b, c] = ms;
  const width = s12 => {
    const r = Math.sqrt(s12), E2 = (s12 - a * a + b * b) / (2 * r), E3 = (S * S - s12 - c * c) / (2 * r);
    return 4 * Math.sqrt(Math.max(0, E2 * E2 - b * b)) * Math.sqrt(Math.max(0, E3 * E3 - c * c));
  };
  const area = adaptiveSimpson(width, (a + b) ** 2, (S - c) ** 2, 1e-12);
  const { mean, sem } = meanSem(ws), exact = area / (128 * Math.PI ** 3 * S * S);
  close(mean, exact, 4 * sem, 'Φ₃ massive');
  assert.ok(sem / exact < 0.01, 'test is sharp');
  assert.throws(() => rambo(3, 2, rand, [1, 1, 1]), /rambo/);
});

test('acceptReject: samples the target, counts tries, throws on a wrong bound', () => {
  const rand = rng(7), xs = [];
  let tries = 0;
  for (let i = 0; i < 4000; i++) { const r = acceptReject(() => rand(), x => 2 * x, 2, rand); xs.push(r.value); tries += r.tries; }
  assert.ok(ksTest(xs, x => x * x).p > 0.01, 'density 2x');
  close(tries / 4000, 2, 0.1, 'mean tries = wmax / ⟨w⟩');
  assert.throws(() => { for (let i = 0; i < 100; i++) acceptReject(() => rand(), x => 2 * x, 1.5, rand); }, /exceeds wmax = 1\.5/);
  assert.throws(() => acceptReject(() => 0, () => 0, 1, rand, { maxTries: 10 }), /no acceptance/);
});
