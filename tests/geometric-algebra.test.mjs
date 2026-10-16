// Geometric algebra checked against algebraic identities, rotation matrices (Rodrigues),
// relativistic velocity addition, the Thomas–Wigner rotation angle for perpendicular boosts, and
// the Lorentz invariants of the electromagnetic field (Doran & Lasenby 2003, ch. 4–7).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  cliffordAlgebra, CL3, STA, rotor3, rotor3Between, rotateVector3, rotorToMatrix3, cross3,
  boostRotor, fourVector, fourVelocity, spatialVelocity, fieldBivector, fieldComponents, fieldInvariants,
} from '../core/geometric-algebra.mjs';
import { rng } from '../core/numerics.mjs';

const close = (a, b, tol = 1e-13) => assert.ok(Math.abs(a - b) <= tol, `${a} vs ${b} (tolerance ${tol})`);
const vclose = (a, b, tol = 1e-13) => { assert.equal(a.length, b.length); Array.from(a).forEach((v, i) => close(v, b[i], tol)); };
const rand = rng(17);
const randomMV = A => Float64Array.from({ length: A.size }, () => 2 * rand() - 1);
const random3 = () => [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1];

test('Clifford algebra axioms: signature, associativity, ab = a·b + a∧b, reverse, inverse', () => {
  for (const [p, q] of [[3, 0], [1, 3], [2, 2], [0, 3], [4, 1]]) {
    const A = cliffordAlgebra(p, q);
    for (let i = 0; i < A.n; i++) { const e = A.blade(1 << i), sq = A.mul(e, e); vclose(sq, A.scalar(i < p ? 1 : -1), 0); }
    const a = randomMV(A), b = randomMV(A), c = randomMV(A);
    vclose(A.mul(A.mul(a, b), c), A.mul(a, A.mul(b, c)), 1e-12);
    vclose(A.mul(a, A.add(b, c)), A.add(A.mul(a, b), A.mul(a, c)), 1e-12);
    vclose(A.reverse(A.mul(a, b)), A.mul(A.reverse(b), A.reverse(a)), 1e-12);
    vclose(A.mul(a, A.inverse(a)), A.scalar(1), 1e-9);
    vclose(A.conjugate(a), A.reverse(A.involute(a)), 0);
    const u = A.grade(a, 1), v = A.grade(b, 1);
    vclose(A.mul(u, v), A.add(A.scalar(A.scalarProduct(u, v)), A.wedge(u, v)), 1e-13);
    vclose(A.wedge(u, v), A.neg(A.wedge(v, u)), 1e-15);
    vclose(A.wedge(u, u), A.zero(), 1e-15);
    // left contraction of a vector into a bivector: u⌋(v∧w) = (u·v)w − (u·w)v
    const w = A.grade(c, 1);
    vclose(A.inner(u, A.wedge(v, w)), A.sub(A.scale(w, A.scalarProduct(u, v)), A.scale(v, A.scalarProduct(u, w))), 1e-13);
  }
  assert.equal(CL3.toString(CL3.add(CL3.blade('e12', 2), CL3.scalar(1))), '1 + 2e12');
  assert.equal(STA.names[0b0011], 'e01');
  assert.throws(() => cliffordAlgebra(5, 4), /p \+ q ≤ 8/);
  assert.throws(() => CL3.inverse(CL3.add(CL3.scalar(1), CL3.basis(1))), /not invertible/);   // 1 + e1 is a null-divisor
});

test('3D rotors reproduce Rodrigues rotation matrices, compose like matrices, and give the cross product', () => {
  for (let k = 0; k < 20; k++) {
    const axis = random3(), angle = 6 * rand() - 3, n = Math.hypot(...axis), u = axis.map(x => x / n);
    const c = Math.cos(angle), s = Math.sin(angle);
    const M = [0, 1, 2].map(i => [0, 1, 2].map(j => c * (i === j) + (1 - c) * u[i] * u[j] + s * [[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]][i][j]));
    const R = rotor3(axis, angle);
    close(CL3.norm2(R), 1, 1e-15);
    rotorToMatrix3(R).forEach((row, i) => vclose(row, M[i], 1e-14));
    const v = random3();
    vclose(rotateVector3(R, v), M.map(row => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]), 1e-14);
    // composition: (R2 R1) ↔ M2 M1
    const R2 = rotor3(random3(), 2 * rand()), M2 = rotorToMatrix3(R2), M12 = rotorToMatrix3(CL3.mul(R2, R));
    M12.forEach((row, i) => vclose(row, [0, 1, 2].map(j => M2[i][0] * M[0][j] + M2[i][1] * M[1][j] + M2[i][2] * M[2][j]), 1e-14));
    // the double cover: rotation by 2π is −1
    vclose(rotor3(axis, 2 * Math.PI), CL3.scalar(-1), 1e-15);
    const a = random3(), b = random3();
    vclose(cross3(a, b), [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]], 1e-15);
    const Rab = rotor3Between(a, b), bn = Math.hypot(...b), an = Math.hypot(...a);
    vclose(rotateVector3(Rab, a).map(x => x / an), b.map(x => x / bn), 1e-14);
  }
  vclose(rotateVector3(rotor3([0, 0, 1], Math.PI / 2), [1, 0, 0]), [0, 1, 0], 1e-15);   // right-handed
});

test('boost rotors: four-velocity, collinear rapidity addition and Einstein velocity addition', () => {
  const v = [0.6, 0, 0];
  vclose(STA.sandwich(boostRotor(v), STA.basis(0)), fourVelocity(v), 1e-15);
  // collinear boosts: (u + v)/(1 + uv)
  const R = STA.mul(boostRotor([0.5, 0, 0]), boostRotor([0.7, 0, 0]));
  vclose(R, boostRotor([1.2 / 1.35, 0, 0]), 1e-15);
  for (let k = 0; k < 10; k++) {
    const u = random3().map(x => 0.55 * x), w = random3().map(x => 0.55 * x);
    // boost by u, then by w (lab frame): velocity w ⊕ u
    const R21 = STA.mul(boostRotor(w), boostRotor(u)), out = spatialVelocity(STA.sandwich(R21, STA.basis(0)));
    const gw = 1 / Math.sqrt(1 - (w[0] ** 2 + w[1] ** 2 + w[2] ** 2)), wu = w[0] * u[0] + w[1] * u[1] + w[2] * u[2];
    const expected = [0, 1, 2].map(i => (w[i] + u[i] / gw + gw / (1 + gw) * wu * w[i]) / (1 + wu));
    vclose(out, expected, 1e-14);
    // boosts preserve the interval x² = t² − r²
    const x = fourVector([1.3, ...random3()]), xp = STA.sandwich(R21, x);
    close(STA.scalarProduct(xp, xp), STA.scalarProduct(x, x), 1e-13);
  }
});

test('perpendicular boosts compose to a boost times the Thomas–Wigner rotation, cos θ = (γ₁+γ₂)/(1+γ₁γ₂)', () => {
  for (const [a, b] of [[0.6, 0.8], [0.3, 0.9], [0.95, 0.5]]) {
    const R21 = STA.mul(boostRotor([0, b, 0]), boostRotor([a, 0, 0]));
    const w = spatialVelocity(STA.sandwich(R21, STA.basis(0)));
    const W = STA.mul(STA.reverse(boostRotor(w)), R21);        // what is left after the pure boost
    vclose(STA.sandwich(W, STA.basis(0)), STA.basis(0), 1e-14);  // a spatial rotation
    for (let i = 0; i < STA.size; i++) if (STA.gradeOf[i] === 2 && (i & 1)) close(W[i], 0, 1e-14);
    const g1 = 1 / Math.sqrt(1 - a * a), g2 = 1 / Math.sqrt(1 - b * b);
    close(2 * W[0] ** 2 - 1, (g1 + g2) / (1 + g1 * g2), 1e-14);  // cos θ = 2cos²(θ/2) − 1
    assert.ok(Math.abs(W[0]) < 1 - 1e-6);                         // non-trivial rotation
  }
});

test('Maxwell bivector: F² = E² − B² + 2I E·B, invariants under boosts, moving-charge field B = v × E', () => {
  for (let k = 0; k < 10; k++) {
    const E = random3(), B = random3(), F = fieldBivector(E, B);
    const c = fieldComponents(F);
    vclose(c.E, E, 1e-15); vclose(c.B, B, 1e-15);
    const inv = fieldInvariants(F), E2 = E.reduce((s, x) => s + x * x, 0), B2 = B.reduce((s, x) => s + x * x, 0);
    close(inv.EsqMinusBsq, E2 - B2, 1e-14); close(inv.EdotB, E[0] * B[0] + E[1] * B[1] + E[2] * B[2], 1e-14);
    assert.deepEqual(STA.grades(STA.mul(F, F), 1e-14).filter(g => g !== 0 && g !== 4), []);
    const R = STA.mul(boostRotor(random3().map(x => 0.5 * x)), STA.exp(STA.scale(STA.blade('e12'), rand())));
    const inv2 = fieldInvariants(STA.sandwich(R, F));
    close(inv2.EsqMinusBsq, inv.EsqMinusBsq, 1e-12); close(inv2.EdotB, inv.EdotB, 1e-12);
  }
  // pure E field (charge at rest) actively boosted to velocity v: E∥ unchanged, E⊥ × γ, B = v × E
  const v = [0.5, 0.2, -0.1], E = [0.3, -0.4, 1.2], { E: E2, B: B2 } = fieldComponents(STA.sandwich(boostRotor(v), fieldBivector(E, [0, 0, 0])));
  vclose(B2, cross3(v, E2), 1e-14);
  const v2 = v.reduce((s, x) => s + x * x, 0), g = 1 / Math.sqrt(1 - v2), par = (E[0] * v[0] + E[1] * v[1] + E[2] * v[2]) / v2;
  vclose(E2, E.map((x, i) => par * v[i] + g * (x - par * v[i])), 1e-14);
  // exp of a non-simple bivector (series path) agrees with the product of commuting simple factors
  const Fs = STA.add(STA.blade('e01', 0.4), STA.blade('e23', 0.7));     // e01 and e23 commute
  vclose(STA.exp(Fs), STA.mul(STA.exp(STA.blade('e01', 0.4)), STA.exp(STA.blade('e23', 0.7))), 1e-14);
});
