import { test } from 'node:test';
import assert from 'node:assert/strict';
import { cmat, identity, fromReal, cadd, cscale, cmul, dagger, kron, trace, expectation, spinOps, embed, eigHermitian, propagator, evolve, thermalState, lindblad, checkDensity } from '../core/quantum.mjs';
import { rng } from '../core/numerics.mjs';

const maxDiff = (A, B) => { let m = 0; for (let i = 0; i < A.re.length; i++) m = Math.max(m, Math.hypot(A.re[i] - B.re[i], A.im[i] - B.im[i])); return m; };
const randomHermitian = (n, seed) => {
  const r = rng(seed), A = cmat(n);
  for (let i = 0; i < n; i++) {for (let j = i; j < n; j++) {
    const a = r() - 0.5, b = i === j ? 0 : r() - 0.5;
    A.re[i * n + j] = A.re[j * n + i] = a; A.im[i * n + j] = b; A.im[j * n + i] = -b;
  }}
  return A;
};

test('spinOps: commutation [Sx, Sy] = i Sz and S² = j(j+1) for j = ½, 1, 3/2', () => {
  for (const j of [0.5, 1, 1.5]) {
    const { x, y, z, dim } = spinOps(j);
    const comm = cadd(cmul(x, y), cmul(y, x), -1);
    assert.ok(maxDiff(comm, cscale(z, [0, 1])) < 1e-14, `j = ${j}`);
    const S2 = cadd(cadd(cmul(x, x), cmul(y, y)), cmul(z, z));
    assert.ok(maxDiff(S2, cscale(identity(dim), j * (j + 1))) < 1e-13);
  }
});

test('kron and embed: dimensions, (A⊗B)(C⊗D) = AC⊗BD, operators on different subsystems commute', () => {
  const a = spinOps(0.5), b = spinOps(1);
  const left = cmul(kron(a.x, b.y), kron(a.z, b.z)), right = kron(cmul(a.x, a.z), cmul(b.y, b.z));
  assert.ok(maxDiff(left, right) < 1e-14);
  const dims = [2, 2, 3], S1 = embed(a.x, 0, dims), I3 = embed(b.z, 2, dims);
  assert.equal(S1.n, 12);
  assert.ok(maxDiff(cmul(S1, I3), cmul(I3, S1)) < 1e-14);
});

test('eigHermitian: H V = V diag(E), V unitary, real symmetric limit, degenerate spectra', () => {
  for (const [n, seed] of [[2, 1], [5, 2], [12, 3]]) {
    const H = randomHermitian(n, seed), { values, vectors } = eigHermitian(H);
    const D = cmat(n); values.forEach((e, k) => { D.re[k * n + k] = e; });
    assert.ok(maxDiff(cmul(H, vectors), cmul(vectors, D)) < 1e-12, `n = ${n}`);
    assert.ok(maxDiff(cmul(dagger(vectors), vectors), identity(n)) < 1e-12);
    for (let k = 1; k < n; k++) assert.ok(values[k] >= values[k - 1]);
  }
  // two decoupled spin-½ in zero field: 4-fold degenerate 0; Heisenberg J S1·S2: triplet J/4 ×3, singlet −3J/4
  const s = spinOps(0.5), dot = cadd(cadd(cmul(embed(s.x, 0, [2, 2]), embed(s.x, 1, [2, 2])), cmul(embed(s.y, 0, [2, 2]), embed(s.y, 1, [2, 2]))), cmul(embed(s.z, 0, [2, 2]), embed(s.z, 1, [2, 2])));
  const e = eigHermitian(dot);
  assert.ok(Math.abs(e.values[0] + 0.75) < 1e-14 && e.values.slice(1).every(v => Math.abs(v - 0.25) < 1e-14));
  assert.ok(maxDiff(cmul(dagger(e.vectors), e.vectors), identity(4)) < 1e-13);
  assert.throws(() => eigHermitian(fromReal([[1, 2], [0, 1]])), /Hermitian/);
});

test('propagator: Larmor precession of a spin-½ (⟨Sx⟩ = ½cos ωt), unitary, composes in time', () => {
  const s = spinOps(0.5), w = 2.3, H = cscale(s.z, w), eig = eigHermitian(H);
  const plusX = fromReal([[0.5, 0.5], [0.5, 0.5]]);                 // |+x⟩⟨+x|
  for (const t of [0, 0.4, 1.7]) assert.ok(Math.abs(expectation(evolve(plusX, propagator(eig, t)), s.x) - 0.5 * Math.cos(w * t)) < 1e-14);
  const U1 = propagator(eig, 0.3), U2 = propagator(eig, 0.5), U3 = propagator(eig, 0.8);
  assert.ok(maxDiff(cmul(U2, U1), U3) < 1e-14);
  assert.ok(maxDiff(cmul(dagger(U3), U3), identity(2)) < 1e-14);
});

test('lindblad: pure dephasing decays the coherence as e^{−γt}·(phase), populations fixed; decay gives e^{−Γt}', () => {
  const s = spinOps(0.5), w = 1.1, g = 0.7;
  const P0 = fromReal([[1, 0], [0, 0]]), P1 = fromReal([[0, 0], [0, 1]]);
  const plusX = fromReal([[0.5, 0.5], [0.5, 0.5]]);
  // dephasing with L = |0⟩⟨0| at rate 2γ gives off-diagonal decay e^{−γt}
  const L = lindblad({ H: cscale(s.z, w), jumps: [{ L: P0, rate: 2 * g }], rho0: plusX });
  L.advanceTo(2, { exact: true });
  const r = L.rho, coh = Math.hypot(r.re[1], r.im[1]);
  assert.ok(Math.abs(coh - 0.5 * Math.exp(-g * 2)) < 1e-8, `${coh}`);
  assert.ok(Math.abs(r.re[0] - 0.5) < 1e-10);
  // amplitude decay |1⟩ → |0⟩ at rate Γ
  const down = fromReal([[0, 1], [0, 0]]), G = 0.9;
  const D = lindblad({ H: cmat(2), jumps: [{ L: down, rate: G }], rho0: P1 });
  D.advanceTo(1.5, { exact: true });
  assert.ok(Math.abs(D.population(1) - Math.exp(-G * 1.5)) < 1e-9);
  const c = checkDensity(D.rho);
  assert.ok(Math.abs(c.trace - 1) < 1e-10 && c.hermiticity < 1e-12 && c.minEigenvalue > -1e-10);
});

test('lindblad with up/down rates in detailed balance relaxes to the thermal state of H', () => {
  const E = 1.3, kT = 0.6, H = fromReal([[0, 0], [0, E]]);
  const down = fromReal([[0, 1], [0, 0]]), up = dagger(down), k = 0.8;
  const L = lindblad({ H, jumps: [{ L: down, rate: k }, { L: up, rate: k * Math.exp(-E / kT) }], rho0: fromReal([[0, 0], [0, 1]]) });
  L.advanceTo(40);
  const th = thermalState(H, kT);
  assert.ok(Math.abs(L.population(1) - th.re[3]) < 1e-8, `${L.population(1)} vs ${th.re[3]}`);
  assert.ok(Math.abs(trace(th).re - 1) < 1e-14);
});
