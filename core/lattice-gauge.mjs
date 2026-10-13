// Lattice gauge theory: pure SU(2) Yang–Mills on a periodic hypercubic lattice in
// d = 2, 3 or 4 dimensions, Wilson action S = β Σ_P (1 − ½ Tr U_P).
//
// This is a first-principles non-abelian gauge theory, the same framework as lattice QCD
// (with SU(2) instead of SU(3), no dynamical quarks, small volumes). Confinement shows up
// as an area law for Wilson loops, W(R,T) ~ e^{−σ R T}, with string tension σ measured
// from Creutz ratios. Updates: Creutz heat bath (Kennedy–Pendleton) + over-relaxation.
//
// SU(2) elements are unit quaternions [a0, a1, a2, a3] ↔ a0·1 + i a·σ.

import { rng as makeRng } from './numerics.mjs';

// ---------------------------------------------------------------- SU(2) algebra (quaternions)

/** out = p·q */
export function qmul(p, q, out = new Float64Array(4)) {
  const [a0, a1, a2, a3] = p, [b0, b1, b2, b3] = q;
  out[0] = a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3;
  out[1] = a0 * b1 + a1 * b0 - a2 * b3 + a3 * b2;
  out[2] = a0 * b2 + a2 * b0 - a3 * b1 + a1 * b3;
  out[3] = a0 * b3 + a3 * b0 - a1 * b2 + a2 * b1;
  return out;
}
export const qdag = (p, out = new Float64Array(4)) => { out[0] = p[0]; out[1] = -p[1]; out[2] = -p[2]; out[3] = -p[3]; return out; };
/** ½ Tr of a quaternion matrix = its scalar part. */
export const halfTr = p => p[0];

// ---------------------------------------------------------------- lattice

export function createLattice({ L = 6, d = 4, beta = 2.3, seed = 1, start = 'cold' } = {}) {
  if (!(Number.isFinite(beta) && beta >= 0)) throw new Error('createLattice: beta must be finite and non-negative');
  const dims = Array.isArray(L) ? L : Array(d).fill(L);
  const D = dims.length;
  const V = dims.reduce((a, b) => a * b, 1);
  const U = new Float64Array(V * D * 4);
  const rand = makeRng(seed);
  for (let i = 0; i < V * D; i++) {
    if (start === 'hot') {
      const q = randomSU2(rand);
      U.set(q, 4 * i);
    } else U[4 * i] = 1;
  }
  // neighbour tables: up[site*D+mu], dn[site*D+mu]
  const up = new Int32Array(V * D), dn = new Int32Array(V * D);
  const stride = dims.map((_, mu) => dims.slice(0, mu).reduce((a, b) => a * b, 1));
  for (let s = 0; s < V; s++) {
    for (let mu = 0; mu < D; mu++) {
      const x = Math.floor(s / stride[mu]) % dims[mu];
      up[s * D + mu] = s + ((x + 1) % dims[mu] - x) * stride[mu];
      dn[s * D + mu] = s + ((x - 1 + dims[mu]) % dims[mu] - x) * stride[mu];
    }
  }
  return { dims, D, V, U, up, dn, beta, rand, sweeps: 0 };
}

function randomSU2(rand) {
  let a, n;
  do { a = [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]; n = Math.hypot(...a); } while (n > 1 || n === 0);
  return a.map(v => v / n);
}

const link = (lat, s, mu) => lat.U.subarray(4 * (s * lat.D + mu), 4 * (s * lat.D + mu) + 4);

/** Sum of the 2(D−1) staples around link (s, μ), so that the local action is −(β/2) Tr(U·A). */
export function staple(lat, s, mu, A = new Float64Array(4)) {
  A.fill(0);
  const { D, up, dn } = lat;
  const t1 = new Float64Array(4), t2 = new Float64Array(4), t3 = new Float64Array(4);
  const sp = up[s * D + mu];
  for (let nu = 0; nu < D; nu++) {
    if (nu === mu) continue;
    // forward: U_ν(s+μ) U_μ(s+ν)† U_ν(s)†
    qmul(link(lat, sp, nu), qdag(link(lat, up[s * D + nu], mu), t1), t2);
    qmul(t2, qdag(link(lat, s, nu), t1), t3);
    for (let k = 0; k < 4; k++) A[k] += t3[k];
    // backward: U_ν(s+μ−ν)† U_μ(s−ν)† U_ν(s−ν)
    const sm = dn[s * D + nu], spm = dn[sp * D + nu];
    qmul(qdag(link(lat, spm, nu), t1), qdag(link(lat, sm, mu), t3), t2);
    qmul(t2, link(lat, sm, nu), t3);
    for (let k = 0; k < 4; k++) A[k] += t3[k];
  }
  return A;
}

/** Kennedy–Pendleton heat-bath sample of a0 with density ∝ √(1−a0²) e^{k a0}. */
function samplea0(k, rand) {
  if (!(Number.isFinite(k) && k >= 0)) throw new Error('heat bath: coupling must be finite and non-negative');
  // Kennedy–Pendleton rejection becomes inefficient as k → 0. Propose from
  // Haar measure instead: accept with exp(k(a0−1)), including the exact k = 0 limit.
  if (k < 1) {
    for (;;) {
      const a0 = randomSU2(rand)[0];
      if (rand() < Math.exp(k * (a0 - 1))) return a0;
    }
  }
  for (;;) {
    const r1 = 1 - rand(), r2 = 1 - rand(), r3 = 1 - rand();
    const lambda2 = -(Math.log(r1) + Math.cos(2 * Math.PI * r2) ** 2 * Math.log(r3)) / (2 * k);
    if (rand() ** 2 <= 1 - lambda2) return 1 - 2 * lambda2;
  }
}

/** One heat-bath update of every link, then `overrelax` over-relaxation sweeps. */
export function sweep(lat, { overrelax = 2, heatbath = true } = {}) {
  const { V, D, beta, rand } = lat;
  if (!(Number.isFinite(beta) && beta >= 0)) throw new Error('sweep: beta must be finite and non-negative');
  const A = new Float64Array(4), Vn = new Float64Array(4), X = new Float64Array(4), tmp = new Float64Array(4);
  if (heatbath) {for (let s = 0; s < V; s++) {for (let mu = 0; mu < D; mu++) {
    staple(lat, s, mu, A);
    const k = Math.hypot(A[0], A[1], A[2], A[3]);
    if (k < 1e-12) { link(lat, s, mu).set(randomSU2(rand)); continue; }
    for (let i = 0; i < 4; i++) Vn[i] = A[i] / k;            // A = k·V with V ∈ SU(2)
    // sample X with weight e^{β k ½Tr X}, then U = X V†
    const a0 = samplea0(beta * k, rand);
    const r = Math.sqrt(Math.max(0, 1 - a0 * a0));
    const ct = 2 * rand() - 1, st = Math.sqrt(1 - ct * ct), ph = 2 * Math.PI * rand();
    X[0] = a0; X[1] = r * st * Math.cos(ph); X[2] = r * st * Math.sin(ph); X[3] = r * ct;
    qmul(X, qdag(Vn, tmp), link(lat, s, mu));
  }}}
  for (let o = 0; o < overrelax; o++) {
    for (let s = 0; s < V; s++) {for (let mu = 0; mu < D; mu++) {
      staple(lat, s, mu, A);
      const k = Math.hypot(A[0], A[1], A[2], A[3]);
      if (k < 1e-12) continue;
      for (let i = 0; i < 4; i++) Vn[i] = A[i] / k;
      // U → V† U† V† : leaves the local action unchanged (microcanonical), decorrelates fast
      const U = link(lat, s, mu);
      qmul(qdag(Vn, tmp), qdag(U, X), A);
      qmul(A, qdag(Vn, tmp), U);
    }}
  }
  // re-unitarise against round-off
  for (let i = 0; i < V * D; i++) {
    const o = 4 * i, n = Math.hypot(lat.U[o], lat.U[o + 1], lat.U[o + 2], lat.U[o + 3]);
    for (let k = 0; k < 4; k++) lat.U[o + k] /= n;
  }
  lat.sweeps++;
}

/** Random local gauge transformation U_μ(s) → g(s) U_μ(s) g(s+μ)†: physics must not change. */
export function gaugeTransform(lat, rand = lat.rand) {
  const g = Array.from({ length: lat.V }, () => randomSU2(rand));
  const t = new Float64Array(4), d = new Float64Array(4);
  for (let s = 0; s < lat.V; s++) {for (let mu = 0; mu < lat.D; mu++) {
    const U = link(lat, s, mu);
    qmul(g[s], U, t);
    qmul(t, qdag(g[lat.up[s * lat.D + mu]], d), U);
  }}
}

/** Total Wilson action S/β = Σ_P (1 − ½ Tr U_P). */
export function action(lat) {
  const nP = lat.V * lat.D * (lat.D - 1) / 2;
  return nP * (1 - plaquette(lat));
}

/** Average plaquette ⟨½ Tr U_P⟩. */
export function plaquette(lat) {
  const { V, D, up } = lat;
  const t1 = new Float64Array(4), t2 = new Float64Array(4), t3 = new Float64Array(4);
  let sum = 0, n = 0;
  for (let s = 0; s < V; s++) {for (let mu = 0; mu < D; mu++) {for (let nu = mu + 1; nu < D; nu++) {
    qmul(link(lat, s, mu), link(lat, up[s * D + mu], nu), t1);
    qmul(t1, qdag(link(lat, up[s * D + nu], mu), t2), t3);
    qmul(t3, qdag(link(lat, s, nu), t2), t1);
    sum += t1[0]; n++;
  }}}
  return sum / n;
}

/**
 * Average R×T Wilson loop ⟨½ Tr W⟩ over all sites and all planes μ < ν, with R along μ
 * (the lower index) and T along ν. On an anisotropic lattice W(R,T) ≠ W(T,R) in general.
 * Loops must fit in the lattice without wrapping: R, T < extent of every direction.
 */
export function wilsonLoop(lat, R, T) {
  if (R >= Math.min(...lat.dims) || T >= Math.min(...lat.dims)) throw new Error(`wilsonLoop: R, T must be < lattice extent ${Math.min(...lat.dims)}`);
  const { V, D, up } = lat;
  let sum = 0, n = 0;
  const P = new Float64Array(4), t = new Float64Array(4), d = new Float64Array(4);
  for (let mu = 0; mu < D; mu++) {for (let nu = 0; nu < D; nu++) {
    if (mu >= nu) continue;                                  // each plane once
    for (let s0 = 0; s0 < V; s0++) {
      P.fill(0); P[0] = 1;
      let s = s0;
      const go = (dir, steps, back) => {
        for (let i = 0; i < steps; i++) {
          if (!back) { qmul(P, link(lat, s, dir), t); P.set(t); s = up[s * D + dir]; }
          else { s = lat.dn[s * D + dir]; qmul(P, qdag(link(lat, s, dir), d), t); P.set(t); }
        }
      };
      go(mu, R, false); go(nu, T, false); go(mu, R, true); go(nu, T, true);
      sum += P[0]; n++;
    }
  }}
  return sum / n;
}

/** Creutz ratio χ(R,T) = −ln[ W(R,T) W(R−1,T−1) / (W(R,T−1) W(R−1,T)) ] → σa² for large loops. */
export const creutzRatio = (W, R, T) => -Math.log(W(R, T) * W(R - 1, T - 1) / (W(R, T - 1) * W(R - 1, T)));
