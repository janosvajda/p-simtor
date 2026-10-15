// Lattice gauge theory: pure SU(2) Yang–Mills on a periodic hypercubic lattice in
// d = 2, 3 or 4 dimensions, Wilson action S = β Σ_P (1 − ½ Tr U_P).
//
// This is a first-principles non-abelian gauge theory, the same framework as lattice QCD
// (with SU(2) instead of SU(3), no dynamical quarks, small volumes). Confinement shows up
// as an area law for Wilson loops, W(R,T) ~ e^{−σ R T}, with string tension σ measured
// from Creutz ratios. Updates: Creutz heat bath (Kennedy–Pendleton) + over-relaxation.
//
// SU(2) elements are unit quaternions [a0, a1, a2, a3] ↔ a0·1 + i a·σ.
//
//   createLattice, sweep                        sampling; sweep(lat, { sites: [s0, s1] }) updates a range (spread a
//                                               slow sweep over frames), { cool: true } minimises the action
//   topology(lat, { improved })                 clover field strength: action density s(x), charge density q(x), Q, S (4D)
//   bpstInstanton({ L, rho, anti })             links of an exact BPST (anti-)instanton: a benchmark for topology
//   gaugeTransform(lat, rand | g)               random or given g(x); returns g. Links are gauge-DEPENDENT:
//                                               unfixed, each link is exactly Haar-distributed and ⟨U⟩ = 0
//                                               (Elitzur's theorem); only traces of closed loops are physical
//   centreTransform(lat, { dir, slice })        multiply the links crossing one slice by −1 (Z₂ centre):
//                                               action unchanged, Polyakov loops along dir flip sign
//   landauGauge(lat, { omega, tol })            lattice Landau gauge fixing (over-relaxed), returns { theta, g }
//   haarSU2, qexp, qlog                         Haar-random SU(2) element; SU(2) exponential and logarithm
//   plaquette, action, wilsonLoop, creutzRatio  averages (wilsonLoop: one loop, link by link)
//   wilsonLoops(lat, { Rmax, Tmax, planes })    every R×T loop at once from line products (~50× faster)
//   polyakovLoop(lat, { dir })                  ⟨½Tr P⟩ and its per-site values: finite-volume/deconfinement diagnostic
//   plaquetteField(lat)                         1 − ½Tr U_P per site and plane: local action density for views
//   apeSmear(lat, { alpha, steps, dirs })       smeared copy of the links (glueball / potential operators)
//   strongCouplingPlaquette(beta, d)            ⟨P⟩ ≈ u + 2(d−2)u⁵, u = I₂(β)/I₁(β): the check in any dimension
//   siteCoords, siteIndex                       site number ↔ coordinates (axis 0 fastest)
//   qmulInto                                    allocation-free quaternion product on flat arrays

import { rng as makeRng, besselI } from './numerics.mjs';

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
/**
 * out[oo..] = A·B on flat arrays without allocating, where A = a[ao..ao+3] (or its dagger if da)
 * and B = b[bo..bo+3] (or its dagger if db). out may not overlap a or b.
 */
export function qmulInto(out, oo, a, ao, b, bo, da = false, db = false) {
  const sa = da ? -1 : 1, sb = db ? -1 : 1;
  const a0 = a[ao], a1 = sa * a[ao + 1], a2 = sa * a[ao + 2], a3 = sa * a[ao + 3];
  const b0 = b[bo], b1 = sb * b[bo + 1], b2 = sb * b[bo + 2], b3 = sb * b[bo + 3];
  out[oo] = a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3;
  out[oo + 1] = a0 * b1 + a1 * b0 - a2 * b3 + a3 * b2;
  out[oo + 2] = a0 * b2 + a2 * b0 - a3 * b1 + a1 * b3;
  out[oo + 3] = a0 * b3 + a3 * b0 - a1 * b2 + a2 * b1;
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

/** Haar-random SU(2) element (uniform on S³). */
export function haarSU2(rand) { return Float64Array.from(randomSU2(rand)); }

/** exp(i θ n̂·σ) for a 3-vector v = θ n̂: the unit quaternion [cos θ, sin θ n̂]. */
export function qexp(v) {
  const th = Math.hypot(v[0], v[1], v[2]), s = th > 1e-300 ? Math.sin(th) / th : 1;
  return Float64Array.from([Math.cos(th), s * v[0], s * v[1], s * v[2]]);
}
/** Inverse of qexp: the 3-vector θ n̂ with θ ∈ [0, π] (the rotation-free branch). */
export function qlog(q) {
  const s = Math.hypot(q[1], q[2], q[3]), th = Math.atan2(s, q[0]), f = s > 1e-300 ? th / s : 1;
  return [f * q[1], f * q[2], f * q[3]];
}

function randomSU2(rand) {
  let a, n;
  do { a = [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]; n = Math.hypot(...a); } while (n > 1 || n === 0);
  return a.map(v => v / n);
}

const link = (lat, s, mu) => lat.U.subarray(4 * (s * lat.D + mu), 4 * (s * lat.D + mu) + 4);

/** Sum of the 2(D−1) staples around link (s, μ), so that the local action is −(β/2) Tr(U·A). */
export function staple(lat, s, mu, A = new Float64Array(4), dirs) {
  A.fill(0);
  const { D, up, dn, U } = lat, t = ST;
  const sp = up[s * D + mu];
  for (let nu = 0; nu < D; nu++) {
    if (nu === mu || (dirs && !dirs.includes(nu))) continue;
    // forward: U_ν(s+μ) U_μ(s+ν)† U_ν(s)†
    qmulInto(t, 0, U, 4 * (sp * D + nu), U, 4 * (up[s * D + nu] * D + mu), false, true);
    qmulInto(t, 4, t, 0, U, 4 * (s * D + nu), false, true);
    for (let k = 0; k < 4; k++) A[k] += t[4 + k];
    // backward: U_ν(s+μ−ν)† U_μ(s−ν)† U_ν(s−ν)
    const sm = dn[s * D + nu], spm = dn[sp * D + nu];
    qmulInto(t, 0, U, 4 * (spm * D + nu), U, 4 * (sm * D + mu), true, true);
    qmulInto(t, 4, t, 0, U, 4 * (sm * D + nu));
    for (let k = 0; k < 4; k++) A[k] += t[4 + k];
  }
  return A;
}
const ST = new Float64Array(8);

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

/**
 * One heat-bath update of every link, then `overrelax` over-relaxation sweeps. sites: [s0, s1)
 * restricts the update to a range of sites (a full sweep spread over several frames; the
 * sweep counter advances when the range ends at V). Each range runs heat bath then over-relaxation,
 * so ranged sweeps are a different but equally valid update order (each link update keeps the
 * Boltzmann weight); for cooling, ranges reproduce a full sweep exactly. cool: true replaces each link by the one that
 * minimises its local action (U = A†/|A|): Wilson cooling, a smoothing that removes ultraviolet
 * noise but also slowly shrinks and destroys instantons; disclose the number of cooling sweeps.
 */
export function sweep(lat, { overrelax = 2, heatbath = true, sites = [0, lat.V], cool = false } = {}) {
  const { D, beta, rand } = lat;
  const [s0, s1] = sites;
  if (!(Number.isFinite(beta) && beta >= 0)) throw new Error('sweep: beta must be finite and non-negative');
  if (!(s0 >= 0 && s1 <= lat.V && s0 <= s1)) throw new Error(`sweep: sites must be a range within [0, ${lat.V}]`);
  const A = new Float64Array(4), Vn = new Float64Array(4), X = new Float64Array(4), tmp = new Float64Array(4);
  if (cool) {
    for (let s = s0; s < s1; s++) {for (let mu = 0; mu < D; mu++) {
      staple(lat, s, mu, A);
      const k = Math.hypot(A[0], A[1], A[2], A[3]);
      if (k < 1e-12) continue;
      const U = link(lat, s, mu);
      U[0] = A[0] / k; U[1] = -A[1] / k; U[2] = -A[2] / k; U[3] = -A[3] / k;
    }}
    if (s1 === lat.V) lat.sweeps++;
    return;
  }
  if (heatbath) {for (let s = s0; s < s1; s++) {for (let mu = 0; mu < D; mu++) {
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
    for (let s = s0; s < s1; s++) {for (let mu = 0; mu < D; mu++) {
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
  for (let i = s0 * D; i < s1 * D; i++) {
    const o = 4 * i, n = Math.hypot(lat.U[o], lat.U[o + 1], lat.U[o + 2], lat.U[o + 3]);
    for (let k = 0; k < 4; k++) lat.U[o + k] /= n;
  }
  if (s1 === lat.V) lat.sweeps++;
}

/**
 * Local gauge transformation U_μ(s) → g(s) U_μ(s) g(s+μ)†: physics must not change. Pass a rand
 * function (default lat.rand; pass your own to leave the Markov chain untouched) for a random g,
 * or an array of V quaternions to apply a given g. Returns g.
 */
export function gaugeTransform(lat, randOrG = lat.rand) {
  const g = typeof randOrG === 'function' ? Array.from({ length: lat.V }, () => randomSU2(randOrG)) : randOrG;
  if (g.length !== lat.V) throw new Error(`gaugeTransform: need ${lat.V} SU(2) elements`);
  const t = new Float64Array(4), d = new Float64Array(4);
  for (let s = 0; s < lat.V; s++) {for (let mu = 0; mu < lat.D; mu++) {
    const U = link(lat, s, mu);
    qmul(g[s], U, t);
    qmul(t, qdag(g[lat.up[s * lat.D + mu]], d), U);
  }}
  return g;
}

/**
 * Centre transformation: multiply every link along `dir` that leaves the slice x_dir = slice by
 * −1 (the Z₂ centre of SU(2)). Plaquettes and the action are unchanged; Polyakov loops along dir
 * change sign. In the confined phase ⟨P⟩ = 0 respects this symmetry; deconfinement breaks it.
 */
export function centreTransform(lat, { dir = lat.D - 1, slice = 0 } = {}) {
  for (let s = 0; s < lat.V; s++) {
    if (siteCoords(lat, s)[dir] !== slice) continue;
    const o = 4 * (s * lat.D + dir);
    for (let k = 0; k < 4; k++) lat.U[o + k] = -lat.U[o + k];
  }
}

/**
 * Lattice Landau gauge: maximise F = (1/(V·D)) Σ ½Tr U_μ(x) by site-wise gauge transformations,
 * with over-relaxation exponent omega (1 = plain, ~1.7 faster). Stops when the gauge-condition
 * violation θ = (1/V) Σ_x |Σ_μ (A_μ(x) − A_μ(x−μ))|² < tol, where A_μ is the vector part of U_μ.
 * Changes lat in place (an exactly gauge-equivalent configuration) and returns
 * { theta, F, iterations, g } with the accumulated g. The result is one Gribov copy, not unique.
 */
export function landauGauge(lat, { omega = 1.7, tol = 1e-12, maxIterations = 5000 } = {}) {
  const { V, D, U, dn } = lat;
  const g = Array.from({ length: V }, () => Float64Array.from([1, 0, 0, 0]));
  const K = new Float64Array(4), h = new Float64Array(4), t = new Float64Array(4);
  const theta = () => {
    let sum = 0;
    for (let s = 0; s < V; s++) {
      const d = [0, 0, 0];
      for (let mu = 0; mu < D; mu++) {
        const a = 4 * (s * D + mu), b = 4 * (dn[s * D + mu] * D + mu);
        for (let k = 0; k < 3; k++) d[k] += U[a + 1 + k] - U[b + 1 + k];
      }
      sum += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    }
    return sum / V;
  };
  const F = () => { let f = 0; for (let i = 0; i < V * D; i++) f += U[4 * i]; return f / (V * D); };
  let it = 0, th = theta();
  for (; it < maxIterations && th > tol; it++) {
    for (let s = 0; s < V; s++) {
      // K = Σ_μ [U_μ(s) + U_μ(s−μ)†]; the best h maximises ½Tr(h K): h = K†/|K|
      K.fill(0);
      for (let mu = 0; mu < D; mu++) {
        const a = 4 * (s * D + mu), b = 4 * (dn[s * D + mu] * D + mu);
        K[0] += U[a] + U[b];
        for (let k = 1; k < 4; k++) K[k] += U[a + k] - U[b + k];
      }
      const n = Math.hypot(K[0], K[1], K[2], K[3]);
      if (!(n > 0)) continue;
      h[0] = K[0] / n; h[1] = -K[1] / n; h[2] = -K[2] / n; h[3] = -K[3] / n;
      if (omega !== 1) h.set(qexp(qlog(h).map(v => v * omega)));
      // apply h at s: U_μ(s) → h U_μ(s), U_μ(s−μ) → U_μ(s−μ) h†
      for (let mu = 0; mu < D; mu++) {
        const a = 4 * (s * D + mu), b = 4 * (dn[s * D + mu] * D + mu);
        qmulInto(t, 0, h, 0, U, a); U.set(t, a);
        qmulInto(t, 0, U, b, h, 0, false, true); U.set(t, b);
      }
      qmulInto(t, 0, h, 0, g[s], 0); g[s].set(t);
    }
    th = theta();
  }
  return { theta: th, F: F(), iterations: it, g, converged: th <= tol };
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

// ---------------------------------------------------------------- fast measurements and fields

/** Coordinates of site s (axis 0 fastest), and the inverse. */
export function siteCoords(lat, site) {
  const x = [];
  let s = site;
  for (const n of lat.dims) { x.push(s % n); s = Math.floor(s / n); }
  return x;
}
export function siteIndex(lat, x) {
  let s = 0, stride = 1;
  lat.dims.forEach((n, mu) => { s += ((x[mu] % n) + n) % n * stride; stride *= n; });
  return s;
}

/** Products of `len` consecutive links along μ from every site: out[(s·len + r−1)·4..] = U_μ(s)…U_μ(s+(r−1)μ). */
function lineProducts(lat, mu, len) {
  const { V, D, U, up } = lat, out = new Float64Array(V * len * 4);
  for (let s = 0; s < V; s++) {
    out.set(U.subarray(4 * (s * D + mu), 4 * (s * D + mu) + 4), s * len * 4);
    let x = up[s * D + mu];
    for (let r = 2; r <= len; r++) { qmulInto(out, (s * len + r - 1) * 4, out, (s * len + r - 2) * 4, U, 4 * (x * D + mu)); x = up[x * D + mu]; }
  }
  return out;
}
/** shift[s·(n+1) + r] = s + rμ. */
function shifts(lat, mu, n) {
  const { V, D, up } = lat, out = new Int32Array(V * (n + 1));
  for (let s = 0; s < V; s++) { let x = s; for (let r = 0; r <= n; r++) { out[s * (n + 1) + r] = x; x = up[x * D + mu]; } }
  return out;
}

/**
 * Every Wilson loop ⟨½Tr W(R,T)⟩ for 1 ≤ R ≤ Rmax, 1 ≤ T ≤ Tmax in one pass, averaged over sites and
 * planes with the same convention as wilsonLoop (planes μ < ν, R along μ, T along ν), or over the
 * given planes [[μ, ν], …] (e.g. [[0, 3], [1, 3], [2, 3]] for spatial R and Euclidean time T).
 * Uses ½Tr(XY†) = X·Y for SU(2): two quaternion products per loop. Returns W(R, T) as a function
 * with W.values (Float64Array, index (R−1)·Tmax + T−1: the vector to resample together), W.Rmax, W.Tmax.
 */
export function wilsonLoops(lat, { Rmax, Tmax = Rmax, planes: chosen } = {}) {
  const { V, D, dims } = lat;
  const planes = chosen ?? Array.from({ length: D }, (_a, mu) => Array.from({ length: D - mu - 1 }, (_b, k) => [mu, mu + k + 1])).flat();
  if (!(Rmax >= 1 && Tmax >= 1)) throw new Error('wilsonLoops: Rmax and Tmax must be ≥ 1');
  for (const [mu, nu] of planes) {
    if (!(mu >= 0 && nu >= 0 && mu < D && nu < D && mu !== nu)) throw new Error(`wilsonLoops: bad plane [${mu}, ${nu}]`);
    if (Rmax >= dims[mu] || Tmax >= dims[nu]) throw new Error(`wilsonLoops: loops must not wrap: Rmax < ${dims[mu]}, Tmax < ${dims[nu]} in plane [${mu}, ${nu}]`);
  }
  const values = new Float64Array(Rmax * Tmax), X = new Float64Array(4), Y = new Float64Array(4);
  const lines = new Map(), shift = new Map();
  const L = mu => { if (!lines.has(mu)) { lines.set(mu, lineProducts(lat, mu, Math.max(Rmax, Tmax))); shift.set(mu, shifts(lat, mu, Math.max(Rmax, Tmax))); } return lines.get(mu); };
  const n = Math.max(Rmax, Tmax);
  for (const [mu, nu] of planes) {
    const Lm = L(mu), Ln = L(nu), Sm = shift.get(mu), Sn = shift.get(nu);
    for (let s = 0; s < V; s++) {
      for (let R = 1; R <= Rmax; R++) {
        const sR = Sm[s * (n + 1) + R];
        for (let T = 1; T <= Tmax; T++) {
          const sT = Sn[s * (n + 1) + T];
          qmulInto(X, 0, Lm, (s * n + R - 1) * 4, Ln, (sR * n + T - 1) * 4);     // along μ, then ν
          qmulInto(Y, 0, Ln, (s * n + T - 1) * 4, Lm, (sT * n + R - 1) * 4);     // along ν, then μ
          values[(R - 1) * Tmax + T - 1] += X[0] * Y[0] + X[1] * Y[1] + X[2] * Y[2] + X[3] * Y[3];
        }
      }
    }
  }
  for (let i = 0; i < values.length; i++) values[i] /= V * planes.length;
  const W = (R, T) => {
    if (!(R >= 1 && R <= Rmax && T >= 1 && T <= Tmax)) throw new Error(`wilsonLoops: W(${R}, ${T}) outside 1..${Rmax} × 1..${Tmax}`);
    return values[(R - 1) * Tmax + T - 1];
  };
  return Object.assign(W, { values, Rmax, Tmax });
}

/**
 * Polyakov loop: the ½Tr of the product of all links along `dir` (default the last axis, Euclidean
 * time), from every site with x_dir = 0. In a finite periodic box ⟨P⟩ ≈ 0 means the confined
 * phase; a large |⟨P⟩| means the box is too short/hot to confine. Returns { mean, values }.
 */
export function polyakovLoop(lat, { dir = lat.D - 1 } = {}) {
  const { V, D, U, up, dims } = lat;
  const values = [], P = new Float64Array(4), t = new Float64Array(4);
  for (let s = 0; s < V; s++) {
    if (siteCoords(lat, s)[dir] !== 0) continue;
    P.set(U.subarray(4 * (s * D + dir), 4 * (s * D + dir) + 4));
    let x = up[s * D + dir];
    for (let k = 1; k < dims[dir]; k++) { qmulInto(t, 0, P, 0, U, 4 * (x * D + dir)); P.set(t); x = up[x * D + dir]; }
    values.push(P[0]);
  }
  return { mean: values.reduce((a, b) => a + b, 0) / values.length, values: Float64Array.from(values) };
}

/**
 * Local action density: values[s·nP + p] = 1 − ½Tr U_P for the plaquette at site s in plane
 * planes[p] = [μ, ν] (μ < ν). Its mean is 1 − plaquette(lat). Gauge invariant site by site.
 */
export function plaquetteField(lat) {
  const { V, D, U, up } = lat;
  const planes = [];
  for (let mu = 0; mu < D; mu++) for (let nu = mu + 1; nu < D; nu++) planes.push([mu, nu]);
  const values = new Float64Array(V * planes.length), t = new Float64Array(8);
  for (let s = 0; s < V; s++) {
    planes.forEach(([mu, nu], p) => {
      qmulInto(t, 0, U, 4 * (s * D + mu), U, 4 * (up[s * D + mu] * D + nu));
      qmulInto(t, 4, U, 4 * (s * D + nu), U, 4 * (up[s * D + nu] * D + mu));
      values[s * planes.length + p] = 1 - (t[0] * t[4] + t[1] * t[5] + t[2] * t[6] + t[3] * t[7]);
    });
  }
  return { values, planes };
}

/**
 * APE smearing (Albanese et al., PLB 192 (1987) 163) of the links in `dirs` (default all but the
 * last axis: spatial smearing for time-slice operators), using only staples within `dirs`:
 *   U'_μ = Proj_SU(2)[ (1 − α) U_μ + α/(2(n − 1)) Σ_{ν ∈ dirs, ν ≠ μ} (staples) ],  n = |dirs|.
 * For SU(2) the projection is exact (normalise the quaternion). Links along other axes are copied.
 * Returns a new lattice object (same geometry) and leaves lat unchanged; smeared links are for
 * measurement only, never for updates.
 */
export function apeSmear(lat, { alpha = 0.5, steps = 1, dirs: chosen } = {}) {
  const dirs = chosen ?? Array.from({ length: lat.D - 1 }, (_, k) => k);
  if (!(alpha >= 0 && alpha <= 1)) throw new Error('apeSmear: alpha must be in [0, 1]');
  if (dirs.length < 2) throw new Error('apeSmear: need at least two smearing directions');
  let cur = { ...lat, U: Float64Array.from(lat.U) };
  const A = new Float64Array(4), w = alpha / (2 * (dirs.length - 1));
  for (let k = 0; k < steps; k++) {
    const next = Float64Array.from(cur.U);
    for (let s = 0; s < lat.V; s++) {
      for (const mu of dirs) {
        staple(cur, s, mu, A, dirs);
        const o = 4 * (s * lat.D + mu);
        // staple() returns the dagger of the usual staple sum: Σ staples = A†
        const q = [(1 - alpha) * cur.U[o] + w * A[0], (1 - alpha) * cur.U[o + 1] - w * A[1], (1 - alpha) * cur.U[o + 2] - w * A[2], (1 - alpha) * cur.U[o + 3] - w * A[3]];
        const n = Math.hypot(...q);
        for (let i = 0; i < 4; i++) next[o + i] = q[i] / n;
      }
    }
    cur = { ...cur, U: next };
  }
  return cur;
}

/** Strong-coupling expansion of the average plaquette in d dimensions: u + 2(d − 2)u⁵, u = I₂(β)/I₁(β). Exact in d = 2 (infinite volume). */
export function strongCouplingPlaquette(beta, d) {
  const u = besselI(2, beta) / besselI(1, beta);
  return u + 2 * (d - 2) * u ** 5;
}

// ---------------------------------------------------------------- topology (4D)

/** Ordered product of links along a path from site s; steps are +(μ+1) forward, −(μ+1) backward. */
function walk(lat, s0, steps, out, tmp) {
  const { D, U, up, dn } = lat;
  let s = s0;
  out[0] = 1; out[1] = out[2] = out[3] = 0;
  for (const d of steps) {
    if (d > 0) { qmulInto(tmp, 0, out, 0, U, 4 * (s * D + d - 1)); s = up[s * D + d - 1]; }
    else { s = dn[s * D - d - 1]; qmulInto(tmp, 0, out, 0, U, 4 * (s * D - d - 1), false, true); }
    out.set(tmp);
  }
}
const leafCache = new Map();
/** The four counter-clockwise n×n leaves of the clover in the (μ, ν) plane. */
function leaves(mu, nu, n) {
  const key = `${mu}${nu}${n}`;
  if (!leafCache.has(key)) {
    const m = mu + 1, v = nu + 1, rep = (d, k) => Array(k).fill(d);
    leafCache.set(key, [
      [...rep(m, n), ...rep(v, n), ...rep(-m, n), ...rep(-v, n)],
      [...rep(v, n), ...rep(-m, n), ...rep(-v, n), ...rep(m, n)],
      [...rep(-m, n), ...rep(-v, n), ...rep(m, n), ...rep(v, n)],
      [...rep(-v, n), ...rep(m, n), ...rep(v, n), ...rep(-m, n)],
    ]);
  }
  return leafCache.get(key);
}

/**
 * Gauge-invariant local densities from the clover field strength (4D only). The clover average of
 * n×n leaves gives F^a_μν = 2·vec(¼Σ leaves)/n²; improved: true combines 1×1 and 2×2 clovers as
 * (4F₁ − F₂)/3, cancelling the O(a²) error (de Forcrand, García Pérez & Stamatescu, NPB 499, 409 (1997)).
 *   s(x) = Σ_{μ<ν} |F_μν|²/(16π²)   action density in units of one instanton, 8π²/g² (∫s = 1 for one)
 *   q(x) = (1/8π²) Σ F·F̃ terms       topological charge density, Q = Σ q ≈ integer after smoothing
 * Sites [s0, s1) (default all). Returns { s, q, Q, S } (S = Σ s, so S ≥ |Q| with equality when self-dual).
 * Raw Monte Carlo configurations are dominated by noise: smooth first (sweep with cool: true).
 * Accuracy (tested on bpstInstanton): improved within 0.5 % for ρ ≥ 3.5a; plain clover ≈ 13 % low at ρ = 2.5a.
 */
export function topology(lat, { improved = true, sites = [0, lat.V] } = {}) {
  if (lat.D !== 4) throw new Error('topology: the topological charge density needs d = 4');
  const planes = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
  const s = new Float64Array(lat.V), q = new Float64Array(lat.V);
  const F = planes.map(() => new Float64Array(3)), leaf = new Float64Array(4), tmp = new Float64Array(4);
  const PI2 = Math.PI * Math.PI;
  let Q = 0, S = 0;
  for (let x = sites[0]; x < sites[1]; x++) {
    planes.forEach(([mu, nu], p) => {
      let f1 = 0, f2 = 0, f3 = 0;
      for (const l of leaves(mu, nu, 1)) { walk(lat, x, l, leaf, tmp); f1 += leaf[1]; f2 += leaf[2]; f3 += leaf[3]; }
      F[p][0] = 0.5 * f1; F[p][1] = 0.5 * f2; F[p][2] = 0.5 * f3;
      if (improved) {
        let g1 = 0, g2 = 0, g3 = 0;
        for (const l of leaves(mu, nu, 2)) { walk(lat, x, l, leaf, tmp); g1 += leaf[1]; g2 += leaf[2]; g3 += leaf[3]; }
        F[p][0] = (4 * F[p][0] - g1 / 8) / 3; F[p][1] = (4 * F[p][1] - g2 / 8) / 3; F[p][2] = (4 * F[p][2] - g3 / 8) / 3;
      }
    });
    let ff = 0, fd = 0;
    for (const f of F) ff += f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
    for (let a = 0; a < 3; a++) fd += F[0][a] * F[5][a] - F[1][a] * F[4][a] + F[2][a] * F[3][a];
    s[x] = ff / (16 * PI2); q[x] = fd / (8 * PI2);
    S += s[x]; Q += q[x];
  }
  return { s, q, Q, S };
}

/**
 * An exact BPST instanton of size rho (lattice units), centred between sites of a periodic L⁴
 * lattice, in singular gauge so the field decays as ρ²/|x|³ towards the edges (Belavin et al.,
 * PLB 59, 85 (1975)). Links are path-ordered exponentials of the continuum A_μ (substeps per link).
 * anti: true gives the anti-instanton. Use it as a manufactured benchmark: topology() should give
 * Q = ±1 and S = 1 (self-dual) up to the stated discretisation error, and smoothing should keep it.
 */
export function bpstInstanton({ L = 16, rho = 3.5, anti = false, substeps = 8, beta = 2.4 } = {}) {
  const lat = createLattice({ L, d: 4, beta });
  const sign = -1, gaugeSign = -1, c = (L - 1) / 2, D = 4, U = lat.U;   // conventions fixed by the Q = +1 test
  // 't Hooft symbols η_aμν (a = 0..2 for σ1..3, μ = 0..3 with 3 the Euclidean time)
  const eta = (a, mu, nu) => {
    if (mu === nu) return 0;
    if (mu < 3 && nu < 3) { const e = [a, mu, nu]; return (e[0] - e[1]) * (e[1] - e[2]) * (e[2] - e[0]) / 2; }
    if (nu === 3) return a === mu ? 1 : 0;
    return a === nu ? -1 : 0;
  };
  const A = (x, mu, out) => {
    const den = x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2 + rho * rho;
    for (let a = 0; a < 3; a++) { let v = 0; for (let nu = 0; nu < 4; nu++) v += eta(a, mu, nu) * x[nu]; out[a] = sign * 2 * v / den; }
  };
  const gAt = (x, out) => {   // singular-gauge rotation g(x) = x̂ as a unit quaternion
    const r = Math.hypot(...x);
    if (r === 0) { out.set([1, 0, 0, 0]); return; }
    out[0] = x[3] / r; for (let i = 0; i < 3; i++) out[i + 1] = gaugeSign * x[i] / r;
  };
  const a = new Float64Array(3), cur = new Float64Array(4), nxt = new Float64Array(4), step = new Float64Array(4), g1 = new Float64Array(4), g2 = new Float64Array(4);
  for (let s = 0; s < lat.V; s++) {
    const xi = siteCoords(lat, s).map(v => v - c);
    for (let mu = 0; mu < D; mu++) {
      cur.set([1, 0, 0, 0]);
      for (let k = 0; k < substeps; k++) {
        const xm = xi.slice(); xm[mu] += (k + 0.5) / substeps;
        A(xm, mu, a);
        step.set(qexp(a.map(v => v / (2 * substeps))));
        qmulInto(nxt, 0, cur, 0, step, 0); cur.set(nxt);
      }
      const xe = xi.slice(); xe[mu] += 1;
      gAt(xi, g1); gAt(xe, g2);
      qmulInto(nxt, 0, g1, 0, cur, 0); qmulInto(cur, 0, nxt, 0, g2, 0, false, true);
      U.set(cur, 4 * (s * D + mu));
    }
  }
  if (anti) {
    // time reflection t → L − 1 − t (about the centre) maps Q → −Q exactly:
    // U'_μ(x) = U_μ(Px) for μ < 3, U'_t(x) = U_t(Px − t̂)†
    const R = new Float64Array(U.length);
    for (let s = 0; s < lat.V; s++) {
      const x = siteCoords(lat, s), px = [...x.slice(0, 3), L - 1 - x[3]], pm = [...x.slice(0, 3), L - 2 - x[3]];
      const sp = siteIndex(lat, px), sm = siteIndex(lat, pm);
      for (let mu = 0; mu < 3; mu++) R.set(U.subarray(4 * (sp * D + mu), 4 * (sp * D + mu) + 4), 4 * (s * D + mu));
      const o = 4 * (sm * D + 3);
      R.set([U[o], -U[o + 1], -U[o + 2], -U[o + 3]], 4 * (s * D + 3));
    }
    U.set(R);
  }
  return lat;
}
