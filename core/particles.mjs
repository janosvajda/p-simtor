// Particle systems with pluggable forces, for any classical point-particle model.
//
//   const sys = createParticles({
//     bodies: [{ m: 1, x: [1, 0, 0], v: [0, 1, 0] }, { m: 1, charge: -1, x: [...], v: [...] }],
//     forces: [gravity({ G: 1 }), uniformField({ B: [0, 0, 1] }), drag({ gamma: 0.1 })],
//   });
//   sys.advanceTo(10); sys.sampleAt(t, outX); sys.diagnostics();
//
// Integrator choice follows the mathematics: position-only forces → symplectic Yoshida-4
// (bounded energy error); any velocity-dependent force (magnetism, drag, relativistic
// corrections) → adaptive Dormand–Prince on the full phase-space state, with the work done by
// dissipative forces integrated alongside so energy balance can always be checked.
//
// Units are the caller's: keep G, k, c, fields and masses consistent.
//
// Periodic boxes: box: [Lx, Ly, Lz] (0 or omitted = open axis; [Lx, Ly] for a 2D box). Positions
// stay unwrapped (smooth integration, displacements for diffusion); short-range forces use the
// minimum image, and wrapped(out) folds positions into the box for drawing. Long-range forces
// (gravity, Coulomb) refuse a periodic box: that needs Ewald summation, which is not provided.

import { mechanicsStepper, odeStepper } from './ode.mjs';

// ---------------------------------------------------------------- forces
// A force adds accelerations: accumulate(x, v, t, bodies, a) for flat arrays (3 per body).
// Optional: potential(x, bodies) (conservative part), power(x, v, bodies) (rate of work done
// on the bodies by a non-conservative force), energyCorrection(x, v, bodies) (extra terms of a
// modified conserved energy), internal: true for pairwise forces (momentum is conserved).
// Short-range forces receive a sixth argument env with env.pairsWithin(x, rc, fn): cell-list
// neighbour search with the minimum image along periodic axes.

/**
 * Call fn(i, j, dx, dy, dz, r2) once for every pair closer than rc, where (dx, dy, dz) = xⱼ − xᵢ
 * (minimum image along periodic axes). Cell lists: cost ∝ n at fixed density.
 * x: flat positions (3 per body); box: [Lx, Ly, Lz] with 0 for open axes.
 */
export function forEachPairWithin(x, n, rc, box, fn) {
  const L = [0, 1, 2].map(k => (box?.[k] > 0 && Number.isFinite(box[k]) ? box[k] : 0));
  if (L.some(l => l && rc > l / 2)) throw new Error('forEachPairWithin: the cutoff must be at most half the box length (minimum image)');
  const lo = [0, 0, 0], nc = [1, 1, 1], cs = [1, 1, 1];
  for (let k = 0; k < 3; k++) {
    let min = 0, span = L[k];
    if (!span) {
      let max = -Infinity; min = Infinity;
      for (let i = 0; i < n; i++) { min = Math.min(min, x[3 * i + k]); max = Math.max(max, x[3 * i + k]); }
      span = Math.max(max - min, 0) * (1 + 1e-12) + 1e-300;
    }
    lo[k] = min; nc[k] = Math.max(1, Math.floor(span / rc)); cs[k] = span / nc[k];
  }
  const cellOf = (i, k) => {
    const c = Math.floor((x[3 * i + k] - lo[k]) / cs[k]);
    return L[k] ? ((c % nc[k]) + nc[k]) % nc[k] : Math.min(nc[k] - 1, Math.max(0, c));
  };
  const head = new Int32Array(nc[0] * nc[1] * nc[2]).fill(-1), next = new Int32Array(n);
  for (let i = 0; i < n; i++) { const c = cellOf(i, 0) + nc[0] * (cellOf(i, 1) + nc[1] * cellOf(i, 2)); next[i] = head[c]; head[c] = i; }
  // distinct neighbouring cells along each axis (1 or 2 cells per axis must not double count)
  const around = (c, k) => [...new Set([c - 1, c, c + 1].filter(v => L[k] || (v >= 0 && v < nc[k])).map(v => ((v % nc[k]) + nc[k]) % nc[k]))];
  const rc2 = rc * rc;
  for (let cz = 0; cz < nc[2]; cz++) {for (let cy = 0; cy < nc[1]; cy++) {for (let cx = 0; cx < nc[0]; cx++) {
    const c = cx + nc[0] * (cy + nc[1] * cz);
    if (head[c] < 0) continue;
    for (const nz of around(cz, 2)) {for (const ny of around(cy, 1)) {for (const nx of around(cx, 0)) {
      const d = nx + nc[0] * (ny + nc[1] * nz);
      for (let i = head[c]; i >= 0; i = next[i]) {for (let j = head[d]; j >= 0; j = next[j]) {
        if (j <= i) continue;
        let dx = x[3 * j] - x[3 * i], dy = x[3 * j + 1] - x[3 * i + 1], dz = x[3 * j + 2] - x[3 * i + 2];
        if (L[0]) dx -= L[0] * Math.round(dx / L[0]);
        if (L[1]) dy -= L[1] * Math.round(dy / L[1]);
        if (L[2]) dz -= L[2] * Math.round(dz / L[2]);
        const r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < rc2) fn(i, j, dx, dy, dz, r2);
      }}
    }}}
  }}}
}

/**
 * Lennard-Jones V = 4ε[(σ/r)¹² − (σ/r)⁶] between all bodies, cut at rc (default 2.5σ).
 * shift: 'force' (default: V and F both vanish at rc, so energy is smooth and Verlet keeps its
 * order), 'energy' (V shifted to 0 at rc; the force jumps there) or 'none'. Per-body sigma and
 * epsilon override the defaults with Lorentz–Berthelot mixing σᵢⱼ = (σᵢ+σⱼ)/2, εᵢⱼ = √(εᵢεⱼ);
 * the cutoff then scales as rc·σᵢⱼ/σ. A cutoff changes the model (pressure, phase boundaries).
 */
export function lennardJones({ epsilon = 1, sigma = 1, cutoff = 2.5 * sigma, shift = 'force' } = {}) {
  if (!(epsilon > 0 && sigma > 0 && cutoff > 0)) throw new Error('lennardJones: epsilon, sigma and cutoff must be positive');
  if (!['force', 'energy', 'none'].includes(shift)) throw new Error('lennardJones: shift must be force, energy or none');
  const params = (b, i, j) => {
    const mixed = b[i].sigma != null || b[j].sigma != null || b[i].epsilon != null || b[j].epsilon != null;
    if (!mixed) return [sigma, epsilon, cutoff];
    const s = ((b[i].sigma ?? sigma) + (b[j].sigma ?? sigma)) / 2;
    return [s, Math.sqrt((b[i].epsilon ?? epsilon) * (b[j].epsilon ?? epsilon)), cutoff * s / sigma];
  };
  const V = (r, s, e) => { const q = (s / r) ** 6; return 4 * e * (q * q - q); };
  const dV = (r, s, e) => { const q = (s / r) ** 6; return 24 * e * (q - 2 * q * q) / r; };    // dV/dr
  const reach = b => cutoff * b.reduce((m, q) => Math.max(m, q.sigma ?? sigma), sigma) / sigma;
  return {
    name: 'Lennard-Jones', internal: true, velocityDependent: false, shortRange: true,
    accumulate(x, _v, _t, b, a, env) {
      env.pairsWithin(x, reach(b), (i, j, dx, dy, dz, r2) => {
        const [s, e, rc] = params(b, i, j), r = Math.sqrt(r2);
        if (r >= rc) return;
        const f = (dV(r, s, e) - (shift === 'force' ? dV(rc, s, e) : 0)) / r;   // force on i = f·(xⱼ − xᵢ)
        a[3 * i] += f * dx / b[i].m; a[3 * i + 1] += f * dy / b[i].m; a[3 * i + 2] += f * dz / b[i].m;
        a[3 * j] -= f * dx / b[j].m; a[3 * j + 1] -= f * dy / b[j].m; a[3 * j + 2] -= f * dz / b[j].m;
      });
    },
    potential(x, b, env) {
      let U = 0;
      env.pairsWithin(x, reach(b), (i, j, _dx, _dy, _dz, r2) => {
        const [s, e, rc] = params(b, i, j), r = Math.sqrt(r2);
        if (r >= rc) return;
        U += V(r, s, e);
        if (shift !== 'none') U -= V(rc, s, e);
        if (shift === 'force') U -= (r - rc) * dV(rc, s, e);
      });
      return U;
    },
  };
}

const pairLoop = (n, fn) => { for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) fn(i, j); };

/** Newtonian gravity between all bodies, optional Plummer softening ε. */
export function gravity({ G = 1, softening = 0 } = {}) {
  const e2 = softening * softening;
  return {
    name: 'gravity', internal: true, velocityDependent: false, longRange: true,
    accumulate(x, _v, _t, b, a) {
      pairLoop(b.length, (i, j) => {
        const dx = x[3 * j] - x[3 * i], dy = x[3 * j + 1] - x[3 * i + 1], dz = x[3 * j + 2] - x[3 * i + 2];
        const r2 = dx * dx + dy * dy + dz * dz + e2, f = G / (r2 * Math.sqrt(r2));
        a[3 * i] += b[j].m * f * dx; a[3 * i + 1] += b[j].m * f * dy; a[3 * i + 2] += b[j].m * f * dz;
        a[3 * j] -= b[i].m * f * dx; a[3 * j + 1] -= b[i].m * f * dy; a[3 * j + 2] -= b[i].m * f * dz;
      });
    },
    potential(x, b) {
      let U = 0;
      pairLoop(b.length, (i, j) => {
        const r = Math.hypot(x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1], x[3 * j + 2] - x[3 * i + 2]);
        U -= G * b[i].m * b[j].m / Math.sqrt(r * r + e2);
      });
      return U;
    },
  };
}

/** Coulomb force k qᵢqⱼ/r² between charged bodies (body.charge), optional softening. */
export function coulomb({ k = 1, softening = 0 } = {}) {
  const e2 = softening * softening;
  return {
    name: 'coulomb', internal: true, velocityDependent: false, longRange: true,
    accumulate(x, _v, _t, b, a) {
      pairLoop(b.length, (i, j) => {
        const qq = (b[i].charge ?? 0) * (b[j].charge ?? 0);
        if (!qq) return;
        const dx = x[3 * j] - x[3 * i], dy = x[3 * j + 1] - x[3 * i + 1], dz = x[3 * j + 2] - x[3 * i + 2];
        const r2 = dx * dx + dy * dy + dz * dz + e2, f = k * qq / (r2 * Math.sqrt(r2));
        a[3 * i] -= f * dx / b[i].m; a[3 * i + 1] -= f * dy / b[i].m; a[3 * i + 2] -= f * dz / b[i].m;
        a[3 * j] += f * dx / b[j].m; a[3 * j + 1] += f * dy / b[j].m; a[3 * j + 2] += f * dz / b[j].m;
      });
    },
    potential(x, b) {
      let U = 0;
      pairLoop(b.length, (i, j) => {
        const r = Math.hypot(x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1], x[3 * j + 2] - x[3 * i + 2]);
        U += k * (b[i].charge ?? 0) * (b[j].charge ?? 0) / Math.sqrt(r * r + e2);
      });
      return U;
    },
  };
}

/** Uniform electric and magnetic fields: a = (q/m)(E + v × B). The magnetic part does no work. */
export function uniformField({ E = [0, 0, 0], B = [0, 0, 0] } = {}) {
  const hasB = B.some(c => c !== 0);
  return {
    name: 'uniform field', internal: false, velocityDependent: hasB,
    accumulate(_x, v, _t, b, a) {
      for (let i = 0; i < b.length; i++) {
        const qm = (b[i].charge ?? 0) / b[i].m;
        if (!qm) continue;
        const vx = v[3 * i], vy = v[3 * i + 1], vz = v[3 * i + 2];
        a[3 * i] += qm * (E[0] + vy * B[2] - vz * B[1]);
        a[3 * i + 1] += qm * (E[1] + vz * B[0] - vx * B[2]);
        a[3 * i + 2] += qm * (E[2] + vx * B[1] - vy * B[0]);
      }
    },
    potential(x, b) {
      let U = 0;
      for (let i = 0; i < b.length; i++) U -= (b[i].charge ?? 0) * (E[0] * x[3 * i] + E[1] * x[3 * i + 1] + E[2] * x[3 * i + 2]);
      return U;
    },
  };
}

/**
 * A fixed central mass GM at `center` acting on every body (test-particle approximation).
 * relativity: true adds the first post-Newtonian (Schwarzschild, harmonic-gauge) correction
 *   a_GR = GM/(c² r³) [ (4GM/r − v²) r + 4 (r·v) v ],
 * which gives the perihelion advance 6πGM/(c² a (1 − e²)) per orbit.
 */
export function centralBody({ GM, center = [0, 0, 0], relativity = false, c = 299792458 } = {}) {
  if (!(GM > 0)) throw new Error('centralBody: GM must be positive');
  const c2 = c * c;
  return {
    name: relativity ? 'central body (Newton + 1PN)' : 'central body', internal: false, velocityDependent: relativity,
    accumulate(x, v, _t, b, a) {
      for (let i = 0; i < b.length; i++) {
        const rx = x[3 * i] - center[0], ry = x[3 * i + 1] - center[1], rz = x[3 * i + 2] - center[2];
        const r = Math.hypot(rx, ry, rz), r3 = r * r * r;
        a[3 * i] -= GM * rx / r3; a[3 * i + 1] -= GM * ry / r3; a[3 * i + 2] -= GM * rz / r3;
        if (relativity) {
          const vx = v[3 * i], vy = v[3 * i + 1], vz = v[3 * i + 2];
          const v2 = vx * vx + vy * vy + vz * vz, rv = rx * vx + ry * vy + rz * vz;
          const k = GM / (c2 * r3), s = 4 * GM / r - v2;
          a[3 * i] += k * (s * rx + 4 * rv * vx); a[3 * i + 1] += k * (s * ry + 4 * rv * vy); a[3 * i + 2] += k * (s * rz + 4 * rv * vz);
        }
      }
    },
    potential(x, b) {
      let U = 0;
      for (let i = 0; i < b.length; i++) U -= GM * b[i].m / Math.hypot(x[3 * i] - center[0], x[3 * i + 1] - center[1], x[3 * i + 2] - center[2]);
      return U;
    },
    // 1PN conserved energy per unit mass adds (1/c²)[3v⁴/8 + 3GMv²/(2r) + (GM/r)²/2]
    energyCorrection: relativity ? (x, v, b) => {
      let E = 0;
      for (let i = 0; i < b.length; i++) {
        const r = Math.hypot(x[3 * i] - center[0], x[3 * i + 1] - center[1], x[3 * i + 2] - center[2]);
        const v2 = v[3 * i] ** 2 + v[3 * i + 1] ** 2 + v[3 * i + 2] ** 2, u = GM / r;
        E += b[i].m * (3 * v2 * v2 / 8 + 1.5 * u * v2 + 0.5 * u * u) / c2;
      }
      return E;
    } : undefined,
  };
}

/** Linear drag a = −γ v (Stokes regime); dissipative. */
export function drag({ gamma } = {}) {
  if (!(gamma >= 0)) throw new Error('drag: gamma must be ≥ 0');
  return {
    name: 'linear drag', internal: false, velocityDependent: true,
    accumulate(_x, v, _t, b, a) { for (let k = 0; k < 3 * b.length; k++) a[k] -= gamma * v[k]; },
    power(_x, v, b) {
      let P = 0;
      for (let i = 0; i < b.length; i++) P -= gamma * b[i].m * (v[3 * i] ** 2 + v[3 * i + 1] ** 2 + v[3 * i + 2] ** 2);
      return P;
    },
  };
}

/** Hookean springs between body pairs: [[i, j, stiffness k, rest length L0], ...]. */
export function springs(pairs) {
  return {
    name: 'springs', internal: true, velocityDependent: false,
    accumulate(x, _v, _t, b, a) {
      for (const [i, j, k, L0] of pairs) {
        const dx = x[3 * j] - x[3 * i], dy = x[3 * j + 1] - x[3 * i + 1], dz = x[3 * j + 2] - x[3 * i + 2];
        const L = Math.hypot(dx, dy, dz) || 1e-300, f = k * (L - L0) / L;
        a[3 * i] += f * dx / b[i].m; a[3 * i + 1] += f * dy / b[i].m; a[3 * i + 2] += f * dz / b[i].m;
        a[3 * j] -= f * dx / b[j].m; a[3 * j + 1] -= f * dy / b[j].m; a[3 * j + 2] -= f * dz / b[j].m;
      }
    },
    potential(x) {
      let U = 0;
      for (const [i, j, k, L0] of pairs) {
        const L = Math.hypot(x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1], x[3 * j + 2] - x[3 * i + 2]);
        U += 0.5 * k * (L - L0) ** 2;
      }
      return U;
    },
  };
}

/** Any external potential energy U(position, body) with its gradient ∇U (per body). */
export function externalPotential({ U, gradU, name = 'external potential' }) {
  const g = [0, 0, 0];
  return {
    name, internal: false, velocityDependent: false,
    accumulate(x, _v, _t, b, a) {
      for (let i = 0; i < b.length; i++) {
        gradU([x[3 * i], x[3 * i + 1], x[3 * i + 2]], b[i], g);
        a[3 * i] -= g[0] / b[i].m; a[3 * i + 1] -= g[1] / b[i].m; a[3 * i + 2] -= g[2] / b[i].m;
      }
    },
    potential(x, b) {
      let s = 0;
      for (let i = 0; i < b.length; i++) s += U([x[3 * i], x[3 * i + 1], x[3 * i + 2]], b[i]);
      return s;
    },
  };
}

// ---------------------------------------------------------------- system

/**
 * bodies: [{ m, charge?, sigma?, epsilon?, x: [x,y,z], v: [vx,vy,vz] }]; forces: from the constructors above;
 * box: optional periodic lengths (see top of file).
 * method: 'auto' | 'yoshida4' | 'leapfrog' | 'dopri5'. h: fixed step (symplectic) or first
 * trial step (adaptive). rtol/atol: adaptive tolerances.
 */
export function createParticles({ bodies, forces, method = 'auto', h, rtol = 1e-10, atol = 1e-12, t0 = 0, box }) {
  if (!Array.isArray(bodies) || !bodies.length) throw new Error('createParticles: bodies are required');
  for (const [i, b] of bodies.entries()) {
    if (!(b.m > 0)) throw new Error(`createParticles: body ${i} needs a positive mass`);
    for (const k of ['x', 'v']) if (!Array.isArray(b[k]) || b[k].length < 2 || !b[k].every(Number.isFinite)) throw new Error(`createParticles: body ${i}.${k} must be finite [x, y] or [x, y, z]`);
  }
  if (!(Number.isFinite(h) && h > 0)) throw new Error('createParticles: a positive finite step h is required');
  const n = bodies.length, N = 3 * n;
  const L = [0, 1, 2].map(k => (box?.[k] > 0 && Number.isFinite(box[k]) ? box[k] : 0));
  const longRange = L.some(l => l) && forces.find(f => f.longRange);
  if (longRange) throw new Error(`createParticles: ${longRange.name} is long-range; a periodic box needs Ewald summation, which is not provided`);
  const env = { box: L, pairsWithin: (x, rc, fn) => forEachPairWithin(x, n, rc, L, fn) };
  const bs = bodies.map(b => ({ ...b }));
  const x0 = new Float64Array(N), v0 = new Float64Array(N);
  bs.forEach((b, i) => { for (let k = 0; k < 3; k++) { x0[3 * i + k] = b.x[k] ?? 0; v0[3 * i + k] = b.v[k] ?? 0; } });
  const velocityDependent = forces.some(f => f.velocityDependent);
  const dissipative = forces.filter(f => f.power);
  const chosen = method === 'auto' ? (velocityDependent ? 'dopri5' : 'yoshida4') : method;
  if (velocityDependent && chosen !== 'dopri5') throw new Error('createParticles: velocity-dependent forces need method "dopri5" (or "auto")');

  const a = new Float64Array(N);
  const totalAccel = (t, x, v) => { a.fill(0); for (const f of forces) f.accumulate(x, v, t, bs, a, env); return a; };

  let stepper, xView, vView, workView;
  if (chosen === 'dopri5') {
    // phase-space state [x, v, W] where W = work done on the system by dissipative forces
    const y0 = new Float64Array(2 * N + 1); y0.set(x0); y0.set(v0, N);
    const dydt = new Float64Array(2 * N + 1);
    stepper = odeStepper({
      method: 'dopri5', y0, t0, h, rtol, atol,
      f(t, y) {
        const x = y.subarray(0, N), v = y.subarray(N, 2 * N);
        totalAccel(t, x, v);
        dydt.set(v, 0); dydt.set(a, N);
        let P = 0;
        for (const f of dissipative) P += f.power(x, v, bs);
        dydt[2 * N] = P;
        return dydt;
      },
    });
    xView = () => stepper.y.subarray(0, N); vView = () => stepper.y.subarray(N, 2 * N); workView = () => stepper.y[2 * N];
  } else {
    stepper = mechanicsStepper({ q0: x0, v0, t0, h, method: chosen, accel: (t, x, out) => { out.set(totalAccel(t, x, null)); } });
    xView = () => stepper.q; vView = () => stepper.v; workView = () => 0;
  }

  const energy = () => {
    const x = xView(), v = vView();
    let T = 0, U = 0, extra = 0;
    for (let i = 0; i < n; i++) T += 0.5 * bs[i].m * (v[3 * i] ** 2 + v[3 * i + 1] ** 2 + v[3 * i + 2] ** 2);
    for (const f of forces) {
      if (f.potential) U += f.potential(x, bs, env);
      if (f.energyCorrection) extra += f.energyCorrection(x, v, bs);
    }
    return { T, U, E: T + U + extra };
  };
  const momentum = () => {
    const v = vView(), p = [0, 0, 0];
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) p[k] += bs[i].m * v[3 * i + k];
    return p;
  };
  const E0 = energy().E, P0 = momentum();
  const conservesMomentum = forces.every(f => f.internal);
  const scaleP = bs.reduce((s, b, i) => s + b.m * Math.hypot(v0[3 * i], v0[3 * i + 1], v0[3 * i + 2]), 0) || 1;

  const sys = {
    bodies: bs, forces, method: chosen, n,
    get t() { return stepper.t; },
    get x() { return xView(); },
    get v() { return vView(); },
    /** Work done on the system by dissipative forces since t0 (negative for drag). */
    get work() { return workView(); },
    box: L,
    /** Positions folded into [0, L) along periodic axes (for drawing); x stays unwrapped. */
    wrapped(out = new Float64Array(N)) {
      const x = xView();
      for (let k = 0; k < N; k++) { const l = L[k % 3]; out[k] = l ? x[k] - l * Math.floor(x[k] / l) : x[k]; }
      return out;
    },
    step: () => stepper.step(),
    advanceTo: (t, maxSteps) => stepper.advanceTo(t, maxSteps),
    get steps() { return stepper.steps; },
    /** Positions at time t within the last step, for drawing. */
    sampleAt(t, out = new Float64Array(N)) {
      if (chosen === 'dopri5') { const y = stepper.sampleAt(t); out.set(y.subarray(0, N)); return out; }
      return stepper.sampleAt(t, out);
    },
    /** Positions and velocities at time t within the last step (e.g. to locate events precisely). */
    sampleStateAt(t) {
      if (chosen === 'dopri5') { const y = stepper.sampleAt(t); return { x: y.slice(0, N), v: y.slice(N, 2 * N) }; }
      const x = new Float64Array(N), v = new Float64Array(N);
      stepper.sampleAt(t, x, v);
      return { x, v };
    },
    energy,
    /**
     * Energy balance |E − E₀ − W| / scale, where W is the work done by dissipative forces:
     * conservation when there is no dissipation, the work–energy theorem when there is.
     */
    energyBalanceError() {
      const { E, T, U } = energy();
      const scale = Math.max(Math.abs(E0), Math.abs(T) + Math.abs(U), 1e-300);
      return Math.abs(E - E0 - workView()) / scale;
    },
    /** Momentum drift normalised by Σm|v₀|, or null when external forces make momentum non-conserved. */
    momentumDrift() {
      if (!conservesMomentum) return null;
      const p = momentum();
      return Math.hypot(p[0] - P0[0], p[1] - P0[1], p[2] - P0[2]) / scaleP;
    },
  };
  return sys;
}
