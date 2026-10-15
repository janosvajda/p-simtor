// Grid evolution that complements core/fields.mjs (same storage and boundary conventions):
//   diffusion1D        ∂u/∂t = D ∂²u/∂x² (+ source)  Crank–Nicolson: walls or periodic, sources, 2nd order
//   waveND             ∂²u/∂t² = c² ∇²u               leapfrog in any dimension, CFL-checked, exact discrete energy
//   stationaryStates1D Hψ = Eψ                        finite differences + symmetric tridiagonal eigensolver
//   transmission1D     scattering T(E), R(E) through any 1D barrier (stationary, integrated from the far side)
//   thermalAverage     (1/kT)∫ f(E) e^{−E/kT} dE  (Bell's tunnelling rate factor, thermionic emission)
//   packetTransmission ∫ |φ(k)|² T(k) dk for a Gaussian packet: the exact benchmark for wave-packet runs
// Source-free periodic diffusion and time-dependent Schrödinger in any dimension: fields.mjs.
//
// Boundary conditions as in laplacianND: 'dirichlet' (zero ghost values), 'neumann' (copied
// edge-cell ghosts: zero face flux, cell-centred) or 'periodic'.

import { solveTridiagonal, lowestTridiagonal, dopri5, gaussLegendre } from './numerics.mjs';
import { laplacianND } from './fields.mjs';

const BCS = ['dirichlet', 'neumann', 'periodic'];
const checkGrid = (u, dx, who) => {
  if (!(u.length >= 3) || !Array.from(u).every(Number.isFinite)) throw new Error(`${who}: need at least 3 finite grid values`);
  if (!(Number.isFinite(dx) && dx > 0)) throw new Error(`${who}: dx must be positive and finite`);
};

/** Solve a cyclic tridiagonal system (periodic boundaries) by Sherman–Morrison. */
function solveCyclic(a, b, c, alpha, beta, d) {
  const n = b.length, gamma = -b[0];
  const bb = Float64Array.from(b);
  bb[0] = b[0] - gamma; bb[n - 1] = b[n - 1] - alpha * beta / gamma;
  const x = solveTridiagonal(a, bb, c, d);
  const uvec = new Float64Array(n); uvec[0] = gamma; uvec[n - 1] = alpha;
  const z = solveTridiagonal(a, bb, c, uvec);
  const fact = (x[0] + beta * x[n - 1] / gamma) / (1 + z[0] + beta * z[n - 1] / gamma);
  for (let i = 0; i < n; i++) x[i] -= fact * z[i];
  return x;
}

/** Normalise a wall: legacy strings, or { type: 'value' | 'flux', value: number | (t => number) }. */
function wall(spec, who) {
  if (spec === 'neumann') return { type: 'flux', value: () => 0 };
  if (spec === 'dirichlet') return { type: 'ghost' };                 // zero ghost value (laplacianND convention)
  if (spec?.type === 'value' || spec?.type === 'flux') {
    const v = spec.value ?? 0, f = typeof v === 'function' ? v : () => v;
    if (!Number.isFinite(f(0))) throw new Error(`${who}: boundary value must be finite`);
    return { type: spec.type, value: f };
  }
  throw new Error(`${who}: a wall is 'dirichlet', 'neumann', { type: 'value', value } or { type: 'flux', value }`);
}

/**
 * Diffusion ∂u/∂t = D ∂²u/∂x² + s(i, t) by Crank–Nicolson on cells x_i = x₀ + (i + ½)dx.
 * Walls: bc 'dirichlet' | 'neumann' | 'periodic' for both, or left/right separately:
 *   { type: 'value', value }  u = value at the wall face (second order; value may be a function of t)
 *   { type: 'flux', value }   inward flux −D∂u/∂n = value at the wall (0 = insulated)
 * setBoundary('left' | 'right', wall) switches a wall mid-run (e.g. a dissolving solid used up).
 * inflow is the amount that entered through the walls since t0 (discrete, consistent with the
 * scheme), so total() − total₀ − ∫source = inflow to round-off.
 * Returns { t, u, step(), advanceTo(t), total(), inflow, setBoundary() }.
 */
export function diffusion1D({ u0, D, dx, dt, bc = 'neumann', left, right, source, t0 = 0 }) {
  checkGrid(u0, dx, 'diffusion1D');
  if (!(D > 0) || !(dt > 0)) throw new Error('diffusion1D: D and dt must be positive');
  if (!BCS.includes(bc)) throw new Error(`diffusion1D: bc must be one of ${BCS.join(', ')}`);
  const periodic = bc === 'periodic';
  if (periodic && (left || right)) throw new Error('diffusion1D: periodic boundaries have no walls');
  const n = u0.length, r = D * dt / (dx * dx);
  const walls = periodic ? null : [wall(left ?? bc, 'diffusion1D'), wall(right ?? bc, 'diffusion1D')];
  const a = new Float64Array(n).fill(-r / 2), b = new Float64Array(n).fill(1 + r), c = new Float64Array(n).fill(-r / 2);
  // diagonal of the discrete Laplacian (×dx²) in the edge cell for each wall type
  const edge = { ghost: -2, value: -3, flux: -1 };
  const setDiag = () => { if (walls) { b[0] = 1 - (r / 2) * edge[walls[0].type]; b[n - 1] = 1 - (r / 2) * edge[walls[1].type]; } };
  setDiag();
  const u = Float64Array.from(u0), rhs = new Float64Array(n), next = new Float64Array(n);
  const lapEdge = (w, inner, own, t) => {                           // ∂²u/∂x² in an edge cell, × dx²
    if (w.type === 'ghost') return inner - 2 * own;
    if (w.type === 'value') return inner - 3 * own + 2 * w.value(t);
    return inner - own + w.value(t) * dx / D;
  };
  // amount entering through a wall during the step, from the same discrete fluxes the scheme uses
  const wallIn = (w, ownOld, ownNew, t) => {
    if (w.type === 'flux') return 0.5 * dt * (w.value(t) + w.value(t + dt));
    const face = (own, tt) => (w.type === 'value' ? 2 * (w.value(tt) - own) : -own) * D / dx;
    return 0.5 * dt * (face(ownOld, t) + face(ownNew, t + dt));
  };
  const s = {
    t: t0, u, steps: 0, inflow: 0,
    step() {
      const t = s.t;
      for (let i = 1; i < n - 1; i++) rhs[i] = u[i] + (r / 2) * (u[i - 1] - 2 * u[i] + u[i + 1]);
      if (periodic) {
        rhs[0] = u[0] + (r / 2) * (u[n - 1] - 2 * u[0] + u[1]);
        rhs[n - 1] = u[n - 1] + (r / 2) * (u[n - 2] - 2 * u[n - 1] + u[0]);
      } else {
        const [L, R] = walls;
        rhs[0] = u[0] + (r / 2) * lapEdge(L, u[1], u[0], t);
        rhs[n - 1] = u[n - 1] + (r / 2) * lapEdge(R, u[n - 2], u[n - 1], t);
        // implicit half: constant parts of the boundary terms at t + dt
        if (L.type === 'value') rhs[0] += r * L.value(t + dt);
        if (L.type === 'flux') rhs[0] += (dt / (2 * dx)) * L.value(t + dt);
        if (R.type === 'value') rhs[n - 1] += r * R.value(t + dt);
        if (R.type === 'flux') rhs[n - 1] += (dt / (2 * dx)) * R.value(t + dt);
      }
      if (source) for (let i = 0; i < n; i++) rhs[i] += dt * 0.5 * (source(i, t) + source(i, t + dt));
      next.set(periodic ? solveCyclic(a, b, c, -r / 2, -r / 2, rhs) : solveTridiagonal(a, b, c, rhs));
      if (walls) s.inflow += wallIn(walls[0], u[0], next[0], t) + wallIn(walls[1], u[n - 1], next[n - 1], t);
      u.set(next); s.t += dt; s.steps++;
    },
    advanceTo(target) { while (s.t < target - 1e-12 * dt) s.step(); },
    /** Σ u dx. Changes only by inflow through the walls and by the source. */
    total() { let sum = 0; for (let i = 0; i < n; i++) sum += u[i]; return sum * dx; },
    setBoundary(side, spec) {
      if (periodic) throw new Error('diffusion1D: periodic boundaries have no walls');
      if (side !== 'left' && side !== 'right') throw new Error("diffusion1D: side must be 'left' or 'right'");
      walls[side === 'left' ? 0 : 1] = wall(spec, 'diffusion1D');
      setDiag();
    },
  };
  return s;
}

/**
 * Wave equation ∂²u/∂t² = c²∇²u on an any-dimensional grid (fields.mjs storage: axis 0 fastest),
 * leapfrog in time with the compact laplacianND. Requires c·dt·√(Σ 1/dxₐ²) ≤ 1 (CFL).
 * Optional mask/maskBc restrict it to an irregular domain (see fields.mjs; values outside stay 0).
 * v0 is the initial ∂u/∂t. energy() is the discrete energy that leapfrog conserves exactly
 * (to round-off); it is positive whenever the CFL condition holds.
 * Returns { t, u, step(), advanceTo(t), energy(), courant }.
 */
export function waveND({ u0, v0, shape, spacing = 1, bc = 'dirichlet', mask, maskBc, c, dt, t0 = 0 }) {
  const grid = { shape, spacing, bc, mask, maskBc }, n = u0.length;
  const dxs = typeof spacing === 'number' ? Array.from(shape, () => spacing) : Array.from(spacing);
  const courant = c * dt * Math.sqrt(dxs.reduce((sum, d) => sum + 1 / (d * d), 0));
  if (!(courant > 0 && courant <= 1)) throw new Error(`waveND: Courant number c·dt·√(Σ1/dx²) = ${courant} must be in (0, 1]`);
  if (v0 && v0.length !== n) throw new Error('waveND: v0 must match u0');
  const cellVolume = dxs.reduce((p, d) => p * d, 1);
  const u = Float64Array.from(u0, (value, i) => (mask && !mask[i] ? 0 : value)), prev = new Float64Array(n), next = new Float64Array(n), lap = new Float64Array(n);
  laplacianND(u, grid, lap);
  // u(−dt) from a Taylor expansion: u − dt·v + ½dt²c²∇²u
  for (let i = 0; i < n; i++) prev[i] = u[i] - dt * (mask && !mask[i] ? 0 : v0?.[i] ?? 0) + 0.5 * dt * dt * c * c * lap[i];
  const s = {
    t: t0, u, steps: 0, courant,
    step() {
      laplacianND(u, grid, lap);
      for (let i = 0; i < n; i++) next[i] = 2 * u[i] - prev[i] + dt * dt * c * c * lap[i];
      prev.set(u); u.set(next); s.t += dt; s.steps++;
    },
    advanceTo(target) { while (s.t < target - 1e-12 * dt) s.step(); },
    /** ½Σ[((uⁿ − uⁿ⁻¹)/dt)² − c²uⁿ·∇²uⁿ⁻¹]dV: exactly conserved by leapfrog because ∇² is symmetric. */
    energy() {
      laplacianND(prev, grid, lap);
      let E = 0;
      for (let i = 0; i < n; i++) { const vt = (u[i] - prev[i]) / dt; E += vt * vt - c * c * u[i] * lap[i]; }
      return 0.5 * E * cellVolume;
    },
  };
  return s;
}

/**
 * Lowest stationary states of H = −ħ²/2m d²/dx² + V(x) on a grid with ψ = 0 beyond both ends
 * (second-order finite differences). Returns { energies, states } with states normalised so that
 * Σ|ψ|²dx = 1. Accuracy O(dx²); confine the domain so ψ is negligible at the walls.
 * The `count` lowest energies by Sturm bisection and their states by shifted inverse iteration
 * (both O(n) per level), orthogonalised within near-degenerate groups (tunnelling doublets):
 * 10⁴ grid points take milliseconds.
 */
export function stationaryStates1D({ V, dx, count = 5, mass = 1, hbar = 1 }) {
  const n = V.length;
  if (!(n >= 3) || !(dx > 0)) throw new Error('stationaryStates1D: need ≥ 3 grid points and dx > 0');
  const t = hbar * hbar / (2 * mass * dx * dx);
  const diag = Float64Array.from(V, v => v + 2 * t), off = new Float64Array(n - 1).fill(-t);
  const energies = lowestTridiagonal(diag, off, count);
  const scale = Math.max(...energies.map(Math.abs), t);
  const states = [];
  energies.forEach((E, k) => {
    // shifted inverse iteration: (H − (E + δ)) x = b converges to the eigenvector of E
    const shift = E + 1e-10 * scale * (1 + k * 1e-3);
    const a = new Float64Array(n).fill(-t), c = new Float64Array(n).fill(-t), b = Float64Array.from(diag, v => v - shift);
    let x = Float64Array.from({ length: n }, (_, i) => Math.sin((k + 1) * Math.PI * (i + 1) / (n + 1)) + 1e-3 * Math.cos(7.3 * i));
    for (let it = 0; it < 4; it++) {
      x = solveTridiagonal(a, b, c, x);
      // remove components along already-found states of nearly equal energy (degenerate or doublet pairs)
      states.forEach((v, j) => {
        if (Math.abs(energies[j] - E) > 1e-6 * scale) return;
        let d = 0; for (let i = 0; i < n; i++) d += v[i] * x[i] * dx;
        for (let i = 0; i < n; i++) x[i] -= d * v[i];
      });
      let s2 = 0; for (const v of x) s2 += v * v;
      const f = 1 / Math.sqrt(s2 * dx);
      for (let i = 0; i < n; i++) x[i] *= f;
    }
    let peak = 0; for (const v of x) peak = Math.max(peak, Math.abs(v));
    const sign = x.find(v => Math.abs(v) > 1e-6 * peak) > 0 ? 1 : -1;   // first lobe positive
    states.push(Float64Array.from(x, v => sign * v));
  });
  return { energies, states };
}

/**
 * Stationary scattering through a 1D potential V(x) (a function; constant outside [x0, x1]).
 * Integrates −(ħ²/2m)ψ'' + Vψ = Eψ from x1 leftwards, starting from the purely transmitted wave
 * e^{ik_R x}, and decomposes ψ at x0 into incoming and reflected waves, so tunnelling probabilities
 * down to ~1e−12 stay accurate (no exponentially growing transfer matrices). E must exceed V at
 * both ends. Returns { T, R } with T + R = 1 to the integration tolerance (a check in itself).
 */
export function transmission1D({ V, E, x0, x1, mass = 1, hbar = 1, rtol = 1e-11, atol = 1e-14 }) {
  if (typeof V !== 'function') throw new Error('transmission1D: V must be a function of x');
  const VL = V(x0), VR = V(x1);
  if (!(E > VL && E > VR)) throw new Error('transmission1D: E must lie above the potential at both ends');
  const c = 2 * mass / (hbar * hbar), kL = Math.sqrt(c * (E - VL)), kR = Math.sqrt(c * (E - VR));
  // state [Re ψ, Im ψ, Re ψ', Im ψ'] at x; integrate in s = x1 − x (forward in s)
  const f = (sv, y) => { const x = x1 - sv, q = c * (V(x) - E); return [-y[2], -y[3], -q * y[0], -q * y[1]]; };
  const y1 = [Math.cos(kR * x1), Math.sin(kR * x1), -kR * Math.sin(kR * x1), kR * Math.cos(kR * x1)];
  const { y } = dopri5(f, 0, y1, x1 - x0, { rtol, atol, h0: Math.min(1e-3, 0.01 / kR), hmax: Math.PI / (4 * Math.max(kL, kR)) });
  // ψ = A e^{ik x} + B e^{−ik x}: A = ½ e^{−ikx}(ψ + ψ'/(ik)), B = ½ e^{ikx}(ψ − ψ'/(ik))
  const [pr, pi, dr, di] = y, cx = Math.cos(kL * x0), sx = Math.sin(kL * x0);
  const qr = di / kL, qi = -dr / kL;                       // ψ'/(ik)
  const Ar = 0.5 * ((pr + qr) * cx + (pi + qi) * sx), Ai = 0.5 * ((pi + qi) * cx - (pr + qr) * sx);
  const Br = 0.5 * ((pr - qr) * cx - (pi - qi) * sx), Bi = 0.5 * ((pi - qi) * cx + (pr - qr) * sx);
  const A2 = Ar * Ar + Ai * Ai;
  return { T: (kR / kL) / A2, R: (Br * Br + Bi * Bi) / A2 };
}

/** (1/kT) ∫_{from}^{∞} f(E) e^{−(E − from)/kT} dE by Gauss–Legendre panels (f bounded, e.g. a transmission). */
export function thermalAverage(f, kT, { from = 0, panels = 40, n = 16, upto = 50 } = {}) {
  if (!(kT > 0)) throw new Error('thermalAverage: kT must be positive');
  let sum = 0;
  const w = upto * kT / panels;
  for (let p = 0; p < panels; p++) sum += gaussLegendre(E => f(E) * Math.exp(-(E - from) / kT), from + p * w, from + (p + 1) * w, n);
  return sum / kT;
}

/** ∫ |φ(k)|² T(k) dk for a Gaussian packet with |φ(k)|² ∝ exp(−(k − k0)²/(2σ_k²)) (σ_k = 1/(2σ_x)). */
export function packetTransmission(T, k0, sigmaK, { width = 8, n = 400 } = {}) {
  const a = k0 - width * sigmaK, b = k0 + width * sigmaK, norm = 1 / (Math.sqrt(2 * Math.PI) * sigmaK);
  let sum = 0;
  const panels = Math.ceil(n / 16), w = (b - a) / panels;
  for (let p = 0; p < panels; p++) sum += gaussLegendre(k => (k > 0 ? T(k) : 0) * norm * Math.exp(-((k - k0) ** 2) / (2 * sigmaK * sigmaK)), a + p * w, a + (p + 1) * w, 16);
  return sum;
}

