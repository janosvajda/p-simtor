// Incompressible (low-Mach) flow in two dimensions, built for live pages and verified against
// exact solutions and published data:
//   latticeBoltzmann2D  D2Q9 lattice Boltzmann (TRT or BGK), half-way bounce-back walls and obstacles
//                       (moving walls by Ladd's correction), velocity inflow, pressure or zero-gradient
//                       outflow, Guo body force, momentum-exchange forces on obstacles
//   vorticity2D         periodic Navier–Stokes in vorticity form, pseudo-spectral with 2/3 dealiasing
//                       and integrating-factor RK4 (viscosity exact), optional forcing and linear drag
//   bilinearSampler / advectParticles   tracers, streaklines and particle paths (RK4)
//   latticeUnits        lattice ↔ physical conversion (Re, Ma, τ)
//   poiseuilleProfile, couetteProfile, taylorGreen, cylinderStrouhal, GHIA_RE100   exact and
//                       published references used by the tests and for on-page comparisons
// Grid storage as in fields.mjs: flat arrays, axis 0 (x) fastest, index = i + nx·j.
// Lattice units: dx = dt = 1, c_s² = 1/3, kinematic viscosity ν = (τ − ½)/3, pressure p = c_s²ρ.

import { fft } from './numerics.mjs';

const CX = Int8Array.of(0, 1, 0, -1, 0, 1, -1, -1, 1);
const CY = Int8Array.of(0, 0, 1, 0, -1, 1, 1, -1, -1);
const W = Float64Array.of(4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36);
const OPP = Int8Array.of(0, 3, 4, 1, 2, 7, 8, 5, 6);
const SIDES = ['left', 'right', 'bottom', 'top'];
const TYPES = ['periodic', 'wall', 'inflow', 'outflow', 'pressure'];

const positiveInt = (v, name, who) => {
  if (!(Number.isSafeInteger(v) && v > 0)) throw new Error(`${who}: ${name} must be a positive integer`);
};
const positive = (v, name, who) => {
  if (!(Number.isFinite(v) && v > 0)) throw new Error(`${who}: ${name} must be positive and finite`);
};
const finite = (v, name, who) => {
  if (!Number.isFinite(v)) throw new Error(`${who}: ${name} must be finite`);
};

// ---------------------------------------------------------------- analytic references

/**
 * Plane Poiseuille flow between fixed walls at y = 0 and y = height, driven by a body force
 * per unit mass `force` (= −(1/ρ)dp/dx) with kinematic viscosity nu: u(y) = force·y(h − y)/(2ν).
 * Any consistent units. Mean velocity = (2/3)·max = force·h²/(12ν). Throws outside 0 ≤ y ≤ h.
 */
export function poiseuilleProfile(y, { height, force, nu }) {
  positive(height, 'height', 'poiseuilleProfile'); positive(nu, 'nu', 'poiseuilleProfile'); finite(force, 'force', 'poiseuilleProfile');
  if (!(y >= 0 && y <= height)) throw new Error('poiseuilleProfile: y must lie between the walls (0 ≤ y ≤ height)');
  return force * y * (height - y) / (2 * nu);
}

/** Plane Couette flow: wall at y = 0 at rest, wall at y = height moving at wallSpeed: u = U·y/h (steady, no pressure gradient). */
export function couetteProfile(y, { height, wallSpeed }) {
  positive(height, 'height', 'couetteProfile'); finite(wallSpeed, 'wallSpeed', 'couetteProfile');
  if (!(y >= 0 && y <= height)) throw new Error('couetteProfile: y must lie between the walls (0 ≤ y ≤ height)');
  return wallSpeed * y / height;
}

/**
 * Decaying Taylor–Green vortex, an exact solution of the 2D incompressible Navier–Stokes equations
 * on the periodic square of side 2π/k: u = U cos kx sin ky·F, v = −U sin kx cos ky·F,
 * ω = ∂v/∂x − ∂u/∂y = −2kU cos kx cos ky·F, p = −(ρU²/4)(cos 2kx + cos 2ky)·F², F = e^{−2νk²t}.
 * Kinetic energy (domain mean) U²/4·F². Returns { ux, uy, omega, p }.
 */
export function taylorGreen(x, y, t, { nu, k = 1, U = 1, rho = 1 }) {
  if (!(nu >= 0 && Number.isFinite(nu))) throw new Error('taylorGreen: nu must be ≥ 0 and finite');
  positive(k, 'k', 'taylorGreen'); finite(U, 'U', 'taylorGreen'); finite(x, 'x', 'taylorGreen'); finite(y, 'y', 'taylorGreen'); finite(t, 't', 'taylorGreen');
  const F = Math.exp(-2 * nu * k * k * t), cx = Math.cos(k * x), sx = Math.sin(k * x), cy = Math.cos(k * y), sy = Math.sin(k * y);
  return { ux: U * cx * sy * F, uy: -U * sx * cy * F, omega: -2 * k * U * cx * cy * F, p: -rho * U * U / 4 * (Math.cos(2 * k * x) + Math.cos(2 * k * y)) * F * F };
}

/**
 * Strouhal number St = f·D/U of laminar parallel vortex shedding behind a circular cylinder,
 * Williamson's fit St = −3.3265/Re + 0.1816 + 1.6×10⁻⁴·Re (C.H.K. Williamson, Annu. Rev. Fluid
 * Mech. 28, 477 (1996), eq. for 49 < Re < 178; St(100) ≈ 0.164). Unconfined cylinder; channel
 * blockage raises St. Throws outside the laminar shedding range.
 */
export function cylinderStrouhal(Re) {
  if (!(Re >= 49 && Re <= 178)) throw new Error('cylinderStrouhal: the fit is valid for laminar shedding, 49 ≤ Re ≤ 178');
  return -3.3265 / Re + 0.1816 + 1.6e-4 * Re;
}

/**
 * Lid-driven square cavity at Re = 100: centreline velocities of Ghia, Ghia & Shin,
 * J. Comput. Phys. 48, 387 (1982), Tables I and II (129×129 multigrid). Coordinates and
 * velocities normalised by the side L and lid speed U (lid at y = 1 moving in +x).
 * u: along the vertical centreline x = ½; v: along the horizontal centreline y = ½.
 */
export const GHIA_RE100 = Object.freeze({
  y: Object.freeze([1, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0]),
  u: Object.freeze([1, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.2109, -0.15662, -0.1015, -0.06434, -0.04775, -0.04192, -0.03717, 0]),
  x: Object.freeze([1, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0]),
  v: Object.freeze([0, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.1089, 0.10091, 0.09233, 0]),
});

// ---------------------------------------------------------------- unit conversion

/**
 * Lattice ↔ physical units for latticeBoltzmann2D. Give either physical scales
 * { length (m), velocity (m/s), viscosity (kinematic, m²/s) } or a Reynolds number { Re },
 * plus `cells` (lattice cells across the length) and `latticeVelocity` (the velocity scale in
 * lattice units, default 0.05). Re = U·L/ν is the same in both systems; Ma = u_lat·√3.
 * Returns { Re, Ma, nu, tau, dx, dt, cells, latticeVelocity, toPhysical: { length, time, velocity } }
 * (dx, dt in physical units, or 1/cells and latticeVelocity/cells of the unit length/velocity when only Re is given).
 * Throws if τ < minTau (default 0.51; τ → ½ is the stability limit of BGK/TRT, so raise cells or
 * latticeVelocity) or Ma ≥ 0.3 (compressibility error ~Ma²).
 */
export function latticeUnits({ length, velocity, viscosity, Re, cells, latticeVelocity = 0.05, minTau = 0.51 }) {
  const who = 'latticeUnits';
  positiveInt(cells, 'cells', who); positive(latticeVelocity, 'latticeVelocity', who);
  const physical = length !== undefined || velocity !== undefined || viscosity !== undefined;
  let L = 1, U = 1, Rey = Re;
  if (physical) {
    if (Re !== undefined) throw new Error(`${who}: give either physical scales or Re, not both`);
    positive(length, 'length', who); positive(velocity, 'velocity', who); positive(viscosity, 'viscosity', who);
    L = length; U = velocity; Rey = U * L / viscosity;
  } else positive(Re, 'Re', who);
  const nu = latticeVelocity * cells / Rey, tau = 3 * nu + 0.5, Ma = latticeVelocity * Math.sqrt(3);
  if (!(tau >= minTau && tau > 0.5)) throw new Error(`${who}: τ = ${tau} < minTau = ${minTau} is unstable for all but very smooth flows; raise cells or latticeVelocity`);
  if (!(Ma < 0.3)) throw new Error(`${who}: Ma = ${Ma.toFixed(3)} ≥ 0.3; lower latticeVelocity (compressibility error grows as Ma²)`);
  const dx = L / cells, dt = dx * latticeVelocity / U;
  return {
    Re: Rey, Ma, nu, tau, dx, dt, cells, latticeVelocity,
    toPhysical: { length: l => l * dx, time: steps => steps * dt, velocity: u => u * dx / dt },
  };
}

// ---------------------------------------------------------------- lattice Boltzmann D2Q9

function normaliseSide(spec, side, who) {
  const s = typeof spec === 'string' ? { type: spec } : { ...spec };
  if (!TYPES.includes(s.type)) throw new Error(`${who}: ${side} boundary type must be one of ${TYPES.join(', ')}`);
  if (s.type === 'inflow' && s.u === undefined) throw new Error(`${who}: ${side} inflow needs u ([ux, uy] or a function of the boundary cell index)`);
  if (s.type === 'pressure') { s.rho ??= 1; positive(s.rho, `${side} rho`, who); }
  return s;
}

/**
 * Two-dimensional lattice Boltzmann method (D2Q9) for weakly compressible, effectively
 * incompressible flow at low Mach number, all in lattice units (dx = dt = 1).
 *   nx, ny        grid cells; fluid cells at x = i, y = j
 *   tau | nu      relaxation time τ > ½, or kinematic viscosity ν = (τ − ½)/3
 *   collision     'trt' (default; two-relaxation-time, magic Λ = 3/16 puts half-way bounce-back
 *                 walls exactly mid-link, so Poiseuille flow is exact) or 'bgk'
 *   solid         Uint8Array(nx·ny): 0 fluid, 1…255 obstacle id (force by momentum exchange per id)
 *   boundaries    { left, right, bottom, top }: 'periodic' (both sides of an axis), 'wall' or
 *                 { type: 'wall', u: [ux, uy] } (half-way bounce-back, Ladd's moving-wall term; walls
 *                 lie half a cell outside the edge cells), { type: 'inflow', u } (velocity bounce-back:
 *                 the same moving-wall rule with a normal component, mass influx ρ·u·n exactly; u is
 *                 [ux, uy] or a function of the boundary cell index), { type: 'pressure', rho } (density
 *                 fixed, default 1, by Guo's non-equilibrium extrapolation: the recommended outlet; it
 *                 absorbs the mass the inflow brings) or 'outflow' (zero-gradient copy; lets density drift).
 *                 Default: periodic everywhere. Corners: a wall wins.
 *                 Impulsive starts send acoustic waves (forces ring for a few D/U): start from potential
 *                 flow or ramp the inflow with setBoundary.
 *   force         [Fx, Fy] uniform body force density (Guo forcing; second order)
 *   density, velocity   initial ρ(i, j) and [ux, uy](i, j) functions (default ρ = 1, u = 0); the
 *                 populations start at equilibrium.
 *   machLimit     step() throws if max |u|·√3 exceeds it (default 0.3) or the state is not finite.
 * Stable for τ > ½ in principle; below τ ≈ 0.51 only smooth, well-resolved flows survive (grid
 * Reynolds number u/ν ≲ 10). Mass Σρ is conserved to round-off with periodic sides, walls (also
 * tangentially moving ones) and obstacles. ~40 M cell updates/s in V8 (a 200×100 grid at ~2000 steps/s). Returns { nx, ny, nu, tau, t, rho, ux, uy, step(n), mass(), force(id), vorticity(out),
 * setBoundary(side, spec), maxMach } (rho, ux, uy are live Float64Arrays of the current state).
 */
export function latticeBoltzmann2D({ nx, ny, tau, nu, collision = 'trt', magic = 3 / 16, solid, boundaries = {}, force = [0, 0], density, velocity, machLimit = 0.3 }) {
  const who = 'latticeBoltzmann2D';
  positiveInt(nx, 'nx', who); positiveInt(ny, 'ny', who);
  if ((tau === undefined) === (nu === undefined)) throw new Error(`${who}: give exactly one of tau and nu`);
  const tauP = tau ?? 3 * nu + 0.5;
  if (!(Number.isFinite(tauP) && tauP > 0.5)) throw new Error(`${who}: τ must exceed ½ (ν > 0); got τ = ${tauP}`);
  if (collision !== 'trt' && collision !== 'bgk') throw new Error(`${who}: collision must be 'trt' or 'bgk'`);
  positive(magic, 'magic', who); positive(machLimit, 'machLimit', who);
  if (!(force?.length === 2 && force.every(Number.isFinite))) throw new Error(`${who}: force must be [Fx, Fy] (finite)`);
  const N = nx * ny;
  if (solid && solid.length !== N) throw new Error(`${who}: solid mask length must equal nx·ny`);
  const mask = solid ? Uint8Array.from(solid) : new Uint8Array(N);
  const sides = SIDES.map(s => normaliseSide(boundaries[s] ?? 'periodic', s, who));
  if ((sides[0].type === 'periodic') !== (sides[1].type === 'periodic') || (sides[2].type === 'periodic') !== (sides[3].type === 'periodic')) {
    throw new Error(`${who}: 'periodic' must be set on both sides of an axis`);
  }
  const omP = 1 / tauP, omM = collision === 'bgk' ? omP : 1 / (magic / (tauP - 0.5) + 0.5);
  const [Fx, Fy] = force;
  // per-side velocity tables (transverse cell index → velocity)
  const sideU = sides.map((s, k) => {
    const len = k < 2 ? ny : nx, ux = new Float64Array(len), uy = new Float64Array(len);
    return { ux, uy, len, s };
  });
  const fillSide = k => {
    const { ux, uy, len, s } = sideU[k];
    for (let p = 0; p < len; p++) {
      const u = typeof s.u === 'function' ? s.u(p) : s.u ?? [0, 0];
      if (!(u?.length === 2 && Number.isFinite(u[0]) && Number.isFinite(u[1]))) throw new Error(`${who}: ${SIDES[k]} boundary velocity must be [ux, uy] (finite)`);
      ux[p] = u[0]; uy[p] = u[1];
    }
  };
  for (let k = 0; k < 4; k++) fillSide(k);

  // pull table: ≥ 0 source cell; −2 − side for a domain boundary; ≤ −10 for a solid neighbour (−10 − index)
  const src = new Int32Array(9 * N), cellKind = new Uint8Array(N);    // cellKind: 0 bulk, 1 boundary, 2 solid
  const buildSources = () => {
    for (let q = 0; q < 9; q++) {
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const x = i + nx * j;
          let si = i - CX[q], sj = j - CY[q], ySide = -1, xSide = -1;
          if (sj < 0 || sj >= ny) {
            const side = sj < 0 ? 2 : 3, type = sides[side].type;
            if (type === 'periodic') sj = (sj + ny) % ny;
            else if (type === 'outflow') sj = j;
            else ySide = side;
          }
          if (si < 0 || si >= nx) {
            const side = si < 0 ? 0 : 1, type = sides[side].type;
            if (type === 'periodic') si = (si + nx) % nx;
            else if (type === 'outflow') si = i;
            else xSide = side;
          }
          let code;
          if (ySide >= 0 && xSide >= 0) code = -2 - (sides[xSide].type === 'wall' && sides[ySide].type !== 'wall' ? xSide : ySide);
          else if (ySide >= 0) code = -2 - ySide;
          else if (xSide >= 0) code = -2 - xSide;
          else { const s = si + nx * sj; code = mask[s] ? -10 - s : s; }
          src[q * N + x] = code;
        }
      }
    }
    for (let x = 0; x < N; x++) {
      let kind = mask[x] ? 2 : 0;
      for (let q = 0; q < 9 && kind === 0; q++) if (src[q * N + x] < 0) kind = 1;
      cellKind[x] = kind;
    }
  };
  buildSources();

  let f = new Float64Array(9 * N), fn = new Float64Array(9 * N);
  const rho = new Float64Array(N), ux = new Float64Array(N), uy = new Float64Array(N);
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const x = i + nx * j;
      if (mask[x]) { rho[x] = 1; continue; }
      const r = density ? density(i, j) : 1, u = velocity ? velocity(i, j) : [0, 0];
      if (!(r > 0 && Number.isFinite(r) && Number.isFinite(u[0]) && Number.isFinite(u[1]))) throw new Error(`${who}: initial density must be positive and velocity finite at (${i}, ${j})`);
      rho[x] = r; ux[x] = u[0]; uy[x] = u[1];
      const usq = 1.5 * (u[0] * u[0] + u[1] * u[1]);
      for (let q = 0; q < 9; q++) { const cu = 3 * (CX[q] * u[0] + CY[q] * u[1]); f[q * N + x] = W[q] * r * (1 + cu + 0.5 * cu * cu - usq); }
    }
  }

  const forces = new Float64Array(2 * 256);
  const g = new Float64Array(9);
  const hp = 1 - 0.5 * omP, hm = 1 - 0.5 * omM, hp9 = hp / 9, hp36 = hp / 36, hm3_9 = hm / 3, hm3_36 = hm / 12;
  // incoming population q at fluid cell x whose pull source is an obstacle (s ≤ −10) or a domain side
  const edge = (fo, q, s, x) => {
    if (s <= -10) {
      const v = fo[OPP[q] * N + x], id = mask[-10 - s];
      forces[2 * id] -= 2 * CX[q] * v; forces[2 * id + 1] -= 2 * CY[q] * v;
      return v;
    }
    const k = -2 - s, b = sideU[k], p = k < 2 ? (x / nx) | 0 : x % nx, type = b.s.type;
    if (type === 'wall' || type === 'inflow') return fo[OPP[q] * N + x] + 6 * W[q] * rho[x] * (CX[q] * b.ux[p] + CY[q] * b.uy[p]);
    // pressure: Guo's non-equilibrium extrapolation, feq(ρ_out, u_x) + [f*_q(x) − feq(ρ_x, u_x)]
    const ax = ux[x], ay = uy[x], ca = 3 * (CX[q] * ax + CY[q] * ay);
    return W[q] * (b.s.rho - rho[x]) * (1 + ca + 0.5 * ca * ca - 1.5 * (ax * ax + ay * ay)) + fo[q * N + x];
  };
  const sim = {
    nx, ny, tau: tauP, nu: (tauP - 0.5) / 3, t: 0, rho, ux, uy, maxMach: 0, solid: mask,
    /** Advance n lattice time steps (stream, boundaries, collide). */
    step(n = 1) {
      positiveInt(n, 'n', who);
      for (let it = 0; it < n; it++) {
        forces.fill(0);
        const fo = f, fw = fn;                                    // locals: closure lets are slow in hot loops
        let maxU2 = 0, sumU2 = 0;
        for (let x = 0; x < N; x++) {
          const kind = cellKind[x];
          if (kind === 2) continue;                                // solid
          let g0, g1, g2, g3, g4, g5, g6, g7, g8;
          if (kind === 0) {                                        // bulk: all sources are fluid cells
            g0 = fo[x]; g1 = fo[N + src[N + x]]; g2 = fo[2 * N + src[2 * N + x]]; g3 = fo[3 * N + src[3 * N + x]];
            g4 = fo[4 * N + src[4 * N + x]]; g5 = fo[5 * N + src[5 * N + x]]; g6 = fo[6 * N + src[6 * N + x]];
            g7 = fo[7 * N + src[7 * N + x]]; g8 = fo[8 * N + src[8 * N + x]];
          } else {
            for (let q = 0; q < 9; q++) { const s = src[q * N + x]; g[q] = s >= 0 ? fo[q * N + s] : edge(fo, q, s, x); }
            g0 = g[0]; g1 = g[1]; g2 = g[2]; g3 = g[3]; g4 = g[4]; g5 = g[5]; g6 = g[6]; g7 = g[7]; g8 = g[8];
          }
          const r = g0 + g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8, mx = g1 - g3 + g5 - g6 - g7 + g8, my = g2 - g4 + g5 + g6 - g7 - g8;
          const vx = (mx + 0.5 * Fx) / r, vy = (my + 0.5 * Fy) / r, u2 = vx * vx + vy * vy, usq = 1.5 * u2;
          if (u2 > maxU2) maxU2 = u2;
          sumU2 += u2;
          rho[x] = r; ux[x] = vx; uy[x] = vy;
          // TRT: symmetric/antisymmetric parts of each pair (q, q̄) relax with ω⁺/ω⁻; Guo source split alike
          const uF = vx * Fx + vy * Fy;
          fw[x] = g0 - omP * (g0 - W[0] * r * (1 - usq)) - hp * 3 * W[0] * uF;
          // pairs (1,3) (2,4) (5,7) (6,8), unrolled: c·u, c·F of the first member
          const ws = r / 9, wd = r / 36;
          let cu, cF, dp, dm, sp, sm;
          cu = vx; cF = Fx;
          dp = 0.5 * (g1 + g3) - ws * (1 + 4.5 * cu * cu - usq); dm = 0.5 * (g1 - g3) - 3 * ws * cu;
          sp = hp9 * (9 * cu * cF - 3 * uF); sm = hm3_9 * cF;
          fw[N + x] = g1 - omP * dp - omM * dm + sp + sm; fw[3 * N + x] = g3 - omP * dp + omM * dm + sp - sm;
          cu = vy; cF = Fy;
          dp = 0.5 * (g2 + g4) - ws * (1 + 4.5 * cu * cu - usq); dm = 0.5 * (g2 - g4) - 3 * ws * cu;
          sp = hp9 * (9 * cu * cF - 3 * uF); sm = hm3_9 * cF;
          fw[2 * N + x] = g2 - omP * dp - omM * dm + sp + sm; fw[4 * N + x] = g4 - omP * dp + omM * dm + sp - sm;
          cu = vx + vy; cF = Fx + Fy;
          dp = 0.5 * (g5 + g7) - wd * (1 + 4.5 * cu * cu - usq); dm = 0.5 * (g5 - g7) - 3 * wd * cu;
          sp = hp36 * (9 * cu * cF - 3 * uF); sm = hm3_36 * cF;
          fw[5 * N + x] = g5 - omP * dp - omM * dm + sp + sm; fw[7 * N + x] = g7 - omP * dp + omM * dm + sp - sm;
          cu = vy - vx; cF = Fy - Fx;
          dp = 0.5 * (g6 + g8) - wd * (1 + 4.5 * cu * cu - usq); dm = 0.5 * (g6 - g8) - 3 * wd * cu;
          sp = hp36 * (9 * cu * cF - 3 * uF); sm = hm3_36 * cF;
          fw[6 * N + x] = g6 - omP * dp - omM * dm + sp + sm; fw[8 * N + x] = g8 - omP * dp + omM * dm + sp - sm;
        }
        const tmp = f; f = fn; fn = tmp;
        sim.t++;
        sim.maxMach = Math.sqrt(3 * maxU2);
        if (!Number.isFinite(sumU2)) throw new Error(`${who}: state became non-finite at t = ${sim.t} (unstable: raise τ or resolution, lower velocities)`);
        if (sim.maxMach > machLimit) throw new Error(`${who}: Mach number ${sim.maxMach.toFixed(3)} exceeds machLimit ${machLimit} at t = ${sim.t}; lower the driving velocity or force`);
      }
      return sim;
    },
    /** Total mass Σρ over fluid cells (lattice units). */
    mass() { let m = 0; for (let x = 0; x < N; x++) if (!mask[x]) m += rho[x]; return m; },
    /** Force [Fx, Fy] on obstacle `id` during the last step (momentum exchange, lattice units per unit depth). */
    force(id = 1) {
      if (!(Number.isInteger(id) && id >= 1 && id <= 255)) throw new Error(`${who}: obstacle id must be an integer 1…255`);
      return [forces[2 * id], forces[2 * id + 1]];
    },
    /** Vorticity ∂uy/∂x − ∂ux/∂y by central differences (solid neighbours count as u = 0, open edges one-sided). */
    vorticity(out = new Float64Array(N)) {
      if (out.length !== N) throw new Error(`${who}: vorticity output must have nx·ny entries`);
      const px = sides[0].type === 'periodic', py = sides[2].type === 'periodic';
      const at = (arr, i, j, x) => {
        if ((i < 0 || i >= nx) && !px) return arr[x];
        if ((j < 0 || j >= ny) && !py) return arr[x];
        const s = (i + nx) % nx + nx * ((j + ny) % ny); return mask[s] ? 0 : arr[s];
      };
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const x = i + nx * j;
          out[x] = mask[x] ? 0 : 0.5 * (at(uy, i + 1, j, x) - at(uy, i - 1, j, x) - at(ux, i, j + 1, x) + at(ux, i, j - 1, x));
        }
      }
      return out;
    },
    /** Replace a boundary (same periodicity), e.g. to ramp an inflow or start a lid. */
    setBoundary(side, spec) {
      const k = SIDES.indexOf(side);
      if (k < 0) throw new Error(`${who}: side must be one of ${SIDES.join(', ')}`);
      const s = normaliseSide(spec, side, who);
      if ((s.type === 'periodic') !== (sides[k].type === 'periodic')) throw new Error(`${who}: setBoundary cannot change periodicity`);
      sides[k] = s; sideU[k].s = s; fillSide(k); buildSources();
    },
  };
  return sim;
}

// ---------------------------------------------------------------- pseudo-spectral vorticity

/** In-place 2D FFT on row-major (x fastest) arrays using scratch columns. */
function fft2(re, im, nx, ny, colRe, colIm, inverse) {
  for (let j = 0; j < ny; j++) fft(re.subarray(j * nx, (j + 1) * nx), im.subarray(j * nx, (j + 1) * nx), inverse);
  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) { colRe[j] = re[i + nx * j]; colIm[j] = im[i + nx * j]; }
    fft(colRe, colIm, inverse);
    for (let j = 0; j < ny; j++) { re[i + nx * j] = colRe[j]; im[i + nx * j] = colIm[j]; }
  }
}

/**
 * Incompressible 2D Navier–Stokes on the periodic box [0, Lx) × [0, Ly) in vorticity form
 *   ∂ω/∂t + u·∇ω = ν∇²ω − μω + forcing,   u = (∂ψ/∂y, −∂ψ/∂x),  ∇²ψ = −ω,
 * pseudo-spectral (power-of-two nx, ny), 2/3-rule dealiasing (also applied to ω₀), Lawson
 * integrating-factor RK4 (viscous and drag terms exact; 4th order in dt for the nonlinear part).
 * The mean flow is zero (ω carries no k = 0 information). forcing: optional steady vorticity source
 * on the grid. With ν = μ = 0 the dealiased (Galerkin) system conserves energy and enstrophy exactly
 * in space; RK4 adds an O(dt⁵) drift per step. step() throws if the advective CFL number
 * max|u|·dt/dx exceeds cflMax (default 1) or the state is not finite.
 * Returns { t, step(n), vorticity(), velocity(), streamfunction(), energy(), enstrophy(), cfl }:
 * energy = ½⟨|u|²⟩ and enstrophy = ½⟨ω²⟩ are domain means.
 */
export function vorticity2D({ omega0, shape, length = [2 * Math.PI, 2 * Math.PI], nu = 0, drag = 0, forcing, dt, dealias = true, cflMax = 1 }) {
  const who = 'vorticity2D';
  if (!(shape?.length === 2)) throw new Error(`${who}: shape must be [nx, ny]`);
  const [nx, ny] = shape;
  for (const n of shape) if (!(Number.isSafeInteger(n) && n >= 4 && !(n & (n - 1)))) throw new Error(`${who}: nx and ny must be powers of two ≥ 4`);
  const [Lx, Ly] = length;
  positive(Lx, 'Lx', who); positive(Ly, 'Ly', who); positive(dt, 'dt', who); positive(cflMax, 'cflMax', who);
  if (!(nu >= 0 && Number.isFinite(nu) && drag >= 0 && Number.isFinite(drag))) throw new Error(`${who}: nu and drag must be ≥ 0 and finite`);
  const N = nx * ny;
  if (!(omega0?.length === N) || !Array.prototype.every.call(omega0, Number.isFinite)) throw new Error(`${who}: omega0 must hold nx·ny finite values`);
  if (forcing && !(forcing.length === N && Array.prototype.every.call(forcing, Number.isFinite))) throw new Error(`${who}: forcing must hold nx·ny finite values`);
  const kx = new Float64Array(N), ky = new Float64Array(N), k2 = new Float64Array(N), keep = new Uint8Array(N), E = new Float64Array(N);
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const x = i + nx * j, mi = i < nx / 2 ? i : i - nx, mj = j < ny / 2 ? j : j - ny;
      const kxv = 2 * Math.PI * mi / Lx, kyv = 2 * Math.PI * mj / Ly;
      k2[x] = kxv * kxv + kyv * kyv;
      kx[x] = i === nx / 2 ? 0 : kxv; ky[x] = j === ny / 2 ? 0 : kyv;     // no derivative of the Nyquist mode
      keep[x] = !dealias || (Math.abs(mi) < nx / 3 && Math.abs(mj) < ny / 3) ? 1 : 0;
      E[x] = Math.exp(-(nu * k2[x] + drag) * dt / 2);
    }
  }
  const colRe = new Float64Array(ny), colIm = new Float64Array(ny);
  const wr = Float64Array.from(omega0), wi = new Float64Array(N);
  fft2(wr, wi, nx, ny, colRe, colIm, false);
  for (let x = 0; x < N; x++) if (!keep[x]) { wr[x] = 0; wi[x] = 0; }
  wr[0] = 0; wi[0] = 0;                                    // mean vorticity is zero on a periodic box
  let fr = null, fi = null;
  if (forcing) {
    fr = Float64Array.from(forcing); fi = new Float64Array(N);
    fft2(fr, fi, nx, ny, colRe, colIm, false);
    for (let x = 0; x < N; x++) if (!keep[x]) { fr[x] = 0; fi[x] = 0; }
    fr[0] = 0; fi[0] = 0;
  }
  const ar = new Float64Array(N), ai = new Float64Array(N), br = new Float64Array(N), bi = new Float64Array(N);
  const cr = new Float64Array(N), ci = new Float64Array(N), dr = new Float64Array(N), di = new Float64Array(N);
  const sr = new Float64Array(N), si = new Float64Array(N);
  const pr = new Float64Array(N), pi = new Float64Array(N), qr = new Float64Array(N), qi = new Float64Array(N);
  let maxU = 0;
  /** Nonlinear term −u·∇ω (+ forcing) of spectral vorticity (hr, hi), written to (outR, outI). */
  const nonlinear = (hr, hi, outR, outI) => {
    for (let x = 0; x < N; x++) {
      const inv = k2[x] > 0 ? 1 / k2[x] : 0, psr = hr[x] * inv, psi = hi[x] * inv;
      // û = i ky ψ̂, v̂ = −i kx ψ̂ packed as û + i v̂; ω̂_x + i ω̂_y likewise
      pr[x] = kx[x] * psr - ky[x] * psi; pi[x] = ky[x] * psr + kx[x] * psi;
      qr[x] = -kx[x] * hi[x] - ky[x] * hr[x]; qi[x] = kx[x] * hr[x] - ky[x] * hi[x];
    }
    fft2(pr, pi, nx, ny, colRe, colIm, true); fft2(qr, qi, nx, ny, colRe, colIm, true);
    let m = 0;
    for (let x = 0; x < N; x++) {
      const u = pr[x], v = pi[x]; m = Math.max(m, Math.abs(u) / (Lx / nx) + Math.abs(v) / (Ly / ny));
      outR[x] = -(u * qr[x] + v * qi[x]); outI[x] = 0;
    }
    maxU = Math.max(maxU, m);
    fft2(outR, outI, nx, ny, colRe, colIm, false);
    for (let x = 0; x < N; x++) {
      if (!keep[x]) { outR[x] = 0; outI[x] = 0; } else if (fr) { outR[x] += fr[x]; outI[x] += fi[x]; }
    }
    outR[0] = 0; outI[0] = 0;
  };
  const toPhysical = (hr, hi, scale) => {
    const re = new Float64Array(N), im = new Float64Array(N);
    for (let x = 0; x < N; x++) { const s = scale(x); re[x] = hr[x] * s; im[x] = hi[x] * s; }
    fft2(re, im, nx, ny, colRe, colIm, true);
    return re;
  };
  const sim = {
    t: 0, cfl: 0, nx, ny, length: [Lx, Ly],
    /** Advance n steps of dt (Lawson integrating-factor RK4). */
    step(n = 1) {
      positiveInt(n, 'n', who);
      for (let it = 0; it < n; it++) {
        maxU = 0;
        nonlinear(wr, wi, ar, ai);
        for (let x = 0; x < N; x++) { sr[x] = E[x] * (wr[x] + 0.5 * dt * ar[x]); si[x] = E[x] * (wi[x] + 0.5 * dt * ai[x]); }
        nonlinear(sr, si, br, bi);
        for (let x = 0; x < N; x++) { sr[x] = E[x] * wr[x] + 0.5 * dt * br[x]; si[x] = E[x] * wi[x] + 0.5 * dt * bi[x]; }
        nonlinear(sr, si, cr, ci);
        for (let x = 0; x < N; x++) { const e2 = E[x] * E[x]; sr[x] = e2 * wr[x] + dt * E[x] * cr[x]; si[x] = e2 * wi[x] + dt * E[x] * ci[x]; }
        nonlinear(sr, si, dr, di);
        for (let x = 0; x < N; x++) {
          const e = E[x], e2 = e * e;
          wr[x] = e2 * wr[x] + dt / 6 * (e2 * ar[x] + 2 * e * (br[x] + cr[x]) + dr[x]);
          wi[x] = e2 * wi[x] + dt / 6 * (e2 * ai[x] + 2 * e * (bi[x] + ci[x]) + di[x]);
        }
        sim.t += dt;
        sim.cfl = maxU * dt;
        if (!Number.isFinite(sim.cfl) || !Number.isFinite(sim.enstrophy())) throw new Error(`${who}: state became non-finite at t = ${sim.t}`);
        if (sim.cfl > cflMax) throw new Error(`${who}: CFL number ${sim.cfl.toFixed(3)} exceeds cflMax ${cflMax}; reduce dt`);
      }
      return sim;
    },
    /** Vorticity ω on the grid (new Float64Array). */
    vorticity() { return toPhysical(wr, wi, () => 1); },
    /** Streamfunction ψ with ∇²ψ = −ω (zero mean). */
    streamfunction() { return toPhysical(wr, wi, x => (k2[x] > 0 ? 1 / k2[x] : 0)); },
    /** Velocity { ux, uy } on the grid (new Float64Arrays). */
    velocity() {
      const re = new Float64Array(N), im = new Float64Array(N);
      for (let x = 0; x < N; x++) {
        const inv = k2[x] > 0 ? 1 / k2[x] : 0, psr = wr[x] * inv, psi = wi[x] * inv;
        re[x] = kx[x] * psr - ky[x] * psi; im[x] = ky[x] * psr + kx[x] * psi;
      }
      fft2(re, im, nx, ny, colRe, colIm, true);
      return { ux: re, uy: im };
    },
    /** Kinetic energy per unit area ½⟨|u|²⟩ (Parseval, exact for the truncated field). */
    energy() { let s = 0; for (let x = 1; x < N; x++) if (k2[x] > 0) s += (wr[x] * wr[x] + wi[x] * wi[x]) / k2[x]; return 0.5 * s / (N * N); },
    /** Enstrophy per unit area ½⟨ω²⟩. */
    enstrophy() { let s = 0; for (let x = 0; x < N; x++) s += wr[x] * wr[x] + wi[x] * wi[x]; return 0.5 * s / (N * N); },
  };
  return sim;
}

// ---------------------------------------------------------------- tracers

/**
 * Bilinear velocity interpolation on a grid with nodes at (x0 + i·dx, y0 + j·dy) (flat arrays,
 * x fastest). periodic: [px, py] wraps that axis with period n·d; otherwise positions are clamped
 * to the node range (zero-gradient). Exact for velocity fields linear in x and y. solid (optional):
 * cells with solid[i] ≠ 0 count as u = 0. Returns velocity(x, y, t, out) for advectParticles
 * (t ignored: the grid is the current frozen field).
 */
export function bilinearSampler({ ux, uy, nx, ny, spacing = 1, origin = [0, 0], periodic = [false, false], solid }) {
  const who = 'bilinearSampler';
  positiveInt(nx, 'nx', who); positiveInt(ny, 'ny', who);
  if (!(ux?.length === nx * ny && uy?.length === nx * ny)) throw new Error(`${who}: ux and uy must have nx·ny entries`);
  if (solid && solid.length !== nx * ny) throw new Error(`${who}: solid mask must have nx·ny entries`);
  const [hx, hy] = typeof spacing === 'number' ? [spacing, spacing] : spacing;
  positive(hx, 'spacing', who); positive(hy, 'spacing', who);
  const [x0, y0] = origin;
  finite(x0, 'origin', who); finite(y0, 'origin', who);
  const axis = (s, n, per) => {
    let i0, w;
    if (per) { const m = ((s % n) + n) % n; i0 = Math.floor(m); w = m - i0; return [i0 % n, (i0 + 1) % n, w]; }
    if (n === 1) return [0, 0, 0];
    const c = Math.min(Math.max(s, 0), n - 1);
    i0 = Math.min(Math.floor(c), n - 2); w = c - i0;
    return [i0, i0 + 1, w];
  };
  return (x, y, _t, out = new Float64Array(2)) => {
    const [i0, i1, wx] = axis((x - x0) / hx, nx, periodic[0]), [j0, j1, wy] = axis((y - y0) / hy, ny, periodic[1]);
    const a = i0 + nx * j0, b = i1 + nx * j0, c = i0 + nx * j1, d = i1 + nx * j1;
    const w00 = (1 - wx) * (1 - wy), w10 = wx * (1 - wy), w01 = (1 - wx) * wy, w11 = wx * wy;
    const val = (arr, k) => (solid && solid[k] ? 0 : arr[k]);
    out[0] = w00 * val(ux, a) + w10 * val(ux, b) + w01 * val(ux, c) + w11 * val(ux, d);
    out[1] = w00 * val(uy, a) + w10 * val(uy, b) + w01 * val(uy, c) + w11 * val(uy, d);
    return out;
  };
}

/**
 * Advance tracer particles dx/dt = u(x, t) by classical RK4 (4th order), in place.
 * xy: flat [x0, y0, x1, y1, …] (Float64Array). velocity(x, y, t, out) writes [ux, uy] into out
 * (e.g. from bilinearSampler or an analytic field). Takes `steps` steps of dt from time t.
 * Throws if a velocity is not finite. Returns xy.
 */
export function advectParticles(xy, velocity, dt, { steps = 1, t = 0 } = {}) {
  const who = 'advectParticles';
  if (!(xy?.length % 2 === 0)) throw new Error(`${who}: xy must be a flat [x, y, …] array of even length`);
  if (typeof velocity !== 'function') throw new Error(`${who}: velocity must be a function (x, y, t, out)`);
  finite(dt, 'dt', who); positiveInt(steps, 'steps', who);
  const k = new Float64Array(2);
  for (let s = 0; s < steps; s++) {
    const t0 = t + s * dt;
    for (let p = 0; p < xy.length; p += 2) {
      const x = xy[p], y = xy[p + 1];
      velocity(x, y, t0, k); const k1x = k[0], k1y = k[1];
      velocity(x + 0.5 * dt * k1x, y + 0.5 * dt * k1y, t0 + 0.5 * dt, k); const k2x = k[0], k2y = k[1];
      velocity(x + 0.5 * dt * k2x, y + 0.5 * dt * k2y, t0 + 0.5 * dt, k); const k3x = k[0], k3y = k[1];
      velocity(x + dt * k3x, y + dt * k3y, t0 + dt, k);
      const nxp = x + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k[0]), nyp = y + dt / 6 * (k1y + 2 * k2y + 2 * k3y + k[1]);
      if (!(Number.isFinite(nxp) && Number.isFinite(nyp))) throw new Error(`${who}: non-finite velocity at particle ${p / 2}`);
      xy[p] = nxp; xy[p + 1] = nyp;
    }
  }
  return xy;
}
