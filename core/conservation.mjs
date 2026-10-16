// Hyperbolic conservation laws q_t + ∇·F(q) = S by Godunov-type finite volumes (Toro, Riemann
// Solvers and Numerical Methods for Fluid Dynamics, 3rd ed. 2009; LeVeque, Finite Volume Methods
// for Hyperbolic Problems 2002):
//   finiteVolume        1D or 2D cell averages: MUSCL reconstruction of primitive variables with a
//                       slope limiter, an approximate or exact Riemann solver at every face, SSP-RK2/3
//                       in time (method of lines; 2D dimension-by-dimension, unsplit), CFL control,
//                       boundary-flux bookkeeping so conservation is checkable to round-off
//   systems             linearAdvection, inviscidBurgers, euler (ideal gas, γ), shallowWater (with
//                       bathymetry: well-balanced hydrostatic reconstruction, Audusse et al. 2004)
//   exact solutions     burgersExactRiemann, eulerExactRiemann (Toro ch. 4), shallowWaterExactRiemann
//                       (Toro, Shock-Capturing Methods for Free-Surface Shallow Flows 2001, ch. 5; wet and dry)
// Storage as fields.mjs: one Float64Array per conserved variable, axis 0 fastest (index i + nx·j).
// Units are the caller's; all formulas are dimensionally consistent (g in length/time², γ dimensionless).

const fail = (who, msg) => { throw new Error(`${who}: ${msg}`); };

// ------------------------------------------------------------------ limiters (undivided slopes)

const LIMITERS = {
  first: () => 0,                                                            // Godunov, first order
  none: (a, b) => 0.5 * (a + b),                                             // unlimited central (Fromm): oscillates at jumps
  minmod: (a, b) => (a * b <= 0 ? 0 : Math.abs(a) < Math.abs(b) ? a : b),
  vanleer: (a, b) => (a * b <= 0 ? 0 : 2 * a * b / (a + b)),
  mc: (a, b) => {                                                            // monotonised central (van Leer 1977)
    if (a * b <= 0) return 0;
    const m = Math.min(2 * Math.abs(a), 2 * Math.abs(b), 0.5 * Math.abs(a + b));
    return a > 0 ? m : -m;
  },
  superbee: (a, b) => {
    if (a * b <= 0) return 0;
    const A = Math.abs(a), B = Math.abs(b), m = Math.max(Math.min(2 * A, B), Math.min(A, 2 * B));
    return a > 0 ? m : -m;
  },
};

// ------------------------------------------------------------------ generic Riemann fluxes

/** HLL flux (Harten, Lax & van Leer 1983) from the system's flux, conversion and wave-speed estimate. */
function makeHLL(sys) {
  const FL = new Float64Array(sys.nvar), FR = new Float64Array(sys.nvar), UL = new Float64Array(sys.nvar), UR = new Float64Array(sys.nvar), S = new Float64Array(2);
  return (pL, pR, ax, out) => {
    sys.waveSpeeds(pL, pR, ax, S);
    const sL = S[0], sR = S[1];
    if (sL >= 0) { sys.flux(pL, ax, out); return; }
    if (sR <= 0) { sys.flux(pR, ax, out); return; }
    sys.flux(pL, ax, FL); sys.flux(pR, ax, FR); sys.primToCons(pL, UL); sys.primToCons(pR, UR);
    const inv = 1 / (sR - sL);
    for (let v = 0; v < sys.nvar; v++) out[v] = (sR * FL[v] - sL * FR[v] + sL * sR * (UR[v] - UL[v])) * inv;
  };
}
/** Rusanov / local Lax–Friedrichs: the most dissipative, most robust choice. */
function makeRusanov(sys) {
  const FL = new Float64Array(sys.nvar), FR = new Float64Array(sys.nvar), UL = new Float64Array(sys.nvar), UR = new Float64Array(sys.nvar);
  return (pL, pR, ax, out) => {
    const s = Math.max(sys.speed(pL, ax), sys.speed(pR, ax));
    sys.flux(pL, ax, FL); sys.flux(pR, ax, FR); sys.primToCons(pL, UL); sys.primToCons(pR, UR);
    for (let v = 0; v < sys.nvar; v++) out[v] = 0.5 * (FL[v] + FR[v]) - 0.5 * s * (UR[v] - UL[v]);
  };
}

// ------------------------------------------------------------------ systems

/**
 * Linear advection q_t + a·∇q = 0 (a: number in 1D, [ax, ay] in 2D). Exact solution: the initial
 * profile translated by a·t. Riemann solvers: 'exact' (upwind), 'hll', 'rusanov' (all identical here).
 */
export function linearAdvection(a) {
  const vel = typeof a === 'number' ? [a] : Array.from(a ?? []);
  if (!(vel.length === 1 || vel.length === 2) || !vel.every(Number.isFinite)) fail('linearAdvection', 'velocity must be a number or [ax, ay]');
  const sys = {
    name: 'advection', nvar: 1, dims: vel.length, names: ['q'], normal: () => -1,
    consToPrim: (c, p) => { p[0] = c[0]; }, primToCons: (p, c) => { c[0] = p[0]; },
    flux: (p, ax, f) => { f[0] = vel[ax] * p[0]; },
    speed: (_p, ax) => Math.abs(vel[ax]),
    waveSpeeds: (_pL, _pR, ax, S) => { S[0] = Math.min(vel[ax], 0); S[1] = Math.max(vel[ax], 0); },
    admissible: p => Number.isFinite(p[0]),
  };
  sys.solvers = { exact: (pL, pR, ax, out) => { out[0] = vel[ax] * (vel[ax] >= 0 ? pL[0] : pR[0]); } };
  sys.solvers.hll = makeHLL(sys); sys.solvers.rusanov = makeRusanov(sys);
  return sys;
}

/**
 * Inviscid Burgers q_t + (q²/2)_x = 0 (1D; in 2D the flux q²/2 along both axes). 'exact' is the
 * Godunov flux of the entropy solution (shock speed (qL + qR)/2, transonic rarefaction → flux 0).
 */
export function inviscidBurgers({ dims = 1 } = {}) {
  if (dims !== 1 && dims !== 2) fail('inviscidBurgers', 'dims must be 1 or 2');
  const sys = {
    name: 'burgers', nvar: 1, dims, names: ['u'], normal: () => -1,
    consToPrim: (c, p) => { p[0] = c[0]; }, primToCons: (p, c) => { c[0] = p[0]; },
    flux: (p, _ax, f) => { f[0] = 0.5 * p[0] * p[0]; },
    speed: p => Math.abs(p[0]),
    waveSpeeds: (pL, pR, _ax, S) => { S[0] = Math.min(pL[0], pR[0], 0); S[1] = Math.max(pL[0], pR[0], 0); },
    admissible: p => Number.isFinite(p[0]),
  };
  sys.solvers = {
    exact: (pL, pR, _ax, out) => {
      const a = pL[0], b = pR[0];
      if (a > b) out[0] = (a + b) / 2 > 0 ? 0.5 * a * a : 0.5 * b * b;
      else out[0] = a > 0 ? 0.5 * a * a : b < 0 ? 0.5 * b * b : 0;
    },
  };
  sys.solvers.hll = makeHLL(sys); sys.solvers.rusanov = makeRusanov(sys);
  return sys;
}

/**
 * Compressible Euler equations, ideal gas p = (γ − 1)(E − ½ρ|u|²). Conserved [ρ, ρu, (ρv), E],
 * primitive [ρ, u, (v), p]. Riemann solvers: 'hllc' (Toro, Spruce & Speares 1994; default),
 * 'hll' (Einfeldt speeds: positivity-preserving), 'roe' (Harten–Hyman entropy fix; can fail
 * in near-vacuum rarefactions such as Toro's 123 problem), 'rusanov', 'exact' (Godunov flux from
 * eulerExactRiemann; slower). Admissible states: ρ > 0, p > 0.
 */
export function euler({ gamma = 1.4, dims = 1 } = {}) {
  if (!(gamma > 1 && Number.isFinite(gamma))) fail('euler', 'gamma must exceed 1');
  if (dims !== 1 && dims !== 2) fail('euler', 'dims must be 1 or 2');
  const nvar = dims + 2, P = dims + 1, g1 = gamma - 1;
  const ke = p => { let k = 0; for (let d = 0; d < dims; d++) k += p[1 + d] * p[1 + d]; return 0.5 * p[0] * k; };
  const sound = p => Math.sqrt(gamma * p[P] / p[0]);
  const sys = {
    name: 'euler', nvar, dims, gamma, names: dims === 1 ? ['rho', 'mx', 'E'] : ['rho', 'mx', 'my', 'E'],
    primNames: dims === 1 ? ['rho', 'u', 'p'] : ['rho', 'u', 'v', 'p'],
    normal: ax => 1 + ax,
    consToPrim: (c, p) => {
      const r = c[0]; p[0] = r; let k = 0;
      for (let d = 0; d < dims; d++) { p[1 + d] = c[1 + d] / r; k += c[1 + d] * p[1 + d]; }
      p[P] = g1 * (c[P] - 0.5 * k);
    },
    primToCons: (p, c) => { c[0] = p[0]; for (let d = 0; d < dims; d++) c[1 + d] = p[0] * p[1 + d]; c[P] = p[P] / g1 + ke(p); },
    flux: (p, ax, f) => {
      const un = p[1 + ax], E = p[P] / g1 + ke(p);
      f[0] = p[0] * un;
      for (let d = 0; d < dims; d++) f[1 + d] = p[0] * un * p[1 + d];
      f[1 + ax] += p[P]; f[P] = (E + p[P]) * un;
    },
    speed: (p, ax) => Math.abs(p[1 + ax]) + sound(p),
    /** Einfeldt (HLLE) estimates with Roe averages. */
    waveSpeeds: (pL, pR, ax, S) => {
      const rl = Math.sqrt(pL[0]), rr = Math.sqrt(pR[0]), w = 1 / (rl + rr);
      let q2 = 0;
      for (let d = 0; d < dims; d++) { const ud = (rl * pL[1 + d] + rr * pR[1 + d]) * w; q2 += ud * ud; }
      const HL = (pL[P] / g1 + ke(pL) + pL[P]) / pL[0], HR = (pR[P] / g1 + ke(pR) + pR[P]) / pR[0];
      const H = (rl * HL + rr * HR) * w, un = (rl * pL[1 + ax] + rr * pR[1 + ax]) * w, a = Math.sqrt(Math.max(g1 * (H - 0.5 * q2), 0));
      S[0] = Math.min(pL[1 + ax] - sound(pL), un - a); S[1] = Math.max(pR[1 + ax] + sound(pR), un + a);
    },
    admissible: p => p[0] > 0 && p[P] > 0 && Number.isFinite(p[0] + p[P] + p[1] + (dims === 2 ? p[2] : 0)),
  };
  const FL = new Float64Array(nvar), FR = new Float64Array(nvar), S = new Float64Array(2), star = new Float64Array(nvar);
  const hllc = (pL, pR, ax, out) => {
    sys.waveSpeeds(pL, pR, ax, S);
    const sL = S[0], sR = S[1];
    if (sL >= 0) { sys.flux(pL, ax, out); return; }
    if (sR <= 0) { sys.flux(pR, ax, out); return; }
    const rL = pL[0], rR = pR[0], uL = pL[1 + ax], uR = pR[1 + ax];
    const sS = (pR[P] - pL[P] + rL * uL * (sL - uL) - rR * uR * (sR - uR)) / (rL * (sL - uL) - rR * (sR - uR));
    const side = sS >= 0 ? pL : pR, sK = sS >= 0 ? sL : sR, F = sS >= 0 ? FL : FR;
    sys.flux(side, ax, F);
    const r = side[0], u = side[1 + ax], fac = r * (sK - u) / (sK - sS), E = side[P] / g1 + ke(side);
    star[0] = fac;
    for (let d = 0; d < dims; d++) star[1 + d] = fac * (d === ax ? sS : side[1 + d]);
    star[P] = fac * (E / r + (sS - u) * (sS + side[P] / (r * (sK - u))));
    // F* = F_K + s_K (U*_K − U_K)
    out[0] = F[0] + sK * (star[0] - r);
    for (let d = 0; d < dims; d++) out[1 + d] = F[1 + d] + sK * (star[1 + d] - r * side[1 + d]);
    out[P] = F[P] + sK * (star[P] - E);
  };
  const roe = (pL, pR, ax, out) => {
    sys.flux(pL, ax, FL); sys.flux(pR, ax, FR);
    const rl = Math.sqrt(pL[0]), rr = Math.sqrt(pR[0]), w = 1 / (rl + rr);
    const HL = (pL[P] / g1 + ke(pL) + pL[P]) / pL[0], HR = (pR[P] / g1 + ke(pR) + pR[P]) / pR[0], H = (rl * HL + rr * HR) * w;
    const u = (rl * pL[1] + rr * pR[1]) * w, v = dims === 2 ? (rl * pL[2] + rr * pR[2]) * w : 0;
    const q2 = u * u + v * v, a2 = g1 * (H - 0.5 * q2);
    if (!(a2 > 0)) fail('euler roe', 'non-physical Roe average (near vacuum); use riemann "hllc" or "hll"');
    const a = Math.sqrt(a2), rho = rl * rr, un = ax === 0 ? u : v, ut = ax === 0 ? v : u;
    const dr = pR[0] - pL[0], dp = pR[P] - pL[P], dun = pR[1 + ax] - pL[1 + ax], dut = dims === 2 ? pR[2 - ax] - pL[2 - ax] : 0;
    const a1 = (dp - rho * a * dun) / (2 * a2), a2w = dr - dp / a2, a3 = (dp + rho * a * dun) / (2 * a2), at = rho * dut;
    // Harten–Hyman entropy fix on the acoustic waves
    const fix = (lam, lamL, lamR) => { const d = Math.max(0, lam - lamL, lamR - lam); return Math.abs(lam) < d ? (lam * lam + d * d) / (2 * d) : Math.abs(lam); };
    const l1 = fix(un - a, pL[1 + ax] - sound(pL), pR[1 + ax] - sound(pR)), l2 = Math.abs(un), l3 = fix(un + a, pL[1 + ax] + sound(pL), pR[1 + ax] + sound(pR));
    // eigenvectors in (ρ, m_normal, m_tangential, E) ordering, mapped to axes below
    const k1 = l1 * a1, k2 = l2 * a2w, k3 = l3 * a3, kt = l2 * at;
    const dMass = k1 + k2 + k3;
    const dMn = k1 * (un - a) + k2 * un + k3 * (un + a);
    const dMt = (k1 + k2 + k3) * ut + kt;
    const dE = k1 * (H - un * a) + k2 * 0.5 * q2 + k3 * (H + un * a) + kt * ut;
    out[0] = 0.5 * (FL[0] + FR[0] - dMass);
    out[1 + ax] = 0.5 * (FL[1 + ax] + FR[1 + ax] - dMn);
    if (dims === 2) out[2 - ax] = 0.5 * (FL[2 - ax] + FR[2 - ax] - dMt);
    out[P] = 0.5 * (FL[P] + FR[P] - dE);
  };
  const sample = new Float64Array(3), prim = new Float64Array(nvar);
  const exact = (pL, pR, ax, out) => {
    const [ps, us] = eulerStar(pL[0], pL[1 + ax], pL[P], pR[0], pR[1 + ax], pR[P], gamma);
    eulerSample(pL[0], pL[1 + ax], pL[P], pR[0], pR[1 + ax], pR[P], gamma, ps, us, 0, sample);
    prim[0] = sample[0]; prim[1 + ax] = sample[1]; prim[P] = sample[2];
    if (dims === 2) prim[2 - ax] = (us >= 0 ? pL : pR)[2 - ax];
    sys.flux(prim, ax, out);
  };
  sys.solvers = { hllc, roe, exact, hll: makeHLL(sys), rusanov: makeRusanov(sys) };
  return sys;
}

/**
 * Shallow-water equations h_t + (hu)_x = 0, (hu)_t + (hu² + ½gh²)_x = −g h b_x (2D adds hv).
 * Conserved [h, hu, (hv)], primitive [h, u, (v)]. Cells with h ≤ dry (default 1e−10 of the length
 * unit) have u = 0. bathymetry: bed elevation b per cell (same size as the grid) → well-balanced
 * hydrostatic reconstruction (Audusse et al., SIAM J. Sci. Comput. 25, 2050 (2004)): a lake at
 * rest η = h + b = const stays at rest to round-off, wet/dry fronts keep h ≥ 0 (CFL ≤ 0.5).
 * Riemann solvers: 'hllc' (default; HLL + contact for hv), 'hll', 'rusanov', 'exact'
 * (shallowWaterExactRiemann; slower). Linear wave speed √(gh).
 */
export function shallowWater({ g = 9.81, dims = 1, bathymetry, dry = 1e-10 } = {}) {
  if (!(g > 0 && Number.isFinite(g))) fail('shallowWater', 'g must be positive');
  if (dims !== 1 && dims !== 2) fail('shallowWater', 'dims must be 1 or 2');
  if (!(dry >= 0)) fail('shallowWater', 'dry threshold must be ≥ 0');
  if (bathymetry && !Array.from(bathymetry).every(Number.isFinite)) fail('shallowWater', 'bathymetry must be finite');
  const nvar = dims + 1;
  const sys = {
    name: 'shallow-water', nvar, dims, g, dry, bathymetry: bathymetry ? Float64Array.from(bathymetry) : null,
    names: dims === 1 ? ['h', 'hu'] : ['h', 'hu', 'hv'], primNames: dims === 1 ? ['h', 'u'] : ['h', 'u', 'v'],
    normal: ax => 1 + ax,
    consToPrim: (c, p) => { const h = c[0]; p[0] = h; for (let d = 0; d < dims; d++) p[1 + d] = h > dry ? c[1 + d] / h : 0; },
    primToCons: (p, c) => { c[0] = p[0]; for (let d = 0; d < dims; d++) c[1 + d] = p[0] * p[1 + d]; },
    flux: (p, ax, f) => {
      const h = p[0], un = p[1 + ax];
      f[0] = h * un;
      for (let d = 0; d < dims; d++) f[1 + d] = h * un * p[1 + d];
      f[1 + ax] += 0.5 * g * h * h;
    },
    speed: (p, ax) => Math.abs(p[1 + ax]) + Math.sqrt(g * Math.max(p[0], 0)),
    /** Toro (2001) estimates: two-rarefaction depth h*, shock-corrected; dry-bed fronts u ± 2√(gh). */
    waveSpeeds: (pL, pR, ax, S) => {
      const hL = pL[0], hR = pR[0], uL = pL[1 + ax], uR = pR[1 + ax], aL = Math.sqrt(g * Math.max(hL, 0)), aR = Math.sqrt(g * Math.max(hR, 0));
      if (hL <= dry) { S[0] = uR - 2 * aR; S[1] = uR + aR; return; }
      if (hR <= dry) { S[0] = uL - aL; S[1] = uL + 2 * aL; return; }
      const as = Math.max(0, 0.5 * (aL + aR) + 0.25 * (uL - uR)), hs = as * as / g;
      const qL = hs > hL ? Math.sqrt(0.5 * (hs + hL) * hs) / hL : 1, qR = hs > hR ? Math.sqrt(0.5 * (hs + hR) * hs) / hR : 1;
      S[0] = uL - aL * qL; S[1] = uR + aR * qR;
    },
    admissible: p => p[0] >= 0 && Number.isFinite(p[0] + p[1] + (dims === 2 ? p[2] : 0)),
  };
  const hll = makeHLL(sys), S = new Float64Array(2);
  const hllc = (pL, pR, ax, out) => {
    hll(pL, pR, ax, out);
    if (dims === 1) return;
    sys.waveSpeeds(pL, pR, ax, S);
    const sL = S[0], sR = S[1], hL = pL[0], hR = pR[0], uL = pL[1 + ax], uR = pR[1 + ax];
    const den = hR * (uR - sR) - hL * (uL - sL);
    const sS = den !== 0 ? (sL * hR * (uR - sR) - sR * hL * (uL - sL)) / den : 0.5 * (uL + uR);
    out[2 - ax] = out[0] * (sS >= 0 ? pL[2 - ax] : pR[2 - ax]);           // tangential velocity carried by the contact
  };
  const sample = new Float64Array(2), prim = new Float64Array(nvar);
  const exact = (pL, pR, ax, out) => {
    swSample(pL[0], pL[1 + ax], pR[0], pR[1 + ax], g, 0, sample, dry);
    prim[0] = sample[0]; prim[1 + ax] = sample[1];
    if (dims === 2) prim[2 - ax] = (sample[1] >= 0 ? pL : pR)[2 - ax];
    sys.flux(prim, ax, out);
  };
  sys.solvers = { hll, hllc, exact, rusanov: makeRusanov(sys) };
  return sys;
}

// ------------------------------------------------------------------ finite-volume engine

const SIDES = ['periodic', 'transmissive', 'reflective'];
function parseBc(bc, dims, who) {
  const side = s => {
    if (SIDES.includes(s)) return { type: s };
    if (s && typeof s === 'object' && s.state) return { type: 'state', state: Float64Array.from(s.state) };
    return fail(who, `boundary must be ${SIDES.join(', ')} or { state: primitive values }`);
  };
  let axes;
  if (typeof bc === 'string') axes = Array.from({ length: dims }, () => [bc, bc]);
  else if (dims === 1 && Array.isArray(bc) && bc.length === 2 && !Array.isArray(bc[0])) axes = [bc];
  else if (Array.isArray(bc) && bc.length === dims) axes = bc.map(e => (Array.isArray(e) ? e : [e, e]));
  else fail(who, 'bc must be a string, [left, right] (1D), or one entry per axis');
  const out = axes.map(([l, r]) => [side(l), side(r)]);
  for (const [l, r] of out) if ((l.type === 'periodic') !== (r.type === 'periodic')) fail(who, 'periodic must be set on both sides of an axis');
  return out;
}

/**
 * Finite-volume solver for q_t + ∇·F(q) = S on a uniform 1D or 2D grid of cell averages.
 *   system      linearAdvection | inviscidBurgers | euler | shallowWater (system.dims = shape.length)
 *   q0          array of conserved-variable arrays (system.names order), each of grid size
 *   shape       [n] or [nx, ny] (n ≥ 2 per axis); spacing: number or [dx, dy]
 *   bc          'periodic' | 'transmissive' (zero gradient) | 'reflective' (wall: normal velocity
 *               mirrored) | { state: primitive values } (fixed inflow); a string for all sides,
 *               [left, right] in 1D, or one entry (string or [low, high]) per axis
 *   limiter     'mc' (default) | 'vanleer' | 'minmod' | 'superbee' | 'none' (unlimited, oscillates)
 *               | 'first' (Godunov); applied to primitive variables
 *   riemann     a key of system.solvers ('hllc' for euler/shallowWater, 'exact' for scalars)
 *   integrator  'ssprk3' (default, Shu–Osher) | 'ssprk2' (Heun)
 *   cfl         Courant number in (0, 1] for the automatic step: dt = cfl / max Σ_a (|λ|_a/dx_a)
 *               (≤ 0.5 recommended for MUSCL; positivity of h with bathymetry needs ≤ 0.5)
 * Faces whose reconstructed state is not admissible (ρ, p ≤ 0, h < 0) fall back to first order.
 * step(dt?) takes the CFL step, or the given dt (throws if its Courant number exceeds 1). After every
 * stage, a non-admissible or non-finite cell throws (reduce cfl, use 'minmod' or 'hll').
 * Conservation: totals()[v] − totals0[v] = inflow[v] to round-off for every conserved variable
 * whose equation has no source (all of them except momentum with bathymetry); inflow accumulates
 * the boundary fluxes with the integrator's weights. Cost per step: stages × faces Riemann solves.
 * Returns { t, q, steps, dt, step(dt?), advanceTo(t, maxSteps?), maxDt(), totals(), inflow,
 *           primitive(v) (copy of primitive variable v), courant }.
 */
export function finiteVolume({ system, q0, shape, spacing = 1, bc = 'transmissive', limiter = 'mc', riemann, integrator = 'ssprk3', cfl = 0.4, t0 = 0 }) {
  const who = 'finiteVolume';
  if (!system?.solvers) fail(who, 'system must come from linearAdvection, inviscidBurgers, euler or shallowWater');
  const dims = shape?.length;
  if (dims !== 1 && dims !== 2) fail(who, 'shape must be [n] or [nx, ny]');
  if (system.dims !== dims) fail(who, `system is ${system.dims}D but shape is ${dims}D`);
  if (!shape.every(n => Number.isInteger(n) && n >= 2)) fail(who, 'each axis needs at least 2 cells');
  const dx = typeof spacing === 'number' ? shape.map(() => spacing) : Array.from(spacing);
  if (dx.length !== dims || !dx.every(h => h > 0 && Number.isFinite(h))) fail(who, 'spacing must be positive (one per axis)');
  const lim = LIMITERS[limiter];
  if (!lim) fail(who, `limiter must be one of ${Object.keys(LIMITERS).join(', ')}`);
  const solverName = riemann ?? (system.solvers.hllc ? 'hllc' : 'exact');
  const solve = system.solvers[solverName];
  if (!solve) fail(who, `riemann must be one of ${Object.keys(system.solvers).join(', ')}`);
  if (!['ssprk2', 'ssprk3'].includes(integrator)) fail(who, "integrator must be 'ssprk2' or 'ssprk3'");
  if (!(cfl > 0 && cfl <= 1)) fail(who, 'cfl must be in (0, 1]');
  const bcs = parseBc(bc, dims, who);
  const nv = system.nvar, size = shape.reduce((a, b) => a * b, 1), strides = [1, shape[0]];
  if (!Array.isArray(q0) || q0.length !== nv) fail(who, `q0 must hold ${nv} arrays (${system.names.join(', ')})`);
  for (const a of q0) { if (a.length !== size) fail(who, 'each q0 array must have the grid size'); for (let i = 0; i < size; i++) if (!Number.isFinite(a[i])) fail(who, 'q0 must be finite'); }
  const bed = system.bathymetry;
  if (bed && bed.length !== size) fail(who, 'bathymetry must have the grid size');
  for (const [l, r] of bcs) for (const s of [l, r]) if (s.type === 'state' && s.state.length !== nv) fail(who, `boundary state needs ${nv} primitive values`);
  const hydro = !!bed, g = system.g;

  const q = q0.map(a => Float64Array.from(a));
  const k = q.map(() => new Float64Array(size)), base = q.map(() => new Float64Array(size)), stage = q.map(() => new Float64Array(size));
  const nmax = Math.max(...shape), pad = nmax + 4;
  const P = Array.from({ length: nv }, () => new Float64Array(pad)), slope = Array.from({ length: nv }, () => new Float64Array(pad));
  const B = new Float64Array(pad), Eta = new Float64Array(pad), slopeEta = new Float64Array(pad);
  const FLf = Array.from({ length: nv }, () => new Float64Array(nmax + 1)), FRf = Array.from({ length: nv }, () => new Float64Array(nmax + 1));
  const hFaceL = new Float64Array(nmax + 1), hFaceR = new Float64Array(nmax + 1), bFaceL = new Float64Array(nmax + 1), bFaceR = new Float64Array(nmax + 1);
  const c = new Float64Array(nv), p = new Float64Array(nv), pL = new Float64Array(nv), pR = new Float64Array(nv), F = new Float64Array(nv);
  const boundaryRate = new Float64Array(nv);
  const inflow = new Float64Array(nv);
  const cellVolume = dx.reduce((a, b) => a * b, 1);

  // fill the padded primitive line for axis ax, line start `start`, then ghosts
  const gather = (Q, ax, start, n) => {
    const st = strides[ax];
    for (let i = 0; i < n; i++) {
      const idx = start + i * st;
      for (let v = 0; v < nv; v++) c[v] = Q[v][idx];
      system.consToPrim(c, p);
      for (let v = 0; v < nv; v++) P[v][i + 2] = p[v];
      if (hydro) B[i + 2] = bed[idx];
    }
    const [lo, hi] = bcs[ax], nrm = system.normal(ax);
    const ghost = (side, gIdx, src) => {
      if (side.type === 'state') { for (let v = 0; v < nv; v++) P[v][gIdx] = side.state[v]; if (hydro) B[gIdx] = B[src]; return; }
      for (let v = 0; v < nv; v++) P[v][gIdx] = P[v][src];
      if (side.type === 'reflective' && nrm >= 0) P[nrm][gIdx] = -P[nrm][src];
      if (hydro) B[gIdx] = B[src];
    };
    if (lo.type === 'periodic') {
      for (let v = 0; v < nv; v++) { P[v][0] = P[v][n]; P[v][1] = P[v][n + 1]; P[v][n + 2] = P[v][2]; P[v][n + 3] = P[v][3]; }
      if (hydro) { B[0] = B[n]; B[1] = B[n + 1]; B[n + 2] = B[2]; B[n + 3] = B[3]; }
    } else {
      const refl = s => s.type === 'reflective';
      ghost(lo, 1, 2); ghost(lo, 0, refl(lo) ? 3 : 2);
      ghost(hi, n + 2, n + 1); ghost(hi, n + 3, refl(hi) ? n : n + 1);
    }
  };

  const rhs = (Q, out) => {
    for (let v = 0; v < nv; v++) out[v].fill(0);
    boundaryRate.fill(0);
    for (let ax = 0; ax < dims; ax++) {
      const n = shape[ax], other = dims === 2 ? shape[1 - ax] : 1, st = strides[ax], ostride = dims === 2 ? strides[1 - ax] : 0;
      const area = cellVolume / dx[ax], invdx = 1 / dx[ax];
      for (let line = 0; line < other; line++) {
        const start = line * ostride;
        gather(Q, ax, start, n);
        for (let v = 0; v < nv; v++) { const a = P[v], s = slope[v]; for (let i = 1; i <= n + 2; i++) s[i] = lim(a[i] - a[i - 1], a[i + 1] - a[i]); }
        if (hydro) {
          for (let i = 0; i <= n + 3; i++) Eta[i] = P[0][i] + B[i];
          for (let i = 1; i <= n + 2; i++) slopeEta[i] = lim(Eta[i] - Eta[i - 1], Eta[i + 1] - Eta[i]);
        }
        for (let f = 0; f <= n; f++) {
          const il = f + 1, ir = f + 2;
          for (let v = 0; v < nv; v++) { pL[v] = P[v][il] + 0.5 * slope[v][il]; pR[v] = P[v][ir] - 0.5 * slope[v][ir]; }
          const firstL = !system.admissible(pL), firstR = !system.admissible(pR);
          if (firstL) for (let v = 0; v < nv; v++) pL[v] = P[v][il];
          if (firstR) for (let v = 0; v < nv; v++) pR[v] = P[v][ir];
          if (!hydro) {
            solve(pL, pR, ax, F);
            for (let v = 0; v < nv; v++) FLf[v][f] = FRf[v][f] = F[v];
            continue;
          }
          // hydrostatic reconstruction: faces of η and h give faces of b (lake at rest preserved)
          const etaL = firstL ? Eta[il] : Eta[il] + 0.5 * slopeEta[il], etaR = firstR ? Eta[ir] : Eta[ir] - 0.5 * slopeEta[ir];
          const hL = pL[0], hR = pR[0], bL = etaL - hL, bR = etaR - hR, bs = Math.max(bL, bR);
          hFaceL[f] = hL; hFaceR[f] = hR; bFaceL[f] = bL; bFaceR[f] = bR;
          const hsL = Math.max(0, hL + bL - bs), hsR = Math.max(0, hR + bR - bs);
          pL[0] = hsL; pR[0] = hsR;
          solve(pL, pR, ax, F);
          for (let v = 0; v < nv; v++) { FLf[v][f] = F[v]; FRf[v][f] = F[v]; }
          FLf[1 + ax][f] += 0.5 * g * (hL * hL - hsL * hsL);
          FRf[1 + ax][f] += 0.5 * g * (hR * hR - hsR * hsR);
        }
        for (let i = 0; i < n; i++) {
          const idx = start + i * st;
          for (let v = 0; v < nv; v++) out[v][idx] -= (FLf[v][i + 1] - FRf[v][i]) * invdx;
          if (hydro) out[1 + ax][idx] -= g * 0.5 * (hFaceR[i] + hFaceL[i + 1]) * (bFaceL[i + 1] - bFaceR[i]) * invdx;
        }
        if (bcs[ax][0].type !== 'periodic') for (let v = 0; v < nv; v++) boundaryRate[v] += (FRf[v][0] - FLf[v][n]) * area;
      }
    }
  };

  const maxRate = Q => {
    let m = 0;
    for (let i = 0; i < size; i++) {
      for (let v = 0; v < nv; v++) c[v] = Q[v][i];
      system.consToPrim(c, p);
      let r = 0;
      for (let ax = 0; ax < dims; ax++) r += system.speed(p, ax) / dx[ax];
      if (r > m) m = r;
    }
    return m;
  };
  const verify = (Q, when) => {
    for (let i = 0; i < size; i++) {
      for (let v = 0; v < nv; v++) c[v] = Q[v][i];
      if (hydro || system.name === 'shallow-water') {
        if (c[0] < 0 && c[0] > -1e-13) { Q[0][i] = 0; c[0] = 0; }
      }
      system.consToPrim(c, p);
      if (!system.admissible(p)) fail(who, `non-admissible state (positivity or finiteness lost) in cell ${i} at t = ${when}; reduce cfl, use limiter 'minmod' or riemann 'hll'`);
    }
  };
  verify(q, t0);

  const weights = integrator === 'ssprk3' ? [1 / 6, 1 / 6, 2 / 3] : [0.5, 0.5];
  const s = {
    t: t0, q, steps: 0, dt: 0, courant: 0, inflow, riemann: solverName,
    maxDt() { const r = maxRate(q); return r > 0 ? cfl / r : Infinity; },
    step(dtIn) {
      const r = maxRate(q);
      const dt = dtIn ?? (r > 0 ? cfl / r : fail(who, 'all wave speeds are zero: pass dt explicitly'));
      if (!(dt > 0 && Number.isFinite(dt))) fail(who, 'dt must be positive and finite');
      s.courant = dt * r;
      if (s.courant > 1) fail(who, `Courant number dt·Σ|λ|/dx = ${s.courant} exceeds 1 (CFL violated)`);
      for (let v = 0; v < nv; v++) base[v].set(q[v]);
      const accumulate = w => { for (let v = 0; v < nv; v++) inflow[v] += w * dt * boundaryRate[v]; };
      // stage 1
      rhs(q, k); accumulate(weights[0]);
      for (let v = 0; v < nv; v++) for (let i = 0; i < size; i++) stage[v][i] = base[v][i] + dt * k[v][i];
      verify(stage, s.t + dt);
      rhs(stage, k); accumulate(weights[1]);
      if (integrator === 'ssprk2') {
        for (let v = 0; v < nv; v++) for (let i = 0; i < size; i++) q[v][i] = 0.5 * base[v][i] + 0.5 * (stage[v][i] + dt * k[v][i]);
      } else {
        for (let v = 0; v < nv; v++) for (let i = 0; i < size; i++) stage[v][i] = 0.75 * base[v][i] + 0.25 * (stage[v][i] + dt * k[v][i]);
        verify(stage, s.t + dt / 2);
        rhs(stage, k); accumulate(weights[2]);
        for (let v = 0; v < nv; v++) for (let i = 0; i < size; i++) q[v][i] = base[v][i] / 3 + (2 / 3) * (stage[v][i] + dt * k[v][i]);
      }
      verify(q, s.t + dt);
      s.t += dt; s.dt = dt; s.steps++;
      return dt;
    },
    advanceTo(target, maxSteps = 1e7) {
      if (!Number.isFinite(target)) fail(who, 'target time must be finite');
      let count = 0;
      while (s.t < target) {
        if (count++ >= maxSteps) return false;
        const r = maxRate(q), dt = r > 0 ? Math.min(cfl / r, target - s.t) : target - s.t;
        if (!(dt > 0)) break;
        s.step(dt);
      }
      return true;
    },
    /** Σ q_v dV per conserved variable (Kahan-free plain sum; round-off ~ size·ε). */
    totals() { return q.map(a => { let t = 0; for (let i = 0; i < size; i++) t += a[i]; return t * cellVolume; }); },
    primitive(v) {
      const out = new Float64Array(size);
      for (let i = 0; i < size; i++) { for (let w = 0; w < nv; w++) c[w] = q[w][i]; system.consToPrim(c, p); out[i] = p[v]; }
      return out;
    },
  };
  return s;
}

// ------------------------------------------------------------------ exact Riemann solvers

/** Entropy solution of the Burgers Riemann problem: sample(ξ = x/t). Shock speed (uL + uR)/2 if uL > uR, else a centred rarefaction. */
export function burgersExactRiemann(uL, uR) {
  if (!(Number.isFinite(uL) && Number.isFinite(uR))) fail('burgersExactRiemann', 'states must be finite');
  const shock = uL > uR, speed = (uL + uR) / 2;
  return { shock, speed, sample: xi => (shock ? (xi < speed ? uL : uR) : xi <= uL ? uL : xi >= uR ? uR : xi) };
}

const STAR = new Float64Array(2);
function eulerStar(rL, uL, pL, rR, uR, pR, g) {
  const aL = Math.sqrt(g * pL / rL), aR = Math.sqrt(g * pR / rR), du = uR - uL;
  if (2 * (aL + aR) / (g - 1) <= du) fail('eulerExactRiemann', 'the states generate a vacuum (2(aL + aR)/(γ − 1) ≤ uR − uL)');
  const z = (g - 1) / (2 * g), AL = 2 / ((g + 1) * rL), AR = 2 / ((g + 1) * rR), BL = (g - 1) / (g + 1) * pL, BR = (g - 1) / (g + 1) * pR;
  const fK = (pp, pK, aK, rK, A, Bk, out, o) => {
    if (pp > pK) { const s = Math.sqrt(A / (pp + Bk)); out[o] = (pp - pK) * s; out[o + 1] = s * (1 - (pp - pK) / (2 * (Bk + pp))); } else { const r = (pp / pK) ** z; out[o] = 2 * aK / (g - 1) * (r - 1); out[o + 1] = 1 / (rK * aK) * (pp / pK) ** (-(g + 1) / (2 * g)); }
  };
  // initial guess: Toro's adaptive PVRS / two-rarefaction / two-shock choice
  const pmin = Math.min(pL, pR), pmax = Math.max(pL, pR);
  const ppv = Math.max(0, 0.5 * (pL + pR) - 0.125 * du * (rL + rR) * (aL + aR));
  let pp;
  if (pmax / pmin <= 2 && ppv >= pmin && ppv <= pmax) pp = ppv;
  else if (ppv < pmin) pp = ((aL + aR - 0.5 * (g - 1) * du) / (aL / pL ** z + aR / pR ** z)) ** (1 / z);
  else { const gL = Math.sqrt(AL / (ppv + BL)), gR = Math.sqrt(AR / (ppv + BR)); pp = (gL * pL + gR * pR - du) / (gL + gR); }
  pp = Math.max(pp, 1e-14 * pmax);
  const f = new Float64Array(4);
  for (let it = 0; it < 100; it++) {
    fK(pp, pL, aL, rL, AL, BL, f, 0); fK(pp, pR, aR, rR, AR, BR, f, 2);
    let next = pp - (f[0] + f[2] + du) / (f[1] + f[3]);
    if (next <= 0) next = 0.1 * pp;
    const change = 2 * Math.abs(next - pp) / (next + pp);
    pp = next;
    if (change < 1e-14) break;
    if (it === 99) fail('eulerExactRiemann', 'pressure iteration did not converge');
  }
  fK(pp, pL, aL, rL, AL, BL, f, 0); fK(pp, pR, aR, rR, AR, BR, f, 2);
  STAR[0] = pp; STAR[1] = 0.5 * (uL + uR) + 0.5 * (f[2] - f[0]);
  return STAR;
}
function eulerSample(rL, uL, pL, rR, uR, pR, g, ps, us, xi, out) {
  const g6 = (g - 1) / (g + 1), z = (g - 1) / (2 * g), aL = Math.sqrt(g * pL / rL), aR = Math.sqrt(g * pR / rR);
  if (xi <= us) {
    if (ps > pL) {
      const sL = uL - aL * Math.sqrt((g + 1) / (2 * g) * ps / pL + z);
      if (xi <= sL) { out[0] = rL; out[1] = uL; out[2] = pL; } else { out[0] = rL * (ps / pL + g6) / (g6 * ps / pL + 1); out[1] = us; out[2] = ps; }
    } else {
      const head = uL - aL, as = aL * (ps / pL) ** z, tail = us - as;
      if (xi <= head) { out[0] = rL; out[1] = uL; out[2] = pL; } else if (xi > tail) { out[0] = rL * (ps / pL) ** (1 / g); out[1] = us; out[2] = ps; } else {
        const c = 2 / (g + 1) + g6 / aL * (uL - xi);
        out[0] = rL * c ** (2 / (g - 1)); out[1] = 2 / (g + 1) * (aL + (g - 1) / 2 * uL + xi); out[2] = pL * c ** (2 * g / (g - 1));
      }
    }
  } else if (ps > pR) {
    const sR = uR + aR * Math.sqrt((g + 1) / (2 * g) * ps / pR + z);
    if (xi >= sR) { out[0] = rR; out[1] = uR; out[2] = pR; } else { out[0] = rR * (ps / pR + g6) / (g6 * ps / pR + 1); out[1] = us; out[2] = ps; }
  } else {
    const head = uR + aR, as = aR * (ps / pR) ** z, tail = us + as;
    if (xi >= head) { out[0] = rR; out[1] = uR; out[2] = pR; } else if (xi < tail) { out[0] = rR * (ps / pR) ** (1 / g); out[1] = us; out[2] = ps; } else {
      const c = 2 / (g + 1) - g6 / aR * (uR - xi);
      out[0] = rR * c ** (2 / (g - 1)); out[1] = 2 / (g + 1) * (-aR + (g - 1) / 2 * uR + xi); out[2] = pR * c ** (2 * g / (g - 1));
    }
  }
}

/**
 * Exact solution of the Euler Riemann problem (Toro ch. 4): Newton iteration on the pressure
 * function to 1e−14 relative, then wave sampling. left/right: { rho, u, p } (> 0). Throws when the data
 * generate a vacuum. Returns { pStar, uStar, rhoStarL, rhoStarR, sample(ξ = x/t) → { rho, u, p } }.
 */
export function eulerExactRiemann(left, right, gamma = 1.4) {
  const { rho: rL, u: uL, p: pL } = left ?? {}, { rho: rR, u: uR, p: pR } = right ?? {};
  if (!(rL > 0 && rR > 0 && pL > 0 && pR > 0 && Number.isFinite(uL) && Number.isFinite(uR))) fail('eulerExactRiemann', 'need rho > 0, p > 0 and finite u on both sides');
  if (!(gamma > 1)) fail('eulerExactRiemann', 'gamma must exceed 1');
  const [ps, us] = eulerStar(rL, uL, pL, rR, uR, pR, gamma);
  const g6 = (gamma - 1) / (gamma + 1);
  const rhoStar = (rK, pK) => (ps > pK ? rK * (ps / pK + g6) / (g6 * ps / pK + 1) : rK * (ps / pK) ** (1 / gamma));
  const out = new Float64Array(3);
  return {
    pStar: ps, uStar: us, rhoStarL: rhoStar(rL, pL), rhoStarR: rhoStar(rR, pR),
    sample(xi) { eulerSample(rL, uL, pL, rR, uR, pR, gamma, ps, us, xi, out); return { rho: out[0], u: out[1], p: out[2] }; },
  };
}

function swStar(hL, uL, hR, uR, g) {
  const aL = Math.sqrt(g * hL), aR = Math.sqrt(g * hR);
  // f_K(h): bore (h > h_K) or rarefaction branch, with its exact derivative
  const fK = (h, hK, aK, o, f) => {
    if (h > hK) { const sq = Math.sqrt(0.5 * g * (h + hK) / (h * hK)); f[o] = (h - hK) * sq; f[o + 1] = sq - (h - hK) * g / (4 * sq * h * h); } else { f[o] = 2 * (Math.sqrt(g * h) - aK); f[o + 1] = Math.sqrt(g / h); }
  };
  const f = new Float64Array(4);
  let h = Math.max((0.5 * (aL + aR) + 0.25 * (uL - uR)) ** 2 / g, 1e-12 * Math.max(hL, hR));
  for (let it = 0; it < 100; it++) {
    fK(h, hL, aL, 0, f); fK(h, hR, aR, 2, f);
    let next = h - (f[0] + f[2] + uR - uL) / (f[1] + f[3]);
    if (next <= 0) next = 0.1 * h;
    const change = 2 * Math.abs(next - h) / (next + h);
    h = next;
    if (change < 1e-14) break;
    if (it === 99) fail('shallowWaterExactRiemann', 'depth iteration did not converge');
  }
  fK(h, hL, aL, 0, f); fK(h, hR, aR, 2, f);
  STAR[0] = h; STAR[1] = 0.5 * (uL + uR) + 0.5 * (f[2] - f[0]);
  return STAR;
}
function swSample(hL, uL, hR, uR, g, xi, out, dry = 0) {
  const aL = Math.sqrt(g * Math.max(hL, 0)), aR = Math.sqrt(g * Math.max(hR, 0));
  const leftFan = () => { const a = (uL + 2 * aL - xi) / 3; out[0] = a * a / g; out[1] = (uL + 2 * aL + 2 * xi) / 3; };
  const rightFan = () => { const a = (-uR + 2 * aR + xi) / 3; out[0] = a * a / g; out[1] = (uR - 2 * aR + 2 * xi) / 3; };
  const set = (h, u) => { out[0] = h; out[1] = u; };
  if (hL <= dry && hR <= dry) return set(0, 0);
  if (hR <= dry) { if (xi <= uL - aL) return set(hL, uL); if (xi < uL + 2 * aL) return leftFan(); return set(0, 0); }
  if (hL <= dry) { if (xi >= uR + aR) return set(hR, uR); if (xi > uR - 2 * aR) return rightFan(); return set(0, 0); }
  if (2 * (aL + aR) <= uR - uL) {                                                // dry region created in the middle
    if (xi <= uL - aL) return set(hL, uL);
    if (xi < uL + 2 * aL) return leftFan();
    if (xi <= uR - 2 * aR) return set(0, 0);
    if (xi < uR + aR) return rightFan();
    return set(hR, uR);
  }
  const [hs, us] = swStar(hL, uL, hR, uR, g), as = Math.sqrt(g * hs);
  if (xi <= us) {
    if (hs > hL) { const sL = uL - aL * Math.sqrt(0.5 * (hs + hL) * hs) / hL; return xi <= sL ? set(hL, uL) : set(hs, us); }
    if (xi <= uL - aL) return set(hL, uL);
    if (xi >= us - as) return set(hs, us);
    return leftFan();
  }
  if (hs > hR) { const sR = uR + aR * Math.sqrt(0.5 * (hs + hR) * hs) / hR; return xi >= sR ? set(hR, uR) : set(hs, us); }
  if (xi >= uR + aR) return set(hR, uR);
  if (xi <= us + as) return set(hs, us);
  return rightFan();
}

/**
 * Exact shallow-water Riemann solution (Toro 2001 ch. 5) for a flat bed, wet or dry: dam break
 * over a dry bed (Ritter 1892: h = (2√(gh₀) − ξ)²/(9g) for −√(gh₀) ≤ ξ ≤ 2√(gh₀)), wet bed
 * (Stoker 1957: rarefaction + bore), and dry-region formation. left/right: { h ≥ 0, u }.
 * Returns { hStar, uStar (null when a dry region forms), sample(ξ = x/t) → { h, u } }.
 */
export function shallowWaterExactRiemann(left, right, g = 9.81) {
  const { h: hL, u: uL } = left ?? {}, { h: hR, u: uR } = right ?? {};
  if (!(hL >= 0 && hR >= 0 && Number.isFinite(uL) && Number.isFinite(uR) && g > 0)) fail('shallowWaterExactRiemann', 'need h ≥ 0, finite u and g > 0');
  const wet = hL > 0 && hR > 0 && 2 * (Math.sqrt(g * hL) + Math.sqrt(g * hR)) > uR - uL;
  const star = wet ? Array.from(swStar(hL, uL, hR, uR, g)) : [null, null];
  const out = new Float64Array(2);
  return { hStar: star[0], uStar: star[1], sample(xi) { swSample(hL, uL, hR, uR, g, xi, out); return { h: out[0], u: out[1] }; } };
}
