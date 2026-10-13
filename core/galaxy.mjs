// Galaxy dynamics: stars orbiting in the gravitational potential of a disk galaxy
// with a dark-matter halo, and its rotation curve.
//
// Mass model: three standard analytic components with exact potentials and forces,
// normalised like MWPotential2014 (Bovy 2015, ApJS 216, 29): circular speed 220 km/s at
// R₀ = 8 kpc, with radial-force fractions 5 % bulge, 60 % disk, 35 % halo.
//   bulge  Hernquist (1990)          Φ = −GM/(r + a)                  (substitutes Bovy's cut-off power law)
//   disk   Miyamoto & Nagai (1975)   Φ = −GM/√(R² + (a + √(z² + b²))²), a = 3 kpc, b = 0.28 kpc
//   halo   Navarro–Frenk–White (1996) Φ = −4πGρ₀r_s³ ln(1 + r/r_s)/r, r_s = 16 kpc
// Units: kpc, km/s, M☉; time unit kpc/(km/s) = 977.8 Myr.

export const G = 4.30091727e-6;                    // kpc (km/s)² / M☉ (IAU/CODATA)
export const TIME_UNIT_MYR = 977.792221;           // 1 kpc / (1 km/s) in Myr
export const R0 = 8, V0 = 220;                     // normalisation point (kpc, km/s)

const HERNQUIST_A = 0.6, MN_A = 3, MN_B = 0.28, NFW_RS = 16;

// Unnormalised component shapes: radial force at (R, z=0) per unit "mass parameter".
const bulgeForceR = (R, a = HERNQUIST_A) => G / ((R + a) ** 2);
const diskForceR = (R, a = MN_A, b = MN_B) => G * R / Math.pow(R * R + (a + b) ** 2, 1.5);
const nfwMassShape = x => x < 1e-3
  ? x * x * (0.5 + x * (-2 / 3 + x * (3 / 4 + x * (-4 / 5 + x * (5 / 6 - 6 / 7 * x)))))
  : Math.log1p(x) - x / (1 + x);
const haloForceR = (R, rs = NFW_RS) => R === 0 ? G / (2 * rs * rs) : G * nfwMassShape(R / rs) / (R * R);   // per 4πρ₀r_s³

/** Component masses / normalisations so that v_c(R₀) = V₀ with the given force fractions. */
export function normalise({ fBulge = 0.05, fDisk = 0.60, fHalo = 0.35, R = R0, v = V0 } = {}) {
  const Ftot = v * v / R;
  return {
    Mb: fBulge * Ftot / bulgeForceR(R),
    Md: fDisk * Ftot / diskForceR(R),
    Mh: fHalo * Ftot / haloForceR(R),             // = 4πρ₀ r_s³
  };
}

export const MW = normalise();

/** Potential Φ(R, z) in (km/s)², with each component switchable. */
export function potential(R, z, { bulge = true, disk = true, halo = true } = {}, P = MW) {
  const r = Math.hypot(R, z);
  let phi = 0;
  if (bulge) phi -= G * P.Mb / (r + HERNQUIST_A);
  if (disk) phi -= G * P.Md / Math.sqrt(R * R + (MN_A + Math.sqrt(z * z + MN_B * MN_B)) ** 2);
  // NFW: ln(1 + x)/r → 1/r_s as r → 0 (finite central potential)
  if (halo) {
    const x = r / NFW_RS;
    const shape = x < 1e-4 ? (1 - x / 2 + x * x / 3 - x * x * x / 4) / NFW_RS : Math.log1p(x) / r;
    phi -= G * P.Mh * shape;
  }
  return phi;
}

/**
 * Acceleration (−∇Φ) at (x, y, z), in (km/s)²/kpc, written into out[0..2].
 * Exact analytic gradients of the three potentials.
 */
export function acceleration(x, y, z, comps, P, out) {
  const R2 = x * x + y * y, r2 = R2 + z * z, r = Math.sqrt(r2);
  let ax = 0, ay = 0, az = 0;
  if (comps.bulge && r > 0) {                       // force is finite and radial; zero exactly at r = 0
    const f = -G * P.Mb / (r * (r + HERNQUIST_A) ** 2);
    ax += f * x; ay += f * y; az += f * z;
  }
  if (comps.disk) {
    const zb = Math.sqrt(z * z + MN_B * MN_B), s = MN_A + zb;
    const D = Math.pow(R2 + s * s, 1.5);
    const f = -G * P.Md / D;
    ax += f * x; ay += f * y; az += f * z * s / zb;
  }
  if (comps.halo && r > 0) {
    const xr = r / NFW_RS;
    const f = -G * P.Mh * nfwMassShape(xr) / (r2 * r);
    ax += f * x; ay += f * y; az += f * z;
  }
  out[0] = ax; out[1] = ay; out[2] = az;
  return out;
}

/** Circular speed (km/s) in the midplane, total and per component. */
export function circularSpeed(R, comps = { bulge: true, disk: true, halo: true }, P = MW) {
  if (!(Number.isFinite(R) && R >= 0)) throw new Error('circularSpeed: radius must be finite and non-negative');
  const vb = comps.bulge ? Math.sqrt(R * P.Mb * bulgeForceR(R)) : 0;
  const vd = comps.disk ? Math.sqrt(R * P.Md * diskForceR(R)) : 0;
  const vh = comps.halo ? Math.sqrt(R * P.Mh * haloForceR(R)) : 0;
  return { total: Math.sqrt(vb * vb + vd * vd + vh * vh), bulge: vb, disk: vd, halo: vh };
}

/** NFW enclosed mass M(<r) = 4πρ₀r_s³ [ln(1 + x) − x/(1 + x)]. */
export const haloMass = (r, P = MW) => P.Mh * nfwMassShape(r / NFW_RS);
export const NFW_SCALE = NFW_RS;
/** NFW density ρ(r) = ρ₀ / [(r/r_s)(1 + r/r_s)²], with ρ₀ = Mh / (4π r_s³). */
export const haloDensity = (r, P = MW) => P.Mh / (4 * Math.PI * NFW_RS ** 3) / ((r / NFW_RS) * (1 + r / NFW_RS) ** 2);

/**
 * Measured Milky Way rotation curve, linear fit 5–25 kpc (Eilers et al. 2019, ApJ 871, 120):
 * v_c(R) = 229.0 − 1.7 (R − 8.122) km/s, systematic uncertainty ≈ 2–5 %.
 */
export const eilers2019 = R => 229.0 - 1.7 * (R - 8.122);

// ---------------------------------------------------------------- star ensemble

/**
 * A disk of test-particle stars integrated with kick-drift-kick leapfrog (symplectic;
 * conserves each star's L_z exactly in an axisymmetric potential, energy to O(Δt²)).
 * velocityModel: 'observed' gives every star the circular speed of the full model with dark
 * matter, regardless of the potential the stars actually feel. These are model reference
 * speeds; the measured rotation curve is compared separately via eilers2019.
 */
export function createDisk({ n = 3000, comps, velocityModel = 'observed', dispersion = 8, Rd = 2.6, Rmax = 25, rand }) {
  const x = new Float64Array(3 * n), v = new Float64Array(3 * n), a = new Float64Array(3 * n);
  const R0s = new Float64Array(n);
  const full = { bulge: true, disk: true, halo: true };
  for (let i = 0; i < n; i++) {
    // exponential surface density Σ ∝ e^{−R/Rd}: R ~ Gamma(2, Rd), truncated to [1, Rmax]
    let R;
    do { R = -Rd * Math.log(rand() * rand()); } while (R < 1 || R > Rmax);
    const ph = 2 * Math.PI * rand();
    const vc = circularSpeed(R, velocityModel === 'observed' ? full : comps).total;
    x[3 * i] = R * Math.cos(ph); x[3 * i + 1] = R * Math.sin(ph); x[3 * i + 2] = 0.1 * rand.normal();
    v[3 * i] = -vc * Math.sin(ph) + dispersion * rand.normal();
    v[3 * i + 1] = vc * Math.cos(ph) + dispersion * rand.normal();
    v[3 * i + 2] = 0.5 * dispersion * rand.normal();
    R0s[i] = R;
  }
  const disk = { n, x, v, a, R0s, comps, t: 0, steps: 0 };
  computeAcc(disk);
  disk.E0 = energies(disk);
  disk.L0 = angularMomenta(disk);
  return disk;
}

function computeAcc(d) {
  const tmp = new Float64Array(3);
  for (let i = 0; i < d.n; i++) {
    acceleration(d.x[3 * i], d.x[3 * i + 1], d.x[3 * i + 2], d.comps, MW, tmp);
    d.a[3 * i] = tmp[0]; d.a[3 * i + 1] = tmp[1]; d.a[3 * i + 2] = tmp[2];
  }
}

/** One KDK leapfrog step of size dt (time unit kpc/(km/s)). */
export function stepDisk(d, dt) {
  const { x, v, a } = d, N = 3 * d.n;
  for (let k = 0; k < N; k++) v[k] += 0.5 * dt * a[k];
  for (let k = 0; k < N; k++) x[k] += dt * v[k];
  computeAcc(d);
  for (let k = 0; k < N; k++) v[k] += 0.5 * dt * a[k];
  d.t += dt; d.steps++;
}

export function energies(d) {
  const E = new Float64Array(d.n);
  for (let i = 0; i < d.n; i++) {
    const X = d.x[3 * i], Y = d.x[3 * i + 1], Z = d.x[3 * i + 2];
    const v2 = d.v[3 * i] ** 2 + d.v[3 * i + 1] ** 2 + d.v[3 * i + 2] ** 2;
    E[i] = 0.5 * v2 + potential(Math.hypot(X, Y), Z, d.comps);
  }
  return E;
}

export function angularMomenta(d) {
  const L = new Float64Array(d.n);
  for (let i = 0; i < d.n; i++) L[i] = d.x[3 * i] * d.v[3 * i + 1] - d.x[3 * i + 1] * d.v[3 * i];
  return L;
}
