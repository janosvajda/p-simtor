// Self-gravitating stellar disk: a particle-mesh (PM) N-body model of a razor-thin disk of stars
// that attract each other, embedded in the rigid dark-matter halo (and optional bulge) of
// core/galaxy.mjs. Spiral arms then form by themselves through swing amplification of particle
// noise (Julian & Toomre 1966; Toomre 1981; Sellwood & Carlberg 1984, 2014; Binney & Tremaine 2008,
// "BT08", §6). Plain ES module, browser + Node, deterministic for a seed, no dependencies.
//
// Units (as core/galaxy.mjs): length kpc, velocity km/s, mass M☉, G = 4.30091727e-6 kpc (km/s)²/M☉,
// time unit 1 kpc/(km/s) = TIME_UNIT_MYR = 977.8 Myr. Disk particles move in the plane z = 0.
//
// Gravity of the disk (isolated, non-periodic): cloud-in-cell (CIC) mass assignment on an N×N mesh,
// convolution with a Plummer-softened Green's function by FFT on the zero-padded 2N×2N mesh
// (Hockney & Eastwood 1988, §6.5: the padding removes every periodic image exactly), CIC
// interpolation of the force. The force kernel g(d) = −G d/(|d|² + ε²)^{3/2} is convolved directly
// (no finite-difference gradient), so the mesh force is antisymmetric: momentum is conserved to
// round-off and there is no self-force. Plummer softening ε stands in for the disk thickness
// (a Plummer-softened thin disk acts in its plane like the same disk seen from height ε).
// Particles outside the mesh contribute no mass and feel the monopole of the mesh mass
// (counted in `outside`; enlarge `box` if it is not 0).
//
// Mesh layout follows core/fields.mjs: axis 0 fastest, index i + N·j; node (i, j) sits at
// x = −box/2 + (i + ½)h, y = −box/2 + (j + ½)h with h = box/N.

import { fft, rng, besselI } from './numerics.mjs';
import { fftND } from './fields.mjs';
import { besselK } from './special.mjs';
import { G, MW, acceleration as galaxyAcceleration, potential as galaxyPotential, circularSpeed } from './galaxy.mjs';

export { G };

const TOOMRE = 3.36;                                  // 2π·0.5348: stellar disk with a Schwarzschild velocity distribution (Toomre 1964)

const need = (cond, message) => { if (!cond) throw new Error(message); };
const positiveFinite = (value, name, who) => need(Number.isFinite(value) && value > 0, `${who}: ${name} must be a positive finite number (got ${value})`);

// ---------------------------------------------------------------- analytic references

/**
 * Kuzmin (1956) disk surface density Σ(R) = M a / [2π (R² + a²)^{3/2}] in M☉/kpc² (BT08 §2.3).
 * M > 0 total mass (M☉), a > 0 scale length (kpc). Throws on invalid input.
 */
export function kuzminSurfaceDensity(R, M, a) {
  need(Number.isFinite(R) && R >= 0, 'kuzminSurfaceDensity: R must be finite and ≥ 0');
  positiveFinite(M, 'M', 'kuzminSurfaceDensity'); positiveFinite(a, 'a', 'kuzminSurfaceDensity');
  return M * a / (2 * Math.PI * (R * R + a * a) ** 1.5);
}

/**
 * In-plane radial acceleration of a Kuzmin disk, a_R = −G M R / (R² + (a + ε)²)^{3/2} in (km/s)²/kpc.
 * Exact also with Plummer softening ε ≥ 0: the Kuzmin potential −GM/√(R² + (a + |z|)²) evaluated at
 * height z = ε is the Plummer-softened in-plane potential. Reference for the PM solver test.
 */
export function kuzminRadialAcceleration(R, M, a, softening = 0) {
  need(Number.isFinite(R) && R >= 0, 'kuzminRadialAcceleration: R must be finite and ≥ 0');
  positiveFinite(M, 'M', 'kuzminRadialAcceleration'); positiveFinite(a, 'a', 'kuzminRadialAcceleration');
  need(Number.isFinite(softening) && softening >= 0, 'kuzminRadialAcceleration: softening must be finite and ≥ 0');
  const s = a + softening;
  return -G * M * R / (R * R + s * s) ** 1.5;
}

/**
 * Circular speed (km/s) of a razor-thin exponential disk Σ = M/(2π R_d²) e^{−R/R_d}, unsoftened
 * (Freeman 1970; BT08 §2.6): v_c² = 4πGΣ₀R_d y² [I₀(y)K₀(y) − I₁(y)K₁(y)], y = R/(2R_d).
 * M: total (untruncated) mass in M☉, Rd in kpc. Not a spherical-mass rotation curve: it peaks
 * at R ≈ 2.2 R_d about 15 % above the spherical estimate.
 */
export function exponentialDiskCircularSpeed(R, M, Rd) {
  need(Number.isFinite(R) && R >= 0, 'exponentialDiskCircularSpeed: R must be finite and ≥ 0');
  positiveFinite(M, 'M', 'exponentialDiskCircularSpeed'); positiveFinite(Rd, 'Rd', 'exponentialDiskCircularSpeed');
  if (R === 0) return 0;
  const y = R / (2 * Rd);
  need(y < 300, 'exponentialDiskCircularSpeed: R/Rd too large for the Bessel products (use the point-mass limit)');
  const Sigma0 = M / (2 * Math.PI * Rd * Rd);
  const v2 = 4 * Math.PI * G * Sigma0 * Rd * y * y * (besselI(0, y) * besselK(0, y) - besselI(1, y) * besselK(1, y));
  return Math.sqrt(Math.max(v2, 0));
}

/**
 * Epicyclic frequency κ (km/s/kpc) from a circular-speed function vc(R) (km/s):
 * κ² = (1/R) d(v_c²)/dR + 2 v_c²/R² (BT08 §3.2.3), derivative by a central difference of width
 * 2·dR (default 1e-3 R). Checks: Kepler κ = Ω, flat rotation curve κ = √2 Ω, solid body κ = 2Ω.
 */
export function epicyclicFrequency(vc, R, dR = 1e-3 * R) {
  need(typeof vc === 'function', 'epicyclicFrequency: vc must be a function of R');
  positiveFinite(R, 'R', 'epicyclicFrequency'); positiveFinite(dR, 'dR', 'epicyclicFrequency');
  need(dR < R, 'epicyclicFrequency: dR must be smaller than R');
  const v2 = vc(R) ** 2, d = (vc(R + dR) ** 2 - vc(R - dR) ** 2) / (2 * dR);
  const k2 = d / R + 2 * v2 / (R * R);
  need(Number.isFinite(k2) && k2 >= 0, 'epicyclicFrequency: κ² is negative or non-finite (unstable rotation curve)');
  return Math.sqrt(k2);
}

/**
 * Toomre (1964) stability parameter of a stellar disk Q = σ_R κ / (3.36 G Σ) (BT08 §6.2.3).
 * sigmaR km/s, kappa km/s/kpc, Sigma M☉/kpc². Q < 1: locally unstable to axisymmetric collapse;
 * 1 ≲ Q ≲ 2: strong swing amplification of non-axisymmetric (spiral) disturbances.
 */
export function toomreQ(sigmaR, kappa, Sigma) {
  need(Number.isFinite(sigmaR) && sigmaR >= 0, 'toomreQ: sigmaR must be finite and ≥ 0');
  need(Number.isFinite(kappa) && kappa >= 0, 'toomreQ: kappa must be finite and ≥ 0');
  positiveFinite(Sigma, 'Sigma', 'toomreQ');
  return sigmaR * kappa / (TOOMRE * G * Sigma);
}

/**
 * Arm number most amplified by swing amplification, m = κ² R / (2π G Σ X) with X = λ_azimuthal/λ_crit
 * ≈ 1–3, peak near X ≈ 2 for a flat rotation curve and Q ≈ 1.5 (Toomre 1981; BT08 §6.3).
 * A rough guide (±1 arm): lighter disks (relative to the halo) make more, thinner arms.
 */
export function swingArmNumber(R, kappa, Sigma, X = 2) {
  positiveFinite(R, 'R', 'swingArmNumber'); positiveFinite(kappa, 'kappa', 'swingArmNumber');
  positiveFinite(Sigma, 'Sigma', 'swingArmNumber'); positiveFinite(X, 'X', 'swingArmNumber');
  return kappa * kappa * R / (2 * Math.PI * G * Sigma * X);
}

// ---------------------------------------------------------------- particle-mesh solver

/**
 * Isolated (non-periodic) 2D particle-mesh gravity solver with Plummer softening.
 * Options: grid N (power of two, 8–2048; the FFT mesh is 2N×2N), box (side length, kpc, centred on
 * the origin), softening ε > 0 (kpc; choose ε ≳ h = box/N), Gconst (default galaxy.mjs G).
 * Methods (all allocation-free after construction):
 *   accelerations(x, y, m, n, ax, ay)  deposit + solve + interpolate; overwrites ax, ay ((km/s)²/kpc)
 *   potentials(x, y, m, n, out)        softened potential Φ_self at the particles ((km/s)²; includes the
 *                                      particle's own smoothed cloud, a near-constant offset)
 *   deposit(x, y, m, n) / loadMesh(massGrid), then solveForces() / solvePotential(), then
 *   fieldAt(px, py, out) → out[0..2] = [a_x, a_y, Φ]; Φ only after solvePotential().
 * Accuracy: a few % beyond ~2 cells from a point mass, exact antisymmetry (Newton's third law)
 * to round-off. Mesh node positions: x_i = −box/2 + (i + ½)h.
 */
export function createPMSolver({ grid = 128, box, softening, Gconst = G } = {}) {
  const who = 'createPMSolver';
  need(Number.isInteger(grid) && grid >= 8 && grid <= 2048 && (grid & (grid - 1)) === 0, `${who}: grid must be a power of two in [8, 2048]`);
  positiveFinite(box, 'box', who); positiveFinite(softening, 'softening', who); positiveFinite(Gconst, 'Gconst', who);
  const N = grid, M = 2 * N, MM = M * M, h = box / N, invh = 1 / h, x0 = -box / 2 + h / 2, eps2 = softening * softening;
  // Green's functions on the padded mesh: displacement index wraps (k < N → k, k > N → k − M);
  // the never-used displacement ±N is set to zero so the kernels stay exactly odd/even.
  const kRe = new Float64Array(MM), kIm = new Float64Array(MM), pRe = new Float64Array(MM), pIm = new Float64Array(MM);
  for (let j = 0; j < M; j++) {
    const dy = (j < N ? j : j - M) * h;
    for (let i = 0; i < M; i++) {
      if (i === N || j === N) continue;
      const dx = (i < N ? i : i - M) * h, s = 1 / Math.sqrt(dx * dx + dy * dy + eps2), k = i + M * j;
      kRe[k] = -Gconst * dx * s * s * s; kIm[k] = -Gconst * dy * s * s * s; pRe[k] = -Gconst * s;
    }
  }
  fftND(kRe, kIm, [M, M]);                                   // FFT(g_x) + i FFT(g_y): one inverse FFT gives a_x + i a_y
  fftND(pRe, pIm, [M, M]);                                   // even kernel → real transform (imaginary part is round-off)
  const re = new Float64Array(MM), im = new Float64Array(MM), rowRe = new Float64Array(M), rowIm = new Float64Array(M);
  const viewsRe = Array.from({ length: M }, (_, j) => re.subarray(j * M, (j + 1) * M));
  const viewsIm = Array.from({ length: M }, (_, j) => im.subarray(j * M, (j + 1) * M));
  const fx = new Float64Array(N * N), fy = new Float64Array(N * N), phi = new Float64Array(N * N);
  const tmp = new Float64Array(3);
  const columns = (count, inverse) => {
    for (let i = 0; i < count; i++) {
      for (let j = 0, k = i; j < M; j++, k += M) { rowRe[j] = re[k]; rowIm[j] = im[k]; }
      fft(rowRe, rowIm, inverse);
      for (let j = 0, k = i; j < M; j++, k += M) { re[k] = rowRe[j]; im[k] = rowIm[j]; }
    }
  };
  // pruned 2D FFTs: only the first N rows hold mass; only the first N×N outputs are needed
  const forward = () => { for (let j = 0; j < N; j++) fft(viewsRe[j], viewsIm[j], false); columns(M, false); };
  const inverse = () => { for (let j = 0; j < M; j++) fft(viewsRe[j], viewsIm[j], true); columns(N, true); };
  const pm = {
    N, M, h, box, softening, Gconst, x0,
    meshMass: 0, comX: 0, comY: 0, outside: 0, potentialReady: false,
    /** CIC mass assignment of n particles; returns the mass on the mesh. */
    deposit(x, y, m, n) {
      re.fill(0); im.fill(0);
      let mass = 0, sx = 0, sy = 0, out = 0;
      for (let p = 0; p < n; p++) {
        const fxp = (x[p] - x0) * invh, fyp = (y[p] - x0) * invh;
        if (!(fxp >= 0 && fxp < N - 1 && fyp >= 0 && fyp < N - 1)) { out++; continue; }
        const i = fxp | 0, j = fyp | 0, tx = fxp - i, ty = fyp - j, mp = m[p], k = i + M * j;
        re[k] += mp * (1 - tx) * (1 - ty); re[k + 1] += mp * tx * (1 - ty);
        re[k + M] += mp * (1 - tx) * ty; re[k + M + 1] += mp * tx * ty;
        mass += mp; sx += mp * x[p]; sy += mp * y[p];
      }
      pm.meshMass = mass; pm.comX = mass > 0 ? sx / mass : 0; pm.comY = mass > 0 ? sy / mass : 0; pm.outside = out;
      pm.potentialReady = false;
      return mass;
    },
    /** Put node masses (N×N, index i + N·j, M☉) on the mesh directly, e.g. an analytic surface density × h². */
    loadMesh(massGrid) {
      need(massGrid && massGrid.length === N * N, `${who}: loadMesh needs N×N node masses`);
      re.fill(0); im.fill(0);
      let mass = 0, sx = 0, sy = 0;
      for (let j = 0; j < N; j++) {for (let i = 0; i < N; i++) {
        const mp = massGrid[i + N * j];
        need(Number.isFinite(mp), `${who}: loadMesh masses must be finite`);
        re[i + M * j] = mp; mass += mp; sx += mp * (x0 + i * h); sy += mp * (x0 + j * h);
      }}
      pm.meshMass = mass; pm.comX = mass ? sx / mass : 0; pm.comY = mass ? sy / mass : 0; pm.outside = 0;
      pm.potentialReady = false;
      return mass;
    },
    /** Mesh accelerations from the deposited mass (one forward and one inverse pruned FFT). */
    solveForces() {
      forward();
      for (let k = 0; k < MM; k++) {
        const a = re[k], b = im[k], c = kRe[k], d = kIm[k];
        re[k] = a * c - b * d; im[k] = a * d + b * c;
      }
      inverse();
      for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) { fx[i + N * j] = re[i + M * j]; fy[i + N * j] = im[i + M * j]; }
    },
    /** Mesh potential from the deposited mass. Call deposit/loadMesh again before solveForces. */
    solvePotential() {
      forward();
      for (let k = 0; k < MM; k++) { re[k] *= pRe[k]; im[k] *= pRe[k]; }
      inverse();
      for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) phi[i + N * j] = re[i + M * j];
      pm.potentialReady = true;
    },
    /** CIC-interpolated [a_x, a_y, Φ] at a point; outside the mesh the monopole of the mesh mass. */
    fieldAt(px, py, out = tmp) {
      const fxp = (px - x0) * invh, fyp = (py - x0) * invh;
      if (fxp >= 0 && fxp < N - 1 && fyp >= 0 && fyp < N - 1) {
        const i = fxp | 0, j = fyp | 0, tx = fxp - i, ty = fyp - j, k = i + N * j;
        const w00 = (1 - tx) * (1 - ty), w10 = tx * (1 - ty), w01 = (1 - tx) * ty, w11 = tx * ty;
        out[0] = w00 * fx[k] + w10 * fx[k + 1] + w01 * fx[k + N] + w11 * fx[k + N + 1];
        out[1] = w00 * fy[k] + w10 * fy[k + 1] + w01 * fy[k + N] + w11 * fy[k + N + 1];
        out[2] = w00 * phi[k] + w10 * phi[k + 1] + w01 * phi[k + N] + w11 * phi[k + N + 1];
      } else {
        const dx = px - pm.comX, dy = py - pm.comY, s = 1 / Math.sqrt(dx * dx + dy * dy + eps2), f = -Gconst * pm.meshMass * s * s * s;
        out[0] = f * dx; out[1] = f * dy; out[2] = -Gconst * pm.meshMass * s;
      }
      return out;
    },
    accelerations(x, y, m, n, ax, ay) {
      pm.deposit(x, y, m, n); pm.solveForces();
      for (let p = 0; p < n; p++) { pm.fieldAt(x[p], y[p], tmp); ax[p] = tmp[0]; ay[p] = tmp[1]; }
    },
    potentials(x, y, m, n, out) {
      pm.deposit(x, y, m, n); pm.solvePotential();
      for (let p = 0; p < n; p++) out[p] = pm.fieldAt(x[p], y[p], tmp)[2];
    },
  };
  return pm;
}

// ---------------------------------------------------------------- the self-gravitating disk

/**
 * A razor-thin exponential stellar disk with self-gravity (PM) in a rigid halo/bulge from galaxy.mjs.
 *
 * Options (units kpc, km/s, M☉):
 *   n (30000) particles of equal mass diskMass/n; Rd (2.6) scale length; Rmax (5 Rd) truncation radius;
 *   diskMass (3e10); Q (1.5) target Toomre Q, constant with radius; halo (true) / bulge (false):
 *   rigid NFW halo / Hernquist bulge of galaxy.mjs with normalisation P (default MW; the galaxy.mjs
 *   Miyamoto–Nagai disk is replaced by the live disk); grid (128), box (2.6 Rmax), softening
 *   (max(h, 0.1 Rd)); selfGravity (true; false = test particles in the rigid potential only); seed.
 *
 * Initial conditions (BT08 §4.8.2): positions from Σ ∝ e^{−R/Rd}; the circular speed v_c(R) is the
 * mesh force of the smooth analytic disk (same softening and mesh as the run) plus the external
 * component; κ(R) from v_c; σ_R = Q·3.36GΣ/κ (Toomre 1964), σ_φ = σ_R κ/(2Ω) (epicyclic ratio), mean
 * streaming from the asymmetric-drift (Jeans) equation v̄_φ² = v_c² + σ_R²(1 − κ²/4Ω² − R/R_d +
 * ∂lnσ_R²/∂lnR), clamped at 0. Gaussian velocities. The centre of mass and the mean velocity are
 * set to zero. The epicyclic approximation is poor where σ_R ≳ v_c/3 (near the centre of a hot disk):
 * the disk relaxes a little at the start; measure Q(R) with toomreQProfile instead of assuming it.
 *
 * Returns a state { n, x, y, vx, vy, ax, ay, m, t, steps, pm, table, ... } with
 *   step(dt)  one kick-drift-kick leapfrog step (dt in kpc/(km/s)); advance(dt, k) = k steps;
 *   energy() { kinetic, self, external, total } in M☉ (km/s)²; momentum() [p_x, p_y];
 *   angularMomentum() L_z; centreOfMass() [x, y]; suggestedDt; outside (particles off the mesh).
 * Every step is the same arithmetic whatever the frame schedule; throws if a position becomes
 * non-finite. table = { R, vc, kappa, Sigma, sigmaR, vphi } is the initial-condition model.
 */
export function createSelfGravDisk({
  n = 30000, Rd = 2.6, Rmax, diskMass = 3e10, Q = 1.5, halo = true, bulge = false, P = MW,
  grid = 128, box, softening, selfGravity = true, seed = 1,
} = {}) {
  const who = 'createSelfGravDisk';
  need(Number.isInteger(n) && n >= 1, `${who}: n must be a positive integer`);
  positiveFinite(Rd, 'Rd', who);
  const rMax = Rmax ?? 5 * Rd;
  positiveFinite(rMax, 'Rmax', who); positiveFinite(diskMass, 'diskMass', who); positiveFinite(Q, 'Q', who);
  const side = box ?? 2.6 * rMax;
  positiveFinite(side, 'box', who);
  need(side > 2.05 * rMax, `${who}: box must exceed 2·Rmax so the whole disk lies on the mesh`);
  const eps = softening ?? Math.max(side / grid, 0.1 * Rd);
  const comps = { bulge: !!bulge, disk: false, halo: !!halo };
  const pm = createPMSolver({ grid, box: side, softening: eps });
  const { N, h, x0 } = pm;

  // analytic truncated exponential disk
  const xm = rMax / Rd, Sigma0 = diskMass / (2 * Math.PI * Rd * Rd * (1 - Math.exp(-xm) * (1 + xm)));
  const SigmaOf = R => (R <= rMax ? Sigma0 * Math.exp(-R / Rd) : 0);

  // v_c(R) of the smooth disk on this mesh: ring averages of the mesh force
  const dR = h / 2, K = Math.ceil(rMax / dR) + 4;
  const T = { R: new Float64Array(K), vc: new Float64Array(K), kappa: new Float64Array(K), Sigma: new Float64Array(K), sigmaR: new Float64Array(K), sigmaPhi: new Float64Array(K), vphi: new Float64Array(K) };
  const vc2 = new Float64Array(K);
  if (selfGravity) {
    const mesh = new Float64Array(N * N);
    let total = 0;
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) { const v = SigmaOf(Math.hypot(x0 + i * h, x0 + j * h)) * h * h; mesh[i + N * j] = v; total += v; }
    for (let k = 0; k < N * N; k++) mesh[k] *= diskMass / total;
    pm.loadMesh(mesh); pm.solveForces();
  }
  const out = new Float64Array(3), nA = 64;
  for (let k = 0; k < K; k++) {
    const R = (k + 0.5) * dR;
    let aR = 0;
    if (selfGravity) {
      for (let q = 0; q < nA; q++) {
        const c = Math.cos(2 * Math.PI * (q + 0.5) / nA), s = Math.sin(2 * Math.PI * (q + 0.5) / nA);
        pm.fieldAt(R * c, R * s, out); aR += out[0] * c + out[1] * s;
      }
      aR /= nA;
    }
    T.R[k] = R;
    vc2[k] = Math.max(0, -R * aR + circularSpeed(R, comps, P).total ** 2);
  }
  for (let k = 0; k < K; k++) {
    const R = T.R[k], a = Math.max(0, k - 2), b = Math.min(K - 1, k + 2);
    const dv2 = (vc2[b] - vc2[a]) / (T.R[b] - T.R[a]);
    T.kappa[k] = Math.sqrt(Math.max(dv2 / R + 2 * vc2[k] / (R * R), 1e-12));
    T.vc[k] = Math.sqrt(vc2[k]);
    T.Sigma[k] = SigmaOf(R);
    T.sigmaR[k] = T.Sigma[k] > 0 ? Q * TOOMRE * G * T.Sigma[k] / T.kappa[k] : 0;
  }
  for (let k = 0; k < K; k++) {
    const R = T.R[k], Omega2 = vc2[k] / (R * R), a = Math.max(0, k - 2), b = Math.min(K - 1, k + 2);
    const s2a = T.sigmaR[a] ** 2, s2b = T.sigmaR[b] ** 2;
    const dlns2 = s2a > 0 && s2b > 0 ? Math.log(s2b / s2a) / Math.log(T.R[b] / T.R[a]) : 0;
    const ratio = Omega2 > 0 ? T.kappa[k] ** 2 / (4 * Omega2) : 1;
    T.sigmaPhi[k] = T.sigmaR[k] * Math.sqrt(ratio);
    T.vphi[k] = Math.sqrt(Math.max(0, vc2[k] + T.sigmaR[k] ** 2 * (1 - ratio - R / Rd + dlns2)));
  }
  const lookup = (arr, R) => {
    const f = Math.min(Math.max(R / dR - 0.5, 0), K - 1.000001), k = Math.floor(f), t = f - k;
    return arr[k] * (1 - t) + arr[k + 1] * t;
  };

  // particles (structure of arrays)
  const rand = rng(seed);
  const x = new Float64Array(n), y = new Float64Array(n), vx = new Float64Array(n), vy = new Float64Array(n);
  const ax = new Float64Array(n), ay = new Float64Array(n), m = new Float64Array(n).fill(diskMass / n);
  for (let p = 0; p < n; p++) {
    let R;
    do { R = -Rd * Math.log((1 - rand()) * (1 - rand())); } while (R > rMax);   // Gamma(2, Rd) ∝ R e^{−R/Rd}
    const ph = 2 * Math.PI * rand(), c = Math.cos(ph), s = Math.sin(ph);
    const vR = lookup(T.sigmaR, R) * rand.normal(), vP = lookup(T.vphi, R) + lookup(T.sigmaPhi, R) * rand.normal();
    x[p] = R * c; y[p] = R * s; vx[p] = vR * c - vP * s; vy[p] = vR * s + vP * c;
  }
  let cx = 0, cy = 0, cvx = 0, cvy = 0;
  for (let p = 0; p < n; p++) { cx += x[p]; cy += y[p]; cvx += vx[p]; cvy += vy[p]; }
  cx /= n; cy /= n; cvx /= n; cvy /= n;
  for (let p = 0; p < n; p++) { x[p] -= cx; y[p] -= cy; vx[p] -= cvx; vy[p] -= cvy; }

  const acc = new Float64Array(3), phiSelf = new Float64Array(n);
  const computeAcc = () => {
    if (selfGravity) pm.accelerations(x, y, m, n, ax, ay);
    else { ax.fill(0); ay.fill(0); }
    if (comps.halo || comps.bulge) {
      for (let p = 0; p < n; p++) { galaxyAcceleration(x[p], y[p], 0, comps, P, acc); ax[p] += acc[0]; ay[p] += acc[1]; }
    }
    s.outside = selfGravity ? pm.outside : 0;
  };
  // orbital period at R = max(Rd/4, 2ε), divided by 50 (energy error ≲ 1e-4 per outer rotation without a bulge)
  const Rin = Math.max(Rd / 4, 2 * eps), vin = Math.max(lookup(T.vc, Rin), 1e-6);
  const s = {
    n, x, y, vx, vy, ax, ay, m, t: 0, steps: 0, pm, table: T, comps, P, Rd, Rmax: rMax, diskMass, Q,
    softening: eps, box: side, grid, selfGravity, outside: 0,
    suggestedDt: 2 * Math.PI * Rin / vin / 50,
    step(dt) {
      need(Number.isFinite(dt) && dt > 0, `${who}: step needs a positive finite dt`);
      const hdt = 0.5 * dt;
      let check = 0;
      for (let p = 0; p < n; p++) {
        vx[p] += hdt * ax[p]; vy[p] += hdt * ay[p];
        x[p] += dt * vx[p]; y[p] += dt * vy[p]; check += x[p] + y[p];
      }
      need(Number.isFinite(check), `${who}: a particle position became non-finite at t = ${s.t}`);
      computeAcc();
      for (let p = 0; p < n; p++) { vx[p] += hdt * ax[p]; vy[p] += hdt * ay[p]; }
      s.t += dt; s.steps++;
      return s;
    },
    advance(dt, count = 1) {
      need(Number.isInteger(count) && count >= 0, `${who}: advance needs a non-negative integer step count`);
      for (let k = 0; k < count; k++) s.step(dt);
      return s;
    },
    /** Energies in M☉ (km/s)²; self = ½ Σ m Φ_self (CIC-interpolated mesh potential). */
    energy() {
      let kin = 0, self = 0, ext = 0;
      for (let p = 0; p < n; p++) kin += 0.5 * m[p] * (vx[p] * vx[p] + vy[p] * vy[p]);
      if (selfGravity) { pm.potentials(x, y, m, n, phiSelf); for (let p = 0; p < n; p++) self += 0.5 * m[p] * phiSelf[p]; }
      if (comps.halo || comps.bulge) for (let p = 0; p < n; p++) ext += m[p] * galaxyPotential(Math.hypot(x[p], y[p]), 0, comps, P);
      return { kinetic: kin, self, external: ext, total: kin + self + ext };
    },
    momentum() { let px = 0, py = 0; for (let p = 0; p < n; p++) { px += m[p] * vx[p]; py += m[p] * vy[p]; } return [px, py]; },
    angularMomentum() { let L = 0; for (let p = 0; p < n; p++) L += m[p] * (x[p] * vy[p] - y[p] * vx[p]); return L; },
    centreOfMass() { let mx = 0, my = 0, mt = 0; for (let p = 0; p < n; p++) { mx += m[p] * x[p]; my += m[p] * y[p]; mt += m[p]; } return [mx / mt, my / mt]; },
  };
  computeAcc();
  return s;
}

// ---------------------------------------------------------------- measurements

const checkSet = (s, who) => need(s && Number.isInteger(s.n) && s.x?.length >= s.n && s.y?.length >= s.n && s.m?.length >= s.n, `${who}: expected a particle set { n, x, y, m, ... }`);
const radialBins = (nR, Rmax, who) => {
  need(Number.isInteger(nR) && nR >= 1, `${who}: nR must be a positive integer`);
  positiveFinite(Rmax, 'Rmax', who);
  return Float64Array.from({ length: nR }, (_, k) => (k + 0.5) * Rmax / nR);
};

/**
 * Surface density Σ (M☉/kpc²) on a polar grid: nR rings of equal width out to Rmax, nPhi sectors
 * from φ = 0 counter-clockwise. Returns { Sigma (index iR + nR·iPhi), R, nR, nPhi, Rmax }.
 */
export function surfaceDensityPolar(s, { nR = 32, nPhi = 64, Rmax = s?.Rmax } = {}) {
  const who = 'surfaceDensityPolar';
  checkSet(s, who); const R = radialBins(nR, Rmax, who);
  need(Number.isInteger(nPhi) && nPhi >= 1, `${who}: nPhi must be a positive integer`);
  const Sigma = new Float64Array(nR * nPhi), dR = Rmax / nR, dphi = 2 * Math.PI / nPhi;
  for (let p = 0; p < s.n; p++) {
    const r = Math.hypot(s.x[p], s.y[p]);
    if (!(r < Rmax)) continue;
    let ph = Math.atan2(s.y[p], s.x[p]); if (ph < 0) ph += 2 * Math.PI;
    Sigma[Math.floor(r / dR) + nR * Math.min(nPhi - 1, Math.floor(ph / dphi))] += s.m[p];
  }
  for (let i = 0; i < nR; i++) { const area = R[i] * dR * dphi; for (let j = 0; j < nPhi; j++) Sigma[i + nR * j] /= area; }
  return { Sigma, R, nR, nPhi, Rmax };
}

/**
 * Fourier amplitudes of the azimuthal mass distribution in each of nR rings out to Rmax,
 * A_m(R) = |Σ_j m_j e^{−imφ_j}| / Σ_j m_j (0 for an axisymmetric disk; shot noise ≈ 1/√N_ring), and
 * phase_m(R) = arg Σ_j m_j e^{−imφ_j} (a spiral cos(mφ + p ln R) gives phase = p ln R + const).
 * annulus[m] = √(mass-weighted mean of A_m² over the rings centred in [Rin, Rout)): a ring-by-ring
 * amplitude, so a wound spiral does not cancel itself as a phase-coherent annulus sum would.
 * noise = the same average for pure shot noise, √(Σ_k w_k/N_k) for equal masses.
 * Returns { R, A (index iR·(mMax+1) + m), phase, annulus, dominant (m ≥ 1 with the largest annulus
 * amplitude), noise, mMax }.
 */
export function fourierModes(s, { mMax = 6, nR = 16, Rmax = s?.Rmax, Rin = 0, Rout = Rmax } = {}) {
  const who = 'fourierModes';
  checkSet(s, who); const R = radialBins(nR, Rmax, who);
  need(Number.isInteger(mMax) && mMax >= 1 && mMax <= 64, `${who}: mMax must be an integer in [1, 64]`);
  need(Rin >= 0 && Rout > Rin, `${who}: need 0 ≤ Rin < Rout`);
  const W = mMax + 1, C = new Float64Array(nR * W), S = new Float64Array(nR * W), mass = new Float64Array(nR), count = new Uint32Array(nR);
  const dR = Rmax / nR;
  for (let p = 0; p < s.n; p++) {
    const r = Math.hypot(s.x[p], s.y[p]);
    if (!(r < Rmax) || r === 0) continue;
    const c1 = s.x[p] / r, s1 = -s.y[p] / r, mp = s.m[p], k = Math.floor(r / dR) * W;
    let c = 1, sn = 0;                                        // e^{−imφ} by recurrence
    for (let mm = 0; mm < W; mm++) {
      C[k + mm] += mp * c; S[k + mm] += mp * sn;
      const cn = c * c1 - sn * s1; sn = c * s1 + sn * c1; c = cn;
    }
    mass[k / W] += mp; count[k / W]++;
  }
  const A = new Float64Array(nR * W), phase = new Float64Array(nR * W), annulus = new Float64Array(W);
  let wsum = 0, noise2 = 0;
  for (let k = 0; k < nR; k++) {
    if (!(mass[k] > 0)) continue;
    const use = R[k] >= Rin && R[k] < Rout;
    if (use) { wsum += mass[k]; noise2 += mass[k] / count[k]; }
    for (let mm = 0; mm < W; mm++) {
      const i = k * W + mm;
      A[i] = Math.hypot(C[i], S[i]) / mass[k]; phase[i] = Math.atan2(S[i], C[i]);
      if (use) annulus[mm] += mass[k] * A[i] * A[i];
    }
  }
  for (let mm = 0; mm < W; mm++) annulus[mm] = wsum > 0 ? Math.sqrt(annulus[mm] / wsum) : 0;
  let dominant = 1;
  for (let mm = 2; mm < W; mm++) if (annulus[mm] > annulus[dominant]) dominant = mm;
  return { R, A, phase, annulus, dominant, noise: wsum > 0 ? Math.sqrt(noise2 / wsum) : Infinity, mMax };
}

/**
 * Pitch angle of the strongest m-armed logarithmic spiral in [Rin, Rout] (Kalnajs 1975;
 * Considère & Athanassoula 1988): A(p) = |Σ m_j e^{−i(mφ_j + p ln R_j)}| / Σ m_j is maximised over
 * p ∈ [−pMax, pMax]; pitch i = atan(m/p). pitchDeg > 0: trailing arms (outer parts lag the rotation,
 * sense taken from the sign of L_z when velocities are given), < 0: leading. Returns
 * { m, p, pitchDeg, amplitude, trailing }. Cost O(N · pMax/dp).
 */
export function spiralPitch(s, { m = 2, Rin, Rout, pMax = 50, dp = 0.25 } = {}) {
  const who = 'spiralPitch';
  checkSet(s, who);
  need(Number.isInteger(m) && m >= 1, `${who}: m must be a positive integer`);
  need(Rin > 0 && Rout > Rin, `${who}: need 0 < Rin < Rout`);
  positiveFinite(pMax, 'pMax', who); positiveFinite(dp, 'dp', who);
  let L = 0;
  if (s.vx && s.vy) for (let q = 0; q < s.n; q++) L += s.m[q] * (s.x[q] * s.vy[q] - s.y[q] * s.vx[q]);
  const sense = L < 0 ? -1 : 1;
  const sel = [];
  for (let q = 0; q < s.n; q++) { const r = Math.hypot(s.x[q], s.y[q]); if (r >= Rin && r < Rout) sel.push(q); }
  need(sel.length > 0, `${who}: no particles in the annulus`);
  const mphi = Float64Array.from(sel, q => m * Math.atan2(s.y[q], s.x[q]) * sense);
  const lnR = Float64Array.from(sel, q => Math.log(Math.hypot(s.x[q], s.y[q])));
  const w = Float64Array.from(sel, q => s.m[q]);
  let wt = 0; for (const v of w) wt += v;
  let best = { p: 0, amplitude: -1 };
  const steps = Math.round(pMax / dp);
  for (let k = -steps; k <= steps; k++) {
    const p = k * dp;
    let c = 0, sn = 0;
    for (let q = 0; q < w.length; q++) { const a = mphi[q] + p * lnR[q]; c += w[q] * Math.cos(a); sn -= w[q] * Math.sin(a); }
    const amp = Math.hypot(c, sn) / wt;
    if (amp > best.amplitude) best = { p, amplitude: amp };
  }
  const pitchDeg = best.p === 0 ? 90 : Math.atan(m / best.p) * 180 / Math.PI;
  return { m, p: best.p, pitchDeg, amplitude: best.amplitude, trailing: best.p > 0 };
}

/**
 * Rotation curve measured from the particles in nR rings out to Rmax: mean streaming speed
 * v̄_φ (km/s), radial velocity dispersion σ_R, and the circular speed from the current total
 * accelerations v_c = √(−R ⟨a_R⟩) (needs s.ax, s.ay). Returns { R, vphi, sigmaR, vc, mass, count }.
 */
export function rotationCurve(s, { nR = 24, Rmax = s?.Rmax } = {}) {
  const who = 'rotationCurve';
  checkSet(s, who); const R = radialBins(nR, Rmax, who);
  need(s.vx && s.vy, `${who}: velocities required`);
  const mass = new Float64Array(nR), sv = new Float64Array(nR), svr = new Float64Array(nR), svr2 = new Float64Array(nR), sa = new Float64Array(nR), count = new Uint32Array(nR);
  const dR = Rmax / nR, hasA = !!(s.ax && s.ay);
  for (let p = 0; p < s.n; p++) {
    const r = Math.hypot(s.x[p], s.y[p]);
    if (!(r < Rmax) || r === 0) continue;
    const k = Math.floor(r / dR), c = s.x[p] / r, sn = s.y[p] / r, mp = s.m[p];
    const vR = s.vx[p] * c + s.vy[p] * sn, vP = -s.vx[p] * sn + s.vy[p] * c;
    mass[k] += mp; sv[k] += mp * vP; svr[k] += mp * vR; svr2[k] += mp * vR * vR; count[k]++;
    if (hasA) sa[k] += mp * (s.ax[p] * c + s.ay[p] * sn);
  }
  const vphi = new Float64Array(nR), sigmaR = new Float64Array(nR), vc = new Float64Array(nR).fill(NaN);
  for (let k = 0; k < nR; k++) {
    if (!(mass[k] > 0)) { vphi[k] = sigmaR[k] = NaN; continue; }
    vphi[k] = sv[k] / mass[k];
    const mvr = svr[k] / mass[k];
    sigmaR[k] = Math.sqrt(Math.max(0, svr2[k] / mass[k] - mvr * mvr));
    if (hasA) vc[k] = Math.sqrt(Math.max(0, -R[k] * sa[k] / mass[k]));
  }
  return { R, vphi, sigmaR, vc, mass, count };
}

/**
 * Toomre Q(R) measured from the particles: Σ and σ_R from rings, κ from the measured v_c(R)
 * (rotationCurve, central differences). Q = σ_R κ/(3.36GΣ). Noisy where rings hold few particles.
 * Returns { R, Q, Sigma, sigmaR, kappa, vc }.
 */
export function toomreQProfile(s, { nR = 24, Rmax = s?.Rmax } = {}) {
  const rc = rotationCurve(s, { nR, Rmax });
  const { R, vc, sigmaR, mass } = rc, dR = Rmax / nR;
  const Sigma = Float64Array.from(mass, (mk, k) => mk / (2 * Math.PI * R[k] * dR));
  const kappa = new Float64Array(nR).fill(NaN), Q = new Float64Array(nR).fill(NaN);
  for (let k = 0; k < nR; k++) {
    const a = Math.max(0, k - 1), b = Math.min(nR - 1, k + 1);
    const d = (vc[b] ** 2 - vc[a] ** 2) / (R[b] - R[a]), k2 = d / R[k] + 2 * vc[k] ** 2 / R[k] ** 2;
    if (Number.isFinite(k2) && k2 > 0) kappa[k] = Math.sqrt(k2);
    if (Sigma[k] > 0 && Number.isFinite(kappa[k])) Q[k] = sigmaR[k] * kappa[k] / (TOOMRE * G * Sigma[k]);
  }
  return { R, Q, Sigma, sigmaR, kappa, vc };
}

/**
 * Local surface density at every particle (M☉/kpc²) by CIC on a cells×cells Cartesian grid of side
 * 2·extent centred on the origin: for colouring a face-on view by density. Writes into out
 * (Float32Array or Float64Array of length ≥ n; allocated if omitted). The work grid can be reused
 * via opts.work (length cells²). Particles off the grid get 0.
 */
export function localDensity(s, { cells = 128, extent = s?.Rmax, out, work } = {}) {
  const who = 'localDensity';
  checkSet(s, who);
  need(Number.isInteger(cells) && cells >= 2, `${who}: cells must be an integer ≥ 2`);
  positiveFinite(extent, 'extent', who);
  const g = work ?? new Float64Array(cells * cells), o = out ?? new Float64Array(s.n);
  need(g.length >= cells * cells && o.length >= s.n, `${who}: out/work too short`);
  g.fill(0);
  const h = 2 * extent / cells, x0 = -extent + h / 2, inv = 1 / h, area = h * h;
  for (let pass = 0; pass < 2; pass++) {
    for (let p = 0; p < s.n; p++) {
      const fx = (s.x[p] - x0) * inv, fy = (s.y[p] - x0) * inv;
      if (!(fx >= 0 && fx < cells - 1 && fy >= 0 && fy < cells - 1)) { if (pass) o[p] = 0; continue; }
      const i = fx | 0, j = fy | 0, tx = fx - i, ty = fy - j, k = i + cells * j;
      const w00 = (1 - tx) * (1 - ty), w10 = tx * (1 - ty), w01 = (1 - tx) * ty, w11 = tx * ty;
      if (pass === 0) {
        const mp = s.m[p] / area;
        g[k] += mp * w00; g[k + 1] += mp * w10; g[k + cells] += mp * w01; g[k + cells + 1] += mp * w11;
      } else o[p] = w00 * g[k] + w10 * g[k + 1] + w01 * g[k + cells] + w11 * g[k + cells + 1];
    }
  }
  return o;
}
