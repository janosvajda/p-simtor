// Self-gravitating disk (core/disk-nbody.mjs): particle-mesh gravity against exact solutions,
// conservation laws, and the spontaneous growth of spiral structure by swing amplification.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  G, createPMSolver, createSelfGravDisk, kuzminSurfaceDensity, kuzminRadialAcceleration,
  exponentialDiskCircularSpeed, epicyclicFrequency, toomreQ, swingArmNumber,
  surfaceDensityPolar, fourierModes, spiralPitch, rotationCurve, toomreQProfile, localDensity,
} from '../core/disk-nbody.mjs';
import { rng } from '../core/numerics.mjs';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const tableAt = (d, R) => Math.min(d.table.R.length - 1, Math.round(R / (d.table.R[1] - d.table.R[0]) - 0.5));

test('analytic references: Freeman disk, epicyclic frequency, Toomre Q, swing arm number', () => {
  // Freeman (1970): v_c peaks at R ≈ 2.15 R_d (BT08 §2.6); Keplerian far outside
  const M = 5e10, Rd = 3;
  let best = 0, Rpk = 0;
  for (let R = 0.5; R < 15; R += 0.005) { const v = exponentialDiskCircularSpeed(R, M, Rd); if (v > best) { best = v; Rpk = R; } }
  assert.ok(Math.abs(Rpk / Rd - 2.15) < 0.02, `peak at ${Rpk / Rd} R_d`);
  assert.ok(rel(best ** 2 * Rd / (G * M), 0.3870) < 0.01, `v_max² R_d/GM = ${best ** 2 * Rd / (G * M)}`);
  assert.ok(rel(exponentialDiskCircularSpeed(40 * Rd, M, Rd) ** 2, G * M / (40 * Rd)) < 0.01);
  // κ: Kepler κ = Ω, flat curve κ = √2 Ω, solid body κ = 2Ω (BT08 §3.2.3)
  const R = 5;
  assert.ok(rel(epicyclicFrequency(r => Math.sqrt(G * 1e11 / r), R), Math.sqrt(G * 1e11 / R) / R) < 1e-6);
  assert.ok(rel(epicyclicFrequency(() => 200, R), Math.SQRT2 * 200 / R) < 1e-9);
  assert.ok(rel(epicyclicFrequency(r => 30 * r, R), 60) < 1e-9);
  // Q = σκ/(3.36GΣ): Σ = 50 M☉/pc², κ = 37 km/s/kpc, σ = 35 km/s → Q ≈ 1.8 (solar-neighbourhood values, BT08 §6.2.3)
  assert.ok(Math.abs(toomreQ(35, 37, 50e6) - 1.79) < 0.02, `${toomreQ(35, 37, 50e6)}`);
  assert.ok(rel(swingArmNumber(8, 37, 50e6, 2), 37 * 37 * 8 / (2 * Math.PI * G * 50e6 * 2)) < 1e-12);
  assert.throws(() => toomreQ(10, 30, 0), /Sigma/);
  assert.throws(() => exponentialDiskCircularSpeed(-1, M, Rd), /R/);
});

test('PM force of a point mass = exact Plummer-softened force; isolated boundary (no images)', () => {
  const box = 64, N = 64, eps = 1, pm = createPMSolver({ grid: N, box, softening: eps });
  const h = pm.h, mass = 1e10;
  // mass near the left edge, field points across the whole box: a periodic solver would feel the image at +box
  for (const [sx, sy] of [[0.3, -0.2], [-28.6, 3.1]]) {
    pm.deposit(new Float64Array([sx]), new Float64Array([sy]), new Float64Array([mass]), 1); pm.solveForces();
    // CIC PM errors (Hockney & Eastwood 1988, §5–7): a few % near 2–4 cells, falling with distance
    let sq = 0, cnt = 0, worst4 = 0, worst10 = 0;
    for (let r = 2 * h; r < 56; r += 0.31) {for (let q = 0; q < 48; q++) {
      const c = Math.cos(q * 0.1309 + 0.05), s = Math.sin(q * 0.1309 + 0.05), px = sx + r * c, py = sy + r * s;
      if (Math.max(Math.abs(px), Math.abs(py)) > box / 2 - h) continue;
      const f = pm.fieldAt(px, py), exact = G * mass * r / (r * r + eps * eps) ** 1.5;
      const err = Math.hypot(f[0] + exact * c, f[1] + exact * s) / exact;
      if (r < 4 * h) { sq += err * err; cnt++; } else worst4 = Math.max(worst4, err);
      if (r >= 10 * h) worst10 = Math.max(worst10, err);
    }}
    assert.ok(Math.sqrt(sq / cnt) < 0.05, `rms 2–4 cells: ${Math.sqrt(sq / cnt)}`);
    assert.ok(worst4 < 0.07, `≥ 4 cells: ${worst4}`);
    assert.ok(worst10 < 0.015, `≥ 10 cells (a periodic image would give O(1)): ${worst10}`);
  }
  // potential too: Φ = −GM/√(r² + ε²)
  pm.deposit(new Float64Array([0.3]), new Float64Array([-0.2]), new Float64Array([mass]), 1); pm.solvePotential();
  for (const r of [3, 8, 20]) assert.ok(rel(pm.fieldAt(0.3 + r, -0.2)[2], -G * mass / Math.hypot(r, eps)) < 0.01);
  // outside the mesh: monopole
  const far = pm.fieldAt(100, -0.2);
  assert.ok(rel(-far[0], G * mass * 99.7 / (99.7 ** 2 + 1) ** 1.5) < 1e-12);
});

test('PM force of a Kuzmin disk = exact formula (softening ε ≡ Kuzmin seen from height ε)', () => {
  const N = 128, box = 40, M = 1e10, a = 2, pm = createPMSolver({ grid: N, box, softening: box / N });
  const h = pm.h, mesh = new Float64Array(N * N);
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) mesh[i + N * j] = kuzminSurfaceDensity(Math.hypot(pm.x0 + i * h, pm.x0 + j * h), M, a) * h * h;
  pm.loadMesh(mesh); pm.solveForces();
  let worst = 0;
  for (let R = 2 * h; R <= 15; R += 0.37) {for (let q = 0; q < 16; q++) {
    const c = Math.cos(q * 0.39 + 0.1), s = Math.sin(q * 0.39 + 0.1), f = pm.fieldAt(R * c, R * s);
    const ex = kuzminRadialAcceleration(R, M, a, pm.softening);
    worst = Math.max(worst, Math.hypot(f[0] - ex * c, f[1] - ex * s) / Math.abs(ex));
  }}
  assert.ok(worst < 0.03, `Kuzmin force error ${worst}`);
});

test('mesh rotation curve of the exponential disk matches Freeman (1970) for small softening', () => {
  const d = createSelfGravDisk({ n: 100, halo: false, grid: 256, softening: 0.02 });
  const Mu = d.diskMass / (1 - Math.exp(-5) * 6);          // untruncated mass for the same Σ₀
  for (const R of [2.6, 5.2, 7.8, 10.4]) {
    const k = tableAt(d, R);
    assert.ok(rel(d.table.vc[k], exponentialDiskCircularSpeed(d.table.R[k], Mu, 2.6)) < 0.02, `R = ${R}`);
  }
});

test("Newton's third law: mesh self-forces sum to zero, momentum conserved to round-off", () => {
  const rand = rng(7), n = 5000, pm = createPMSolver({ grid: 64, box: 40, softening: 0.5 });
  const x = new Float64Array(n), y = new Float64Array(n), m = new Float64Array(n), ax = new Float64Array(n), ay = new Float64Array(n);
  for (let p = 0; p < n; p++) { x[p] = 3 * rand.normal() + 3; y[p] = 2 * rand.normal(); m[p] = 1e6 * (0.5 + rand()); }
  pm.accelerations(x, y, m, n, ax, ay);
  let fx = 0, fy = 0, scale = 0;
  for (let p = 0; p < n; p++) { fx += m[p] * ax[p]; fy += m[p] * ay[p]; scale += m[p] * Math.hypot(ax[p], ay[p]); }
  assert.equal(pm.outside, 0);
  assert.ok(Math.hypot(fx, fy) / scale < 1e-12, `net force ${Math.hypot(fx, fy) / scale}`);
  // a self-gravitating disk without halo (unstable, but on the mesh): total momentum stays zero
  const d = createSelfGravDisk({ n: 4000, halo: false, grid: 64, Q: 1, box: 60 });
  let pscale = 0; for (let p = 0; p < d.n; p++) pscale += d.m[p] * Math.hypot(d.vx[p], d.vy[p]);
  d.advance(0.002, 60);
  assert.equal(d.outside, 0);
  const [px, py] = d.momentum();
  assert.ok(Math.hypot(px, py) / pscale < 1e-12, `momentum ${Math.hypot(px, py) / pscale}`);
  const [cx, cy] = d.centreOfMass();
  assert.ok(Math.hypot(cx, cy) < 1e-9, `centre of mass drift ${Math.hypot(cx, cy)} kpc`);
});

test('rigid halo + bulge: L_z exact (10 rotations), energy error bounded and O(Δt²)', () => {
  // KDK leapfrog with central forces: kicks exert no torque, drifts keep x × v
  const run = f => {
    const d = createSelfGravDisk({ n: 600, selfGravity: false, bulge: true, Q: 1.5 });
    const L0 = d.angularMomentum(), E0 = d.energy().total, dt = d.suggestedDt * f;
    const T = 2 * Math.PI * 2 * d.Rd / d.table.vc[tableAt(d, 2 * d.Rd)], chunk = Math.ceil(T / dt);
    let worst = 0;
    for (let k = 0; k < 10; k++) { d.advance(dt, chunk); worst = Math.max(worst, rel(d.energy().total, E0)); }
    return { dL: rel(d.angularMomentum(), L0), worst };
  };
  const a = run(1), b = run(0.5);
  assert.ok(a.dL < 1e-12 && b.dL < 1e-12, `ΔL/L ${a.dL}, ${b.dL}`);
  assert.ok(a.worst < 5e-3, `energy ${a.worst}`);
  assert.ok(a.worst / b.worst > 3, `energy error ratio for Δt/2: ${a.worst / b.worst}`);
});

// Equal set-up for the hot and cool disks: only Q differs (Sellwood & Carlberg 2014 use the same contrast).
const SETUP = { n: 20000, grid: 128, diskMass: 3e10, Rd: 2.6, seed: 3 };
const ANNULUS = { Rin: 1.5 * 2.6, Rout: 3.5 * 2.6, Rmax: 13, nR: 20 };
function evolve(Q, { pitch = false } = {}) {
  const d = createSelfGravDisk({ ...SETUP, Q }), dt = 0.002;
  const E0 = d.energy().total, L0 = d.angularMomentum();
  const series = [];
  for (let k = 0; k <= 8; k++) {
    if (k) d.advance(dt, 25);
    const f = fourierModes(d, ANNULUS);
    const entry = { t: d.t, A: Array.from(f.annulus), dominant: f.dominant, noise: f.noise };
    if (pitch) entry.pitch = spiralPitch(d, { m: f.dominant, Rin: ANNULUS.Rin, Rout: ANNULUS.Rout, pMax: 30, dp: 0.5 });
    series.push(entry);
  }
  return { d, series, dE: rel(d.energy().total, E0), dL: rel(d.angularMomentum(), L0) };
}

test('hot disk (Q = 2.5) stays axisymmetric over ~2 rotations at 2 R_d (4 at R_d)', t => {
  const { d, series, dE, dL } = evolve(2.5);
  const noise = series[0].noise;
  const maxA2 = Math.max(...series.map(s => s.A[2])), maxA = Math.max(...series.flatMap(s => s.A.slice(1)));
  t.diagnostic(`hot: noise 1/√N = ${noise.toFixed(4)}, max A2 = ${maxA2.toFixed(4)}, max A_m = ${maxA.toFixed(4)}, ΔE/E = ${dE.toExponential(1)}, ΔL/L = ${dL.toExponential(1)}`);
  assert.ok(maxA2 < 2.5 * noise, `A2 ${maxA2}`);
  assert.ok(maxA < 2.5 * noise, `A_m ${maxA}`);
  assert.ok(dL < 1e-4 && dE < 1e-3, `ΔL ${dL}, ΔE ${dE}`);
  assert.ok(d.t > 0.39);
});

test('cool disk (Q = 1.2) grows trailing multi-armed spirals; arm number ≈ swing-amplification estimate', t => {
  const { d, series, dE, dL } = evolve(1.2, { pitch: true });
  const noise = series[0].noise;
  let peak = series[0];
  for (const s of series) if (s.A[s.dominant] > peak.A[peak.dominant]) peak = s;
  const m = peak.dominant, growth = peak.A[m] / noise, fromStart = peak.A[m] / series[0].A[m];
  // swing estimate at 2 R_d from the initial model (Toomre 1981, X ≈ 2)
  const k = tableAt(d, 2 * d.Rd), mX = swingArmNumber(d.table.R[k], d.table.kappa[k], d.table.Sigma[k], 2);
  t.diagnostic(`cool: A_${m} peaks at ${peak.A[m].toFixed(3)} = ${growth.toFixed(1)}× noise, ${fromStart.toFixed(1)}× its t = 0 value, at t = ${peak.t.toFixed(3)} (${(peak.t * 977.8).toFixed(0)} Myr); m_X = ${mX.toFixed(2)}; pitch ${peak.pitch.pitchDeg.toFixed(1)}° (trailing ${peak.pitch.trailing}); ΔE/E = ${dE.toExponential(1)}, ΔL/L = ${dL.toExponential(1)}`);
  assert.ok(growth > 6 && fromStart > 5, `growth ${growth}, ${fromStart}`);
  assert.ok(m >= 2 && Math.abs(m - mX) <= 1.5, `m = ${m}, swing estimate ${mX}`);
  assert.ok(peak.pitch.trailing && peak.pitch.pitchDeg > 10 && peak.pitch.pitchDeg < 70, `pitch ${peak.pitch.pitchDeg}`);
  assert.ok(dL < 1e-4 && dE < 1e-3, `ΔL ${dL}, ΔE ${dE}`);
  // the spirals heat the disk: Q rises (the reason arms are transient and recur)
  const q = toomreQProfile(d, { nR: 10 });
  const qMid = (q.Q[3] + q.Q[4] + q.Q[5]) / 3;
  assert.ok(qMid > 1.3, `Q after ${qMid}`);
});

test('initial conditions: measured rotation curve and Toomre Q(R) match the model', () => {
  const d = createSelfGravDisk({ n: 40000, grid: 64, Q: 1.3 });
  const rc = rotationCurve(d, { nR: 10 }), q = toomreQProfile(d, { nR: 10 });
  for (let k = 1; k < 9; k++) {
    assert.ok(rel(rc.vc[k], d.table.vc[tableAt(d, rc.R[k])]) < 0.03, `v_c ring ${k}`);
    assert.ok(Math.abs(q.Q[k] - 1.3) < 0.2, `Q ring ${k}: ${q.Q[k]}`);
  }
  const sd = surfaceDensityPolar(d, { nR: 13, nPhi: 8 });
  let ringMass = 0; for (let j = 0; j < 8; j++) ringMass += sd.Sigma[2 + 13 * j];
  const R = sd.R[2], Sigma = ringMass / 8;
  assert.ok(rel(Sigma, d.table.Sigma[tableAt(d, R)]) < 0.05, `Σ(${R}) = ${Sigma}`);
  const rho = localDensity(d, { cells: 64 });
  assert.equal(rho.length, d.n);
  assert.ok(rho.every(v => v >= 0) && rho[0] > 0);
});

test('fourierModes and spiralPitch recover an imposed logarithmic spiral', () => {
  const n = 40000, rand = rng(11), x = new Float64Array(n), y = new Float64Array(n), m = new Float64Array(n).fill(1);
  const arms = 3, pitch = 25 * Math.PI / 180, p = arms / Math.tan(pitch);
  for (let k = 0; k < n; k++) {
    const R = 2 + 8 * rand();
    let ph;                                                       // density ∝ 1 + 0.4 cos(3φ + p ln R)
    do { ph = 2 * Math.PI * rand(); } while (rand() * 1.4 > 1 + 0.4 * Math.cos(arms * ph + p * Math.log(R)));
    x[k] = R * Math.cos(ph); y[k] = R * Math.sin(ph);
  }
  const s = { n, x, y, m };
  const f = fourierModes(s, { Rin: 2, Rout: 10, Rmax: 10, nR: 16 });
  assert.equal(f.dominant, 3);
  assert.ok(Math.abs(f.annulus[3] - 0.2) < 0.02, `A3 = ${f.annulus[3]}`);
  assert.ok(f.annulus[2] < 3 * f.noise, `A2 = ${f.annulus[2]}`);
  // the phase of Σ e^{−3iφ} winds with R as +p ln R
  const k1 = 5, k2 = 13, W = f.mMax + 1;
  let dphase = f.phase[k2 * W + 3] - f.phase[k1 * W + 3];
  dphase -= 2 * Math.PI * Math.round((dphase - p * Math.log(f.R[k2] / f.R[k1])) / (2 * Math.PI));
  assert.ok(Math.abs(dphase - p * Math.log(f.R[k2] / f.R[k1])) < 0.3, `phase winding ${dphase}`);
  const sp = spiralPitch(s, { m: 3, Rin: 2, Rout: 10, dp: 0.1 });
  assert.ok(sp.trailing && Math.abs(sp.pitchDeg - 25) < 1.5, `pitch ${sp.pitchDeg}`);
});

test('reproducible with a seed; independent of how many steps run per frame', () => {
  const opts = { n: 3000, grid: 64, Q: 1.2, seed: 5 };
  const a = createSelfGravDisk(opts), b = createSelfGravDisk(opts), c = createSelfGravDisk({ ...opts, seed: 6 });
  for (let k = 0; k < 12; k++) a.step(0.002);
  for (let k = 0; k < 4; k++) b.advance(0.002, 3);
  assert.deepEqual(a.x, b.x); assert.deepEqual(a.vy, b.vy); assert.equal(a.t, b.t);
  assert.notDeepEqual(a.x.slice(0, 10), c.x.slice(0, 10));
});

test('invalid input throws', () => {
  assert.throws(() => createPMSolver({ grid: 100, box: 10, softening: 1 }), /power of two/);
  assert.throws(() => createPMSolver({ grid: 64, box: 10, softening: 0 }), /softening/);
  assert.throws(() => createSelfGravDisk({ n: 10, Q: -1 }), /Q/);
  assert.throws(() => createSelfGravDisk({ n: 10, Rmax: 20, box: 30 }), /box/);
  const d = createSelfGravDisk({ n: 10, grid: 32 });
  assert.throws(() => d.step(Number.NaN), /dt/);
  assert.throws(() => fourierModes({ n: 1 }), /particle set/);
});

test('performance: 30k particles on a 128² mesh (256² FFT)', t => {
  const d = createSelfGravDisk({ n: 30000, grid: 128 });
  d.advance(0.001, 5);
  const t0 = performance.now(); d.advance(0.001, 40); const rate = 40 / ((performance.now() - t0) / 1000);
  t.diagnostic(`${rate.toFixed(0)} steps/s for 30000 particles, 128² mesh`);
  assert.ok(rate > 5, `${rate} steps/s`);                // sanity only: timing depends on machine load
});
