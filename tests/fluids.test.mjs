import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  poiseuilleProfile, couetteProfile, taylorGreen, cylinderStrouhal, GHIA_RE100, latticeUnits,
  latticeBoltzmann2D, vorticity2D, bilinearSampler, advectParticles,
} from '../core/fluids.mjs';
import { rng } from '../core/numerics.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const walls = { bottom: 'wall', top: 'wall' };

// ---------------------------------------------------------------- analytic references

test('analytic profiles: Poiseuille solves ν u″ = −F with no slip and mean F h²/12ν; Couette is linear', () => {
  const p = { height: 2, force: 0.3, nu: 0.05 }, h = 1e-3;
  for (const y of [0.1, 0.7, 1, 1.9]) {
    const d2 = (poiseuilleProfile(y + h, p) - 2 * poiseuilleProfile(y, p) + poiseuilleProfile(y - h, p)) / (h * h);
    close(p.nu * d2, -p.force, 1e-6);
  }
  close(poiseuilleProfile(0, p), 0, 0); close(poiseuilleProfile(2, p), 0, 0);
  let mean = 0; const n = 2000;
  for (let i = 0; i < n; i++) mean += poiseuilleProfile((i + 0.5) * 2 / n, p) / n;
  close(mean, p.force * 4 / (12 * p.nu), 1e-6);
  close(couetteProfile(0.5, { height: 2, wallSpeed: 3 }), 0.75, 1e-15);
  assert.throws(() => poiseuilleProfile(2.1, p), /between the walls/);
  assert.throws(() => poiseuilleProfile(1, { ...p, nu: 0 }), /nu/);
});

test('taylorGreen satisfies the incompressible Navier–Stokes equations (finite-difference residuals)', () => {
  const o = { nu: 0.07, k: 1.3, U: 0.8, rho: 1.2 }, h = 1e-4;
  const at = (x, y, t) => taylorGreen(x, y, t, o);
  for (const [x, y, t] of [[0.3, 1.1, 0.2], [2.0, -0.4, 1.5], [-1.2, 0.9, 3]]) {
    const c = at(x, y, t), xp = at(x + h, y, t), xm = at(x - h, y, t), yp = at(x, y + h, t), ym = at(x, y - h, t);
    const tp = at(x, y, t + h), tm = at(x, y, t - h);
    close((xp.ux - xm.ux + yp.uy - ym.uy) / (2 * h), 0, 1e-8, 'divergence');
    close((xp.uy - xm.uy - yp.ux + ym.ux) / (2 * h), c.omega, 1e-7, 'vorticity');
    for (const comp of ['ux', 'uy']) {
      const dt = (tp[comp] - tm[comp]) / (2 * h), dx = (xp[comp] - xm[comp]) / (2 * h), dy = (yp[comp] - ym[comp]) / (2 * h);
      const lap = (xp[comp] + xm[comp] + yp[comp] + ym[comp] - 4 * c[comp]) / (h * h);
      const gp = comp === 'ux' ? (xp.p - xm.p) / (2 * h) : (yp.p - ym.p) / (2 * h);
      close(dt + c.ux * dx + c.uy * dy + gp / o.rho - o.nu * lap, 0, 2e-6, `momentum ${comp}`);
    }
  }
});

test('cylinderStrouhal: Williamson (1996) fit gives St(100) ≈ 0.164; Ghia table is consistent', () => {
  close(cylinderStrouhal(100), 0.1643, 1e-4);
  assert.throws(() => cylinderStrouhal(300), /laminar/);
  assert.equal(GHIA_RE100.y.length, GHIA_RE100.u.length);
  assert.equal(GHIA_RE100.x.length, GHIA_RE100.v.length);
  close(GHIA_RE100.u[0], 1, 0); close(Math.min(...GHIA_RE100.u), -0.2109, 0);
});

test('latticeUnits: Re and Ma are preserved, τ = 3ν + ½, conversions are consistent; unstable choices throw', () => {
  // water (ν = 1e-6 m²/s) past a 1 cm object at 1 cm/s: Re = 100
  const u = latticeUnits({ length: 0.01, velocity: 0.01, viscosity: 1e-6, cells: 20, latticeVelocity: 0.05 });
  close(u.Re, 100, 1e-9); close(u.nu, 0.05 * 20 / 100, 1e-15); close(u.tau, 3 * u.nu + 0.5, 1e-15);
  close(u.Ma, 0.05 * Math.sqrt(3), 1e-15);
  close(u.toPhysical.length(20), 0.01, 1e-15); close(u.toPhysical.velocity(0.05), 0.01, 1e-15);
  close(u.toPhysical.time(1) * 1e-6 / u.dx ** 2, u.nu, 1e-15, 'ν_lat = ν dt/dx²');
  close(latticeUnits({ Re: 100, cells: 20 }).tau, u.tau, 1e-15);
  assert.throws(() => latticeUnits({ Re: 1e6, cells: 10, latticeVelocity: 0.01 }), /unstable/);
  assert.throws(() => latticeUnits({ Re: 100, cells: 10, latticeVelocity: 0.2 }), /Ma/);
});

// ---------------------------------------------------------------- lattice Boltzmann

const poiseuilleError = (ny, collision) => {
  const tau = 0.8, nu = (tau - 0.5) / 3, F = 8 * nu * 0.02 / (ny * ny);
  const s = latticeBoltzmann2D({ nx: 2, ny, tau, collision, boundaries: walls, force: [F, 0] });
  s.step(Math.ceil(25 * ny * ny / (Math.PI ** 2 * nu)));       // 25 e-folds of the slowest mode
  let e = 0, n = 0;
  for (let j = 0; j < ny; j++) { const ex = poiseuilleProfile(j + 0.5, { height: ny, force: F, nu }); e += (s.ux[2 * j] - ex) ** 2; n += ex * ex; }
  return Math.sqrt(e / n);
};

test('LBM Poiseuille: BGK converges at 2nd order; TRT with Λ = 3/16 is exact (Ginzburg)', () => {
  const e = [8, 16, 32].map(ny => poiseuilleError(ny, 'bgk'));
  close(Math.log2(e[0] / e[1]), 2, 0.05, 'order 8→16');
  close(Math.log2(e[1] / e[2]), 2, 0.05, 'order 16→32');
  assert.ok(e[2] < 1e-3);
  assert.ok(poiseuilleError(12, 'trt') < 1e-9, 'half-way bounce-back exactly at mid-link');
});

test('LBM Couette: a moving wall (Ladd) gives the exact linear profile', () => {
  const ny = 12, U = 0.05;
  const s = latticeBoltzmann2D({ nx: 1, ny, tau: 1, boundaries: { bottom: 'wall', top: { type: 'wall', u: [U, 0] } } });
  s.step(6000);
  for (let j = 0; j < ny; j++) close(s.ux[j], couetteProfile(j + 0.5, { height: ny, wallSpeed: U }), 1e-10);
});

test('LBM channel: parabolic velocity inflow + pressure outlet develop Poiseuille flow; flux in = flux out; vorticity = −du/dy', () => {
  const nx = 48, ny = 16, nu = 0.1, umax = 0.04, prof = j => 4 * umax * (j + 0.5) * (ny - j - 0.5) / (ny * ny);
  const s = latticeBoltzmann2D({ nx, ny, nu, boundaries: { ...walls, left: { type: 'inflow', u: [0, 0] }, right: { type: 'pressure' } } });
  s.step(200);
  s.setBoundary('left', { type: 'inflow', u: j => [prof(j), 0] });      // impulsive start of the inflow
  s.step(6000);
  const flux = i => { let q = 0; for (let j = 0; j < ny; j++) q += s.rho[i + nx * j] * s.ux[i + nx * j]; return q; };
  close(flux(nx - 1), flux(2), 1e-6 * flux(2), 'steady mass flux');
  // weak compressibility: ρ falls along the pressure gradient, so u = flux/ρ grows by ~1 %; compare the shape
  const mid = nx / 2, scale = s.ux[mid + nx * (ny / 2)] / prof(ny / 2), w = s.vorticity(), F = 8 * nu * umax / (ny * ny);
  close(scale, 1, 0.02, 'amplitude');
  for (let j = 1; j < ny - 1; j++) {
    close(s.ux[mid + nx * j], scale * prof(j), 1e-3 * umax, `u at y = ${j + 0.5}`);
    close(w[mid + nx * j], -scale * F * (ny - 2 * (j + 0.5)) / (2 * nu), 1e-3 * umax, 'ω = −∂u/∂y');
  }
});

test('LBM: mass is conserved to round-off with walls, a moving lid, an obstacle and a body force', () => {
  const nx = 40, ny = 24, rand = rng(3), solid = new Uint8Array(nx * ny);
  for (let j = 8; j < 14; j++) for (let i = 10; i < 15; i++) solid[i + nx * j] = 1;
  for (const collision of ['trt', 'bgk']) {
    const s = latticeBoltzmann2D({
      nx, ny, tau: 0.6, collision, solid, force: [2e-5, -1e-5],
      boundaries: { bottom: 'wall', top: { type: 'wall', u: [0.05, 0] } },
      density: () => 1 + 0.02 * (rand() - 0.5), velocity: () => [0.02 * (rand() - 0.5), 0.02 * (rand() - 0.5)],
    });
    const m0 = s.mass();
    s.step(500);
    close(s.mass(), m0, 1e-12 * m0, collision);
  }
});

test('LBM momentum exchange: steady drag on an obstacle in a periodic box equals the total body force', () => {
  const nx = 32, ny = 24, F = 1e-6, solid = new Uint8Array(nx * ny);
  let fluid = 0;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) if ((i - 15.5) ** 2 + (j - 11.5) ** 2 < 25) solid[i + nx * j] = 2; else fluid++;
  const s = latticeBoltzmann2D({ nx, ny, tau: 1, solid, force: [F, 0] });
  s.step(8000);
  const [fx, fy] = s.force(2);
  close(fx, F * fluid, 1e-6 * F * fluid, 'drag = Σ body force (no walls)');
  close(fy, 0, 1e-9 * F * fluid, 'no lift by symmetry');
  assert.deepEqual(s.force(1), [0, 0]);
});

test('LBM Taylor–Green vortex decays at the rate 2νk²', () => {
  const n = 32, k = 2 * Math.PI / n, U = 0.02, tau = 0.65, nu = (tau - 0.5) / 3;
  const tg = (i, j) => taylorGreen(i, j, 0, { nu, k, U });
  const s = latticeBoltzmann2D({ nx: n, ny: n, tau, velocity: (i, j) => { const v = tg(i, j); return [v.ux, v.uy]; }, density: (i, j) => 1 + 3 * tg(i, j).p });
  const amp = () => { let a = 0, b = 0; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const v = tg(i, j); a += s.ux[i + n * j] * v.ux; b += v.ux * v.ux; } return a / b; };
  s.step(50); const a1 = amp(), t1 = s.t;
  s.step(1000); const a2 = amp(), t2 = s.t;
  const rate = -Math.log(a2 / a1) / (t2 - t1);
  close(rate / (2 * nu * k * k), 1, 5e-3);
});

test('LBM lid-driven cavity at Re = 100 matches Ghia, Ghia & Shin (1982) centrelines', () => {
  const n = 48, U = 0.1;
  const s = latticeBoltzmann2D({ nx: n, ny: n, nu: U * n / 100, boundaries: { left: 'wall', right: 'wall', bottom: 'wall', top: { type: 'wall', u: [U, 0] } } });
  s.step(Math.round(16 * n / U));
  // cell centres at (j + ½)/n; walls at 0 and 1 (u = 0, lid u = 1)
  const uc = j => 0.5 * (s.ux[n / 2 - 1 + n * j] + s.ux[n / 2 + n * j]) / U, vc = i => 0.5 * (s.uy[i + n * (n / 2 - 1)] + s.uy[i + n * (n / 2)]) / U;
  const interp = (f, z, end) => {
    const y = z * n - 0.5;
    if (y < 0) return f(0) * z * n / 0.5;
    if (y > n - 1) return f(n - 1) + (end - f(n - 1)) * (y - n + 1) / 0.5;
    const j = Math.min(Math.floor(y), n - 2), w = y - j;
    return f(j) * (1 - w) + f(j + 1) * w;
  };
  GHIA_RE100.y.forEach((y, k) => close(interp(uc, y, 1), GHIA_RE100.u[k], 0.01, `u(0.5, ${y})`));
  GHIA_RE100.x.forEach((x, k) => close(interp(vc, x, 0), GHIA_RE100.v[k], 0.01, `v(${x}, 0.5)`));
});

test('LBM Kármán vortex street at Re = 100: Strouhal number near Williamson (1996)', () => {
  // coarse (D = 10 cells) and confined (blockage 1/6, far-field velocity on the sides), which raises St
  // above the unconfined 0.164; started from potential flow with a small asymmetric kick
  const D = 10, U = 0.07, nx = 16 * D, ny = 6 * D, cx = 3 * D, cy = ny / 2 + 0.5, R2 = D * D / 4;
  const solid = new Uint8Array(nx * ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) if ((i - cx) ** 2 + (j - cy) ** 2 <= R2) solid[i + nx * j] = 1;
  const potential = (i, j) => {
    const x = i - cx, y = j - cy, r2 = x * x + y * y;
    if (r2 < R2) return [0, 0];
    return [U * (1 - R2 * (x * x - y * y) / (r2 * r2)), -2 * U * R2 * x * y / (r2 * r2) + 0.2 * U * Math.exp(-((x - D) ** 2 + y * y) / R2)];
  };
  const far = { type: 'inflow', u: [U, 0] };
  const s = latticeBoltzmann2D({
    nx, ny, nu: U * D / 100, solid, boundaries: { left: far, right: { type: 'pressure' }, bottom: far, top: far },
    velocity: potential, density: (i, j) => { const [a, b] = potential(i, j); return 1 + 1.5 * (U * U - a * a - b * b); },
  });
  const T = D / U, steps = Math.round(45 * T), window = Math.round(20 * T), lift = new Float64Array(window);
  let drag = 0;
  for (let k = 0; k < steps; k++) {
    s.step();
    if (k >= steps - window) { const [fx, fy] = s.force(1); lift[k - steps + window] = fy; drag += fx; }
  }
  // dominant frequency of the lift by a periodogram scan (robust to acoustic noise)
  let best = 0, St = 0;
  for (let st = 0.05; st < 0.5; st += 5e-4) {
    const om = 2 * Math.PI * st * U / D;
    let c = 0, si = 0;
    for (let k = 0; k < window; k++) { c += lift[k] * Math.cos(om * k); si += lift[k] * Math.sin(om * k); }
    if (c * c + si * si > best) { best = c * c + si * si; St = st; }
  }
  assert.ok(St > 0.16 && St < 0.21, `St = ${St} (Williamson unconfined ${cylinderStrouhal(100).toFixed(3)})`);
  const Cd = 2 * drag / window / (U * U * D), Clmax = 2 * Math.max(...lift.map(Math.abs)) / (U * U * D);
  assert.ok(Cd > 1.2 && Cd < 2.3, `Cd = ${Cd}`);
  assert.ok(Clmax > 0.2, `shedding developed: |Cl| = ${Clmax}`);
  assert.ok(s.maxMach < 0.3);
});

test('LBM guards: τ ≤ ½, Mach limit, inconsistent periodicity, bad masks throw', () => {
  assert.throws(() => latticeBoltzmann2D({ nx: 4, ny: 4, tau: 0.5 }), /τ must exceed/);
  assert.throws(() => latticeBoltzmann2D({ nx: 4, ny: 4, tau: 0.8, nu: 0.1 }), /exactly one/);
  assert.throws(() => latticeBoltzmann2D({ nx: 4, ny: 4, tau: 0.8, boundaries: { left: 'wall' } }), /both sides/);
  assert.throws(() => latticeBoltzmann2D({ nx: 4, ny: 4, tau: 0.8, solid: new Uint8Array(3) }), /solid/);
  assert.throws(() => latticeBoltzmann2D({ nx: 4, ny: 4, tau: 0.8, boundaries: { ...walls, left: 'inflow', right: 'outflow' } }), /inflow needs u/);
  const s = latticeBoltzmann2D({ nx: 4, ny: 4, tau: 0.8, velocity: () => [0.25, 0] });
  assert.throws(() => s.step(), /Mach number/);
});

// ---------------------------------------------------------------- pseudo-spectral vorticity

const grid2 = (n, f) => { const a = new Float64Array(n * n), h = 2 * Math.PI / n; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) a[i + n * j] = f(i * h, j * h); return a; };
const randomVorticity = (n, seed) => {
  const r = rng(seed), modes = [];
  for (let m = 0; m < 12; m++) modes.push([1 + Math.floor(r() * 5), Math.floor(r() * 9) - 4, r() * 6.3, r() - 0.5]);
  return grid2(n, (x, y) => modes.reduce((s, [a, b, ph, amp]) => s + amp * Math.cos(a * x + b * y + ph), 0));
};

test('vorticity2D: Taylor–Green decays exactly as e^{−2νk²t}; velocity matches the exact field', () => {
  const n = 32, nu = 0.05, k = 2;
  const s = vorticity2D({ omega0: grid2(n, (x, y) => taylorGreen(x, y, 0, { nu, k }).omega), shape: [n, n], nu, dt: 0.05 });
  s.step(40);
  const w = s.vorticity(), { ux, uy } = s.velocity(), h = 2 * Math.PI / n;
  let e = 0;
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < n; i++) {
      const ex = taylorGreen(i * h, j * h, s.t, { nu, k });
      e = Math.max(e, Math.abs(w[i + n * j] - ex.omega), Math.abs(ux[i + n * j] - ex.ux), Math.abs(uy[i + n * j] - ex.uy));
    }
  }
  assert.ok(e < 1e-12, `max error ${e}`);
  close(s.energy(), 0.25 * Math.exp(-4 * nu * k * k * s.t), 1e-13, 'energy U²/4·e^{−4νk²t}');
});

test('vorticity2D: inviscid energy and enstrophy are conserved; the time error is 4th order', () => {
  const n = 32, omega0 = randomVorticity(n, 5);
  const s = vorticity2D({ omega0, shape: [n, n], dt: 0.01 });
  const E0 = s.energy(), Z0 = s.enstrophy();
  s.step(200);
  close(s.energy(), E0, 1e-7 * E0, 'energy'); close(s.enstrophy(), Z0, 1e-7 * Z0, 'enstrophy');
  const run = dt => { const r = vorticity2D({ omega0, shape: [n, n], nu: 0.01, dt }); r.step(Math.round(1 / dt)); return r.vorticity(); };
  const a = run(0.04), b = run(0.02), c = run(0.01);
  const d1 = Math.max(...a.map((v, i) => Math.abs(v - b[i]))), d2 = Math.max(...b.map((v, i) => Math.abs(v - c[i])));
  close(Math.log2(d1 / d2), 4, 0.25, 'order in dt');
});

test('vorticity2D: forced Kolmogorov flow with drag relaxes to ω = f/(νn² + μ)', () => {
  const n = 16, nu = 0.1, mu = 0.05, m = 2, F = 0.3;
  const s = vorticity2D({ omega0: new Float64Array(n * n), forcing: grid2(n, (_x, y) => F * Math.cos(m * y)), shape: [n, n], nu, drag: mu, dt: 0.1 });
  s.step(1000);
  const w = s.vorticity(), h = 2 * Math.PI / n;
  for (let j = 0; j < n; j++) close(w[3 + n * j], F * Math.cos(m * j * h) / (nu * m * m + mu), 1e-7);   // Lawson RK4 fixed point: O((λdt)⁴) relative
});

test('vorticity2D: streamfunction satisfies ∇²ψ = −ω; CFL and shape guards throw', () => {
  const n = 16, omega0 = randomVorticity(n, 9);
  const s = vorticity2D({ omega0, shape: [n, n], dt: 0.01, dealias: false });
  const psi = s.streamfunction(), w = s.vorticity();
  // check a single Fourier mode: ω = cos(x + 2y) ⇒ ψ = ω/5
  const one = vorticity2D({ omega0: grid2(n, (x, y) => Math.cos(x + 2 * y)), shape: [n, n], dt: 0.01 });
  const p1 = one.streamfunction(), w1 = one.vorticity();
  p1.forEach((v, i) => close(v, w1[i] / 5, 1e-14));
  assert.equal(psi.length, w.length);
  assert.throws(() => vorticity2D({ omega0, shape: [12, 16], dt: 0.01 }), /powers of two/);
  const fast = vorticity2D({ omega0: Float64Array.from(omega0, v => 100 * v), shape: [n, n], dt: 0.5 });
  assert.throws(() => fast.step(), /CFL/);
});

// ---------------------------------------------------------------- tracers

test('bilinearSampler is exact for linear fields (clamped and periodic); advectParticles is 4th-order RK4', () => {
  const nx = 7, ny = 5, ux = new Float64Array(nx * ny), uy = new Float64Array(nx * ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const x = 0.5 + i * 0.3, y = -1 + j * 0.3; ux[i + nx * j] = 2 * x - y + 1; uy[i + nx * j] = 0.5 * x + 3 * y; }
  const v = bilinearSampler({ ux, uy, nx, ny, spacing: 0.3, origin: [0.5, -1] });
  const out = new Float64Array(2);
  for (const [x, y] of [[0.61, -0.83], [1.9, -0.05], [2.2, 0.2]]) {
    v(x, y, 0, out); close(out[0], 2 * x - y + 1, 1e-13); close(out[1], 0.5 * x + 3 * y, 1e-13);
  }
  v(-5, 9, 0, out); close(out[0], 2 * 0.5 - 0.2 + 1, 1e-13, 'clamped to the last node');
  const per = bilinearSampler({ ux: Float64Array.from({ length: 16 }, (_, i) => i % 4), uy: new Float64Array(16), nx: 4, ny: 4, periodic: [true, true] });
  close(per(3.5, 0, 0, out)[0], 1.5, 1e-15, 'wraps between the last and first node');
  close(per(-0.5, 7, 0, out)[0], 1.5, 1e-15);
  // solid-body rotation u = (−y, x): exact orbit is a circle
  const rot = (x, y, _t, o) => { o[0] = -y; o[1] = x; return o; };
  const err = dt => { const p = Float64Array.of(1, 0, 0, 2); advectParticles(p, rot, dt, { steps: Math.round(Math.PI / dt) }); return Math.hypot(p[0] + 1, p[1], p[2], p[3] + 2); };
  close(Math.log2(err(0.1 * Math.PI) / err(0.05 * Math.PI)), 4, 0.1, 'RK4 order');
  assert.throws(() => advectParticles(new Float64Array(3), rot, 0.1), /even length/);
});

test('tracers in a Taylor–Green vortex follow streamlines (ψ conserved along steady paths)', () => {
  const n = 64, nu = 0, h = 2 * Math.PI / n;
  const s = vorticity2D({ omega0: grid2(n, (x, y) => taylorGreen(x, y, 0, { nu }).omega), shape: [n, n], dt: 0.05 });
  const { ux, uy } = s.velocity();
  const v = bilinearSampler({ ux, uy, nx: n, ny: n, spacing: h, periodic: [true, true] });
  const psi = (x, y) => Math.cos(x) * Math.cos(y);                // u = ∂ψ/∂y, v = −∂ψ/∂x up to sign convention
  const p = Float64Array.of(1.0, 0.7, 2.1, 1.4, 4.0, 4.5);
  const before = [psi(p[0], p[1]), psi(p[2], p[3]), psi(p[4], p[5])];
  advectParticles(p, v, 0.02, { steps: 300 });
  [0, 1, 2].forEach(k => close(psi(p[2 * k], p[2 * k + 1]), before[k], 2e-3, `tracer ${k}`));
});
