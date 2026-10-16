import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  finiteVolume, linearAdvection, inviscidBurgers, euler, shallowWater,
  burgersExactRiemann, eulerExactRiemann, shallowWaterExactRiemann,
} from '../core/conservation.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const rel = (a, b, tol, msg) => close(a, b, tol * Math.abs(b), msg);
const conserved = (fv, before, tol = 1e-13) => fv.totals().forEach((v, k) => close(v - before[k] - fv.inflow[k], 0, tol * Math.max(1, Math.abs(before[k]), Math.abs(fv.inflow[k])), `variable ${k}`));
/** Agreement with a printed table value: half a unit in its last digit (or 1e−5 relative, inputs being rounded too). */
const printed = (a, text, msg) => close(a, Number(text), Math.max(0.5 * 10 ** -((text.split('.')[1] ?? '').length), 1e-5 * Math.abs(Number(text))), msg);

/** 1D Euler Riemann problem on [0, 1], diaphragm at x0; returns the solver and the L1 density error against the exact solution. */
function eulerShockTube(n, L, R, T, opts = {}, x0 = 0.5, gamma = 1.4) {
  const dx = 1 / n, q0 = [new Float64Array(n), new Float64Array(n), new Float64Array(n)];
  for (let i = 0; i < n; i++) { const s = (i + 0.5) * dx < x0 ? L : R; q0[0][i] = s.rho; q0[1][i] = s.rho * s.u; q0[2][i] = s.p / (gamma - 1) + 0.5 * s.rho * s.u * s.u; }
  const fv = finiteVolume({ system: euler({ gamma }), q0, shape: [n], spacing: dx, ...opts });
  const before = fv.totals();
  fv.advanceTo(T);
  const ex = eulerExactRiemann(L, R, gamma), rho = fv.primitive(0), p = fv.primitive(2);
  let e = 0; for (let i = 0; i < n; i++) e += Math.abs(rho[i] - ex.sample(((i + 0.5) * dx - x0) / fv.t).rho) * dx;
  return { fv, e, before, minRho: Math.min(...rho), minP: Math.min(...p) };
}

test('exact Euler Riemann solver reproduces Toro (2009) Table 4.3 star states', () => {
  const cases = [
    [{ rho: 1, u: 0, p: 1 }, { rho: 0.125, u: 0, p: 0.1 }, ['0.30313', '0.92745', '0.42632', '0.26557']],             // Sod
    [{ rho: 1, u: -2, p: 0.4 }, { rho: 1, u: 2, p: 0.4 }, ['0.00189', '0.00000', '0.02185', '0.02185']],              // 123 problem
    [{ rho: 1, u: 0, p: 1000 }, { rho: 1, u: 0, p: 0.01 }, ['460.894', '19.5975', '0.57506', '5.99924']],
    [{ rho: 1, u: 0, p: 0.01 }, { rho: 1, u: 0, p: 100 }, ['46.0950', '-6.19633', '5.99242', '0.57511']],
    [{ rho: 5.99924, u: 19.5975, p: 460.894 }, { rho: 5.99242, u: -6.19633, p: 46.0950 }, ['1691.64', '8.68975', '14.2823', '31.0426']],
  ];
  for (const [L, R, [p, u, rl, rr]] of cases) {
    const s = eulerExactRiemann(L, R);
    printed(s.pStar, p, 'p*'); printed(s.uStar, u, 'u*'); printed(s.rhoStarL, rl, 'ρ*L'); printed(s.rhoStarR, rr, 'ρ*R');
    const mid = s.sample(s.uStar);
    rel(mid.p, s.pStar, 1e-14, 'contact carries p*');
  }
  // Rankine–Hugoniot across the Sod shock: mass and momentum fluxes continuous in the shock frame
  const sod = eulerExactRiemann({ rho: 1, u: 0, p: 1 }, { rho: 0.125, u: 0, p: 0.1 }), g = 1.4;
  const S = 0 + Math.sqrt(1.4 * 0.1 / 0.125) * Math.sqrt((g + 1) / (2 * g) * sod.pStar / 0.1 + (g - 1) / (2 * g));
  const a = sod.sample(S - 1e-9), b = sod.sample(S + 1e-9);
  close(a.rho * (a.u - S), b.rho * (b.u - S), 1e-12); close(a.rho * (a.u - S) ** 2 + a.p, b.rho * (b.u - S) ** 2 + b.p, 1e-12);
  assert.throws(() => eulerExactRiemann({ rho: 1, u: -20, p: 1 }, { rho: 1, u: 20, p: 1 }), /vacuum/);
});

test('Sod shock tube: every Riemann solver converges to the exact solution and conserves to round-off', () => {
  const L = { rho: 1, u: 0, p: 1 }, R = { rho: 0.125, u: 0, p: 0.1 };
  for (const riemann of ['hllc', 'hll', 'roe', 'exact', 'rusanov']) {
    const { fv, e, before } = eulerShockTube(200, L, R, 0.2, { riemann });
    assert.ok(e < 0.004, `${riemann}: L1(ρ) = ${e}`);
    conserved(fv, before);
  }
  const errs = [100, 200, 400].map(n => eulerShockTube(n, L, R, 0.2).e);
  const order = Math.log2(errs[0] / errs[2]) / 2;
  assert.ok(order > 0.75 && order < 1.2, `L1 order ${order} (discontinuous solution: ≈ 1)`);
  // a first-order scheme is measurably worse than MUSCL
  assert.ok(eulerShockTube(200, L, R, 0.2, { limiter: 'first' }).e > 1.5 * errs[1]);
});

test('Lax, Toro 123 (near-vacuum) and strong-shock (Toro test 3) problems: accurate and positive', () => {
  const lax = eulerShockTube(200, { rho: 0.445, u: 0.698, p: 3.528 }, { rho: 0.5, u: 0, p: 0.571 }, 0.13);
  assert.ok(lax.e < 0.015, `Lax L1 ${lax.e}`); conserved(lax.fv, lax.before);
  const t123 = eulerShockTube(200, { rho: 1, u: -2, p: 0.4 }, { rho: 1, u: 2, p: 0.4 }, 0.15);
  assert.ok(t123.e < 0.005 && t123.minRho > 0 && t123.minP > 0, `123 L1 ${t123.e}`);
  const strong = eulerShockTube(400, { rho: 1, u: 0, p: 1000 }, { rho: 1, u: 0, p: 0.01 }, 0.012, { limiter: 'minmod' });
  assert.ok(strong.e < 0.08 && strong.minP > 0, `test 3 L1 ${strong.e}`); conserved(strong.fv, strong.before);
});

test('Riemann solvers are consistent: F(q, q) equals the physical flux, in both directions in 2D', () => {
  const systems = [[euler({ dims: 2 }), [1.3, 0.4, -0.7, 2.1]], [shallowWater({ dims: 2 }), [1.7, -0.3, 0.8]], [euler(), [0.6, -1.2, 0.9]]];
  for (const [sys, p] of systems) {
    const want = new Float64Array(sys.nvar), got = new Float64Array(sys.nvar);
    for (let ax = 0; ax < sys.dims; ax++) {
      sys.flux(p, ax, want);
      for (const [name, solve] of Object.entries(sys.solvers)) { solve(p, p, ax, got); want.forEach((w, v) => close(got[v], w, 1e-13 * (1 + Math.abs(w)), `${sys.name} ${name} axis ${ax}`)); }
    }
  }
});

test('linear advection: measured L1 order 2 (unlimited, van Leer), ≈ 1 for Godunov; periodic total to round-off', () => {
  const order = (limiter, integrator = 'ssprk3') => {
    const errs = [64, 128, 256].map(n => {
      const dx = 1 / n, avg = (i, t) => (Math.cos(2 * Math.PI * (i * dx - t)) - Math.cos(2 * Math.PI * ((i + 1) * dx - t))) / (2 * Math.PI * dx);
      const fv = finiteVolume({ system: linearAdvection(1), q0: [Float64Array.from({ length: n }, (_, i) => avg(i, 0))], shape: [n], spacing: dx, bc: 'periodic', limiter, integrator });
      fv.advanceTo(1);
      close(fv.totals()[0], 0, 1e-14);
      let e = 0; for (let i = 0; i < n; i++) e += Math.abs(fv.q[0][i] - avg(i, fv.t)) * dx;
      return e;
    });
    return Math.log2(errs[1] / errs[2]);
  };
  close(order('none'), 2, 0.05, 'unlimited MUSCL'); close(order('none', 'ssprk2'), 2, 0.1, 'SSP-RK2');
  close(order('vanleer'), 2, 0.1, 'van Leer'); assert.ok(order('mc') > 1.85); assert.ok(order('minmod') > 1.7);
  close(order('first'), 1, 0.1, 'Godunov');
});

test('Euler smooth density wave: exact translation, second-order convergence (MC limiter)', () => {
  const errs = [64, 128, 256].map(n => {
    const dx = 1 / n, avg = (i, t) => 1 + 0.2 * (Math.cos(2 * Math.PI * (i * dx - t)) - Math.cos(2 * Math.PI * ((i + 1) * dx - t))) / (2 * Math.PI * dx);
    const rho = Float64Array.from({ length: n }, (_, i) => avg(i, 0));
    const fv = finiteVolume({ system: euler(), q0: [rho, Float64Array.from(rho), Float64Array.from(rho, r => 1 / 0.4 + 0.5 * r)], shape: [n], spacing: dx, bc: 'periodic' });
    const before = fv.totals();
    fv.advanceTo(1);
    conserved(fv, before);
    let e = 0; for (let i = 0; i < n; i++) e += Math.abs(fv.q[0][i] - avg(i, fv.t)) * dx;
    return e;
  });
  assert.ok(Math.log2(errs[1] / errs[2]) > 1.8, `order ${Math.log2(errs[1] / errs[2])}`);
});

test('inviscid Burgers: entropy shock at speed (uL + uR)/2, rarefaction fan, transonic fan', () => {
  for (const [uL, uR] of [[1, 0], [0, 1], [-1, 1], [2, -1]]) {
    const n = 200, dx = 2 / n, x = i => -1 + (i + 0.5) * dx;
    const fv = finiteVolume({ system: inviscidBurgers(), q0: [Float64Array.from({ length: n }, (_, i) => (x(i) < 0 ? uL : uR))], shape: [n], spacing: dx });
    const before = fv.totals();
    fv.advanceTo(0.4);
    const ex = burgersExactRiemann(uL, uR);
    let e = 0; for (let i = 0; i < n; i++) e += Math.abs(fv.q[0][i] - ex.sample(x(i) / fv.t)) * dx;
    assert.ok(e < 0.01 * Math.abs(uL - uR) * 2, `${uL}|${uR}: L1 ${e}`);
    conserved(fv, before);
  }
  // periodic sine steepens into a shock; total conserved to round-off, maximum never exceeded (TVD)
  const n = 128, dx = 2 * Math.PI / n, fv = finiteVolume({ system: inviscidBurgers(), q0: [Float64Array.from({ length: n }, (_, i) => 0.5 + Math.sin((i + 0.5) * dx))], shape: [n], spacing: dx, bc: 'periodic', limiter: 'minmod' });
  const T0 = fv.totals()[0];
  fv.advanceTo(3);
  close(fv.totals()[0], T0, 1e-13 * T0); assert.ok(Math.max(...fv.q[0]) <= 1.5 + 1e-12 && Math.min(...fv.q[0]) >= -0.5 - 1e-12);
});

test('exact shallow-water Riemann solver: Ritter formula, Rankine–Hugoniot bore, Riemann invariant', () => {
  const g = 9.81, c0 = Math.sqrt(g);
  const ritter = shallowWaterExactRiemann({ h: 1, u: 0 }, { h: 0, u: 0 }, g);
  for (const xi of [-c0 + 1e-9, -1, 0, 1.5, 2 * c0 - 1e-3]) { const s = ritter.sample(xi); close(s.h, (2 * c0 - xi) ** 2 / (9 * g), 1e-14); close(s.u, 2 * (c0 + xi) / 3, 1e-14); }
  assert.equal(ritter.sample(2 * c0 + 1e-9).h, 0); assert.equal(ritter.sample(-c0 - 1e-9).h, 1);
  // Stoker: hL = 1, hR = 0.1 → rarefaction + bore
  const st = shallowWaterExactRiemann({ h: 1, u: 0 }, { h: 0.1, u: 0 }, g), hs = st.hStar, us = st.uStar;
  close(us + 2 * Math.sqrt(g * hs), 2 * c0, 1e-12, 'u + 2√(gh) constant across the left rarefaction');
  const S = hs * us / (hs - 0.1);                                                 // bore speed from mass conservation
  close(S * hs * us, (hs * us * us + 0.5 * g * hs * hs) - 0.5 * g * 0.01, 1e-12, 'momentum jump condition (still water ahead)');
  close(st.sample(S - 1e-9).h, hs, 1e-14); close(st.sample(S + 1e-9).h, 0.1, 0);
  // a dry region opens between two receding streams
  const apart = shallowWaterExactRiemann({ h: 1, u: -8 }, { h: 1, u: 8 }, g);
  assert.equal(apart.uStar, null); assert.equal(apart.sample(0).h, 0);
});

test('dam break: Ritter (dry bed) and Stoker (wet bed) converge to the exact solution; h ≥ 0; mass to round-off', () => {
  const g = 9.81;
  const dam = (n, hR, opts = {}) => {
    const dx = 2 / n, x = i => -1 + (i + 0.5) * dx;
    const fv = finiteVolume({ system: shallowWater({ g }), q0: [Float64Array.from({ length: n }, (_, i) => (x(i) < 0 ? 1 : hR)), new Float64Array(n)], shape: [n], spacing: dx, ...opts });
    const m0 = fv.totals()[0];
    fv.advanceTo(0.15);
    close(fv.totals()[0] - m0 - fv.inflow[0], 0, 1e-14);
    assert.ok(Math.min(...fv.q[0]) >= 0);
    const ex = shallowWaterExactRiemann({ h: 1, u: 0 }, { h: hR, u: 0 }, g);
    let e = 0; for (let i = 0; i < n; i++) e += Math.abs(fv.q[0][i] - ex.sample(x(i) / fv.t).h) * dx;
    return e;
  };
  for (const hR of [0, 0.1]) {
    const e = [100, 200, 400].map(n => dam(n, hR));
    assert.ok(e[2] < 0.002 && Math.log2(e[0] / e[2]) / 2 > 0.85, `hR = ${hR}: ${e}`);
  }
  for (const riemann of ['hll', 'rusanov', 'exact']) assert.ok(dam(200, 0.1, { riemann }) < 0.005, riemann);
});

test('shallow water: lake at rest over a bump and a dry island stays at rest (well balanced); waves travel at √(gh)', () => {
  const g = 9.81, n = 100, dx = 1 / n, x = i => (i + 0.5) * dx;
  for (const top of [0.8, 1.5]) {                                              // submerged bump; emerging island (dry cells)
    const b = Float64Array.from({ length: n }, (_, i) => top * Math.exp(-50 * (x(i) - 0.5) ** 2));
    const fv = finiteVolume({ system: shallowWater({ g, bathymetry: b }), q0: [Float64Array.from(b, v => Math.max(0, 1 - v)), new Float64Array(n)], shape: [n], spacing: dx, bc: 'reflective' });
    const m0 = fv.totals()[0];
    fv.advanceTo(0.5);
    assert.ok(Math.max(...fv.q[1].map(Math.abs)) < 1e-13, `momentum stays zero (top ${top})`);
    for (let i = 0; i < n; i++) if (fv.q[0][i] > 0) close(fv.q[0][i] + b[i], 1, 1e-13, 'flat surface');
    close(fv.totals()[0], m0, 1e-14);
  }
  // small Gaussian hump on depth h0 splits into two pulses moving at ±√(g h0)
  const m = 800, h = 40 / m, y = i => -20 + (i + 0.5) * h, h0 = 2, amp = 1e-4;
  const fv = finiteVolume({ system: shallowWater({ g }), q0: [Float64Array.from({ length: m }, (_, i) => h0 + amp * Math.exp(-(y(i) ** 2))), new Float64Array(m)], shape: [m], spacing: h, bc: 'periodic' });
  fv.advanceTo(2);
  let mass = 0, moment = 0; for (let i = 0; i < m; i++) if (y(i) > 0) { mass += fv.q[0][i] - h0; moment += (fv.q[0][i] - h0) * y(i); }
  rel(moment / mass / fv.t, Math.sqrt(g * h0), 1e-3, 'pulse speed');
});

test('2D Euler blast in a closed box: x↔y and mirror symmetry, mass/momentum/energy to round-off, positivity', () => {
  const n = 40, dx = 1 / n, size = n * n, q0 = [0, 1, 2, 3].map(() => new Float64Array(size));
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const r = Math.hypot((i + 0.5) * dx - 0.5, (j + 0.5) * dx - 0.5); q0[0][i + n * j] = 1; q0[3][i + n * j] = (r < 0.1 ? 10 : 0.1) / 0.4; }
  const fv = finiteVolume({ system: euler({ dims: 2 }), q0, shape: [n, n], spacing: dx, bc: 'reflective' });
  const T0 = fv.totals();
  fv.advanceTo(0.15);
  const T1 = fv.totals();
  close(T1[0], T0[0], 1e-14 * T0[0]); close(T1[1], 0, 1e-14); close(T1[2], 0, 1e-14); close(T1[3], T0[3], 1e-14 * T0[3]);
  let transpose = 0, mirror = 0;
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < n; i++) {
      transpose = Math.max(transpose, Math.abs(fv.q[0][i + n * j] - fv.q[0][j + n * i]));
      mirror = Math.max(mirror, Math.abs(fv.q[0][i + n * j] - fv.q[0][n - 1 - i + n * j]));
    }
  }
  assert.ok(transpose < 1e-13 && mirror < 1e-12, `symmetry ${transpose} ${mirror}`);
  assert.ok(Math.min(...fv.primitive(3)) > 0 && Math.min(...fv.q[0]) > 0);
});

test('2D Riemann quadrant (Lax–Liu configuration 3) and circular dam break: symmetric and conservative', () => {
  const n = 32, dx = 1 / n, size = n * n, q0 = [0, 1, 2, 3].map(() => new Float64Array(size));
  const st = (x, y) => (x >= 0.5 ? (y >= 0.5 ? [1.5, 0, 0, 1.5] : [0.5323, 0, 1.206, 0.3]) : (y >= 0.5 ? [0.5323, 1.206, 0, 0.3] : [0.138, 1.206, 1.206, 0.029]));
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const [r, u, v, p] = st((i + 0.5) * dx, (j + 0.5) * dx), k = i + n * j; q0[0][k] = r; q0[1][k] = r * u; q0[2][k] = r * v; q0[3][k] = p / 0.4 + 0.5 * r * (u * u + v * v); }
  const fv = finiteVolume({ system: euler({ dims: 2 }), q0, shape: [n, n], spacing: dx, bc: 'transmissive' });
  const before = fv.totals();
  fv.advanceTo(0.2);
  conserved(fv, before);
  let sym = 0; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) sym = Math.max(sym, Math.abs(fv.q[0][i + n * j] - fv.q[0][j + n * i]), Math.abs(fv.q[1][i + n * j] - fv.q[2][j + n * i]));
  assert.ok(sym < 1e-13, `diagonal symmetry ${sym}`);
  const m = 32, h = 1 / m, w0 = [new Float64Array(m * m), new Float64Array(m * m), new Float64Array(m * m)];
  for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) w0[0][i + m * j] = Math.hypot((i + 0.5) * h - 0.5, (j + 0.5) * h - 0.5) < 0.2 ? 2 : 1;
  const sw = finiteVolume({ system: shallowWater({ dims: 2 }), q0: w0, shape: [m, m], spacing: h, bc: 'reflective' });
  const M0 = sw.totals()[0];
  sw.advanceTo(0.08);
  close(sw.totals()[0], M0, 1e-14 * M0);
  let s2 = 0; for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) s2 = Math.max(s2, Math.abs(sw.q[0][i + m * j] - sw.q[0][j + m * i]));
  assert.ok(s2 < 1e-13);
  // well-balanced in 2D: lake at rest over a Gaussian seamount
  const b = Float64Array.from({ length: m * m }, (_, k) => 0.5 * Math.exp(-40 * (((k % m) + 0.5) * h - 0.5) ** 2 - 40 * ((Math.floor(k / m) + 0.5) * h - 0.5) ** 2));
  const lake = finiteVolume({ system: shallowWater({ dims: 2, bathymetry: b }), q0: [Float64Array.from(b, v => 1 - v), new Float64Array(m * m), new Float64Array(m * m)], shape: [m, m], spacing: h, bc: 'reflective' });
  lake.advanceTo(0.1);
  assert.ok(Math.max(...lake.q[1].map(Math.abs), ...lake.q[2].map(Math.abs)) < 1e-13);
});

test('input validation: CFL violation, bad options and inadmissible states throw', () => {
  const q0 = [new Float64Array(16).fill(1), new Float64Array(16), new Float64Array(16).fill(2.5)];
  const fv = finiteVolume({ system: euler(), q0, shape: [16], spacing: 1 / 16 });
  assert.throws(() => fv.step(10 * fv.maxDt()), /Courant/);
  assert.throws(() => finiteVolume({ system: euler(), q0, shape: [16], limiter: 'nope' }), /limiter/);
  assert.throws(() => finiteVolume({ system: euler(), q0, shape: [16], riemann: 'nope' }), /riemann/);
  assert.throws(() => finiteVolume({ system: euler(), q0, shape: [16], cfl: 1.5 }), /cfl/);
  assert.throws(() => finiteVolume({ system: euler({ dims: 2 }), q0, shape: [16] }), /2D/);
  assert.throws(() => finiteVolume({ system: euler(), q0: [q0[0], q0[1], new Float64Array(16).fill(-1)], shape: [16] }), /non-admissible/);
  assert.throws(() => finiteVolume({ system: euler(), q0, shape: [16], bc: ['periodic', 'reflective'] }), /periodic/);
});
