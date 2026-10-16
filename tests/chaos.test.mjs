import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  logisticMap, tentMap, henonMap, standardMap, circleMap, iterateMap, bifurcationDiagram,
  superstableParameters, feigenbaumDelta, rotationNumber, mapLyapunov,
  lorenz, rossler, duffing, henonHeiles, doublePendulum, lyapunovSpectrum, poincareSection,
  correlationDimension, boxCountingDimension, smallEigenvalues, classifyFixedPoint, findFixedPoint, findPeriodicOrbit,
} from '../core/chaos.mjs';
import { odeStepper } from '../core/ode.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const pointsOf = (flat, d) => Array.from({ length: flat.length / d }, (_, i) => flat.subarray(d * i, d * i + d));

// ---------------------------------------------------------------- maps

test('maps: definitions, analytic Jacobians against finite differences, determinants', () => {
  const cases = [[logisticMap(3.7), [0.31]], [tentMap(1.9), [0.3]], [henonMap(), [0.4, -0.2]],
    [standardMap(1.3, { wrap: false }), [0.7, 0.4]], [circleMap(0.3, 0.8), [0.21]]];
  for (const [m, x] of cases) {
    const J = m.jacobian(x), e = 1e-6;
    for (let j = 0; j < m.dim; j++) {
      const p = Float64Array.from(x), q = Float64Array.from(x); p[j] += e; q[j] -= e;
      const fp = m.map(p), fq = m.map(q);
      for (let i = 0; i < m.dim; i++) close(J[i][j], (fp[i] - fq[i]) / (2 * e), 1e-8, `${m.name} J[${i}][${j}]`);
    }
  }
  const Jh = henonMap().jacobian([0.4, 0]);
  close(Jh[0][0] * Jh[1][1] - Jh[0][1] * Jh[1][0], -0.3, 1e-16, 'Hénon det J = −b');
  const Js = standardMap(2).jacobian([1.1, 0]);
  close(Js[0][0] * Js[1][1] - Js[0][1] * Js[1][0], 1, 1e-15, 'standard map area-preserving');
  const w = standardMap(5).map([6, 3]);
  assert.ok(w[0] >= 0 && w[0] < 2 * Math.PI && w[1] >= -Math.PI && w[1] < Math.PI, 'wrapped to the torus');
  assert.equal(iterateMap(logisticMap(2.5), [0.1], 5, { transient: 200 })[4], 0.6, 'stable fixed point 1 − 1/r');
  assert.throws(() => iterateMap(logisticMap(5), [0.5], 100), /diverged/);
});

test('bifurcation diagram: period 1, 2, 4 and chaos of the logistic map', () => {
  const rows = bifurcationDiagram({ family: logisticMap, params: [2.8, 3.2, 3.5, 3.9], x0: [0.5], transient: 2000, samples: 64 });
  const distinct = row => new Set(Array.from(row.values, v => v.toFixed(6))).size;
  assert.deepEqual(rows.map(distinct).slice(0, 3), [1, 2, 4]);
  assert.ok(distinct(rows[3]) > 50);
  // exact period-2 orbit x = (r + 1 ± √((r − 3)(r + 1)))/(2r)
  const r = 3.2, x2 = (r + 1 + Math.sqrt((r - 3) * (r + 1))) / (2 * r);
  close(Math.max(...rows[1].values), x2, 1e-12);
});

test('Feigenbaum: superstable parameters of the logistic map and δ → 4.669201609…', () => {
  const R = superstableParameters({ g: (r, x) => r * x * (1 - x), xc: 0.5, r0: 2, r1: 1 + Math.sqrt(5), count: 12 });
  close(R[2], 3.4985616993277, 1e-12, 'R_2: superstable 4-cycle');
  close(R[11], 3.5699456718, 2e-7, 'converging to r∞ = 3.5699456718 (Feigenbaum 1978)');
  const d = feigenbaumDelta(R);
  close(d.at(-1), 4.669201609102990, 2e-6, 'δ');
  assert.ok(Math.abs(d.at(-1) - 4.6692016091) < Math.abs(d[3] - 4.6692016091) / 100, 'converges');
  assert.throws(() => feigenbaumDelta([1, 2]), /three/);
});

test('Lyapunov exponents of maps: logistic r = 4 → ln 2, tent → ln μ exactly, Hénon spectrum, area preservation', () => {
  close(mapLyapunov(logisticMap(4), [0.3], { n: 400000 }).exponents[0], Math.LN2, 2e-3, 'logistic r = 4 (Ulam–von Neumann)');
  close(mapLyapunov(tentMap(1.9), [0.3], { n: 1000 }).exponents[0], Math.log(1.9), 1e-14, 'tent');
  assert.ok(mapLyapunov(logisticMap(3.2), [0.3]).exponents[0] < -0.5, 'periodic window: negative');
  const h = mapLyapunov(henonMap(), [0.1, 0.1], { n: 200000 });
  close(h.exponents[0], 0.41922, 5e-3, 'Hénon λ1 (Sprott 2003)');
  close(h.sum, Math.log(0.3), 1e-12, 'Hénon λ1 + λ2 = ln b exactly');
  const s = mapLyapunov(standardMap(5, { wrap: false }), [0.5, 0.2], { n: 20000 });
  close(s.sum, 0, 1e-12, 'standard map: λ1 + λ2 = 0');
  close(s.exponents[0], Math.log(5 / 2), 0.15, 'standard map K = 5: λ ≈ ln(K/2) (Chirikov 1979)');
  assert.equal(h.history.length, 10);
});

test('circle map: rotation number Ω at K = 0, locked to 1/2 in the tongue at K = 1', () => {
  close(rotationNumber(circleMap(0.3141, 0), [0.1]), 0.3141, 1e-12);
  close(rotationNumber(circleMap(0.5, 1), [0.1], 20000), 0.5, 1e-4);
  close(rotationNumber(circleMap(0.49, 0.95), [0.1], 20000), 0.5, 1e-4, 'mode-locked plateau');
});

// ---------------------------------------------------------------- flows

test('flows: Jacobians against finite differences; Lorenz divergence; Hénon–Heiles and double pendulum conserve energy', () => {
  const y = [0.3, -0.2, 0.25, 0.1];
  for (const [s, x, t] of [[lorenz(), [1, 2, 3], 0], [rossler(), [1, -2, 0.5], 0], [duffing(), [0.4, -0.3], 1.3], [henonHeiles(), y, 0], [doublePendulum(), [0.5, -0.3, 0.2, 1.1], 0]]) {
    const J = s.jacobian(t, x), e = 1e-6;
    for (let j = 0; j < s.dim; j++) {
      const p = Float64Array.from(x), q = Float64Array.from(x); p[j] += e; q[j] -= e;
      const fp = s.f(t, p), fq = s.f(t, q);
      for (let i = 0; i < s.dim; i++) close(J[i][j], (fp[i] - fq[i]) / (2 * e), 1e-7, `${s.name} J[${i}][${j}]`);
    }
  }
  const J = lorenz().jacobian(0, [3, -1, 7]);
  close(J[0][0] + J[1][1] + J[2][2], -(10 + 1 + 8 / 3), 1e-15);
  for (const [s, y0] of [[henonHeiles(), [0, 0.1, 0.35, 0.1]], [doublePendulum(), [2, 2.5, 0, 0]]]) {
    const st = odeStepper({ f: s.f, y0, h: 0.01, rtol: 1e-12, atol: 1e-14 });
    st.advanceTo(20);
    close(s.energy(st.y), s.energy(y0), 1e-9 * Math.max(1, Math.abs(s.energy(y0))), `${s.name} energy`);
  }
});

test('Lorenz Lyapunov spectrum: λ ≈ (0.906, 0, −14.57), Σλ = −(σ + 1 + β) exactly', () => {
  const r = lyapunovSpectrum({ system: lorenz(), y0: [1, 1, 20], duration: 300, transient: 10, rtol: 1e-9, atol: 1e-11 });
  close(r.exponents[0], 0.9056, 0.02, 'λ1 (Sprott 2003)');
  close(r.exponents[1], 0, 0.01, 'λ2: flow direction');
  close(r.exponents[2], -14.5723, 0.03, 'λ3');
  close(r.sum, -(10 + 1 + 8 / 3), 1e-5, 'trace of J');
  // a linear flow: exponents are the real parts of the eigenvalues, exactly
  const lin = { name: 'linear', dim: 2, f: (_t, [x, v]) => Float64Array.of(0.2 * x + v, -0.5 * v), jacobian: () => [Float64Array.of(0.2, 1), Float64Array.of(0, -0.5)] };
  const l = lyapunovSpectrum({ system: lin, y0: [1, 1], duration: 400, interval: 2 });
  close(l.exponents[0], 0.2, 5e-3); close(l.exponents[1], -0.5, 5e-3); close(l.sum, -0.3, 1e-8);
});

test('Poincaré section: decoupled oscillators give the exact rotation; Hénon–Heiles points stay on the energy shell', () => {
  const w = Math.SQRT2;
  const sys = { name: 'osc', dim: 4, f: (_t, [x, y, px, py]) => Float64Array.of(px, py, -x, -w * w * y) };
  const { points, times } = poincareSection({ system: sys, y0: [0.5, 0.3, 0, 0.2], section: (_t, y) => y[0], direction: 1, count: 6 });
  // x = 0.5 cos t crosses 0 upward at t = 3π/2 + 2πk; y(t) = 0.3 cos ωt + (0.2/ω) sin ωt
  points.forEach((p, k) => {
    const t = 1.5 * Math.PI + 2 * Math.PI * k;
    close(times[k], t, 1e-8, 'crossing time');
    close(p[1], 0.3 * Math.cos(w * t) + 0.2 / w * Math.sin(w * t), 1e-8, 'y');
    close(p[3], -0.3 * w * Math.sin(w * t) + 0.2 * Math.cos(w * t), 1e-8, 'py');
    assert.ok(Math.abs(p[0]) < 1e-10 && p[2] > 0);
  });
  const hh = henonHeiles(), E = 1 / 8, y0 = [0, 0.1, 0, 0.3];
  y0[2] = Math.sqrt(2 * (E - hh.energy(y0)));
  const sec = poincareSection({ system: hh, y0, section: (_t, y) => y[0], direction: 1, count: 50 });
  assert.equal(sec.points.length, 50);
  for (const p of sec.points) { close(hh.energy(p), E, 1e-9, 'H'); close(p[0], 0, 1e-9, 'on the section'); }
});

// ---------------------------------------------------------------- dimensions

test('fractal dimensions: Cantor set exact; Hénon and Lorenz against Grassberger–Procaccia and Russell et al.', () => {
  // midpoints of the 2¹² level-12 middle-thirds intervals: N(3⁻ᵏ) = 2ᵏ exactly
  let pts = [0.5];
  for (let lvl = 1; lvl <= 12; lvl++) pts = pts.flatMap(x => [x / 3, 2 / 3 + x / 3]);
  const cantor = boxCountingDimension(pts.map(x => [x]), { sizes: [2, 3, 4, 5, 6, 7, 8].map(k => 3 ** -k) });
  assert.deepEqual(Array.from(cantor.counts), [4, 8, 16, 32, 64, 128, 256]);
  close(cantor.dimension, Math.log(2) / Math.log(3), 1e-12);
  close(correlationDimension(pts.map(x => [x]), { rMin: 1e-3, rMax: 0.1 }).dimension, Math.log(2) / Math.log(3), 0.05, 'Cantor D2');
  // uniform random points in a square → 2 (a smooth set)
  const grid = [];
  for (let i = 0; i < 50; i++) for (let j = 0; j < 50; j++) grid.push([(i + 0.5 * ((j * 7) % 3) / 3) / 50, (j + 0.37) / 50]);
  close(correlationDimension(grid, { rMin: 0.05, rMax: 0.2 }).dimension, 2, 0.1, 'square');
  const henon = pointsOf(iterateMap(henonMap(), [0.1, 0.1], 3000, { transient: 100 }), 2);
  close(correlationDimension(henon, { rMin: 1e-3, rMax: 0.1 }).dimension, 1.21, 0.05, 'Hénon D2');
  const many = pointsOf(iterateMap(henonMap(), [0.1, 0.1], 100000, { transient: 100 }), 2);
  close(boxCountingDimension(many, { sizes: [0.05, 0.02, 0.01, 0.005] }).dimension, 1.261, 0.06, 'Hénon D0');
  const s = odeStepper({ f: lorenz().f, y0: [1, 1, 20], h: 0.01, rtol: 1e-9, atol: 1e-11 });
  s.advanceTo(20);
  const lp = [];
  for (let k = 1; k <= 2500; k++) { s.advanceTo(20 + 0.1 * k); lp.push(s.sampleAt(20 + 0.1 * k)); }
  close(correlationDimension(lp, { rMin: 0.3, rMax: 3, theiler: 5 }).dimension, 2.05, 0.08, 'Lorenz D2');
  assert.throws(() => correlationDimension(lp, { rMin: 3, rMax: 1 }), /rMin/);
});

// ---------------------------------------------------------------- fixed points and periodic orbits

test('eigenvalues of small matrices: closed forms, a rotation, companion matrices', () => {
  const e2 = smallEigenvalues([[0, 1], [-4, 0]]);
  close(e2[0].re, 0, 1e-15); close(Math.abs(e2[0].im), 2, 1e-15);
  const e3 = smallEigenvalues([[2, 0, 0], [0, -1, 0], [1, 0, 5]]);
  assert.deepEqual(e3.map(e => Math.round(e.re * 1e12) / 1e12), [5, 2, -1]);
  // companion of (λ − 1)(λ − 2)(λ + 3)(λ² + 1) = λ⁵ − 2λ³ ... built from its roots
  const roots = [1, 2, -3], poly = [1];                    // coefficients high → low
  const mulBy = c => { const out = Array(poly.length + 1).fill(0); poly.forEach((p, i) => { out[i] += p; out[i + 1] -= c * p; }); poly.splice(0, poly.length, ...out); };
  roots.forEach(mulBy);
  const quad = [1, 0, 1], full = Array(poly.length + 2).fill(0);
  poly.forEach((p, i) => quad.forEach((q, j) => { full[i + j] += p * q; }));
  const n = full.length - 1, C = Array.from({ length: n }, (_, i) => Float64Array.from({ length: n }, (_c, j) => (i === 0 ? -full[j + 1] : i === j + 1 ? 1 : 0)));
  const ev = smallEigenvalues(C);
  close(ev[0].re, 2, 1e-10); close(ev[1].re, 1, 1e-10); close(ev.at(-1).re, -3, 1e-10);
  close(Math.abs(ev[2].im), 1, 1e-10); close(ev[2].re, 0, 1e-10);
  assert.throws(() => smallEigenvalues([]), /size/);
});

test('Lorenz fixed points: exact positions, origin saddle eigenvalues, Hopf at ρ_H = σ(σ+β+3)/(σ−β−1)', () => {
  const sys = lorenz(), c = Math.sqrt(8 / 3 * 27);
  const o = findFixedPoint(sys, [0.1, -0.1, 0.2]);
  assert.ok(o.converged); assert.equal(o.type, 'saddle'); assert.equal(o.unstableDimension, 1);
  const disc = Math.sqrt(11 ** 2 + 4 * 10 * 27);
  close(o.eigenvalues[0].re, (-11 + disc) / 2, 1e-10); close(o.eigenvalues[1].re, -8 / 3, 1e-10); close(o.eigenvalues[2].re, (-11 - disc) / 2, 1e-10);
  const p = findFixedPoint(sys, [8, 8, 26]);
  [c, c, 27].forEach((v, i) => close(p.x[i], v, 1e-10));
  assert.equal(p.type, 'saddle-focus');                     // ρ = 28 > ρ_H ≈ 24.74: unstable complex pair
  const rhoH = 10 * (10 + 8 / 3 + 3) / (10 - 8 / 3 - 1);
  for (const [rho, type] of [[rhoH - 0.5, 'stable focus'], [rhoH + 0.5, 'saddle-focus']]) {
    const cc = Math.sqrt(8 / 3 * (rho - 1));
    assert.equal(findFixedPoint(lorenz({ rho }), [cc, cc, rho - 1]).type, type, `ρ = ${rho}`);
  }
  // double pendulum: small-angle normal modes ω² = (g/l)(2 ∓ √2) for equal masses and lengths
  const dp = findFixedPoint(doublePendulum(), [0.01, -0.01, 0, 0]);
  assert.equal(dp.type, 'center');
  const om = dp.eigenvalues.map(e => Math.abs(e.im)).sort((a, b) => a - b);
  close(om[0], Math.sqrt(9.81 * (2 - Math.SQRT2)), 1e-6); close(om[2], Math.sqrt(9.81 * (2 + Math.SQRT2)), 1e-6);
});

test('periodic orbits of maps: Hénon fixed points, logistic 2-cycle multiplier, standard-map elliptic point', () => {
  const a = 1.4, b = 0.3, xs = (-(1 - b) + Math.sqrt((1 - b) ** 2 + 4 * a)) / (2 * a);
  const h = findPeriodicOrbit(henonMap(), [0.6, 0.2], 1);
  close(h.x[0], xs, 1e-12); close(h.x[1], b * xs, 1e-12); assert.equal(h.type, 'saddle');
  const r = 3.3, l = findPeriodicOrbit(logisticMap(r), [0.8], 2);
  close(l.multipliers[0].re, -r * r + 2 * r + 4, 1e-9, 'multiplier −r² + 2r + 4');
  assert.equal(l.type, 'stable node');
  close(l.x[0], (r + 1 + Math.sqrt((r - 3) * (r + 1))) / (2 * r), 1e-12);
  // standard map K = 0.5: (π, 0) is elliptic (|tr J| = |2 − K| < 2), (0, 0) hyperbolic
  const e = findPeriodicOrbit(standardMap(0.5, { wrap: false }), [3, 0.01], 1);
  close(e.x[0], Math.PI, 1e-12); assert.equal(e.type, 'elliptic');
  assert.equal(findPeriodicOrbit(standardMap(0.5, { wrap: false }), [0.1, 0.01], 1).type, 'saddle');
  assert.equal(classifyFixedPoint([[0.5, 0], [0, 0.2]], { kind: 'map' }).type, 'stable node');
  assert.throws(() => classifyFixedPoint([[1]], { kind: 'orbit' }), /kind/);
});
