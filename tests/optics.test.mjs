import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  fresnelCS, singleSlit, doubleSlit, grating, circularAperture, fresnelEdge, fraunhofer, angularSpectrum,
  gaussianBeam, beamFromQ, propagateQ, abcd, abcdSystem, cavityMode,
  reflect, refract, fresnelCoefficients, brewsterAngle, criticalAngle, traceRay, paraxialSystem,
  sellmeier, GLASS, abbeNumber, rainbow, thinFilm, jones, fdtd2D, yeeDispersion,
} from '../core/optics.mjs';
import { gaussLegendre, besselJ, besselJZeros } from '../core/numerics.mjs';
import { c } from '../core/constants.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const deg = Math.PI / 180;

// ---------------------------------------------------------------- Fresnel integrals

test('fresnelCS: Abramowitz & Stegun Table 7.7 values, odd symmetry, limits ½', () => {
  // A&S Table 7.7 (C(x), S(x) to 7 decimals)
  for (const [x, C, S] of [[0.5, 0.4923442, 0.0647324], [1, 0.7798934, 0.4382591], [1.5, 0.4452612, 0.6975050], [2, 0.4882534, 0.3434157]]) {
    const r = fresnelCS(x);
    close(r.C, C, 6e-8, `C(${x})`); close(r.S, S, 6e-8, `S(${x})`);
    const m = fresnelCS(-x);
    assert.equal(m.C, -r.C); assert.equal(m.S, -r.S);
  }
  assert.deepEqual(fresnelCS(Infinity), { C: 0.5, S: 0.5 });
  assert.deepEqual(fresnelCS(0), { C: 0, S: 0 });
  // asymptotic C ≈ ½ + sin(πx²/2)/(πx) − cos(πx²/2)/(π²x³) (A&S 7.3.9–10 leading terms)
  for (const x of [20, 37.3]) {
    const a = Math.PI * x * x / 2, r = fresnelCS(x);
    close(r.C, 0.5 + Math.sin(a) / (Math.PI * x) - Math.cos(a) / (Math.PI ** 2 * x ** 3), 3 / (Math.PI ** 3 * x ** 5), 'asymptotic C');
    close(r.S, 0.5 - Math.cos(a) / (Math.PI * x) - Math.sin(a) / (Math.PI ** 2 * x ** 3), 3 / (Math.PI ** 3 * x ** 5), 'asymptotic S');
  }
  assert.throws(() => fresnelCS(NaN));
});

test('fresnelCS: agrees with Gauss–Legendre quadrature across the series/continued-fraction switch', () => {
  const quad = (f, x) => { let s = 0; const p = 40; for (let k = 0; k < p; k++) s += gaussLegendre(f, k * x / p, (k + 1) * x / p, 32); return s; };
  for (const x of [0.1, 0.9, 1.49, 1.5, 1.51, 2.7, 4.2, 6]) {
    const r = fresnelCS(x);
    close(r.C, quad(t => Math.cos(Math.PI * t * t / 2), x), 2e-14, `C(${x})`);
    close(r.S, quad(t => Math.sin(Math.PI * t * t / 2), x), 2e-14, `S(${x})`);
  }
});

// ---------------------------------------------------------------- closed-form diffraction

test('slits and gratings: zeros at a sinθ = mλ, N = 1 and N = 2 limits, N − 1 zeros between orders', () => {
  const lambda = 500e-9, a = 2e-6, d = 7e-6;
  for (const m of [1, 2, 3]) close(singleSlit(Math.asin(m * lambda / a), { a, lambda }), 0, 1e-28);
  close(singleSlit(0, { a, lambda }), 1, 0);
  for (const th of [0.01, 0.05, 0.13, 0.2]) {
    close(grating(th, { a, d, N: 1, lambda }), singleSlit(th, { a, lambda }), 1e-15);
    const g = Math.PI * d * Math.sin(th) / lambda;
    close(doubleSlit(th, { a, d, lambda }), singleSlit(th, { a, lambda }) * Math.cos(g) ** 2, 1e-14);
  }
  const N = 6;
  for (let m = 1; m <= 2; m++) {
    const th = Math.asin(m * lambda / d);                      // principal maximum: envelope value
    close(grating(th, { a, d, N, lambda }), singleSlit(th, { a, lambda }), 1e-12);
    for (let p = 1; p < N; p++) close(grating(Math.asin((m + p / N) * lambda / d), { a, d, N, lambda }), 0, 1e-25);
  }
  assert.throws(() => grating(0, { a: 3e-6, d: 2e-6, N: 2, lambda }));
});

test('Airy pattern: first zero at 1.21967 λ/D, 83.8 % encircled energy (Born & Wolf §8.5.2)', () => {
  const lambda = 633e-9, D = 1e-3, j11 = besselJZeros(1, 1)[0];
  close(j11 / Math.PI, 1.219670, 1e-6);
  close(circularAperture(Math.asin(j11 / Math.PI * lambda / D), { D, lambda }), 0, 1e-28);
  close(circularAperture(0, { D, lambda }), 1, 0);
  const I = v => circularAperture(Math.asin(v * lambda / (Math.PI * D)), { D, lambda });
  let enc = 0;
  for (let k = 0; k < 8; k++) enc += gaussLegendre(v => I(v) * v / 2, k * j11 / 8, (k + 1) * j11 / 8, 24);
  close(enc, 1 - besselJ(0, j11) ** 2 - besselJ(1, j11) ** 2, 1e-12, 'Rayleigh encircled energy');
  close(enc, 0.838, 5e-4);
});

test('knife edge: I/I0 = 1/4 at the shadow edge, first maximum 1.3704 at v = 1.2172 (Born & Wolf §8.7)', () => {
  const lambda = 500e-9, z = 1, s = Math.sqrt(lambda * z / 2);   // x = v·s
  close(fresnelEdge(0, { z, lambda }), 0.25, 1e-15);
  let best = 0, vBest = 0;
  for (let v = 1.1; v < 1.35; v += 1e-5) { const I = fresnelEdge(v * s, { z, lambda }); if (I > best) { best = I; vBest = v; } }
  close(best, 1.3704, 1e-4); close(vBest, 1.2172, 2e-4);
  close(fresnelEdge(1e3 * s, { z, lambda }), 1, 1e-3);
  const v = 8; close(fresnelEdge(-v * s, { z, lambda }), 1 / (2 * Math.PI ** 2 * v * v), 1e-5, 'deep shadow asymptote');
});

// ---------------------------------------------------------------- FFT optics

test('fraunhofer: exact Dirichlet kernel of a sampled slit, sinc² near axis, separable 2D rectangle', () => {
  const n = 1024, M = 33, dx = 1e-6, lambda = 1e-6;
  const re = new Float64Array(n);
  for (let i = n / 2 - 16; i <= n / 2 + 16; i++) re[i] = 1;
  const r = fraunhofer({ re, shape: [n], dx, lambda });
  for (let k = 0; k < n; k++) {
    const q = k - n / 2, exact = q === 0 ? M * dx : dx * Math.sin(Math.PI * q * M / n) / Math.sin(Math.PI * q / n);
    close(r.re[k], exact, 1e-18, `bin ${q}`); close(r.im[k], 0, 1e-18);
  }
  for (const q of [1, 5, 12, 30]) close(r.normalized[n / 2 + q], singleSlit(Math.asin(r.sinTheta[n / 2 + q]), { a: M * dx, lambda }), 2e-3);
  // 2D rectangle = product of the 1D transforms
  const nx = 32, ny = 16, re2 = new Float64Array(nx * ny), ax = new Float64Array(nx), ay = new Float64Array(ny);
  for (let i = 13; i <= 19; i++) ax[i] = 1;
  for (let j = 6; j <= 9; j++) ay[j] = 1;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) re2[i + nx * j] = ax[i] * ay[j];
  const fx = fraunhofer({ re: ax, shape: [nx], dx, lambda }), fy = fraunhofer({ re: ay, shape: [ny], dx, lambda });
  const f2 = fraunhofer({ re: re2, shape: [nx, ny], dx, lambda });
  for (let j = 0; j < ny; j++) { for (let i = 0; i < nx; i++) {
    const pr = fx.re[i] * fy.re[j] - fx.im[i] * fy.im[j], pi = fx.re[i] * fy.im[j] + fx.im[i] * fy.re[j];
    close(f2.re[i + nx * j], pr, 1e-24); close(f2.im[i + nx * j], pi, 1e-24);
  } }
  assert.throws(() => fraunhofer({ re: new Float64Array(100), shape: [100], dx, lambda }));
});

test('angularSpectrum: Gaussian beam w(z) = w0√(1+(z/zR)²) and exact power conservation', () => {
  const lambda = 1e-6, w0 = 8 * lambda, n = 128, dx = lambda / 2, re = new Float64Array(n * n);
  const x = i => (i - n / 2) * dx;
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) re[i + n * j] = Math.exp(-(x(i) ** 2 + x(j) ** 2) / w0 ** 2);
  const beam = gaussianBeam({ w0, lambda }), power = (a, b) => a.reduce((s, v, k) => s + v * v + b[k] * b[k], 0);
  const p0 = power(re, new Float64Array(n * n));
  for (const z of [beam.zR / 2, beam.zR]) {
    const u = angularSpectrum({ re, shape: [n, n], dx, lambda, z });
    const P = power(u.re, u.im);
    close(P, p0, 1e-12 * p0, 'power');
    let m2 = 0;
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const k = i + n * j; m2 += x(i) ** 2 * (u.re[k] ** 2 + u.im[k] ** 2); }
    const w = 2 * Math.sqrt(m2 / P);                     // ⟨x²⟩ = w²/4
    close(w / beam.w(z), 1, 2e-3, `w at z = ${z}`);         // non-paraxial correction O((λ/w0)²)
  }
});

test('angularSpectrum: Talbot self-imaging of a cosine grating (exact distance and 2d²/λ paraxially)', () => {
  const lambda = 1e-6, d = 8 * lambda, n = 64, dx = d / 8, re = new Float64Array(n);
  for (let i = 0; i < n; i++) re[i] = 1 + 0.5 * Math.cos(2 * Math.PI * i * dx / d);
  const zT = lambda / (1 - Math.sqrt(1 - (lambda / d) ** 2));     // exact revival for orders 0, ±1
  close(zT / (2 * d * d / lambda), 1, 5e-3, 'paraxial Talbot length');
  const u = angularSpectrum({ re, shape: [n], dx, lambda, z: zT });
  const ph = Math.atan2(u.im[0], u.re[0]);                   // global phase e^{ikz}
  const g = { re: Math.cos(-ph), im: Math.sin(-ph) };
  for (let i = 0; i < n; i++) {
    close(u.re[i] * g.re - u.im[i] * g.im, re[i], 1e-11); close(u.re[i] * g.im + u.im[i] * g.re, 0, 1e-11);
  }
  const h = angularSpectrum({ re, shape: [n], dx, lambda, z: zT / 2 });  // half Talbot: shifted by d/2
  for (let i = 0; i < n; i++) close(Math.hypot(h.re[i], h.im[i]), 1 - 0.5 * Math.cos(2 * Math.PI * i * dx / d), 1e-11);
});

// ---------------------------------------------------------------- Gaussian beams and resonators

test('Gaussian beam: q through free space and a thin lens (waist λf/(πw0) at the back focal plane)', () => {
  const lambda = 1064e-9, w0 = 0.5e-3, f = 0.2, b = gaussianBeam({ w0, lambda });
  const q = propagateQ(b.q(0), abcd.freeSpace(0.7));
  close(q.re, 0.7, 1e-15); close(q.im, b.zR, 1e-15);
  const p = beamFromQ(q, lambda);
  close(p.w, b.w(0.7), 1e-15); close(p.R, b.R(0.7), 1e-9); close(p.w0, w0, 1e-15);
  const M = abcdSystem(abcd.freeSpace(f), abcd.thinLens(f), abcd.freeSpace(f));
  const out = beamFromQ(propagateQ(b.q(0), M), lambda);
  close(out.z, 0, 1e-15, 'waist at back focal plane'); close(out.w0, lambda * f / (Math.PI * w0), 1e-15);
  close(out.w, out.w0, 1e-15);
  assert.throws(() => beamFromQ({ re: 1, im: -1 }, lambda));
});

test('cavityMode: confocal waist √(λL/2π), round-trip self-consistency, mirrors match wavefronts, Gouy = acos√(g1g2)', () => {
  const lambda = 633e-9, L = 0.3;
  close(cavityMode({ L, R1: L, R2: L, lambda }).w0, Math.sqrt(lambda * L / (2 * Math.PI)), 1e-15);
  const R1 = 1.0, R2 = 0.5, m = cavityMode({ L, R1, R2, lambda });
  const q1 = { re: -m.waistFromMirror1, im: m.zR };              // just after mirror 1, heading to mirror 2
  const rt = abcdSystem(abcd.freeSpace(L), abcd.curvedMirror(R2), abcd.freeSpace(L), abcd.curvedMirror(R1));
  const back = propagateQ(q1, rt);
  close(back.re, q1.re, 1e-12); close(back.im, q1.im, 1e-12);
  const at2 = beamFromQ(propagateQ(q1, abcd.freeSpace(L)), lambda);
  close(at2.R, R2, 1e-12, 'wavefront radius equals mirror 2 radius'); close(at2.w, m.w2, 1e-15);
  close(beamFromQ(q1, lambda).R, -R1, 1e-12, 'mirror 1 (diverging away from it)');
  close(m.gouy, Math.acos(Math.sqrt(m.g1 * m.g2)), 1e-13);
  const plano = cavityMode({ L, R1: Infinity, R2: 1, lambda });
  close(plano.waistFromMirror1, 0, 0);
  assert.throws(() => cavityMode({ L: 1, R1: 0.4, R2: 0.4, lambda }), /not stable/);
  assert.throws(() => cavityMode({ L, R1: Infinity, R2: Infinity, lambda }));
});

// ---------------------------------------------------------------- interfaces

test('refract/reflect: Snell law in 3D, TIR returns null, mirror law', () => {
  const n1 = 1, n2 = 1.5, d = [0.3, -0.2, 0.9], N = [0, 0, -1];
  const t = refract(d, N, n1, n2), u = d.map(v => v / Math.hypot(...d));
  const sin1 = Math.hypot(u[0], u[1]), sin2 = Math.hypot(t[0], t[1]);
  close(n1 * sin1, n2 * sin2, 1e-15); close(t[0] * u[1] - t[1] * u[0], 0, 1e-15, 'plane of incidence');
  assert.equal(refract([Math.sin(50 * deg), Math.cos(50 * deg)], [0, 1], 1.5, 1), null);
  const r = reflect([1, -1], [0, 1]);
  close(r[0], Math.SQRT1_2, 1e-15); close(r[1], Math.SQRT1_2, 1e-15);
});

test('fresnelCoefficients: R + T = 1, normal incidence, Brewster, TIR phases (Fresnel rhomb, Hecht §8.6)', () => {
  for (const [n1, n2] of [[1, 1.5], [1.5, 1], [1.33, 2.4]]) {
    for (let th = 0; th < 89.5 * deg; th += 3.7 * deg) {
      const f = fresnelCoefficients(th, n1, n2);
      close(f.Rs + f.Ts, 1, 1e-14); close(f.Rp + f.Tp, 1, 1e-14);
      if (f.tir) { close(f.Rs, 1, 1e-14); close(f.Rp, 1, 1e-14); }
    }
    const f0 = fresnelCoefficients(0, n1, n2);
    close(f0.Rs, ((n1 - n2) / (n1 + n2)) ** 2, 1e-15); close(f0.rp.re, -f0.rs.re, 1e-15);
    close(fresnelCoefficients(brewsterAngle(n1, n2), n1, n2).Rp, 0, 1e-28);
  }
  close(criticalAngle(1.5, 1), Math.asin(1 / 1.5), 1e-15);
  assert.throws(() => criticalAngle(1, 1.5));
  const f = fresnelCoefficients(54.6 * deg, 1.51, 1);
  const delta = Math.atan2(f.rp.im, f.rp.re) - Math.atan2(f.rs.im, f.rs.re);
  close(Math.abs(Math.atan2(Math.sin(delta), Math.cos(delta))) / deg, 45, 0.1, 'rhomb retardance per reflection');
});

// ---------------------------------------------------------------- ray tracing

test('traceRay/paraxialSystem: thick-lens lensmaker focal length, paraxial focus, undercorrected spherical aberration', () => {
  const n = 1.5168, R1 = 0.05, R2 = -0.05, t = 0.005;
  const lens = [{ z: 0, R: R1, n }, { z: t, R: R2, n: 1 }];
  const { efl, bfd } = paraxialSystem(lens);
  close(1 / efl, (n - 1) * (1 / R1 - 1 / R2 + (n - 1) * t / (n * R1 * R2)), 1e-12, 'lensmaker');
  const cross = h => { const r = traceRay({ p: [0, h, -0.01], d: [0, 0, 1] }, lens); assert.equal(r.status, 'ok'); return r.p[2] - r.p[1] * r.d[2] / r.d[1]; };
  close((cross(1e-5) - t) / bfd, 1, 1e-7, 'paraxial ray meets the back focal point');
  const lsa = h => cross(h) - cross(1e-5);
  assert.ok(lsa(0.01) < 0 && lsa(0.015) < lsa(0.01), 'positive lens: marginal rays focus short, growing with height');
  close(lsa(0.004) / lsa(0.002), 4, 0.05, 'third-order longitudinal SA ∝ h²');
  // a flat surface reproduces Snell, apertures clip, a strongly curved surface can be missed
  const s = traceRay({ p: [0, 0, -1], d: [0, Math.sin(0.4), Math.cos(0.4)] }, [{ z: 0, R: Infinity, n: 1.5 }]);
  close(Math.hypot(s.d[0], s.d[1]), Math.sin(0.4) / 1.5, 1e-15);
  assert.equal(traceRay({ p: [0, 0.02, 0], d: [0, 0, 1] }, [{ z: 0.01, R: 0.05, n, aperture: 0.01 }]).status, 'clipped');
  assert.equal(traceRay({ p: [0, 0.2, 0], d: [0, 0, 1] }, [{ z: 0.01, R: 0.05, n }]).status, 'missed');
  assert.equal(traceRay({ p: [0, 0, 0], d: [0, 0.8, 0.6] }, [{ z: 0.01, R: Infinity, n: 1 }], 1.5).status, 'tir');   // inside glass at 53° > 41.8°
});

// ---------------------------------------------------------------- dispersion and rainbow

test('sellmeier: Schott N-BK7 n_d = 1.5168, ν_d = 64.17; fused silica n_d = 1.4585 (Malitson 1965)', () => {
  close(sellmeier(587.5618e-9, GLASS.BK7), 1.5168, 5e-5);
  close(abbeNumber(GLASS.BK7), 64.17, 0.05);
  close(sellmeier(587.5618e-9, GLASS.fusedSilica), 1.4585, 5e-5);
  close(abbeNumber(GLASS.fusedSilica), 67.8, 0.2);
  assert.ok(sellmeier(400e-9, GLASS.BK7) > sellmeier(700e-9, GLASS.BK7), 'normal dispersion');
  assert.throws(() => sellmeier(10e-6, GLASS.BK7), /range/);
});

test('rainbow: Descartes minimum deviation, primary ≈ 42.0° and secondary ≈ 50.9° for n = 1.333', () => {
  const n = 1.333, D = (i, k) => { const r = Math.asin(Math.sin(i) / n); return 2 * (i - r) + k * (Math.PI - 2 * r); };
  for (const [k, angle] of [[1, 42.0], [2, 50.9]]) {
    const rb = rainbow(n, k);
    close(rb.angle / deg, angle, 0.1);
    let min = Infinity;
    for (let i = 0; i < Math.PI / 2; i += 1e-5) min = Math.min(min, D(i, k));
    close(rb.deviation, min, 1e-9, 'stationary point is the minimum');
  }
  close(rainbow(n).deviation / deg, 138.0, 0.1);
  assert.throws(() => rainbow(0.9));
});

// ---------------------------------------------------------------- multilayers

test('thinFilm: bare interface equals Fresnel; lossless R + T = 1; TIR substrate gives T = 0', () => {
  for (const pol of ['s', 'p']) {
    for (let th = 0; th < 85 * deg; th += 7 * deg) {
      const f = fresnelCoefficients(th, 1.2, 1.7), m = thinFilm({ lambda: 1e-6, n0: 1.2, ns: 1.7, theta0: th, polarization: pol });
      const r = pol === 's' ? f.rs : f.rp, t = pol === 's' ? f.ts : f.tp;
      close(m.r.re, r.re, 1e-14); close(m.r.im, r.im, 1e-14); close(m.t.re, t.re, 1e-14); close(m.t.im, t.im, 1e-14);
      close(m.T, pol === 's' ? f.Ts : f.Tp, 1e-14);
    }
    const g = fresnelCoefficients(60 * deg, 1.5, 1), m = thinFilm({ lambda: 1e-6, n0: 1.5, ns: 1, theta0: 60 * deg, polarization: pol });
    const r = pol === 's' ? g.rs : g.rp;
    close(m.R, 1, 1e-14); close(m.T, 0, 0); close(m.r.re, r.re, 1e-14); close(m.r.im, r.im, 1e-14);
    const layers = [{ n: 2.1, d: 137e-9 }, { n: 1.38, d: 211e-9 }, { n: 1.9, d: 50e-9 }, { n: 1.46, d: 400e-9 }];
    for (const th of [0, 0.3, 1.1]) { const s = thinFilm({ lambda: 633e-9, layers, ns: 1.52, theta0: th, polarization: pol }); close(s.R + s.T, 1, 1e-13); }
  }
  close(thinFilm({ lambda: 1e-6, n0: 1, ns: 1.5, theta0: brewsterAngle(1, 1.5), polarization: 'p' }).R, 0, 1e-28);
});

test('thinFilm: quarter-wave AR coating R = 0, Bragg mirror closed form, single film Airy formula', () => {
  const lambda = 550e-9, ns = 1.52, nAR = Math.sqrt(ns);
  close(thinFilm({ lambda, layers: [{ n: nAR, d: lambda / (4 * nAR) }], ns }).R, 0, 1e-30);
  const nH = 2.35, nL = 1.38, N = 6, layers = [];
  for (let k = 0; k < N; k++) layers.push({ n: nH, d: lambda / (4 * nH) }, { n: nL, d: lambda / (4 * nL) });
  const Y = ns * (nH / nL) ** (2 * N);                       // Macleod eq. for (HL)^N on a substrate
  close(thinFilm({ lambda, layers, ns }).R, ((1 - Y) / (1 + Y)) ** 2, 1e-12);
  // single film n1, thickness d, normal incidence: Airy summation (Born & Wolf §1.6.4)
  const n1 = 1.9, d = 313e-9, r01 = (1 - n1) / (1 + n1), r12 = (n1 - ns) / (n1 + ns), b = 2 * Math.PI * n1 * d / lambda;
  const R = (r01 ** 2 + r12 ** 2 + 2 * r01 * r12 * Math.cos(2 * b)) / (1 + r01 ** 2 * r12 ** 2 + 2 * r01 * r12 * Math.cos(2 * b));
  close(thinFilm({ lambda, layers: [{ n: n1, d }], ns }).R, R, 1e-14);
});

test('thinFilm: oblique incidence on an absorbing substrate matches Fresnel with complex N (both polarisations)', () => {
  const C = (re, im = 0) => ({ re, im }), mul = (a, b) => C(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
  const div = (a, b) => { const d = b.re ** 2 + b.im ** 2; return C((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
  const sqrt = a => { const r = Math.hypot(a.re, a.im); return C(Math.sqrt((r + a.re) / 2), Math.sign(a.im || 1) * Math.sqrt((r - a.re) / 2)); };
  const abs2 = a => a.re ** 2 + a.im ** 2, N = C(1.6, 0.9), n1 = 1.33;
  for (const th of [0.2, 0.7, 1.2]) {
    const s1 = n1 * Math.sin(th), c1 = Math.cos(th), Nc2 = sqrt(C(N.re ** 2 - N.im ** 2 - s1 * s1, 2 * N.re * N.im));   // N cosθ2
    const Rs = abs2(div(C(n1 * c1 - Nc2.re, -Nc2.im), C(n1 * c1 + Nc2.re, Nc2.im)));
    const NNc1 = mul(mul(N, N), C(c1)), n1Nc2 = mul(C(n1), Nc2);         // N²cosθ1 vs n1·N cosθ2 (p, scaled by N)
    const Rp = abs2(div(C(NNc1.re - n1Nc2.re, NNc1.im - n1Nc2.im), C(NNc1.re + n1Nc2.re, NNc1.im + n1Nc2.im)));
    close(thinFilm({ lambda: 1e-6, n0: n1, ns: N, theta0: th, polarization: 's' }).R, Rs, 1e-14);
    close(thinFilm({ lambda: 1e-6, n0: n1, ns: N, theta0: th, polarization: 'p' }).R, Rp, 1e-14);
  }
});

test('thinFilm: absorbing media (metal R formula), frustrated TIR conserves energy and decays with the gap', () => {
  const N = { re: 0.2, im: 3.4 };                             // silver-like at visible
  close(thinFilm({ lambda: 600e-9, ns: N }).R, ((N.re - 1) ** 2 + N.im ** 2) / ((N.re + 1) ** 2 + N.im ** 2), 1e-14);
  const film = thinFilm({ lambda: 600e-9, layers: [{ n: N, d: 20e-9 }], ns: 1.5 });
  assert.ok(film.A > 0.01 && film.A < 1 && film.T > 0, 'thin metal absorbs and transmits');
  let prev = 1;
  for (const gap of [50e-9, 200e-9, 400e-9]) {
    const f = thinFilm({ lambda: 600e-9, n0: 1.5, layers: [{ n: 1, d: gap }], ns: 1.5, theta0: 60 * deg, polarization: 's' });
    close(f.R + f.T, 1, 1e-13); assert.ok(f.T < prev); prev = f.T;
  }
});

// ---------------------------------------------------------------- polarisation

test('jones: Malus law, three polarisers give 1/8, waveplates (QWP → circular, HWP rotates by 2θ), Stokes', () => {
  for (const th of [0, 0.3, 1.1]) close(jones.intensity(jones.apply(jones.polarizer(th), jones.linear(0))), Math.cos(th) ** 2, 1e-15);
  // 0°, 45°, 90° polarisers pass 1/4 of light polarised along the first (1/8 of unpolarised light)
  const three = jones.system(jones.polarizer(0), jones.polarizer(Math.PI / 4), jones.polarizer(Math.PI / 2));
  close(jones.intensity(jones.apply(three, jones.linear(0))), 1 / 4, 1e-15);
  close(jones.intensity(jones.apply(jones.system(jones.polarizer(0), jones.polarizer(Math.PI / 2)), jones.linear(0))), 0, 1e-30);
  const circ = jones.stokes(jones.apply(jones.quarterWave(Math.PI / 4), jones.linear(0)));
  close(Math.abs(circ[3]), 1, 1e-15); close(circ[1], 0, 1e-15);
  const rot = jones.stokes(jones.apply(jones.halfWave(0.2), jones.linear(0)));
  close(rot[1], Math.cos(0.8), 1e-15); close(rot[2], Math.sin(0.8), 1e-15);
  close(jones.stokes(jones.circular(1))[3], 1, 1e-15);
  close(jones.stokes(jones.apply(jones.quarterWave(0), jones.linear(Math.PI / 4)))[3], 1, 1e-15, 'slow y axis delayed by π/2 gives [1, i]/√2');
  const v = jones.apply(jones.waveplate(0.7, 0.4), [{ re: 0.6, im: 0.1 }, { re: -0.2, im: 0.5 }]), S = jones.stokes(v);
  close(S[0] ** 2, S[1] ** 2 + S[2] ** 2 + S[3] ** 2, 1e-15); close(S[0], 0.66, 1e-15, 'waveplates are unitary');
  const r = jones.stokes(jones.apply(jones.rotator(0.5), jones.linear(0)));
  close(r[2], Math.sin(1), 1e-15);
});

// ---------------------------------------------------------------- FDTD

test('yeeDispersion: second-order approach to ω = ck; FDTD modes obey it exactly (TM and TE)', () => {
  const dx = 1e-7, dt = 0.5 * dx / (c * Math.SQRT2), err = h => yeeDispersion({ kx: 2e6, ky: 1e6, dx: h, dt: 0.5 * h / (c * Math.SQRT2) }) / (c * Math.hypot(2e6, 1e6)) - 1;
  close(err(dx) / err(dx / 2), 4, 0.01);
  assert.ok(err(dx) < 0, 'Yee waves lag');
  for (const mode of ['TM', 'TE']) {
    const nx = 32, ny = 16, sim = fdtd2D({ nx, ny, dx, dt, mode, boundary: 'periodic' });
    const kx = 2 * Math.PI * 3 / (nx * dx), ky = 2 * Math.PI * 2 / (ny * dx), F = mode === 'TM' ? sim.Ez : sim.Hz;
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) F[i + nx * j] = Math.cos(kx * i * dx + ky * j * dx);
    const cosw = Math.cos(yeeDispersion({ kx, ky, dx, dt }) * dt), hist = [F[5]];
    for (let s = 0; s < 60; s++) { sim.step(); hist.push(F[5]); }
    for (let s = 1; s < 60; s++) if (Math.abs(hist[s]) > 0.2) close((hist[s + 1] + hist[s - 1]) / (2 * hist[s]), cosw, 1e-12, mode);
  }
  assert.throws(() => fdtd2D({ nx: 8, ny: 8, dx, dt: 1.01 * dx / (c * Math.SQRT2) }), /Courant/);
});

test('fdtd2D: discrete energy exactly conserved in closed PEC and periodic boxes with a dielectric (TM, TE)', () => {
  const nx = 48, ny = 40, dx = 1e-7, epsR = new Float64Array(nx * ny).fill(1);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) if ((i - 30) ** 2 + (j - 20) ** 2 < 64) epsR[i + nx * j] = 4;
  for (const mode of ['TM', 'TE']) { for (const boundary of ['pec', 'periodic']) {
    const sim = fdtd2D({ nx, ny, dx, epsR, mode, boundary, courant: 0.9 }), F = mode === 'TM' ? sim.Ez : sim.Hz;
    for (let j = 1; j < ny - 1; j++) for (let i = 1; i < nx - 1; i++) F[i + nx * j] = Math.exp(-((i - 15) ** 2 + (j - 18) ** 2) / 10);
    const W0 = sim.energy();
    sim.step(800);
    close(sim.energy() / W0, 1, 1e-12, `${mode} ${boundary}`);
  } }
});

/** Quasi-1D run (periodic in y, Mur in x): soft line source, returns probe time series. */
function pulseRun({ mode, epsOf, nx, src, probes, steps, width = 20 }) {
  const dx = 1e-7, ny = 2, epsR = new Float64Array(nx * ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) epsR[i + nx * j] = epsOf(i);
  const tau = width * dx * Math.sqrt(epsOf(src)) / c, t0 = 5 * tau;
  const sim = fdtd2D({ nx, ny, dx, epsR, mode, boundary: { x: 'mur', y: 'periodic' }, sources: [{ i: src, j: 'all', f: t => Math.exp(-(((t - t0) / tau) ** 2)) }] });
  const F = mode === 'TM' ? sim.Ez : sim.Hz, out = probes.map(() => []);
  for (let s = 0; s < steps; s++) { sim.step(); probes.forEach((p, k) => out[k].push(F[p])); }
  return { out, dt: sim.dt, dx, sim };
}
const peak = (series, sign = 1, from = 0, to = series.length) => {
  let k = from;
  for (let i = from; i < to; i++) if (sign * series[i] > sign * series[k]) k = i;
  const [a, b, d] = [series[k - 1], series[k], series[k + 1]], off = (a - d) / (2 * (a - 2 * b + d));
  return { t: k + off, v: b - (a - d) * off / 4 };
};

test('fdtd2D: plane pulse travels at c/n; Mur walls absorb it', () => {
  const n = 1.5, { out, dt, dx, sim } = pulseRun({ mode: 'TM', epsOf: () => n * n, nx: 1600, src: 200, probes: [500, 1300], steps: 3700, width: 30 });
  const v = 800 * dx / ((peak(out[1]).t - peak(out[0]).t) * dt);
  close(v / (c / n), 1, 2e-3);
  const Wpulse = sim.energy();                                 // both pulses still inside
  sim.step(3000);                                              // ... and now both have left
  const W = sim.energy();
  assert.ok(W < 1e-5 * Wpulse, `energy left after the pulse exits: ${W / Wpulse}`);
});

test('fdtd2D: reflection and transmission at a dielectric step match Fresnel (r = −1/3 for E, +1/3 for H)', () => {
  const epsOf = i => (i >= 900 ? 4 : 1);
  for (const mode of ['TM', 'TE']) {
    const { out } = pulseRun({ mode, epsOf, nx: 1500, src: 250, probes: [550, 1150], steps: 3600 });
    const inc = peak(out[0], 1, 0, 1800), sgn = mode === 'TM' ? -1 : 1;
    const ref = peak(out[0], sgn, Math.round(inc.t) + 300, 3600), tr = peak(out[1], 1);
    close(ref.v / inc.v, sgn / 3, 0.01, `${mode} reflection`);
    close(tr.v / inc.v, mode === 'TM' ? 2 / 3 : 4 / 3, 0.02, `${mode} transmission`);
  }
});
