import { test } from 'node:test';
import assert from 'node:assert/strict';
import { laplacianND, gradientND, divergenceND, fftND, solvePoissonPeriodic, diffusionPeriodic, schrodingerND } from '../core/fields.mjs';
import { conjugateGradient } from '../core/linalg.mjs';

const close = (a, b, tol = 1e-12) => assert.ok(Math.abs(a - b) <= tol, `${a} vs ${b} (tolerance ${tol})`);
const vectorClose = (a, b, tol = 1e-12) => { assert.equal(a.length, b.length); a.forEach((value, i) => close(value, b[i], tol)); };
const mode = (nx, ny, mx = 1, my = 1) => Float64Array.from({ length: nx * ny }, (_, index) => Math.cos(2 * Math.PI * (mx * (index % nx) / nx + my * Math.floor(index / nx) / ny)));

test('ND FFT matches direct 3D DFT, round trips and obeys Parseval', () => {
  const shape = [4, 2, 4], n = 32;
  const originalRe = Float64Array.from({ length: n }, (_, i) => Math.sin(0.7 * i) + 0.2 * i);
  const originalIm = Float64Array.from({ length: n }, (_, i) => Math.cos(0.3 * i));
  const re = originalRe.slice(), im = originalIm.slice();
  fftND(re, im, shape);
  for (let index = 0; index < n; index++) {
    const kx = index % 4, ky = Math.floor(index / 4) % 2, kz = Math.floor(index / 8);
    let expectedRe = 0, expectedIm = 0;
    for (let j = 0; j < n; j++) {
      const angle = -2 * Math.PI * (kx * (j % 4) / 4 + ky * (Math.floor(j / 4) % 2) / 2 + kz * Math.floor(j / 8) / 4);
      expectedRe += originalRe[j] * Math.cos(angle) - originalIm[j] * Math.sin(angle);
      expectedIm += originalRe[j] * Math.sin(angle) + originalIm[j] * Math.cos(angle);
    }
    close(re[index], expectedRe, 2e-13); close(im[index], expectedIm, 2e-13);
  }
  const energy = (a, b) => a.reduce((sum, value, i) => sum + value * value + b[i] * b[i], 0);
  close(energy(re, im) / n, energy(originalRe, originalIm), 1e-12);
  fftND(re, im, shape, true); vectorClose(re, originalRe, 5e-15); vectorClose(im, originalIm, 5e-15);
  assert.throws(() => fftND(re, im, [3, 2]), /power of two/);
  assert.throws(() => fftND(re, re, shape), /overlap/);
});

test('anisotropic ND derivatives match polynomial identities and known discrete Fourier symbols', () => {
  const nx = 9, ny = 7, spacing = [0.2, 0.3], shape = [nx, ny];
  const quadratic = Float64Array.from({ length: nx * ny }, (_, i) => (i % nx * spacing[0]) ** 2 + 3 * (Math.floor(i / nx) * spacing[1]) ** 2);
  const options = { shape, spacing, bc: 'dirichlet' };
  const L = laplacianND(quadratic, options), [gx, gy] = gradientND(quadratic, options);
  for (let j = 1; j < ny - 1; j++) { for (let i = 1; i < nx - 1; i++) {
    const index = i + nx * j;
    close(L[index], 8, 1e-12); close(gx[index], 2 * i * spacing[0]); close(gy[index], 6 * j * spacing[1]);
  } }
  const wave = mode(16, 8, 2, 1), grid = { shape: [16, 8], spacing: [0.1, 0.2] };
  const lambda = 4 * Math.sin(2 * Math.PI / 16) ** 2 / 0.1 ** 2 + 4 * Math.sin(Math.PI / 8) ** 2 / 0.2 ** 2;
  vectorClose(laplacianND(wave, grid), Float64Array.from(wave, value => -lambda * value), 2e-12);
  const components = gradientND(wave, grid), divGrad = divergenceND(components, grid);
  const widerSymbol = Math.sin(4 * Math.PI / 16) ** 2 / 0.1 ** 2 + Math.sin(2 * Math.PI / 8) ** 2 / 0.2 ** 2;
  vectorClose(divGrad, Float64Array.from(wave, value => -widerSymbol * value), 1e-12);
});

test('3D zero-flux stencil conserves the cell sum; mixed boundaries and alias guards work', () => {
  const shape = [4, 3, 5], u = Float64Array.from({ length: 60 }, (_, i) => Math.sin(0.2 * i) + 2);
  close(laplacianND(u, { shape, spacing: [1, 0.3, 2], bc: 'neumann' }).reduce((s, v) => s + v, 0), 0, 1e-13);
  const constant = new Float64Array(12).fill(1);
  vectorClose(laplacianND(constant, { shape: [4, 3], bc: ['periodic', 'neumann'] }), new Float64Array(12), 0);
  assert.throws(() => laplacianND(u, { shape }, u), /overlap/);
  assert.throws(() => laplacianND(u, { shape, spacing: [1, 2] }), /spacing/);
  assert.throws(() => laplacianND(u, { shape, bc: 'invalid' }), /boundary/);
  assert.throws(() => divergenceND([u], { shape }), /component/);
});

test('Poisson reproduces continuous Fourier modes and exact discrete residuals with the documented sign', () => {
  const nx = 32, ny = 16, spacing = [2 / nx, 3 / ny], shape = [nx, ny], u = mode(nx, ny, 2, 1);
  const lambda = (2 * Math.PI) ** 2 + (2 * Math.PI / 3) ** 2;
  const f = Float64Array.from(u, value => lambda * value);
  const spectral = solvePoissonPeriodic(f, { shape, spacing });
  vectorClose(spectral.u, u, 3e-14);
  const discreteSource = Float64Array.from(laplacianND(u, { shape, spacing }), value => -value);
  const discrete = solvePoissonPeriodic(discreteSource, { shape, spacing, discretization: 'finite-difference' });
  vectorClose(discrete.u, u, 3e-14);
  vectorClose(laplacianND(discrete.u, { shape, spacing }), Float64Array.from(discreteSource, value => -value), 2e-12);
  assert.throws(() => solvePoissonPeriodic(new Float64Array(8).fill(1), { shape: [8] }), /zero mean/);
  const screened = solvePoissonPeriodic(new Float64Array(8).fill(3), { shape: [8], screening: 2 });
  vectorClose(screened.u, new Float64Array(8).fill(0.75), 0);
});

test('matrix-free Dirichlet Poisson combines the shared field stencil and CG solver', () => {
  const shape = [15, 11], spacing = [0.1, 0.2], n = shape[0] * shape[1];
  const exact = Float64Array.from({ length: n }, (_, index) => Math.sin(Math.PI * ((index % 15) + 1) / 16) * Math.sin(Math.PI * (Math.floor(index / 15) + 1) / 12));
  const lambda = 4 * Math.sin(Math.PI / 32) ** 2 / 0.1 ** 2 + 4 * Math.sin(Math.PI / 24) ** 2 / 0.2 ** 2;
  const b = Float64Array.from(exact, value => lambda * value);
  const operator = (x, out) => { laplacianND(x, { shape, spacing, bc: 'dirichlet' }, out); for (let i = 0; i < n; i++) out[i] = -out[i]; };
  const result = conjugateGradient(operator, b, { rtol: 1e-11 });
  assert.equal(result.converged, true); vectorClose(result.x, exact, 2e-13);
});

test('finite-difference Poisson converges spatially at second order against a continuum mode', () => {
  const error = n => {
    const shape = [n, n], spacing = [1 / n, 1 / n], exact = mode(n, n);
    const source = Float64Array.from(exact, value => 8 * Math.PI ** 2 * value);
    const { u } = solvePoissonPeriodic(source, { shape, spacing, discretization: 'finite-difference' });
    return Math.max(...Float64Array.from(u, (value, i) => Math.abs(value - exact[i])));
  };
  const coarse = error(16), medium = error(32), fine = error(64);
  close(Math.log2(coarse / medium), 2, 0.02); close(Math.log2(medium / fine), 2, 0.02);
});

test('periodic diffusion has the exact decay of 2D modes, conserves mass, and is frame-independent', () => {
  const shape = [32, 16], spacing = [1 / 32, 2 / 16], u0 = Float64Array.from(mode(...shape), value => 2 + value);
  const config = { u0, shape, spacing, D: 0.1, dt: 0.01 };
  const direct = diffusionPeriodic(config), framed = diffusionPeriodic(config), total0 = direct.total();
  direct.advanceTo(0.2);
  for (const time of [0.01, 0.07, 0.13, 0.2]) framed.advanceTo(time);
  assert.deepEqual(direct.u, framed.u);
  const decay = Math.exp(-0.1 * ((2 * Math.PI) ** 2 + Math.PI ** 2) * 0.2);
  vectorClose(direct.u, Float64Array.from(u0, value => 2 + (value - 2) * decay), 4e-15);
  close(direct.total(), total0, 1e-14);
  assert.throws(() => direct.advanceTo(NaN), /target/);
  assert.throws(() => diffusionPeriodic({ ...config, t0: 1e20 }).step(), /underflow/);
});

test('ND Schrodinger plane wave has the exact phase including potential and conserved norm', () => {
  const shape = [8, 16], spacing = [0.25, 0.125], n = 128, mass = 2, hbar = 0.7, dt = 0.003;
  const phase = Float64Array.from({ length: n }, (_, i) => 2 * Math.PI * ((i % 8) / 8 + 2 * Math.floor(i / 8) / 16));
  const re0 = Float64Array.from(phase, value => Math.cos(value)), im0 = Float64Array.from(phase, value => Math.sin(value)), V = new Float64Array(n).fill(1.3);
  const s = schrodingerND({ re0, im0, V, shape, spacing, dt, mass, hbar });
  const N0 = s.norm(); for (let i = 0; i < 100; i++) s.step();
  const omega = hbar / (2 * mass) * (Math.PI ** 2 + (2 * Math.PI) ** 2) + 1.3 / hbar;
  vectorClose(s.re, Float64Array.from(phase, value => Math.cos(value - omega * s.t)), 2e-13);
  vectorClose(s.im, Float64Array.from(phase, value => Math.sin(value - omega * s.t)), 2e-13);
  close(s.norm(), N0, 5e-13);
  assert.deepEqual(re0, Float64Array.from(phase, value => Math.cos(value)));
  vectorClose(s.probability(), new Float64Array(n).fill(1), 3e-13);
});

test('2D split-operator has second-order time convergence for a noncommuting potential', () => {
  const shape = [16, 16], n = 256, spacing = [0.25, 0.25];
  const re0 = Float64Array.from({ length: n }, (_, i) => Math.exp(-0.6 * ((i % 16 - 8) ** 2 + (Math.floor(i / 16) - 8) ** 2)));
  const V = Float64Array.from({ length: n }, (_, i) => 0.7 * Math.cos(2 * Math.PI * (i % 16) / 16) + 0.3 * Math.sin(2 * Math.PI * Math.floor(i / 16) / 16));
  const at = dt => { const s = schrodingerND({ re0, V, shape, spacing, dt }); for (let i = 0; i < Math.round(0.16 / dt); i++) s.step(); return s; };
  const a = at(0.02), b = at(0.01), c = at(0.005);
  const difference = (u, v) => Math.hypot(...Float64Array.from(u.re, (value, i) => value - v.re[i]), ...Float64Array.from(u.im, (value, i) => value - v.im[i]));
  close(Math.log2(difference(a, b) / difference(b, c)), 2, 0.08);
});

test('field propagation rejects invalid scales and reports a real work limit', () => {
  const re0 = new Float64Array(8).fill(1), shape = [8];
  assert.throws(() => schrodingerND({ re0, shape, dt: 0.1, hbar: 0 }), /positive/);
  assert.throws(() => schrodingerND({ re0, V: [1], shape, dt: 0.1 }), /length/);
  assert.throws(() => schrodingerND({ re0, shape, dt: 0.1, t0: 1e20 }).step(), /underflow/);
  const s = schrodingerND({ re0, shape, dt: 0.01 });
  assert.equal(s.advanceTo(1, 2), false); close(s.t, 0.02); assert.equal(s.steps, 2);
  assert.throws(() => s.advanceTo(Infinity), /target/);
});
