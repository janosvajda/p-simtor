import { test } from 'node:test';
import assert from 'node:assert/strict';
import { phaseWinding, fftAxis, fftND } from '../core/fields.mjs';

test('phaseWinding finds a vortex (+1) and an antivortex (−1) at their positions, nothing elsewhere', () => {
  const nx = 64, ny = 48, re = new Float64Array(nx * ny), im = new Float64Array(nx * ny);
  const v1 = [20.3, 24.4], v2 = [44.6, 20.2];
  for (let j = 0; j < ny; j++) {for (let i = 0; i < nx; i++) {
    const th = Math.atan2(j - v1[1], i - v1[0]) - Math.atan2(j - v2[1], i - v2[0]);
    const amp = Math.exp(-(((i - 32) / 30) ** 2 + ((j - 24) / 22) ** 2));
    re[i + nx * j] = amp * Math.cos(th); im[i + nx * j] = amp * Math.sin(th);
  }}
  const { vortices } = phaseWinding(re, im, { shape: [nx, ny] });
  assert.equal(vortices.length, 2);
  const plus = vortices.find(v => v.charge === 1), minus = vortices.find(v => v.charge === -1);
  assert.ok(Math.hypot(plus.x - v1[0], plus.y - v1[1]) < 1 && Math.hypot(minus.x - v2[0], minus.y - v2[1]) < 1);
});

test('fftAxis along each axis composes to fftND', () => {
  const shape = [8, 4], n = 32, re = Float64Array.from({ length: n }, (_, i) => Math.sin(i * 0.7)), im = Float64Array.from({ length: n }, (_, i) => Math.cos(i * 0.3));
  const a = [Float64Array.from(re), Float64Array.from(im)], b = [Float64Array.from(re), Float64Array.from(im)];
  fftND(a[0], a[1], shape);
  fftAxis(b[0], b[1], shape, 0); fftAxis(b[0], b[1], shape, 1);
  for (let i = 0; i < n; i++) assert.ok(Math.abs(a[0][i] - b[0][i]) < 1e-12 && Math.abs(a[1][i] - b[1][i]) < 1e-12);
});
