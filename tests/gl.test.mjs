import { test } from 'node:test';
import assert from 'node:assert/strict';
import { perspective, lookAt, multiply, invert, orbitCamera, colormap, normalize, COLORMAPS } from '../ui/gl.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const apply = (m, p) => [0, 1, 2, 3].map(r => m[r] * p[0] + m[4 + r] * p[1] + m[8 + r] * p[2] + m[12 + r] * (p[3] ?? 1));

test('camera maths: lookAt places the target ahead, perspective maps the frustum, invert is exact', () => {
  const v = lookAt([0, 0, 5], [0, 0, 0], [0, 1, 0]);
  apply(v, [0, 0, 0]).forEach((x, i) => close(x, [0, 0, -5, 1][i], 1e-6));
  const P = perspective(Math.PI / 2, 2, 1, 100);
  const near = apply(P, [0, 0, -1]), far = apply(P, [0, 0, -100]);
  close(near[2] / near[3], -1, 1e-6, 'near plane → −1'); close(far[2] / far[3], 1, 1e-5, 'far plane → +1');
  const edge = apply(P, [2, 1, -1]);
  close(edge[0] / edge[3], 1, 1e-6, 'right edge at 45°·aspect'); close(edge[1] / edge[3], 1, 1e-6, 'top edge at 45°');
  const M = multiply(P, v), I = multiply(M, invert(M));
  for (let i = 0; i < 16; i++) close(I[i], i % 5 === 0 ? 1 : 0, 1e-5);
});

test('orbitCamera: target projects to the centre; z is up; the right-hand side appears on the right', () => {
  const cam = orbitCamera({ target: [1, 2, 3], distance: 10, yaw: 0, pitch: 0 });
  const c = cam.project([1, 2, 3], 800, 600);
  close(c.x, 400, 1e-3); close(c.y, 300, 1e-3); close(c.depth, 10, 1e-4);
  assert.ok(cam.project([1, 2, 4], 800, 600).y < 300, '+z is up on screen');
  // camera on +x looking back at the target: +y is to the viewer's right
  assert.ok(cam.project([1, 3, 3], 800, 600).x > 400);
  close(Math.hypot(...cam.eye().map((e, i) => e - [1, 2, 3][i])), 10, 1e-12);
});

test('colour maps: published endpoints, monotonic lightness for sequential maps, neutral midpoint for coolwarm', () => {
  const h = s => [1, 3, 5].map(i => parseInt(s.slice(i, i + 2), 16) / 255);
  colormap('viridis', 0).forEach((x, i) => close(x, h('#440154')[i], 1e-12));
  colormap('viridis', 1).forEach((x, i) => close(x, h('#fde725')[i], 1e-12));
  colormap('inferno', 2).forEach((x, i) => close(x, h('#fcffa4')[i], 1e-12, 'clamped above 1'));
  const lum = c => { const lin = c.map(x => (x <= 0.04045 ? x / 12.92 : ((x + 0.055) / 1.055) ** 2.4)); return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]; };
  for (const name of ['viridis', 'inferno', 'magma']) {
    let prev = -1;
    for (let t = 0; t <= 1.0001; t += 0.01) { const L = lum(colormap(name, t)); assert.ok(L >= prev - 1e-9, `${name} lightness rises at ${t}`); prev = L; }
  }
  const mid = colormap('coolwarm', 0.5);
  assert.ok(Math.max(...mid) - Math.min(...mid) < 0.01, 'coolwarm midpoint is grey');
  assert.deepEqual(COLORMAPS, ['viridis', 'inferno', 'magma', 'coolwarm']);
  close(normalize(100, [1, 1e4], true), 0.5, 1e-15); close(normalize(3, [2, 4]), 0.5, 1e-15);
  assert.throws(() => colormap('rainbow', 0.5), /unknown/);
});
