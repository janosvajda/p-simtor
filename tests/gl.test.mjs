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
  assert.deepEqual(COLORMAPS, ['viridis', 'inferno', 'magma', 'coolwarm', 'phase']);
  close(normalize(100, [1, 1e4], true), 0.5, 1e-15); close(normalize(3, [2, 4]), 0.5, 1e-15);
  assert.throws(() => colormap('rainbow', 0.5), /unknown/);
});

test('relative camera matrices give the same projection as absolute ones; phase map is cyclic with constant lightness', () => {
  const cam = orbitCamera({ target: [3, -2, 1], distance: 7, yaw: 1.1, pitch: 0.3 });
  const abs = cam.matrices(1.5).viewProj, rel = cam.matrices(1.5, { relative: true }).viewProj;
  for (const p of [[3, -2, 1], [4, 0, -1], [-2, 5, 3]]) {
    const a = apply(abs, p), r = apply(rel, [p[0] - 3, p[1] + 2, p[2] - 1]);
    close(a[0] / a[3], r[0] / r[3], 1e-5); close(a[1] / a[3], r[1] / r[3], 1e-5);
  }
  colormap('phase', 0).forEach((x, i) => close(x, colormap('phase', 1)[i], 1e-12, 'cyclic'));
  const toLin = x => (x <= 0.04045 ? x / 12.92 : ((x + 0.055) / 1.055) ** 2.4);
  const okL = c => { const [r, g, b] = c.map(toLin); const l = Math.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b), m = Math.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b), s = Math.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b); return 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s; };
  for (let t = 0; t < 1; t += 0.05) close(okL(colormap('phase', t)), 0.72, 0.02, `OKLab lightness at ${t}`);
  assert.ok(COLORMAPS.includes('phase'));
});

test('orbitCamera.fit frames a box on wide and tall screens; limits clamp the distance', () => {
  const box = [[-2, -1, -0.5], [2, 1, 0.5]], r = Math.hypot(4, 2, 1) / 2;
  for (const aspect of [16 / 9, 0.46]) {
    const cam = orbitCamera({ yaw: 0.7, pitch: 0.4 }).fit(box, aspect, 1);
    // every corner projects inside the screen
    for (const x of [-2, 2]) {for (const y of [-1, 1]) {for (const z of [-0.5, 0.5]) {
      const p = cam.project([x, y, z], 1000 * aspect, 1000);
      assert.ok(p.x >= -1 && p.x <= 1000 * aspect + 1 && p.y >= -1 && p.y <= 1001, `corner inside at aspect ${aspect}`);
    }}}
    assert.ok(Math.abs(Math.sin(Math.min(cam.fov / 2, Math.atan(Math.tan(cam.fov / 2) * aspect))) * cam.distance - r) < 1e-9, 'bounding sphere touches the narrower edge');
  }
  const limited = orbitCamera({ minDistance: 10, maxDistance: 20 }).fit(box, 1);
  assert.equal(limited.distance, 10, 'clamped to the minimum distance');
});

test('mesh geometry and transforms: outward unit normals, exact radii, rotation order, composition', async () => {
  const { geometry, transform, compose } = await import('../ui/gl.mjs');
  const box = geometry.box(2, 3, 4);
  for (let f = 0; f < 6; f++) {                       // each face: 4 vertices, normal points away from the centre
    const i = 4 * f, n = [box.normals[3 * i], box.normals[3 * i + 1], box.normals[3 * i + 2]];
    const c = [0, 1, 2].map(k => (box.positions[3 * i + k] + box.positions[3 * (i + 2) + k]) / 2);
    close(Math.hypot(...n), 1, 1e-7); assert.ok(n[0] * c[0] + n[1] * c[1] + n[2] * c[2] > 0, `face ${f} outward`);
  }
  const sph = geometry.sphere(1.5, 12);
  for (let i = 0; i < sph.positions.length; i += 3) close(Math.hypot(sph.positions[i], sph.positions[i + 1], sph.positions[i + 2]), 1.5, 1e-6);
  const cyl = geometry.cylinder(0.5, 2, 16);
  assert.ok(Math.max(...cyl.positions.filter((_, i) => i % 3 === 2)) === 2 && Math.min(...cyl.positions.filter((_, i) => i % 3 === 2)) === 0);
  const T = transform({ translate: [1, 2, 3], rotate: [0, 0, Math.PI / 2], scale: 2 });
  const map = (m, p) => [0, 1, 2].map(r => m[r] * p[0] + m[4 + r] * p[1] + m[8 + r] * p[2] + m[12 + r]);
  map(T, [1, 0, 0]).forEach((v, i) => close(v, [1, 4, 3][i], 1e-12, 'x → y, doubled, then shifted'));
  const R = transform({ rotate: [Math.PI / 2, 0, Math.PI / 2] });  // about x first, then z
  map(R, [0, 1, 0]).forEach((v, i) => close(v, [0, 0, 1][i], 1e-12, 'y →(x-rot) z →(z-rot) z'));
  map(compose(transform({ translate: [5, 0, 0] }), T), [0, 0, 0]).forEach((v, i) => close(v, [6, 2, 3][i], 1e-12));
});
