import { test } from 'node:test';
import assert from 'node:assert/strict';
import { hydrogenRadial, complexYlm, realYlm, gaussLegendre, integrateToInfinity } from '../core/numerics.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const radial = (f) => integrateToInfinity(f, { panels: 40, n: 48 });

test('hydrogenRadial: closed forms, normalisation, ⟨r⟩ = (3n² − l(l+1))/2Z, orthogonality', () => {
  for (const r of [0.1, 1, 3.7]) {
    close(hydrogenRadial(1, 0, r), 2 * Math.exp(-r), 1e-15);
    close(hydrogenRadial(2, 1, r), r * Math.exp(-r / 2) / (2 * Math.sqrt(6)), 1e-15);
    close(hydrogenRadial(2, 0, r, { Z: 3 }), 2 * (1.5) ** 1.5 * (1 - 1.5 * r) * Math.exp(-1.5 * r), 1e-14, 'Z = 3');
  }
  for (const [n, l] of [[1, 0], [2, 1], [3, 0], [3, 2], [5, 3], [8, 1]]) {
    close(radial(r => (hydrogenRadial(n, l, r) * r) ** 2), 1, 1e-12, `norm ${n}${l}`);
    close(radial(r => hydrogenRadial(n, l, r) ** 2 * r ** 3), (3 * n * n - l * (l + 1)) / 2, 1e-10, `⟨r⟩ ${n}${l}`);
  }
  close(radial(r => hydrogenRadial(2, 0, r) * hydrogenRadial(5, 0, r) * r * r), 0, 1e-13, 'orthogonal');
});

test('complexYlm: standard forms, orthonormality on the sphere, relation to the real harmonics', () => {
  const th = 0.9, ph = 2.1;
  close(complexYlm(1, 0, th, ph).re, Math.sqrt(3 / (4 * Math.PI)) * Math.cos(th), 1e-15);
  const y11 = complexYlm(1, 1, th, ph), c = -Math.sqrt(3 / (8 * Math.PI)) * Math.sin(th);
  close(y11.re, c * Math.cos(ph), 1e-15); close(y11.im, c * Math.sin(ph), 1e-15);
  const y1m = complexYlm(1, -1, th, ph);
  close(y1m.re, -y11.re, 1e-15); close(y1m.im, y11.im, 1e-15, 'Y₁⁻¹ = −conj(Y₁¹)');
  close(realYlm(3, 2, th, ph), Math.SQRT2 * complexYlm(3, 2, th, ph).re, 1e-15, 'real = √2 (−1)^m Re Y');
  const inner = (a, b) => gaussLegendre(t => gaussLegendre(p => {
    const A = complexYlm(...a, t, p), B = complexYlm(...b, t, p);
    return (A.re * B.re + A.im * B.im) * Math.sin(t);
  }, 0, 2 * Math.PI, 32), 0, Math.PI, 32);
  close(inner([4, -3], [4, -3]), 1, 1e-13); close(inner([4, 2], [3, 2]), 0, 1e-13); close(inner([2, 1], [2, -1]), 0, 1e-13);
});
