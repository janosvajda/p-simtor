// Test fixture for ui/gl.mjs: a static ring of points, a 2D field and a 3D Gaussian blob.
// Exists only so `npm run test:browser` can check the GPU renderers pixel by pixel.

export const PARAMS = [{ key: 'n', label: 'Points on the ring', min: 3, max: 20000, value: 2000, integer: true, log: true }];

export function createSim({ n }) {
  const ring = new Float32Array(3 * n), vals = new Float32Array(n);
  for (let i = 0; i < n; i++) { const a = 2 * Math.PI * i / n; ring.set([Math.cos(a), Math.sin(a), 0], 3 * i); vals[i] = i / n; }
  const N = 32, blob = new Float32Array(N ** 3);
  for (let k = 0; k < N; k++) {for (let j = 0; j < N; j++) {for (let i = 0; i < N; i++) {
    const x = (i + 0.5) / N * 2 - 1, y = (j + 0.5) / N * 2 - 1, z = (k + 0.5) / N * 2 - 1;
    blob[i + N * (j + N * k)] = Math.exp(-(x * x + y * y + z * z) / 0.1);
  }}}
  return { t: 0, ring, vals, blob, N, step() { this.t += 0.01; }, diagnostics: () => [], state: () => ({ n }) };
}

export function verify() {
  const { blob, N } = createSim({ n: 10 });
  const peak = Math.max(...blob);
  return { summary: 'GPU renderer fixture', checks: [{ label: 'Gaussian blob peaks near the centre cell (exp(−3/(0.1·N²)) for N = 32, float32 storage)', value: peak.toFixed(6), pass: Math.abs(peak - Math.exp(-3 / (0.1 * N * N))) < 1e-7 }] };
}
