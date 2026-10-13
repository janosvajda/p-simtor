import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  planckLambda, planckNu, rayleighJeans, WIEN_B, WIEN_X_NU, SIGMA_SB, exitance, bandRadiance,
  MEAN_PHOTON_ENERGY_KT, spectrumToXYZ, chromaticity, blackbodyColor, c2,
} from '../core/radiation.mjs';
import { h, kB, c } from '../core/constants.mjs';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('Wien displacement constant derived from the exact root matches CODATA 2018', () => {
  assert.ok(rel(WIEN_B, 2.897771955e-3) < 1e-9, `${WIEN_B}`);
});

test('Stefan–Boltzmann constant from first principles matches CODATA 2018', () => {
  assert.ok(rel(SIGMA_SB, 5.670374419e-8) < 1e-9, `${SIGMA_SB}`);
});

test('integrating Planck’s law reproduces σT⁴/π', () => {
  for (const T of [300, 5772, 30000]) {
    const peak = WIEN_B / T;
    assert.ok(rel(bandRadiance(T, peak / 100, peak * 1e4), exitance(T) / Math.PI) < 1e-9, `T=${T}`);
  }
});

test('B_λ peaks at λ_max = b/T and B_ν peaks at ν = x k T / h', () => {
  const T = 5772, lm = WIEN_B / T;
  assert.ok(planckLambda(lm, T) > planckLambda(lm * 1.001, T) && planckLambda(lm, T) > planckLambda(lm * 0.999, T));
  const nm = WIEN_X_NU * kB * T / h;
  assert.ok(planckNu(nm, T) > planckNu(nm * 1.001, T) && planckNu(nm, T) > planckNu(nm * 0.999, T));
  // ν_max ≠ c/λ_max: the two peaks describe different densities
  assert.ok(rel(c / nm, lm) > 0.5);
});

test('Rayleigh–Jeans is the long-wavelength limit of Planck', () => {
  assert.ok(rel(planckLambda(1, 300), rayleighJeans(1, 300)) < 1e-4);
  assert.ok(planckLambda(1e-7, 300) < 1e-30 * rayleighJeans(1e-7, 300), 'no ultraviolet catastrophe');
});

test('mean photon energy is 2.701 kT', () => {
  assert.ok(Math.abs(MEAN_PHOTON_ENERGY_KT - 2.7011780) < 1e-7);
});

test('colour: Planckian locus points (limited by the CMF fit, ~0.002)', () => {
  // CIE illuminant A: Planckian at 2856 K with c2 = 1.435e-2 m K → x = 0.44757, y = 0.40745
  const TA = 2856 * c2 / 1.435e-2;
  const [xA, yA] = chromaticity(spectrumToXYZ(nm => planckLambda(nm * 1e-9, TA)));
  assert.ok(Math.abs(xA - 0.44757) < 0.003 && Math.abs(yA - 0.40745) < 0.003, `A: ${xA}, ${yA}`);
  // Planckian locus at 6500 K: (0.3135, 0.3236)
  const [x65, y65] = blackbodyColor(6500).xy;
  assert.ok(Math.abs(x65 - 0.3135) < 0.002 && Math.abs(y65 - 0.3236) < 0.002, `6500 K: ${x65}, ${y65}`);
  // hot bodies are blue, cool bodies are red
  const [r1, , b1] = blackbodyColor(2000).rgb, [r2, , b2] = blackbodyColor(20000).rgb;
  assert.ok(r1 > b1 && b2 > r2);
});

test('Planckian locus from 1000 K to 10000 K (CIE tables), and cold bodies tend to the red end', () => {
  const locus = { 1000: [0.6528, 0.3444, 0.009], 1500: [0.5857, 0.3931, 0.004], 2000: [0.5267, 0.4133, 0.0025],
    2856: [0.4476, 0.4074, 0.0012], 4000: [0.3805, 0.3768, 0.0006], 6500: [0.3135, 0.3236, 0.0003], 10000: [0.2807, 0.2884, 0.0003] };
  for (const [T, [x, y, tol]] of Object.entries(locus)) {
    const [cx, cy] = blackbodyColor(+T).xy;
    assert.ok(Math.abs(cx - x) < tol && Math.abs(cy - y) < tol, `${T} K: (${cx}, ${cy}) vs (${x}, ${y})`);
  }
  // T → 0: spectral locus red end (CIE 1931: x = 0.7347, y = 0.2653 for λ ≥ 700 nm)
  for (const T of [2.7255, 1, 0.01]) {
    const [cx, cy] = blackbodyColor(T).xy;
    assert.ok(Math.abs(cx - 0.7347) < 1e-4 && Math.abs(cy - 0.2653) < 1e-4, `${T} K`);
  }
});

test('blackbody colour is finite for every temperature (no underflow to NaN)', () => {
  for (const T of [1e-3, 0.5, 2.7255, 10, 310, 1e4, 1e6, 1e8]) {
    const { rgb, xy } = blackbodyColor(T);
    assert.ok(rgb.every(Number.isFinite) && xy.every(Number.isFinite), `${T} K: ${rgb}`);
    assert.ok(rgb.every(v => v >= 0 && v <= 255));
  }
  assert.deepEqual(blackbodyColor(2.7255).rgb.slice(0, 2), [255, 0], 'deep red, not pink');
});
