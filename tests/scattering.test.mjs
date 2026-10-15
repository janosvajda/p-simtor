import { test } from 'node:test';
import assert from 'node:assert/strict';
import { transmission1D, thermalAverage, packetTransmission } from '../core/pde.mjs';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
// Eckart barrier V0/cosh²(x/a): T = sinh²(πka) / (sinh²(πka) + cosh²((π/2)√(8mV0a²/ħ² − 1)))  (ħ = m = 1)
const eckartT = (E, V0, a) => { const k = Math.sqrt(2 * E), s = Math.sinh(Math.PI * k * a) ** 2; return s / (s + Math.cosh((Math.PI / 2) * Math.sqrt(8 * V0 * a * a - 1)) ** 2); };

test('transmission1D reproduces the exact Eckart transmission from deep tunnelling to above the barrier', () => {
  const V0 = 5, a = 1, V = x => V0 / Math.cosh(x / a) ** 2;
  for (const E of [0.5, 2, 4.9, 5.5, 9]) {
    const { T, R } = transmission1D({ V, E, x0: -25, x1: 25 });
    assert.ok(rel(T, eckartT(E, V0, a)) < 1e-6, `E = ${E}: ${T} vs ${eckartT(E, V0, a)}`);
    assert.ok(Math.abs(T + R - 1) < 1e-8);
  }
});

test('transmission1D: rectangular barrier vs 1/(1 + V0² sinh²(κw)/(4E(V0 − E))), and the heavier particle tunnels less', () => {
  const V0 = 3, w = 2, V = x => (Math.abs(x) < w / 2 ? V0 : 0);
  for (const E of [0.4, 1.5, 2.8]) {
    const kap = Math.sqrt(2 * (V0 - E)), exact = 1 / (1 + V0 * V0 * Math.sinh(kap * w) ** 2 / (4 * E * (V0 - E)));
    assert.ok(rel(transmission1D({ V, E, x0: -3, x1: 3 }).T, exact) < 1e-5, `E = ${E}`);
  }
  assert.ok(transmission1D({ V, E: 1, x0: -3, x1: 3, mass: 2 }).T < transmission1D({ V, E: 1, x0: -3, x1: 3 }).T / 10);
});

test('thermalAverage: ∫e^{−E/kT} = kT; a step T(E) = θ(E − V) gives e^{−V/kT} (classical rate)', () => {
  assert.ok(Math.abs(thermalAverage(() => 1, 0.7) - 1) < 1e-12);
  const Vb = 2, kT = 0.5;
  assert.ok(rel(thermalAverage(E => (E > Vb ? 1 : 0), kT, { panels: 200 }), Math.exp(-Vb / kT)) < 2e-3);
});

test('packetTransmission: constant T gives T; a narrow packet gives T(k0)', () => {
  assert.ok(Math.abs(packetTransmission(() => 0.3, 5, 0.4) - 0.3) < 1e-10);
  const T = k => 1 / (1 + Math.exp(-(k - 4)));
  assert.ok(Math.abs(packetTransmission(T, 4.5, 1e-3) - T(4.5)) < 1e-6);
});
