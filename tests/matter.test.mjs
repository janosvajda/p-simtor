import { test } from 'node:test';
import assert from 'node:assert/strict';
import { MATERIALS, betheStoppingPower, csdaRange, densityEffect } from '../core/matter.mjs';

const MMU = 105.6583755;

test('Bethe minimum ionisation reproduces the PDG table for every material (within 0.2 %)', () => {
  for (const [key, mat] of Object.entries(MATERIALS)) {
    let best = Infinity;
    for (let l = -1; l <= 3; l += 0.002) best = Math.min(best, betheStoppingPower(10 ** l, { mass: MMU, material: key }));
    assert.ok(Math.abs(best / mat.minIon - 1) < 2e-3, `${key}: ${best} vs ${mat.minIon}`);
  }
});

test('stopping power scales as z², depends on βγ only (not on mass) up to the small Wmax term; δ is continuous', () => {
  const s1 = betheStoppingPower(1, { mass: MMU, material: 'water' }), s2 = betheStoppingPower(1, { mass: MMU, z: 2, material: 'water' });
  assert.ok(Math.abs(s2 / s1 - 4) < 1e-12);
  const p = betheStoppingPower(1, { mass: 938.272, material: 'water' });
  assert.ok(Math.abs(p / s1 - 1) < 0.01);
  for (const mat of Object.values(MATERIALS)) {for (const x of [mat.x0, mat.x1]) {
    const bg = 10 ** x;
    assert.ok(Math.abs(densityEffect(bg * (1 + 1e-9), mat) - densityEffect(bg * (1 - 1e-9), mat)) < 1e-2, `${mat.name} at x = ${x}`);
  }}
});

test('CSDA range: monotonic; a 4.12 MeV muon (from pion decay at rest) stops within about 1 cm of liquid hydrogen', () => {
  const R = T => csdaRange(T, { mass: MMU, material: 'liquidHydrogen' }) / MATERIALS.liquidHydrogen.rho;
  assert.ok(R(2) < R(4.12) && R(4.12) < R(10));
  const r = R(4.11987);
  // PDG muon table: 0.884 cm; the uncorrected Bethe formula is ~6 % long at these energies (documented)
  assert.ok(r > 0.88 && r < 0.96, `range ${r} cm`);
  assert.throws(() => betheStoppingPower(1, { mass: MMU, material: 'unobtainium' }), /unknown material/);
});
