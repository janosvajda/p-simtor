// The starter (psim new --starter) must itself be a correct, verified page source.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { PARAMS, createSim, verify } from '../templates/starter/physics.mjs';
import { validateModuleReport } from '../core/verification.mjs';

test('starter physics verifies against the exact decay law', () => {
  const r = validateModuleReport(verify());
  assert.ok(r.checks.every(c => c.pass), JSON.stringify(r.checks));
  const sim = createSim(Object.fromEntries(PARAMS.map(p => [p.key, p.value])));
  for (let i = 0; i < 200; i++) sim.step();
  assert.ok(Math.abs(sim.N / sim.exact(sim.t) - 1) < 1e-9);
});

test('starter page has every required part', () => {
  const html = readFileSync(new URL('../templates/starter/page.html', import.meta.url), 'utf8');
  for (const part of ['intro:', 'actions:', 'reveal:', 'timeline:', 'seconds:', 'legend:', 'view.steady(', 'bench', 'createPlot(', 'name="generator" content="p-simtor"', 'createInstrument({']) assert.ok(html.includes(part), part);
  assert.match(readFileSync(new URL('../templates/starter/physics.mjs', import.meta.url), 'utf8'), /ONE IDEA[\s\S]*Learn:[\s\S]*Picture:[\s\S]*Action:/);
});

test('starter view draws the atoms themselves: the sample follows the law it illustrates', () => {
  const sim = createSim({ halfLife: 10 });
  for (const [k, expect] of [[1, 200], [2, 100], [3, 50]]) assert.ok(Math.abs(sim.intact(10 * k) - expect) <= 0.03 * 400, `${k} half-lives: ${sim.intact(10 * k)}`);
});
