// The recipe catalog must point at real, buildable examples, and typical requests must
// rank the right recipe first (that is what keeps agents fast and cheap).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const { recipes } = JSON.parse(readFileSync(join(ROOT, 'examples/catalog.json'), 'utf8'));
const psim = (...args) => execFileSync('node', [join(ROOT, 'scripts/psim.mjs'), ...args], { encoding: 'utf8' });

test('every compose recipe has physics.mjs, page.html and a built page', () => {
  for (const r of recipes.filter(x => x.kind === 'compose')) {
    for (const f of ['physics.mjs', 'page.html', `${r.id}.html`]) assert.ok(existsSync(join(ROOT, r.dir, f)), `${r.id}: ${f}`);
  }
});

test('every example prompt in the README ranks its own recipe first', () => {
  const readme = readFileSync(join(ROOT, 'README.md'), 'utf8');
  const rows = [...readme.matchAll(/^\| `\/p-simtor ([^`]+)` \| \[[^\]]+\]\(examples\/([\w-]+)\//gm)];
  assert.ok(rows.length >= 7, `found ${rows.length} README example rows`);
  for (const [, prompt, id] of rows) {
    const first = psim('recipes', ...prompt.split(' ')).match(/^ {2}([\w-]+) {2}\(match/m)?.[1];
    assert.equal(first, id, prompt);
  }
});

test('typical requests rank the right recipe first', () => {
  const cases = {
    'simulate a chaotic three-body gravitational system': 'three-body',
    'simulate three stars orbiting each other': 'three-body',
    'model planets orbiting the sun': 'solar-system',
    'simulate a proton moving through a uniform magnetic field': 'proton-magnetic-field',
    'simulate black body radiation': 'blackbody',
    'black hole event horizon with an accretion disk': 'black-hole',
    'gravitational lensing einstein ring': 'black-hole',
    'dark matter galaxy simulation': 'galaxy',
    'how gps works with relativity time dilation': 'gps-relativity',
  };
  for (const [request, expected] of Object.entries(cases)) {
    const first = psim('recipes', ...request.split(' ')).match(/^ {2}([\w-]+) {2}\(match/m)?.[1];
    assert.equal(first, expected, request);
  }
});
