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

test('README lists every gallery demo with a live link and never links to source-view pages', () => {
  const readme = readFileSync(join(ROOT, 'README.md'), 'utf8');
  const { demos } = JSON.parse(readFileSync(join(ROOT, 'docs', 'demos.json'), 'utf8'));
  for (const d of demos) assert.ok(readme.includes(`https://janosvajda.github.io/p-simtor/${d.id}.html`), `README links the live ${d.id} demo`);
  assert.equal(readme.match(/\]\(examples\/[^)]*\.html\)/g), null, 'no links to examples/*.html (GitHub shows them as source)');
});

test('recipe prompts rank their own recipe first', () => {
  const cases = {
    'model the planets orbiting the Sun': 'solar-system', 'simulate three stars orbiting each other': 'three-body',
    'show why galaxies need dark matter': 'galaxy', 'how does relativity affect GPS clocks?': 'gps-relativity',
    'simulate a black hole event horizon': 'black-hole', 'simulate black body radiation': 'blackbody',
    'simulate a proton moving through a magnetic field': 'proton-magnetic-field',
  };
  for (const [prompt, id] of Object.entries(cases)) {
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
