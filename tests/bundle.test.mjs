import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import vm from 'node:vm';
import { bundleHtml, bundleModules, auditHtml } from '../scripts/lib/bundle.mjs';

const SKILL_DIR = join(dirname(fileURLToPath(import.meta.url)), '..');

test('bundled modules run and give the same physics as the source modules', async () => {
  const code = bundleModules(`
import { integrateTo } from 'psim:core/nbody.mjs';
import { THREE_BODY_PRESETS as P } from 'psim:core/presets.mjs';
result = integrateTo(P['figure-eight'].bodies, { integrator: 'logh4' }, 1).x[0];
`, SKILL_DIR);
  assert.ok(!/^\s*import\s/m.test(code), 'no import statements left');
  const ctx = { result: null };
  vm.runInNewContext(code, ctx);
  const { integrateTo } = await import('../core/nbody.mjs');
  const { THREE_BODY_PRESETS } = await import('../core/presets.mjs');
  assert.equal(ctx.result, integrateTo(THREE_BODY_PRESETS['figure-eight'].bodies, { integrator: 'logh4' }, 1).x[0]);
});

test('gravity template bundles into a self-contained page', () => {
  const html = bundleHtml(readFileSync(join(SKILL_DIR, 'templates/gravity.html'), 'utf8'), {
    skillDir: SKILL_DIR, config: { title: 'x' }, verification: { a: 1 },
  });
  assert.deepEqual(auditHtml(html), []);
  assert.ok(html.includes('/* psim:ui/instrument.css */'));
  assert.ok(!/<link[^>]+psim:/.test(html));
});

test('JSON payloads cannot close the script tag', () => {
  const html = bundleHtml('<script type="application/json" id="psim-config">null</script>', {
    skillDir: SKILL_DIR, config: { note: '</script><script>alert(1)</script>' },
  });
  assert.equal(html.match(/<\/script>/g).length, 1);
});

test('audit flags external resources and missing verification', () => {
  const problems = auditHtml('<script src="https://cdn.example/x.js"></script><script type="application/json" id="psim-verification">null</script>');
  assert.equal(problems.length, 2);
});
