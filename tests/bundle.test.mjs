import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, mkdtempSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
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

test('bundled modules export exactly what Node sees, for every core and ui module', async () => {
  // The bundle must never silently drop an export (it once missed `export const R0 = 8, V0 = 220`).
  const { readdirSync } = await import('node:fs');
  const { register } = await import('node:module');
  register(join(SKILL_DIR, 'scripts/lib/psim-hooks.mjs'), import.meta.url);   // psim: imports in example modules
  const modules = [
    ...['core', 'ui'].flatMap(dir => readdirSync(join(SKILL_DIR, dir)).filter(n => n.endsWith('.mjs')).map(f => `${dir}/${f}`)),
    ...readdirSync(join(SKILL_DIR, 'examples')).filter(d => !d.includes('.')).map(d => `examples/${d}/physics.mjs`)
      .filter(k => { try { readFileSync(join(SKILL_DIR, k)); return true; } catch { return false; } }),
  ];
  {
    for (const key of modules) {
      const code = key.startsWith('examples/')
        ? bundleModules(`import * as M from './${key.split('/').pop()}';\nresult = Object.keys(M).sort();`, SKILL_DIR, join(SKILL_DIR, dirname(key)))
        : bundleModules(`import * as M from 'psim:${key}';\nresult = Object.keys(M).sort();`, SKILL_DIR);
      const ctx = { result: null, document: undefined, devicePixelRatio: 1 };
      vm.runInNewContext(code, ctx);
      const real = Object.keys(await import(join(SKILL_DIR, key))).sort();
      assert.deepEqual([...ctx.result], real, key);
    }
  }
});

test('multi-declarator exports are bundled, including nested commas and comments', () => {
  const dir = mkdtempSync(join(tmpdir(), 'psim-bundle-'));
  writeFileSync(join(dir, 'm.mjs'), "export const a = 1, b = Math.max(2, 3), c = [4, 5] /* , d */;\nexport let e = 'x, y', f = { g: 1, h: 2 };\n");
  const code = bundleModules("import { a, b, c, e, f } from './m.mjs';\nresult = [a, b, c[1], e, f.h];", SKILL_DIR, dir);
  const ctx = { result: null };
  vm.runInNewContext(code, ctx);
  assert.deepEqual([...ctx.result], [1, 3, 5, 'x, y', 2]);
});

test('all new maths modules compute the same results in a self-contained bundle as in Node', async () => {
  const exercise = (L, F, P) => {
    const A = [[2, 1], [1, 3]], b = [1, 2], x = L.solveLinear(A, b);
    const eigen = L.eigenSymmetric(A), cg = L.conjugateGradient(A, b);
    const root = L.newtonSolve(v => [v[0] * v[0] - 2], [1]);
    const fit = L.leastSquares([[1, -1], [1, 0], [1, 1]], [1, 0, 3]);
    const re = new Float64Array([1, 0, -1, 0]), im = new Float64Array(4), shape = [4];
    F.fftND(re, im, shape); F.fftND(re, im, shape, true);
    const source = new Float64Array([2, 0, -2, 0]), grid = { shape, discretization: 'finite-difference' };
    const field = F.solvePoissonPeriodic(source, grid).u;
    const grad = F.gradientND(field, grid), lap = F.laplacianND(field, grid), div = F.divergenceND(grad, grid);
    const diffusion = F.diffusionPeriodic({ u0: field, ...grid, D: 0.1, dt: 0.1 }); diffusion.step();
    const quantum = F.schrodingerND({ re0: field, shape, dt: 0.01 }); quantum.step();
    const p = P.normalizeWeights([1, 2]), q = P.normalizeLogWeights([-3, -4]);
    const covariance = P.covarianceMatrix([[0, 1], [1, 3], [2, 5]]);
    const ensemble = P.canonicalEnsemble([0, 1, 2], { kT: 0.7, degeneracies: [1, 2, 1] });
    return JSON.stringify({ x, product: L.matrixVector(A, x), eigen, cg, root, fit, field, grad, lap, div,
      re, im, diffusion: diffusion.u, quantum: quantum.probability(), norm: quantum.norm(),
      p, q, H: P.shannonEntropy(p), KL: P.relativeEntropy(p, q), logZ: P.logSumExp([-3, -4]), covariance, ensemble });
  };
  const entry = `
import * as L from 'psim:core/linalg.mjs';
import * as F from 'psim:core/fields.mjs';
import * as P from 'psim:core/probability.mjs';
result = (${exercise.toString()})(L, F, P);
`;
  const source = `<script type="application/json" id="psim-verification">null</script><script type="module">${entry}</script>`;
  const html = bundleHtml(source, { skillDir: SKILL_DIR, verification: { checks: [] } });
  assert.deepEqual(auditHtml(html), []);
  const ctx = { result: null };
  vm.runInNewContext(html.match(/<script type="module">([\s\S]*?)<\/script>/)[1], ctx);
  const [L, F, P] = await Promise.all([import('../core/linalg.mjs'), import('../core/fields.mjs'), import('../core/probability.mjs')]);
  assert.equal(ctx.result, exercise(L, F, P));
});
