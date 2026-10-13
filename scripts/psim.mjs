#!/usr/bin/env node
// p-simtor command-line tool. Node >= 20, built-in modules only.
//
//   node scripts/psim.mjs gravity  [--preset NAME | --scenario FILE] [--out FILE] [--title T] [--subtitle S]
//   node scripts/psim.mjs verify   --scenario FILE [--out report.json]
//   node scripts/psim.mjs build    SOURCE.html --out FILE [--verify-module physics.mjs] [--config FILE]
//   node scripts/psim.mjs verify-module physics.mjs [--out report.json]
//   node scripts/psim.mjs check    PAGE.html [--shots DIR] [--wait MS]
//   node scripts/psim.mjs presets
//
// Run with no arguments for help.

import { readFileSync, writeFileSync, mkdirSync, existsSync } from 'node:fs';
import { dirname, resolve, join, basename } from 'node:path';
import { tmpdir } from 'node:os';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { bundleHtml, auditHtml } from './lib/bundle.mjs';
import './register.mjs';
import { checkPage } from './lib/chrome.mjs';
import { THREE_BODY_PRESETS } from '../core/presets.mjs';
import { verifyScenario } from '../core/nbody-verify.mjs';
import { INTEGRATORS } from '../core/nbody.mjs';

const SKILL_DIR = resolve(dirname(fileURLToPath(import.meta.url)), '..');

function parseArgs(argv) {
  const out = { _: [] };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (a.startsWith('--')) {
      const key = a.slice(2);
      const next = argv[i + 1];
      if (next == null || next.startsWith('--')) out[key] = true;
      else { out[key] = next; i++; }
    } else out._.push(a);
  }
  return out;
}

function fail(msg) {
  console.error(`psim: ${msg}`);
  process.exit(1);
}

function writeOut(file, text) {
  mkdirSync(dirname(resolve(file)), { recursive: true });
  writeFileSync(file, text);
}

const pct = x => x.toExponential(2);

// ------------------------------------------------------------------ scenario validation

function validateGravityScenario(s, where) {
  const errs = [];
  if (!Array.isArray(s.bodies) || s.bodies.length < 2) errs.push('bodies: need at least 2 bodies');
  else s.bodies.forEach((b, i) => {
    if (!(b.m > 0)) errs.push(`bodies[${i}].m must be > 0`);
    for (const k of ['x', 'v']) {
      if (!Array.isArray(b[k]) || b[k].length < 2 || b[k].some(q => !Number.isFinite(q))) errs.push(`bodies[${i}].${k} must be [x, y] or [x, y, z] numbers`);
    }
  });
  if (s.bodies?.length > 64) errs.push('bodies: the gravity template supports at most 64 bodies (direct O(N²) summation)');
  if (!(s.duration > 0)) errs.push('duration must be > 0 (physical time to verify and to scale playback)');
  if (s.integrator && !INTEGRATORS[s.integrator]) errs.push(`integrator must be one of ${Object.keys(INTEGRATORS).join(', ')}`);
  if (s.softening != null && !(s.softening >= 0)) errs.push('softening must be >= 0');
  if (s.units?.G != null && !(s.units.G > 0)) errs.push('units.G must be > 0');
  if (errs.length) fail(`invalid scenario${where ? ` in ${where}` : ''}:\n  - ${errs.join('\n  - ')}`);
  for (const b of s.bodies) { b.x = [...b.x, 0, 0].slice(0, 3); b.v = [...b.v, 0, 0].slice(0, 3); }
  return s;
}

function summarise(key, rep) {
  const c = rep.conservation, cv = rep.convergence;
  const lines = [
    `  ${key}: ${rep.integratorLabel}, t = 0 → ${rep.duration}, ${c.steps.toLocaleString()} steps`,
    `    |ΔE/E0| = ${pct(c.dE)}   ΔP = ${pct(c.dP)}   ΔL = ${pct(c.dL)}`,
    `    order ${cv.observedOrder == null ? 'unresolved' : cv.observedOrder.toFixed(2)} (expected ${cv.expectedOrder}),`
      + ` converged ${cv.converged ? 'over whole run' : `until t ≈ ${cv.horizon.toFixed(2)} (chaotic beyond)`}`,
    `    deterministic: ${rep.determinism.identical}   frame-rate independent: ${rep.frameIndependence.identical}`,
  ];
  const problems = [];
  if (c.dE > 1e-5) problems.push(`energy error ${pct(c.dE)} is large: increase stepsPerCrossing or use logh4`);
  if (cv.observedOrder != null && Math.abs(cv.observedOrder - cv.expectedOrder) > 0.6) problems.push('observed convergence order differs from theory');
  if (!rep.determinism.identical || !rep.frameIndependence.identical) problems.push('reproducibility check failed');
  for (const b of rep.benchmarks) if (b.matchesLiterature === false) problems.push(`benchmark failed: ${b.name}`);
  return { text: lines.join('\n'), problems };
}

// ------------------------------------------------------------------ commands

function cmdPresets() {
  console.log('Gravity presets (Hénon units, G = 1):');
  for (const [k, p] of Object.entries(THREE_BODY_PRESETS)) console.log(`  ${k.padEnd(14)} ${p.title}  —  ${p.source}`);
}

function loadScenario(args) {
  if (args.scenario) {
    const s = JSON.parse(readFileSync(args.scenario, 'utf8'));
    return validateGravityScenario(s, args.scenario);
  }
  return null;
}

function cmdGravity(args) {
  const custom = loadScenario(args);
  const presetKey = args.preset ?? (custom ? null : 'pythagorean');
  if (presetKey && !THREE_BODY_PRESETS[presetKey]) fail(`unknown preset "${presetKey}". Run: psim presets`);

  const scenarios = {};
  let initial;
  if (custom) {
    initial = custom.key ?? 'scenario';
    scenarios[initial] = {
      title: custom.title ?? 'Custom scenario', note: custom.note ?? '', source: custom.source ?? '',
      duration: custom.duration, speed: custom.speed, bodies: custom.bodies,
    };
  }
  const includePresets = custom ? args['with-presets'] === true || custom.includePresets === true : true;
  if (includePresets) {
    for (const [k, p] of Object.entries(THREE_BODY_PRESETS)) {
      scenarios[k] = { title: p.title, note: p.note, source: p.source, duration: p.duration, bodies: p.bodies };
    }
  }
  initial ??= presetKey;
  const base = custom ?? {};
  const integrator = base.integrator ?? 'logh4';
  const stepsPerCrossing = base.stepsPerCrossing ?? 2000;
  const softening = base.softening ?? 0;
  const units = base.units ?? { G: 1, system: 'Hénon N-body units (G = 1)', length: 'L', time: 'T', mass: 'M' };

  console.log('Verifying (conservation, convergence, determinism, frame independence, benchmarks)…');
  const verification = {};
  let problems = [];
  for (const [k, s] of Object.entries(scenarios)) {
    const rep = verifyScenario({ bodies: s.bodies, duration: s.duration, integrator, stepsPerCrossing, softening, G: units.G ?? 1 });
    verification[k] = rep;
    const sum = summarise(k, rep);
    console.log(sum.text);
    problems = problems.concat(sum.problems.map(p => `${k}: ${p}`));
  }

  const initialTitle = scenarios[initial].title;
  const config = {
    title: args.title ?? base.pageTitle ?? (custom ? initialTitle : 'Three-body problem'),
    subtitle: args.subtitle ?? base.subtitle ?? `${initialTitle} — Newtonian point masses, ${INTEGRATORS[integrator].label}`,
    tier: base.tier ?? (softening > 0 ? 'approximate' : 'validated'),
    tierReason: base.tierReason,
    units, integrator, stepsPerCrossing, softening,
    trailTime: base.trailTime ?? Math.max(1, scenarios[initial].duration / 4),
    notes: base.notes ?? [],
    initial, scenarios,
  };

  const out = args.out ?? 'simulations/three-body.html';
  const html = bundleHtml(readFileSync(join(SKILL_DIR, 'templates/gravity.html'), 'utf8'), { skillDir: SKILL_DIR, config, verification });
  const audit = auditHtml(html);
  writeOut(out, html);
  console.log(`\nWrote ${out} (${(html.length / 1024).toFixed(0)} KB, self-contained)`);
  if (audit.length) console.log(`Audit problems:\n  - ${audit.join('\n  - ')}`);
  if (problems.length) {
    console.log(`Verification warnings:\n  - ${problems.join('\n  - ')}`);
    process.exitCode = 2;
  }
}

function cmdVerify(args) {
  const s = loadScenario(args);
  if (!s) fail('verify needs --scenario FILE');
  const rep = verifyScenario({ ...s, G: s.units?.G ?? 1 });
  const sum = summarise(s.key ?? 'scenario', rep);
  console.log(sum.text);
  if (args.out) { writeOut(args.out, JSON.stringify(rep, null, 2)); console.log(`Wrote ${args.out}`); }
  if (sum.problems.length) { console.log(`Warnings:\n  - ${sum.problems.join('\n  - ')}`); process.exitCode = 2; }
}

/**
 * Run the verify() export of a physics module (the agent's own model) in Node.
 * Contract: verify() returns (or resolves to) { checks: [{ label, value, pass, state? }], notes?, summary? }.
 */
async function runModuleVerification(file) {
  const mod = await import(pathToFileURL(resolve(file)).href);
  if (typeof mod.verify !== 'function') fail(`${file} has no exported verify() function`);
  const t0 = Date.now();
  const rep = await mod.verify();
  if (!rep || !Array.isArray(rep.checks) || !rep.checks.length) fail(`${file}: verify() must return { checks: [...] } with at least one check`);
  rep.generated = new Date().toISOString();
  rep.runtimeMs = Date.now() - t0;
  const failed = rep.checks.filter(c => c.pass === false);
  for (const c of rep.checks) console.log(`  ${c.pass === false ? 'FAIL' : c.state === 'warn' ? 'warn' : ' ok '}  ${c.label}: ${c.value}`);
  if (failed.length) { console.log(`${failed.length} verification check(s) FAILED. Fix the model or method before delivering.`); process.exitCode = 2; }
  return rep;
}

async function cmdVerifyModule(args) {
  const file = args._[1];
  if (!file) fail('verify-module needs a physics module path');
  const rep = await runModuleVerification(file);
  if (args.out) { writeOut(args.out, JSON.stringify(rep, null, 2)); console.log(`Wrote ${args.out}`); }
}

async function cmdBuild(args) {
  const src = args._[1];
  if (!src) fail('build needs a SOURCE.html');
  if (!args.out) fail('build needs --out FILE');
  const config = args.config ? JSON.parse(readFileSync(args.config, 'utf8')) : undefined;
  let verification = args.verification ? JSON.parse(readFileSync(args.verification, 'utf8')) : undefined;
  if (args['verify-module']) {
    console.log(`Verifying ${args['verify-module']}…`);
    verification = await runModuleVerification(args['verify-module']);
  }
  const html = bundleHtml(readFileSync(src, 'utf8'), { skillDir: SKILL_DIR, srcDir: dirname(resolve(src)), config, verification });
  writeOut(args.out, html);
  const audit = auditHtml(html);
  console.log(`Wrote ${args.out} (${(html.length / 1024).toFixed(0)} KB)`);
  if (audit.length) { console.log(`Audit problems:\n  - ${audit.join('\n  - ')}`); process.exitCode = 2; }
}

async function cmdCheck(args) {
  const page = args._[1];
  if (!page || !existsSync(page)) fail('check needs an existing PAGE.html');
  const shots = args.shots ?? join(tmpdir(), 'psim-shots', basename(page, '.html'));
  const res = await checkPage(resolve(page), { shotsDir: shots, waitMs: +(args.wait ?? 3000) });
  for (const r of res.runs) {
    console.log(`  ${r.label.padEnd(18)} ${r.ok ? 'ok ' : 'FAIL'}  ${r.summary}`);
    if (r.screenshot) console.log(`  ${''.padEnd(18)} screenshot: ${r.screenshot}`);
  }
  if (res.problems.length) { console.log(`Problems:\n  - ${res.problems.join('\n  - ')}`); process.exitCode = 2; }
  else console.log('All runtime checks passed. Now look at every screenshot before delivering.');
}

const HELP = `p-simtor — scientifically verified interactive physics simulations

Commands:
  presets                                  list built-in gravity presets
  gravity [--preset NAME|--scenario FILE]  verify + build a gravity instrument page
          [--out FILE] [--title T] [--subtitle S] [--with-presets]
  verify  --scenario FILE [--out FILE]     run the gravity verification suite only
  build   SOURCE.html --out FILE           bundle a composed page (imports psim: and ./local modules)
          [--verify-module physics.mjs]    run the module's verify() and embed the results
          [--config FILE] [--verification FILE]
  verify-module physics.mjs [--out FILE]   run a physics module's own verification checks
  check   PAGE.html [--shots DIR]          headless Chrome: errors, liveness, overflow, UI→physics tests,
                                           screenshots (default: system temp dir)

Skill directory: ${SKILL_DIR}`;

const args = parseArgs(process.argv.slice(2));
const cmd = args._[0];
try {
  if (cmd === 'presets') cmdPresets();
  else if (cmd === 'gravity') cmdGravity(args);
  else if (cmd === 'verify') cmdVerify(args);
  else if (cmd === 'build') await cmdBuild(args);
  else if (cmd === 'verify-module') await cmdVerifyModule(args);
  else if (cmd === 'check') await cmdCheck(args);
  else console.log(HELP);
} catch (e) {
  fail(e.stack ?? e.message);
}
