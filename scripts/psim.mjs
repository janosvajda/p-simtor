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

import { readFileSync, writeFileSync, mkdirSync, existsSync, copyFileSync } from 'node:fs';
import { spawn } from 'node:child_process';
import { dirname, resolve, join, basename } from 'node:path';
import { tmpdir } from 'node:os';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { bundleHtml, auditHtml } from './lib/bundle.mjs';
import './register.mjs';
import { checkPage } from './lib/chrome.mjs';
import { THREE_BODY_PRESETS } from '../core/presets.mjs';
import { solarSystemScenario, solarSystemScenarios, PLANETS } from '../core/solar-system.mjs';
import { verifyScenario } from '../core/nbody-verify.mjs';
import { INTEGRATORS } from '../core/nbody.mjs';
import { validateModuleReport, gravityReportProblems } from '../core/verification.mjs';

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

/** The platform's command for opening a file in the default browser. */
function openCommand(file) {
  if (process.platform === 'darwin') return ['open', [file]];
  if (process.platform === 'win32') return ['cmd', ['/c', 'start', '', file]];
  return ['xdg-open', [file]];
}

/** How to open a built page: clickable link plus a terminal command. Printed after every build. */
function printOpenHint(file) {
  const abs = resolve(file);
  const [cmd, cmdArgs] = openCommand(abs);
  const shown = [cmd, ...cmdArgs.map(a => (a.includes(' ') ? `"${a}"` : a))].join(' ');
  console.log(`\n━━ Open it (end your reply with these lines) ━━\n  link:     ${pathToFileURL(abs).href}\n  terminal: ${shown}\n  or:       node "${join(SKILL_DIR, 'scripts', 'psim.mjs')}" open "${abs}"`);
}

// ------------------------------------------------------------------ recipes (start from a verified example)

function loadCatalog() {
  return JSON.parse(readFileSync(join(SKILL_DIR, 'examples', 'catalog.json'), 'utf8')).recipes;
}

/** Rank recipes against a free-text request by keyword hits. */
function matchRecipes(query) {
  const q = ` ${query.toLowerCase().replace(/[^a-z0-9]+/g, ' ')} `;
  return loadCatalog()
    .map(r => ({ r, score: r.keywords.reduce((s, k) => s + (q.includes(` ${k.replace(/[^a-z0-9]+/g, ' ')} `) ? k.split(' ').length : 0), 0) }))
    .sort((a, b) => b.score - a.score);
}

function cmdRecipes(args) {
  const query = args._.slice(1).join(' ');
  const list = query ? matchRecipes(query) : loadCatalog().map(r => ({ r, score: null }));
  console.log(query ? `Recipes ranked for "${query}":` : 'Verified recipes (start from the closest one):');
  for (const { r, score } of list) {
    if (query && score === 0) continue;
    console.log(`\n  ${r.id}${score ? `  (match ${score})` : ''}\n    ${r.title}`);
    console.log(`    start: ${r.kind === 'gravity' ? `node psim.mjs ${r.command}` : `node psim.mjs new ${r.id} simulations/${r.id}`}`);
    console.log(`    adapt: ${r.adapt}`);
  }
  if (query && list.every(x => x.score === 0)) console.log('  No close recipe: compose from the toolkit (references/compose.md), starting from the most similar layout.');
}

function cmdNew(args) {
  const id = args._[1];
  const r = loadCatalog().find(x => x.id === id);
  if (!r) fail(`unknown recipe "${id}". Run: psim recipes`);
  if (r.kind === 'gravity') { console.log(`This recipe needs no code. Run:\n  node "${join(SKILL_DIR, 'scripts', 'psim.mjs')}" ${r.command}`); return; }
  const dest = args._[2] ?? join('simulations', id);
  mkdirSync(dest, { recursive: true });
  for (const f of ['physics.mjs', 'page.html']) {
    const to = join(dest, f);
    if (existsSync(to) && !args.force) fail(`${to} exists (use --force to overwrite)`);
    copyFileSync(join(SKILL_DIR, r.dir, f), to);
  }
  const psim = `node "${join(SKILL_DIR, 'scripts', 'psim.mjs')}"`;
  console.log(`Copied the verified "${id}" recipe to ${dest}/ (physics.mjs, page.html).`);
  console.log(`It already builds and passes its checks. Change only what the request needs, then:`);
  console.log(`  ${psim} build ${join(dest, 'page.html')} --verify-module ${join(dest, 'physics.mjs')} --out ${dest}.html`);
  console.log(`  ${psim} check ${dest}.html`);
  console.log(`Adapt: ${r.adapt}`);
}

function cmdOpen(args) {
  const file = args._[1];
  if (!file || !existsSync(file)) fail('open needs an existing PAGE.html');
  const [cmd, cmdArgs] = openCommand(resolve(file));
  const child = spawn(cmd, cmdArgs, { stdio: 'ignore', detached: true });
  child.on('error', () => fail(`could not launch a browser; open ${pathToFileURL(resolve(file)).href} manually`));
  child.unref();
  console.log(`Opened ${pathToFileURL(resolve(file)).href}`);
}

// ------------------------------------------------------------------ scenario validation

function validateGravityScenario(s, where) {
  const errs = [];
  if (!Array.isArray(s.bodies) || s.bodies.length < 2) errs.push('bodies: need at least 2 bodies');
  else {s.bodies.forEach((b, i) => {
    if (!(Number.isFinite(b.m) && b.m > 0)) errs.push(`bodies[${i}].m must be finite and > 0`);
    for (const k of ['x', 'v']) {
      if (!Array.isArray(b[k]) || b[k].length < 2 || b[k].some(q => !Number.isFinite(q))) errs.push(`bodies[${i}].${k} must be [x, y] or [x, y, z] numbers`);
    }
  });}
  if (s.bodies?.length > 64) errs.push('bodies: the gravity template supports at most 64 bodies (direct O(N²) summation)');
  if (!(Number.isFinite(s.duration) && s.duration > 0)) errs.push('duration must be finite and > 0 (physical time to verify and to scale playback)');
  if (s.integrator && !INTEGRATORS[s.integrator]) errs.push(`integrator must be one of ${Object.keys(INTEGRATORS).join(', ')}`);
  if (s.softening != null && !(Number.isFinite(s.softening) && s.softening >= 0)) errs.push('softening must be finite and >= 0');
  if (s.units?.G != null && !(Number.isFinite(s.units.G) && s.units.G > 0)) errs.push('units.G must be finite and > 0');
  if (s.stepsPerCrossing != null && !(Number.isFinite(s.stepsPerCrossing) && s.stepsPerCrossing > 0)) errs.push('stepsPerCrossing must be finite and > 0');
  if (errs.length) fail(`invalid scenario${where ? ` in ${where}` : ''}:\n  - ${errs.join('\n  - ')}`);
  for (const b of s.bodies) { b.x = [...b.x, 0, 0].slice(0, 3); b.v = [...b.v, 0, 0].slice(0, 3); }
  return s;
}

function summarise(key, rep) {
  const c = rep.conservation, cv = rep.convergence;
  const orderNote = cv.orderH > c.h * (1 + 1e-12) ? ` measured at coarser h = ${pct(cv.orderH)},` : '';
  const lines = [
    `  ${key}: ${rep.integratorLabel}, t = 0 → ${rep.duration}, ${c.steps.toLocaleString()} steps`,
    `    max |ΔE/E0| = ${pct(c.dE)}   max ΔP = ${pct(c.dP)}   max ΔL = ${pct(c.dL)}`,
    `    order ${cv.observedOrder == null ? 'unresolved' : cv.observedOrder.toFixed(2)} (expected ${cv.expectedOrder}),${orderNote} converged ${cv.converged ? 'over whole run' : `until t ≈ ${cv.horizon.toFixed(2)} (chaotic beyond)`}`,
    `    deterministic: ${rep.determinism.identical}   frame-rate independent: ${rep.frameIndependence.identical}`,
  ];
  const problems = gravityReportProblems(rep);
  return { text: lines.join('\n'), problems };
}

// ------------------------------------------------------------------ commands

function cmdPresets() {
  console.log('Gravity presets (Hénon units, G = 1):');
  for (const [k, p] of Object.entries(THREE_BODY_PRESETS)) console.log(`  ${k.padEnd(14)} ${p.title}  —  ${p.source}`);
  console.log(`  ${'solar-system'.padEnd(14)} Sun + planets at J2000 from JPL elements (AU, yr)  [--planets ${PLANETS.map(p => p.name).join(',')}] [--years N]`);
}

function loadScenario(args) {
  if (args.scenario) {
    const s = JSON.parse(readFileSync(args.scenario, 'utf8'));
    return validateGravityScenario(s, args.scenario);
  }
  return null;
}

function cmdGravity(args) {
  let custom = loadScenario(args);
  let extra = [];
  if (args.preset === 'solar-system') {
    if (typeof args.planets === 'string' || args.years) {
      const planets = typeof args.planets === 'string' ? args.planets.split(',') : undefined;
      custom = validateGravityScenario(solarSystemScenario({ planets, years: args.years ? +args.years : undefined }));
    } else {
      [custom, ...extra] = solarSystemScenarios().map(sc => validateGravityScenario(sc));
    }
    delete args.preset;
  }
  const presetKey = args.preset ?? (custom ? null : 'pythagorean');
  if (presetKey && !THREE_BODY_PRESETS[presetKey]) fail(`unknown preset "${presetKey}". Run: psim presets`);

  const scenarios = {};
  let initial;
  if (custom) {
    initial = custom.key ?? 'scenario';
    scenarios[initial] = {
      title: custom.title ?? 'Custom scenario', note: custom.note ?? '', source: custom.source ?? '',
      duration: custom.duration, speed: custom.speed, bodies: custom.bodies, trailTime: custom.trailTime, fitRadius: custom.fitRadius, stepsPerCrossing: custom.stepsPerCrossing,
    };
    for (const sc of extra) {
      scenarios[sc.key] = { title: sc.title, note: sc.note ?? '', source: sc.source ?? '', duration: sc.duration, speed: sc.speed, bodies: sc.bodies, trailTime: sc.trailTime, fitRadius: sc.fitRadius, stepsPerCrossing: sc.stepsPerCrossing };
    }
  }
  const includePresets = custom ? args['with-presets'] === true || custom.includePresets === true : true;
  if (includePresets) {
    for (const [k, p] of Object.entries(THREE_BODY_PRESETS)) {
      scenarios[k] = { title: p.title, note: p.note, source: p.source, duration: p.duration, bodies: p.bodies, stepsPerCrossing: p.stepsPerCrossing };
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
    const rep = verifyScenario({ bodies: s.bodies, duration: s.duration, integrator, stepsPerCrossing: s.stepsPerCrossing ?? stepsPerCrossing, softening, G: units.G ?? 1 });
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
  printOpenHint(out);
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
  validateModuleReport(rep);
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
  if (verification != null) {
    validateModuleReport(verification);
    if (verification.checks.some(c => c.pass === false)) process.exitCode = 2;
  }
  const html = bundleHtml(readFileSync(src, 'utf8'), { skillDir: SKILL_DIR, srcDir: dirname(resolve(src)), config, verification });
  writeOut(args.out, html);
  const audit = auditHtml(html);
  console.log(`Wrote ${args.out} (${(html.length / 1024).toFixed(0)} KB)`);
  printOpenHint(args.out);
  if (audit.length) { console.log(`Audit problems:\n  - ${audit.join('\n  - ')}`); process.exitCode = 2; }
}

async function cmdCheck(args) {
  const page = args._[1];
  if (!page || !existsSync(page)) fail('check needs an existing PAGE.html');
  const shots = args.shots ?? join(tmpdir(), 'psim-shots', basename(page, '.html'));
  const res = await checkPage(resolve(page), { shotsDir: shots, waitMs: +(args.wait ?? 3000), advance: args.advance });
  for (const r of res.runs) {
    console.log(`  ${r.label.padEnd(18)} ${r.ok ? 'ok ' : 'FAIL'}  ${r.summary}`);
    if (r.screenshot) console.log(`  ${''.padEnd(18)} screenshot: ${r.screenshot}`);
    for (const f of r.more ?? []) console.log(`  ${''.padEnd(18)} scrolled:   ${f}`);
  }
  if (res.problems.length) { console.log(`Problems:\n  - ${res.problems.join('\n  - ')}`); process.exitCode = 2; }
  else console.log('All runtime checks passed. Now look at every screenshot, including the scrolled ones, before delivering.');
}

const HELP = `p-simtor — scientifically verified interactive physics simulations

Commands:
  recipes [request words]                  list verified starting points, ranked for a request
  new     RECIPE [DIR]                     copy a verified example into the project to adapt
  presets                                  list built-in gravity presets
  gravity [--preset NAME|--scenario FILE]  verify + build a gravity instrument page
          [--out FILE] [--title T] [--subtitle S] [--with-presets]
  verify  --scenario FILE [--out FILE]     run the gravity verification suite only
  build   SOURCE.html --out FILE           bundle a composed page (imports psim: and ./local modules)
          [--verify-module physics.mjs]    run the module's verify() and embed the results
          [--config FILE] [--verification FILE]
  verify-module physics.mjs [--out FILE]   run a physics module's own verification checks
  open    PAGE.html                        open a built page in the default browser
  check   PAGE.html [--shots DIR] [--advance T]   headless Chrome: errors, liveness, overflow, UI→physics tests,
                                           extreme inputs, screenshots (default: system temp dir);
                                           --advance T fast-forwards to time T first

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
  else if (cmd === 'open') cmdOpen(args);
  else if (cmd === 'recipes') cmdRecipes(args);
  else if (cmd === 'new') cmdNew(args);
  else console.log(HELP);
} catch (e) {
  fail(e.stack ?? e.message);
}
