#!/usr/bin/env node
// Installs this skill for Claude Code, Codex and/or OpenCode. Node built-ins only.
//
//   node scripts/install.mjs                      all agents, for the current user (global)
//   node scripts/install.mjs --agent claude-code  one agent (claude-code | codex | opencode)
//   node scripts/install.mjs --project DIR        into one project instead of globally
//   node scripts/install.mjs --uninstall          remove what this script installed
//
// A skill is a folder containing SKILL.md; installing copies it to where each agent looks.

import { cpSync, existsSync, mkdirSync, readFileSync, rmSync, statSync } from 'node:fs';
import { homedir } from 'node:os';
import { basename, dirname, join, relative, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const SKILL_DIR = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const NAME = readFileSync(join(SKILL_DIR, 'SKILL.md'), 'utf8').match(/^name:\s*(\S+)/m)?.[1];
if (!NAME) { console.error('install: SKILL.md has no name'); process.exit(1); }

// Where each agent discovers skills.
const TARGETS = {
  'claude-code': { global: join(homedir(), '.claude', 'skills'), project: join('.claude', 'skills'), invoke: `/${NAME} <what to simulate>` },
  codex: { global: join(homedir(), '.agents', 'skills'), project: join('.agents', 'skills'), invoke: `$${NAME} <what to simulate>` },
  opencode: { global: join(homedir(), '.config', 'opencode', 'skills'), project: join('.agents', 'skills'), invoke: `/${NAME} <what to simulate>` },
};

// Repository files that are not part of the installed skill.
const EXCLUDE = new Set(['.git', '.github', 'node_modules', '.DS_Store', 'specification.md']);

const args = process.argv.slice(2);
const opt = k => { const i = args.indexOf(k); return i < 0 ? null : (args[i + 1] && !args[i + 1].startsWith('--') ? args[i + 1] : true); };
const agentArg = opt('--agent');
const projectArg = opt('--project');
const uninstall = args.includes('--uninstall');

if (args.includes('--help') || args.includes('-h')) {
  console.log(readFileSync(fileURLToPath(import.meta.url), 'utf8').split('\n').slice(1, 8).map(l => l.replace(/^\/\/ ?/, '')).join('\n'));
  process.exit(0);
}

const agents = agentArg && agentArg !== true ? agentArg.split(',') : Object.keys(TARGETS);
for (const a of agents) if (!TARGETS[a]) { console.error(`install: unknown agent "${a}" (use ${Object.keys(TARGETS).join(', ')})`); process.exit(1); }

const [major] = process.versions.node.split('.').map(Number);
if (major < 22) console.warn(`warning: Node ${process.versions.node} found; p-simtor needs Node 22+ for its browser checks.`);

const projectRoot = projectArg && projectArg !== true ? resolve(projectArg) : null;
const done = new Set();
for (const a of agents) {
  const base = projectRoot ? join(projectRoot, TARGETS[a].project) : TARGETS[a].global;
  const dest = join(base, NAME);
  if (done.has(dest)) { console.log(`  ${a.padEnd(12)} ${dest} (shared)`); continue; }
  done.add(dest);
  if (resolve(dest) === SKILL_DIR) { console.log(`  ${a.padEnd(12)} ${dest} (already the source)`); continue; }
  if (uninstall) {
    if (existsSync(join(dest, 'SKILL.md'))) { rmSync(dest, { recursive: true, force: true }); console.log(`  ${a.padEnd(12)} removed ${dest}`); }
    else console.log(`  ${a.padEnd(12)} not installed at ${dest}`);
    continue;
  }
  if (existsSync(dest) && !existsSync(join(dest, 'SKILL.md'))) {
    console.error(`install: ${dest} exists and is not a skill; leaving it alone`);
    process.exitCode = 1;
    continue;
  }
  rmSync(dest, { recursive: true, force: true });   // replace an older version cleanly
  mkdirSync(base, { recursive: true });
  cpSync(SKILL_DIR, dest, {
    recursive: true,
    filter: src => {
      const rel = relative(SKILL_DIR, src);
      return !rel.split(/[\\/]/).some(part => EXCLUDE.has(part)) && !/-shots$/.test(basename(src));
    },
  });
  const ok = existsSync(join(dest, 'SKILL.md')) && statSync(join(dest, 'scripts', 'psim.mjs')).isFile();
  console.log(`  ${a.padEnd(12)} ${ok ? 'installed' : 'FAILED'} → ${dest}`);
  if (!ok) process.exitCode = 1;
}

if (!uninstall && !process.exitCode) {
  console.log(`\nStart a new agent session${projectRoot ? ` in ${projectRoot}` : ''}, then:`);
  for (const a of agents) console.log(`  ${a.padEnd(12)} ${TARGETS[a].invoke}`);
}
