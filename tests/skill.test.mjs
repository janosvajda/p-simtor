// Validates SKILL.md against the Agent Skills specification (agentskills.io/specification)
// and checks that every file the skill tells an agent to use actually exists.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, existsSync, readdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join, basename, resolve } from 'node:path';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const text = readFileSync(join(ROOT, 'SKILL.md'), 'utf8');
const m = text.match(/^---\n([\s\S]*?)\n---\n([\s\S]*)$/);

function frontmatter() {
  assert.ok(m, 'SKILL.md must start with YAML frontmatter');
  const fm = {};
  let parent = null;
  for (const line of m[1].split('\n')) {
    const top = line.match(/^([\w-]+):\s*(.*)$/);
    const nested = line.match(/^\s+([\w/-]+):\s*(.*)$/);
    if (top) { parent = top[1]; fm[parent] = top[2] === '' ? {} : top[2].replace(/^"(.*)"$/, '$1'); }
    else if (nested && parent) fm[parent][nested[1]] = nested[2].replace(/^"(.*)"$/, '$1');
  }
  return fm;
}

test('frontmatter follows the Agent Skills spec', () => {
  const fm = frontmatter();
  const allowed = ['name', 'description', 'license', 'compatibility', 'metadata', 'allowed-tools'];
  for (const k of Object.keys(fm)) assert.ok(allowed.includes(k), `unexpected frontmatter field "${k}"`);
  assert.match(fm.name, /^[a-z0-9]+(-[a-z0-9]+)*$/, 'name: lowercase letters, digits, single hyphens');
  assert.ok(fm.name.length <= 64);
  assert.ok(fm.description.length > 0 && fm.description.length <= 1024, `description length ${fm.description.length}`);
  assert.ok(!fm.compatibility || fm.compatibility.length <= 500);
  for (const v of Object.values(fm.metadata ?? {})) assert.equal(typeof v, 'string');
});

test('skill name matches the directory it is installed under', () => {
  // The repository root is the skill; installers name the folder after the skill.
  const fm = frontmatter();
  assert.equal(fm.name, 'p-simtor');
  const dir = basename(ROOT);
  if (dir !== 'p-simtor') console.log(`  note: checked out as "${dir}"; installers create "p-simtor"`);
});

test('SKILL.md body stays within the recommended size', () => {
  assert.ok(text.split('\n').length < 500, 'keep SKILL.md under 500 lines');
  assert.ok(m[2].split(/\s+/).length < 3500, 'keep instructions compact (< ~5000 tokens)');
});

test('every file referenced from SKILL.md and references/ exists', () => {
  const docs = [join(ROOT, 'SKILL.md'), ...readdirSync(join(ROOT, 'references')).map(f => join(ROOT, 'references', f))];
  for (const doc of docs) {
    const body = readFileSync(doc, 'utf8');
    for (const [, link] of body.matchAll(/\]\(([^)\s]+)\)/g)) {
      if (link.startsWith('#') || /^[a-z][a-z\d+.-]*:/i.test(link)) continue;
      const path = decodeURIComponent(link.split(/[?#]/)[0]);
      assert.ok(existsSync(resolve(dirname(doc), path)), `${basename(doc)} links to missing ${link}`);
    }
    for (const [, path] of body.matchAll(/`((?:core|ui|scripts|templates|examples)\/[\w./-]+\.(?:mjs|html|css))`/g)) {
      assert.ok(existsSync(join(ROOT, path)), `${basename(doc)} mentions missing ${path}`);
    }
  }
});

test('evals file is well-formed', () => {
  const ev = JSON.parse(readFileSync(join(ROOT, 'evals/evals.json'), 'utf8'));
  assert.equal(ev.skill_name, 'p-simtor');
  assert.ok(ev.evals.length >= 5);
  for (const e of ev.evals) {
    assert.equal(typeof e.prompt, 'string');
    assert.equal(typeof e.expected_output, 'string');
    assert.ok(Array.isArray(e.files));
  }
});

test('fact-check step is part of the workflow and its reference exists', () => {
  assert.match(text, /Fact check:/);
  assert.match(text, /references\/facts\.md/);
  const facts = readFileSync(join(ROOT, 'references/facts.md'), 'utf8');
  assert.match(facts, /Mercury or Venus has moons/);
});
