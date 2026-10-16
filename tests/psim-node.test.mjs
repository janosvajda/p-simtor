// `psim node SCRIPT.mjs [args]` forwards every argument to the script as process.argv.slice(2).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { tmpdir } from 'node:os';
import { fileURLToPath } from 'node:url';

const psim = fileURLToPath(new URL('../scripts/psim.mjs', import.meta.url));

test('psim node forwards arguments, including ones that look like psim words', () => {
  const dir = mkdtempSync(join(tmpdir(), 'node-dir-'));
  const script = join(dir, 'node.mjs');
  writeFileSync(script, 'console.log(JSON.stringify(process.argv.slice(2)));');
  const out = execFileSync(process.execPath, [psim, 'node', script, '5', '--n', '3', 'node', 'check'], { encoding: 'utf8' });
  assert.deepEqual(JSON.parse(out.trim().split('\n').pop()), ['5', '--n', '3', 'node', 'check']);
});
