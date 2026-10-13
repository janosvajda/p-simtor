// Node module resolution hook: maps "psim:core/x.mjs" to the skill's own files, so a
// physics module written in the user's project imports p-simtor code the same way in
// Node (tests, verify-module) as in the bundled browser page.
import { fileURLToPath, pathToFileURL } from 'node:url';
import { dirname, join } from 'node:path';

const SKILL_DIR = join(dirname(fileURLToPath(import.meta.url)), '..', '..');

export async function resolve(specifier, context, next) {
  if (specifier.startsWith('psim:')) {
    return { url: pathToFileURL(join(SKILL_DIR, specifier.slice(5))).href, shortCircuit: true };
  }
  return next(specifier, context);
}
