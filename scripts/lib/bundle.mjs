// Turns an HTML source that references p-simtor modules into one self-contained file.
//
//   <link rel="stylesheet" href="psim:ui/instrument.css">        -> inlined <style>
//   <script type="module"> import { x } from 'psim:core/nbody.mjs' -> modules inlined,
//                                                                    each in its own scope
//   <script type="application/json" id="psim-config">               -> filled from --config
//   <script type="application/json" id="psim-verification">         -> filled from --verification
//
// Modules are wrapped individually (no shared top-level scope), so names never collide.
// Only the module syntax p-simtor itself uses is supported: named/namespace imports and
// `export function|const|let|class` declarations.

import { readFileSync } from 'node:fs';
import { join, dirname, normalize, resolve, relative } from 'node:path';

const IMPORT_RE = /^\s*import\s+(?:\*\s+as\s+(\w+)|\{([^}]*)\})\s+from\s+['"]([^'"]+)['"];?\s*$/gm;
const EXPORT_RE = /^export\s+(?:async\s+)?(function\*?|const|let|class)\s+([A-Za-z_$][\w$]*)/gm;

function parseImports(code) {
  const out = [];
  for (const m of code.matchAll(IMPORT_RE)) {
    out.push({ text: m[0], ns: m[1], names: m[2], spec: m[3] });
  }
  return out;
}

function bindingFor(imp, key) {
  const ref = `__psim[${JSON.stringify(key)}]`;
  if (imp.ns) return `const ${imp.ns} = ${ref};`;
  const parts = imp.names.split(',').map(s => s.trim()).filter(Boolean)
    .map(s => s.replace(/^(\w+)\s+as\s+(\w+)$/, '$1: $2'));
  return `const { ${parts.join(', ')} } = ${ref};`;
}

/**
 * Resolve a specifier to { key, file }. "psim:core/x.mjs" -> skill file; "./x.mjs" from the
 * page -> file next to the page (the agent's own physics module); relative imports inside
 * a module resolve against that module.
 */
function resolveSpec(spec, from, { skillDir, srcDir }) {
  if (spec.startsWith('psim:')) {
    const rel = normalize(spec.slice(5));
    return { key: rel, file: join(skillDir, rel) };
  }
  if (spec.startsWith('./') || spec.startsWith('../')) {
    if (!from) {
      if (!srcDir) throw new Error(`Relative import "${spec}" needs the page's directory`);
      const file = resolve(srcDir, spec);
      return { key: `local:${relative(srcDir, file)}`, file };
    }
    const file = resolve(dirname(from.file), spec);
    const key = from.key.startsWith("local:") ? `local:${relative(srcDir, file)}` : normalize(relative(skillDir, file));
    return { key, file };
  }
  throw new Error(`Unsupported import "${spec}" (only psim: and relative imports are bundled)`);
}

export function bundleModules(entryCode, skillDir, srcDir) {
  const order = [];
  const seen = new Set();
  const ctx = { skillDir, srcDir };
  const visit = (mod, stack) => {
    if (seen.has(mod.file)) return;
    if (stack.includes(mod.file)) throw new Error(`Import cycle: ${[...stack, mod.file].join(' -> ')}`);
    const code = readFileSync(mod.file, 'utf8');
    const imports = parseImports(code);
    let body = code;
    for (const imp of imports) {
      const dep = resolveSpec(imp.spec, mod, ctx);
      visit(dep, [...stack, mod.file]);
      body = body.replace(imp.text, bindingFor(imp, dep.key));
    }
    const names = [...body.matchAll(EXPORT_RE)].map(m => m[2]);
    body = body.replace(/^export\s+(?=(?:async\s+)?(?:function|const|let|class)\b)/gm, '');
    if (/^export\b/m.test(body)) throw new Error(`${mod.key}: unsupported export form (use export function/const/let/class)`);
    seen.add(mod.file);
    order.push(`// ---- ${mod.key.startsWith('local:') ? mod.key : `psim:${mod.key}`}\n__psim[${JSON.stringify(mod.key)}] = (() => {\n${body}\nreturn { ${names.join(', ')} };\n})();`);
  };
  let entry = entryCode;
  for (const imp of parseImports(entryCode)) {
    const dep = resolveSpec(imp.spec, null, ctx);
    visit(dep, []);
    entry = entry.replace(imp.text, bindingFor(imp, dep.key));
  }
  return `const __psim = {};\n${order.join('\n\n')}\n\n// ---- page\n${entry}`;
}

function jsonForScript(value) {
  return JSON.stringify(value, null, 1).replace(/<\//g, '<\\/');
}

export function bundleHtml(html, { skillDir, srcDir, config, verification } = {}) {
  html = html.replace(/<link\s+rel="stylesheet"\s+href="psim:([^"]+)"\s*\/?>/g, (_, p) =>
    `<style>/* psim:${p} */\n${readFileSync(join(skillDir, p), 'utf8')}</style>`);
  html = html.replace(/<script type="module">([\s\S]*?)<\/script>/g, (whole, code) =>
    /from\s+['"](?:psim:|\.\.?\/)/.test(code) ? `<script type="module">\n${bundleModules(code, skillDir, srcDir)}\n</script>` : whole);
  if (config !== undefined) {
    html = html.replace(/(<script type="application\/json" id="psim-config">)[\s\S]*?(<\/script>)/,
      (_, a, b) => `${a}${jsonForScript(config)}${b}`);
  }
  if (verification !== undefined) {
    html = html.replace(/(<script type="application\/json" id="psim-verification">)[\s\S]*?(<\/script>)/,
      (_, a, b) => `${a}${jsonForScript(verification)}${b}`);
  }
  return html;
}

/** Static checks on a finished page: self-contained, no unresolved references. */
export function auditHtml(html) {
  const problems = [];
  const unresolved = [...html.matchAll(/(?:href|src)="psim:[^"]+"|^\s*import\s.*from\s+['"][^'"]+['"]/gm)].map(m => m[0].trim());
  if (unresolved.length) problems.push(`Unresolved p-simtor references: ${unresolved.join(', ')}`);
  const external = [...html.matchAll(/<(?:script|link)[^>]+(?:src|href)="(https?:[^"]+)"/g)].map(m => m[1]);
  if (external.length) problems.push(`External resources (page will not work offline): ${external.join(', ')}`);
  if (/<script type="application\/json" id="psim-verification">\s*(null)?\s*<\/script>/.test(html)) {
    problems.push('No verification report embedded (psim-verification is empty)');
  }
  return problems;
}
