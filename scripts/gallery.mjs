#!/usr/bin/env node
// The docs/ gallery (GitHub Pages) shows REAL demos: pages that the installed skill produced
// in a separate project from a real /p-simtor prompt. Nothing here builds simulations.
//
//   node scripts/gallery.mjs add PAGE.html --prompt "/p-simtor show why galaxies need dark matter" [--title T] [--id ID]
//   node scripts/gallery.mjs remove ID
//   node scripts/gallery.mjs list
//   node scripts/gallery.mjs order ID…      put these demos first, in this order
//   node scripts/gallery.mjs build           regenerate docs/index.html (add/remove do this too)
//
// `add` copies the page into docs/, takes a thumbnail (if Chrome is available) and records the
// prompt, the source path and the date in docs/demos.json.

import { copyFileSync, existsSync, mkdirSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { basename, dirname, join, resolve } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const DOCS = join(ROOT, 'docs');
const MANIFEST = join(DOCS, 'demos.json');
const { version } = JSON.parse(readFileSync(join(ROOT, 'package.json'), 'utf8'));

const args = process.argv.slice(2);
const opt = k => { const i = args.indexOf(`--${k}`); return i >= 0 ? args[i + 1] : undefined; };
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const fail = m => { console.error(`gallery: ${m}`); process.exit(1); };

mkdirSync(join(DOCS, 'thumbs'), { recursive: true });
const load = () => (existsSync(MANIFEST) ? JSON.parse(readFileSync(MANIFEST, 'utf8')) : { demos: [] });
const save = m => writeFileSync(MANIFEST, `${JSON.stringify(m, null, 2)}\n`);

/** Facts read from the page itself: title, verification checks, generator tag. */
function inspect(html) {
  const title = html.match(/<title>([^<]*)<\/title>/)?.[1]?.trim();
  const rep = html.match(/<script type="application\/json" id="psim-verification">([\s\S]*?)<\/script>/)?.[1];
  let checks = null;
  try {
    const r = JSON.parse(rep);
    const list = Array.isArray(r?.checks) ? r.checks : r && typeof r === 'object' ? Object.values(r).flatMap(x => x?.benchmarks ?? []) : [];
    if (list.length) checks = { total: list.length, passed: list.filter(c => c.pass !== false && c.matchesLiterature !== false).length };
  } catch { /* no report */ }
  return { title, checks, generator: /<meta name="generator" content="p-simtor"/.test(html) };
}

async function thumbnail(id) {
  try {
    const { launch } = await import('./lib/chrome.mjs');
    const b = await launch();
    await b.goto(`${pathToFileURL(join(DOCS, `${id}.html`)).href}?theme=dark&autoplay=1`, { width: 1280, height: 800, scale: 1 });
    await new Promise(r => setTimeout(r, 3500));
    const shot = await b.send('Page.captureScreenshot', { format: 'jpeg', quality: 82, clip: { x: 0, y: 64, width: 1280 - 372, height: 736, scale: 0.6 } });
    writeFileSync(join(DOCS, 'thumbs', `${id}.jpg`), Buffer.from(shot.data, 'base64'));
    await b.close();
    return true;
  } catch (e) {
    console.log(`  no thumbnail (${e.message.split('\n')[0]})`);
    return false;
  }
}

function build() {
  const { demos } = load();
  const cards = demos.map(d => `
    <a class="card" href="${esc(d.id)}.html">
      ${d.thumb ? `<img src="thumbs/${esc(d.id)}.jpg" alt="${esc(d.title)}" loading="lazy" width="545" height="442">` : '<div class="ph"></div>'}
      <div class="txt">
        <h2>${esc(d.title)}</h2>
        <div class="cmd"><span>Prompt</span><code>${esc(d.prompt)}</code></div>
        ${d.checks ? `<p class="chk"><span>Verification</span>${d.checks.passed} of ${d.checks.total} checks passed (details on the page)</p>` : ''}
        <p class="meta">Generated ${esc(d.date)} by an agent using the installed skill${d.skillVersion ? ` (p-simtor ${esc(d.skillVersion)})` : ''}</p>
      </div>
    </a>`).join('');
  const empty = '<p class="empty">No demos yet. Generate one in a test project with <code>/p-simtor …</code>, then add it with <code>node scripts/gallery.mjs add &lt;page.html&gt; --prompt "…"</code>.</p>';
  const html = `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>p-simtor demos</title>
<meta name="description" content="Interactive physics simulations created by the p-simtor agent skill from real prompts.">
<style>
:root { --bg: #07090d; --panel: #0e1219; --line: rgba(196,214,255,.10); --ink: #e6e9ef; --ink2: #b6bdcb; --muted: #7b8396; --accent: #8cc8ff; color-scheme: dark; }
@media (prefers-color-scheme: light) { :root { --bg: #efede7; --panel: #fbfaf7; --line: rgba(30,34,46,.10); --ink: #171a21; --ink2: #3a3f4b; --muted: #6b7180; --accent: #1c6fd1; color-scheme: light; } }
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink); font: 15px/1.55 ui-sans-serif, -apple-system, "Segoe UI", system-ui, sans-serif; }
main { max-width: 1180px; margin: 0 auto; padding: 56px 20px 72px; }
header h1 { margin: 0 0 10px; font-size: clamp(28px, 4vw, 40px); letter-spacing: -.02em; }
header p { margin: 0; max-width: 760px; color: var(--ink2); }
.tag { display: inline-block; margin-bottom: 14px; font: 600 11px/1 ui-monospace, Menlo, monospace; letter-spacing: .1em; text-transform: uppercase; color: var(--accent); }
.grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 20px; margin-top: 40px; }
.card { display: flex; flex-direction: column; background: var(--panel); border: 1px solid var(--line); border-radius: 14px; overflow: hidden; color: inherit; text-decoration: none; transition: transform .15s, border-color .15s; }
.card:hover { transform: translateY(-2px); border-color: var(--accent); }
.card img, .card .ph { width: 100%; height: auto; aspect-ratio: 545 / 442; object-fit: cover; display: block; background: #0a0d13; border-bottom: 1px solid var(--line); }
.txt { padding: 16px 18px 18px; display: grid; gap: 10px; }
.txt h2 { margin: 0; font-size: 17px; }
.cmd, .chk { display: grid; gap: 2px; margin: 0; font-size: 13px; color: var(--ink2); }
.cmd span, .chk span { font: 600 10.5px/1.6 ui-monospace, Menlo, monospace; letter-spacing: .08em; text-transform: uppercase; color: var(--accent); }
code { font: 12.5px ui-monospace, Menlo, monospace; color: var(--ink); overflow-wrap: anywhere; }
.meta { margin: 0; font-size: 12px; color: var(--muted); }
.empty { margin-top: 40px; color: var(--muted); }
</style>
</head>
<body>
<main>
  <header>
    <span class="tag">p-simtor ${esc(version)} · demos</span>
    <h1>Physics simulations created by p-simtor</h1>
    <p>Each demo was generated by a coding agent using the installed p-simtor skill, from the prompt shown on its card, in a separate project. Every page computes its motion from the real equations, shows its own verification, and runs entirely in your browser.</p>
  </header>
  ${demos.length ? `<section class="grid">${cards}\n  </section>` : empty}
</main>
</body>
</html>
`;
  writeFileSync(join(DOCS, 'index.html'), html);
  writeFileSync(join(DOCS, '.nojekyll'), '');
  console.log(`docs/index.html: ${demos.length} demo(s).  Open: ${pathToFileURL(join(DOCS, 'index.html')).href}`);
}

const cmd = args[0];
if (cmd === 'add') {
  const src = args[1];
  if (!src || !existsSync(src)) fail('add needs an existing PAGE.html (the page the skill built in your test project)');
  const prompt = opt('prompt');
  if (!prompt) fail('add needs --prompt "the exact /p-simtor prompt you typed"');
  const html = readFileSync(src, 'utf8');
  const info = inspect(html);
  if (!info.generator) fail(`${src} was not produced by p-simtor (no generator meta tag)`);
  const id = opt('id') ?? basename(src, '.html').toLowerCase().replace(/[^a-z0-9-]+/g, '-');
  const m = load();
  m.demos = m.demos.filter(d => d.id !== id);
  copyFileSync(src, join(DOCS, `${id}.html`));
  const thumb = !args.includes('--no-thumb') && await thumbnail(id);
  // version of the skill that built it, if the page sits in a project with the skill installed
  let skillVersion;
  for (let dir = dirname(resolve(src)); dir !== dirname(dir); dir = dirname(dir)) {
    const pkg = join(dir, '.claude', 'skills', 'p-simtor', 'package.json');
    if (existsSync(pkg)) { skillVersion = JSON.parse(readFileSync(pkg, 'utf8')).version; break; }
  }
  m.demos.push({
    id, title: opt('title') ?? info.title ?? id, prompt: prompt.startsWith('/p-simtor') ? prompt : `/p-simtor ${prompt}`,
    date: new Date().toISOString().slice(0, 10), source: resolve(src), checks: info.checks, thumb, skillVersion,
  });
  save(m);
  console.log(`Added "${id}"${info.checks ? ` (${info.checks.passed}/${info.checks.total} checks passed)` : ''}.`);
  build();
} else if (cmd === 'remove') {
  const id = args[1];
  const m = load();
  if (!m.demos.some(d => d.id === id)) fail(`no demo "${id}"`);
  m.demos = m.demos.filter(d => d.id !== id);
  save(m);
  rmSync(join(DOCS, `${id}.html`), { force: true });
  rmSync(join(DOCS, 'thumbs', `${id}.jpg`), { force: true });
  console.log(`Removed "${id}".`);
  build();
} else if (cmd === 'list') {
  const { demos } = load();
  if (!demos.length) console.log('No demos yet.');
  for (const d of demos) console.log(`  ${d.id.padEnd(28)} ${d.date}  ${d.prompt}`);
} else if (cmd === 'order') {
  // move the named demos to the front, in the given order; the rest keep their order after them
  const ids = args.slice(1), m = load();
  const unknown = ids.filter(id => !m.demos.some(d => d.id === id));
  if (!ids.length || unknown.length) fail(unknown.length ? `no demo "${unknown.join('", "')}"` : 'order needs demo ids (see list)');
  m.demos = [...ids.map(id => m.demos.find(d => d.id === id)), ...m.demos.filter(d => !ids.includes(d.id))];
  save(m);
  build();
} else if (cmd === 'build') {
  build();
} else {
  console.log(readFileSync(fileURLToPath(import.meta.url), 'utf8').split('\n').slice(1, 12).map(l => l.replace(/^\/\/ ?/, '')).join('\n'));
}
