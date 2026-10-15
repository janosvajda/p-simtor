# Composing a simulation for any domain

Use this when no built-in model fits the scientific question. First define that question, state and observables using [method selection](method-selection.md). The reference implementation is `examples/proton-magnetic-field/` (`physics.mjs`, `page.html` and the built page). Reuse its structure, replacing its physics and verification as needed.

## Layout in the user's project

```
simulations/
  hydrogen-2p/
    physics.mjs          # the model + verify(); pure JS, no DOM
    physics.test.mjs     # optional node:test unit tests
    page.html            # createInstrument({...}) + draw()
  hydrogen-2p.html       # built, self-contained deliverable
```

## physics.mjs: the simulation contract

```js
import { me, hbar, a0 } from 'psim:core/constants.mjs';   // resolves in Node and in the bundle

export const PARAMS = [   // generates the controls; any change rebuilds the sim via createSim
  { key: 'B', label: 'Magnetic field B', min: 0.01, max: 10, value: 1, log: true, unit: 'T', hint: '…' },
  { key: 'n', label: 'Particles', min: 100, max: 5000, value: 1000, log: true, integer: true },   // integer counts
  { key: 'state', label: 'Orbital', value: '2p_z', options: [{ value: '2p_z', label: '2p (m = 0)' }] },
];
// Typed values are clamped to [min, max] (or to `limits: [lo, hi]` when a wider range is
// physically meaningful). Set max to what the page can compute smoothly: `psim check`
// types 1e15, −1e15 and 0 into every input and fails if the page freezes or errors.

export function createSim(params) {
  return {
    t,                       // getter: physical time, or clearly labelled sampling/scan coordinate
    step(),                  // one dynamics or sampling step for evolving models
    diagnostics(),           // [{ key, label, value, kind: 'error'|'value', unit?, good?, warn?, text?, pending? }]
    state(),                 // JSON-serialisable state (export + automated checks)
    record?(),               // { header, rows } for CSV export
    // plus anything draw() needs (trail, grids, analytic reference values…)
  };
}

export function verify({ log } = {}) {   // run by psim verify-module / build --verify-module
  log?.('2D exact check');                // optional progress lines for slow (Monte Carlo) checks
  return { summary, checks: [{ label, value, pass }], notes };
}
```

Rules:
- Keep calculation, observation and representation distinct. The calculation updates/samples the model state; observation functions derive the requested quantities; `draw()` maps these results to the declared visual channels. Document units, conditions, provenance and the progress axis alongside the model. A new subject uses this same contract without needing a new recipe or domain guide.
- `createSim` must be deterministic for given params. Use a seeded RNG (e.g. mulberry32) for stochastic models and expose the seed as a param.
- Diagnostics of `kind: 'error'` are relative errors (plotted on a log axis; the first one also appears in the HUD and in `psim check`'s summary under its own label). Use `good`/`warn` thresholds from `verification.md`. A model with no conserved quantity to drift (Monte Carlo, sampling, data) needs no error diagnostic: show statistical results as `kind: 'value'` with `text: '0.314 ± 0.004'`, or a pull against an exact value as the error.
- A quantity that is not measured yet (too few samples) is `{ …, pending: true }`: the row shows "measuring…" (`texts.pending`) and the checks accept it. Never show "—", NaN or ∞ as a value: `psim check` fails it and names the row.
- `verify()` compares against **exact or published** results, not against the code itself. Use refinement for numerical dynamics and sampling/equilibration checks for statistical methods; see [verification](verification.md).
- Static (time-independent) visualisations: pass `timeDependent: false` to `createInstrument`. Then `step()` is not needed.
- Set the title/subtitle, `tierReason`, `model.assumptions`, `model.conventions` and references from the selected specification. Derive each legend/convention from the renderer's mapping of quantities to visual channels. Use `timeUnit` and `timeLabel` (what the clock measures, e.g. "time since strike"; default "simulated time") for the actual progress coordinate, so the readout never shows a bare number; a static parameter scan needs no artificial playback.

## The toolkit: don't rewrite it

All general mathematics (numerical primitives, ODEs and particles, matrices, grid fields and PDEs, probability, stochastic dynamics, data) is catalogued in [math-toolkit.md](math-toolkit.md), with APIs, conventions, conditions and benchmarks. The table below lists the remaining components.

| Need | Use |
|---|---|
| Plots: log/linear axes with 1-2-5 ticks, function and data series, error bars, shaded validity bands, reference lines, markers, draggable cursor, highlight linked to equations | `psim:ui/plot.mjs` → `createPlot` |
| Equations split into physical terms with live values and colour links; a "where" block; numbered derivations | `psim:ui/equation.mjs` → `createEquation`, `createDerivation`, `M` helpers |
| Black-body, colour | `psim:core/radiation.mjs` |
| Geodesics, orbits, redshift | `psim:core/schwarzschild.mjs` |
| Lattice gauge theory (SU(2) sampling, all Wilson loops, Polyakov loop, plaquette field, APE smearing) | `psim:core/lattice-gauge.mjs` |
| QCD running coupling α_s(Q) (MS-bar, 1–4 loops, flavour thresholds) | `psim:core/qcd.mjs` |
| Particle data with uncertainties and edition (masses, lifetimes, G_F, α_s, Michel parameters) | `PDG` in `psim:core/constants.mjs` |
| Relativistic kinematics: 4-vectors, boosts, two-body decay, RAMBO phase space, accept–reject | `psim:core/kinematics.mjs` |

## page.html

```html
<link rel="stylesheet" href="psim:ui/instrument.css">
<script type="application/json" id="psim-verification">null</script>
<script type="module">
import { createInstrument } from 'psim:ui/instrument-app.mjs';
import { PARAMS, createSim } from './physics.mjs';
createInstrument({
  title, subtitle, tier: 'validated' | 'approximate' | 'analytical' | 'conceptual', tierReason,
  params: PARAMS, create: createSim,
  draw(ctx, sim, view, t) { /* view.W, view.H (CSS px), view.color('--body-1'), view.color('--glow') */ },
  onPointer(canvas, view, getSim) { /* camera drag etc. */ },
  timeUnit: 'ns', timeLabel: 'time since injection', speed: 30, duration: 1000,   // HUD shows "time since injection  12.4 ns"
  hudLine: sim => '…',
  model: { legend: [], equations: [{ caption, mathml }], assumptions: [], conventions: [], references: [] },
});
</script>
```

- `draw` is called every frame with a cleared, DPR-scaled 2D context. It must only **read** `sim`.
- **Workbench layout** (`layout: 'workbench'`): the stage becomes a card at the top, followed by equations and figures. Pass `bench({ eqs, figs, createPlot, createEquation, createDerivation, M, setParam, getSim })`, which builds them and returns `{ update(sim, t, { rebuilt }) }`. It's called after every rebuild and about 10×/s. See `examples/blackbody/page.html`.
- **Equations**: give every term a `value: () => …` that reads live state, a `color` matching its plot series id, and `link: plot` so hovering highlights the curve.
- **Fast-forward support**: `psim check --advance T` calls `window.psim.advanceTo(T)` and then `window.psim.fit?.()`. If your view must re-frame after a long run, define `window.psim.fit`.
- **GPU or other self-tests**: define `window.psim.selfTest = () => ({ label, pass, …measurements })` (or an array of them) after `createInstrument`. `psim check` runs it on a fresh page loaded with `?autoplay=0`, so the sim is at t = 0: advance it yourself (`window.psim.advanceTo(t)`) before testing anything that needs samples. It may be `async` (return a Promise). A self-test that reads GPU pixels must draw first (call your renderers, as `drawGL` does) and read with `gl.readPixels` before the next frame clears the canvas. See `tests/fixtures/gl-check/page.html`.
- **Several calculations on one page** (e.g. an analytical curve and a Monte Carlo run): keep the expensive one in `PARAMS`/`createSim` (a change rebuilds it), and give the cheap one bench-local controls (`slider`/`select` from `ui/instrument.mjs` inside `bench`) that only redraw its figure. Label each figure's tier in its caption ("analytical", "Monte Carlo, SU(2) 6⁴"); the header tier describes the page as a whole (the least certain part).
- **Slow steps** (a Monte Carlo sweep of 10–30 ms): the shell checks its per-frame budget (`frameBudget`, default 9 ms) after every step, so at most one slow step runs past it, and the HUD says "compute-limited" when playback cannot keep up. Set `speed` to what one frame can compute and `speedRange: [lo, hi]` so the speed slider cannot ask for 100× that. Measure a step with `psim node timing.mjs` (below) before choosing.
- **Discrete progress axes** (sweeps, generations, events): `timeUnit: 'sweeps', timeDigits: 0` shows "122 sweeps", not "122.800".
- **GPU scenes** (`psim:ui/gl.mjs`): pass `gl: true` and `drawGL(gl, sim, view, t)` to `createInstrument` (works with both layouts, including `workbench`). The shell stacks a WebGL2 canvas under the 2D one; `draw` still runs on top for labels, axes and scale bars, and PNG export combines both. Renderers, by what is drawn:
  - `createPointRenderer(gl)`: `update(positions, { values | colors })` or `append(…)` for points that arrive during the run (detector hits, decays; `clear()` empties), `draw(camera, { size, sizeUnit: 'px' | 'world', colormap, range, log, blend: 'additive' | 'alpha', soft, count })`: 10⁵–10⁶ glowing bodies, stars, molecules or sampled positions. `opacity` scales every point. In the light theme use `blend: 'alpha', soft: true` so points keep the same Gaussian look and apparent size as the dark theme's glow. A growing cloud (accumulated events): allocate the full buffer once with `update()`, then `draw({ count: n })`; never re-upload every frame.
  - `createFieldRenderer(gl)`: `update(values, nx, ny, { range, log })`, `draw({ rect, colormap, smooth })`: a 2D field as a colour-mapped image (NaN cells transparent; `smooth: false` shows the grid honestly).
  - `createVolumeRenderer(gl)`: `update(values, nx, ny, nz)`, `draw(camera, { box, colormap, opacity, steps, mode: 'emission' | 'absorption', smooth })`: ray-marched 3D densities (orbital clouds, gas); `smooth: false` shows each cell as a block (a coarse lattice). Opacity follows the normalised value and values ≤ 0 are transparent; for a field of either sign (charge density, vorticity, Re ψ) use `update(…, { signed: true })` with `colormap: 'coolwarm'`: hue gives the sign, opacity |value|. In `'emission'` mode high `opacity` adds colours up to white; use `'absorption'` when hue must stay readable.
  - `createSphereRenderer(gl)`: `draw(camera, [{ center, radius, color, emissive }], { light, ambient })`: planets, moons, balls, atoms as exact ray-cast spheres; `light` is a point (each sphere lit from its own direction: correct day/night line) or `{ direction }`; `emissive` for a glowing Sun.
  - `createLineRenderer(gl)`: `update([[x,y,z,…], …])`, `draw(camera, { width, color, opacity, dash: [on, off], depthTest })`: orbits, trajectories, field lines, axes; constant pixel width; dashed for reference curves. A line may be `{ points, color, opacity }`, so differently coloured or age-faded tracks need one renderer. `polylines.circle(center, normal, r)`, `.boxEdges([min, max])`, `.grid(dims, { spacing, origin })` (lattice lines) and `.arrow(from, to)` build common shapes.
  - `drawColorBar(ctx, { x, y, width, height, colormap, range, log, label })` in `draw()`: the colour-bar legend for any colour-mapped GPU quantity (ticks at nice values or decades). Every colour mapping needs one, with the quantity and unit as its label.
  - Large coordinates (astronomical distances): spheres and lines are drawn relative to `camera.target` in float64, so put the target near what you look at; `camera.matrices(aspect, { relative: true })` gives the matching matrices for custom shaders.
  - `createMeshRenderer(gl)` with `geometry.box/cylinder/cone/sphere/disc/plane/arrow(...)`, `transform({ translate, rotate, scale })` (`rotate` is Euler angles or a 3×3 matrix such as `alignZ(direction)`, which points an arrow, cone or cylinder along a vector) and `compose(a, b)`: solid lit objects (apparatus, playground, buildings, planets as solids). `draw(camera, [{ mesh, matrix, color, opacity, noShadow }], { light: { direction }, ambient, shadows: { z, opacity } })`; shadows are exact projections onto the ground plane for parallel light (e.g. `sunDirection()` from `core/sun.mjs`), overlapping shadows do not darken twice. Draw 2D labels on top with `camera.project`.
  - `orbitCamera({ target, distance, yaw, pitch })` (z up): `attach(canvas)` in `onPointer` for drag and zoom (on a scrolling page the plain wheel scrolls; Ctrl/⌘ + wheel or pinch zooms); `fit(box, view.W / view.H, margin = 1.1)` frames a bounding box on any screen shape (call it in `onFit` and when the aspect changes, so phones see the whole scene). It fits the box's circumscribed sphere: for a ball of radius r pass `[[-r, -r, -r], [r, r, r]].map(c => c.map(v => v / Math.sqrt(3)))`, or the box will look √3 too small; `minDistance`/`maxDistance`/`minPitch`/`maxPitch` limit the camera; `project([x, y, z], view.W, view.H)` places 2D labels on GPU objects.
  - `colormap(name, t)` (`viridis`, `inferno`, `magma`, `coolwarm`, and cyclic `phase` for angles such as a wavefunction's phase) gives the same colours in 2D drawing and legends.
  If WebGL2 is unavailable the shell shows the reason; keep the scientific readouts in 2D so the page still informs. A page with custom shaders can still create its own canvas (see `examples/black-hole/`). Pixel self-tests: see `tests/fixtures/gl-check/page.html`.
- MathML equations: use Unicode bold letters (𝐫 𝐯 𝐄 𝐁 𝐮) for vectors, because MathML Core ignores `mathvariant="bold"`. Leading and trailing spaces inside `<mtext>` collapse: write `<mtext>&#160;on ℝ⁴</mtext>`.
- Only `psim:` and relative `./` imports are bundled. **No CDN scripts**: the deliverable must work offline. Check `core/` first (FFTs, ODE and PDE solvers, linear algebra and statistics are there); anything else, write it or add it to the project as a local module.

## Shell options and UI helpers (no need to read the ui/ source)

`createInstrument(spec)` options beyond the example above:

| Option | Meaning |
|---|---|
| `layout: 'workbench'`, `bench({ eqs, figs, createPlot, createEquation, createDerivation, M, setParam, getSim })` | Stage card on top, equations and figures below; return `{ update(sim, t, { rebuilt }) }` |
| `gl: true`, `drawGL(gl, sim, view, t)` | GPU layer under the 2D canvas (see GPU scenes above) |
| `timeLabel`, `timeUnit`, `speed`, `duration`, `timeDependent: false` | What the clock measures and its unit; `timeDependent: false` for static models (no transport bar) |
| `controls(container, { values, rebuild, setParam, setParams, play, pause, getSim })` | Add your own controls to the Parameters section. Use `setParams({ … })` to change parameters: the sliders and dropdowns show the new values and the model rebuilds once (setting `values[key]` + `rebuild()` leaves the sliders stale) |
| `extraSections(side, { section, el, row, slider, toggle, button, values, rebuild, view })` | Add whole sections to the side panel |
| `onPointer(canvas, view, getSim)`, `onFit(view)` | Pointer handling (e.g. `camera.attach(canvas)`); the transport's fit button |
| `hudLine(sim)` | Second HUD line under the time readout |
| `presets: [{ label, values: { key: value, … }, play }]` | "Try this" buttons: set several parameters at once (sliders and dropdowns show the new values; one rebuild) |
| `formatTime(t)`, `timeDigits` | HUD time text (e.g. a clock "14:30", "day 3"); default `t unit` with `timeDigits` (3) decimals; 0 for counts |
| `glResolution` (`'auto'` default, or a number ≤ 1) | GL layer pixel ratio as a fraction of the device's: `'auto'` lowers it while the GPU cannot keep up (large or fullscreen windows) and raises it again; `view.glScale` reports it. The 2D overlay (labels, text) always stays sharp |
| `frameBudget` (ms, default 9), `speedRange: [lo, hi]` | Computing time per frame; the playback-speed slider's range (default speed/100 … speed×100) |
| `hudError: false` | Keep the error diagnostic in the side panel only (pages for non-specialists) |
| `viewControls(container, { el, slider, select, toggle, button, view, getSim })` | A "View" section after Parameters for presentation-only controls (projection angle, colour map, slice): they redraw and never rebuild the model. `slider(…)` returns an element with `setValue(v)` for controls that follow an animation |
| `testSteps` (default 300) | Steps per run in `psim check`'s parameter-edit test; lower it for slow steps (Monte Carlo sweeps) |
| `texts: { parameters, paramsNote, diagnostics, presets, pending, view, viewNote }` | Section titles and the parameter note in the audience's words (`paramsNote: ''` hides it); `pending` replaces "measuring…" |

**Experimenter actions that are not parameters** (swing the masses now, fire one particle, open a shutter): call a method on `getSim()` from a `controls` or `viewControls` button, and record each action with its time in the simulation state (a schedule), so `state()`, CSV export and `advanceTo` replay it. Keep `createSim(params)` deterministic: the schedule starts empty or comes from a parameter.

**Build-time data**: expensive results that are not part of the live model (a parameter scan taking minutes) can be returned as extra fields of `verify()`; they are embedded with the report and available as `report` in `bench`, `view.report` and `window.psim.report`. Label them "computed at build time" beside the figure. Keep `verify()` to a few minutes; layout iterations reuse the cached report, and an unverified build (`build` without `--verify-module`) is fine for an early layout check, never for delivery.

**Overlays on a GPU field** whose colour map has a fixed dark end (inferno, viridis) need fixed light overlay colours in both themes; theme tokens are for overlays on the page background.

`view.theme` is `'dark'` or `'light'` (pick additive or alpha blending, mesh colours). `view.setLegend([…])` replaces the legend beside the view when the representation changes (e.g. a different projection); `model.legend` is only the initial one.

**Full screen** is built in: a button at the right end of the playback bar (a corner button on static pages) and the F key put the scene, its HUD, legend and playback bar on the whole screen (Esc leaves). Where the Fullscreen API is missing (iPhone Safari) the scene fills the window instead. Nothing to add per page; draw from `view.W`/`view.H` every frame and the view adapts (call `camera.fit` in `onFit`, which runs after the switch).

It returns `{ values, rebuild, view, sim }`. `controls` and `bench` callbacks also receive `setParam(key, value)`, `setParams({ … }, { play })`, `play()` and `pause()`; `window.psim` has the same (`window.psim.setParam`, `setParams`, `play`, `pause`, `getSim`). `view.speed` and `view.playing` let `draw` disclose the current time compression.

**Colours**: plots, `view.color()` and every GPU renderer take a theme token (`'--body-0'`, preferred: it follows the theme) or any CSS colour (`'#ff8800'`, `'rgb(…)'`, `'rgba(…)'`, names); `cssColor(c)` from `ui/gl.mjs` gives `[r, g, b, a]` in 0–1. Per-point `colors` for the point renderer are Float32 RGB in 0–1, three per point.

`psim:ui/plot.mjs` `createPlot({ title, x: { label, unit, min, max, log }, y: { … }, height, onCursor })` returns `{ el, … }`; append `el` to a figure slot. Every series has an id (adding the same id again replaces it): `fn(id, x => y, { color: '--body-0', label, dash, fill })`, `data(id, [[x, y], …], { color, label, dots: true, r, errors: [σ, …] })` (error bars from a parallel array), `band(id, x0, x1, { label })` (shaded validity range), `hline(id, y, { label })`, `vline(id, x, { label })`, `marker(id, x, y, { label })`, `remove(id)`, `setX({ min, max })`/`setY(…)`, `cursor(x)` (draggable; drags call `onCursor(x)`), `highlight(id)`, then `draw()` (`invalidate()` forces a redraw after changing data in place). Put each plot's `el` and its caption in one wrapper element: the figure grid places children in cells, so loose captions and controls end up in separate cells. A tabulated curve: `fn(id, interpolate(xs, ys))` with `interpolate` from `core/data.mjs`.

`psim:ui/instrument.mjs` helpers for those callbacks: `el(tag, attrs, ...children)` (DOM element; `attrs.html` sets innerHTML, `on…` attrs add listeners), `button(label, onClick)`, `toggle(label, checked, onChange)`, `select(options, value, onChange)` with options `[{ value, label }]`, `slider({ label, min, max, value, log, unit, digits, integer, onInput })`, `section(title, { open })`, `row(label, value)`, number formatting `sci(x, digits)` (2.4×10⁻¹¹) and `fix(x, d)`, `cssVar('--token')` for theme colours, and `download(filename, text, mime)`.

## Build and check

```bash
node "$PSIM/scripts/psim.mjs" verify-module simulations/hydrogen-2p/physics.mjs
node "$PSIM/scripts/psim.mjs" build simulations/hydrogen-2p/page.html \
     --verify-module simulations/hydrogen-2p/physics.mjs --out simulations/hydrogen-2p.html
node "$PSIM/scripts/psim.mjs" check simulations/hydrogen-2p.html
```

`check` also runs the generic interaction test: it changes the first numeric parameter through the real input (×1.3, or towards the middle of its range when that would leave [min, max]) and requires the recomputed state to differ, so make the first numeric parameter one that always acts on the model (not a value used only in one mode), and make `state()` include an observable that this parameter changes within `testSteps` steps (counts that need longer to differ make the test fail). It takes the first numeric input in the side panel: put presentation sliders in `viewControls` (after Parameters), not in an `extraSections` block above them. Do not wrap `step()` to skip work: `window.psim.steps(n)` and `advanceTo` must run the real steps (the shell's frame budget already paces playback). A deterministic calculation computed once and replayed (the same wave for every detected electron, a precomputed run) is fine: then `step()` advances the replay and the sampling, and the page says what is replayed. Look at every screenshot.

`build --verify-module` reuses the last passing report while the module, the `.mjs` files beside it and the skill's `core/` are unchanged (layout fixes then rebuild in seconds); `--reverify` forces a new run. `verify-module` always runs. If a check returns a non-finite value, every check is still printed and the bad ones are marked `BAD`.

Non-default modes (a dropdown value, a different apparatus): `psim check page.html --params mode=thomson,B=2` sets them through `window.psim.setParams` before the screenshots and tests (its screenshots go to their own folder). Check every mode the page offers, not only the default.

Timing or exploration scripts that import `psim:` modules: `node "$PSIM/scripts/psim.mjs" node scratch.mjs [args]` (measure a step's cost before choosing lattice sizes or particle counts; delete the script afterwards).
