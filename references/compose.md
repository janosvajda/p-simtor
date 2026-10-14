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
    diagnostics(),           // [{ key, label, value, kind: 'error'|'value', unit?, good?, warn?, text? }]
    state(),                 // JSON-serialisable state (export + automated checks)
    record?(),               // { header, rows } for CSV export
    // plus anything draw() needs (trail, grids, analytic reference values…)
  };
}

export function verify() {   // run by psim verify-module / build --verify-module
  return { summary, checks: [{ label, value, pass }], notes };
}
```

Rules:
- Keep calculation, observation and representation distinct. The calculation updates/samples the model state; observation functions derive the requested quantities; `draw()` maps these results to the declared visual channels. Document units, conditions, provenance and the progress axis alongside the model. A new subject uses this same contract without needing a new recipe or domain guide.
- `createSim` must be deterministic for given params. Use a seeded RNG (e.g. mulberry32) for stochastic models and expose the seed as a param.
- Diagnostics of `kind: 'error'` are relative errors (plotted on a log axis; the first one also appears in the HUD). Use `good`/`warn` thresholds from `verification.md`.
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
| Lattice gauge theory | `psim:core/lattice-gauge.mjs` |

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
- **GPU or other self-tests**: define `window.psim.selfTest = () => ({ label, pass, …measurements })` after `createInstrument`. `psim check` runs it (see `examples/black-hole/page.html`).
- For WebGL (volumetric densities, large particle counts), create your own canvas inside `draw`'s first call, or use `extraSections`. Keep a 2D-canvas fallback message if WebGL is unavailable.
- MathML equations: use Unicode bold letters (𝐫 𝐯 𝐄 𝐁 𝐮) for vectors, because MathML Core ignores `mathvariant="bold"`.
- Only `psim:` and relative `./` imports are bundled. **No CDN scripts**: the deliverable must work offline. If you need a library (for example an FFT), write it or add it to the project as a local module.

## Build and check

```bash
node "$PSIM/scripts/psim.mjs" verify-module simulations/hydrogen-2p/physics.mjs
node "$PSIM/scripts/psim.mjs" build simulations/hydrogen-2p/page.html \
     --verify-module simulations/hydrogen-2p/physics.mjs --out simulations/hydrogen-2p.html
node "$PSIM/scripts/psim.mjs" check simulations/hydrogen-2p.html
```

`check` also runs the generic interaction test: it edits the first numeric parameter through the real input and requires the recomputed state to differ. Look at every screenshot.
