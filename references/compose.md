# Composing a simulation for any domain

Use this when no built-in template fits. The reference implementation is `examples/proton-magnetic-field/` (`physics.mjs`, `page.html` and the built page). Copy its structure.

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
  { key: 'state', label: 'Orbital', value: '2p_z', options: [{ value: '2p_z', label: '2p (m = 0)' }] },
];

export function createSim(params) {
  return {
    t,                       // getter: physical time in the page's timeUnit
    step(),                  // one integrator step (time-dependent models)
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
- `createSim` must be deterministic for given params. Use a seeded RNG (e.g. mulberry32) for stochastic models and expose the seed as a param.
- Diagnostics of `kind: 'error'` are relative errors (plotted on a log axis; the first one also appears in the HUD). Use `good`/`warn` thresholds from `verification.md`.
- `verify()` compares against **exact or published** results, not against the code itself. Include at least one convergence check for numerical methods.
- Static (time-independent) visualisations: pass `timeDependent: false` to `createInstrument`. Then `step()` is not needed.

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
  timeUnit: 'ns', speed: 30, duration: 1000,
  hudLine: sim => '…',
  model: { legend: [], equations: [{ caption, mathml }], assumptions: [], conventions: [], references: [] },
});
</script>
```

- `draw` is called every frame with a cleared, DPR-scaled 2D context. It must only **read** `sim`.
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
