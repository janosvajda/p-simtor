---
name: p-simtor
description: Build scientifically verified, beautifully crafted interactive physics simulations as self-contained HTML instruments - planets and the Solar System, three-body chaos, galaxies and dark matter, black holes and gravitational lensing, relativity and GPS clocks, charged particles in fields, black-body radiation, lattice gauge theory (quark confinement), molecular dynamics and quantum wavefunctions. Use whenever someone asks to simulate, model, visualise, animate or explain a physical system or phenomenon, even if they never say "simulation". The skill computes the evolution from the real equations, measures accuracy against exact results and real measurements, labels each result's trust tier, and never animates fabricated physics.
license: MIT
compatibility: Written for Claude Code; also follows Agent Skills conventions for Codex and OpenCode. Needs Node.js 22+ (built-in modules only). Chrome or Chromium for the browser checks.
metadata:
  version: "0.3.0"
---

# p-simtor

Turns a description of a physical system into an **interactive scientific instrument**: one offline HTML file whose motion is computed from the governing equations, whose accuracy is measured and shown, and whose visuals reveal the physics. **The mathematics determines the behaviour; the visuals reveal it.**

`PSIM` = this skill's directory (Claude Code: `${CLAUDE_SKILL_DIR}`; otherwise the folder containing this file, e.g. `.claude/skills/p-simtor`). Run tools as `node "$PSIM/scripts/psim.mjs" <command>`. Write output into the **user's project** (`simulations/`), never into the skill.

## Fast path: start from a verified recipe (do this first)

```bash
node "$PSIM/scripts/psim.mjs" recipes <words from the request>   # ranked verified starting points
node "$PSIM/scripts/psim.mjs" new <recipe> simulations/<name>     # copy a working, verified example
```
Recipes: `three-body`, `solar-system` (no code: one `gravity` command), `proton-magnetic-field`, `blackbody`, `black-hole` (incl. lensing), `galaxy` (dark matter), `gps-relativity`. A copied recipe **already builds and passes its checks**. Change only what the request needs (parameters, presets, one extra plot or term), rebuild, check, deliver. Most requests need tens of lines of edits, not a new 900-line page.

No close recipe? Compose (below), starting from the recipe with the most similar *layout*.

## Work economically: fewer tokens, fewer errors

1. **Reuse before writing.** Before writing any formula or solver, check `core/` (its exports are listed under Files). Never re-derive what exists: deflection, shadow, Planck, Kepler, N-body, orbits, clock rates, galaxy potentials, RNG, FFT, ODE solvers, quadrature, roots, special functions.
2. **No throwaway prototypes.** Put experiments directly into `verify()` and run `psim verify-module`. It's fast and it stays as evidence.
3. **Precision with purpose.** Verify to the precision that matters physically and for the claim (see thresholds in `references/verification.md`). 1e-8 relative is excellent for almost everything. Don't chase 1e-15.
4. **Two strikes rule.** If a check fails twice, stop tuning and re-examine the expectation: a sign convention, a known discretisation error, an approximation's real limit. Fix the cause or report the limit honestly.
5. **Browser check once,** at the end (`psim check`), then look at the dark-desktop and phone screenshots. To see long-term behaviour, use `psim check --advance <time>` (in the page's time unit); don't write your own browser scripts. Use `verify-module` while developing.
6. **Small edits over rewrites:** edit the copied files in place. Keep `page.html` under ~400 lines by using `createPlot`, `createEquation`, `createDerivation` and the `workbench` layout.

## First: fact-check the request

Before building, read [references/facts.md](references/facts.md) and check every object, number and claim in the request (for example, Mercury has **no moons**; electrons don't orbit like planets). Use checked data from `core/` where it exists (planets and their moon counts: `core/solar-system.mjs`).

**Every final reply starts with a line `Fact check: …`**, before anything else:
- if the request contains a wrong premise: correct it in one kind sentence and say what you built instead ("Mercury has no moons, so I simulated Mercury alone; ask if you want a hypothetical moon");
- otherwise: state the key facts the simulation relies on, with sources ("Mercury's orbit from JPL J2000 elements; relativistic advance 42.98″/century").

Never silently drop a wrong part of a request, and never invent objects as if they were real.

## Non-negotiable rules

1. **Computed, not choreographed.** Every displayed value comes from the model's equations or a cited dataset. No keyframes, splines or easing passed off as physics.
2. **Rendering never touches physics;** simulated time is independent of frame rate (the shells handle this).
3. **Parameters change the calculation,** never just speed or appearance.
4. **Measure before you claim:** numbers for conservation, convergence and benchmarks, not adjectives.
5. **Label the trust tier** with the reason; put assumptions, limits and visual conventions in the page.
6. **Right mathematics for the scale:** quantum results are probability densities, never planetary electrons; QCD is never computed quark trajectories.
7. **Disclose visual conventions** (exaggerated sizes, symbolic colours, time compression).
8. **Compute the strongest real model that fits** (lattice Monte Carlo over a fitted curve, ray-traced geodesics over a drawn black hole). State limits in one line at the end, never as the opening.

| Tier | Use when | Evidence |
|---|---|---|
| `validated` | Numerical solution of an established model | Conservation + convergence + an exact or published benchmark |
| `approximate` | Numerical with modelling approximations (softening, test particles, empirical potentials) | Same, plus the approximation named |
| `analytical` | Closed form evaluated exactly | Identity / normalisation checks |
| `conceptual` | Nothing computable here | Every element cited, marked illustrative |

## Compose (no close recipe)

1. Read [references/method-selection.md](references/method-selection.md) and the domain guide: [gravity](references/gravity.md), [electromagnetism](references/electromagnetism.md), [radiation & relativity](references/radiation-relativity.md), [particle physics](references/particle-physics.md), [molecular](references/molecular.md), [quantum](references/quantum.md).
2. Copy the most similar recipe with `psim new`, then replace its physics. The contract, the toolkit table and the workbench/self-test APIs are in [references/compose.md](references/compose.md).
3. `physics.mjs`: `PARAMS`, `createSim(params)`, and `verify()` against exact results or real measurements. Constants come from `core/constants.mjs`.
4. `page.html`: `createInstrument({ … })` with a domain-appropriate `draw()`; equations split into terms (`createEquation`). Visual rules: [references/visual-language.md](references/visual-language.md).

```bash
node "$PSIM/scripts/psim.mjs" verify-module simulations/x/physics.mjs
node "$PSIM/scripts/psim.mjs" build simulations/x/page.html --verify-module simulations/x/physics.mjs --out simulations/x.html
node "$PSIM/scripts/psim.mjs" check simulations/x.html
```
Gravity needs no code: `psim gravity --preset <name> | --scenario file.json` (format in `references/gravity.md`). Exit code 2 from any command means a check failed: fix it, don't ship it. If Chrome is blocked (sandboxed agents), say the browser check was not run.

## Deliver

Start with the `Fact check:` line (see above). Keep the rest short: the model and method, the trust tier, the measured verification numbers, and one line of limits.

**The final message must END with an "Open it" block. These are its last lines, every time, even after follow-up edits:**
```
**Open it:** [<name>.html](file:///absolute/path/to/simulations/<name>.html)
or run: open "/absolute/path/to/simulations/<name>.html"
Try: <two or three things to try>. Want me to open it now?
```
Copy the link and command from the build's "Open it:" output (it prints the right command for macOS, Linux or Windows). The page works offline and needs no server. If the user says yes, run `psim open <page.html>`.

## Files

- `scripts/psim.mjs`: `recipes`, `new`, `gravity`, `presets`, `build`, `verify-module`, `verify`, `check`, `open`
- `core/`:
  - `numerics` (rng, fft, rk4Step, dopri5, gaussLegendre, adaptiveSimpson, integrateToInfinity, brent, lgamma, erf, besselI, laguerre, legendreP, realYlm, solveTridiagonal, linearFit, meanWithError)
  - `constants` (CODATA 2018, IAU)
  - `nbody`, `nbody-verify`, `kepler`, `presets`, `solar-system`
  - `radiation` (Planck, Wien, σ, colour)
  - `schwarzschild` (light: tracePhoton, deflectionExact, shadowAngle; orbits with Einstein's correction: precessionPerOrbit, precessionSI, relativisticOrbit; ISCO, diskRedshift)
  - `relativity` (clock rates, GPS)
  - `galaxy` (bulge/disk/NFW potentials, rotation curves, star disks)
  - `lattice-gauge` (SU(2) Monte Carlo, Wilson loops)
- `ui/`: `instrument-app.mjs` (shell; `stage` or `workbench` layout), `plot.mjs`, `equation.mjs`, `instrument.mjs`, `orbit-view.mjs`, `instrument.css`
- `examples/`: the recipes (`catalog.json`), each with `physics.mjs`, `page.html` and a built page
