---
name: p-simtor
description: Build interactive scientific simulations and visual explanations as self-contained HTML instruments. Use when a user wants to simulate or visualise a physical, statistical or biophysical system, including unfamiliar subjects. Translate the scientific question into a sourced mathematical model, computed observables, a suitable numerical method and an explicit visual representation; verify claims within the model's stated scope.
license: MIT
compatibility: Written for Claude Code; also follows Agent Skills conventions for Codex and OpenCode. Needs Node.js 22+ (built-in modules only). Chrome or Chromium for the browser checks.
metadata:
  version: "0.3.0"
---

# p-simtor

Turns a description of a physical system into an **interactive scientific instrument**: one offline HTML file whose behaviour is computed from the governing equations or a cited dataset, whose accuracy is measured and shown, and whose visuals reveal the physics. **The mathematics determines the behaviour; the visuals reveal it.**

`PSIM` = this skill's directory (Claude Code: `${CLAUDE_SKILL_DIR}`; otherwise the folder containing this file, e.g. `.claude/skills/p-simtor`). Run tools as `node "$PSIM/scripts/psim.mjs" <command>`. Write output into the **user's project** (`simulations/`), never into the skill.

## Shared workflow for every subject

Use this process even when the subject has no guide, recipe or solver in this repository. [references/method-selection.md](references/method-selection.md) explains the decisions in detail.

1. **Define the scientific aim.** Identify the effect to explain or quantity to predict, scale and requested fidelity. Ask one focused question if missing information changes the scientific answer. Preserve choices already made by the user.
2. **Formulate the model.** Establish the state variables, equations or dataset, units, parameters, initial/boundary conditions and assumptions from checked sources. Classify its mathematical structure, not just the subject's name. If those foundations are missing, research them before writing a solver.
3. **Define observables.** Specify how each requested quantity is calculated from the state or data, which parts are directly measurable, and the applicable uncertainty. Do not substitute an easy-to-draw quantity for the quantity the user requested.
4. **Define the time axis precisely.** For anything that evolves, decide the run's span from the question and declare it with `timeline` (see compose.md): the run covers the process asked about with a small margin and stops there. If the request states a timescale ("in 200 fs", "over a day", "10 years"), the span covers it (≲ 1.5× unless a longer window is needed; then that extra window is a labelled part of the timeline with its reason). Measure the characteristic times in the model (half-life, arrival, crossing, period) in `verify()`, mark them on the timeline as computed key times next to the stated or measured value, and report the comparison in the reply. A clock that keeps counting, a run that ends at an arbitrary time, or a stated timescale the page never shows are errors.
5. **Define the representation.** Map each visual channel to an observable or a disclosed symbol: position, height, colour, arrows, markers and playback axis. State whether the axes describe physical space, internal coordinates, probabilities or sampling. This determines the legend and conventions for any subject.
6. **Choose a feasible method.** Select an analytical, numerical, stochastic or data-based calculation from the problem's structure and accuracy needs. Reuse suitable tools; implement missing mathematics when feasible. If a reduction would change the requested scientific scope, obtain agreement before substituting it. **When you cannot ask** (a non-interactive run, or the user asked you not to): build the most representative computable version, put the scope decision in the `Fact check:` line, the page title or subtitle and the model panel, and list the alternatives you did not build.
7. **Establish evidence, then build.** Choose independent benchmarks, applicable invariants and resolution/sampling checks for the actual claims. Implement the calculation, observables and renderer separately; verify, inspect the page, and report only what the evidence supports.

Keep these decisions in the physics module's documentation and the page's model panel. A topic's absence from the repository is a reason to research its model, not to choose an unrelated recipe or stop at an unsupported animation.

## Optional shortcut: reuse a matching recipe

```bash
node "$PSIM/scripts/psim.mjs" recipes <words from the request>   # ranked verified starting points
node "$PSIM/scripts/psim.mjs" new <recipe> simulations/<name>     # copy a working, verified example
```
Recipes: `three-body`, `solar-system` (no code: one `gravity` command), `proton-magnetic-field`, `blackbody`, `black-hole` (incl. lensing), `galaxy` (dark matter), `gps-relativity`. A copied recipe **already builds and passes its checks for its original model**. For a matching request, change the needed parameters, presets, plots or terms, rebuild and check. Replacing its physics requires new verification; the original checks do not validate the replacement.

Without a matching model, use the compose contract below. Reuse a layout if helpful; copying a physics recipe is not required.

## Work economically: fewer tokens, fewer errors

1. **Reuse suitable implementations.** Check `core/` before writing a solver or numerical primitive. Verify that its variables, equations and regime fit the chosen model; a similar name or visual is insufficient. Independent derivations or reference solutions remain useful verification.
   For any general mathematics, read [references/math-toolkit.md](references/math-toolkit.md) for working imports, conditions and failure handling.
2. **No throwaway prototypes.** Put experiments directly into `verify()` and run `psim verify-module`. It stays as evidence. For a timing measurement (the cost of one step before choosing sizes), run a short script with `psim node script.mjs`, then delete it.
3. **Precision with purpose.** Choose tolerances from the observable, benchmark and relevant uncertainty (see `references/verification.md`). Numerical drift, a Monte Carlo error bar and a measured constant need different criteria. Don't chase extra digits that do not support the claim.
4. **Two strikes rule.** If a check fails twice, stop tuning and re-examine the expectation: a sign convention, a known discretisation error, an approximation's real limit. Fix the cause or report the limit honestly.
5. **Browser check once,** at the end (`psim check`), then look at the dark-desktop and phone screenshots, including the `scrolled:` ones (figures below the stage, the control column, the phone page). To see long-term behaviour, use `psim check --advance <time>` (in the page's time unit); don't write your own browser scripts. Use `verify-module` while developing.
6. **Keep implementation compact.** For a matching recipe, make focused edits. For a new model, use the shared shell, plots and equation components instead of duplicating an application. Keep physics, observables and checks in `physics.mjs` and only presentation in `page.html`; a rich 3D page may be longer, but repeated drawing code belongs in a shared renderer, not copied per page.

## Check premises and sources

Before building, check the request's physical claims, equations, units and interpretation. [references/facts.md](references/facts.md) contains examples, not a complete list of possible errors. Use checked repository data where applicable and primary sources for missing or unfamiliar claims. An ambiguous request is not automatically a false premise.

**Every reply delivering a simulation starts with a line `Fact check: …`**, before anything else:
- if the request contains a wrong premise: correct it briefly and explain the agreed scope of the result;
- otherwise: state the key facts or modelling assumptions the simulation relies on, with sources.

Never silently drop a wrong part of a request, and never invent objects as if they were real.

## Non-negotiable rules

1. **Computed, not choreographed.** Every quantitative claim comes from the model's equations or a cited dataset. Symbolic drawings and analogies must be identified as such. No keyframes, splines or easing passed off as computed physical motion.
2. **Rendering never touches physics;** simulated time is independent of frame rate (the shells handle this).
3. **Physical parameters change the calculation.** Keep presentation controls such as palette or playback speed separate from changes to the model's state, equations or conditions.
4. **Measure before you claim:** numbers for conservation, convergence and benchmarks, not adjectives.
5. **Label the trust tier for the stated model.** Numerical verification establishes solver accuracy; physical validation requires an appropriate independent benchmark. A checked toy model is not validation of the full theory.
6. **Respect the mathematical meaning of the state.** Positions, fields, amplitudes, distributions and statistical samples are different objects. Do not turn one into another without a stated physical approximation or measurement operation.
7. **Explain what the viewer is seeing.** Put the meaning of colour, geometry and the playback axis beside the view, with details in the model panel. Optical colour, a symbolic charge label and a heat-map value are different encodings. Monte Carlo sweeps and parameter scans are not physical time.
8. **Answer the scientific question with a defensible model.** Prefer a verified computation when it serves that question; use a cited dataset or a labelled analogy when appropriate. Name any reduction that changes the meaning in the title/subtitle or adjacent caption, and explain its limits in the model panel and final reply.

| Tier | Use when | Evidence |
|---|---|---|
| `validated` | Verified numerical solution of a specified established model | Appropriate invariants + refinement + an exact or published benchmark; name the validated scope |
| `approximate` | Reduced, effective or deliberately simplified model | Model-specific checks, with omissions and regime of validity named |
| `analytical` | Closed form of the stated model | Identity / normalisation / limiting-case checks; exact evaluation does not establish physical applicability |
| `conceptual` | Explanatory diagram or analogy | Sources and explicit symbolic meaning; no numerical prediction claimed |

Mixed pages must label computed results, imported data and illustrative overlays individually. A page-wide badge must not lend numerical credibility to an uncomputed overlay.

## Implement the selected model

1. Use the contract, toolkit and workbench APIs in [references/compose.md](references/compose.md). Optional domain notes are linked from method selection; they accelerate research but do not define the supported subjects.
2. `physics.mjs`: implement `PARAMS`, `createSim(params)` and `verify()` for the selected model and observables. Reuse a matching recipe or create the module directly. Use `core/constants.mjs` for available constants; cite the values, units and provenance of additional inputs.
3. `page.html`: use `createInstrument({ … })`, a renderer appropriate to the state, and equations split into terms (`createEquation`). Derive legends and conventions from the representation decisions, using [references/visual-language.md](references/visual-language.md).

```bash
node "$PSIM/scripts/psim.mjs" verify-module simulations/x/physics.mjs
node "$PSIM/scripts/psim.mjs" build simulations/x/page.html --verify-module simulations/x/physics.mjs --out simulations/x.html
node "$PSIM/scripts/psim.mjs" check simulations/x.html
```
Gravity needs no code: `psim gravity --preset <name> | --scenario file.json` (format in `references/gravity.md`). Exit code 2 from any command means a check failed: fix it, don't ship it. If Chrome is blocked (sandboxed agents), say the browser check was not run.

## Deliver

Start with the `Fact check:` line (see above). Keep the rest short: the model and method, the trust tier, the measured verification numbers, and one line of limits.

**A final message delivering or updating a runnable page must END with an "Open it" block:**
```
**Open it:** [<name>.html](file:///absolute/path/to/simulations/<name>.html)
or run: open "/absolute/path/to/simulations/<name>.html"
Try: <two or three things to try>. Want me to open it now?
```
Copy the link and command from the build's "Open it:" output (it prints the right command for macOS, Linux or Windows). The page works offline and needs no server. If the user says yes, run `psim open <page.html>`.

Clarifying questions and discussions before a build need no delivery block. Never invent an artifact, verification result or open link.

## Files

- `scripts/psim.mjs`: `recipes`, `new`, `gravity`, `presets`, `build`, `verify-module`, `verify`, `check`, `open`, `node` (scratch scripts with `psim:` imports)
- `core/`:
  - `numerics` (rng, fft, rk4Step, dopri5, gaussLegendre, adaptiveSimpson, integrateToInfinity, brent, lgamma, gammaP, gammaQ, erf, besselI, besselJ, besselJZeros, laguerre, complexYlm, hydrogenRadial, legendreP, realYlm, solveTridiagonal, eigenTridiagonal, linearFit, meanWithError)
  - `constants` (CODATA 2018, IAU, μ_B, g_e, α mass; `PDG` particle data with uncertainties and edition; IUPAC `ATOMIC_WEIGHT`)
  - `matter` (Bethe stopping power, ranges, PDG material table), `spin` (spin ½ / qubits: rotation, Born rule, measurement)
  - `quantum` (complex matrices, spin operators for any j, tensor products, Hermitian eigensolver, exact propagators, Lindblad open-system dynamics)
  - `linalg` (matrixVector, solveLinear, leastSquares, eigenSymmetric, cholesky, eigenGeneralized, expm, conjugateGradient, newtonSolve)
  - `fields` (laplacianND, gradientND, divergenceND, fftND, solvePoissonPeriodic, diffusionPeriodic, schrodingerND with absorbing edges, probabilityCurrent, contourLines, densityLevel; optional domain masks)
  - `probability` (maxwellSpeed, effusiveSpeed, logSumExp, normalizeWeights, normalizeLogWeights, shannonEntropy, relativeEntropy, covarianceMatrix, canonicalEnsemble)
  - `ode` (odeStepper: rk4/dopri5 with dense output; mechanicsStepper: leapfrog/yoshida4; refinementOrder, lyapunovExponent, linearChain)
  - `particles` (createParticles with optional periodic box; gravity, coulomb, lennardJones, uniformField, centralBody, drag, springs, externalPotential; forEachPairWithin)
  - `pde` (diffusion1D with held/flux walls and sources, waveND, stationaryStates1D)
  - `stochastic` (gillespie, eulerMaruyama, langevinBAOAB with heat, metropolis, metropolisSites, autocorrelationTime, effectiveSampleSize)
  - `data` (interpolate, histogram, createHistogram, ksTest, chiSquareTest, jackknife, jackknifeRefit, bootstrap, curveFit, meanSem, pull)
  - `kinematics` (4-vectors, boosts, random directions and rotations, twoBodyDecay, rambo phase space, acceptReject)
  - `nbody`, `nbody-verify`, `kepler`, `presets`, `solar-system`
  - `radiation` (Planck, Wien, σ, colour)
  - `schwarzschild` (light: tracePhoton, deflectionExact, shadowAngle; orbits with Einstein's correction: precessionPerOrbit, precessionSI, relativisticOrbit; ISCO, diskRedshift)
  - `relativity` (clock rates, GPS)
  - `sun` (the Sun's position in the sky, sunrise/sunset, equation of time)
  - `galaxy` (bulge/disk/NFW potentials, rotation curves, star disks)
  - `lattice-gauge` (SU(2) Monte Carlo, wilsonLoops, polyakovLoop, plaquetteField, apeSmear, strongCouplingPlaquette, gauge transforms, centre symmetry, Landau gauge, topology and cooling, BPST instanton benchmark)
  - `projection` (4D rotations and projections, S³ maps: stereographic, axis–angle, Hopf)
  - `qcd` (running coupling α_s: β coefficients, 1–4-loop MS-bar running with flavour thresholds)
- `ui/`: `instrument-app.mjs` (shell; `stage` or `workbench` layout; `gl: true` adds a GPU layer), `gl.mjs` (WebGL2 points, 2D fields, 3D volumes, spheres, lines, solid meshes with shadows and arrows, orbit camera, colour maps and colour bars, line shapes), `plot.mjs`, `equation.mjs`, `instrument.mjs`, `orbit-view.mjs`, `instrument.css`
- `examples/`: the recipes (`catalog.json`), each with `physics.mjs`, `page.html` and a built page
