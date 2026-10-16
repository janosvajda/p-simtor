# Changelog

## Unreleased — lessons from a generated black-hole page

- **Fixes for problems three newly generated astrophysics pages ran into:**
  - **Action results are now shown.** Each action can declare `reveal`: the seconds of playback until its result is plainly visible (default 8). `psim check` saves `-action-N.png` 2 s after the press and `-action-N-result.png` at `reveal`, and the newcomer test looks at both. Actions also get `seek(t)`, to jump to just before the moment they are about; pages had been relying on the undocumented `window.psim.advanceTo`.
  - **The starter's main view is a picture.** It shows 400 atoms going dark one by one, and the decay curve sits below it as evidence. Before, the starter was a chart, which taught the opposite of the rule. Its unused "atoms at the start" setting is gone.
  - **Key times outside the run are allowed.** A measured event after the point where the model stops being valid (for example, impacts after a computed break-up) is listed in the side panel as "after this run". Before, it threw an error, which pushed pages to extend the run beyond what the model covers.
  - **`psim node SCRIPT.mjs [args]`** forwards every argument after `node`, and the script reads them as `process.argv.slice(2)`. This is now documented and tested.

- **Two bugs found by newly generated pages:**
  - The shell accepts a verification `notes` given as a string. The starter returned one, which would have printed one paragraph per character.
  - `psim check` treats a page whose clock starts before zero (a negative `timeline.start`) as running once the clock passes its start. Before, such a page was reported as not running, with no reason given; the summary now gives one.

A page generated from Kip Thorne's colloquium "Geometrodynamics: The Nonlinear Dynamics of Curved Spacetime" exposed four gaps in the skill:

- **`psim check` now presses every main action button** (`spec.actions`). For each button it:
  - starts from 40 % through the run;
  - plays for 2 s and saves a screenshot (`-action-N.png`);
  - fails on an error or a frozen page;
  - flags a press that changes nothing.

  Until now the page's main demonstration was never exercised, and no screenshot showed its effect. The newcomer test now includes these shots.
- **SKILL.md step 7: what "most representative" means** when the full request cannot be computed and the user cannot be asked, in order:
  1. The same system and headline phenomenon.
  2. The largest stage of it that can be computed from the real equations and verified.
  3. A run that stops where the model stops being valid, with a marker naming what is not computed.

  Any missing stage may be shown only as labelled published data. Evaluating an exact closed-form solution counts as computed; a reduced model solved exactly is still `approximate`.
- **The convergence-order rule no longer misfires** (references/verification.md):
  - Orders are measured where errors are well above round-off.
  - Too low fails.
  - Higher than theory is accepted when the reason is stated.
  - Residuals of discretised equations are judged by their order under refinement, not by an absolute threshold.
- **`core/gravitational-waves.mjs`** ([gravitational-waves.md](references/gravitational-waves.md)):
  - circular inspiral driven by the radiated energy (exact solution of the leading-order equations);
  - strain h₊ and h×;
  - the linearised tidal field E_ij of the binary in the near, transition and wave zones;
  - the energy flux of that field;
  - chirp and binary-pulsar formulas.

  The guide says which stages of a black-hole collision a browser can compute (inspiral, field) and which it cannot (merger, ring-down). Tested against:
  - the Newtonian limit;
  - the vacuum field equations, with the residual's order under refinement;
  - the quadrupole formula, to 5e-9;
  - the Hulse–Taylor period decay;
  - GW150914's chirp mass.

## Unreleased — one idea per page, designed before the code

These skill changes answer the owner's review of two reworked pages. In one, the main button sat under a long text in the side panel. In the other, a label jumped from frame to frame, and the page as a whole explained nothing a newcomer could follow.

- **SKILL.md step 5 (new): design the page around one idea, before any code.** It names:
  - the one thing a newcomer should learn;
  - the one picture that shows it;
  - the one action that demonstrates it.

  Everything else goes below the view or into the side panel.
- **Step 8: the newcomer test.** It runs on the agent's own first-screen screenshots, desktop and phone:
  - Does the first screen show the one idea?
  - Is the main action visible?
  - Are the labels steady?
  - Is anything unexplained?

  Any "no" means the page is not finished.
- **`psim new --starter <id>`.** A blank page with every required part already in place, which builds and passes every check (`templates/starter/`): the one-idea header, intro, main action, timeline with `seconds`, legend, and a full-view scene with a steady label. The agent replaces the model and the drawing instead of assembling the page from the docs.
- **`actions`.** The page's main action ("Heat it up") is a large button directly under the view. `psim check` fails a page whose action is below the fold at 1440×900.
- **`view.steady(id, x, y)`.** A label attached to something moving follows it smoothly and ignores frame-to-frame jitter.
- **Fixes from the first two fresh generations' audits:**
  - "As it heats up" questions are scanned over the control parameter; the skill no longer invents a heating time.
  - Labels must be true at every moment and use everyday words; a symbol appears only with its plain name.
  - Colour bars use plain words, sit small at an edge of the view, and don't count towards the label limit.
  - quantum.md gives grid-resolution advice for sharp walls (≥ 15–20 points per wavelength, checked by refinement).
  - `psim check` saves a `-midrun` screenshot, which the newcomer test also judges.
  - The phone view leaves room for the player and action rows.
  - visual-language.md describes the current layout.
  - `--preset` appears in the help.
- **Fixes from the third fresh generation's audit:**
  - The starter's folder option is explained where the starter is introduced, and equations or the workbench layout are clearly optional.
  - Scans mark characteristic parameter values as their key points.
  - The newcomer test requires the intro and the "look for" items to be true in every preset.
  - The diagnostics plot names the real progress axis (not "t").
  - `hudLine(sim, t)` gets the displayed time, so the readout and the picture show the same moment.
- **The intro shows two lines with "more".** The settings below it stay in reach.

## Unreleased — the animation first, explained in plain words

The project owner's review found these problems:
- Playback controls sat far below the animation.
- Text covered the scenes.
- There was no explanation a non-physicist could follow.
- Runs were too short.
- Some pages were slow.

The shell layout is redesigned for every page:

- **One scene card:**
  - The view has nothing from the shell drawn over it, apart from a small "Finished · ↺ Replay" pill at the end of a run.
  - Directly under the view sits a video-player bar: restart, play, step, the clock, the timeline with its key moments, speed, fit and full screen.
  - Below that, a key row holds the legend and one short live readout.
  - On a phone, the first screen holds the view and the whole player bar.
- **`intro: { text, watch }` (required).** "What you are seeing" explains in everyday words what is shown and what to look for. It is the first item in the side panel on wide screens, and comes right under the player on phones. It is the same on every page, including the gravity template.
- **Error figures** no longer appear under the animation (`hudError` defaults to false); they stay in Diagnostics.
- **The gravity template** uses the same layout and gains full screen. F now toggles full screen there too, as in the app shell.
- **New `psim check` rows:**
  - Plain-language intro: required, 12–90 words.
  - Playback length at default speed: fails under 30 s; 45–120 s is recommended.
  - Smooth playback in frames per second at 1440×900: fails under 24 fps.
- **SKILL.md rules 4 and 7, compose.md and visual-language.md** cover the layout, the intro, watchable run lengths, at most four short labels on the canvas, and a large centred scene.
- **The animation view is a dark screen in both page themes.** Panels, text and figures follow the light/dark toggle. `view.color` and the gravity view read the view's dark tokens; `view.theme` is always `'dark'`, and `view.pageTheme` gives the page theme. Before this change, pages that painted their own dark scenes showed black boxes in light pages.
- **`core/disk-nbody.mjs`** ([spiral-galaxies.md](references/spiral-galaxies.md)): a self-gravitating star disk in a rigid halo, using an isolated particle-mesh solver (zero-padded FFT, CIC).
  - Initial conditions are set by Toomre Q.
  - Measurements: Fourier arm amplitudes, pitch angle, Q(R) and the rotation curve.
  - Checked against exact Kuzmin and Freeman forces, momentum and L_z conservation, and a hot disk that stays smooth.
  - A cool disk (Q = 1.2) grows trailing three-armed spirals, against a swing-amplification estimate of 2.7.
- **The GL test fixture** pins `glResolution: 1`, so its pixel checks are exact.

## Unreleased — seven new demos and shell fixes they exposed

- **Gallery:** seven demos in a new "Waves, flows and the cosmos" section. They were generated from prompts by the installed skill in the test project, using the new modules:
  - Kármán vortex street
  - Dzhanibekov effect
  - age of the universe
  - soliton collision
  - laser cavity
  - Chandrasekhar limit
  - dam break
- **`--advance T` holds at T.** `psim check` and the gallery thumbnail now pause at T. Before, playback continued, and screenshots showed a later time.
- **Timeline labels are placed by measured width** in up to three rows. A label that fits nowhere is hidden; its tick keeps the tooltip, and the side panel lists every key time. The start label is dropped when it would touch the end label. This fixes label collisions on phones.
- **`view.insets` and `view.overlays`.** Pages can now see how much of the stage the HUD, legend, timeline and playback bar cover, instead of guessing margins.
- **`drawColorBar`** restores the canvas state and keeps its label clear of the top tick.
- **Docs.**
  - Equation `link` belongs on `createEquation`, and term values may be strings.
  - Timeline label rules, and how to put run-measured times on the timeline.
  - The `linearFit` return shape.
  - Vortex-street lessons in fluids.md: inlet distance, side boundaries, TRT vs BGK, ramped start, browser throughput.
  - Source terms by operator splitting, and KdV step size for accuracy, in nonlinear-waves.md.

## Unreleased — scientific accuracy review

- Fixed cancelling hypergeometric polynomials using exact rational summation of the input doubles, with a correctly rounded binary64 result. Added scaled Bessel K recurrence, protected exponential tails, a zeta expansion near zero and overflow-safe Student-t tails; unreliable near-integer polylog expansions now reject explicitly.
- Replaced characteristic-polynomial eigenvalue estimates with the shared scaled Hessenberg/Francis QR solver. Stiff integration now preserves reused callback outputs, forms meaningful finite-difference Jacobians at large state scales and respects fixed-step `hmax`.
- Corrected Nyquist contributions to wave energies/momentum, real-field odd derivatives, and initial-state projection for quadratic dealiasing.
- Stabilised opaque-film optics, large grating intensities, Fresnel phases, small robust-fit residuals, and constrained fitting. Corrected Saha underflow and dilute-gas root brackets, small-redshift cosmological integrals, symmetric-top limits and nearly null geometric-algebra inverses.
- Added independent regression checks and a bundled-calculation comparison. Updated numerical accuracy, uncertainty and stability guidance; physical model approximations remain explicitly identified. See [the review](docs/scientific-audit.md).

## Unreleased — new physics and mathematics

Seven new subject areas in `core/`, each with tests against exact solutions or published values and a reference guide linked from method selection (514 Node tests in total).

- **Waves and optics** (`optics`, [optics.md](references/optics.md)):
  - Fresnel integrals and closed-form diffraction (slits, gratings, Airy disc, knife edge).
  - FFT far field and exact angular-spectrum propagation (checked against the Talbot effect and Gaussian beam widths).
  - Gaussian beams, ABCD matrices and cavity modes.
  - 3D ray tracing, Fresnel coefficients, Sellmeier glasses and rainbows.
  - Thin-film stacks, Jones calculus.
  - A 2D Yee FDTD solver with exact discrete dispersion and energy checks.
- **Incompressible fluids** (`fluids`, [fluids.md](references/fluids.md)):
  - D2Q9 lattice Boltzmann with obstacles, inflow/outflow and momentum-exchange forces, and a periodic spectral vorticity solver.
  - Benchmarks: Poiseuille and Couette exact, Taylor–Green decay, the Ghia lid-driven cavity, and the Kármán street Strouhal number.
- **Nonlinear waves and conservation laws** (`spectral`, `conservation`, [nonlinear-waves.md](references/nonlinear-waves.md)):
  - Chebyshev collocation, ETDRK4 and split-step solvers for KdV, NLS, Kuramoto–Sivashinsky, Burgers and sine-Gordon, with exact solitons.
  - Finite-volume Euler, shallow-water and Burgers solvers, checked against exact Riemann solutions (Toro) and the Ritter/Stoker dam breaks.
- **Rigid bodies and chaos** (`quaternion`, `rigid`, `chaos`, [rigid-bodies.md](references/rigid-bodies.md), [chaos.md](references/chaos.md)):
  - Exact torque-free motion, a Lie–Poisson free body and heavy top, and the intermediate-axis flip period.
  - Feigenbaum δ, Lyapunov spectra (Lorenz sum exact), Poincaré sections, correlation and box dimensions.
- **Cosmology and stars** (`cosmology`, `stars`, [cosmology.md](references/cosmology.md), [stars.md](references/stars.md)):
  - Planck 2018 ΛCDM: age 13.786 Gyr, distances matching astropy, growth factor, sound horizon.
  - Lane–Emden polytropes, Chandrasekhar white dwarfs, the mass formula against AME2020, Saha ionisation, stellar timescales.
- **Stiff ODEs and special functions** (`stiff`, `special`, [stiff-and-special.md](references/stiff-and-special.md)):
  - RADAU5, Rodas3 and variable-order BDF with the `odeStepper` interface, checked on Robertson, Van der Pol, HIRES and Oregonator.
  - Airy, Bessel Y/K and spherical, Hermite, hypergeometric, elliptic, zeta and polylog, incomplete beta and Student t, and exact Wigner 3j/6j and Clebsch–Gordan.
- **Optimisation and geometry** (`optimize`, `autodiff`, `qmc`, `geometric-algebra`, [optimization.md](references/optimization.md)):
  - BFGS/L-BFGS, Nelder–Mead, bounded and robust Levenberg–Marquardt, differential evolution.
  - Forward and reverse automatic differentiation.
  - Sobol and Halton sequences with scrambling.
  - Clifford algebras with rotors, boosts and field invariants.
- `rigid.jacobiElliptic` now re-exports the general `special.jacobiElliptic`, which also covers m > 1.

## Unreleased — precise timers

- **`timeline` in the shell.** Every evolving page declares its run span and key times (computed, measured, requested). A bar under the time readout shows progress and labelled markers, and clicking it jumps to a time. The run stops exactly at its declared end.
- **`psim check`.**
  - It advances to 3× the declared duration and requires the clock to stop exactly at the end, with every marker inside the run.
  - Time-dependent pages without a timeline get a visible warning.
- **SKILL.md workflow step 4: define the time axis.** The span follows the question, any longer window is labelled with its reason, and characteristic times are measured in `verify()` and compared with the stated value.
- **Checklists.** `verification.md` adds timescale claims and `visual-language.md` adds a timeline question.
- **Shell robustness.**
  - Timelines may start at a non-zero time, such as a logarithmic clock starting at −1.
  - The end time snaps to exactly the declared value.
  - Labels near the edges of the bar are anchored, so they are not clipped.
  - The bar has a backing on every page, and on phones it shares the bottom row with the legend without overlapping.
  - Counted progress (sweeps, steps, generations) shows whole numbers.
  - The "finished" text is shorter.
- **GPU self-tests.**
  - `window.psim.nextFrame()` resolves in the task that drew a frame, so `gl.readPixels` sees the drawn buffer.
  - `psim check` calls the self-test after a drawn frame and repeats a failed one once. A pass on the second attempt is reported as such, with the first reading.
  - The Rutherford, Stern–Gerlach and weak-force gallery pages now use `nextFrame`. Rutherford had read the buffer after compositing, which gave blank reads.
- **Self-tests.** The enzyme page's rectangular-barrier mode and the dark-matter page's halo-off mode now have correct self-tests. Both were found by `--params` checks.
- **Recipes.** The galaxy, GPS-relativity and proton-in-a-field recipes declare their timelines, with computed key times: an orbit at R₀ is 223 Myr, a GPS orbit 11.967 h, a gyroperiod 65.66 ns.
- **Gallery.** Every time-dependent gallery demo declares its timeline. The retinal (vision) demo now runs to 300 fs by default, with S₁'s half-life of 138 fs (computed) and photorhodopsin at about 200 fs (measured) marked; a View option extends the run to 1000 fs for the long-time yield.

## Unreleased — quantum biology builds

- **New `core/quantum.mjs`** for few-level and open quantum systems:
  - complex matrices, spin operators for any j, tensor products;
  - a degenerate-safe complex Hermitian eigensolver, and exact propagators;
  - Lindblad dynamics, thermal states and density-matrix checks.
  - Tested against exact Larmor precession, Heisenberg spectra, e^{−γt} dephasing, e^{−Γt} decay and detailed balance.
- **GPU:**
  - meshes coloured per vertex by a computed quantity (`values`, `colormap`, `range`), with a pixel test;
  - `geometry.surface` for parametric surfaces (polar plots, energy surfaces);
  - `camera.ray` for picking and `camera.fitPoints` for tight framing;
  - per-point sizes in the point renderer.
- **`psim check`:**
  - an option and its value passed as one argument is now an error, not silently ignored;
  - `--params` also applies to the self-test and the edit test;
  - `--advance` and `--params` screenshots never overwrite the default run's;
  - "steps —" is gone from the summary.
- **Phone fix:** derivation lines scroll on phones instead of being clipped.
- **From the vision build:**
  - **Verification cache fix.** The cache no longer stores a stale key when files change during a long `verify()`. This caused repeated 10-minute re-verifications.
  - `psim check --preset "label"` tests a "Try this" preset.
  - New `coupledSchrodinger2`, for two coupled surfaces (Rabi oscillations and the uncoupled limit tested).
  - `schrodingerND`: per-axis masses and `energy()`.
  - `cmapTexture` is exported for custom shaders.
  - Guidance: refinement runs must end at the same time; split-operator energy is judged by its order; the plot `alpha` and `width` options are documented.
- **From the DNA tunnelling build:**
  - `stationaryStates1D` is about 10⁴× faster on fine grids (Sturm bisection plus inverse iteration; doublets stay orthogonal), so 10⁴ points take 37 ms.
  - New `lowestTridiagonal`, and a values-only mode for `eigenTridiagonal`.
  - New `stationaryDistribution` (subtraction-free, populations down to 10⁻³⁰).
  - `psim check` now:
    - names the elements that overflow on phones;
    - says when the later tests were skipped.
  - Equation values wrap on phones.
  - Guidance on the secular approximation, double-well benchmarks, logarithmic clocks, phone legend placement, and open biological questions.
- **From the enzyme-tunnelling build:**
  - `transmission1D`: stationary 1D scattering for any barrier, reproducing the exact Eckart and rectangular results.
  - `thermalAverage` and `packetTransmission`.
  - Chemistry energy units and AME2020 nuclide masses.
  - `meshRenderer.remove` frees GPU buffers.
  - The hidden lin/log redraw check now looks only at parameter controls.
  - `verification.md` warns to refine Δt with Δx.
  - `quantum.md` uses the packet-averaged transmission benchmark.
- **Docs:**
  - open quantum systems in `quantum.md`;
  - MathML needs `<mo>` for `|`;
  - citing a checked reproduction when the primary source is paywalled.

## Unreleased — runs that end

- **Clock stops at the end.** A simulation that ends sets `sim.done` (or `sim.ended`). The shell stops the clock there, shows "finished" in the HUD, and restarts the run when you press play. `advanceTo` and the frame loop stop too. Before this, the binary black-hole merger's clock kept counting after the computed inspiral ended. The same applied to the wave-packet tunnelling, black-hole inspiral, Stern–Gerlach, bubble chamber and decay chain demos; they now declare where their run ends (the only change to those pages). The pre-delivery checklist in `visual-language.md` now asks about it.

## Unreleased — README

- The README lists all demos, grouped by subject and marked when they use WebGL, with the correct count. The list is generated by `gallery.mjs`, which updates it on every add, remove and order. A test fails if the README and `docs/demos.json` disagree.
- New sections: WebGL 3D, full screen and speed; what the agent can start from (updated for the new modules). The install section is unchanged.

## Unreleased — full screen and faster GPU scenes

- **Full screen** for every page: a button on the playback bar (a corner button on static pages) and the F key. It shows the scene, HUD, legend and controls on the whole screen; Esc leaves. Where the Fullscreen API is unavailable (iPhone Safari), the scene fills the window instead. Tested in the browser fixture. The old "fit view" button looked like a fullscreen button; it now has its own icon.
- **Adaptive GL resolution.** When frames fall behind while JavaScript is cheap (the GPU is the bottleneck), the GL layer renders at a lower pixel ratio and recovers when there is headroom. `glResolution` fixes it; `view.glScale` reports it. The layer also never exceeds 8 M device pixels.
- **Volume renderer.** Smooth volumes are ray-marched at ≤ 1.5 M rays into an offscreen texture and upscaled. Hydrogen orbitals on a 2560×1440 Retina screen went from 25 to 58 fps, and hold 52 fps in fullscreen.
- **Renderer sizes.** Point and line sizes and field rectangles follow the GL canvas's real pixel ratio (`pixelRatio(gl)`).
- **From the BEC build:**
  - `psim check --advance` no longer dies after 10 s (the inner timeout ignored the outer one).
  - The build-time report is available to `bench`, `view` and `window.psim`, so precomputed scans can be shown.
  - New `fftAxis` and `phaseWinding` (vortex detection).
  - `quantum.md` covers Gross–Pitaevskii split-step stability, imaginary-time Δτ bias, vortices and a unitarity budget that scales with the run.
  - `compose.md` no longer says to write your own FFT.
- **Faster FFT** with cached twiddle tables, so `schrodingerND` is 1.75× faster (see below).
- **Gallery.** Every page is rebuilt from its unchanged source with this runtime. `gallery add` takes `--advance T` and `--click "button"` so thumbnails show the developed state.

## Unreleased — gaps reported by the famous-experiment builds

The Rutherford, Stern–Gerlach, bubble chamber and Cavendish builds were audited; these gaps are now filled and tested.

- **`psim check --params mode=…`** checks and screenshots non-default modes. Two builders had to hand-edit defaults to see theirs. The check's step-running evaluations no longer hit a fixed 10 s limit on slow Monte Carlo pages, and a timeout now names the expression.
- **Colours.** Plots, `view.color()` and the GPU renderers accept theme tokens and any CSS colour. Before this, `rgba()` tokens rendered as garbage in GPU lines and hex colours vanished from plots. `cssColor()` is exported.
- **New `core/matter.mjs`.** Bethe stopping power with the density effect, CSDA range, and PDG data for liquid hydrogen, water, silicon, lead and air. Each material's minimum ionisation is within 0.1 % of the PDG value.
- **New `core/spin.mjs`.** Spin ½ / qubit states, precession, Born probabilities and measurement. `quantum.md` now describes the semiclassical beam model.
- **Statistics (`core/data.mjs`).**
  - `curveFit` (Levenberg–Marquardt with parameter covariance; pull rms ≈ 1 on synthetic data).
  - `jackknifeRefit` (re-runs any fit with each block left out).
- **Other core additions.**
  - `maxwellSpeed` and `effusiveSpeed`, with exact samplers.
  - `rngFor(seed, i)`: an independent random stream per item.
  - `odeStepper`: `advanceTo(t, …, { exact: true })` and a mutable `hmax`; its documentation now warns that events within a single step can be missed.
- **Double slit build:**
  - **Faster solver.** `schrodingerND` is 1.75× faster (cached FFT twiddles, precomputed phase factors, and one finiteness check per step).
  - **New `probabilityCurrent`.** A spectral probability current, so what absorbing layers take can be checked against the flux through a line (agreement within 2 % in tests).
  - **`quantum.md`.** It now points to `schrodingerND` instead of telling builders to write an FFT. It also gives the right points-per-wavelength for spectral methods, says to refine Δt with Δx², calls the norm-plus-absorbed check a tautology, and says to sample detections from the flux.
- **Constants.** New: μ_B, g_e, the alpha mass, π± and π⁰ lifetimes, the K± mass and lifetime, and IUPAC `ATOMIC_WEIGHT`.
- **Shell.**
  - `getSim` is available in `controls` and on `window.psim`.
  - `view.theme`.
  - Point renderer `append()`/`clear()` for events that arrive during a run (pixel-tested).
- **Guidance.**
  - Experimenter actions as a recorded schedule.
  - `state()` must change within `testSteps`.
  - Captions wrapped with their figures.
  - Energy *balance* normalisation.
  - Meaning of `state: 'warn'` (a disclosed known limit).
  - Importance-sampling disclosure.
  - Position-space detector tracks with invisible neutrinos.

## Unreleased — gaps reported by the particle-physics builds

The asymptotic-freedom, Yang–Mills mass-gap and weak-force builds audited the skill; their gaps are now filled and tested.

- **Guidance.**
  - Open problems: a `facts.md` row for "prove/solve an open problem" (Millennium problems). It says to build evidence for a regularised version and never write "proved".
  - Non-interactive runs: a rule in SKILL.md step 5 for when the user can't be asked. Build the most representative computable version, state the scope in the Fact check line, title and model panel, and list the alternatives.
  - `particle-physics.md`:
    - a mass-gap section (smeared loops, time-slice correlators, GEVP, the bias of variational estimates, continuum evidence);
    - a strong-coupling check for any dimension;
    - measured sweep and measurement costs;
    - jackknife errors;
    - running α_s via `core/qcd.mjs`.
  - New `weak-interaction.md`: candidate processes, Michel-parameter benchmarks, why recomputing τ_μ from G_F is circular, and momentum-space views with undetectable neutrinos.
- **Lattice (`core/lattice-gauge.mjs`).**
  - `wilsonLoops` measures every R×T loop in one pass, about 16× faster than repeated `wilsonLoop` calls. It is identical to `wilsonLoop` to 10⁻¹⁶ and can be limited to chosen planes.
  - New functions: `polyakovLoop`, `plaquetteField` (local action density), `apeSmear`, `strongCouplingPlaquette(β, d)`, `siteCoords`/`siteIndex` and `qmulInto`.
  - `staple` no longer allocates, so sweeps are about 2.5× faster (6⁴: 6 ms).
- **QCD (new `core/qcd.mjs`).**
  - β coefficients to four loops in Casimir form; for SU(3) they reproduce the ζ₃ polynomial.
  - MS-bar running of α_s with flavour thresholds, which agrees with the exact 1-loop and 2-loop solutions to 10⁻¹⁰.
- **Particle data.** `PDG` in `core/constants.mjs`: masses, lifetimes, G_F, m_W, m_Z, m_H, α_s, quark masses and Michel parameters, each with its uncertainty and edition.
- **Statistics (`core/data.mjs`).**
  - New: `jackknife` and `bootstrap` (blocked, for vector observables and nonlinear estimators), the incremental `createHistogram`, `chiSquareTest`, `meanSem` and `pull`.
  - New in `core/numerics.mjs`: `gammaP` and `gammaQ`.
- **Kinematics (new `core/kinematics.mjs`).** 4-vectors and boosts, random directions and rotations, `twoBodyDecay`, RAMBO phase space (exact weights), and `acceptReject`, which refuses a wrong bound instead of biasing the sample.
- **Linear algebra (`core/linalg.mjs`).** New: `cholesky` and `eigenGeneralized` (the variational method for correlation matrices).
- **GPU (`ui/gl.mjs`).**
  - New helpers: `drawColorBar`, `geometry.arrow`, `alignZ`, `polylines` (circle, box edges, lattice grid, arrow); `transform` now accepts a rotation matrix, and `compile` is exported.
  - Renderer changes:
    - lines can carry their own colour and opacity in one draw;
    - points accept `count` (a growing cloud without re-uploading) and `soft` (the same look in light and dark themes);
    - volumes accept `smooth: false`.
- **Shell.**
  - The frame budget is now checked after every step, so a slow Monte Carlo step can't stall a frame. It is set with `frameBudget`.
  - New options: `speedRange`, `timeDigits` (integer sweep counts) and `pending: true` diagnostics ("measuring…").
  - `window.psim.setParam` now exists, as the documentation already said.
- **CLI.**
  - `psim check`:
    - names the page's own error label, not "|ΔE/E0|";
    - its parameter-edit test stays inside [min, max];
    - a value-less diagnostic names its row.
  - `verify-module`:
    - prints every check and marks malformed ones `BAD`, instead of aborting on the first;
    - `verify({ log })` prints progress.
  - `build --verify-module` reuses the last passing report while the physics files are unchanged (`--reverify` forces a new run).
  - New `psim node script.mjs` for timing scripts that use `psim:` imports.
- **Docs (`compose.md`).** The plot API; several calculations on one page; that self-tests run at t = 0 and may be async; and MathML spaces.

Also from the 4D gauge-theory builds (4D lattice projection, SU(2) as a 3-sphere):

- **Topology (from the QCD-vacuum build).**
  - `topology()`: clover field strength, giving the action density s(x), the charge density q(x), Q and S. The improved clover is within 1.5 % of Q = 1 on an exact instanton; the plain clover is 7–13 % low.
  - `bpstInstanton()`: an exact BPST instanton as a benchmark; its anti-instanton is built by time reflection and gives exactly −Q.
  - Cooling (`sweep({ cool: true })`): the action never rises and the instanton survives, becoming self-dual.
  - Ranged sweeps (`sweep({ sites })`) spread a slow sweep over frames.
- **Signed volumes.** `createVolumeRenderer` `update(…, { signed: true })` shows hue as the sign and opacity as |value| (pixel-tested).
- **`verification.md`.** Node and Chrome can differ in the last bit of exp, log and cos, so chaotic or Monte Carlo pages should verify statistics, not one seed's exact trajectory.
- **New `core/projection.mjs`. 4D rotations and 4D → 3D projections, and maps of S³: stereographic, axis–angle, and the Hopf map and fibres. All are tested.
- **Lattice.**
  - `gaugeTransform` accepts a given g (or a separate rand) and returns g.
  - New: `centreTransform` (Z₂ centre symmetry), `landauGauge` (over-relaxed, converges to θ < 10⁻¹² with every plaquette unchanged), `haarSU2`, `qexp` and `qlog`.
- **Shell.**
  - `viewControls`: a "View" section for presentation-only controls that never rebuild the model.
  - `view.setLegend()`: the legend can follow the representation.
  - `testSteps`: lowers the cost of the check's parameter-edit test for slow steps.
- **Guidance.**
  - `visual-language.md`: how to show more than three dimensions (slice, projection, manifold section, small multiples, each with a caption), periodic boundaries, gauge-dependent quantities, and additive-blend saturation.
  - `particle-physics.md`: Elitzur's theorem (links are Haar-distributed; only loops are physical), finite-temperature lattices and the Polyakov loop, and the check's warm-up.
  - `facts.md`: a row on gauge-dependent links; a rule to say when a value was quoted from memory.
  - `compose.md`: the √3 trap in `camera.fit`, and that GPU self-tests must draw before reading pixels.

## Unreleased — gaps reported by four audited builds

All four builds (lightning, sundial, see-saw, hammer and feather) read SKILL.md and six references and reused the toolkit; their audits named these gaps, now filled and tested:

- Solid objects: `createMeshRenderer`, `geometry` (box, cylinder, cone, sphere, disc, plane), `transform`/`compose`, Lambert lighting and exact planar shadows without double darkening (pixel tests: top/side shading 1.0/0.25, shadow edges at the ray-predicted x = −3 and −1, depth order).
- Shell: `presets` ("Try this", several parameters at once, every control in sync, one rebuild), `setParam`/`setParams`/`play`/`pause` for `controls` and `bench`, dropdowns now follow `setParam`, `formatTime`, `hudError: false`, `texts`, `view.speed`/`view.playing`, a readable HUD backing over GPU scenes.
- Camera: `fit(box, aspect)`, distance and pitch limits, and the wheel no longer hijacks page scrolling.
- `odeStepper` events (located on the dense output, optionally terminal): ground contact matches √(2h/g) to 10⁻¹².
- `core/sun.mjs`: the Sun's position, sunrise/sunset and equation of time (Meeus ch. 25).
- Constants: standard gravity, Earth/Moon GM and radii, gas constant, dry air and the ISA atmosphere.
- Recipe search: object names (earth, sun, electron, temperature …) are now weak mentions, so unrelated requests get "No close recipe" instead of the Solar System recipe.
- `visual-language.md`: writing for children, students and the public without lowering the standards.

## Unreleased — gaps found by auditing seven primary-school builds

- An independent audit of the builders' tool calls showed every build read SKILL.md and 6–7 references, and reused core/ (Kepler, planet data, ODE steppers, roots, quadrature, contours, radiation colours, GPU points, camera). The gaps it found:
- `ui/gl.mjs`: `createSphereRenderer` (exact ray-cast lit spheres with correct day/night line and overlap; three builds had hand-drawn their own) and `createLineRenderer` (thick 3D polylines, solid or dashed). Both draw relative to the camera target in float64, so astronomical coordinates stay precise; `camera.matrices(aspect, { relative: true })` for custom shaders. Cyclic `phase` colour map (constant lightness in OKLab).
- `compose.md` now documents every `createInstrument` option and the `ui/instrument.mjs` helpers (`el`, `button`, `toggle`, `select`, `slider`, `section`, `sci`, …): every builder had been reading the ui/ source to find them.
- SKILL.md: the page-length rule now separates physics from presentation instead of a fixed line count.
- Tests: sphere silhouette equals the exact angular size, shading equals n·l from an exact ray–sphere intersection (130 vs 129 of 255), a sphere 10⁸ units from the origin renders identically, lines land at `project()` with the requested width; relative matrices equal absolute ones; the phase map is cyclic with constant OKLab lightness.

## Unreleased — GPU-rendered scenes

- `ui/gl.mjs`: WebGL2 renderers chosen by what is drawn: glowing point sprites (10⁵–10⁶ bodies or samples, colour-mapped by a computed value), 2D fields as colour-mapped images, ray-marched 3D volumes (emission or emission–absorption), an orbit camera with `project()` for 2D labels, and perceptually uniform colour maps (viridis, inferno, magma, coolwarm).
- `createInstrument({ gl: true, drawGL })` stacks the GPU layer under the 2D overlay; PNG export composites both; pages explain when WebGL2 is unavailable.
- Tests: exact camera/projection/inverse maths and colour-map endpoints and monotonic lightness in Node; a browser fixture checks pixels (a point lands where `project()` says, field cells take the colour-map colours exactly, a centred volume is symmetric and bounded).
- References: when to use the GPU (by element count and state type) in `visual-language.md`; the API in `compose.md`.

## Unreleased — gaps found by ten real builds

Ten simulations were generated in a test project from single prompts through the installed skill. Their reports listed the following gaps; each fix is tested against exact results.

- `psim check`: pages publish their first frame at once (slow pages were reported as "frames 0"), the check waits for a live page, and it now screenshots everything below the first screen (figures, control column, phone page). SKILL.md asks for those screenshots to be inspected.
- Plots: error bars on log axes no longer break when value − error ≤ 0; reference-line labels near the top edge are no longer clipped; band and vertical-line labels that would overlap are stacked.
- Equations: `M.sup`/`M.sub` wrap both parts, so a multi-part exponent such as e⁻ˣ renders (it showed as "e − x" in the black-body derivation, found by the new scrolled screenshots).
- Particles: optional periodic `box` with the minimum image, cell-list neighbour search (`forEachPairWithin`), and a `lennardJones` force (shifted-force cutoff, Lorentz–Berthelot mixing). Gravity and Coulomb refuse a periodic box (Ewald summation is not provided).
- `metropolisSites` for discrete states (checked against exact enumeration of the 4×4 Ising model); `langevinBAOAB` reports the heat exchanged with the bath.
- `diffusion1D`: walls held at a value or given a flux (constant or time-dependent), switchable mid-run, with exact `inflow` bookkeeping.
- `schrodingerND`: optional absorbing layers; `absorbed()` reports the probability removed.
- `lyapunovExponent` (linearised equations or twin trajectory; Lorenz λ₁ ≈ 0.906 reproduced), `linearChain` (Bateman with branching) and `expm` (matrix exponential).
- `complexYlm` and `hydrogenRadial`; `contourLines` (marching squares) and `densityLevel` (probability-enclosing level).

## Unreleased — reusable matrix, field and probability maths

- Added 20 importable functions across `core/linalg.mjs`, `core/fields.mjs` and `core/probability.mjs`: real dense/sparse equations, nonlinear roots, eigenmodes, multidimensional derivatives/FFT and periodic field solves, stable probabilities, entropy/covariance and finite-spectrum equilibrium.
- Solvers reject invalid dimensions/domains and report unresolved iteration budgets. Matrix methods handle dimensional rescaling; probabilities use shifted weights and preserve tiny thermodynamic/relative-entropy signals.
- Added 28 numerical and bundling tests against independent solutions, identities, residuals, orthogonality, conservation and measured spatial/time convergence. All new modules compute identical results in Node and a self-contained bundle.
- Added a mathematical API guide and routed skill/model-selection/compose instructions to the implementations and their conditions. No new topic-specific recipes are required to use them.

## Unreleased — general building blocks by mathematical structure

- `core/ode.mjs`: general ODE stepping (RK4, adaptive Dormand–Prince with Hermite dense output) and symplectic mechanics (leapfrog, Yoshida-4); `refinementOrder` measures the convergence order.
- `core/particles.mjs`: interacting bodies built from force terms (gravity, Coulomb, uniform E/B, central body with optional 1PN, drag, springs, external potential). Chooses a symplectic method for position-only forces, otherwise an adaptive one. Reports energy balance including dissipated work, and momentum drift.
- `core/pde.mjs`: Crank–Nicolson diffusion with walls and sources, `waveND` (leapfrog in any dimension on `laplacianND`, CFL-checked, exactly conserved discrete energy), stationary states. `eigenTridiagonal` added to `core/numerics.mjs`. Merged with `core/fields.mjs`: `schrodinger1D` and `laplacian1D` removed in favour of `schrodingerND` and `laplacianND`; `wave1D`/`wave2D` replaced by `waveND`; Neumann now means the same (cell-centred zero flux) in every module.
- One catalogue for general mathematics: `references/math-toolkit.md`.
- `besselJ` and `besselJZeros` in `core/numerics.mjs` (checked against published values, an independent integral representation, the power series and McMahon asymptotics).
- Irregular domains: optional `mask`/`maskBc` for `laplacianND`, `gradientND`, `divergenceND` and `waveND`. Symmetric, so energy is still conserved to round-off. A curved boundary is a staircase (first order, measured against the exact Bessel mode of a clamped disk).
- `core/stochastic.mjs`: Gillespie, Euler–Maruyama, Langevin BAOAB, Metropolis, autocorrelation time and effective sample size.
- `core/data.mjs`: monotone or linear interpolation, histograms with Poisson errors, Kolmogorov–Smirnov test.
- Every module is tested against exact results: analytic solutions, conservation laws, stationary distributions and measured convergence orders. No topic modules were added.

## Unreleased — general model selection

- A shared workflow formulates any subject as a mathematical specification: scientific aim, state, governing law/data, conditions, observables, representation, method and evidence. Unfamiliar subjects require researched foundations, not another mandatory topic guide.
- Method selection uses mathematical structure rather than a table of topic keywords. Reuse of a solver and reuse of a layout are independent; recipes and domain notes are optional shortcuts.
- Calculation, observation and representation are explicit separate responsibilities. Legends derive from the mapping of quantities to visual channels, and tests derive from each scientific claim. Ask about material ambiguity or a reduction that changes scientific scope.
- Domain notes remain optional worked guidance. Clarified pure SU(2) versus QCD and the meaning of symbolic visuals; general rules apply the same distinctions to unfamiliar representations.
- Verification and trust tiers describe the selected model's scope. Clarifying replies do not require invented artifact links or verification reports.
- Added a shared generalisation rubric and seven cross-cutting evaluation specifications covering unfamiliar supplied laws, missing foundations, observables, representations and scope changes. These are not recorded as completed agent evaluations. No new simulation recipes or solver code were added.

## Unreleased — review fixes

- Removed `scripts/install.mjs`: install with `npx skills add janosvajda/p-simtor` (see README).

- Dormand–Prince enforces local error tolerance at any time scale and throws on step underflow or an exhausted budget.
- Schwarzschild deflection uses a scaled turning-point root and cancellation-free quadrature; relativistic orbits complete the requested radial periods with refined periapsis events, retaining coordinate-time animation support.
- SU(2) heat-bath updates sample Haar measure at zero coupling and use an efficient weak-coupling rejection sampler. Negative and non-finite coupling are rejected.
- Verification requires explicit boolean results and finite measurements, including imported reports. Gravity checks gate quantitative benchmarks, measure maximum conservation errors over the whole run, and attach solver settings. Pythagorean and Lagrange recipes use finer steps to resolve error spikes hidden by endpoint-only checks.
- Order calibration reports its step separately when production refinements reach round-off; production-step conservation and predictability measurements remain separate.
- Changing gravity's solver or accuracy withdraws stale reports and marks verification needed, including in exported data. Browser input handlers and CDP requests have deadlines; frozen renderers are closed without further interaction.
- Browser stress tests exposed fast Solar System flybys: timestep selection now includes pair crossing speeds. Non-finite or backward numerical steps restore the last valid state and pause playback with an explanation.
- Galaxy circular speeds have a finite central limit. Galaxy text distinguishes model speeds, expansion, unbound stars and distant stars; Solar System notes identify JPL fitted approximate elements.
- Blackbody's axis selector works independently of temperature; plot ranges include both cold and hot spectral peaks.
- Added 17 regression tests and browser interaction checks for stale reports and spectrum axes (144 Node tests total).

## 0.3.1 — 2026-10-13: math audit

Independent line-by-line review of every core module and example physics, with a regression test for each finding (tests: 69 → 116; every exported function in core/ is now tested).

Fixed (wrong results):
- `solveKepler` diverged for negative mean anomaly at e ≥ 0.8 (now range-reduced; guards e ∉ [0, 1)).
- Solar System time units: G = 4π² is the Gaussian year, not the Julian year used everywhere else; now G = GM☉·yr²/AU³ = 39.476926 (Earth's phase drift over 30 yr: 10′ → 2′).
- `Msun` mixed CODATA 2014 G with 2018; now GM☉/G. `hbar` = h/2π exactly.
- Proton example: orbit drawn off-centre by r·ΩΔt/2 (Boris half-step start); now centred to 1e-8 r, with a new check.
- `besselI` for negative odd arguments; `linearFit` cancellation with offset x; `realYlm` underflow at high l (new stable `normalizedLegendre`).

Fixed (edge cases): galaxy potential/force at r = 0; `dopri5` hang on NaN; `erf(±∞)`; Γ at poles; `meanWithError` on short series; `integrateToInfinity` tail cut; N-body step for coincident bodies; energy drift for E₀ = 0; `turningPoint`, `deflectionExact` and `shadowAngle` outside their domains (shadow angle at r₀ = 3M returned NaN from round-off); Wilson loops larger than the lattice; black-hole colour lookup half-texel offset; `lagrangeTriangle` with G ≠ 1.

Fixed after the first user test:
- Bundler dropped all but the first name in `export const a = …, b = …` (V0 in the galaxy page became undefined → NaN diagnostic). Now parsed properly; a test compares every bundled module's exports with Node's.
- `psim check` now fails when any diagnostic shows "—" or NaN.
- Final reply must end with the "Open it" link block; the build output marks it.
- Phone layout: legend moved bottom-right, clear of scale bars.

Robustness (found by user testing and the new stress test):
- Inputs are clamped to their allowed range (and rounded for counts); typing 100 000 000 stars no longer freezes the galaxy page. Gravity edits are limited to 100× the scenario's scales; an impossible setup keeps the last valid one and says why.
- `psim check` now types 10¹⁵, −10¹⁵ and 0 into every numeric input and fails if the page freezes or errors.
- Black-body colour underflowed to NaN below ~200 K (the CMB preset crashed the page); the colour is now computed from a rescaled spectrum, the CIE red-end ratio is enforced beyond 700 nm (cold bodies tend to deep red, not green), and out-of-gamut colours are clipped. Planckian-locus accuracy is now ≤ 0.001 above 2500 K.
- Black-body stage: cool bodies say they glow in the infrared or microwaves; the spectrum strip always contains the peak.
- Demo gallery: `docs/` for GitHub Pages holds only real demos (pages the installed skill generated in a separate project). `node scripts/gallery.mjs add <page> --prompt "…"` imports one with its prompt, date and verification count; `list` / `remove` manage them.

From the Mercury test run:
- Kepler solver rejected round-off eccentricity (e = −1e−16 for a circular orbit), a bug introduced by the audit guard; |e| < 1e−12 now counts as 0.
- SKILL.md: check the request's facts first (e.g. Mercury has no moons) and say so in the first sentence.
- `psim check --advance T`: fast-forward before screenshots (trails kept, view re-framed).
- Orbits with Einstein's correction: `precessionPerOrbit`, `precessionSI`, `relativisticOrbit` (exact, Cutler–Kennefick–Poisson 1994), `ellipticK`; checked against orbit integration, the weak-field series and Mercury's 42.98″/century.

Fact checking and orbits in time (from the second Mercury run):
- Every final reply starts with a `Fact check:` line; `references/facts.md` lists common wrong premises; planet moon counts (NASA, dated) in `core/solar-system.mjs`.
- `relativisticOrbit` integrates coordinate time too: points [φ, r, t], periapsis times and `at(t)` for animation (Newtonian period = Kepler to 2e−10).

Clarified: misleading comments (TCG vs TT, sRGB matrix provenance, LUT encoding, frame-independence check), galaxy energy diagnostic normalisation, stronger 4D strong-coupling lattice test.

## 0.3.0 — 2026-10-13

- Faster, cheaper generation: verified recipe catalog (`examples/catalog.json`), `psim recipes <request>` ranks starting points, `psim new <recipe>` copies a working example into the project. SKILL.md rewritten around this fast path plus economy rules (reuse core, no prototypes, precision with purpose, two-strikes rule, browser check once); 25 % shorter.
- Solar System: `core/solar-system.mjs` (JPL J2000 elements → state vectors), `psim gravity --preset solar-system` with inner / outer / whole-system views.
- Galaxy and dark matter: `core/galaxy.mjs` (Hernquist bulge, Miyamoto–Nagai disk, NFW halo normalised like MWPotential2014), `examples/galaxy`.
- Relativity and GPS: `core/relativity.mjs`, `examples/gps-relativity` (orbit and clocks integrated; GPS ICD and IAU constants reproduced).
- Real-world tests: Earth's J2000 position, Mars' orbital period.
- UI: small numbers in inputs use scientific notation; verification values wrap under their labels.

## 0.2.0 — 2026-10-13

- Rules: compute the strongest real model instead of leading with limitations.
- `core/numerics.mjs`: seeded RNG, FFT, RK4, Dormand–Prince 5(4), Gauss–Legendre / adaptive / semi-infinite quadrature, Brent, Γ, erf, Bessel I, Laguerre, Legendre, real spherical harmonics, tridiagonal solver, linear fit, blocking error analysis.
- `core/radiation.mjs` and `examples/blackbody`: Planck's law, Wien and Stefan–Boltzmann from first principles, CIE colour; verified against CODATA, the solar constant and the CMB.
- `core/schwarzschild.mjs` and `examples/black-hole`: GPU ray-traced null geodesics, accretion disk with exact gravitational and Doppler shifts; GPU self-test against the exact shadow.
- `core/lattice-gauge.mjs`: SU(2) lattice gauge theory (heat bath, over-relaxation, Wilson loops, Creutz ratios), checked against the exact 2D solution.
- `ui/plot.mjs`, `ui/equation.mjs`, workbench layout: scientific plots, equations split into live colour-linked terms, derivations.
- Bundler: `export { … }` lists. Checker: page self-tests (`window.psim.selfTest`).
- Built-in installer `scripts/install.mjs`; Anatomy references removed from the spec.
- Delivery: builds print a clickable file:// link and the OS open command; new `psim open`; SKILL.md requires the agent to end with how to open the page.
- Code quality: ESLint (dev-only) with strict correctness rules; codebase lint-clean; `npm run check`.
- `tests/real-world.test.mjs`: 13 validations against measurements (Mercury's perihelion, solar light bending, M87*, GPS, Earth's year, CMB, Balmer-α, …).
- Fix: default N-body step now respects the shortest pair orbital time; extreme mass ratios (Sun–Earth) previously got steps far too large.

## 0.1.0 — 2026-10-13

First working version.

- Agent Skill packaging (`SKILL.md`, Agent Skills spec-compliant frontmatter), installed with its own `scripts/install.mjs` (Node built-ins only; no dependencies).
- Gravity instrument: Newtonian N-body (≤ 64 bodies) with Leapfrog, Yoshida-4, LogH and LogH+Yoshida-4 integrators, Hermite dense output, frame-rate-independent stepping, five cited three-body presets, custom scenarios, chaos twin, in-page convergence study, CSV/JSON/PNG export.
- Build-time verification: conservation, convergence order, predictability horizon, determinism, frame independence, Kepler / figure-eight / Pythagorean benchmarks.
- Compose path for any domain: `ui/instrument-app.mjs` shell, simulation contract, `verify-module`, bundler with local-module support, Node `psim:` resolver.
- Worked example: relativistic proton in uniform E and B fields (Boris pusher) with six exact checks.
- Headless-Chrome checks over the DevTools protocol: runtime errors, liveness, overflow, screenshots, UI→physics interaction tests.
- References: method selection, verification, gravity, electromagnetism, molecular dynamics, quantum, particle physics, compose contract, visual language.
