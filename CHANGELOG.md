# Changelog

## Unreleased — review fixes

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
