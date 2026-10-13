# p-simtor

**An Agent Skill that turns a description of a physical system into a verified, interactive simulation.**

Ask your coding agent to *"simulate a chaotic three-body system"* and p-simtor has it choose the right theory and numerical method, integrate the real equations, measure how accurate the result is, and deliver a single self-contained HTML page: a scientific instrument you can play, pause, step, edit and export, labelled with how far its numbers can be trusted.

> The mathematics determines the behaviour. The visuals reveal it.

## Install

```bash
git clone https://github.com/<owner>/p-simtor.git
node p-simtor/scripts/install.mjs
```

This installs the skill for Claude Code, Codex and OpenCode for your user account. Options:

```bash
node p-simtor/scripts/install.mjs --agent claude-code     # one agent only
node p-simtor/scripts/install.mjs --project ~/my-project  # one project instead of all
node p-simtor/scripts/install.mjs --uninstall             # remove it again
```

Start a new agent session afterwards. Nothing else is needed: no API keys, no `npm install`, no installer packages. The installer and the skill use only Node.js built-ins.

**Requirements:** Node.js 22 or later; Chrome or Chromium for browser checks (set `CHROME_PATH` if it isn't in a standard location).

## Use

| Agent | Invocation |
|---|---|
| Claude Code | `/p-simtor simulate a gravitational three-body system` |
| OpenCode | `/p-simtor …`, or just describe the system |
| Codex | `$p-simtor simulate a proton in a 2 T magnetic field` (also listed under `/skills`) |

The skill also triggers on its own when you ask to simulate, model or visualise a physical system.

## What you get

A single `.html` file in your project (default `simulations/`) that opens offline by double-clicking, with:

- **Computed motion.** The page runs the same integrator code that was verified in Node. Changing a mass, charge, field or energy re-solves the model; it never just changes the animation speed.
- **A trust-tier badge:** *Validated numerical*, *Approximate numerical*, *Analytical*, or *Conceptual illustration*, with the reason given.
- **Live diagnostics:** energy, momentum and angular-momentum error, plotted on log axes while it runs.
- **An embedded verification report** measured at build time: conservation, step-refinement convergence order, determinism, frame-rate independence, and published benchmarks.
- **Model panel:** equations (MathML), assumptions, validity limits, visual conventions, references.
- **Export:** trajectory CSV, run JSON, PNG.
- **Dark and light themes**, phone layout, keyboard controls (Space, `.` to step, R, F).

## Capabilities, honestly

| Domain | Status | What exists |
|---|---|---|
| Gravity: two-, three- and N-body (≤ 64 bodies), orbits, chaos | **Implemented and verified** | Template, CLI, 5 cited presets, custom scenarios, chaos twin, in-page convergence study |
| Charged particles in E and B fields | **Worked example, verified** | `examples/proton-magnetic-field`: relativistic Boris pusher, 6 exact checks |
| Molecular dynamics (Lennard-Jones, thermal motion) | Guidance only | `references/molecular.md` + compose contract; no template yet |
| Atomic orbitals, wavefunctions, wavepackets, tunnelling | Guidance only | `references/quantum.md` + compose contract; no template yet |
| Quarks, gluons, confinement | Guidance only (analytical or conceptual tiers) | `references/particle-physics.md`: Cornell potential, running coupling, cited lattice results; **no QCD dynamics, by design** |

"Guidance only" means the agent gets the physics, the method choices and the required checks, and builds the model through the compose path (`physics.mjs` + `verify()` + `page.html`). Those outputs have **not** been generated and verified as part of this release.

## Measured results (v0.1.0)

Three-body problem, Hénon units, LogH + Yoshida-4 integrator, defaults:

| Scenario | Steps | \|ΔE/E₀\| | Convergence order | Converged until |
|---|---|---|---|---|
| Figure-eight (Chenciner–Montgomery), 4 periods | 13,808 | 1.4×10⁻¹⁴ | 4.01 | whole run |
| Pythagorean (Burrau), t = 0 → 80 | 92,426 | 2.4×10⁻¹¹ | 4.00 | t ≈ 59 (chaotic) |
| Hierarchical triple, t = 0 → 40 | 405,989 | 1.1×10⁻¹³ | 4.00 | whole run |
| Broucke A2, t = 0 → 20 | 16,681 | 2.3×10⁻¹² | 4.00 | whole run |

- Kepler orbit with e = 0.9 over 10 orbits vs the analytic solution: position error 9.4×10⁻⁵.
- Figure-eight returns to its starting state after one published period within 4.1×10⁻⁸.
- Pythagorean run reproduces the published outcome: the mass-3 body is ejected and masses 4 and 5 form a binary.
- A fixed-step 4th-order integrator on the same Pythagorean problem gives an energy error of about 2×10². That's why the time-transformed integrator is the default.
- In headless Chrome, the page's integration is **bit-identical** to an independent Node run, before and after a mass is edited through the real input.

Proton in a uniform field (relativistic Boris):
- energy is conserved to round-off in pure B;
- the orbit radius matches the exact discrete-Boris radius r√(1+(ΩΔt/2)²) to 3×10⁻¹¹;
- phase-error convergence order is 2.00;
- the 1 GeV gyroperiod is 135.5044 ns vs 135.5043 ns exact;
- the E×B drift matches E/B.

## How it works

```
SKILL.md                 workflow, honesty rules, trust tiers (what the agent reads)
references/              loaded on demand: method selection, verification, per-domain
                         physics, compose contract, visual language
core/                    physics, plain ES modules that run in Node and the browser
  nbody.mjs              forces, invariants, Leapfrog / Yoshida-4 / LogH / LogH-4,
                         Hermite dense output, frame-independent stepper
  nbody-verify.mjs       convergence study, determinism, frame independence, benchmarks
  kepler.mjs presets.mjs constants.mjs (CODATA 2018, IAU)
ui/                      instrument.css (visual system), instrument.mjs (controls, plots,
                         clock, export), instrument-app.mjs (shell for any domain),
                         orbit-view.mjs (orthographic 3D orbit renderer)
templates/gravity.html   the gravity instrument
scripts/psim.mjs         CLI: presets · gravity · verify · build · verify-module · check
scripts/install.mjs      installer (Node built-ins only)
scripts/lib/             bundler (inlines modules into one file), Chrome DevTools checks
examples/                three-body and proton-in-B, sources and built pages
tests/                   node:test suites for physics, bundler, and skill-spec compliance
evals/evals.json         agent-level evaluation prompts with expected behaviour
```

### Design decisions

- **Plain JavaScript, no build toolchain, no dependencies.** The same module computes the physics in Node (verification and tests) and in the page (interaction), so the verified code is the shipped code. The check proves this: page and Node results are bit-identical. Python or WASM would split the physics across two runtimes, and a bundler or npm install would make every user's environment a failure point.
- **Time-transformed symplectic integration for few-body gravity.** LogH (Mikkola–Tanikawa, Preto–Tremaine) gives exact Kepler orbit shapes and automatic step shrinking through close encounters, at leapfrog cost. Composed to 4th order with Yoshida coefficients. RK4 and fixed-step methods were rejected on evidence (see measured results).
- **Frame-rate independence by construction.** The step sequence depends only on the initial value problem. Frames only choose *when* to look, using Hermite interpolation between computed states. Tested with two very different frame schedules.
- **Single self-contained HTML deliverables.** The bundler inlines each module in its own scope and embeds the verification report. Nothing loads from a CDN, so the page works offline and can be emailed or archived.
- **Canvas 2D first, WebGL when the physics needs it.** Orbits and helices are line art and look sharper in 2D canvas. Volumetric quantum densities are where GPU ray-marching earns its cost (documented in `references/quantum.md`).
- **Compose contract over a pile of templates.** New domains need a model with `createSim`, `diagnostics` and `verify`, plus a `draw` function. The shell supplies controls, the clock, verification display, export and testing hooks, so the skill can grow without inventing fake physics to fill a catalogue.

## Verification status of this release

| Claim | Status |
|---|---|
| `scripts/install.mjs` installs, reinstalls and uninstalls for Claude Code, Codex and OpenCode locations | **Verified** (fresh empty project) |
| Full workflow from the *installed* copy (build, verify, browser check, compose path, self-tests) | **Verified** in a fresh project |
| `git clone` from GitHub | Not yet tested (repository not published) |
| A fresh Claude Code session discovers `/p-simtor` and follows the workflow end to end | Not yet tested by a separate agent session |
| Codex `$p-simtor` and OpenCode `/p-simtor` invocation | Documented per each agent's conventions; **untested** (agents not available) |
| Molecular, quantum, QCD outputs | Guidance only; not generated or verified |

## Development

```bash
npm test               # physics, bundler and SKILL.md spec tests (node:test)
npm run test:browser   # rebuild examples and run headless-Chrome checks
npm run examples       # rebuild example pages
```

## License

MIT
