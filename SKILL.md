---
name: p-simtor
description: Build scientifically verified, beautifully crafted interactive physics simulations as self-contained HTML instruments - gravitational two-, three- and N-body orbits and chaos, charged particles in electric and magnetic fields, molecular dynamics, atomic orbitals and quantum wavepackets, and honestly qualified particle-physics visualisations (quarks, gluons, confinement). Use whenever someone asks to simulate, model, visualise, animate or explain a physical system or phenomenon - orbits, planets, atoms, molecules, electrons, fields, waves, particles, quantum states - even if they never say "simulation". The skill picks the governing theory and numerical method, computes the evolution from the real equations, measures accuracy (conservation laws, convergence, published benchmarks), labels each result's trust tier, and refuses to animate fabricated physics.
license: MIT
compatibility: Written for Claude Code; also follows Agent Skills conventions for Codex and OpenCode. Needs Node.js 22+ (built-in modules only). Chrome or Chromium for the browser checks.
metadata:
  version: "0.1.0"
---

# p-simtor

p-simtor turns a description of a physical system into an **interactive scientific instrument**: a single HTML file whose motion is computed from the governing equations, whose accuracy is measured and shown, and whose visuals are designed to reveal the physics.

**The mathematics determines the behaviour. The visuals reveal it.** Nothing on screen moves unless a calculation moved it.

## Locate the skill directory

All commands below run from the skill directory, `$PSIM`:

- Claude Code: `PSIM="${CLAUDE_SKILL_DIR}"`
- Otherwise use the directory containing this `SKILL.md`, typically `~/.claude/skills/p-simtor`, `.claude/skills/p-simtor`, `.agents/skills/p-simtor`, `~/.agents/skills/p-simtor` or `~/.config/opencode/skills/p-simtor`.

Check it with `node "$PSIM/scripts/psim.mjs"`, which prints help. Generated files go in the **user's project** (default `./simulations/`), never inside the skill directory.

## Non-negotiable rules

1. **Computed, not choreographed.** Every displayed position, field value or probability comes from the model's equations, or from a cited published dataset. Never animate keyframes, splines or easing curves and call them physics.
2. **Rendering never touches physics.** Drawing reads state; it never writes it. Simulated time advances on its own clock and is independent of frame rate (the provided shells do this).
3. **Parameters change the calculation.** A mass, charge, field or energy control re-solves the model. It never just rescales speed or appearance.
4. **Measure before you claim.** Report conservation errors, convergence order and benchmark agreement as numbers. Never write "accurate" without a measurement behind it.
5. **Label the trust tier** (see below) and give the reason. Put assumptions, limits and visual conventions in the page.
6. **Right mathematics for the scale.** Quantum systems are shown as probability densities and wavefunctions, never as planets orbiting a nucleus. QCD is never shown as computed quark trajectories. Do not force every scale into classical point particles.
7. **Disclose visual conventions**: exaggerated sizes, symbolic colours, non-physical camera motion, time compression.
8. **When the requested fidelity is infeasible, say so** and build the best defensible alternative at the right tier. Never fake it.

## Trust tiers

| Tier | Use when | Evidence required |
|---|---|---|
| `validated` | Numerical solution of an established model, verified | Conservation and convergence measured, plus at least one exact or published benchmark |
| `approximate` | Numerical, with modelling approximations that change results (softening, thermostats, truncated basis, empirical potentials) | Same checks, and the approximation and its effect named |
| `analytical` | Closed-form solution evaluated exactly (hydrogen orbitals, free wavepacket, Kepler orbit, Cornell potential) | Normalisation or identity checks on the evaluated formula |
| `conceptual` | Schematic of something that cannot be computed here (QCD vacuum, string breaking dynamics) | Every element tied to an established result and cited; marked illustrative |

The tier labels the **model as computed in the page**, not the real world. A validated Newtonian three-body run is still Newtonian.

## Workflow

### 1. Classify the request
Read [references/method-selection.md](references/method-selection.md). Decide which of these the request is:
- an exact analytic solution exists (evaluate it; tier `analytical`);
- ODE/PDE integration is needed (`validated` or `approximate`);
- statistical or stochastic (ensembles, thermostats, Monte Carlo);
- quantum (wavefunctions and probabilities, not trajectories);
- needs a specialised package, or is infeasible in a browser (say so and offer the best defensible alternative).

### 2. Establish the theory and the numbers
Write down the governing equations, the units, the constants (import them from `core/constants.mjs`: CODATA 2018 and IAU values, never typed from memory), the initial and boundary conditions, and the regime of validity. Read the domain reference:

| Domain | Reference | Status in this version |
|---|---|---|
| Gravity: 2-, 3-, N-body, orbits, chaos | [references/gravity.md](references/gravity.md) | **Ready**: built-in template, `psim gravity` |
| Charged particles, E/B fields, Lorentz force | [references/electromagnetism.md](references/electromagnetism.md) | **Worked example**: `examples/proton-magnetic-field/` |
| Molecular dynamics, Lennard-Jones, thermal motion | [references/molecular.md](references/molecular.md) | Compose (guidance and checks, no template yet) |
| Atomic orbitals, wavefunctions, wavepackets, tunnelling | [references/quantum.md](references/quantum.md) | Compose (guidance and checks, no template yet) |
| Quarks, gluons, confinement, QCD | [references/particle-physics.md](references/particle-physics.md) | Compose, tier `analytical`/`conceptual` only |

### 3. Choose the method
Pick it by accuracy, stability and cost, not habit. Use symplectic or time-reversible integrators for long Hamiltonian runs, adaptive or regularised stepping for close encounters, exact rotations (Boris) for magnetic gyration, and split-operator or Crank–Nicolson for the Schrödinger equation. Record why in the page's model panel.

### 4. Build the simulation
**Path A: gravity (ready).** Write a scenario JSON (format in [references/gravity.md](references/gravity.md)) or use a preset:
```bash
node "$PSIM/scripts/psim.mjs" presets
node "$PSIM/scripts/psim.mjs" gravity --preset figure-eight --out simulations/three-body.html
node "$PSIM/scripts/psim.mjs" gravity --scenario simulations/my-system.json --out simulations/my-system.html
```
This verifies every scenario (conservation, step-refinement convergence, determinism, frame-rate independence, Kepler, figure-eight and Pythagorean benchmarks) and bakes the measured report into the page. Exit code 2 means a check failed: fix it, don't ship it.

**Path B: compose (any other domain).** Read [references/compose.md](references/compose.md). In the user's project, write:
1. `physics.mjs`: the model, following the simulation contract, plus a `verify()` function that checks it against exact or published results;
2. `page.html`: calls `createInstrument()` from `psim:ui/instrument-app.mjs` with a domain-appropriate `draw()`.

Start by copying `examples/proton-magnetic-field/` (the reference implementation). Then:
```bash
node "$PSIM/scripts/psim.mjs" verify-module simulations/foo/physics.mjs
node "$PSIM/scripts/psim.mjs" build simulations/foo/page.html --verify-module simulations/foo/physics.mjs --out simulations/foo.html
```
Unit-test the physics in Node: `node --import "$PSIM/scripts/register.mjs" --test simulations/foo/*.test.mjs`.

### 5. Verify in the browser
```bash
node "$PSIM/scripts/psim.mjs" check simulations/three-body.html
```
This loads the page in headless Chrome (dark, light and phone widths). It checks for runtime errors, live integration and horizontal overflow, and runs interaction tests: for gravity it drives the real UI and confirms the page's integration is **bit-identical** to an independent Node run, before and after a mass edit. **Open and look at every screenshot** it prints, and fix anything that looks wrong: overlaps, unreadable text, empty stage, misleading scale. In sandboxed agents (Codex), Chrome may be blocked; ask to approve the command or have the user run it. Report the check as not run rather than assume it passed.

### 6. Visual craft
Read [references/visual-language.md](references/visual-language.md) before drawing anything new. Each phenomenon gets its own representation: trails for trajectories, field lines and arrows for fields, density volumes or isosurfaces for probability, energy-level diagrams for spectra. Never use one generic particle animation for everything.

### 7. Deliver
Tell the user, briefly:
- the file(s) produced and how to open them (double-click the HTML file; works offline);
- what model was computed, by which method, and the trust tier;
- the measured verification numbers (energy error, convergence order, benchmarks) and the predictability horizon for chaotic systems;
- limitations, and anything not verified (for example, browser checks not run).

## What p-simtor will not do

- Present a classical planetary electron, a hand-drawn "quark orbit", or an eased tween as computed physics.
- Claim quantitative molecular accuracy without a validated force field and its citation.
- Claim to solve QCD, nuclear structure or turbulence when it has not. Offer the defensible alternative (see `particle-physics.md`).
- Silently fall back to a fake when something is out of scope. State the limitation in the reply and in the page.

## Files

- `scripts/psim.mjs`: CLI (`presets`, `gravity`, `verify`, `build`, `verify-module`, `check`)
- `core/`: physics, browser- and Node-compatible: `nbody.mjs` (forces, invariants, Leapfrog/Yoshida/LogH integrators, Hermite dense output), `nbody-verify.mjs`, `kepler.mjs`, `presets.mjs`, `constants.mjs`
- `ui/`: `instrument.css` (visual system), `instrument.mjs` (controls, plots, clock, export), `instrument-app.mjs` (shell for composed pages), `orbit-view.mjs`
- `templates/gravity.html`: the gravity instrument
- `examples/`: reference implementations with built pages
- `references/`: method selection, verification, per-domain guides, compose contract, visual language
