# Incompressible fluids

Status: **compose** with verified tools in `core/fluids.mjs` (D2Q9 lattice Boltzmann with walls, obstacles, inflow and outflow; periodic pseudo-spectral Navier–Stokes in vorticity form; tracers; exact and published references). No template yet: follow `compose.md`. Two dimensions, low Mach number, Newtonian fluid of constant viscosity. Compressible flow, shocks, shallow water and KdV solitons are a different model: see `core/spectral.mjs`.

A fluid page has to say **which flow**, at **which Reynolds number**, in **which geometry**. Re = UL/ν decides everything visible (creeping flow, steady eddies, vortex shedding, turbulence). Give it on the page with the length and velocity it uses. A 2D simulation is a model of a 3D flow: above Re ≈ 190 the real cylinder wake turns three-dimensional (Williamson 1996), and 2D turbulence has an inverse energy cascade that 3D turbulence does not.

## Candidate models (pick by the question)

| Question | Model | Observables and benchmark |
|---|---|---|
| Flow in a pipe or channel; why is the profile parabolic? | Plane Poiseuille (body force or pressure drop between walls), `latticeBoltzmann2D` with walls | Profile against `poiseuilleProfile`; mean = (2/3) u_max; wall shear; BGK 2nd-order convergence, TRT exact |
| Flow past an obstacle; vortex street, drag | LBM with a `solid` mask, `inflow` + `pressure` outlet | Strouhal number vs `cylinderStrouhal(Re)` (Williamson 1996, 49 ≤ Re ≤ 178); C_d, C_l from `force(id)`; onset of shedding near Re ≈ 47 |
| Driven cavity, recirculation | LBM with a moving top wall | Centreline velocities vs `GHIA_RE100` (Ghia, Ghia & Shin 1982) |
| Decay of eddies, viscosity | Taylor–Green vortex in `vorticity2D` or LBM | Amplitude ∝ e^{−2νk²t} (`taylorGreen`), energy ∝ e^{−4νk²t} |
| 2D turbulence, vortex merger, shear layers (Kelvin–Helmholtz) | `vorticity2D` on a periodic box, small ν, optional forcing and drag | Energy and enstrophy conserved at ν = 0; spectra; inverse cascade; KH growth rate ≈ ΔU·k/2 for a thin layer (vortex sheet limit) |
| Forced steady flow | Kolmogorov forcing in `vorticity2D` | Laminar state ω = f/(νk² + μ), unstable above a critical Re |
| Streaklines, dye, particle paths | `bilinearSampler` + `advectParticles` on either solver | Tracers stay on streamlines of steady flows (ψ constant) |
| Couette flow, sliding plates | LBM with one moving wall | Linear profile (`couetteProfile`), exact with bounce-back |

## Method

- **Lattice Boltzmann** (`latticeBoltzmann2D`, lattice units dx = dt = 1): choose resolution N across the body and a lattice velocity u ≤ 0.1, then ν = uN/Re and τ = 3ν + ½. `latticeUnits` does this and converts back to metres and seconds; it refuses Ma ≥ 0.3 and τ < 0.51. Below τ ≈ 0.52 the grid Reynolds number u/ν must stay ≲ 10, so high Re needs more cells, not a smaller τ. TRT (default, Λ = 3/16) is more stable and puts bounce-back walls exactly mid-link; BGK is the textbook scheme.
- **Boundaries:** walls lie half a cell outside the edge cells (a channel of ny cells has height ny). Inflow is a velocity bounce-back (exact mass influx), the outlet a fixed-density `pressure` side. Impulsive starts launch sound waves that ring in forces for a few D/U; start from potential flow (with ρ from Bernoulli) or ramp the inflow by `setBoundary`.
- **Obstacles:** `solid` cells with an id 1…255; `force(id)` is the momentum-exchange force of the last step. C_d = 2F_x/(ρU²D), C_l = 2F_y/(ρU²D). Staircase circles are first-order in shape: D ≥ 10 cells for a few-percent St.
- **Vorticity solver** (`vorticity2D`): periodic only, power-of-two grid, 2/3 dealiasing, viscosity integrated exactly. Choose dt with CFL max|u|dt/dx ≲ 0.5; `step` throws above `cflMax`. Resolve the Kolmogorov scale: kmax·η ≳ 1, or add drag/hyperviscosity knowingly.
- **Performance:** LBM ≈ 40 M cell updates per second in V8, so a 200 × 100 grid runs ~2000 steps/s; call `step(n)` once per frame with n chosen for the frame budget. A 128² spectral grid takes a few ms per step.

**Flow past a body (vortex street), lessons from a real build:**
- Put the inlet about 8 D upstream of the body; at 4 D the Strouhal number came out 6–11 % high at Re = 100.
- Do not use `'outflow'` (zero gradient) on the side walls: it made the cylinder shed below the onset (Re 40–45) and the Strouhal number depended on the domain. Use `'periodic'` or free-stream inflow on the sides and report the blockage D/H.
- Near τ ≈ 0.515 the default TRT collision blew up at the inflow edge (Re 150–180) while BGK stayed stable: if TRT fails there, try BGK or more cells before lowering τ.
- Start from rest with a ramped inflow (or a potential-flow start at a lower lattice velocity); an impulsive start trips the Mach guard and sends sound waves that ring in the lift force for hundreds of D/U. Measure the shedding frequency from a wake velocity probe (lock-in or zero crossings after the transient), not from the lift force alone.
- Throughput measured in a browser page: about 22–25 M cell updates per second (the Node figure of 40 M is an upper bound); size the grid for ≥ 30 frames per second.

## Verification

- **Exact solutions:** Poiseuille (TRT to 1e-9, BGK at measured order 2), Couette, Taylor–Green decay rate (LBM within 0.5 %; spectral to round-off).
- **Conservation:** LBM mass to round-off with walls, moving lids, obstacles and body forces; momentum balance in a periodic box (steady obstacle drag = total body force); spectral energy and enstrophy at ν = 0 (1e-7 over 200 steps); measured 4th-order time convergence.
- **Published data:** Ghia et al. (1982) cavity centrelines at Re = 100 (|Δu| < 0.01 U at 48²); Williamson (1996) Strouhal number at Re = 100 (0.164 unconfined; the coarse confined test case gives 0.16–0.21: say which).
- **Guards:** τ, Mach and CFL limits throw with a message; non-finite states throw. Show max Mach (LBM) or CFL (spectral) on the page.
- **Compressibility:** LBM is weakly compressible: density varies by ~Ma² and along pressure gradients. Compare profile shapes or ρu fluxes, not raw u, when the pressure drop is large.

## Representation

- Colour vorticity with a diverging map centred on zero (sign = rotation sense); speed with a sequential map. Draw obstacles solid and walls as lines at the half-cell position.
- Streaklines (continuously released dye) differ from streamlines (instantaneous ψ contours) and pathlines in unsteady flow: label which one is drawn.
- State the frame: a vortex street seen in the cylinder frame; switch to the fluid frame only with a label.
- Show Re, the resolution (cells per D), and Ma or CFL in the model panel; for a converted physical scenario, give the physical U, L, ν.

**Visible scope:** "2D incompressible (low-Mach lattice Boltzmann / periodic spectral) Navier–Stokes at the stated Re and resolution; no turbulence model, no 3D instabilities, staircase obstacle boundaries." Say which of these limits matters for the question.

## API (`core/fluids.mjs`)

| Export | Signature | Result and conditions |
|---|---|---|
| `latticeBoltzmann2D` | `({ nx, ny, tau \| nu, collision = 'trt', magic = 3/16, solid, boundaries, force = [0, 0], density, velocity, machLimit = 0.3 })` | `{ nx, ny, tau, nu, t, rho, ux, uy, maxMach, step(n), mass(), force(id), vorticity(out?), setBoundary(side, spec) }`; lattice units; sides `left/right/bottom/top`: `'periodic'` (pairwise), `'wall'`, `{ type: 'wall', u }`, `{ type: 'inflow', u }` (u array or function of edge index), `{ type: 'pressure', rho = 1 }`, `'outflow'`; throws for τ ≤ ½, Mach > machLimit, non-finite state |
| `latticeUnits` | `({ length, velocity, viscosity } \| { Re }, cells, latticeVelocity = 0.05, minTau = 0.51)` | `{ Re, Ma, nu, tau, dx, dt, toPhysical: { length, time, velocity } }`; throws for τ < minTau or Ma ≥ 0.3 |
| `vorticity2D` | `({ omega0, shape: [nx, ny], length = [2π, 2π], nu = 0, drag = 0, forcing, dt, dealias = true, cflMax = 1 })` | `{ t, cfl, step(n), vorticity(), velocity(), streamfunction(), energy(), enstrophy() }`; periodic, powers of two ≥ 4; zero mean flow; energy/enstrophy are domain means |
| `bilinearSampler` | `({ ux, uy, nx, ny, spacing = 1, origin = [0, 0], periodic = [false, false], solid })` | `(x, y, t, out) → out = [ux, uy]`; exact for linear fields; clamped outside non-periodic axes; solid cells count as 0 |
| `advectParticles` | `(xy, velocity, dt, { steps = 1, t = 0 })` | RK4 in place on flat `[x0, y0, …]`; 4th order; throws on non-finite velocity |
| `poiseuilleProfile` | `(y, { height, force, nu })` | force·y(h − y)/(2ν), 0 ≤ y ≤ h |
| `couetteProfile` | `(y, { height, wallSpeed })` | U·y/h |
| `taylorGreen` | `(x, y, t, { nu, k = 1, U = 1, rho = 1 })` | `{ ux, uy, omega, p }`, exact NS solution, period 2π/k |
| `cylinderStrouhal` | `(Re)` | Williamson (1996) fit −3.3265/Re + 0.1816 + 1.6e-4 Re; 49 ≤ Re ≤ 178 |
| `GHIA_RE100` | `{ y, u, x, v }` | Ghia, Ghia & Shin (1982) Re = 100 centreline velocities (normalised by L and lid speed) |
