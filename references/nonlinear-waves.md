# Nonlinear waves, spectral methods and conservation laws

Status: generic tools with exact benchmarks: `core/spectral.mjs` (Chebyshev and Fourier spectral methods; integrators for NLS, KdV, Kuramoto–Sivashinsky, viscous Burgers and sine-Gordon) and `core/conservation.mjs` (Godunov-type finite volumes for advection, Burgers, compressible Euler and shallow water, with exact Riemann solvers). There is no ready-made recipe. Read [method selection](method-selection.md) and [verification](verification.md) first. The tools in [math toolkit](math-toolkit.md) cover linear PDEs. This guide does not cover incompressible flow (lattice Boltzmann, vorticity Navier–Stokes).

Choose the method by the shape of the solution. **Smooth waves** (solitons, breathers, dispersive or chaotic patterns on a periodic domain) need spectral methods: exponential accuracy and invariants close to round-off. **Shocks, bores and contact discontinuities** need finite volumes: a spectral method rings (Gibbs) at a jump. A finite-volume scheme converges only at first order there, but it conserves exactly and keeps ρ, p and h positive.

## Candidate models (pick by the question)

| Question | Model | Observables and benchmark |
|---|---|---|
| What is a soliton? Why do solitons pass through each other? | KdV u_t + 6uu_x + u_xxx = 0, two-soliton collision | Exact Hirota solution (`kdvTwoSoliton`); after the collision the faster soliton is shifted forward by Δ₁ = (2/k₁)ln((k₁+k₂)/(k₁−k₂)) and the slower one back by Δ₂ (`kdvPhaseShifts`). Invariants ∫u, ∫u², ∫(½u_x² − u³) |
| Optical pulses in a fibre, Bose–Einstein condensates, deep-water wave groups | Focusing NLS iψ_t + ½ψ_xx + κ\|ψ\|²ψ = 0 | Bright soliton (`nlsSoliton`): N = 2A, P = 2Av, H = Av² − A³/3. Norm conserved to round-off. Rogue wave: the Peregrine breather (`nlsPeregrine`) reaches peak \|ψ\| = 3 on a unit background |
| Topological defects (domain walls, dislocations, Josephson fluxons) | Sine-Gordon u_tt = u_xx − sin u | Kink (`sineGordonKink`, energy 8γ, Lorentz contraction) and breather (`sineGordonBreather`, energy 16√(1−ω²)). The winding number is fixed |
| Spatiotemporal chaos (flame fronts, thin films) | Kuramoto–Sivashinsky u_t + uu_x + u_xx + u_xxxx = 0 | Linear growth rate k² − k⁴ (maximum ¼ at k = 1/√2). Energy budget dE/dt = ∫(u_x² − u_xx²). Kassam & Trefethen (2005) run with L = 32π. Compare statistics, not trajectories |
| How does a shock form? What sets its width? | Burgers u_t + uu_x = νu_xx; inviscid limit | Exact Cole–Hopf solution (`burgersColeHopf`). Front width ~ν/Δu. Inviscid shock speed (u_L + u_R)/2 (`burgersExactRiemann`) |
| Shock tubes, blast waves, supersonic flow | Compressible Euler equations, ideal gas | Exact Riemann solution (`eulerExactRiemann`, Toro ch. 4): Sod, Lax, 123, strong-shock problems. 2D: Lax–Liu quadrants, blast in a box |
| Dam break, tsunami run-up, tidal bores | Shallow-water equations with bathymetry | Ritter (dry bed) and Stoker (wet bed) dam break (`shallowWaterExactRiemann`). Lake at rest preserved exactly. Linear wave speed √(gh) |
| Eigenmodes and boundary-value problems on an interval (stability, quantum wells, beams) | Chebyshev collocation | Spectral convergence. −u″ → (kπ/2)². Harmonic oscillator → 2k + 1. Verify by the coefficient tail or by increasing n |

## Method

- **Fourier pseudo-spectral (periodic):** n must be a power of two. `fourierDerivative` multiplies by (ik)^m. The stiff linear part (dispersion, hyperdiffusion) is integrated exactly. Use ETDRK4 (`etdrk4Periodic` and the `kdv`, `kuramotoSivashinsky` and `viscousBurgers` wrappers; Kassam & Trefethen, SIAM J. Sci. Comput. 26, 1214 (2005)), which is fourth order with no linear stability limit. For Hamiltonian waves use split-step Fourier (`nls`, `sineGordon`): it is symplectic or unitary, and `order: 4` gives Yoshida's composition. A new semilinear periodic PDE needs only `linear(k)`, `nonlinear(...)` and `multiplier(k)` in `etdrk4Periodic`.
- **Explicit limits:** ETDRK4 enforces its configured advective Courant limit. NLS split steps preserve the norm, but this alone does not guarantee accuracy or stability of perturbations: numerical sideband instabilities are possible. Sine-Gordon integrates the free wave exactly but treats sin(u) by explicit kicks; even its spatially uniform small-amplitude mode requires dt < 2 with Strang splitting. Refine dt and the spatial grid together and check the actual observables.
- **Aliasing:** `dealias: true` projects the initial field and each nonlinear term onto \|mode\| < n/3. This is the 2/3 rule for quadratic nonlinearities; higher-degree products need stricter truncation or padding. The projected initial condition can differ from the supplied samples, so disclose that truncation. Without dealiasing, invariant errors also depend on the spectral tail.
- **Nyquist mode:** odd derivatives of real fields set the Nyquist mode to zero. Energy terms involving second derivatives retain its contribution. Complex NLS kinetic energy and momentum use all FFT modes, with the negative Nyquist frequency convention of `fftFreq`. A mode at the grid limit is a numerical edge case, not evidence that the physical solution is spatially resolved.

For split-step stability conditions see [Faou, Gauckler & Lubich (2014)](https://arxiv.org/abs/1306.0656) and the [soliton instability analysis by Lakoba](https://arxiv.org/abs/1008.4974). Conservation of norm alone is insufficient evidence for a resolved physical solution.
- **Non-periodic domains:** for the sine-Gordon kink, `winding: Q` evolves u − 2πQx/L, which is periodic. For a localised packet, make the domain large enough that the tails are below the tolerance at the edges.
- **Chebyshev:** `chebyshevBVP` (Dirichlet or Robin conditions) and `chebyshevEigen` (Dirichlet) build dense collocation matrices with O(n³) cost; keep n ≲ 200. Only about the lowest 2n/π eigenvalues are resolved. The others are spurious, so confirm each reported value at a larger n.
- **Source terms with finite volumes** (bed friction, gravity on a slope, heating): `finiteVolume` advances the conservative flux part only. Add a source by operator splitting: after each `step(dt)` update `fv.q` in place with the exact or implicit source solution over the same dt (Strang: half source, step, half source for second order), and check the balance it changes (momentum lost to friction = ∫ source dt).
- **KdV accuracy vs stability:** the Courant limit (0.85) only keeps `kdv` stable. In one gallery build, a 10⁻⁶ match with the exact two-soliton solution needed a step about 0.2 of that limit (the default gave 1.8×10⁻⁵): choose dt by halving until the observable stops changing, and report the measured error.
- **Finite volumes:** `finiteVolume({ system, q0, shape, spacing, bc, limiter, riemann, integrator, cfl })` uses MUSCL reconstruction of primitive variables with limiters (minmod, MC, van Leer, superbee), a Riemann solver at every face (HLLC, HLL, Roe, Rusanov, exact), and SSP-RK2 or SSP-RK3. In 2D the sweeps are dimension-by-dimension but unsplit. Reach for **HLLC + MC** first. Use **HLL + minmod** for very strong shocks or near-vacuum. Roe needs an entropy fix and fails near vacuum (it throws there). 'exact' is the reference flux but is slower.
- **Well balancing:** shallow water with `bathymetry` uses hydrostatic reconstruction (Audusse et al. 2004), so a lake at rest stays at rest to round-off, including dry islands. A naive source term would generate spurious currents of order g·Δb.
- **Boundaries:** the choices are periodic, transmissive (zero gradient, approximately non-reflecting for outgoing waves), reflective (a wall) and `{ state }` (fixed inflow). Transmissive boundaries reflect a little when the flow is subsonic. Keep the waves away from them, or report the effect.

## Verification

- Exact solutions at the same time: soliton, two-soliton, breather, Cole–Hopf, exact Riemann solutions. Report the max or L1 error, not "looks right".
- Measured order: halve dt (ETDRK4 ≈ 4, Strang ≈ 2, Yoshida ≈ 4) and dx (smooth finite-volume flow ≈ 2; shocks ≈ 1 in L1). Report the measured number.
- Invariants: the NLS norm is conserved to round-off. The finite-volume `totals() − totals₀ − inflow` is zero to round-off for every conserved variable without a source. H, momentum and the KdV invariants drift at O(dt^order), so check the order, not only the size. The KS energy budget closes as dt⁴ when dealiased.
- Positivity: ρ > 0, p > 0, h ≥ 0 are checked after every stage; the solver throws rather than continuing with negative states.
- Symmetry: a symmetric 2D initial state must remain symmetric (x↔y transpose, mirror) to about 1e−13.
- Published numbers: Toro (2009) Table 4.3 star states (p* = 0.30313, u* = 0.92745 for Sod). Hirota phase shifts. Ritter's parabola h = (2√(gh₀) − x/t)²/(9g).
- Chaos: beyond the Lyapunov time, compare statistics (energy, spectrum, amplitude range), never pointwise trajectories. A run at a different dt is a different trajectory.

## Representation

- Draw the wave profile with the exact solution as a thin overlay curve and show the live error number. Space-time plots (x horizontal, t upward, colour = u) show soliton collisions, phase shifts and KS chaos at a glance.
- Solitons: mark each peak and its free trajectory x₀ + ct, so the phase shift is visible as a jump.
- Shock tubes: show ρ, u and p in three stacked panels with the exact Riemann solution, and label the rarefaction, contact and shock.
- Shallow water: draw the bed as a filled solid and the water surface η = h + b above it. Colour dry cells differently. Mark the √(gh) wave fronts.
- 2D Euler: show density as a heat map, plus Schlieren (\|∇ρ\|) for shocks. Show the conserved totals and their drift as numbers.

**Visible scope:** for example, "1D dimensionless KdV on a periodic domain; dissipation and higher-order dispersion neglected"; "Inviscid ideal-gas Euler, second-order finite volumes: shocks smeared over 2–3 cells; no viscosity, heat conduction or turbulence"; "Shallow-water (hydrostatic, depth-averaged) equations: no dispersion, so no undular bores; no bottom friction."

## API

### `core/spectral.mjs`

| Export | Result | Conditions |
|---|---|---|
| `chebyshevNodes(n, [a, b])` | n + 1 Chebyshev–Lobatto points, ascending | n ≥ 1, a < b |
| `chebyshevDiffMatrix(n, [a, b])` | `{ x, D, D2 }`: exact derivatives of the degree-n interpolant | Entries grow like n⁴; n ≲ 200 |
| `chebyshevInterpolant(values, [a, b])` | Barycentric interpolant `x => p(x)` | values at the ascending nodes |
| `chebyshevCoefficients(values)` | a_k with p = Σa_k T_k | Resolution check: \|a_n\| ≈ ε·max\|a_k\| |
| `clenshawCurtisWeights(n, [a, b])`, `clenshawCurtis(f, a, b, n = 64)` | Nodes and positive weights; the integral | Exact for degree ≤ n; geometric convergence for analytic f |
| `chebyshevBVP({ n, interval, a, b, c, f, left, right })` | Solves au″ + bu′ + cu = f: `{ x, u, interpolant, tail }` | left and right: `{ value }` or `{ alpha, beta, value }`; throws if singular |
| `chebyshevEigen({ n, interval, a, b, c, vectors })` | Eigenvalues of au″ + bu′ + cu = λu with u = 0 at both ends: `{ x, values: [{re, im}], vectors: [{ value, u }] }` | About 2n/π modes resolved; confirm at a larger n |
| `eigenNonsymmetric(A)` | Eigenvalues `[{ re, im }]` of a dense real matrix (balance, Hessenberg, Francis QR) | O(n³); throws if the iteration does not converge |
| `fourierDerivative(u, length, order = 1, out?)` | d^m u/dx^m of periodic samples | n a power of two; Nyquist mode zeroed for odd m |
| `etdrk4Periodic({ re0, im0?, dx, dt, linear, nonlinear, multiplier?, speed?, maxCourant = 0.85, dealias? })` | ETDRK4 stepper: `{ t, re, im, step(), advanceTo(t), spectrum(), courant }` | u_t = Lu + N(u) with diagonal L(k); throws on a Courant violation or a non-finite state |
| `kdv({ u0, dx, dt, dealias? })` | ETDRK4 stepper with `u`, `mass()`, `momentum()`, `energy()` | 6·max\|u\|·dt/dx ≤ 0.85 |
| `kuramotoSivashinsky({ u0, dx, dt, dealias?, maxCourant = 2 })` | Stepper with `u`, `mean()`, `energy()`, `dissipation()` | Chaotic: use statistics |
| `viscousBurgers({ u0, nu, dx, dt, dealias? })` | Stepper with `u`, `mean()`, `energy()` | Resolve the front: dx ≲ ν/max\|u\| |
| `nls({ re0, im0?, dx, dt, kappa = 1, V?, order = 2 })` | Split-step stepper with `norm()`, `hamiltonian()`, `momentum()` | Unitary; order 2 or 4 |
| `sineGordon({ u0, v0?, dx, dt, winding = 0, order = 2 })` | Symplectic stepper with `u`, `v`, `energy()` | u0 includes the 2πQ ramp |
| `nlsSoliton(x, t, { amplitude, velocity, x0, phase })`, `nlsPeregrine(x, t)` | ψ as `{ re, im }` | κ = 1 |
| `kdvSoliton(x, t, { k, x0 })`, `kdvTwoSoliton(x, t, { k1, k2, x1, x2 })`, `kdvPhaseShifts(k1, k2)` | Exact KdV solutions; `{ fast, slow }` shifts | k₁ > k₂ > 0 for the shifts |
| `burgersColeHopf({ primitive, nu, umax }, x, t)` | Exact u(x, t) by quadrature | primitive F(y) = ∫₀^y u₀ |
| `sineGordonKink(x, t, { velocity, x0, charge })`, `sineGordonBreather(x, t, { omega, x0 })` | `{ u, ut }` | \|v\| < 1, 0 < ω < 1 |

### `core/conservation.mjs`

| Export | Result | Conditions |
|---|---|---|
| `finiteVolume({ system, q0, shape, spacing, bc, limiter = 'mc', riemann, integrator = 'ssprk3', cfl = 0.4 })` | `{ t, q, step(dt?), advanceTo(t), maxDt(), totals(), inflow, primitive(v), courant }` | 1D or 2D; n ≥ 2 per axis; throws if Courant > 1, on a non-admissible state, or on bad options |
| `linearAdvection(a)` | Scalar system (a or [ax, ay]) | — |
| `inviscidBurgers({ dims })` | Scalar system; 'exact' is the Godunov entropy flux | — |
| `euler({ gamma = 1.4, dims })` | Conserved [ρ, ρu, (ρv), E], primitive [ρ, u, (v), p]; solvers hllc, hll, roe, rusanov, exact | ρ, p > 0 |
| `shallowWater({ g = 9.81, dims, bathymetry?, dry = 1e−10 })` | Conserved [h, hu, (hv)]; well-balanced with bathymetry; solvers hllc, hll, rusanov, exact | h ≥ 0; cfl ≤ 0.5 for wet/dry fronts |
| `burgersExactRiemann(uL, uR)` | `{ shock, speed, sample(ξ) }` | Entropy solution |
| `eulerExactRiemann({ rho, u, p }, { rho, u, p }, gamma)` | `{ pStar, uStar, rhoStarL, rhoStarR, sample(ξ) → { rho, u, p } }` | Throws when the data generate a vacuum |
| `shallowWaterExactRiemann({ h, u }, { h, u }, g)` | `{ hStar, uStar, sample(ξ) → { h, u } }`; star values are null when a dry region opens | Flat bed; wet or dry |

Tests (`tests/spectral.test.mjs`, `tests/conservation.test.mjs`) check each export against exact solutions, published values (Toro Table 4.3, the Kassam–Trefethen setup, Hirota phase shifts, Ritter), measured convergence orders, invariants to round-off, positivity and symmetry.

**Limitations:** periodic Fourier grids only (power-of-two n). No Orr–Sommerfeld or generalised (A − λB) eigenproblems. Uniform Cartesian finite-volume grids without adaptive refinement. No viscous terms in the Euler equations. Roe is not implemented for shallow water. The 2D sweeps are dimension-by-dimension, so oblique shocks are slightly more diffused than with a genuinely multidimensional solver.
