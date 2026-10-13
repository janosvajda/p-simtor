# Method selection

How to go from "simulate X" to a defensible mathematical model and numerical method. Work through the questions in order and write the answers down. They become the page's model panel.

## 1. What kind of physics is this?

| Signal in the request | Nature | Representation | Typical tier |
|---|---|---|---|
| Orbits, planets, stars, gravity, "three-body" | Classical Hamiltonian ODE | Trajectories, phase space | validated |
| Charged particle, cyclotron, magnetic bottle, beam | Classical (relativistic) ODE in fields | Trajectories plus field lines | validated |
| Gas, liquid, "atoms bouncing", temperature, Lennard-Jones | Many-body ODE plus statistical mechanics | Particles plus distributions (speed histogram, g(r)) | approximate or validated within the potential |
| Orbital shapes, electron cloud, "where is the electron" | Stationary quantum states | Probability density \|ψ\|², phase colour, nodes | analytical |
| Wavepacket, tunnelling, double slit, time evolution | Time-dependent Schrödinger PDE | ψ(x,t) as amplitude plus phase, \|ψ\|², T/R coefficients | validated |
| Energy levels, spectra, transitions | Eigenvalue problem | Level diagram, spectrum | analytical or validated |
| Quarks, gluons, confinement, proton structure | Quantum field theory (non-perturbative) | Potentials, flux-tube results, parton distributions (cited) | analytical or conceptual |
| Nuclear decay, fission | Stochastic and quantum | Decay statistics (exact exponential law), not nucleon trajectories | analytical, or approximate (Monte Carlo) |
| Fluids, turbulence, weather | PDE, often infeasible at full fidelity | Reduced model (2D, low Re, lattice Boltzmann) with stated limits | approximate |

## 2. Is there an exact solution?

Use it whenever one exists, even inside a numerical page, as the benchmark:

- Kepler two-body orbit (`core/kepler.mjs`), figure-eight and Lagrange solutions of three-body problem
- Uniform-field charged-particle motion: helix, gyroradius, E×B drift
- Hydrogen-like eigenfunctions ψ_nlm and energies E_n = −Z²·13.6057 eV / n²
- Free Gaussian wavepacket spreading; square-barrier transmission T(E)
- Harmonic oscillator (classical and quantum)
- Radioactive decay N(t) = N₀ e^(−λt)

If the whole request is answered by the exact solution, evaluate it (tier `analytical`). Don't integrate numerically just to look busier.

## 3. Choose the numerical method

| Problem | Method | Why | Avoid |
|---|---|---|---|
| Few-body gravity with close encounters | LogH leapfrog with Yoshida-4 composition (`logh4`) | Time-transformed symplectic: steps shrink near encounters; exact Kepler shape | Fixed-step anything (energy error ~1e2 on Pythagorean problem, see tests) |
| Softened N-body (N ≳ 10), long runs | Yoshida-4 or KDK leapfrog, fixed step | Symplectic: bounded energy error | RK4 (secular energy drift) |
| Charged particle in B (± E) | Boris pusher (relativistic form) | Exact energy conservation in pure B; phase-space volume preserving | RK4 (spirals out or in) |
| Molecular dynamics (NVE) | Velocity Verlet | Symplectic, time-reversible, standard | Euler; changing dt mid-run |
| Thermostatted MD (NVT) | Langevin (BAOAB) or Nosé–Hoover | Samples the canonical ensemble correctly | Velocity rescaling presented as physical dynamics |
| 1D/2D time-dependent Schrödinger | Split-operator FFT (periodic), Crank–Nicolson (boxed) | Unitary: norm conserved to round-off | Explicit Euler (unstable) |
| Stationary states without closed form | Finite differences plus tridiagonal eigensolver, or shooting | Converges at known order | Guessing shapes |
| Smooth non-Hamiltonian ODEs | Adaptive RK45/DOP853 with error control | Error estimate per step | Fixed RK4 without a convergence test |

## 4. Units, scales and precision

- Pick a unit system where numbers are O(1): Hénon N-body units (G = 1), AU–M☉–yr (G ≈ 4π², `G_AU_MSUN_YR`), reduced Lennard-Jones units (σ = ε = m = 1), atomic units (ħ = mₑ = e = a₀ = 1), SI for single particles in lab fields. Show the unit system in the page and convert to physical units on display when it helps.
- Constants come from `core/constants.mjs` (CODATA 2018 / IAU). Never type them from memory.
- Doubles have ~16 digits. Relative errors below ~1e-14 are round-off. Chaotic systems amplify round-off exponentially. Report the predictability horizon, not a false long-term trajectory.

## 5. Time stepping and the display clock

- The step size is set by the physics (crossing time, gyroperiod, vibrational period, CFL condition), never by the frame rate.
- The shells advance the integrator to a target **physical** time each frame, within a compute budget. If the budget is hit, playback slows ("compute-limited") rather than taking larger steps.
- Positions between steps come from interpolation of computed states (cubic Hermite for ODEs). Interpolation is drawing, not physics.

## 6. Hard requests: compute the strongest real thing

Never open with "this can't be simulated". Find the most first-principles calculation that fits in a page:

| Request | Build this (real computation) | Then add |
|---|---|---|
| QCD, confinement | SU(2) lattice gauge Monte Carlo (`core/lattice-gauge.mjs`): Wilson loops, area law, measured string tension | Cornell fit to the measured potential; one-line limits (SU(2), quenched, finite lattice) |
| Black holes | GPU null-geodesic ray tracing (`core/schwarzschild.mjs`) | Accretion disk, black-body colours, redshift, shadow check |
| Molecules | Classical MD with a cited force field; small-molecule normal modes | Hydrogen-like orbitals for single-electron intuition |
| Many-electron atoms | Numerical radial Schrödinger (Numerov) in a screened potential | Hydrogen-like comparison |
| Nuclear structure | Woods–Saxon single-particle levels by numerical diagonalisation; shell gaps at magic numbers | Semi-empirical mass formula |
| Turbulence | 2D Navier–Stokes (spectral, FFT in `core/numerics.mjs`) at moderate Re | Energy spectrum vs k^(−5/3) / k^(−3) |
| Cosmology | Friedmann equations with Planck 2018 parameters (ODE) | Ages, distances, horizon |

State the remaining limits in one line, in the model panel and at the end of the reply.
