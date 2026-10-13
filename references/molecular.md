# Molecular dynamics and statistical mechanics

Status: **compose**. No template yet. Follow `compose.md`, and use this guide for the physics and checks.

## Lennard-Jones fluid (the standard first model)

V(r) = 4ε[(σ/r)¹² − (σ/r)⁶]. Work in reduced units: σ = ε = m = k_B = 1, time unit τ = σ√(m/ε). For argon (a validated LJ fit): σ = 3.405 Å, ε/k_B = 119.8 K, m = 39.948 u, so τ ≈ 2.15 ps (Rahman 1964; standard textbook values, Allen & Tildesley).

- **Box and boundaries:** periodic boundary conditions with the minimum-image convention. Box side L ≥ 2 r_c.
- **Cutoff:** r_c = 2.5σ, with the potential shifted (V(r) − V(r_c)) so energy is continuous. State the shift: it changes pressure and energy slightly versus the full potential.
- **Integrator:** velocity Verlet with Δt = 0.002–0.005 τ. Energy must fluctuate without drift (|ΔE/E| ≲ 1e-4 over 1e4 steps at Δt = 0.005).
- **Initial state:** lattice (FCC or square in 2D) plus Maxwell–Boltzmann velocities from a **seeded** RNG, total momentum removed, rescaled to target T only during equilibration. Don't rescale during production.
- **Neighbour lists or cells** keep it O(N). N = 100–2000 is fine in a browser (2D is more legible for teaching; say which).
- **Thermostats:** NVE is the honest default. For NVT use Langevin (BAOAB) or Nosé–Hoover and label it. Berendsen does not sample the canonical ensemble; if used, say so.

## Verification for MD `verify()`

- Energy drift in NVE over a long run (threshold above), and second-order convergence of the energy error with Δt.
- Speed distribution vs Maxwell–Boltzmann at the measured temperature (2D: f(v) = (m v/kT) e^(−mv²/2kT)); report a KS statistic and the sample size.
- Momentum conserved to round-off.
- Optional benchmark: LJ equation of state or g(r) at a published state point (for example T* = 0.728, ρ* = 0.8442, near the triple point; Verlet 1967). Compare the first g(r) peak position (~1.1σ) and height qualitatively. Claim quantitative agreement only if you actually computed it with enough sampling.

## Molecules (bonded)

- Diatomic vibration: harmonic or Morse potential V = D(1 − e^(−a(r−r₀)))². Benchmark against the exact Morse levels if quantised, or the classical period if classical.
- Rotation: rigid rotor. Molecular orientation via quaternions.
- Realistic molecules (water, proteins) need a validated force field (TIP3P/SPC/E, AMBER, CHARMM, OPLS) with **all** parameters cited. A browser page can run small rigid-water boxes with a cited force field (tier `approximate`). For anything quantitative beyond that, recommend OpenMM/GROMACS and keep the page illustrative.

## Visual notes

Atoms are drawn as spheres at van der Waals scale (σ/2) or smaller (disclose it). Colour encodes kinetic energy with a perceptual sequential map, or species identity, never both without a legend. Show distributions (speed histogram with the MB curve, g(r), T(t), P(t)) beside the particles: statistical mechanics is about the distributions. Periodic images can be ghosted at the box edge to make periodicity visible.
