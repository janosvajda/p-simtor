# Quantum and atomic systems

Status: **compose**. No template yet. Follow `compose.md`. Quantum results are **probabilities and amplitudes**. Electrons do not move on trajectories. Never draw an electron circling a nucleus and call it quantum mechanics. A Bohr-model picture may appear only as a labelled historical or conceptual comparison.

## Spin, two-level systems and measurement (`core/spin.mjs`)

Spin ½ (and any qubit) is an internal state, a spinor: `ketAlong(n̂, ±1)`, `rotate(ψ, n̂, φ)` (precession, exp(−iφσ·n̂/2)), `bloch(ψ)`, `probability(ψ, n̂, s)` and `measure(ψ, n̂, rand)` (Born rule, collapse). A particle with spin moving through apparatus (Stern–Gerlach, spin rotators, neutron or atom beams) is modelled **semiclassically**: a classical centre-of-mass path per branch plus the spinor, with the branch chosen by the Born rule where the beams separate. Say so in the model panel; it is not a solution of the full spatial Schrödinger equation (no recombination or interference between branches, unless you build that). Checks: P(+) = cos²(θ/2), sequential analysers (z → x → z gives ½), and the classical-moment alternative as a labelled contrast. Constants: `muB`, `ge` in `core/constants.mjs`; thermal beams: `maxwellSpeed`, `effusiveSpeed` in `core/probability.mjs`.

## Hydrogen-like orbitals (exact; tier `analytical`)

ψ_nlm(r,θ,φ) = R_nl(r) Y_lm(θ,φ), with energy E_n = −Z² Ry / n² (Ry = 13.605693 eV, `core/constants.mjs`).

- Radial: R_nl(r) = N ρ^l e^(−ρ/2) L^(2l+1)_(n−l−1)(ρ), with ρ = 2Zr/(n a₀). Use **generalised Laguerre** L^(α)_k by its three-term recurrence. Normalise by N² = (2Z/na₀)³ (n−l−1)! / (2n (n+l)!).
- Angular: for the familiar lobes use **real** spherical harmonics (p_x, p_y, p_z, d_xy, …). State that real orbitals are superpositions of m = ±|m| states. Complex Y_lm give rings of |ψ|² with phase winding: show the phase as hue if you use them.
- Examples: 1s: R = 2(Z/a₀)^(3/2) e^(−Zr/a₀). 2p: R = (1/(2√6))(Z/a₀)^(3/2)(Zr/a₀) e^(−Zr/2a₀); p_z ∝ R₂₁ cos θ.
- Reduced-mass correction (μ = mₑm_p/(mₑ+m_p)) shifts energies by ~5e-4. Mention it if quoting energies to 4+ digits.

**`verify()` checks:** ∫|ψ|² dV = 1 by quadrature (≤1e-6); orthogonality of states you display; n−l−1 radial nodes and l angular nodes, located numerically; ⟨r⟩ = (a₀/2Z)[3n² − l(l+1)] compared with quadrature; energy formula vs Ry.

**Representation options:**
- Point cloud sampled from |ψ|² (exact sampling: inverse-CDF in r from the radial distribution r²R², then rejection in angle). Use a seeded RNG and state the sample count.
- Volumetric density by ray-marching |ψ|² on the GPU (WebGL fragment shader). Evaluate ψ analytically per sample, with no texture approximation, and use a log or gamma transfer function (disclose it).
- Isosurface enclosing 90% probability (find the level by sorting sampled densities). Colour the lobes by the **sign or phase of ψ**, with a legend.
- Radial distribution P(r) = r²|R|² plot alongside, with nodes marked.

## Time-dependent problems (tier `validated`)

iħ ∂ψ/∂t = [−ħ²/2m ∇² + V]ψ. Use atomic units, or ħ = m = 1 with stated scales.

- **Split-operator FFT** (periodic box, any dimension): use `schrodingerND` from `core/fields.mjs` (ψ ← e^(−iVΔt/2ħ) F⁻¹[e^(−iħk²Δt/2m) F[e^(−iVΔt/2ħ) ψ]], unitary to round-off, second order in Δt; ~5 ms per step at 256² in Node). `absorb: { width, strength }` adds absorbing edges when waves must leave the box (say so: norm is then not conserved by design).
- **Crank–Nicolson** (hard-wall box, 1D): (1 + iHΔt/2ħ)ψⁿ⁺¹ = (1 − iHΔt/2ħ)ψⁿ, solved with the Thomas tridiagonal algorithm (`solveTridiagonal`). Unitary, second order.
- **Grid:** the spectral split-operator is exact for every mode below Nyquist, so ~4–5 points per shortest wavelength (k_max Δx ≲ 1.3) suffice; finite differences (Crank–Nicolson) need ~10 (k_max Δx ≲ 0.6). Sharp walls and steps carry high-k content, so confirm by refinement either way.
- **Refinement:** the splitting error grows with the largest kinetic energy on the grid, ∝ (ħk_max²/2m)Δt. Halving Δx at fixed Δt quadruples k_max² and can make results *worse*; refine Δt with Δx (Δt ∝ Δx²), and report the convergence order from the asymptotic pair (e.g. Δt/2 → Δt/4 → Δt/8).
- **What leaves through absorbers** must be measured independently: integrate the flux of `probabilityCurrent(re, im, grid)` through a line in front of the layer over time. `absorbed()` is just the norm deficit, so "norm + absorbed = 1" is a tautology, not a check.
- **Nonlinear and rotating problems** (Gross–Pitaevskii, i∂ψ/∂t = [−½∇² + V + g|ψ|² − ΩL_z]ψ): `schrodingerND` handles only a static linear V. Write the split-step yourself with `fft`/`fftAxis` (single-axis transforms make the rotation term exact by alternating directions: Bao, Wang & Markowich) and `fftND`; evaluate g|ψ|² at each potential half-step. Split-step on nonlinear equations has a resonance instability at large Δt (energy jumps by orders of magnitude): find the stable Δt by halving and say so.
- **Ground states by imaginary time** (t → −iτ, renormalise each step): the fixed point carries an O(Δτ) bias in μ and E; extrapolate Δτ → 0 (two step sizes) before comparing with Thomas–Fermi or exact results. Benchmarks: μ = ħω at g = 0; 2D Thomas–Fermi μ = √(g/π) (ω = 1, ∫|ψ|² = 1).
- **Vortices**: `phaseWinding(re, im, { shape })` gives the ±1 phase windings per plaquette and their positions; a single imprinted vortex must give exactly one core and ⟨L_z⟩ = ±1.
- **Unitarity budget**: round-off grows with steps and grid size; judge norm drift against √(steps)·10⁻¹⁵·(grid points)^½ rather than a fixed 10⁻¹².
- **Detections** (single particles on a screen, Tonomura-style build-up): sample events from the computed flux into the screen (time and position), not from |ψ|² at one instant, and label the dots as draws from the computed distribution, not computed trajectories.

**`verify()` checks:** norm conservation (≤1e-12 without absorbers; with absorbers a flux budget, above); free Gaussian width σ(t) = σ₀√(1 + (ħt/2mσ₀²)²) and centre ⟨x⟩ = x₀ + ħk₀t/m vs exact; for a square barrier, the transmitted probability for a narrow-in-k packet vs the exact T(E) (agreement to the packet's k-spread; state it); energy expectation ⟨H⟩ conserved; second-order convergence.

**Representation:** plot Re ψ, Im ψ and |ψ|², or |ψ| with phase as hue (a standard domain-colouring convention, which must be disclosed). Draw V(x) on the same axes (scaled, disclosed) with the mean energy as a line. For tunnelling, show the transmitted and reflected probabilities as numbers over time.

## Superposition and energy levels

- Superposition of eigenstates evolves exactly: ψ(t) = Σ cₙ ψₙ e^(−iEₙt/ħ) (tier `analytical`). Show the beat frequency (E₂ − E₁)/h as a number.
- Level diagrams: draw energies to scale (the vertical axis is energy, with units). Transitions are drawn as arrows with wavelength λ = hc/ΔE, computed.
- Multi-electron atoms have no closed-form orbitals. Hydrogen-like shapes with an effective Z are a labelled approximation. Hartree–Fock or DFT is out of scope for a page: say so.
