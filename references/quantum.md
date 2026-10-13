# Quantum and atomic systems

Status: **compose**. No template yet. Follow `compose.md`. Quantum results are **probabilities and amplitudes**. Electrons do not move on trajectories. Never draw an electron circling a nucleus and call it quantum mechanics. A Bohr-model picture may appear only as a labelled historical or conceptual comparison.

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

- **Split-operator FFT** (periodic box): ψ ← e^(−iVΔt/2ħ) F⁻¹[e^(−iħk²Δt/2m) F[e^(−iVΔt/2ħ) ψ]]. Unitary to round-off, second order in Δt. Needs a radix-2 FFT (write it, about 40 lines, and test it against a direct DFT). Use a complex absorbing potential at the edges if waves must leave the box (and say so: norm is then not conserved by design).
- **Crank–Nicolson** (hard-wall box): (1 + iHΔt/2ħ)ψⁿ⁺¹ = (1 − iHΔt/2ħ)ψⁿ, solved with the Thomas tridiagonal algorithm. Unitary, second order.
- Grid: at least ~10 points per shortest wavelength (k_max Δx ≲ 0.6). Check convergence by halving Δx and Δt.

**`verify()` checks:** norm conservation (≤1e-12 without absorbers); free Gaussian width σ(t) = σ₀√(1 + (ħt/2mσ₀²)²) and centre ⟨x⟩ = x₀ + ħk₀t/m vs exact; for a square barrier, the transmitted probability for a narrow-in-k packet vs the exact T(E) (agreement to the packet's k-spread; state it); energy expectation ⟨H⟩ conserved; second-order convergence.

**Representation:** plot Re ψ, Im ψ and |ψ|², or |ψ| with phase as hue (a standard domain-colouring convention, which must be disclosed). Draw V(x) on the same axes (scaled, disclosed) with the mean energy as a line. For tunnelling, show the transmitted and reflected probabilities as numbers over time.

## Superposition and energy levels

- Superposition of eigenstates evolves exactly: ψ(t) = Σ cₙ ψₙ e^(−iEₙt/ħ) (tier `analytical`). Show the beat frequency (E₂ − E₁)/h as a number.
- Level diagrams: draw energies to scale (the vertical axis is energy, with units). Transitions are drawn as arrows with wavelength λ = hc/ΔE, computed.
- Multi-electron atoms have no closed-form orbitals. Hydrogen-like shapes with an effective Z are a labelled approximation. Hartree–Fock or DFT is out of scope for a page: say so.
