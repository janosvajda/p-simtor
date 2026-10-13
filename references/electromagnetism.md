# Electromagnetism, charged particles, nuclear dynamics

## Single particles in prescribed fields (worked example: `examples/proton-magnetic-field/`)

Equation of motion: du/dt = (q/m)(E + u×B/γ), with u = γv and γ = √(1 + u²/c²). Constants come from `core/constants.mjs` (`mp`, `me`, `e`, `c`, `MeV`).

**Method: relativistic Boris pusher** (half electric kick, magnetic rotation, half electric kick, drift). It conserves energy exactly in pure B and preserves phase-space volume. Its known errors: rotation angle 2·atan(ΩΔt/2) instead of ΩΔt (phase error O(Δt²)), and orbit radius r·√(1+(ΩΔt/2)²). Use ≥ 64 steps per gyration for visuals and ≥ 256 for quantitative work.

**Exact references for `verify()`:**
- r = γ m v⊥ / (|q| B), T = 2π γ m / (|q| B), pitch = v∥ T
- E×B drift v_d = E×B/B² (guiding centre; non-relativistic for E ≪ cB). Measure it from the guiding centre R = x + (γm/qB)(w × B̂), not from raw positions.
- Uniform E only: hyperbolic motion, x(t) = (mc²/qE)(√(1 + (qEt/mc)²) − 1).
- Magnetic mirror: magnetic moment μ = γ²m v⊥²/(2B) is an adiabatic invariant (approximately conserved, so check it to a stated tolerance).

Sign conventions are a classic bug. A positive charge in B = +ẑ gyrates **clockwise** seen from +z. Check that the orbit centre is where you think before comparing radii.

## Fields

- Static point charges or dipoles: Coulomb / Biot–Savart. Draw field lines by integrating dr/ds = F/|F| with RK4 and adaptive ds, seeded uniformly in solid angle around sources (seed count ∝ charge, so line density means field strength). Stop at sinks or the boundary.
- Vector-field glyphs: arrow length ∝ log|F| or normalised, with colour for magnitude. Disclose which one.
- Time-dependent fields (radiation, waves): 1D/2D FDTD (Yee scheme), Courant number ≤ 1/√d. Verify against plane-wave propagation speed and a reflection coefficient.

## Many charged particles

Direct Coulomb summation for N ≲ 1000. Softening changes physics, so label it `approximate`. Plasma: particle-in-cell is beyond a quick page. Offer a reduced model and say so.

## Nuclear interactions

Nucleon–nucleon forces are not a simple classical potential, and nuclear structure requires many-body quantum methods.
- Rutherford scattering of α particles off a nucleus (Coulomb only, outside the nuclear radius): validated classical trajectories plus the exact cross-section dσ/dΩ ∝ 1/sin⁴(θ/2). Good benchmark.
- Nuclear potential well (Woods–Saxon: V(r) = −V₀/(1 + e^((r−R)/a)), R = 1.25 A^(1/3) fm, a ≈ 0.65 fm, V₀ ≈ 50 MeV): an empirical model, tier `approximate`. Cite it.
- Decay: exponential law exactly (analytic), or Monte Carlo with a seeded RNG (approximate, show the statistical error).
- Fission and fusion dynamics: conceptual only, unless a published model is implemented and cited.
