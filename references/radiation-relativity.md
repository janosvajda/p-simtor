# Thermal radiation and general relativity

## Black-body radiation (`core/radiation.mjs`, example `examples/blackbody/`)

Everything follows from Planck's law, B_λ(λ,T) = (2hc²/λ⁵) · 1/(e^{hc/λkT} − 1):
- `planckLambda`, `planckNu`: spectral radiance per wavelength and per frequency. These have different peaks; say which one you show.
- `WIEN_B` (from the exact root of x = 5(1 − e⁻ˣ)), `SIGMA_SB` (from 2π⁵k⁴/15h³c²), `exitance` σT⁴, `photonDensity`, `MEAN_PHOTON_ENERGY_KT`.
- `bandRadiance(T, λ0, λ1)`: power in a band, for example the visible fraction.
- `rayleighJeans`, `wienApprox`: the classical and short-wavelength limits. Plotting Rayleigh–Jeans next to Planck shows the ultraviolet catastrophe.
- Colour: `blackbodyColor(T)` gives the exact chromaticity (CIE 1931 analytic fit, ±0.002) mapped to sRGB with normalised brightness (disclose this). `wavelengthColor(nm)` gives spectrum strips.

**Real-world checks to use in `verify()`:** solar constant 1361 W/m² from T_eff = 5772 K; solar luminosity 3.828×10²⁶ W; CMB peak 160.2 GHz and 411 photons/cm³ at 2.7255 K; Planckian locus at 6500 K (0.3135, 0.3236); CIE illuminant A.

**Natural extensions** (all computable): cavity modes and the derivation of Planck's law by counting modes; greybody emissivity; a planet's equilibrium temperature; a star's colour–magnitude placement; the cooling of a hot object (dT/dt = −σAT⁴/C, an ODE).

## Black holes (`core/schwarzschild.mjs`, example `examples/black-hole/`)

Geometric units G = c = 1, lengths in M (r_s = 2M). These are exact general-relativity results for test particles and light:

| Quantity | Value | Function |
|---|---|---|
| Event horizon | r = 2M | `R_S` |
| Photon sphere | r = 3M | `R_PHOTON` |
| Shadow edge (critical impact parameter) | b = 3√3 M | `B_CRIT`, `shadowAngle(r₀)` |
| ISCO | r = 6M, efficiency 1 − √(8/9) = 5.72 % | `R_ISCO`, `orbitalEnergy` |
| Photon orbits | u'' = −u + 3Mu² | `tracePhoton(b)` (RK4 in φ) |
| Exact deflection | quadrature | `deflectionExact(b)` |
| Disk redshift | g = √(1−3M/r) / ((1 − Ωλ)√(1−2M/r₀)) | `diskRedshift` |

**Rendering:** ray-trace on the GPU (WebGL2 fragment shader). Each pixel's photon moves in a plane, so a 1D ODE in that plane is all a pixel needs; it's exactly `tracePhoton`. Copy the shader in `examples/black-hole/page.html`. Thin disk: first equatorial crossing with 6M ≤ r ≤ r_out; emitted black-body at T(r), observed at gT with intensity ∝ g⁴T⁴. **Always add a GPU self-test**: render the shadow on a white sky, read pixels and compare the edge with `shadowAngle` (see `window.psim.selfTest` in the example; `psim check` runs it).

**Approximations to disclose:** the Shakura–Sunyaev temperature profile (rather than the full Page–Thorne profile), the opaque thin disk, the static observer, and the tone-mapped display brightness.

**Extensions:**
- *Orbits with Einstein's correction (Mercury, stars near Sgr A\*, binary pulsars as test-mass approximations):* use `relativisticOrbit(p, e, { orbits, newtonian })` for the rosette and the motion: it returns points [φ, r, t] and `at(t)` for the position at any time (times in GM/c³), plus a Newtonian twin with `newtonian: true`. Don't write your own orbit integrator for this. Use and `precessionPerOrbit(p, e)` or `precessionSI({ GM, a, e })` for the exact numbers (Cutler, Kennefick & Poisson 1994). Mercury: 42.98″ per century. The real effect is tiny: offer a clearly labelled exaggeration, and quote the real value.
- *Rotating (Kerr) black holes:* separable geodesics with Carter's constant. Doable, but needs a new solver and new checks.
- *Gravitational waves, inspiralling binaries, colliding black holes:* see [gravitational-waves.md](gravitational-waves.md) (`core/gravitational-waves.mjs`).
- *Gravitational time dilation:* dτ/dt = √(1−2M/r). GPS clocks run fast by ≈ +38.5 µs/day net.
- *Event Horizon Telescope:* the shadow diameter 2√27 GM/(c²D) gives about 40 µas for M87* (M = 6.5×10⁹ M☉, D = 16.8 Mpc). The measured ring is 42 ± 3 µas.
