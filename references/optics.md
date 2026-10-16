# Waves and optics

Status: **compose** with verified tools in `core/optics.mjs` (diffraction, Fresnel integrals, angular-spectrum propagation, Gaussian beams and resonators, rays and lenses, dispersion, thin films, polarisation, 2D FDTD). No template yet: follow `compose.md`. SI units throughout (metres, seconds, radians, V/m, A/m). Complex numbers are `{ re, im }`; the time factor is e^{−iωt}, so an absorbing medium has N = n + iκ with κ ≥ 0.

Light is a wave first. Pick the **lowest model that answers the question** and say which one it is: rays (geometric optics) when apertures are ≫ λ and no fringes are asked for; scalar waves (Fresnel/Fraunhofer, angular spectrum) for diffraction and interference; full Maxwell (FDTD, transfer matrices) when polarisation, near fields or sub-wavelength structure matter.

## Candidate models (pick by the question)

| Question | Model | Observables and benchmark |
|---|---|---|
| Why do slits make fringes? How sharp are grating orders? | Fraunhofer closed forms `singleSlit`, `doubleSlit`, `grating` | Zeros at a sinθ = mλ; principal maxima at d sinθ = mλ with N − 1 zeros between; N = 1, 2 limits |
| How small can a telescope or eye resolve? | Airy pattern `circularAperture` | First dark ring at 1.21967 λ/D (j₁,₁/π); 83.8 % encircled energy (Born & Wolf §8.5.2) |
| Fringes at a shadow edge, Cornu spiral | `fresnelEdge`, `fresnelCS` | I/I₀ = ¼ at the edge, peak 1.3704 at v = 1.2172; A&S Table 7.7 |
| Arbitrary aperture (hexagonal mirror, mask, hologram) far field | `fraunhofer` (FFT) | Exact DFT of the samples (Dirichlet kernel), sinc² near the axis |
| How does a beam or pattern change as it travels? Talbot carpet | `angularSpectrum` (exact scalar propagation) | Gaussian w(z) = w₀√(1 + (z/z_R)²); power conserved; Talbot revival at λ/(1 − √(1 − λ²/d²)) ≈ 2d²/λ |
| Laser beams, focusing, resonators | `gaussianBeam`, `abcd`, `propagateQ`, `cavityMode` | Waist λf/(πw₀) at the back focal plane; confocal w₀ = √(λL/2π); stability 0 < g₁g₂ < 1 |
| Lenses, aberrations, image formation | `traceRay` (exact) with `paraxialSystem` | Lensmaker thick-lens EFL; real rays converge on the paraxial focus as h → 0; spherical aberration ∝ h² and short for a positive lens |
| Reflection at glass/water, Brewster, TIR, fibres | `fresnelCoefficients`, `refract`, `criticalAngle`, `brewsterAngle` | R + T = 1; R_p = 0 at Brewster; Fresnel rhomb 45° per reflection at n = 1.51, 54.6° (Hecht §8.6) |
| Prisms, chromatic aberration, rainbows | `sellmeier` with `GLASS`, `abbeNumber`, `rainbow` | N-BK7 n_d = 1.5168, ν_d = 64.17 (Schott); rainbow 42.0° / 50.9° for n = 1.333 |
| Soap films, AR coatings, dielectric mirrors, frustrated TIR | `thinFilm` (stable scattering recursion, equivalent to transfer matrices) | Quarter-wave AR R = 0; (HL)ᴺ mirror closed form; Airy single-film formula; R + T = 1 lossless; opaque layers use decaying factors |
| Polarisers, wave plates, 3D glasses, optical activity | `jones` | Malus cos²θ; QWP at 45° → circular; HWP rotates by 2θ; S₀² = S₁² + S₂² + S₃² |
| Light in structures: waveguides, photonic crystals, scattering off objects, radiation from an antenna | `fdtd2D` (Yee) | Exact Yee dispersion `yeeDispersion`; speed c/n; Fresnel r = −1/3 at a 1→4 permittivity step; discrete energy conserved in closed boxes |

## Method

- **Closed forms first.** Every pattern with a textbook formula is drawn from the formula, not from a simulation; use simulations to show what has no formula, and check them against the formula where both exist.
- **Fraunhofer by FFT:** pad the aperture into a power-of-two window (≥ 4× its size for smooth patterns); the frequency grid is f = (k − n/2)/(n·dx) and sinθ = λf. Only |f| ≪ 1/(2dx) reproduces the continuous transform; |sinθ| > 1 is evanescent and must not be drawn as a far-field direction.
- **Angular spectrum:** exact for scalar fields in a homogeneous medium, any z. The window is periodic: pad so the field stays inside, or light wraps round. Sample at dx ≤ λ/2 when wide angles matter. For z beyond the window's far field use `fraunhofer` (or Fresnel scaling) instead.
- **Gaussian beams:** q = z + i z_R transforms as (Aq + B)/(Cq + D) with real-angle ABCD matrices (`abcd`); q is the in-medium value (z_R = πnw₀²/λ₀). Paraxial: w₀ ≳ 2λ.
- **Ray tracing:** `traceRay` is exact (sphere intersection + vector Snell); `paraxialSystem` is the first-order limit of the same surfaces. Rays report `status` ('ok', 'missed', 'clipped', 'tir'): draw only 'ok' rays as transmitted.
- **Thin films:** characteristic matrices (Born & Wolf §1.6, Macleod ch. 2); s and p separately, oblique incidence and absorbing layers allowed; unpolarised light averages R and T.
- **FDTD:** Yee grid, leapfrog in time, Courant limit Δt ≤ √(min εr)/(c√(1/Δx² + 1/Δy²)) enforced (throws). ≥ 20 cells per shortest wavelength in the densest medium for percent-level phase. Walls: 'pec', 'periodic' or first-order Mur (exact only at normal incidence with cΔt = Δx; a few % reflection at oblique incidence: keep scatterers away from walls, or say so). Sources are soft (added) unless `hard`.

## Verification

- Closed forms: zeros and maxima at the stated positions; limiting cases (N = 1 grating = single slit; bare interface in `thinFilm` = `fresnelCoefficients`).
- Conservation: R + T = 1 (lossless), Σ|U|² under `angularSpectrum`, the FDTD discrete energy in closed boxes (round-off), unitarity of wave plates.
- Published values: Airy 1.21967 λ/D, knife edge 1.3704, BK7 n_d and ν_d, rainbow 42°, Fresnel rhomb 45°.
- Measured convergence: Yee dispersion error ∝ Δx² (ratio 4 per halving); third-order spherical aberration ∝ h²; Gaussian-beam width within O((λ/w₀)²) of the paraxial formula.
- FDTD against Fresnel: pulse reflection −1/3 (E) and +1/3 (H) at a 1 → 4 permittivity step to ~1 %; speed c/n to 0.2 %; nothing left after the pulse exits through Mur walls.
- Determinism: same inputs give the same arrays; no random seeds are involved.

## Representation

- Intensity patterns: plot I/I₀ against sinθ or screen position with the zeros marked from the formula; show a log scale option for weak side lobes (disclose it). Fringes on a "screen" image use a perceptual colour map or the wavelength colour, stated.
- Wave fields: show Re E (signed, diverging map) for phase fronts, |E|² for energy; animate Re[E e^{−iωt}] from one complex solution rather than re-solving per frame.
- Rays: draw only to the last valid surface; colour by wavelength when dispersion is modelled; mark the paraxial focus and show the caustic for spherical aberration.
- Polarisation: polarisation ellipse plus the Poincaré sphere point (S₁, S₂, S₃)/S₀; state the handedness convention ([1, i]/√2 is S₃ = +1, counter-clockwise facing the source).
- FDTD: show the grid scale, the walls and their type, and the dielectric map as an outline.

**Visible scope:** "Scalar wave optics (or ray optics / 2D Maxwell, as chosen) in linear, isotropic, non-magnetic media; [paraxial / first-order Mur walls / sampling] approximations as stated." Name what is not modelled: vector diffraction near sub-wavelength apertures, coherence and bandwidth (unless averaged over λ), nonlinear and anisotropic media, and the third dimension in 2D FDTD.

## API (`core/optics.mjs`)

| Export | Signature | Result and conditions |
|---|---|---|
| `fresnelCS` | `(x)` | `{ C, S }`, A&S normalisation ∫cos(πt²/2); odd; ±∞ → ±½; ~1e-15 absolute; NaN throws |
| `singleSlit` | `(theta, { a, lambda })` | I/I₀ = sinc²(πa sinθ/λ), normal incidence |
| `doubleSlit` | `(theta, { a, d, lambda })` | sinc²β cos²(πd sinθ/λ) |
| `grating` | `(theta, { a, d, N, lambda })` | sinc²β [sin Nγ/(N sin γ)]², exact limit at principal maxima; d ≥ a, integer N ≥ 1 |
| `circularAperture` | `(theta, { D, lambda })` | Airy [2J₁(v)/v]², v = πD sinθ/λ |
| `fresnelEdge` | `(x, { z, lambda })` | Knife-edge I/I₀ at screen distance z; x > 0 illuminated side |
| `fraunhofer` | `({ re, im?, shape, dx, lambda })` | `{ re, im, intensity, normalized, fx, sinTheta }`, centred (zero frequency at n/2); shape `[n]` or `[nx, ny]`, powers of two |
| `angularSpectrum` | `({ re, im?, shape, dx, lambda, z, n = 1 })` | Propagated `{ re, im }`; z ≥ 0; periodic window; evanescent parts decay |
| `gaussianBeam` | `({ w0, lambda, n = 1 })` | `{ w0, zR, divergence, w(z), R(z), gouy(z), q(z) }`; paraxial |
| `beamFromQ` | `(q, lambda, n = 1)` | `{ w, R, z, zR, w0 }`; throws if Im q ≤ 0 |
| `propagateQ` | `(q, M)` | (Aq + B)/(Cq + D) |
| `abcd` | `.freeSpace(d)`, `.thinLens(f)`, `.curvedMirror(R)`, `.flatInterface(n1, n2)`, `.curvedInterface(n1, n2, R)` | 2×2 real-angle matrices; R > 0 centre downstream (concave mirror R > 0) |
| `abcdSystem` | `(...matrices)` | Product in the order light meets them |
| `cavityMode` | `({ L, R1, R2, lambda, n = 1 })` | `{ g1, g2, zR, w0, waistFromMirror1, w1, w2, gouy, transverseModeSpacing }`; throws unless 0 < g₁g₂ < 1 or g₁ = g₂ = 0 |
| `reflect` | `(d, normal)` | Unit reflected direction, 2D or 3D |
| `refract` | `(d, normal, n1, n2)` | Unit refracted direction, or `null` on TIR |
| `fresnelCoefficients` | `(theta1, n1, n2)` | `{ rs, rp, ts, tp (complex), Rs, Rp, Ts, Tp, theta2, tir }`; Hecht signs (r_p = −r_s at normal incidence) |
| `brewsterAngle`, `criticalAngle` | `(n1, n2)` | atan(n₂/n₁); asin(n₂/n₁), throws if n₁ ≤ n₂ |
| `traceRay` | `({ p, d }, surfaces, n0 = 1)` | `{ points, p, d, n, status }`; surfaces `[{ z, R, n, aperture? }]` ordered by z, ray towards +z |
| `paraxialSystem` | `(surfaces, n0 = 1)` | `{ matrix, efl, bfd }` from first to last vertex |
| `GLASS` | `.BK7`, `.fusedSilica` | Sellmeier `{ B, C (µm²), range (m), source }` |
| `sellmeier` | `(lambda, glass)` | n(λ), λ in metres; throws outside `range` |
| `abbeNumber` | `(glass)` | ν_d at the d, F, C lines |
| `rainbow` | `(n, k = 1)` | `{ incidence, refraction, deviation, angle }` for k internal reflections; geometric only |
| `thinFilm` | `({ lambda, n0 = 1, layers, ns, theta0 = 0, polarization = 's' })` | `{ r, t, R, T, A }` ('unpolarized': R, T, A only); layer `{ n (number or { re, im }), d }` |
| `jones` | `.linear(θ)`, `.circular(±1)`, `.polarizer(θ)`, `.waveplate(Γ, θ)`, `.halfWave(θ)`, `.quarterWave(θ)`, `.rotator(φ)`, `.system(...M)`, `.apply(M, v)`, `.intensity(v)`, `.stokes(v)` | Complex 2-vectors and 2×2 matrices; polariser θ = transmission axis, wave-plate θ = fast axis (slow axis delayed by Γ); Stokes `[S0, S1, S2, S3]` |
| `yeeDispersion` | `({ kx, ky = 0, dx, dy = dx, dt, n = 1 })` | Numerical ω (rad/s) of the Yee scheme; throws beyond the Courant limit |
| `fdtd2D` | `({ nx, ny, dx, dy = dx, dt?, courant = 0.5, mode = 'TM', epsR?, boundary = 'mur', sources = [] })` | Live fields (`Ez, Hx, Hy` or `Hz, Ex, Ey`), `step(k)`, `advanceTo(t)`, `energy()` (J/m, exact invariant for PEC/periodic without sources), `dt`, `dtMax`, `t`, `steps`; throws if Δt > Δt_max |
