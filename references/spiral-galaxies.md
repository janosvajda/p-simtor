# Spiral galaxies: how spiral arms form

Status: verified tools in `core/disk-nbody.mjs` (self-gravitating particle-mesh disk in the rigid halo and bulge of `core/galaxy.mjs`, tested in `tests/disk-nbody.test.mjs`). There is no ready-made recipe. For a disk of test stars in a fixed potential (rotation curves, dark matter) use `core/galaxy.mjs` and the `galaxy` recipe.

Spiral arms are **not material arms** made of the same stars. Material arms would wind up within a few rotations, because inner stars orbit faster (the winding problem), yet most disk galaxies show spirals. Arms are **density waves**: regions of higher density that stars move through. In an N-body disk they form **by themselves** from the stars' mutual gravity. Shear turns any small leading disturbance into a trailing one. While it swings round, self-gravity and the epicyclic motion of the stars amplify it about tenfold (**swing amplification**: Julian & Toomre 1966; Goldreich & Lynden-Bell 1965; Toomre 1981). Particle noise is enough to start the process. The arms heat the disk, the heating raises Q, and they fade. New arms then grow somewhere else, so the spirals are **transient and recurrent** (Sellwood & Carlberg 1984, 2014). Say this in the Fact check line when the request says "stars form the arms" or treats the arms as rigid objects. Long-lived, quasi-stationary density waves (Lin & Shu 1964) and arms driven by a bar or a companion galaxy (M51) are other explanations, and this model does not include them.

The essential ingredients are **self-gravity of the disk** (stars attract each other) and **differential rotation**. Without self-gravity, as in `galaxy.mjs` stars that feel only the fixed potential, no arms form. Use this as the control experiment.

## Candidate models (pick by the question)

| Question | Model | Observables and benchmark |
|---|---|---|
| How do spiral arms form? | `createSelfGravDisk({ Q: 1.2 })`: cool disk (Q ≈ 1.0–1.3), the default MW-like halo, `diskMass` 3e10 M☉ | A_m(t) grows from the noise level ~1/√N_ring by ≳ 8× within about one rotation; the arms are trailing; the arm number is near the swing estimate m ≈ κ²R/(2πGΣX), X ≈ 2 |
| Why do some galaxies have two arms and others many? | The same model with a heavier or lighter disk relative to the halo (`diskMass`, or the halo normalisation `P`) | Heavier disk → fewer, broader arms (m ≈ 2, a bar becomes likely); lighter disk → more, thinner arms. Compare the dominant m with `swingArmNumber` at 2 R_d |
| Why is the disk not full of arms? (stability) | The same model with Q = 2.5–3 (a hot disk) | It stays axisymmetric: A_m stays within ~2× the noise. Toomre Q from `toomreQ` and `toomreQProfile` |
| Material arms vs density waves | Two runs: `selfGravity: false` (test stars, a pattern painted on at t = 0 winds up) and `selfGravity: true` | A painted pattern winds up with pitch angle tan i ≈ 1/(R t |dΩ/dR|) → 0; self-gravitating arms keep a pitch of 15–40° |
| Why do arms fade and recur? | Long cool run; plot A_m(t) and Q(R, t) | Arms grow and decay over a few rotations while Q(R) rises. The disk heats and the arms weaken unless cooling (gas, star formation) is added, which this model does not include |
| Rotation curve of the disk + halo | `rotationCurve` of the particles, `exponentialDiskCircularSpeed` (Freeman 1970) for the disk alone, `galaxy.mjs` halo | Eilers et al. (2019) data via `eilers2019` for a Milky Way comparison (the default disk is lighter than the MW's, so this is only approximate) |

## Method

- **Gravity** (`createPMSolver`): particle-mesh. Cloud-in-cell mass assignment, convolution with a Plummer-softened Green's function by FFT on a zero-padded 2N×2N mesh (Hockney & Eastwood 1988: isolated boundaries, no periodic images), and CIC force interpolation. The force kernel is convolved directly, so pair forces are exactly antisymmetric: momentum is conserved to round-off and there is no self-force. The cost is O(N_particles + N_mesh² log N_mesh) per step, and no memory is allocated per step.
- **Softening as thickness:** the disk is razor-thin (2D). Plummer softening ε stands in for the vertical thickness: in its plane, a Plummer-softened thin disk acts like the same disk seen from height ε. Use ε ≈ the mesh cell to 0.1 R_d. Larger ε suppresses the short-wavelength (many-armed) response, so a 64² mesh grows weaker spirals than a 128² mesh.
- **External potential:** the NFW halo and the optional Hernquist bulge of `galaxy.mjs` act as rigid analytic forces. Its Miyamoto–Nagai disk is replaced by the live disk. A rigid halo cannot absorb angular momentum, which suppresses bar growth compared with a live halo (Athanassoula 2002).
- **Initial conditions:** Σ ∝ e^{−R/R_d}. v_c(R) is the mesh force of the smooth disk plus the external force. σ_R = Q·3.36GΣ/κ, σ_φ = σ_R κ/2Ω, and the mean rotation follows from the asymmetric-drift equation. In a hot disk with a light centre, σ_R approaches v_c near the centre: there the epicyclic approximation fails, a few stars leave the mesh (`outside`), and the disk relaxes a little. Measure Q(R) rather than assume it.
- **Integration:** kick-drift-kick leapfrog with a fixed Δt (`suggestedDt` = inner orbital period / 50). `step` and `advance` do the same arithmetic whatever the frame schedule.
- **Time:** 1 time unit = 977.8 Myr (`TIME_UNIT_MYR` in `galaxy.mjs`). One rotation at 2 R_d ≈ 0.18 units ≈ 180 Myr in the default model.

## Verification

- PM force of a point mass vs the exact softened force: rms 3–5 % at 2–4 cells, < 7 % beyond 4 cells, and < 1.5 % beyond 10 cells, including across the whole box from a mass near one edge (a periodic solver would be wrong by O(1) there).
- Kuzmin disk on the mesh vs the exact formula with a → a + ε: < 3 % beyond 2 cells. Mesh rotation curve of the exponential disk vs Freeman (1970): < 2 % only with small softening ε and no truncation (the test's setup); with the default softening (≈ 0.1 R_d) and truncation at 5 R_d a demo measured 4.6 % — compare with a softened, truncated reference or report the difference.
- Newton's third law: Σ m a / Σ m|a| < 1e-12. Total momentum of an isolated disk is conserved to 1e-15 and the centre of mass stays put.
- In the rigid axisymmetric potential, L_z is conserved to 1e-12 over 10 rotations. Energy is bounded, and the error falls ≥ 3× when Δt is halved (O(Δt²)).
- Self-gravitating runs: ΔL_z/L_z ≈ 4e-6, and ΔE/E ≈ 1e-5 (cool) to 1e-4 (hot, a few stars leave the mesh) over 2 rotations.
- Hot disk (Q = 2.5): max A_2 = 0.047 against a noise level of 0.030. Cool disk (Q = 1.2): A_3 grows to 0.26 (8.5× noise) in ≈ 200 Myr. The swing estimate gives m_X = 2.7, the measured arms are trailing, and the pitch angle is ≈ 31°.
- The analysis tools recover an imposed 3-armed logarithmic spiral (amplitude, arm number, 25° pitch, phase winding).
- Reproducible with a seed; identical results for any steps-per-frame split.
- Things to check in a new page: A_m(t) at t = 0 is at the noise level (`fourierModes(...).noise`); `outside` stays ≈ 0; the measured Q(R) at t = 0 matches the requested Q; ΔE/E and ΔL/L are shown.

## Representation

- **Face-on view** of the particles, with the halo drawn as a faint circle or not at all, since it is rigid. Draw about 30–60k points. `createPointRenderer` from the UI kit handles more than a 2D canvas.
- **Colour by local surface density** (`localDensity`). Arms then show up as density enhancements, which is what they are. Do **not** colour by stellar age or draw blue arms and red inter-arm regions: this model forms no stars and has no gas. In real galaxies the blue arms come from young stars born in compressed gas. Say so in the legend if you mention it.
- **Tracer stars:** highlight a handful of stars and draw their trails. They pass **through** the arms, which shows that the arms are waves, not material.
- **A_m(t) plot** for m = 1–6 on a log axis, with the noise level as a dashed line: growth, saturation and recurrence. Mark the time in Myr and in rotations at 2 R_d.
- **Controls:** Q (stability), the disk-to-halo mass ratio (arm number), self-gravity on/off (control experiment), and N (noise). Show the swing estimate m_X next to the measured dominant m.
- A polar (R, φ) density map (`surfaceDensityPolar`) shows the arms as straight inclined stripes in ln R: the slope gives the pitch angle (`spiralPitch`).

**Visible scope:** "Collisionless N-body disk with a rigid halo; no gas, no star formation; arms are transient recurrent spirals in this model." Add if relevant: "2D disk, thickness approximated by softening; no bar-driven or tidally driven spirals."

## API (`core/disk-nbody.mjs`)

Units: kpc, km/s, M☉, G = 4.30091727e-6 kpc (km/s)²/M☉. A particle set is `{ n, x, y, vx, vy, m, ax?, ay? }` (Float64Arrays).

| Export | Signature → result | Conditions |
|---|---|---|
| `createSelfGravDisk` | `({ n = 30000, Rd = 2.6, Rmax = 5Rd, diskMass = 3e10, Q = 1.5, halo = true, bulge = false, P = MW, grid = 128, box = 2.6Rmax, softening = max(h, 0.1Rd), selfGravity = true, seed = 1 })` → state `{ n, x, y, vx, vy, ax, ay, m, t, steps, step(dt), advance(dt, k), energy() → { kinetic, self, external, total }, momentum(), angularMomentum(), centreOfMass(), suggestedDt, outside, table: { R, vc, kappa, Sigma, sigmaR, sigmaPhi, vphi }, pm }` | grid a power of two; box > 2.05 Rmax; throws on a non-finite position |
| `createPMSolver` | `({ grid, box, softening, Gconst = G })` → `{ accelerations(x, y, m, n, ax, ay), potentials(x, y, m, n, out), deposit, loadMesh(massGrid), solveForces(), solvePotential(), fieldAt(x, y, out) → [a_x, a_y, Φ], h, x0, outside, meshMass }` | isolated boundaries; off-mesh particles feel the mesh monopole |
| `surfaceDensityPolar` | `(s, { nR = 32, nPhi = 64, Rmax })` → `{ Sigma (iR + nR·iPhi), R, nR, nPhi }` | M☉/kpc² |
| `fourierModes` | `(s, { mMax = 6, nR = 16, Rmax, Rin, Rout })` → `{ R, A, phase (iR·(mMax+1) + m), annulus, dominant, noise }` | annulus = ring-rms amplitude, so a wound spiral does not cancel |
| `spiralPitch` | `(s, { m, Rin, Rout, pMax = 50, dp = 0.25 })` → `{ m, p, pitchDeg, amplitude, trailing }` | log-spiral transform; > 0 = trailing for the sense of L_z |
| `rotationCurve` | `(s, { nR = 24, Rmax })` → `{ R, vphi, sigmaR, vc, mass, count }` | v_c from the current accelerations |
| `toomreQProfile` | `(s, { nR = 24, Rmax })` → `{ R, Q, Sigma, sigmaR, kappa, vc }` | noisy in sparse rings |
| `localDensity` | `(s, { cells = 128, extent = Rmax, out, work })` → Σ at each particle | for colouring |
| `toomreQ`, `swingArmNumber` | `(sigmaR, kappa, Sigma)`, `(R, kappa, Sigma, X = 2)` | stellar disk (3.36); X ≈ 1–3 |
| `epicyclicFrequency` | `(vc(R), R, dR)` → κ | numerical derivative |
| `exponentialDiskCircularSpeed` | `(R, M, Rd)` → km/s | Freeman 1970, unsoftened, untruncated M |
| `kuzminSurfaceDensity`, `kuzminRadialAcceleration` | `(R, M, a)`, `(R, M, a, softening = 0)` | exact also with Plummer softening |

Performance (Node 24, plain JS): 30k particles on a 128² mesh (256² FFT) run at about 100–130 steps/s, roughly 7 ms per FFT solve and 1–2 ms for the particles. 60k particles also exceed 30 steps/s (measured ≈ 100+). At `suggestedDt` (≈ 0.0013 in the default model), one rotation at 2 R_d takes about 150 steps. In a page, run 4–8 steps per frame.
