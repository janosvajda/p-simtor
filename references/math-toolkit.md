# Reusable mathematical tools

This is the single catalogue of general mathematics in `core/`: read it before writing a solver or numerical primitive. These are working Node/browser ES modules, independent of the scientific subject. Import them in a composed physics module; run that module through `psim verify-module` and build with the existing bundler.

Choose tools from the model's mathematical structure. Their numerical checks do not establish that your chosen physical model is applicable. Use suitable units, bound UI inputs by measured cost, and verify the new model's own observables.

## Numerical primitives: `core/numerics.mjs`

Constants (`core/constants.mjs`): CODATA 2018, IAU, standard gravity `gn`, Earth/Moon `GM` and radii, gas constant `R`, dry-air `Mair`, `gammaAir` and ISA sea-level values.

Seeded `rng(seed)` (with `.normal()`) and `rngFor(seed, i)` (an independent reproducible stream per item, to re-trace one particle alone), 1D `fft`/`fftFreq`, regularised incomplete gamma `gammaP`/`gammaQ` (χ² p-values), `rk4Step`, `dopri5`, quadrature (`gaussLegendre`, `adaptiveSimpson`, `integrateToInfinity`), `brent` roots, special functions (`lgamma`, `gamma`, `erf`, `besselI`, `besselJ`, `besselJZeros`, `laguerre`, `complexYlm`, `hydrogenRadial`, `legendreP`, `normalizedLegendre`, `realYlm`, `ellipticK`), `solveTridiagonal`, `eigenTridiagonal`, `linearFit` and `meanWithError` (blocking error bars for correlated series).

**Bessel functions.** `besselJ(n, x)` (integer n ≥ 0, any real x; relative accuracy even where J_n is tiny) and `besselJZeros(n, count)` give the exact modes of circular and cylindrical problems, for example ω_{n,s} = c·j_{n,s}/R for a clamped membrane, and serve as benchmarks for grid solvers on a disk.

**Atomic functions.** `complexYlm(l, m, θ, φ)` → `{ re, im }` (Condon–Shortley phase; `realYlm` gives the real p_x/p_y/d… combinations) and `hydrogenRadial(n, l, r, { Z, a0 })` (normalised so ∫R²r²dr = 1; point nucleus, non-relativistic) give exact one-electron orbitals ψ = R·Y. Use the reduced-mass a₀ when spectroscopic precision matters.

**Two eigensolvers, by matrix structure.** `eigenTridiagonal(diag, off)` (implicit QL, O(n²), unsorted values with unit vectors) handles 1D operators with thousands of points; `eigenSymmetric` (below) handles any dense symmetric matrix of modest size.

**Two error-bar estimators for correlated samples.** `meanWithError` (blocking) and `effectiveSampleSize` (integrated autocorrelation, below) estimate the same uncertainty by independent methods. Use either; agreement between them is useful evidence. Use `autocorrelationTime` when you must thin a chain (for example before `ksTest`).

## Evolution of finite-dimensional systems: `core/ode.mjs`, `core/particles.mjs`

| Export | Use and result | Conditions |
|---|---|---|
| `odeStepper({ f, y0, t0, method, h, rtol, atol, hmax, events })` | `'rk4'` fixed step or `'dopri5'` adaptive; `step`, `advanceTo`, `sampleAt(t)` (cubic Hermite dense output), step/reject counts | Non-stiff ODEs; dense output lets frames sample any time without changing the step sequence. `events: [{ g(t, y), direction, terminal }]` locates sign changes of g (ground contact, collisions, thresholds) and can stop there; results in `events`, `stopped` |
| `mechanicsStepper({ accel, q0, v0, t0, h, method })` | Symplectic `'leapfrog'` or `'yoshida4'` for ẍ = a(t, q); `sampleAt(t, outQ, outV)`; negative h runs backwards | Accelerations independent of velocity; fixed h (symplectic only at fixed step) |
| `refinementOrder(run, h)` | Measured convergence order from runs at h, h/2, h/4 | Asymptotic regime; report the measured order, not the nominal one |
| `lyapunovExponent({ f, jvp, y0, duration, interval, transient })` | Largest Lyapunov exponent: linearised equations when `jvp(t, y, w)` = J·w is given, otherwise a twin trajectory; `{ lambda, history }` | A finite-time estimate for one start: report how `history` settles and compare both methods |
| `linearChain({ rates, links })` | Exact (Bateman) solution of first-order chains with branching: decay series, sequential reactions, compartments; `at(t, x0)`, `matrix` | Distinct rates (equal rates: `expm`); links go from lower to higher index |
| `createParticles({ bodies, forces, method, h, rtol, atol, box })` | Interacting bodies from force terms: `gravity`, `coulomb` (uses `body.charge`), `lennardJones` (shifted-force cutoff by default; per-body `sigma`/`epsilon` with Lorentz–Berthelot mixing), `uniformField({E, B})`, `centralBody` (optional 1PN with `energyCorrection`), `drag`, `springs`, `externalPotential`. `method: 'auto'` picks Yoshida-4 for position-only forces, otherwise Dormand–Prince with the dissipated work integrated alongside. `box: [Lx, Ly, Lz]` makes axes periodic (0 = open; `[Lx, Ly]` for 2D): positions stay unwrapped, short-range forces use the minimum image, `wrapped()` folds positions for drawing | `energy()`, `energyBalanceError()` (includes work done by drag/fields), `momentumDrift()` (null when external forces act). Short-range forces use cell lists (cost ∝ N); gravity and Coulomb are all-pairs and refuse a periodic box (that needs Ewald summation) |
| `forEachPairWithin(x, n, rc, box, fn)` | Every pair closer than rc once, with the minimum-image displacement, by cell lists | For your own short-range interactions; rc ≤ half the box |

## The Sun in the sky: `core/sun.mjs`

`sunEquatorial(date)` (declination, right ascension, equation of time), `sunPosition(date, lat, lon)` (altitude, azimuth from north), `sunDirection(date, lat, lon)` ([east, north, up] unit vector, ready for `createMeshRenderer` shadows) and `sunriseSunset(date, lat, lon)` (with polar day/night). Meeus ch. 25 low-precision theory, ~0.01°; checked against Meeus' example, 2024 equinox/solstice, equation-of-time extremes and published London sunrise/sunset. No refraction except the standard −0.833° horizon.

## Grid evolution with walls, sources and waves: `core/pde.mjs`

Uses the storage and boundary conventions of `fields.mjs` below (axis 0 fastest; cell-centred Neumann).

| Export | Calculation | Limits and conventions |
|---|---|---|
| `diffusion1D({ u0, D, dx, dt, bc, left, right, source })` | Crank–Nicolson, second order, unconditionally stable; `total()`, `inflow`, `setBoundary(side, wall)` | Walls: `bc` for both, or `left`/`right` = `{ type: 'value', value }` (u held at the wall face, e.g. a saturated surface) or `{ type: 'flux', value }` (inward flux; 0 = insulated); values may be functions of t; switch a wall mid-run (a solid used up). `total() − total₀ = inflow` (+ source) to round-off. With D·dt/dx² ≫ 1 the sharpest grid modes decay slowly (CN's amplification → −1): use smaller dt when early transients or exact steady states matter. Source-free periodic diffusion in any dimension: `diffusionPeriodic` |
| `waveND({ u0, v0, shape, spacing, bc, mask, maskBc, c, dt })` | Leapfrog with `laplacianND` in any dimension; `energy()` is the discrete energy leapfrog conserves to round-off | Throws unless c·dt·√(Σ1/dx²) ≤ 1; constant c; optional irregular domain (see masks below) |
| `stationaryStates1D({ V, dx, count, mass, hbar })` | Lowest levels and normalised states of −ħ²/2m ψ'' + Vψ = Eψ | ψ = 0 beyond both ends; O(dx²); domain must contain the states |

## Matrices and equations: `core/linalg.mjs`

Dense matrices are arrays of finite rows (ordinary arrays or typed rows); vectors may be arrays or `Float64Array`. Inputs are preserved. These routines use real arithmetic.

| Export | Use and result | Conditions |
|---|---|---|
| `matrixVector(A, x)` | Rectangular matrix-vector product, compensated row sums | Consistent finite dimensions |
| `solveLinear(A, b, options?)` | Square Ax = b, row equilibration and partial pivoting; returns x | Nonsingular at the specified relative `pivotTolerance` |
| `leastSquares(A, b, options?)` | Column-scaled, pivoted Householder QR; returns `{ x, residualNorm, rank, permutation }` | Rows ≥ columns, full column rank; `rankTolerance` controls numerical rank rejection |
| `eigenSymmetric(A, options?)` | Cyclic Jacobi; returns ascending `values`, matching unit `vectors` (`vectors[k]` belongs to `values[k]`), `sweeps`, `offDiagonalNorm` | Real symmetric, modest dense size; relative Frobenius `rtol`, symmetry tolerance and sweep budget |
| `cholesky(A)` | Lower-triangular L with A = LLᵀ (arrays of rows) | Symmetric positive definite; throws with the failing pivot otherwise (e.g. a noisy correlation matrix) |
| `eigenGeneralized(A, B)` | A v = λ B v: ascending `values`, B-normalised `vectors` (vᵀBv = 1) | A symmetric, B SPD. The variational (GEVP) method for correlation matrices: A = C(t), B = C(t₀), largest λ → ground-state energy −ln λ/(t − t₀) |
| `expm(A)` | Matrix exponential (scaling and squaring, degree-13 Padé): x(t) = e^{At}x₀ for any linear x' = Ax | Normwise accuracy; check tiny components against `linearChain` or a structured solution |
| `conjugateGradient(A, b, options?)` | Matrix-free or dense SPD solve; returns `{ x, converged, iterations, residualNorm }` | Symmetric positive-definite operator; optional SPD inverse preconditioner |
| `newtonSolve(f, x0, options?)` | Damped nonlinear root search; returns x, residual, residual norm, convergence flag, iterations and reason | Scaled equations; pure callbacks; analytic `jacobian` or central differences using `xScale` |

For iterative solves, **check `converged` before using the result**. CG checks the true residual against max(atol, rtol·‖b‖), and returns false at its work limit. Positive curvature checks detect some invalid operators; they do not prove that an arbitrary callback is SPD. A periodic Laplacian's constant null mode needs a gauge/constraint or a compatible restricted space.

Newton success bounds the residual, not distance to a root, uniqueness or physical validity. Provide an analytical Jacobian when finite differences would cross a model's domain boundary. Singular systems throw; QR does not manufacture a minimum-norm answer for rank-deficient data. Dense work is cubic; use a matrix-free operator for larger sparse problems.

```js
import { solveLinear, conjugateGradient } from 'psim:core/linalg.mjs';
const x = solveLinear([[2, 1], [1, 3]], [1, 2]); // [0.2, 0.6]
const result = conjugateGradient([[2, 1], [1, 3]], [1, 2]);
if (!result.converged) throw new Error(`Unresolved linear solve: ${result.residualNorm}`);
```

Independent tests cover dimensional rescaling, near-dependent QR columns, exact discrete eigenmodes, eigenvector residuals/orthogonality, sparse manufactured solutions, preconditioning and nonlinear roots. Algorithms and conditions: [LAPACK factorisations](https://www.netlib.org/lapack/lug/node38.html), [QR least squares](https://www.netlib.org/lapack/lug/node40.html), [CG operator and stopping conditions](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.cg.html).

## Multidimensional fields: `core/fields.mjs`

Flat field storage uses **axis 0 fastest**: shape `[nx, ny, nz]` means index `i + nx*(j + ny*k)`. `spacing` is a positive scalar or one value per axis. Operators accept any-dimensional shapes; FFT routines need a power of two on every axis. Tests exercise one, two and three dimensions.

| Export | Calculation | Limits and conventions |
|---|---|---|
| `laplacianND(u, grid, out?)` | Compact second-order finite-difference ∇²u | Per-axis periodic, zero-exterior Dirichlet, or copied-edge-cell Neumann ghosts; output must not overlap input |
| `gradientND(u, grid)` | Central derivative per axis | Returns flat components; same ghost rules |
| `divergenceND(components, grid)` | Central divergence | Composition with `gradientND` has a wider stencil than `laplacianND` |
| `fftND(re, im, shape, inverse?)` | In-place complex separable radix-2 transform | Separate non-overlapping storage; inverse divides by total grid points |
| `solvePoissonPeriodic(source, options)` | (-∇² + screening²)u = source | Periodic FFT; spectral or finite-difference symbol; unscreened source must have zero mean within declared tolerance, solution mean fixed to zero |
| `diffusionPeriodic(config)` | ∂u/∂t = D∇²u with no source | Exact Fourier propagation of the selected spatial operator, constant non-negative D; periodic boundaries |
| `schrodingerND(config)` | iħ∂ψ/∂t = [−ħ²∇²/(2m) + V]ψ | Periodic Strang FFT, static real V, isotropic positive mass; optional `absorb: { width, strength }` (complex absorbing layers at the edges so outgoing waves leave instead of wrapping; `absorbed()` = probability removed). Without `absorb` no automatic normalisation or boundaries |
| `contourLines(values, { shape, spacing, origin }, level)` | Marching-squares level sets of a 2D field as polylines `[{ points, closed }]` | O(dx²) for smooth fields; a drawing aid |
| `densityLevel(values, fraction)` | The density level whose super-level set holds `fraction` of the total ("the 90 % surface") | Uniform grid, non-negative values; cell resolution: refine to check |

For `laplacianND`, `gradientND` and `divergenceND`, `grid` is `{ shape, spacing, bc }`; `bc` may be a string or per-axis array.

**Irregular domains.** Add `mask` (one entry per cell, 1 = inside) and `maskBc` (`'dirichlet'`: zero outside, such as a clamped edge; `'neumann'`: zero flux, such as a free or insulated edge) to the grid. The operators then act only on masked cells and return 0 outside, and remain symmetric, so `waveND` energy stays exactly conserved and `conjugateGradient` applies. A curved boundary becomes a staircase of cells, which is **first order in dx**: the clamped disk's fundamental is low by 3.9 %, 2.1 %, 0.9 % and 0.46 % at 20, 40, 80 and 160 cells across. Report the measured error against an exact mode (for a disk: `besselJZeros`), not the nominal second order of the stencil. Neumann means cell-centred zero face flux and conserves the plain cell sum. This differs from a node-centred mirror stencil with half-weight endpoints. Dirichlet means zero values at the ghost grid points, not that the supplied edge samples are fixed walls.

Poisson options include `screening` (default 0), `discretization` (`spectral` default or `finite-difference`), and source-compatibility `rtol`/`atol`. It returns `{ u, sourceMean, discretization }`. Do not remove a nonzero source mean as an unreported physical change. Positive screening removes the constant-mode null space.

Diffusion takes `{ u0, shape, spacing, D, dt, t0?, discretization? }`, and returns t, u, `step()`, `advanceTo(t)` and `total()`. It evaluates each target from the initial coefficients; a change of frame schedule cannot change the field at a given time. There is no temporal truncation error for this linear problem, but spatial resolution still matters. Spectral truncation alone does not guarantee pointwise positivity.

Schrödinger takes `{ re0, im0?, V?, shape, spacing, dt, mass?, hbar?, t0? }` and returns t/re/im, `step()`, `advanceTo(target, maxSteps?)`, `probability()` and `norm()`. Norm includes cell volume. `advanceTo` returns false when its fixed-step budget is exhausted; it can finish just after the requested time. Check spatial resolution, wraparound and second-order splitting error independently; norm conservation alone does not establish physical accuracy.

```js
import { laplacianND } from 'psim:core/fields.mjs';
import { conjugateGradient } from 'psim:core/linalg.mjs';
const grid = { shape: [16, 16], spacing: [0.1, 0.1], bc: 'dirichlet' };
const A = (x, out) => {
  laplacianND(x, grid, out);
  for (let i = 0; i < out.length; i++) out[i] = -out[i];
};
// b is your finite source on the interior nodes; the zero ghosts make -∇² SPD.
const solve = b => conjugateGradient(A, b);
```

Tests compare the ND transform to a direct 3D DFT, verify Parseval, polynomial derivatives, flux conservation, Fourier-mode solutions, Poisson residuals/spatial order, diffusion decay and quantum phase/norm/splitting order. Transform conventions: [NumPy's DFT definitions](https://numpy.org/doc/stable/reference/routines.fft.html). The evolution benchmarks follow independently from the displayed equations' Fourier modes.

## Probability and equilibrium: `core/probability.mjs`

| Export | Use | Conditions/interpretation |
|---|---|---|
| `logSumExp(logWeights)` | Stable ln Σ exp(wᵢ), preserving small corrections with `log1p` | Finite or −∞ inputs; empty/all-impossible support returns −∞ |
| `normalizeWeights(weights)` | Stable probabilities from non-negative weights | Positive total required, scaled to avoid sum overflow |
| `normalizeLogWeights(logWeights)` | Stable normalized probabilities from log weights | At least one finite weight; max shifting and division |
| `shannonEntropy(p, {base}?)` | Discrete information entropy, nats by default or bits at base 2 | Normalized p, 0 ln 0 = 0; not automatically thermodynamic entropy |
| `relativeEntropy(p, q)` | KL divergence in nats | Matching normalized supports; p > 0 at q = 0 gives +∞ |
| `covarianceMatrix(samples, {ddof}?)` | Mean and full joint covariance | Rows are joint samples; ddof=1 sample/default or 0 population; not a standard error for correlated draws |
| `canonicalEnsemble(energies, options)` | Finite-spectrum equilibrium populations and moments | Positive energy-scale `kT`, temperature-independent levels, optional positive integer degeneracies |

For canonical equilibrium, `kT` means **k_B T in the same energy units as the levels**. Outputs are `probabilities`, `meanEnergy`, `energyVariance`, `freeEnergy`, `entropy` (S/k_B) and `heatCapacity` (C/k_B). Level populations include degeneracy: Z = Σgᵢexp(−Eᵢ/kT), so treating a degenerate level as one microstate would give the wrong entropy.

The routine reports `referenceEnergy` and `logShiftedPartition` for Z_shift = Σgᵢexp(−(Eᵢ−E_ref)/kT). This preserves probabilities and free energy even when direct exponentials all overflow/underflow. If needed, ln Z = ln Z_shift − E_ref/kT, provided that unshifted number is representable. Zero temperature is a separate ground-state limit; do not pass kT=0. Missing degrees of freedom and an infinite spectrum are not supplied by this routine.

```js
import { canonicalEnsemble } from 'psim:core/probability.mjs';
const equilibrium = canonicalEnsemble([0, 1, 2], {
  kT: 0.7, degeneracies: [1, 2, 1],
});
// Entropy is per modeled system in units of k_B, not an arbitrary visual disorder score.
```

The probability tests use exact distributions, covariance identities, binomial-level partition functions, thermodynamic derivatives, energy-offset invariance and cold/hot limits. They retain tiny entropy and nearly equal-distribution KL signals. Sources: [Blanchard–Higham–Higham on stable log-sum-exp/softmax](https://eprints.maths.manchester.ac.uk/2765/3/paper_nomarks.pdf), [MIT canonical partition functions and degeneracy](https://ocw.mit.edu/courses/5-62-physical-chemistry-ii-spring-2008/bda82786b6014d3f468798339f658629_03_562ln08.pdf).

Infinite KL or impossible-support limits are mathematically meaningful, but cannot be a finite measured value in a passing verification report. Identify the support problem or limit explicitly. All three modules reject invalid dimensions, domains and unresolved work rather than manufacture a prediction; rescale problems whose outputs or intermediates exceed floating-point range.

## Stochastic dynamics and sampling: `core/stochastic.mjs`

Every routine takes `rand = rng(seed)`, so runs are reproducible.

| Export | Use | Conditions/interpretation |
|---|---|---|
| `gillespie({ x0, reactions, rand })` | Exact event-by-event simulation of counts; reactions `{ rate: x => propensity, change }`; `advanceTo(t)` gives the state at t | Well-mixed populations; elapsed time is physical |
| `eulerMaruyama({ f, g, x0, dt, rand })` | Itô SDE dX = f dt + g dW, diagonal noise | Strong order ½, weak order 1: refine dt |
| `langevinBAOAB({ x0, v0, force, mass, gamma, kT, dt, rand })` | Underdamped Langevin; samples the canonical distribution with O(dt²) configurational error; `heat` = energy given by the bath, so H − H₀ − heat = O(dt²) checks a thermostatted run | Physical time when γ is a physical friction rate; otherwise a sampler |
| `metropolis({ logP, x0, rand, stepSize, propose })` | Random-walk sampling of any density known up to a constant; `acceptance()` | Symmetric proposals; steps are not physical time; discard warm-up |
| `metropolisSites({ state, propose, deltaLogP, rand, order })` | Single-site Metropolis for discrete states (Ising/Potts spins, occupations); `deltaLogP(i, value, state)` = ln P change (−ΔE/kT); one `step()` = one sweep | Symmetric proposals; sweeps are not physical time; verify on a small lattice by exact enumeration |
| `autocorrelationTime(xs)`, `effectiveSampleSize(xs)` | τ (1 for independent draws), n/τ, mean and its standard error | Needs a series much longer than τ |

## Working with data: `core/data.mjs`

| Export | Use | Conditions/interpretation |
|---|---|---|
| `interpolate(xs, ys, { method, extrapolate })` | `'monotone'` (Fritsch–Carlson: no invented extrema) or `'linear'`; NaN outside the data unless `extrapolate` | Strictly increasing xs; interpolation is not new data, so label it |
| `histogram(samples, { bins, range, density })` | Counts, density values and Poisson errors; `inRange` | Errors assume independent samples |
| `ksTest(samples, cdf)` | Kolmogorov–Smirnov D and asymptotic p against a model CDF | Independent samples: thin correlated chains first |
| `createHistogram({ bins, range })` | Incremental: `add(x, w)` per event; `counts`, `errors()` (√Σw²), `density()`, `underflow`/`overflow`/`nan`, `edges`, `centers` | For live Monte Carlo pages (events arrive one by one); `histogram` is the batch version |
| `chiSquareTest(observed, expected, { ddof, errors })` | Pearson χ² (or Σ((o−e)/σ)² with `errors`) → `{ chi2, dof, p, z, skipped }` | Expected counts ≳ 5 per bin; bins with expected ≤ 0 are skipped and counted |
| `jackknife(samples, estimator, { bins })` | Blocked leave-one-out → `{ value, error, bias, blocks }`; samples are numbers or vectors measured together; the estimator may return a vector | Nonlinear functions of means (ratios, logs, effective masses, fits). Blocks ≫ autocorrelation time. A NaN on any block gives a NaN error: report "unresolved" |
| `bootstrap(samples, estimator, { resamples, bins, rand })` | Block bootstrap → `{ value, error, quantiles(p) }` | Same conditions; asymmetric intervals from `quantiles` |
| `jackknifeRefit(samples, fit, { bins })` | Jackknife of an estimator that needs the samples themselves (a refit, a median): fit re-run with each block left out | B + 1 fits; blocks ≫ autocorrelation time |
| `curveFit(model, xs, ys, p0, { sigma, jacobian })` | Levenberg–Marquardt nonlinear least squares → `{ p, errors, cov, chi2, dof, p_value, converged }`; without sigma, errors from the scatter | Start values within the basin (frequency of an oscillation within ~20 %); check pulls on synthetic data |
| `meanSem(xs)`, `pull(measured, error, expected)` | Mean ± standard error of independent samples; (measured − expected)/error | For correlated series use `meanWithError` |

## Relativistic kinematics and event sampling: `core/kinematics.mjs`

4-vectors are `[E, px, py, pz]` with c = 1 and metric (+, −, −, −); units are the caller's (state them).

| Export | Use | Conditions |
|---|---|---|
| `dot4(p, q)`, `mass(p)`, `add4(...ps)` | Invariant products, invariant mass, sums | `mass` returns 0 for slightly negative p² from round-off |
| `boost(p, β)`, `restFrameBoost(P)` | Active boost by velocity β; the β taking P to its rest frame (`boost(P, restFrameBoost(P))` has zero momentum) | \|β\| < 1 |
| `randomUnitVector(rand)`, `randomRotation(rand)`, `orthonormalFrame(n)` | Uniform direction on S²; uniform SO(3) rotation (3×3 rows); right-handed frame [e₁, e₂, n] | |
| `twoBodyDecay(M, m1, m2, rand)` | Isotropic two-body decay in the rest frame, momentum from the Källén function | M ≥ m1 + m2 |
| `rambo(n, √s, rand, masses?)` | Uniform n-body phase space (Kleiss, Stirling & Ellis 1986) → `{ momenta, weight }`; weight is the phase-space volume in the (2π) convention (constant when massless) | n ≥ 2; exact 4-momentum conservation |
| `acceptReject(sample, weight, wmax, rand)` | Unweighted events from weighted ones → `{ value, tries }` | Throws if a weight exceeds `wmax` (a wrong bound never biases silently) |

Benchmarks: massless 3-body phase space gives x = 2E/√s distributed as 2x on [0, 1]; 2-body phase-space volume 1/(8π); muon-decay spectra in [weak interaction](weak-interaction.md).

## Higher dimensions and the 3-sphere: `core/projection.mjs`

| Export | Use | Conditions |
|---|---|---|
| `rotation4({ xy, xz, xw, yz, yw, zw })`, `apply4(R, p)` | 4D rotation from plane angles (row-major 4×4) | Orthogonal, det +1 (tested) |
| `project4(p, { mode: 'perspective' \| 'orthographic', eye })` | 4D → 3D for ui/gl.mjs; perspective scales by eye/(eye − w) | Perspective throws for points behind the 4D eye; connectivity exact, lengths not |
| `stereographic(q)`, `inverseStereographic(x)` | S³ ↔ ℝ³ from the pole −1 (conformal) | Points near −1 go far out: clip or use `axisAngle` |
| `axisAngle(q)` | S³ → ball of radius π (distance from the identity) | Every point finite; antipodal shell points are identified |
| `hopfMap(q)`, `hopfFibre(p, n)` | S³ → S² and the U(1) circle above p (left multiplication by exp(itσ₃)) | Unit quaternions |

## Lattice gauge fields and the QCD coupling: `core/lattice-gauge.mjs`, `core/qcd.mjs`

Subject-specific but reusable: SU(2) sampling, `wilsonLoops` (all R×T at once), `polyakovLoop`, `plaquetteField`, `apeSmear`, `strongCouplingPlaquette(β, d)`, `siteCoords`/`siteIndex`, allocation-free `qmulInto`, gauge tools (`gaugeTransform` with a given g, `centreTransform`, `landauGauge`, `haarSU2`, `qexp`/`qlog`), topology (`topology`: clover s(x), q(x), Q; `bpstInstanton` benchmark; cooling and ranged sweeps via `sweep` options); the MS-bar running coupling `runAlphaS({ loops, thresholds })`, `betaCoefficients(nf, { group })`, `alphaS1Loop`, `exactTwoLoop`. Usage, costs and required checks: [particle physics](particle-physics.md).

## Matter, spin and gas: `core/matter.mjs`, `core/spin.mjs`, `core/probability.mjs`

| Export | Use | Conditions |
|---|---|---|
| `betheStoppingPower(βγ, { mass, z, material })`, `csdaRange(T, …)`, `MATERIALS` | Mean energy loss (MeV cm²/g) and range (g/cm²) of heavy charged particles in liquid hydrogen, water, silicon, lead, air (PDG data; minimum ionisation within 0.1 %) | 0.1 ≲ βγ ≲ 1000; mean loss only (no straggling); not for e± |
| `ketAlong`, `rotate`, `bloch`, `probability`, `measure` | Spin ½ / qubit states, precession, Born probabilities, projective measurement | Pure states; see [quantum](quantum.md) for the semiclassical beam model |
| `maxwellSpeed(T, m)`, `effusiveSpeed(T, m)` | Gas speed laws (∝ v² and ∝ v³ e^{−mv²/2kT}) with pdf, cdf, exact sampler, mean | m in kg (`ATOMIC_WEIGHT[el]·u·1000`) |

Fields: `probabilityCurrent(re, im, { shape, spacing, mass, hbar })` from `core/fields.mjs` gives j = (ħ/m)Im(ψ*∇ψ) per axis, spectrally (the flux check for absorbing layers and transmission).

ODE notes: `odeStepper(...).advanceTo(t, Infinity, { exact: true })` ends exactly at t; `s.hmax` is mutable (cap steps inside a region so events there are not stepped over).

