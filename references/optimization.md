# Optimisation, automatic differentiation, quasi-Monte Carlo and geometric algebra

Four subject-independent modules: `core/optimize.mjs`, `core/autodiff.mjs`, `core/qmc.mjs` and `core/geometric-algebra.mjs`. Like the other core modules, they are plain ES modules that run in Node and the browser and have no dependencies. With fixed seeds they are deterministic. Each one is tested against exact answers in `tests/<module>.test.mjs`.

## Choosing a tool

| Problem | Use | Why / condition |
|---|---|---|
| Minimum of f(x) on an interval | `brentMinimize` (fast) or `goldenSection` (simplest) | f unimodal on [a, b]. Accuracy limit √ε·\|x\|, because f is flat within that distance |
| A few parameters; f noisy, non-smooth or without derivatives | `nelderMead` | n ≲ 10. There is no convergence theorem, so restart from the result to confirm |
| A smooth f with a gradient | `bfgs` (n ≲ 300) or `lbfgs` (any n) | Superlinear convergence. Get gradients from `autodiff.gradient` or supply analytic ones |
| Simple bounds l ≤ x ≤ u | `lbfgs` with `lower`/`upper` (projected), or `nelderMead` with bounds | Variables pinned at a bound are reported in `active` |
| Fit a model to data (residual vector) | `levenbergMarquardt`; `curveFit` (`core/data.mjs`) for y = model(x, p) | Bounds, robust losses and covariance. Check pulls on synthetic data |
| Outliers in the data | `levenbergMarquardt` with `loss: 'soft_l1'`/`'huber'`/`'cauchy'` and `fScale` ≈ the inlier σ | Robust covariances are only indicative |
| Many local minima (cluster energies, multimodal fits) | `differentialEvolution` (box) or `simulatedAnnealing`, then polish | Heuristics. Repeat with several seeds and report the spread |
| Derivatives of code you wrote | `autodiff`: `gradient` (reverse, cost independent of n), `derivative`/`jvp` (forward), `hessian`, `hessianVectorProduct` | Exact to rounding. Write the function using the module's operations |
| Derivatives of a black box | `finiteGradient` / `finiteHessian` | Truncation and rounding trade-off; steps are chosen automatically |
| Is my analytic gradient right? | `checkGradient` at several random points | Normwise relative error against 4th-order differences |
| Integrals over 2 to 21 dimensions (smooth integrand) | `qmcIntegrate` (scrambled Sobol) | Error ~N⁻¹ or better, against N⁻¹ᐟ² for Monte Carlo, with replicate error bars |
| Quality of a point set | `starDiscrepancy` (exact, small sets), `l2StarDiscrepancy` (any size) | |
| Rotations, boosts, EM field invariants | `geometric-algebra`: `rotor3`, `boostRotor`, `fieldBivector` | No gimbal lock. Rotors compose by multiplication |

## Conditions and verification advice

- **Convergence reports a stopping criterion.** Check `converged` and `reason`, the gradient or constrained optimality residual, and the local curvature. A small step or cost change alone does not certify stationarity or a minimum. For non-convex problems, run from several starts or seeds; agreement still does not prove global optimality.
- **Choose tolerances for the quantity.** `gtol` is absolute and bounds a gradient norm; step and cost tolerances have the method-specific relative/absolute meanings in the API. Rescale the objective and parameters. A large constant offset can obscure cost changes, so inspect residuals and gradients as well as the reported cost.
- **Finite-difference gradients** (the default in `bfgs`/`lbfgs` when no `gradient` is given) are accurate to about 1e-10 relative. Use `gtol` ≳ 1e-7 with them, or switch to `autodiff.gradient`. If f is undefined outside its domain, give an analytic or AD gradient.
- **Levenberg–Marquardt covariance.** With `absoluteSigma: true` the residuals must be (y − model)/σ with trusted σ. The default scales the covariance by χ²/dof, as `curveFit` does without `sigma`. Parameters pinned at a bound get NaN errors. Report the bound, not an error bar.
- **Stochastic global search** (`differentialEvolution`, `simulatedAnnealing`) finds the global minimum with high probability on benchmark functions, but it gives no proof. Use `polish: true` for full precision.
- **Autodiff** differentiates the branch your code takes. At kinks no classical derivative need exist: `abs(0)` uses a zero subgradient convention, while tied `min`/`max` choose a branch. Do not use `Math.*` on AD values; use the module's functions. Do not close over a `Dual` from an outer `derivative` call inside an inner one; use `secondDerivative`, `hessian` or `hessianVectorProduct` instead.
- **QMC.** A single unscrambled estimate has no sampling error bar. `qmcIntegrate` requires randomised replicates and rejects `scramble: 'none'`; use the sequence generators directly for deterministic points. Prefer R ≳ 10 replicates and n = 2^m for Sobol. The sample-variance error estimate is not a bound and does not include bias from finite coordinate resolution (Sobol uses 32 bits). Smoothness controls the convergence rate; measure it for the actual integrand.
- **Geometric algebra** conventions follow Doran & Lasenby (2003): STA signature (+ − − −), σ_k = e_k e0, F = E + IB with I = e0e1e2e3, active transformations x ↦ R x R̃, and c = 1. The rotors R and −R give the same transformation.

## Verified by

- **optimize.** Exact minima of Moré, Garbow & Hillstrom (1981, *ACM TOMS* 7, 17) problems: Rosenbrock #1, Beale #5, helical valley #7, Powell singular #13, and extended Rosenbrock #21 with n = 100. Also Himmelblau's four minima, the global minimum of Rastrigin in 5D (DE) and 2D (SA, three seeds), and bound-constrained minima (Rosenbrock with x ≤ ½ → (½, ¼); a box-clipped quadratic). The tests measure BFGS superlinear convergence, the golden-section ratio 1/φ per evaluation, the Brent/golden evaluation ratio and the finite-difference error orders 1, 2 and 4. The strong-Wolfe conditions are asserted. The covariance of a linear model is exactly σ²(AᵀA)⁻¹. On exponential data the results agree with `curveFit`. A robust fit with 20 % gross outliers is also checked.
- **autodiff.** Every elementary function is checked against analytic first and second derivatives, in forward, reverse and nested modes. The gradient of a 13-atom Lennard-Jones cluster equals minus the analytic forces to 1e-13. The tests also check the Rosenbrock Hessian, Hessian-vector products (H·v), a 6D Hessian against finite differences and the polar-map Jacobian.
- **qmc.** The tests reproduce the published first Sobol points. They check that every Joe–Kuo polynomial is primitive and every mᵢ is odd with mᵢ < 2ⁱ (Joe & Kuo 2008, *SIAM J. Sci. Comput.* 30, 2635). All 21 coordinates are stratified, and dimensions 1–2 form a (0, m, 2)-net, both unscrambled and scrambled. The star discrepancy is checked in closed form (1D and single points), and Warnock's L2 formula against quadrature. On Genz test integrands (Genz 1984) the measured RMSE slopes are about −1.2 (Sobol + Owen) and −1.0 (Halton + permutations), against −0.5 for Monte Carlo. Error bars cover the exact values.
- **geometric-algebra.** The tests check the axioms in five signatures, rotors against Rodrigues matrices and their composition, and the cross product. Boosts are checked against Einstein velocity addition, the collinear rapidity sum and invariance of the interval. The Thomas–Wigner angle for perpendicular boosts satisfies cos θ = (γ₁+γ₂)/(1+γ₁γ₂). The field tests cover F² = E² − B² + 2I E·B, its boost invariance, and the field of a moving charge (B = v × E, E⊥ × γ).

```js
import { lbfgs, levenbergMarquardt } from 'psim:core/optimize.mjs';
import { gradient, add, sub, mul, pow, square } from 'psim:core/autodiff.mjs';
import { qmcIntegrate } from 'psim:core/qmc.mjs';
import { rotor3, rotateVector3 } from 'psim:core/geometric-algebra.mjs';

// LJ cluster relaxation with exact AD gradients
const energy = x => { let E = 0; for (let i = 0; i < x.length; i += 3) for (let j = i + 3; j < x.length; j += 3) {
  const r2 = add(add(square(sub(x[i], x[j])), square(sub(x[i + 1], x[j + 1]))), square(sub(x[i + 2], x[j + 2])));
  const s6 = pow(r2, -3); E = add(E, mul(4, sub(square(s6), s6))); } return E; };
const x0 = [0, 0, 0, 1.2, 0, 0, 0.5, 1, 0.1];               // three atoms
const fg = x => gradient(energy, x);
const relaxed = lbfgs(x => fg(x).value, x0, { gradient: x => fg(x).gradient, gtol: 1e-8 });
if (!relaxed.converged) throw new Error(relaxed.reason);

const I = qmcIntegrate(x => Math.exp(-(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)), 3, { n: 4096 }); // I.value ± I.error
const v = rotateVector3(rotor3([0, 0, 1], Math.PI / 2), [1, 0, 0]); // [0, 1, 0]
```

## API

### `core/optimize.mjs`

Every minimiser returns `{ x, f, iterations, evaluations, converged, reason, … }`. Invalid input throws.

| Export | Result | Conditions |
|---|---|---|
| `goldenSection(f, a, b, { xtol, atol, maxIterations })` | `{ x, f, iterations, evaluations, converged, reason, bracket }` | f unimodal on [a, b]. Stops at bracket ≤ 2(xtol·\|x\| + atol), with defaults √ε and 1e-10·(b − a) |
| `brentMinimize(f, a, b, options)` | Same as `goldenSection` | Same conditions. Parabolic steps, typically 2–4× fewer evaluations |
| `nelderMead(f, x0, { step, xtol, ftol, maxIterations, maxEvaluations, adaptive, lower, upper })` | Adds `simplex` | Gao–Han adaptive coefficients. Converged when the simplex diameter ≤ xtol and the f spread ≤ ftol |
| `wolfeLineSearch(f, gradient, x, d, { f0, g0, alpha0, c1, c2, alphaMax, maxEvaluations })` | `{ alpha, x, f, g, evaluations, ok }` | d must be a descent direction (throws otherwise). 0 < c1 < c2 < 1. Once f stops changing beyond its rounding error (within 1e-12·\|f\|), the Hager–Zhang approximate Wolfe condition replaces Armijo, so ∇f can still be reduced |
| `bfgs(f, x0, { gradient, gtol, maxIterations, c1, c2, onIteration, gradientOptions })` | Adds `gradient`, `gradientEvaluations`, `inverseHessian` | f is C². Converged when ‖∇f‖∞ ≤ gtol. Dense n×n |
| `lbfgs(f, x0, { gradient, memory, gtol, maxIterations, lower, upper, c1, c2, onIteration })` | Adds `gradient`, `gradientEvaluations`, `active` | Strong Wolfe without bounds. Projected search with bounds. Converged when the projected gradient ≤ gtol |
| `levenbergMarquardt(residuals, p0, { jacobian, lower, upper, loss, fScale, absoluteSigma, maxIterations, xtol, ftol, gtol, lambda0 })` | `{ p, cost, residuals, chi2, dof, jacobian, cov, errors, active, iterations, evaluations, jacobianEvaluations, converged, reason }` | m ≥ number of free parameters. `cov` is null if JᵀWJ is singular. `loss`: `'linear'`, `'soft_l1'`, `'huber'` or `'cauchy'` |
| `simulatedAnnealing(f, x0, { lower, upper, seed, maxEvaluations, T0, Tend, step, polish })` | `{ x, f, evaluations, acceptance, converged, reason, polished? }` | Heuristic: `converged` refers only to the Nelder–Mead polish |
| `differentialEvolution(f, lower, upper, { seed, popSize, F, CR, strategy, maxGenerations, ftol, rtol, polish })` | `{ x, f, generations, evaluations, converged, reason, population, values, polished? }` | Finite box. Converged when the population's f spread ≤ ftol + rtol·\|f\| |
| `finiteGradient(f, x, { order, h, scale, f0 })` | Float64Array | order 1, 2 (default) or 4. Throws if f is not finite at a stencil point |
| `finiteHessian(f, x, { gradient, scale, h })` | Symmetric rows | With `gradient`: ~1e-10. From f alone: ~1e-7 |
| `checkGradient(f, gradient, x, { rtol, scale, h })` | `{ ok, maxError, worst, analytic, numeric, errors }` | Normwise relative error. Check at several points |

### `core/autodiff.mjs`

| Export | Result | Conditions |
|---|---|---|
| `Dual`, `dual(x, dx)`, `Var` | Forward-mode dual numbers (nestable) and reverse-mode tape nodes | Build them through the drivers below |
| `add sub mul div neg square pow sqrt cbrt exp expm1 log log1p sin cos tan asin acos atan atan2 sinh cosh tanh abs hypot min max sum dot` | The same type as the inputs (number, Dual or Var) | `pow(a, b)` with a differentiated b needs a > 0. `abs`, `min` and `max` are one-sided at ties |
| `value(x)` | The primal number | Use it for branching and output |
| `derivative(f, x)` | `{ value, derivative }` | One forward pass |
| `secondDerivative(f, x)` | `{ value, first, second }` | Hyper-dual pass |
| `jvp(f, x, v)` | `{ value, jvp }` | Directional derivative. Scalar or vector f |
| `gradientForward(f, x)` | `{ value, gradient }` | n passes |
| `jacobian(F, x)` | `{ value, jacobian }` | n passes. m rows |
| `gradient(f, x)` | `{ value, gradient }` | Reverse mode: f must return a scalar. Memory ∝ number of operations |
| `hessianVectorProduct(f, x, v)` | `{ value, gradient, hvp }` | Forward-over-reverse, one pass |
| `hessian(f, x)` | `{ value, gradient, hessian }` | n passes |

### `core/qmc.mjs`

| Export | Result | Conditions |
|---|---|---|
| `sobol(dim, { scramble, seed, skip })` | `{ dim, point(i), next(), points(n) }`, coordinates in [0, 1) | dim ≤ 21, use n = 2^m. `scramble`: `'none'`, `'shift'` or `'owen'` (hash-based nested uniform). At most 2³² points |
| `halton(dim, { scramble, seed, skip })` | Same interface | `scramble`: `'none'`, `'permute'` or `'shift'`. Unscrambled point 0 is the origin |
| `qmcIntegrate(f, dim, { n, replicates, sequence, scramble, seed, lower, upper })` | `{ value, error, n, replicates, estimates }` | `sequence`: `'sobol'` (n = 2^m), `'halton'` or `'random'`. At least 2 replicates (≳ 10 recommended). Throws on a non-finite integrand |
| `starDiscrepancy(points, { maxWork })` | Exact D* | Cost ~(N+1)^d·N·d. Throws above `maxWork` |
| `l2StarDiscrepancy(points)` | Exact L2 star discrepancy (Warnock) | O(N²d) |
| `SOBOL_DIRECTIONS`, `SOBOL_MAX_DIM` | The Joe–Kuo table for dimensions 2–21, and 21 | |

### `core/geometric-algebra.mjs`

| Export | Result | Conditions |
|---|---|---|
| `cliffordAlgebra(p, q, { firstIndex })` | An algebra object: `zero scalar vector blade basis pseudoscalar add sub scale neg mul wedge inner scalarProduct commutator grade grades reverse involute conjugate dual norm2 inverse exp sandwich equals toString`, plus `n p q size names gradeOf` | p + q ≤ 8. Multivectors are Float64Array(2ⁿ) indexed by blade bitmask. `inner` is the left contraction. `inverse` throws for null divisors. `exp` is closed-form when a² is a scalar, otherwise a series |
| `CL3`, `STA` | Cl(3,0) with e1–e3, and Cl(1,3) with e0–e3 (+ − − −) | |
| `rotor3(axis, angle)` | CL3 rotor exp(−I n̂ θ/2) | Right-handed. Compose as R₂R₁ |
| `rotor3Between(a, b)` | Smallest rotation taking â to b̂ | Throws for antiparallel vectors |
| `rotateVector3(R, v)`, `rotorToMatrix3(R)`, `cross3(a, b)` | [x, y, z]; 3×3 rows; a × b = −I(a∧b) | |
| `boostRotor(v)` | STA rotor with R e0 R̃ = γ(e0 + v) | \|v\| < 1, c = 1 |
| `fourVector([t, x, y, z])`, `fourVelocity(v)`, `spatialVelocity(u)` | STA vectors; 3-velocity | |
| `fieldBivector(E, B)`, `fieldComponents(F)`, `fieldInvariants(F)` | F = E + IB; `{ E, B }`; `{ EsqMinusBsq, EdotB }` | Natural units. Use F ↦ R F R̃ to transform the field |
