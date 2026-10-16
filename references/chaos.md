# Chaos and nonlinear dynamics

Status: generic tools (`core/chaos.mjs`, built on `odeStepper` in `core/ode.mjs` and `newtonSolve` in `core/linalg.mjs`). No ready-made recipe. Read [method selection](method-selection.md) for integrator choice and [rigid bodies](rigid-bodies.md) for the Dzhanibekov flip (instability, not chaos).

"Chaos" is a measurable property, not a look: a positive Lyapunov exponent on a bounded attractor or energy shell. A page that only shows a tangled trajectory has not shown chaos. Pick one quantitative diagnostic per page and verify it against an exact identity or a published value.

## Candidate models (pick by the question)

| Question | Model | Observables and benchmark |
|---|---|---|
| What is chaos? Route from order to chaos | Logistic map, bifurcation diagram | Period doubling at superstable R_n (R_2 = 3.4985616993), δ → 4.669201609… (Feigenbaum 1978), r∞ = 3.5699456718; λ = ln 2 at r = 4 |
| Butterfly effect, strange attractor | Lorenz (10, 28, 8/3) | λ = (0.9056, 0, −14.5723) (Sprott, Chaos and Time-Series Analysis, 2003); Σλ = −(σ+1+β) exactly; D2 ≈ 2.05 (Grassberger & Procaccia, PRL 50, 346 (1983)) |
| A fractal from a simple rule | Hénon map (1.4, 0.3) | λ1 = 0.41922, λ1 + λ2 = ln 0.3 exactly; D2 = 1.21, D0 = 1.261 (Russell, Hanson & Ott, PRL 45, 1175 (1980)) |
| Chaos in Hamiltonian systems; KAM tori | Standard map, Hénon–Heiles Poincaré section | Islands and chaotic sea; λ1 + λ2 = 0; λ ≈ ln(K/2) for large K; regular at H = 1/12, mostly chaotic at H = 1/6 |
| Double pendulum, driven oscillators | `doublePendulum`, `duffing` | Normal modes ω² = (g/l)(2 ∓ √2) at small angles; energy conserved (pendulum); strobed Duffing section |
| Mode locking, devil's staircase | Sine circle map | Rotation number = Ω at K = 0; plateaus (Arnold tongues) for K > 0 |
| Where does the system settle, and is it stable? | `findFixedPoint`, `findPeriodicOrbit`, `classifyFixedPoint` | Lorenz C± = (±√(β(ρ−1)), ±√(β(ρ−1)), ρ−1), Hopf at ρ_H = σ(σ+β+3)/(σ−β−1) ≈ 24.74; logistic 2-cycle multiplier −r² + 2r + 4 |

## Method

- **Maps:** objects `{ dim, map, jacobian }`; iterate with `iterateMap`, sweep with `bifurcationDiagram`. Feigenbaum: compute superstable parameters (`superstableParameters`), then `feigenbaumDelta`; they converge faster and more reliably than bifurcation points read off a diagram.
- **Flows:** objects `{ dim, f, jacobian }` that plug straight into `odeStepper` (adaptive Dormand–Prince, tight rtol). Hamiltonian chaos over long times: keep rtol ≤ 1e-10 and report the energy drift, or use `mechanicsStepper` when the force is position-only.
- **Lyapunov exponents:** maps `mapLyapunov` (QR of tangent products); flows `lyapunovSpectrum` (Benettin QR, all exponents) or `lyapunovExponent` in `core/ode.mjs` (largest only, also without a Jacobian). Discard a transient; show the running estimate (`history`).
- **Poincaré sections:** `poincareSection` with an event function (Brent on the dense output). Hénon–Heiles: section x = 0, direction +1, plot (y, py); start all orbits on one energy shell by solving for px.
- **Dimensions:** `correlationDimension` (pairs, O(N²), N ≤ a few thousand; Theiler window for flows) and `boxCountingDimension` (grid cells). The scaling range is a judgement: plot ln C against ln r and show the fit range.

## Verification

- Exact identities: Σλ = ⟨ln|det J|⟩ for maps (Hénon ln b, area-preserving 0), Σλ = mean divergence for flows (Lorenz −(σ+1+β)); tent map λ = ln μ to round-off; a linear flow gives the eigenvalue real parts.
- Published values with tolerances that reflect finite-time or finite-sample error (λ1 to ~0.02, D2 to ~0.05–0.1): cite the source in the page.
- Feigenbaum ratios converge geometrically; show successive δ_n, not one number.
- Fixed points against closed forms; eigenvalues against the characteristic polynomial; classification changes at the analytic bifurcation (Lorenz Hopf at ρ_H).
- Poincaré points on the section (|g| ≈ 0) and on the energy shell; an integrable test case (decoupled oscillators) gives the exact rotation.
- Sensitivity, shown honestly: two nearby initial conditions separate at ≈ e^{λ1 t} until saturation; the shadowing argument means individual long chaotic trajectories are not accurate, but their statistics are. Say so in the page.

## Representation

- Separate the space shown: phase space (attractor), parameter space (bifurcation diagram), section plane (Poincaré), or time series. Label axes with the variables.
- Bifurcation diagrams: tiny translucent points, many samples per parameter; mark R_n and r∞.
- Attractors: fade old trajectory segments; colour by time or by local divergence, not arbitrarily.
- Show λ1 converging (running estimate with a band), and two twin trajectories diverging on a log scale.

**Visible scope:** "Deterministic model, finite-time and finite-sample estimates; individual long chaotic trajectories are sensitive to rounding, their statistics are not." Name the parameters used.

## API

All in `core/chaos.mjs`. Maps return `{ name, dim, params, map(x), jacobian(x) }`; flows `{ name, dim, params, f(t, y), jacobian(t, y), energy? }`.

- Maps: `logisticMap(r)`, `tentMap(mu)` (μ < 2: at 2 floating-point orbits collapse), `henonMap({ a, b })`, `standardMap(K, { wrap })` (state [θ, p]; use `wrap: false` for orbit searches), `circleMap(Omega, K)` (lift; take mod 1 to draw).
- `iterateMap(m, x0, n, { transient })` → Float64Array(n·dim); throws on divergence.
- `bifurcationDiagram({ family, params, x0, transient, samples, coordinate, follow })` → `[{ param, values }]`.
- `superstableParameters({ g(r, x), xc, r0, r1, count ≤ 16 })` → Float64Array R_n. `feigenbaumDelta(points)` → ratios δ_n.
- `rotationNumber(m, x0, n, { transient })` → mean advance per iterate of a 1-D lift.
- `mapLyapunov(m, x0, { n, transient, record })` → `{ exponents, sum, history }`.
- Flows: `lorenz({ sigma, rho, beta })`, `rossler({ a, b, c })`, `duffing({ delta, alpha, beta, gamma, omega })` (state [x, v], non-autonomous), `henonHeiles()` (state [x, y, px, py], with `energy`), `doublePendulum({ m1, m2, l1, l2, g })` (state [θ1, θ2, ω1, ω2], with `energy`; Jacobian by central differences).
- `lyapunovSpectrum({ system, y0, t0, duration, interval, transient, h, rtol, atol })` → `{ exponents, sum, history }`; needs `system.jacobian`.
- `poincareSection({ system, y0, t0, section(t, y), direction, count, tmax, h, rtol, atol })` → `{ points, times, stepper }`.
- `correlationDimension(points, { rMin, rMax, bins, theiler })` → `{ dimension, r2, curve }`. `boxCountingDimension(points, { sizes })` → `{ dimension, r2, counts }`. Points are arrays of coordinates.
- `smallEigenvalues(J)` → `[{ re, im }]` for n ≤ 8, descending real part. Uses the shared balanced Hessenberg/Francis QR solver; forward accuracy depends on eigenvalue conditioning. It does not form a characteristic polynomial, which can badly corrupt repeated roots.
- `classifyFixedPoint(J, { kind: 'flow' | 'map', tol })` → `{ eigenvalues, type, stable, hyperbolic, unstableDimension }`.
- `findFixedPoint(flow, guess, { t, tol })` → `{ x, converged, residualNorm, ...classification }`.
- `findPeriodicOrbit(map, x0, period, { tol })` → `{ x, orbit, converged, multipliers, ...classification }`.
