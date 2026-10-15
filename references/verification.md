# Verification

Every page carries evidence for its trust tier. Derive checks from the chosen equations, observable definitions and scientific claims. The catalogue below gives examples, not a closed list or a universal test suite.

## Check catalogue

| Check | How | Applies to | Typical pass threshold |
|---|---|---|---|
| Energy conservation | \|E(t) − E₀\| / \|E₀\| over the run | Hamiltonian ODEs, NVE MD, Boris in pure B | ≤1e-8 few-body; ≤1e-4 MD (fluctuates, must not drift) |
| Momentum / angular momentum | Normalised by Σm\|v\| and Σm\|r\|\|v\| | Isolated systems | ≤1e-10 (round-off level) |
| Convergence order | Errors at h, h/2, h/4: p = log₂(e₁/e₂) | Every numerical integrator | p within ±0.3 of theoretical order |
| Predictability horizon | Time until h/2 and h/4 runs disagree by 1e-3 × system size | Chaotic systems | Report it; don't hide it |
| Exact solution | Compare against closed form | Kepler, helix, harmonic oscillator, free wavepacket | Set by the method's known error |
| Published benchmark | Reproduce a literature outcome | Figure-eight period, Pythagorean ejection, LJ argon equation of state | Matches literature within stated precision |
| Known discretisation error | Compare with the method's exact discrete behaviour (Boris radius r√(1+(ΩΔt/2)²)) | Methods with analysable errors | ~round-off |
| Normalisation | ∫\|ψ\|² dV = 1 | Quantum states | ≤1e-6 (quadrature-limited) |
| Unitarity | Norm of ψ(t) over time | TDSE solvers | ≤1e-12 split-operator / CN |
| Orthogonality, nodes | ⟨ψ_a\|ψ_b⟩ = 0; n−l−1 radial nodes | Eigenstates | ≤1e-6 |
| Statistical distribution | Speed histogram vs Maxwell–Boltzmann; χ² or KS test | MD and Monte Carlo at equilibrium | p-value not tiny; state the sample size |
| Determinism | Two runs bit-identical | All (use seeded RNG for stochastic models) | identical |
| Frame-rate independence | Advance with two different frame schedules | All time-dependent pages | identical |
| UI → physics | Change a parameter in the real UI; state must change and match an independent run | All interactive pages (`psim check`) | identical to Node run |

## Precision with purpose (stop here)

| Kind of check | Good enough |
|---|---|
| Identity of an evaluated formula (F = −dV/dr, analytic gradient vs finite difference) | 1e-8 relative |
| Conservation in a long symplectic run | 1e-8 (few-body), 1e-4 (many particles), bounded, no drift |
| Energy *balance* of a driven or damped system (ΔE = W_in − W_lost) | Normalise by the largest energy flow in the run (e.g. the work of one drive), not by a small stored energy; then require the expected convergence order under step refinement |
| Comparison with an exact solution at the method's resolution | the method's known error at that step size |
| Comparison with a measurement | its stated uncertainty (often 0.1–5 %) |
| Convergence order | within ±0.3 of theory |

Going further costs tokens without changing any claim. A **differential** measurement (subtract a reference run that shares the systematic error) beats brute-force precision when the effect is tiny (see the Mercury test in `tests/real-world.test.mjs`).

## Rules

- A check that never fails is not a check. Make sure each test would catch the bug it targets. p-simtor's own convergence check caught a bug where the step size was silently ignored.
- Put measured values in the page (Verification panel), not adjectives. Gravity reports measure the maximum conservation error over every computed step, and also retain the final error.
- A report applies only to its tested setup. Gravity pages withdraw the report and mark verification needed after solver settings or initial conditions change; exported JSON follows the same rule.
- At very fine timesteps, floating-point noise can dominate refinement ratios. Gravity measures order at an explicitly reported calibration step based on the initial dynamical scales (at most 2000 steps per crossing); conservation and the predictability horizon still use the actual finer production step.
- When a check fails, fix the model, method or step size. Don't loosen the threshold to make it pass. If a limit is real (chaos, approximation), report it as a limit.
- Separate the **model's** validity from the **numerics'** accuracy. A perfectly converged softened N-body run is still softened.

## Verification follows the scientific claim

Choose checks after defining the model and observables. Deterministic replay or a finite output is useful engineering evidence, but does not establish physical accuracy. Benchmark against a solution, dataset or independently justified limiting case outside the implementation being checked. Disclose whether a comparison uses measured data, an exact result of the same reduced model, or an empirical fit.

- **Fields/PDEs:** refine spatial resolution and time step independently. With spectral or explicit schemes the stable/accurate step shrinks with the grid (split-operator error ∝ k_max²Δt; explicit diffusion Δt ∝ Δx²): halving Δx at a fixed Δt can make results worse, so refine Δt along with Δx. Test boundary conditions, conservation or the stated loss/drive balance, and a linear mode or other independent limit. A homogeneous oscillator check alone does not verify spatial propagation.
- **Monte Carlo:** check equilibration and correlated-sample uncertainty, varying sample length and independent seeds. Study volume and discretisation for claims requiring those limits. Sweeps have no integrator convergence order and do not conserve energy per update; they sample the specified distribution.
- **Node and the browser are different engines:** `Math.exp`, `log`, `cos` and `hypot` may differ in the last bit, and a chaotic or Monte Carlo model turns that into a different sample within a few steps. The embedded report is computed in Node, so for such models verify statistical claims (means, distributions, rates) rather than one seed's exact trajectory, and do not claim the page shows the exact configuration the report measured unless a check proves it (e.g. `core/numerics.mjs` `rng` itself is portable).
- **Known, disclosed limits:** a check that compares with a better reference and finds the documented shortfall of the model (e.g. the plain Bethe formula 6 % long at a few MeV) is `{ pass: true, state: 'warn' }` with the deviation in its value, and the limit stated on the page. Never use `warn` to ship a failure of the implementation.
- **Importance sampling and restricted sampling** (aiming only at small impact parameters, oversampling a rare region): compare with the exact expectation of the restricted problem, and state the restriction and weights beside any count or rate.
- **Ratios, logarithms and fits of Monte Carlo averages:** resample the measured vectors together with `jackknife` or `bootstrap` (`core/data.mjs`); propagating independent errors ignores their correlation. Variational estimates (an effective mass at small separation, a trial wavefunction energy) are upper bounds: compare with a published fitted value only with that systematic stated.
- **Derived observables and estimators:** check the definition, units, normalisation and an independently known limit. Study binning, state truncation, quadrature and sample size separately from dynamics where applicable. Uncertainty must account for dependencies in the input samples. A property of an ensemble does not automatically hold for every sampled frame.
- **Reduced or fitted models:** verify the calculation within that model and state where it is applicable. A benchmark for one model, state space or estimator cannot validate a different theory or a broader physical claim.
- **Representation:** verify that visual values match their declared source quantity and mapping. Changing presentation alone must leave the model's state and observables unchanged. A colour, projection or symbolic label is not itself evidence for a physical prediction.

Use uncertainty for noisy observables, not a universal relative-error threshold. If a claimed signal is unresolved, report the missing evidence and withhold that claim rather than force its expected shape. Symbolic colour and illustrative geometry are explanation conventions, not physical verification results.

## Reporting format for composed models

`verify()` in a physics module returns:
```js
{ summary: 'One line on what was checked against what',
  checks: [{ label: 'Energy conserved in pure B (200 gyrations)', value: '|Δγ/γ| = 0', pass: true, state?: 'warn' }],
  notes: ['Neglected effects and why they are negligible'] }
```
`pass: false` makes `psim build --verify-module` exit with code 2, and the page shows the failed check in red.

Every check requires an explicit boolean `pass`, a label, and a finite measured value (number or formatted text). Missing results, NaN and infinity are rejected; imported report files undergo the same checks. Adaptive ODE tolerances are local error controls, not a guarantee on accumulated global error. `dopri5` throws on exhausted budgets or an unresolvable step; its `stopped` flag denotes only an intentional `onStep` stop.
