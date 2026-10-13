# Verification

Every page carries evidence for its trust tier. Choose checks that match the mathematics: there is no universal test.

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

## Reporting format for composed models

`verify()` in a physics module returns:
```js
{ summary: 'One line on what was checked against what',
  checks: [{ label: 'Energy conserved in pure B (200 gyrations)', value: '|Δγ/γ| = 0', pass: true, state?: 'warn' }],
  notes: ['Neglected effects and why they are negligible'] }
```
`pass: false` makes `psim build --verify-module` exit with code 2, and the page shows the failed check in red.

Every check requires an explicit boolean `pass`, a label, and a finite measured value (number or formatted text). Missing results, NaN and infinity are rejected; imported report files undergo the same checks. Adaptive ODE tolerances are local error controls, not a guarantee on accumulated global error. `dopri5` throws on exhausted budgets or an unresolvable step; its `stopped` flag denotes only an intentional `onStep` stop.
