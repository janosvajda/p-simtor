# Stiff ODEs and special functions

Two modules extend `core/ode.mjs` and `core/numerics.mjs`. Both are plain ES modules that run in Node and the browser and depend only on other `core/` modules:

- `core/stiff.mjs` provides implicit integrators for stiff systems M y′ = f(t, y), a Jacobian-spectrum diagnostic and an explicit-or-implicit chooser.
- `core/special.mjs` provides special functions of mathematical physics that `numerics.mjs` does not already have.

`numerics.mjs` keeps Γ, ln Γ, erf, P/Q (incomplete Γ), J_n, I_n, K(m), L_n^(α), P_l^m and Y_lm. Use them from there; `special.mjs` does not duplicate them.

## Stiff integration: when and which

A system is **stiff** when its Jacobian has decaying modes much faster than the motion you want to follow. Examples are chemical kinetics with rates spread over many decades, relaxation oscillators, heat conduction after spatial discretisation, and circuits with tiny time constants. In that case an explicit method such as `odeStepper` with `dopri5` is limited by stability (h·|λ_max| ≲ 3.3), not by accuracy. It then takes millions of tiny steps, or it oscillates and blows up.

Symptoms: `dopri5` needs far more steps than the smoothness of the solution suggests, its step size barely changes after a transient, or its rejection count is high.

1. **Diagnose.**
   - `stiffnessRatio({ f, t, y })` gives the Jacobian eigenvalues at one state. Stiffness is local: Robertson kinetics is not stiff at t = 0 but becomes stiff after about 1e-3 s.
   - `odeAuto({ f, y0, t0, tEnd })` samples the spectrum at t0 and at the end of a short explicit trial run. It decides on the stiffness index S = max(−Re λ)·(tEnd − t0) > 1e3 and returns a ready stepper.
2. **Choose a method** (all of them adapt h; for very large n the dense LU costs O(n³) per factorisation):

| method | order | stability | use when |
|---|---|---|---|
| `'radau5'` (default) | 5 (dense output 3) | L-stable, also for eigenvalues near the imaginary axis | Tight tolerances, oscillatory stiff problems, DAEs with a mass matrix (index 1). The most robust choice |
| `'rosenbrock'` (Rodas3) | 3 (embedded 2) | L-stable, stiffly accurate | Loose tolerances (rtol ≳ 1e-6). No Newton iteration, so it is cheap per step and predictable. It needs an accurate Jacobian every step: pass `jac` |
| `'bdf'` | variable 1–5 | BDF1–2 A-stable; BDF3–5 only A(α) | Large, smooth, strongly damped problems; fewest Jacobians. Avoid when there are eigenvalues with large imaginary parts (lightly damped oscillations): use radau5 |

3. **Verify**, as for any solver:
   - Run at two tolerances (rtol and rtol/100) and compare with an independent benchmark. Their difference is a convergence estimate, not a certified error bound. Local tolerances do not guarantee the same bound on the accumulated global error.
   - Compare two methods. Robertson at t = 4e5 agrees between all three methods to 5e-9.
   - Check invariants. Linear invariants such as the Robertson total mass are kept to round-off by all three methods when the Jacobian is exact.
   - Check the equilibria, or an asymptotic result such as the Van der Pol period.
   - `adaptive: false, h` gives fixed steps for `refinementOrder` studies.

**Conditions and failure handling.**
- Use a forward integration direction. The mass matrix must be constant, and only radau5 accepts one. A singular M requires consistent initial values.
- Without `jac`, forward differences are used. This costs n extra f-evaluations per Jacobian; accuracy depends on state scales, smoothness and cancellation. At tight tolerances give an analytic Jacobian, especially for Rodas3. Derivative and Jacobian callbacks may reuse output work arrays.
- The stepper throws on non-finite states or derivatives, step-size underflow, a singular iteration matrix, or 60 consecutive rejections. A fixed-step Newton failure also throws.
- Events follow `odeStepper` semantics: sign changes of g are located on the dense output with Brent's method, and terminal events stop exactly at the event. A sign change that reverses within one step is missed, so cap `s.hmax` near thin regions.

## Special functions: which tool

| Need | Tool | Notes |
|---|---|---|
| Quantum tunnelling and turning points, WKB connection, Stark effect, rainbow caustics, the triangular well | `airy(x)` → `{ ai, aip, bi, bip }`, `airyAiZeros(k)` | Linear potential V = Fx: ψ ∝ Ai((x − E/F)/ℓ); triangular-well levels E_s = −a_s·(ħ²F²/2m)^{1/3} |
| Spherical waves, partial-wave scattering, the spherical box (zeros of j_l) | `sphericalBesselJ`, `sphericalBesselY` | Small kr is accurate (series); Σ(2l+1)j_l² = 1 checks truncation of a partial-wave sum |
| Cylindrical exterior problems, waveguides | `besselY(n, x)`; J_n is in numerics | Wronskian J_{n+1}Y_n − J_nY_{n+1} = 2/(πx) |
| Screened/Yukawa fields in 2D, heat loss from cylinders, the Debye-shielded line charge | `besselK(n, x)`; I_n is in numerics | Exponentially scaled recurrence retains high-order values even when K₀ underflows |
| Harmonic oscillator states of any n (coherent states, the classical limit) | `oscillatorWavefunction(n, ξ)`; `hermiteH(n, x)` for the polynomial | ψ in physical units: (mω/ħ)^{1/4} ψ_n(x√(mω/ħ)); no overflow for n up to thousands |
| Closed forms that reduce to hypergeometrics (Legendre functions, the Pöschl–Teller and Morse potentials, the Coulomb continuum, Rabi/Landau–Zener) | `hyp2f1`, `hyp1f1` | Real arguments only; throws outside the supported region rather than guessing |
| Pendulum at large amplitude, nonlinear oscillators, cnoidal waves, elliptic orbits in GR | `jacobiElliptic(u, m)`, `ellipticE(m)`, ellipticK (numerics) | Parameter m = k² throughout. Pendulum: sin(θ/2) = k·sn(K − ω₀t \| k²), period 4K(k²)/ω₀ |
| Statistics (Student-t p-values, binomial and F tails), Beta priors | `incompleteBeta`, `studentTCdf`, `studentTPValue`, `beta`, `logBeta` | χ² uses gammaQ from numerics |
| Γ-function derivatives, entropy of Dirichlet/Beta distributions, harmonic sums | `digamma` | ψ(n+1) = H_n − γ |
| Bose–Einstein and Fermi–Dirac integrals, Debye/Stefan constants, Casimir sums | `zeta`, `polylog` | g_s(z) = Li_s(z); BEC at nλ³ = ζ(3/2); f_s(z) = −Li_s(−z) (z ≤ 1, or integer s for z < −1) |
| Radiative transfer (E_n slabs), the Theis well function, line-source heat, Coulomb logarithms | `expintEn`, `expintE1`, `expintEi` | E_1 = Theis W(u) |
| Coupling of angular momenta, selection rules and line strengths, recoupling | `wigner3j`, `clebschGordan`, `wigner6j` | Exact BigInt arithmetic; `{ exact: true }` returns the rational square |
| Plasma dispersion (Re Z), Voigt profiles, the ion-acoustic response | `dawson(x)` | F(x) = (√π/2)e^{−x²}erfi x |

**Verification advice.** Check special-function values with an identity that uses a different method, not with the same formula:
- Wronskians: Ai/Bi, j/y, J/Y, I/K.
- The ODE the function satisfies, e.g. Ai″ = xAi or F′ = 1 − 2xF.
- An integral representation evaluated by quadrature: K_ν, Rayleigh's j_n, E_n, F(φ|m).
- Orthogonality: Hermite functions, 3j, 6j, Clebsch–Gordan.
- A published table: A&S, DLMF.

Oscillatory functions are only accurate in the *absolute* sense near their zeros, so compare there with tolerances scaled by the envelope.

## API: `core/stiff.mjs`

| Export | Result | Accuracy and conditions |
|---|---|---|
| `stiffStepper({ f, jac?, y0, t0 = 0, method = 'radau5', rtol = 1e-6, atol = 1e-9, h?, hmax = ∞, events = [], mass?, adaptive = true, maxOrder = 5 })` | `{ t, y, hmax, method, step(hLimit?), advanceTo(t, maxSteps, { exact }), sampleAt(t, out?), steps, rejected, events, stopped, stats: { fevals, jacobians, decompositions, solves, newton, order } }`, the same interface as odeStepper | Local error per step ≈ rtol·\|y\| + atol (atol may be an array). `jac(t, y)` returns n rows. The mass matrix is constant (radau5 only). `adaptive: false` uses a fixed h. Throws on invalid input, non-finite states, underflow, a singular iteration matrix, or 60 rejections. Measured orders: radau5 5.0, rosenbrock 3.0, BDF1 1.0, BDF2 2.0 |
| `stiffnessRatio(J \| { f, t, y, jac? }, { zeroTolerance = 1e-9 })` | `{ eigenvalues: [{ re, im }], spectralRadius, fastest, slowest, ratio }` | Eigenvalues by balanced Hessenberg QR (~1e-12·‖J‖). Neutral modes (\|Re λ\| ≤ zeroTolerance·ρ) are excluded from slowest and ratio. The result is local to one state |
| `odeAuto({ f, jac?, y0, t0 = 0, tEnd, rtol, atol, probes = [], trialSteps = 200, threshold = 1e3, events? })` | `{ stiff, method, stiffnessIndex, fastest, ratio, reason, trialSpan, stepper }` | A heuristic: S = max(−Re λ)·(tEnd − t0) over t0, the probes and the end of a dopri5 trial run. It returns a radau5 stiffStepper or a dopri5 odeStepper |

## API: `core/special.mjs`

| Export | Result | Accuracy and conditions |
|---|---|---|
| `airy(x)` | `{ ai, aip, bi, bip }` | Any finite x. Maclaurin series (\|x\| ≤ 1), asymptotic expansions (\|x\| ≥ 10) and Taylor integration of y″ = xy in between. About 1e-14 relative for x ≥ 0; for x < 0, about 1e-14 relative to the envelope \|x\|^{∓1/4}/√π. For x < −10 the phase error is about eps·(2/3)\|x\|^{3/2}. Ai underflows and Bi overflows past x ≈ 104 |
| `airyAiZeros(count)` | `[a_1, a_2, …]` (negative) | Asymptotic estimate plus Newton, about 1e-15. Checked against A&S Table 10.13 |
| `sphericalBesselJ(n, x)` | j_n(x) | Integer n ≥ 0, real x. Series, upward or Miller recurrence; ~1e-15 relative away from zeros |
| `sphericalBesselY(n, x)` | y_n(x) | Upward recurrence; −Infinity at x = 0; overflows for x ≪ n |
| `besselY(n, x)` | Y_n(x) | Integer n, x ≥ 0. Neumann series (x < 25) or Hankel (x ≥ 25), then recurrence; ~1e-15 relative to the envelope. Throws for x < 0 |
| `besselK(n, x)` | K_n(x) | Integer n, x ≥ 0. Series (x ≤ 2) or Steed's continued fraction, then scaled recurrence; error grows with order, approximately n·eps away from under/overflow |
| `hermiteH(n, x)` | Physicists' H_n(x) | Recurrence; overflows for large n·\|x\| |
| `oscillatorWavefunction(n, ξ)` | Normalised ψ_n(ξ) | Normalised recurrence with running exponent; any n (tested to 2000); ~√n·eps of max\|ψ\| |
| `hyp2f1(a, b, c, x)` | ₂F₁ | Real parameters, x ≤ 1. Series, the 1 − x connection and Pfaff transformations; log-scaled Gauss sum at x = 1. Terminating polynomials of degree ≤ 2048 use exact rational summation for the input doubles, rounded once; larger degrees throw. Nonterminating paths check cancellation and may reject unsupported regions |
| `hyp1f1(a, b, x)` | ₁F₁ = M(a, b, x) | Real arguments. Exact rational polynomial summation through degree 2048; otherwise series, Kummer transformation or large-\|x\| asymptotics with cancellation checks. A denominator pole throws unless the polynomial terminates before it |
| `ellipticE(m)` | E(m) | m ≤ 1 (m = k²), AGM; ~1e-15 |
| `jacobiElliptic(u, m)` | `{ sn, cn, dn, am }` | Any real u, m (A&S 16.10–11 transformations). am is continuous for m < 1. About 1e-15 absolute for \|u\| ≲ 1e3·K |
| `digamma(x)` | ψ(x) | x ≠ 0, −1, …; ~1e-15 absolute·(1 + \|ln x\|) |
| `beta(a, b)`, `logBeta(a, b)` | B, ln B | a, b > 0; Stirling differences avoid cancellation for a ≫ b |
| `incompleteBeta(x, a, b)` | I_x(a, b) | 0 ≤ x ≤ 1. Continued fraction; ~1e-14 for a, b ≲ 100, degrading as eps·min(a, b) when both are large |
| `studentTCdf(t, ν)`, `studentTPValue(t, ν, tails = 2)` | P(T ≤ t), the p-value | ν > 0 (real). Tail p-values keep relative accuracy for ν ≲ 1e6 (absolute ~1e-10 at ν = 1e8) |
| `zeta(s)` | ζ(s) | Real s ≠ 1. Borwein's η algorithm for s > 0 and the functional equation for s < 0; ~1e-15 |
| `polylog(s, z)` | Li_s(z) | Real s, z ≤ 1 (z < −1 only for integer s); ~1e-14. Throws for z > 1, divergent z = 1, and orders within about 1e-4 of an integer that cancel |
| `expintEn(n, x)`, `expintE1(x)`, `expintEi(x)` | E_n, E_1, Ei | x ≥ 0 for E_n (Ei for any x ≠ 0); ~1e-15; Ei overflows past x ≈ 710 |
| `wigner3j(j1, j2, j3, m1, m2, m3, { exact })` | Number, or `{ sign, num, den }` with value sign·√(num/den) | Integer or half-integer arguments; exact BigInt Racah sum; selection rules give 0; practical to j ~ a few hundred |
| `clebschGordan(j1, m1, j2, m2, J, M, { exact })` | ⟨j1 m1; j2 m2 \| J M⟩ | Condon–Shortley phase |
| `wigner6j(j1, …, j6, { exact })` | {j1 j2 j3; j4 j5 j6} | Exact Racah formula; 0 unless all four triads close |
| `dawson(x)` | F(x) | Series, Rybicki sampling (h = 0.2) or asymptotic; ~1e-15 |

**Tests.** `tests/stiff.test.mjs` and `tests/special.test.mjs` check all of the above against:
- Hairer & Wanner reference solutions: Robertson at t = 40 and 1e11, HIRES, and Oregonator.
- The Van der Pol asymptotic period.
- L- and A-stability on y′ = λy.
- Measured convergence orders.
- An exact DAE solution.
- A&S and DLMF tables.
- Wronskians, quadrature representations, orthogonality relations and closed forms.
