# Stars and stellar remnants

Status: verified module `core/stars.mjs` (tests in `tests/stars.test.mjs`). No ready-made recipe. It covers polytropes (Lane–Emden), zero-temperature white dwarfs with the exact degenerate-electron equation of state, the nuclear binding curve, Saha ionisation, power-law energy generation and global stellar scalings. It is not a stellar-evolution code: no opacity tables, convection or composition change. Read [radiation and relativity](radiation-relativity.md) for black-body spectra and colours and [cosmology](cosmology.md) for recombination redshifts.

Common misreadings: the Sun does not "burn" chemically. It is not on fire, and its centre is not 5772 K: that is the surface temperature, and the centre is ≈ 1.57×10⁷ K. Fusion in the core proceeds at about 10⁻³ W/kg, far less power per kilogram than a human body (≈ 1 W/kg). The Sun shines for billions of years because it is enormous. A white dwarf is held up by electron degeneracy pressure, not by heat, and a heavier white dwarf is *smaller*. Say this in the Fact check line when the request assumes otherwise.

## Candidate models (pick by the question)

| Question | Model | Observables and benchmark |
|---|---|---|
| What is inside a star? Why is the centre so hot? | `polytrope({ n: 3, M, R, mu })` (Eddington's standard model) with `profile(r)` | ρc/ρ̄ = 54.18; Pc = 11.05 GM²/R⁴. Sun: Tc ≈ 1.2×10⁷ K, ρc ≈ 76 g cm⁻³ vs the Standard Solar Model's 1.57×10⁷ K and 152 g cm⁻³ (`SSM_CENTRE`): APPROXIMATE, label it |
| Lane–Emden equation itself | `laneEmden(n)` | Exact n = 0, 1, 5; tabulated ξ₁, −ξ₁²θ′(ξ₁) for n = 1.5, 2, 3, 4 (Chandrasekhar 1939) |
| Why is there a maximum white-dwarf mass? | `whiteDwarf({ xc })` along a sequence, `chandrasekharMass(μe)` | M → M_Ch = 1.4563 M☉ (μe = 2, mass unit u) as R → 0; Sirius B 1.018 M☉ → 0.0080 R☉ (observed 0.0081, Bond et al. 2017) |
| Why is iron the end of fusion? | `semfBinding(Z, A)` vs `AME2020_BA` | B/A peaks at Ni-62 (8.7946 MeV) > Fe-58 > Fe-56; SEMF within 1 % for A ≥ 40, −19 % for He-4 |
| Valley of stability | `mostStableZ(A)` | within 1.5 of the stable isobar (Fe-56, Pb-208, U-238) |
| Spectral lines vs temperature; when does hydrogen ionise? | `sahaHydrogen(T, nH)`, `halfIonisationTemperature(nH)`, `sahaRatio` | 50 % at ≈ 9550 K for P_e = 20 Pa (Carroll & Ostlie §8.1: 9600 K) |
| pp chain vs CNO cycle | `epsPP`, `epsCNO` (power laws T⁴ and T^19.9) | crossover ≈ 1.8×10⁷ K for X = 0.7, X_CNO = 0.01: APPROXIMATE, label it |
| HR diagram, luminosity, temperature | `stefanBoltzmannLuminosity`, `effectiveTemperature`, `mainSequenceLuminosity`, `mainSequenceRadius`, `mainSequenceLifetime` | IAU 2015 B3 L☉, R☉ give T_eff = 5772 K (exact check); the main-sequence relations are APPROXIMATE (factor-2 scatter) |
| How long can a star shine? | `kelvinHelmholtzTime`, `nuclearTime` | Sun: 31.4 Myr (contraction, Kelvin's argument for a young Earth) vs ≈ 10 Gyr (fusion) |
| Maximum luminosity, massive stars | `eddingtonLuminosity(M, { X })` | 1.257×10³¹ W per M☉ for pure hydrogen = 3.28×10⁴ L☉ |

## Method

- **Lane–Emden** θ″ + 2θ′/ξ + θⁿ = 0 is integrated with dopri5 (rtol 1e-12) from a series start. The first zero ξ₁ is located by event detection on dense output. `theta(ξ)` and `dtheta(ξ)` interpolate the stored steps (cubic Hermite, 1e-9) and return NaN outside the star. Results are cached per n.
- **Polytropes** scale the dimensionless solution to (M, R). The temperature uses the ideal gas with mean molecular weight μ (0.62 for solar composition when fully ionised). Radiation pressure is neglected (below 0.1 % of the pressure at the solar centre).
- **White dwarfs:** P(x) = (πm_e⁴c⁵/3h³)[x(2x² − 3)√(1+x²) + 3 asinh x], ρ = μe m_u (8π/3)(m_e c/h)³ x³. The integration uses w = √(1+x²) − 1, which keeps full precision in the non-relativistic limit. A whole mass–radius curve (xc = 10⁻³ … 10⁴) takes a few ms. `whiteDwarfOfMass` inverts it.
- **Chandrasekhar mass convention:** with m = u (μe = 2 exactly for He, C and O) M_Ch = 5.825/μe² M☉ = 1.4563 M☉; with m = m_p one gets 1.4353 M☉. Older texts round to 1.44. State the convention. The effective limit of real white dwarfs is a few per cent lower because of Coulomb corrections, inverse β decay and general relativity, which are not modelled.
- **SEMF:** Rohlf coefficients (aV 15.75, aS 17.8, aC 0.711, aA 23.7, aP 11.18 MeV). Swap in another set as the third argument.
- **Units:** SI throughout. `R_SUN`, `L_SUN`, `TEFF_SUN` are IAU 2015 B3 nominal values; `M_SUN` = GM☉/G.

## Verification

- Exact: Lane–Emden n = 0, 1, 5 profiles and constants; n = 1 polytrope density sin(πr/R)/(πr/R); the white-dwarf low-density limit equals the n = 3/2 polytrope (2e-6); the high-density limit equals M_Ch (1e-6).
- Measured: the mass integral of a profile returns M (1e-8); finite-difference hydrostatic equilibrium dP/dr = −Gmρ/r² (1e-6); convergence of ξ₁ and ω under rtol refinement.
- Published: Chandrasekhar's Lane–Emden table; the AME2020 binding energies (consistent with the `NUCLIDE_MASS` atomic masses to 0.01 keV for H-2, C-12, O-16); σ_T (CODATA); Sirius B; T_eff = 5772 K from IAU nominal values.
- Label every approximate model on screen: polytrope Sun, power-law ε, main-sequence fits, SEMF for light nuclei.

## Representation

- Interior profiles: plot ρ, P, T and m(r) normalised to central values against r/R, on log axes for ρ and P. Mark the region where fusion happens (ε ∝ ρT⁴ is concentrated inside ≈ 0.2 R☉).
- Mass–radius of white dwarfs: plot R against M with the M_Ch asymptote. Mark Sirius B and 40 Eri B as observed points with error bars.
- Binding curve: B/A against A with the SEMF curve and the AME points, the peak at Ni-62 marked, and arrows for "fusion releases energy" (left) and "fission releases energy" (right).
- HR diagram: log L against log T_eff with temperature increasing to the *left*, which is the astronomical convention; say so on the axis.

**Visible scope:** "Polytropic and zero-temperature models; no stellar evolution, opacity or convection; nuclear rates are power-law fits near 15 MK; white dwarfs without Coulomb or GR corrections."

## API

| Export | Signature → result | Conditions |
|---|---|---|
| `laneEmden` | `(n, { rtol = 1e-12, xiMax })` → `{ n, xi1, dtheta1, omega, densityRatio, Wn, xiEnd, theta(ξ), dtheta(ξ) }` | 0 ≤ n ≤ 10; n ≥ 5 has no zero (xi1 = Infinity, values at xiMax = 50) |
| `polytrope` | `({ n, M, R, mu = 0.62 })` → `{ rhoMean, rhoC, Pc, Tc, K, alpha, profile(r) → { rho, P, T, m } }` | 0 ≤ n < 5; ideal gas, no radiation pressure |
| `polytropeMass` | `(star, n = 64)` → ∫4πr²ρ dr, kg | consistency check |
| `chandrasekharMass` | `(mue = 2, { mUnit = u })` → kg | |
| `whiteDwarf` | `({ xc } \| { rhoC }, mue = 2, rtol)` → `{ mass, radius, xc, rhoC, mue }` | 0 < xc < 1e6; T = 0, no Coulomb/GR |
| `whiteDwarfOfMass` | `(M, { mue })` → same as `whiteDwarf` | 0 < M < 0.999 M_Ch |
| `semfBinding` | `(Z, A, k = SEMF_ROHLF)` → MeV | integers, 0 ≤ Z ≤ A |
| `mostStableZ` | `(A, k)` → continuous Z | |
| `bindingFromMass` | `(Z, A, massU)` → MeV from the atomic mass | |
| `AME2020_BA`, `SEMF_ROHLF` | data: `{ Z, A, BA (keV) }`; coefficient set | Wang et al. 2021; Rohlf 1994 |
| `sahaRatio` | `(T, ne, { chi = 13.598 eV, uRatio = 0.5 })` → n_{i+1}/n_i | LTE, no pressure ionisation |
| `sahaHydrogen` | `(T, nH)` → `{ x, ne }` | pure hydrogen |
| `halfIonisationTemperature` | `(nH)` → K where x = 0.5 | 0 < nH < 1e30 m⁻³ |
| `epsPP`, `epsCNO` | `(rho, X, T)`, `(rho, X, XCNO, T)` → W/kg | APPROXIMATE power laws near T₆ ≈ 15 |
| `stefanBoltzmannLuminosity`, `effectiveTemperature` | `(R, T)` → W; `(L, R)` → K | positive inputs |
| `eddingtonLuminosity` | `(M, { X = 1, kappa })` → W | electron scattering by default |
| `kelvinHelmholtzTime`, `nuclearTime` | `(M, R, L, { q = 1 })`, `(M, L, { efficiency = 0.007, coreFraction = 0.1 })` → s | |
| `mainSequenceLuminosity`, `mainSequenceRadius`, `mainSequenceLifetime` | `(M/M☉)` → L/L☉, R/R☉, Gyr | APPROXIMATE, 0.08–300 M☉ |
| `R_SUN`, `L_SUN`, `TEFF_SUN`, `M_SUN`, `SIGMA_T`, `CHI_H`, `NEUTRON_MASS_U`, `U_MEV`, `YEAR`, `SSM_CENTRE` | constants and reference values | sources in the module header |
