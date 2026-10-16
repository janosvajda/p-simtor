# Cosmology (the expanding universe)

Status: verified module `core/cosmology.mjs` (tests in `tests/cosmology.test.mjs`). No ready-made recipe. The module covers the homogeneous FLRW background and linear growth. It does not cover perturbation spectra, the CMB power spectrum, recombination physics or structure formation beyond linear theory. Read [radiation and relativity](radiation-relativity.md) for black-body radiation and [stars](stars.md) for Saha ionisation.

"The universe expands into space" and "galaxies fly apart through space" are common misreadings. Space itself stretches, so comoving distances stay fixed and physical distances scale with a(t). Recession speeds above c are allowed (D > c/H₀). The edge of the observable universe is the particle horizon (≈ 14.2 Gpc comoving for Planck 2018), not c × age (13.8 Gly). Say this in the Fact check line when the request assumes otherwise.

## Candidate models (pick by the question)

| Question | Model | Observables and benchmark |
|---|---|---|
| How old is the universe? Why is the age not 1/H₀? | `age(p)` for ΛCDM and for comparison models | Planck 2018: 13.787 ± 0.020 Gyr (`PLANCK18_DERIVED.age`); EdS gives exactly 2/(3H₀), Milne exactly 1/H₀ |
| How does a(t) evolve? Will it recollapse? | `friedmannTrajectory` (dopri5 on the acceleration equation, started from the Friedmann constraint) | Closed matter-only universe: cycloid a = A(1 − cos η) with crunch at 2πB/H₀; de Sitter e^{H₀t}; EdS t^{2/3} |
| How far away is a galaxy at redshift z? | `comovingDistance`, `transverseComovingDistance`, `luminosityDistance`, `angularDiameterDistance`, `distanceModulus` (Hogg 1999 definitions) | astropy FlatLambdaCDM(70, 0.3) values to 1e-9; Mattig formula for ΩΛ = 0; Milne D_L = D_H z(1 + z/2) |
| Why does the angular size of galaxies reach a minimum? | `angularDiameterDistance(z)` peaks near z ≈ 1.6 (Planck 2018) | Etherington D_L = (1+z)² D_A |
| When did the expansion start to accelerate? | `zAcceleration`, `decelerationParameter` | (2ΩΛ/Ωm)^{1/3} − 1 without radiation; ≈ 0.64 for Planck 2018 |
| Radiation, matter and Λ eras | `zMatterRadiation`, `zMatterLambda`, `E(z)` term by term | z_eq = 3387 ± 21 (Planck); model 3403 (massive ν counted as matter, disclose) |
| Growth of structure | `growthFactor` (Heath integral), `growthRate`, `growthFactorCPT` (labelled fit) | EdS D = a; open ΩΛ = 0 closed form; f ≈ Ωm(a)^0.55 |
| Baryon acoustic oscillations | `soundHorizon(z, p)` (integral), `soundHorizonDragFit` (Aubourg et al. 2015) | r_drag = 147.21 ± 0.23 Mpc and r* = 144.57 Mpc (Planck) to 0.2 % |
| Cooling CMB | `cmbTemperature(z)` = T₀(1+z) | 2973 K at z* = 1089.8 |

## Method

- **Model:** `cosmology({ H0, Om, OL | flat: true, Or?, Tcmb, Neff, Obh2 })`. Curvature is derived as Ωk = 1 − Ωm − Ωr − ΩΛ, so a curved model is never inconsistent. Radiation defaults to photons at T_cmb plus Neff = 3.046 massless neutrinos. Pass `Or: 0` to compare with textbook dust models.
- **Preset:** `PLANCK18` (arXiv:1807.06209 Table 2, TT,TE,EE+lowE+lensing+BAO; H₀ = 67.66, Ωm = 0.3111, flat). Its ΩΛ is 0.68881 because radiation is included (Table 2 rounds to 0.6889).
- **Integrals** are taken in s = √a, where every integrand is smooth, with Gauss–Legendre on geometric panels. Ages and distances are accurate to about 1e-12. Don't reimplement them with a fixed grid in z.
- **Dynamics:** `friedmannTrajectory(p, { tMax, n, origin })` integrates ä = H₀²(−Ωr/a³ − Ωm/2a² + ΩΛa) through turnaround and stops at a Big Crunch (`crunch`). Use `origin: 'today'` for models without a Big Bang (de Sitter). The error falls with `rtol` (1e-6 in a at the default 1e-10, 1e-7 at 1e-12).
- **Units:** Gyr, Mpc, km s⁻¹ Mpc⁻¹, K. The constants are `MPC` (IAU 2015 B2), `GYR` (Julian) and `HUBBLE_TIME_UNIT` = 977.79 Gyr.

## Verification

- Exact limits in `verify()`: EdS age 2/(3H₀) and D_C = 2D_H(1 − 1/√(1+z)); Milne age 1/H₀; de Sitter lookback ln(1+z)/H₀; radiation-only age 1/(2H₀); flat ΛCDM age (2/(3√ΩΛ)) asinh√(ΩΛ/Ωm)/H₀.
- Curvature: the Mattig formula (ΩΛ = 0, open and closed), the closed cycloid for age and for a(t) including recollapse.
- Published values: Planck 2018 age (1σ), z_eq (1σ), r_drag and r* (0.2 %), and the astropy documentation distances.
- Consistency: lookback(z) = t₀ − t(z); q(z_acc) = 0; D_L = (1+z)²D_A; growth from the Heath integral equals the direct growth ODE (1e-7).
- **Circularity:** Planck's age and distances come from fitting this same model. Reproducing them checks the code, not ΛCDM. Independent tests are local H₀ (SH0ES 73.0 ± 1.0 differs from Planck by about 5σ: the "Hubble tension"; disclose it if H₀ is a slider), supernova distance moduli and BAO.

## Representation

- Draw comoving and physical coordinates as separate, labelled views. A comoving grid stays fixed while a(t) rescales the physical one.
- Use log axes for a(t) and z: the eras span 10 decades. Mark z_eq, z_acc and today.
- Show distances as several curves (D_C, D_L, D_A, light-travel distance c·lookback) on one plot. Users routinely mix them up.
- Animated "expansion" shows galaxies at fixed comoving positions scaled by a(t), with no motion through space. Photons are drawn with stretching wavelength (λ ∝ a).
- Sliders for Ωm and ΩΛ: show where the model has no Big Bang (functions throw "no Big Bang" or "E² ≤ 0") and say so on screen instead of drawing a curve.

**Visible scope:** "Homogeneous FLRW background with matter, radiation (massless neutrinos), Λ and curvature; linear growth without radiation; massive neutrinos and dark-energy evolution (w ≠ −1) not modelled."

## API

| Export | Signature → result | Conditions |
|---|---|---|
| `cosmology` | `({ H0, Om, OL \| flat, Or?, Tcmb = 2.7255, Neff = 3.046, Obh2?, name })` → frozen `{ H0, h, Om, Or, OL, Ok, Tcmb, Neff, Obh2, Ogamma }` | throws on H0 ≤ 0, Om < 0, missing OL without `flat` |
| `PLANCK18`, `PLANCK18_DERIVED` | preset model; published age, z_eq, z*, r*, z_drag, r_drag, Ωch² with 1σ errors | arXiv:1807.06209 Table 2 |
| `omegaGammaH2` | `(Tcmb)` → Ωγh² (2.4730e-5 at 2.7255 K) | Tcmb ≥ 0 |
| `E`, `hubbleParameter` | `(z, p)` → H/H₀; H in km/s/Mpc | throws where E² ≤ 0 or z ≤ −1 |
| `hubbleTime`, `hubbleDistance` | `(p)` → 1/H₀ in Gyr; c/H₀ in Mpc | |
| `age` | `(p, z = 0)` → Gyr since the Big Bang | throws without a Big Bang |
| `lookbackTime` | `(z, p)` → Gyr | any z the model reached |
| `comovingDistance` | `(z, p)` → line-of-sight D_C, Mpc | |
| `transverseComovingDistance` | `(z, p)` → D_M with sinh/sin curvature, Mpc | |
| `angularDiameterDistance`, `luminosityDistance` | `(z, p)` → D_A, D_L, Mpc | |
| `distanceModulus` | `(z, p)` → μ = 5 log₁₀(D_L/10 pc) | z > 0 |
| `particleHorizon` | `(p, z = 0)` → comoving horizon, Mpc | Infinity when Ωr = Ωm = 0 |
| `zMatterRadiation`, `zMatterLambda` | `(p)` → equality redshifts | need Ωr > 0; ΩΛ, Ωm > 0 |
| `decelerationParameter`, `zAcceleration` | `(z, p)` → q; `(p)` → z where q = 0 | throws if ΩΛ ≤ 0 |
| `cmbTemperature` | `(z, T0 = 2.7255)` → K | |
| `growthFactor`, `growthRate` | `(a, p)` → D (→ a early), f = dln D/dln a | radiation ignored; Ωm > 0 |
| `growthFactorCPT` | `(a, p)` → Carroll–Press–Turner fit | APPROXIMATE (~1 %) |
| `soundHorizon` | `(z, p)` → comoving r_s, Mpc (integral) | needs `Obh2` |
| `soundHorizonDragFit` | `({ obh2, ocbh2, onuh2 })` → r_drag, Mpc | fitting formula, 0.02 % near Planck |
| `friedmannTrajectory` | `(p, { tMax (Gyr), n = 200, origin = 'bang' \| 'today', aStart = 1e-4, rtol = 1e-10 })` → `{ t, a, H, crunch }` | stops at a Big Crunch |
| `MPC`, `GYR`, `C_KMS`, `HUBBLE_TIME_UNIT`, `T_CMB`, `NEUTRINO_FACTOR` | constants | |
