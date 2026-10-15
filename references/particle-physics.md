# Quarks, gluons and QCD

Status: a tested **pure SU(2) gauge module** (`core/lattice-gauge.mjs`) and a **QCD running-coupling module** (`core/qcd.mjs`), not a full QCD solver or ready-made recipe. Determine the question first: confinement, the mass gap, running coupling, hadron structure and quark–gluon plasma need different models. Read [Higgs guidance](higgs.md) for electroweak mass generation and [weak interaction](weak-interaction.md) for decays and the weak force. Particle data with uncertainties: `PDG` in `core/constants.mjs` (cite its `source` field).

QCD has an SU(3) colour gauge group and quark fields. The existing SU(2) module supports a related reduced gauge-theory demonstration; obtain agreement before substituting it for a request requiring actual SU(3) QCD or dynamical quarks. See [PDG, QCD](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-qcd.pdf) and [PDG, lattice QCD](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-lattice-qcd.pdf).

## 1. Pure SU(2) lattice gauge theory (tier `approximate`)

`core/lattice-gauge.mjs` samples pure SU(2) Yang–Mills theory (no dynamical quarks) on a periodic Euclidean lattice in 2, 3 or 4 dimensions with the Wilson action, using heat-bath plus over-relaxation updates. Name the gauge group, dimension, action and lattice size. Euclidean configurations and Monte Carlo sweeps are not a real-time movie of gluons.

```js
import { createLattice, sweep, wilsonLoops, polyakovLoop, plaquetteField, creutzRatio } from 'psim:core/lattice-gauge.mjs';
import { jackknife } from 'psim:core/data.mjs';
const lat = createLattice({ L: 6, d: 4, beta: 2.3, seed: 1 });
for (let i = 0; i < 50; i++) sweep(lat);                          // initial warm-up; diagnose equilibration
const samples = [];                                               // one vector of loops per measurement
for (let i = 0; i < 400; i++) { sweep(lat); samples.push(wilsonLoops(lat, { Rmax: 3, Tmax: 3 }).values.slice()); }
const W = m => (R, T) => m[(R - 1) * 3 + T - 1];
const chi22 = jackknife(samples, m => creutzRatio(W(m), 2, 2));    // { value, error }: correlated resampling
```

Measured cost in Node (heat bath + 2 over-relaxation; browsers are similar):

| Lattice | `sweep` | `wilsonLoops` (all 3×3) | `apeSmear` ×5 |
|---|---|---|---|
| 6⁴ | 6 ms | 3 ms | 5 ms |
| 8⁴ | 15 ms | 6 ms | 10 ms |
| 12³ | 4 ms | 1 ms | 2 ms |
| 16³ | 11 ms | 2 ms | 6 ms |

Use `wilsonLoops`, not repeated `wilsonLoop` calls (one loop at a time, link by link: ~20× slower). `polyakovLoop(lat)` gives the finite-volume diagnostic: |⟨P⟩| well above its noise means the box is too short to confine at that β. `plaquetteField(lat)` gives the local action density 1 − ½Tr U_P per site and plane for GPU views; `siteCoords`/`siteIndex` convert site numbers.

**What to show:**
- Averaged Wilson loops ⟨W(R,T)⟩ versus loop area and perimeter in a specified confined regime; show whether the data support an area law.
- Creutz ratios as finite-loop estimators of σa². Small loops on a small lattice do not automatically resolve an asymptotic string tension.
- The effective static potential aV_eff(R,T) = −ln[⟨W(R,T)⟩/⟨W(R,T−1)⟩], with its T dependence checked before claiming a ground-state potential. Do not take the logarithm of each individual loop and then average it.
- A computed plaquette-action slice with a labelled scalar colour bar. It shows an observable on a sampled SU(2) lattice configuration, not a literal coloured picture of quarks or gluons.
- The plaquette history, warm-up interval, number of measurements and statistical uncertainty. Use `meanWithError` for correlated scalar series; for ratios, logarithms and fits use `jackknife(samples, estimator)` (or `bootstrap`) from `core/data.mjs` on the measured loop vectors, so correlations are kept. A fixed warm-up count alone does not establish equilibration.

If loop averages are non-positive or the signal is unresolved within uncertainty, withhold the logarithmic potential/string-tension claim and state the sampling limitation. Never take absolute values or clamp the data just to obtain the expected curve.

**Required checks in `verify()`** (all in `tests/lattice-gauge.test.mjs` and `tests/lattice-measure.test.mjs`):
- 2D infinite-volume solution ⟨P⟩ = I₂(β)/I₁(β), ⟨W(R,T)⟩ = ⟨P⟩^{RT}, in a regime where the periodic finite-volume correction is negligible or accounted for;
- strong coupling in d dimensions, `strongCouplingPlaquette(β, d)` = u + 2(d − 2)u⁵ (u + 4u⁵ in 4D, u + 2u⁵ in 3D), at small β (≲ 0.5–1);
- gauge invariance of observables;
- over-relaxation preserves the action;
- Creutz ratios consistent with the expected confined regime, with enough sampling to resolve the claim.

Also study sampling length, seed dependence and finite volume where the claim needs them. The checks above verify the SU(2) model, not SU(3) QCD or measured proton properties.

Run one sweep (plus its measurement) per `step()`; the shell's frame budget then keeps the page responsive. Set `timeUnit: 'sweeps', timeDigits: 0`, a `speed` of what one frame can compute from the table above, and `speedRange` so the slider cannot request 100× that. Keep lattices at the sizes in the table; the shell has no Web Worker support. Set `testSteps` (e.g. 40) so `psim check`'s parameter-edit test does not run 600 slow sweeps, and check with `--advance` past the warm-up so the figures have data.

**Gauge-dependent links.** Without gauge fixing, each single link is exactly Haar-distributed (uniform on S³) at every β and ⟨½Tr U⟩ = 0 (Elitzur, PRD 12, 3978 (1975)): links do **not** cluster near 𝟙 as β grows; plaquettes and other closed loops do. Links cluster only in a fixed gauge: `landauGauge(lat)` (a Gribov copy, not unique; it changes `lat` in place to an exactly gauge-equivalent configuration and returns the accumulated g). `gaugeTransform(lat, g)` applies a given g (or pass a separate `rand` to leave the Markov chain untouched); `centreTransform(lat, { dir, slice })` flips the Polyakov loop but not the action (Z₂ centre symmetry, the order parameter's symmetry). Test each with its invariants.

**Finite temperature** (L³ × N_t with small N_t): the Polyakov loop along t is the deconfinement order parameter. For SU(2) with N_t = 4 the transition is at β_c = 2.2986(6) (Fingberg, Heller & Karsch, NPB 392, 493 (1993)). Space–time loops on such a lattice: `wilsonLoops(lat, { Rmax, Tmax, planes: [[0, 3], [1, 3], [2, 3]] })` (each plane is checked against its own extents).

**Topology and instantons (4D).** `topology(lat, { improved: true })` gives the gauge-invariant action density s(x) (in units of one instanton, 8π²/g²) and topological charge density q(x) from clover field strengths, with Q = Σq and S = Σs. Raw Monte Carlo samples are pure noise in q(x) (Q non-integer, S hundreds of instantons): smooth first with `sweep(lat, { cool: true })` and disclose the number of cooling sweeps. Cooling also shrinks instantons; lumps narrower than ~2.5a are unresolved (Q no longer near an integer), so report them as such. Benchmark: `bpstInstanton({ L, rho, anti })` builds an exact BPST (anti-)instanton; on 12⁴ with ρ = 3a the improved clover gives Q = 0.986 (plain clover 7–13 % low), and after a few cooling sweeps S/|Q| ≈ 1.04 (self-dual). Show q(x) with a signed volume (`createVolumeRenderer` `signed: true`). Spread an expensive sweep over frames with `sweep(lat, { sites: [s0, s1] })`. Sources: Belavin et al., PLB 59, 85 (1975); de Forcrand, García Pérez & Stamatescu, NPB 499, 409 (1997).

Views of 4D configurations and of SU(2) as S³: `core/projection.mjs` and [visual language](visual-language.md) (slices, projections, group space).

**Visible scope:** “Pure SU(2) gauge model; colour shows plaquette action; playback counts Monte Carlo sweeps.” Explain absent dynamical quarks (hence no quark-pair string breaking), finite spacing/volume and the need for continuum/volume extrapolation in the model panel. Physical-unit conversion requires a cited scale-setting procedure for the same model.

## 2. The mass gap: glueball masses (tier `approximate`)

"Prove the Yang–Mills mass gap" asks for an open Millennium Prize problem (Jaffe & Witten, Clay Mathematics Institute, 2000): no computation proves it (see [facts](facts.md)). What can be computed is the lightest glueball mass in a regularised theory, with evidence that the ratio to the string tension stays finite as the lattice spacing shrinks.

- **Model:** pure SU(2) in 2+1 dimensions (cheap, with a same-action benchmark: Teper, PRD 59, 014512 (1999)) or 3+1 dimensions (Lucini, Teper & Wenger, JHEP 06 (2004) 012). Name the group, dimension, β, volume.
- **Operators:** zero-momentum sums over each time slice of spatial Wilson loops built from **smeared** links: `apeSmear(lat, { alpha: 0.5, steps })` smooths the spatial links, then `wilsonLoops(smeared, { Rmax, planes: spatial planes })` per time slice. Unsmeared plaquettes have almost no overlap with the ground state; larger smeared loops (R = 2–4 at several smearing levels) do much better.
- **Correlators and the variational method:** C_ij(Δ) = ⟨O_i(t) O_j(t+Δ)⟩ − ⟨O_i⟩⟨O_j⟩ averaged over t; solve C(1)v = λC(0)v with `eigenGeneralized(C(1), C(0))` from `core/linalg.mjs`; the largest λ gives a·m = −ln λ. Errors: `jackknife` over binned correlators.
- **Bias:** an effective mass at Δ = 1 is a variational upper estimate; compare with a published fitted mass only with that systematic stated. Effective masses at larger Δ quickly lose signal: show them with error bars and withhold unresolved points.
- **Continuum:** one β shows a resolved gap on that lattice only (any finite lattice has a discrete spectrum). The evidence that matters is m/√σ approaching a finite value as a → 0, from several β; an imported continuum extrapolation must be labelled as imported.

## 3. Other QCD questions and analytical context

| Topic | What to compute | Source to cite |
|---|---|---|
| Static potential fit | Fit V(r) = −A/r + σr + C to resolved data. If interpreting A = C_F α_s ħc, use C_F = 3/4 for SU(2), 4/3 for SU(3); do not impose SU(3) numbers on the SU(2) run | Eichten et al., PRD 17, 3090 (1978); PDG for colour factors |
| Physical scale | Sommer definition r₀²F(r₀) = 1.65; determine r₀/a and state the model-specific calibration rather than borrowing physical QCD units | Sommer, NPB 411, 839 (1994) |
| String breaking | Use a calculation/dataset with dynamical quarks and its masses; the breaking distance is not universal | Bali et al., PRD 71, 114513 (2005) |
| Asymptotic freedom | `runAlphaS({ loops })` from `core/qcd.mjs`: MS-bar, 1–4 loops, flavour thresholds at the PDG quark masses, from `PDG.alphaS_MZ`; `betaCoefficients(nf, { group })` for other gauge groups. Shade Q ≲ 1 GeV (α_s ≳ 0.5, perturbation theory fails). Check: 1-loop closed form `alphaS1Loop`, exact 2-loop `exactTwoLoop`. The thresholds use continuous matching (no decoupling constants): say so | PDG QCD review; van Ritbergen, Vermaseren & Larin, PLB 400, 379 (1997) for β₃ |
| Proton structure | Parton distributions from published fit tables (cite set and version) | NNPDF / CT18 / MSHT |

A running-coupling formula is an analytical parameter scan, not confinement dynamics. A Cornell fit is phenomenological context. Published hadron/plasma data must be identified as imported data, including provenance and uncertainty; do not attribute it to the local SU(2) calculation.

## Rules

- Put an adjacent legend such as **“Red/green/blue label colour charge; they are not optical colours.”** Colour charge is a physical quantum property; the names and screen pigments are conventions. Gauge-invariant observables do not depend on this choice of colour basis. See [DOE, quarks and gluons](https://www.energy.gov/science/doe-explainsquarks-and-gluons) for the optical-colour distinction and PDG for the gauge description.
- Do not draw three coloured quarks and classical springs and call that a numerical QCD solution. A labelled pedagogical diagram can accompany the computed observables.
- A separation slider on a formula is a "parameter sweep", not time evolution; label it.
- Quote numbers with source and uncertainty.
- Put model identity and any important reduction beside the view; keep fuller limits in the model panel and final reply.
