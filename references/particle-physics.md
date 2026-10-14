# Quarks, gluons and QCD

Status: a tested **pure SU(2) gauge module**, not a full QCD solver or ready-made recipe. Determine the question first: confinement, running coupling, hadron structure and quark–gluon plasma need different models. Read [Higgs guidance](higgs.md) for electroweak mass generation.

QCD has an SU(3) colour gauge group and quark fields. The existing SU(2) module supports a related reduced gauge-theory demonstration; obtain agreement before substituting it for a request requiring actual SU(3) QCD or dynamical quarks. See [PDG, QCD](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-qcd.pdf) and [PDG, lattice QCD](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-lattice-qcd.pdf).

## 1. Pure SU(2) lattice gauge theory (tier `approximate`)

`core/lattice-gauge.mjs` samples pure SU(2) Yang–Mills theory (no dynamical quarks) on a periodic Euclidean lattice in 2, 3 or 4 dimensions with the Wilson action, using heat-bath plus over-relaxation updates. Name the gauge group, dimension, action and lattice size. Euclidean configurations and Monte Carlo sweeps are not a real-time movie of gluons.

```js
import { createLattice, sweep, plaquette, wilsonLoop, creutzRatio } from 'psim:core/lattice-gauge.mjs';
const lat = createLattice({ L: 6, d: 4, beta: 2.3, seed: 1 });   // 6⁴: ~15 ms per sweep in Node
for (let i = 0; i < 50; i++) sweep(lat);                          // initial warm-up; diagnose equilibration
// measure W(R,T) over many sweeps; σa² ≈ χ(R,T) = −ln[W(R,T)W(R−1,T−1)/(W(R,T−1)W(R−1,T))]
```

**What to show:**
- Averaged Wilson loops ⟨W(R,T)⟩ versus loop area and perimeter in a specified confined regime; show whether the data support an area law.
- Creutz ratios as finite-loop estimators of σa². Small loops on a small lattice do not automatically resolve an asymptotic string tension.
- The effective static potential aV_eff(R,T) = −ln[⟨W(R,T)⟩/⟨W(R,T−1)⟩], with its T dependence checked before claiming a ground-state potential. Do not take the logarithm of each individual loop and then average it.
- A computed plaquette-action slice with a labelled scalar colour bar. It shows an observable on a sampled SU(2) lattice configuration, not a literal coloured picture of quarks or gluons.
- The plaquette history, warm-up interval, number of measurements and statistical uncertainty. Use `meanWithError` for correlated scalar series; for ratios/fits, block and resample the measured loop vectors together to retain correlations. A fixed warm-up count alone does not establish equilibration.

If loop averages are non-positive or the signal is unresolved within uncertainty, withhold the logarithmic potential/string-tension claim and state the sampling limitation. Never take absolute values or clamp the data just to obtain the expected curve.

**Required checks in `verify()`** (all in `tests/lattice-gauge.test.mjs`):
- 2D infinite-volume solution ⟨P⟩ = I₂(β)/I₁(β), ⟨W(R,T)⟩ = ⟨P⟩^{RT}, in a regime where the periodic finite-volume correction is negligible or accounted for;
- 4D strong coupling ⟨P⟩ ≈ u + 4u⁵;
- gauge invariance of observables;
- over-relaxation preserves the action;
- Creutz ratios consistent with the expected confined regime, with enough sampling to resolve the claim.

Also study sampling length, seed dependence and finite volume where the claim needs them. The checks above verify the SU(2) model, not SU(3) QCD or measured proton properties.

Run the Monte Carlo in chunks per frame (one sweep per frame for 6⁴) so the page stays responsive. Use a Web Worker if you go larger.

**Visible scope:** “Pure SU(2) gauge model; colour shows plaquette action; playback counts Monte Carlo sweeps.” Explain absent dynamical quarks (hence no quark-pair string breaking), finite spacing/volume and the need for continuum/volume extrapolation in the model panel. Physical-unit conversion requires a cited scale-setting procedure for the same model.

## 2. Other QCD questions and analytical context

| Topic | What to compute | Source to cite |
|---|---|---|
| Static potential fit | Fit V(r) = −A/r + σr + C to resolved data. If interpreting A = C_F α_s ħc, use C_F = 3/4 for SU(2), 4/3 for SU(3); do not impose SU(3) numbers on the SU(2) run | Eichten et al., PRD 17, 3090 (1978); PDG for colour factors |
| Physical scale | Sommer definition r₀²F(r₀) = 1.65; determine r₀/a and state the model-specific calibration rather than borrowing physical QCD units | Sommer, NPB 411, 839 (1994) |
| String breaking | Use a calculation/dataset with dynamical quarks and its masses; the breaking distance is not universal | Bali et al., PRD 71, 114513 (2005) |
| Asymptotic freedom | SU(3) one-loop α_s(Q²) = 12π/((33−2n_f) ln(Q²/Λ²)) within its perturbative regime; specify n_f, scheme, thresholds and input scale | PDG QCD review; quote its version and uncertainties for numerical inputs |
| Proton structure | Parton distributions from published fit tables (cite set and version) | NNPDF / CT18 / MSHT |

A running-coupling formula is an analytical parameter scan, not confinement dynamics. A Cornell fit is phenomenological context. Published hadron/plasma data must be identified as imported data, including provenance and uncertainty; do not attribute it to the local SU(2) calculation.

## Rules

- Put an adjacent legend such as **“Red/green/blue label colour charge; they are not optical colours.”** Colour charge is a physical quantum property; the names and screen pigments are conventions. Gauge-invariant observables do not depend on this choice of colour basis. See [DOE, quarks and gluons](https://www.energy.gov/science/doe-explainsquarks-and-gluons) for the optical-colour distinction and PDG for the gauge description.
- Do not draw three coloured quarks and classical springs and call that a numerical QCD solution. A labelled pedagogical diagram can accompany the computed observables.
- A separation slider on a formula is a "parameter sweep", not time evolution; label it.
- Quote numbers with source and uncertainty.
- Put model identity and any important reduction beside the view; keep fuller limits in the model panel and final reply.
