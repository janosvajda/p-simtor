# Quarks, gluons and QCD

**Default: compute it.** Confinement can be computed from first principles on a lattice in the browser. Build that first, then add analytical context next to it.

## 1. Lattice gauge theory (tier `approximate`; the real computation)

`core/lattice-gauge.mjs` simulates pure SU(2) Yang–Mills theory (gluons without quarks) on a periodic lattice in 2, 3 or 4 dimensions with the Wilson action, using heat-bath plus over-relaxation updates. This is the same method as lattice QCD, with SU(2) instead of SU(3) and modest volumes.

```js
import { createLattice, sweep, plaquette, wilsonLoop, creutzRatio } from 'psim:core/lattice-gauge.mjs';
const lat = createLattice({ L: 6, d: 4, beta: 2.3, seed: 1 });   // 6⁴: ~15 ms per sweep in Node
for (let i = 0; i < 50; i++) sweep(lat);                          // thermalise
// measure W(R,T) over many sweeps; σa² ≈ χ(R,T) = −ln[W(R,T)W(R−1,T−1)/(W(R,T−1)W(R−1,T))]
```

**What to show:**
- Wilson loops W(R,T) falling with the **area** RT (confinement), not the perimeter, with an exponential fit.
- Creutz ratios approaching the string tension σa².
- The static potential V(R) = −ln[W(R,T)/W(R,T−1)] rising linearly. Overlay a Cornell fit to the **measured** points.
- A live view of the gauge field: the plaquette action density on a 2D slice, coloured, updating each sweep. This is the actual gluon field configuration, and it's a far better visual than drawn springs.
- The plaquette history, showing thermalisation; statistical error bars from `meanWithError` (blocking).

**Required checks in `verify()`** (all in `tests/lattice-gauge.test.mjs`):
- 2D exact solution ⟨P⟩ = I₂(β)/I₁(β), W(R,T) = ⟨P⟩^{RT};
- 4D strong coupling ⟨P⟩ ≈ u + 4u⁵;
- gauge invariance of observables;
- over-relaxation preserves the action;
- positive Creutz ratio (area law).

Run the Monte Carlo in chunks per frame (one sweep per frame for 6⁴) so the page stays responsive. Use a Web Worker if you go larger.

**One-line limits for the page:** SU(2) not SU(3), no dynamical quarks (so no string breaking), finite lattice spacing and volume. The trend is real; continuum numbers need extrapolation.

## 2. Analytical context (tier `analytical`), next to the lattice

| Topic | What to compute | Source to cite |
|---|---|---|
| Static potential fit | Cornell V(r) = −(4/3)α_s ħc/r + σr, σ ≈ 0.18 GeV² ≈ 0.9 GeV/fm | Eichten et al., PRD 17, 3090 (1978) |
| Physical scale | Sommer scale r₀ ≈ 0.5 fm: r²F(r) = 1.65 sets the lattice spacing a in fm | Sommer, NPB 411, 839 (1994) |
| String breaking (needs quarks) | Breaks near 1.2–1.4 fm | Bali et al., PRD 71, 114513 (2005) |
| Asymptotic freedom | One-loop α_s(Q²) = 12π/((33−2n_f) ln(Q²/Λ²)); α_s(M_Z) = 0.1180 | PDG |
| Proton structure | Parton distributions from published fit tables (cite set and version) | NNPDF / CT18 / MSHT |

## Rules

- Colour charge drawn as red, green or blue is a symbol; say so once.
- A separation slider on a formula is a "parameter sweep", not time evolution; label it.
- Quote numbers with source and uncertainty.
- State the one-line limits at the end, not as the opening.
