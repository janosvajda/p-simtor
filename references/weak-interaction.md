# The weak interaction

Status: guidance plus generic tools (`core/kinematics.mjs` for 4-vectors, decays and phase space; `PDG` in `core/constants.mjs` for inputs). No ready-made recipe. Read [particle physics](particle-physics.md) for QCD and [Higgs](higgs.md) for electroweak mass relations.

"The weak force" is not a push or pull you can watch between objects. It is carried by the W and Z bosons, its range is ħ/(m_W c) ≈ 2.5×10⁻¹⁸ m, and it is seen through **decays, flavour change and parity violation**. Say this in the Fact check line when the request reads like a force field to animate, then pick one computable process.

## Candidate models (pick by the question)

| Question | Model | Observables and benchmark |
|---|---|---|
| What does the weak force do? Why does it break mirror symmetry? | Polarised muon decay μ⁺ → e⁺ν_eν̄_μ in V−A (Fermi) theory: purely leptonic, no hadronic uncertainty | Michel spectrum dΓ/dx ∝ x²(3 − 2x) (massless e, ρ = 3/4); angular asymmetry x²[(3 − 2x) + P cosθ (2x − 1)] (ξ = 1, δ = 3/4), so ⟨cosθ⟩ = P/9 for μ⁺ and −P/9 for μ⁻; Michel parameters against `PDG.michelRho`, `michelDelta`, `michelXi`, `michelEta` |
| Radioactive β decay | Fermi theory with the Fermi function for the nuclear charge | Electron spectrum and Kurie plot; endpoint energy from mass differences (cite the nuclide data) |
| How strong is it, and why "weak"? | G_F ≈ √2 g²/(8 m_W²) at tree level; the amplitude ratio of W exchange to photon exchange versus momentum transfer | Show it is weak only at low Q² (≪ m_W²) and comparable to electromagnetism near Q² ~ m_W². The tree-level G_F from m_W misses the radiative correction Δr (a few %): label it |
| Neutrino scattering | Four-fermion cross-sections at low energy, W/Z propagators above | σ ∝ G_F² s at low energy; growth stops at s ~ m_W² |
| The Z boson | Breit–Wigner line shape of e⁺e⁻ → Z → f f̄ | Peak at m_Z with width Γ_Z; number of light neutrinos N_ν = 2.9840 ± 0.0082 (ALEPH, DELPHI, L3, OPAL, SLD, Phys. Rept. 427, 257 (2006)) |

## Method

- **Event sampling:** draw each decay or scattering event from phase space and weight it by the squared matrix element: `rambo(n, √s, rand, masses)` (uniform n-body phase space, exact 4-momentum conservation) or `twoBodyDecay(M, m1, m2, rand)`, then `acceptReject(sample, weight, wmax, rand)`. `acceptReject` throws if a weight exceeds `wmax`, so a wrong bound cannot bias the sample silently. Derive `wmax` analytically and test it.
- **Frames:** work in the rest frame where possible; `boost(p, β)` and `restFrameBoost(P)` move between frames, and `dot4` is invariant under both.
- **Decay times:** exponential with the measured lifetime (`PDG.tauMu`, `PDG.tauPiCharged`, `PDG.tauKCharged`).
- **Tracks in a detector** (bubble or cloud chamber, emulsion): position-space paths are computable (Lorentz force plus the mean energy loss `betheStoppingPower` from `core/matter.mjs`); neutrinos are straight lines that leave no trace, drawn faint and labelled not observable. Spin precession in a field B at ω = g e B/(2m) uses the muon anomaly `PDG.amu`.
- **Accumulating histograms:** `createHistogram({ bins, range })` from `core/data.mjs` grows event by event with Poisson errors; compare with the closed form by `chiSquareTest` and the spectrum shape by `ksTest`.

## Verification

- Exact kinematics: 4-momentum conserved to round-off in every event; final-state masses as requested.
- Sampled spectra against closed forms: `ksTest` per particle species (p-values, not just a figure).
- Asymmetries against ±1/9 (or the energy-resolved value), with `meanSem` errors and pulls.
- Michel parameters refitted from the sample, compared with the V−A values and with `PDG` measurements.
- **Circularity:** G_F is itself extracted from the muon lifetime (MuLan: Tishchenko et al., PRD 87, 052003 (2013)). Recomputing τ_μ from G_F is a consistency check of the formula and corrections, not an independent prediction: say so.
- Determinism: the same seed with different step schedules gives the same state.

## Representation

- A decay is drawn in **momentum space** (arrows from the decay point, length ∝ |p|), not as particles flying through the room, unless you compute positions. Say which space the view shows.
- Neutrinos are not detected: draw them dashed or faint and say "not observable; reconstructed from conservation" in the legend.
- Show the spin (or polarisation) direction and the measured asymmetry together, so the mirror asymmetry is visible: under a parity mirror the momentum flips but the spin does not.
- Historical context when relevant: parity violation was observed in ⁶⁰Co β decay (Wu et al., Phys. Rev. 105, 1413 (1957)) and in muon decay (Garwin, Lederman & Weinrich, Phys. Rev. 105, 1415 (1957)).

**Visible scope:** "Tree-level V−A theory of one process; electron mass and radiative corrections neglected (percent level near the spectrum endpoint)." Say which parts of the electroweak theory are not modelled.
