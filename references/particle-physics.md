# Quarks, gluons and QCD

Status: **compose**, tiers `analytical` or `conceptual` only.

## What is and isn't computable here

Quantum chromodynamics at the hadron scale is non-perturbative. Quantitative results (hadron masses, the static potential, flux-tube profiles, parton distributions) come from lattice QCD on supercomputers or from global fits to experimental data. A browser page cannot solve QCD. Quarks inside a proton have **no classical trajectories**. "Three little balls connected by springs bouncing around" is not a simulation of QCD and must never be presented as one.

Say this to the user plainly, then build one of the defensible alternatives below.

## Defensible visualisations

| Topic | What to compute | Tier | Source to cite |
|---|---|---|---|
| Confinement | Cornell static potential V(r) = −(4/3) α_s ħc / r + σ r, with σ ≈ 0.18 GeV² (≈ 0.9 GeV/fm), α_s ≈ 0.3 (static-potential fits) | analytical (phenomenological fit to lattice and quarkonium data) | Eichten et al., Phys. Rev. D 17, 3090 (1978); lattice determinations of σ |
| Force between quarks | F = −dV/dr: tends to a constant ≈ 0.9 GeV/fm ≈ 14 tonnes-force at large r | analytical | Same |
| String breaking | Energy stored σr reaches ~2 m_light-meson threshold near r ≈ 1.2–1.4 fm, where creating a q̄q pair is favoured | conceptual, with the threshold computed from cited numbers | Bali et al., Phys. Rev. D 71, 114513 (2005) |
| Flux tube | Chromoelectric field concentrated in a tube of width ~0.3–0.5 fm between static sources; draw a profile from published lattice data | conceptual (shape from cited lattice results, not computed) | Bali, Schlichter & Schilling (1995); Cardoso et al. (2013) |
| Asymptotic freedom | One-loop running coupling α_s(Q²) = 12π / ((33 − 2n_f) ln(Q²/Λ²)); compare with α_s(M_Z) = 0.1180 | analytical (perturbative; invalid below ~1 GeV, so shade that region) | PDG Review "Quantum chromodynamics" |
| Proton structure | Parton distribution functions x f(x, Q²) for u, d, g, sea, from published fit tables | analytical display of published data (cite the set and version; never invent shapes) | CTEQ/NNPDF/MSHT; PDG |
| Colour charge | SU(3) colour algebra: 8 gluons, colour-singlet combinations (qqq, qq̄) | conceptual (symbolic colours, disclosed: colour charge is not visual colour) | Any QCD text |
| Scattering | Tree-level cross-sections (for example e⁺e⁻ → qq̄, R-ratio = N_c Σ Q_q²) | analytical | PDG |

## Rules for these pages

- The tier badge says `conceptual` or `analytical`, and the model panel says what was **not** computed ("No QCD dynamics are simulated; the potential is a fit to lattice results").
- Colour charge shown as red, green or blue is a symbol. Say so.
- If you animate anything (for example a stretched flux tube as separation increases), the animated parameter is r, and every frame's energy and field profile is evaluated from the cited formula or data at that r. The motion itself has no dynamical meaning, so label it "parameter sweep, not time evolution".
- Numbers quoted (σ, α_s, Λ, string-breaking distance) carry their source and uncertainty.
- Invite the user to specialised tools for real calculations: lattice QCD codes (openQCD, Chroma), PDF libraries (LHAPDF), event generators (Pythia), which can't run in a page.
