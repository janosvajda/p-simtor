# Higgs fields and mass generation

Status: **compose guidance**, not a built-in solver or recipe. Read this for Higgs fields, electroweak symmetry breaking or elementary-particle mass generation. First determine whether the user wants the potential and mass mechanism, propagation of disturbances, or collider production/decay. These are different models.

## Potential and mass relations

For the Standard Model, explain the complex scalar doublet and its potential V(Φ) = μ²Φ†Φ + λ(Φ†Φ)². For λ > 0 and μ² < 0, the vacuum has Φ†Φ = v²/2 with v² = −μ²/λ. In a specified gauge, write the neutral component as (v + h)/√2. At tree level, m_h² = 2λv²; charged-fermion masses satisfy m_f = y_f v/√2, and m_W = gv/2, m_Z = √(g² + g′²)v/2, in natural units. Couplings are inputs, not predictions of this display. See [PDG, Higgs review, §11.2](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-higgs-boson.pdf).

An analytical potential and mass-relation view can answer the mechanism question without inventing particle trajectories. If the user asks for collider rates, lifetimes or events, use the relevant quantum calculation or cited experimental/generator data; a potential plot cannot predict them.

## Classical scalar disturbances, when that reduction is accepted

A useful reduced model with a canonically normalised real scalar φ is:

```text
V(φ) = λ(φ² − v²)² / 4,       λ > 0, v > 0
∂²φ/∂t² − ∇²φ + λφ(φ² − v²) = 0       (c = ħ = 1)
E = ∫ [ (∂φ/∂t)²/2 + |∇φ|²/2 + V(φ) ] d^d x
```

This real-scalar model has two minima, φ = ±v. Linearising φ = v + h gives a curvature mass m_eff² = V″(v) = 2λv² and dispersion ω² = k² + m_eff². These follow from the stated potential and provide independent small-amplitude benchmarks.

Name it **classical scalar model inspired by the Higgs potential**, with its spatial dimension and units. It omits the full electroweak doublet, gauge fields and quantum fluctuations. Its classical waves are not simulated individual Higgs bosons, and its two-vacuum structure is not the full Standard Model vacuum structure. If damping is added, identify the bath/phenomenological damping and account for lost energy.

## Explain the geometry and colour

- A potential surface plots **field values** and potential energy, not a landscape located in ordinary space. A moving marker shows a computed homogeneous field amplitude only if its equation was integrated; otherwise it is an illustration.
- A spatial field plot uses position on the horizontal axes and amplitude/energy density as height or colour. Label the quantity, scale, units and whether the field value depends on a gauge choice.
- Use an adjacent caption such as “Colour encodes field amplitude; this is a classical scalar reduction.” Do not imply that the Higgs field has a visible colour or behaves as a frictional medium.
- Explain that field–mass coupling changes the equations. Do not demonstrate mass generation by arbitrarily slowing drawn particles. The Higgs boson is an excitation associated with the field, not the field itself; see [CERN's Higgs explanation](https://home.cern/science/physics/higgs-boson/what/).

## Checks matched to the selected model

For an analytical view, check stationary points, potential derivatives and the stated mass relations. For scalar dynamics, verify the homogeneous small-amplitude frequency and a spatial mode's dispersion, energy behaviour under the specified boundaries, and independent spatial/time refinement. Keep undamped energy conservation separate from a damped run's energy balance.

Use `analytical` for evaluated relations within the named model, and `approximate` for the classical reduction. Label any pedagogical analogy individually. Neither label means a full quantum Standard Model simulation has been built.
