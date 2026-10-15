# Fact check before building

Requests often contain a wrong premise. Building it silently spreads the error; correcting it is part of the job. Check every object, number and claim in the request against what is known. When something is wrong, say so in the reply's "Fact check:" line, then build the closest real thing or a clearly labelled hypothetical.

Use data from `core/` where it exists: planets' orbits and moon counts (`core/solar-system.mjs`), constants (`core/constants.mjs`). Don't type facts from memory when a checked value is available.

## Common wrong premises (and what to do)

| The request assumes | The fact | Build instead |
|---|---|---|
| Mercury or Venus has moons | Neither has a moon | The planet alone, or a clearly labelled hypothetical moon |
| Pluto is the ninth planet | Dwarf planet since 2006 (IAU) | Include it only as a dwarf planet, if asked |
| Electrons orbit the nucleus like planets | Electrons are probability clouds (quantum mechanics) | Orbital probability densities; a Bohr model only as a labelled historical picture |
| Black holes suck things in | Far away, a black hole pulls like any mass of the same size; orbits stay stable outside about 3 r_s | Orbits plus the true capture region (photon sphere, ISCO) |
| The Sun is burning (fire) | Nuclear fusion in the core | Fusion, if energy production is the topic |
| Seasons come from Earth's distance to the Sun | Seasons come from the 23.4° axial tilt (Earth is closest in January) | Tilt and sunlight angle |
| The Moon's far side is always dark | It gets as much sunlight as the near side; it always faces away from Earth | Tidal locking, sunlit phases |
| There is no gravity on the ISS | Gravity there is about 90 % of surface gravity; astronauts are in free fall | Orbit as continuous free fall |
| Heavier objects fall faster | Without air, all fall equally (equivalence principle) | Free fall with optional air drag |
| Light from distant stars arrives instantly | Light travels at c: 8.3 min from the Sun, 4.2 years from Proxima Centauri | Light-travel delay |
| Sound travels in space | Sound needs a medium | Light or other signals |
| Galaxies are held together by their visible stars | Visible matter is far too little: dark matter (or modified gravity) is needed | The rotation-curve comparison |
| Relativity matters only near light speed | GPS clocks need it daily (38 µs/day) | Clock-rate effects at everyday speeds and gravity |
| Quarks are literally red, green or blue | Colour charge is a physical quantum property; the colour names and screen pigments are symbolic, not optical colours | Label the charge diagram or measured field quantity beside the view; see [QCD](particle-physics.md) |
| The Higgs field gives mass by dragging particles to a stop | Mass terms arise through field couplings; a friction animation does not model the Higgs mechanism | Potential/mass relations or an explicitly reduced scalar-field calculation; see [Higgs](higgs.md) |
| Prove or solve an open problem (a Millennium Prize problem such as Yang–Mills existence and mass gap, Navier–Stokes regularity; an unsolved equation) | It is open; no simulation or finite computation proves a statement about the continuum or infinite volume | Numerical evidence for a regularised, named version (e.g. a lattice mass gap at stated spacing and volume) compared with a published calculation; say "evidence, not a proof" on the page and never write "proved" |
| Gauge-field link values (or potentials) show the field, and they settle near the identity as the coupling weakens | Links are gauge-dependent: unfixed, each is exactly Haar-distributed with ⟨U⟩ = 0 at every coupling (Elitzur's theorem); only closed loops are physical | Colour plaquettes or Wilson loops; show links only as gauge-dependent geometry, or in a named fixed gauge (e.g. Landau); see [particle physics](particle-physics.md) |
| The weak force is a push or pull you can watch acting between objects | It acts over ~10⁻¹⁸ m via W and Z bosons; it is seen through decays, flavour change and parity violation | One computable weak process (muon decay, β decay, a cross-section) with its observables; see [weak interaction](weak-interaction.md) |
| More visually chaotic molecules always mean more entropy | Entropy depends on the stated ensemble, state space and coarse-graining; fluctuations and exported entropy matter | A defined probability/partition-function calculation with scope and uncertainty; see [entropy](entropy.md) |

## Rules

- **Correct kindly and briefly**, one sentence, then move on to what can be built.
- **Clarify ambiguity separately from errors.** A broad topic name ("the weak force", "quantum gravity") is not a false premise, but name the common misconception it invites if the page would otherwise feed it, as context rather than as a correction. Ask about the scientific aim instead of correcting a claim the user never made. Check unfamiliar claims against sources even when they have no entry in the example table.
- **Never invent objects** as if they were real. A hypothetical is fine when it is labelled in the page and in the reply.
- **Uncertain or changing facts** (moon counts, new measurements): say "as of <year>, per <source>".
- **A value not in `core/` and not checked against its source in this session:** cite the paper, and say in the reply that the value was quoted from memory and not re-checked, so the user can verify it.
- **When the user insists on a counterfactual** ("what if Mercury had a moon?"), build it gladly. Label it as hypothetical, and point out which physics is real (the orbital mechanics) and which part is invented (the moon).
