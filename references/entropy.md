# Entropy, molecules and ensembles

Status: **compose guidance**, not a built-in recipe or entropy solver. Read this for entropy, mixing, irreversibility or molecular state populations. “Entropy of molecules” needs a defined state space: gas expansion, internal rotational/vibrational states and interacting-fluid entropy are different questions. Ask which one matters when the request is ambiguous.

## Define what the entropy measures

Name the system, ensemble/bath, accessible states, entropy definition, units and estimator before plotting a number. Distinguish dimensionless information entropy H from physical entropy S in J/K or units of k_B. Particle speed, visual messiness and entropy are not interchangeable.

For isolated Hamiltonian dynamics, fine-grained Gibbs entropy of the full phase-space ensemble is conserved. Coarse-grained or macroscopic entropy can increase as the distribution spreads and correlations become unresolved; finite systems can fluctuate. Do not force every frame's entropy to increase. See [Chakraborti et al., entropy during ideal-gas free expansion](https://arxiv.org/abs/2109.07742).

## Gas expansion or mixing

For an accepted ideal-gas demonstration, evolve particles with specified boundaries, then compute occupancies from their actual positions. For M equal-volume spatial bins, p_i = n_i/N and:

```text
H_position = −Σ_i p_i ln p_i       (natural logarithm; 0 ln 0 = 0)
0 ≤ H_position ≤ ln M
```

This is a coarse-grained **positional** information entropy. Under an independent-particle approximation, N k_B H_position is a positional contribution with a bin-dependent reference; it omits momenta, internal states and correlations. Never call it the total thermodynamic entropy of an interacting fluid. Keep the bin grid and physical domain fixed when comparing times.

For equilibrium endpoints of an ideal gas at unchanged temperature and particle number, the thermodynamic change is ΔS = N k_B ln(V_f/V_i). Compare that endpoint prediction with the positional estimator only under its assumptions, with finite-sampling uncertainty. Do not insert the endpoint formula as a time-dependent result before the gas equilibrates.

For mixing distinct ideal species at common temperature and pressure, ΔS_mix = −k_B Σ_s N_s ln x_s, where N_s is the particle count of species s and x_s its mole fraction. Artificially painting otherwise identical molecules red and blue does not make them distinct species or create thermodynamic entropy of mixing. Specify what makes species physically distinguishable; see [MIT, ideal gas, mixing and Gibbs' paradox](https://ocw.mit.edu/courses/8-333-statistical-mechanics-i-statistical-mechanics-of-particles-fall-2013/0cdb670a825c6aaea7e229ac40776cbf_MIT8_333F13_Lec13.pdf).

## Internal molecular states

For a specified molecule and cited discrete levels E_i with degeneracies g_i, a canonical internal-state calculation uses:

```text
β = 1/(k_B T)
Z = Σ_i g_i exp(−β E_i)
P_i = g_i exp(−β E_i)/Z
U = Σ_i P_i E_i
S_internal = k_B (ln Z + βU)
```

P_i is the probability of an entire energy level. If calculating entropy as a probability sum, include its degeneracy: S_internal = −k_B Σ_i P_i ln(P_i/g_i). A plain −Σ_i P_i ln P_i would discard states within degenerate levels. See [MIT, canonical partition function and degeneracy](https://ocw.mit.edu/courses/5-62-physical-chemistry-ii-spring-2008/bda82786b6014d3f468798339f658629_03_562ln08.pdf).

Identify which translational, rotational, vibrational or electronic degrees of freedom are included and how many levels were retained. For liquids or flexible molecules, an occupancy or speed histogram alone is inadequate: use a sourced method appropriate to the requested entropy, such as free-energy/thermodynamic integration, or explicitly present a limited estimator.

## Visuals and checks

- Show the sampled distribution beside its entropy curve. Label colour as species, speed, energy or an artificial tag; a colour key must not silently define physical distinguishability.
- For expansion, benchmark particle number and energy under the specified boundaries, histogram normalisation, a one-bin distribution (H = 0) and uniform occupancy (H = ln M). Study bin count and sample size separately from integration-step convergence; report fluctuations across independent seeds and temporal correlations where relevant.
- For molecular states, check Σ_i P_i = 1, agreement of the probability and partition-function entropy formulas, stability as more levels are included, and low-temperature behaviour for the specified ground-state degeneracy.
- Identify system–bath heat/work accounting for driven or thermostatted examples. The system's entropy need not increase when it exports entropy to its surroundings.

Use `analytical` for exact relations within the named model and `approximate` for finite-sample or reduced estimates. State the measurement scope next to the plot. No algorithm should smooth, clamp or monotonise the entropy curve to imitate the second law.
