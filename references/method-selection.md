# Model and method selection

Apply this workflow to any requested scientific system. Select a mathematical problem and measurable outputs before looking for a recipe. Optional domain notes at the end provide shortcuts; they are not a catalogue of what the agent may simulate.

## 1. Establish the scientific aim

Identify what the user wants to explain, compare or predict, at what scale and fidelity. Determine which question an interaction should answer. A phenomenon's name does not specify its observables, experiment or interpretation.

Ask when alternatives change the scientific result, relevant theory or scope. Preserve given choices and resolve routine implementation details yourself. If the equations, quantities or experimental conditions are unknown, obtain or research them; do not fill the gap with whatever existing recipe looks similar.

Check factual premises separately from ambiguity. A hypothetical system is valid when its assumptions are identified. Correct an impossible premise and agree on the resulting scientific scope before building a different system.

## 2. Write the mathematical specification

Record these decisions in the physics module and page's model panel; a separate planning document is optional.

| Decision | What to establish |
|---|---|
| State | Variables actually calculated: scalar/vector/tensor fields, finite coordinates, amplitudes, populations, distributions, graph states or other justified objects |
| Law or data | Governing equations, constraints, transition rules, target probability measure or sourced dataset; signs, coordinate/basis conventions and approximations |
| Inputs | Parameters, units, dimensionality, spatial domain, initial/boundary conditions, external drives, environment or ensemble where relevant |
| Observables | Functions or estimators of the state that answer the scientific question; normalization, units, measurement interpretation and uncertainty |
| Progress axis | Physical time, sample/sweep number, iteration, parameter value, or a static view; iteration is not automatically physical evolution |
| Representation | The meaning of each axis, geometric object, colour, glyph and displayed number; physical quantities versus symbols |
| Validity and evidence | Applicable scale/regime, neglected processes, reference results, required accuracy, resource constraints and checks supporting each claim |

Separate the three operations:

```text
calculation:     equations, samples or source data -> state
observation:     state -> requested quantities and uncertainty
representation: state/observables -> pixels, labels and interaction
```

The renderer reads results; it does not invent or alter the state. A model may have no physical trajectories or no physical time at all. For a data view, identify interpolation/filtering separately from a newly simulated result. Choose what the mathematics means before deciding what it looks like.

## 3. Establish sources and validity

Use primary scientific literature, authoritative documentation or user-supplied equations/data. Verify units, normalization, sign conventions, assumptions and the regime in which the model answers the question. Locate a benchmark or independently derivable limit before coding.

Use checked repository constants and algorithms when they apply. For missing inputs, record the source, units, uncertainty and version/date when relevant. Clearly identify illustrative parameter values. Do not invent empirical coefficients, measurements or a new physical law to complete an animation.

If sources conflict or the model does not answer the requested question, resolve the issue or narrow the claim explicitly. A convenient effective model, fitted curve or dataset is not an automatic substitute for a requested theory. Ask before a material scope change.

## 4. Select a method by mathematical structure

These criteria apply across subjects. A name in the request is not a solver choice.

| Mathematical structure | Selection criteria | Evidence to seek |
|---|---|---|
| Closed-form relationship | Evaluate directly when it answers the question; identify its assumptions and domain | Identities, derivatives, normalization and limiting cases |
| Finite-dimensional evolution | Identify stiffness, Hamiltonian structure, constraints, dissipation and characteristic times; choose a compatible integrator | Independent trajectory/solution, applicable invariants or balances, timestep refinement |
| Spatial evolution | Identify transport/wave/diffusion structure, geometry, boundary conditions and stability restrictions; select spatial and time discretizations together | Known mode or manufactured solution, boundary checks, independent spatial/time refinement |
| Stationary/eigenvalue/constraint problem | Match the operator and constraints; control conditioning, normalization and residuals | Residuals plus independent limits/eigenvalues, resolution and constraint checks |
| Stochastic physical dynamics | Derive transition rates or noise from the model; preserve the meaning of elapsed time and the correct measure | Rate/distribution benchmarks, ensemble statistics, applicable time/space refinement |
| Equilibrium sampling | Specify the target distribution and updates that sample it; diagnose warm-up and correlations | Distribution benchmarks, effective independent information, sample-length and seed studies |
| Statistical observable or integral | Specify the state space, estimator, normalization and approximation/truncation | Independent identity/limit, uncertainty, sample/bin/truncation or quadrature refinement |
| Empirical data or fitted model | Preserve provenance and uncertainties; distinguish observations, interpolation, fit and extrapolation | Independent validation data and stated fit/residual uncertainty; no solver-accuracy claim for imported data |

Hybrid systems may combine rows. Identify which controls and diagnostics apply to each part, and check the coupling/interface as well as individual components. A method being available in `core/` is useful only if its assumptions match the specification.

Concrete reusable implementations of general mathematics are catalogued in [math-toolkit.md](math-toolkit.md). Check their conditions and reported convergence before writing another implementation or using a result.

## 5. Check feasibility without changing the question

Estimate degrees of freedom, memory, spatial/temporal resolution, equilibration/sampling work and rendering cost. Bound user inputs by physical validity and a measured compute budget. Test demanding allowed cases, not just a visually pleasing initial state.

Choose a unit system that keeps relevant scales numerically manageable and display its physical meaning. Set resolution from the phenomenon and error requirement, not from frame rate. When a compute budget is reached, reduce playback throughput rather than silently changing the equations or timestep.

If the requested fidelity needs another backend or more resources, present the concrete choices: the requested calculation with available resources, a sourced data view, or a reduced model with specific lost capabilities. Get agreement before changing scientific scope; do not claim the full request was solved by the reduction.

For dynamics, label physical time. For sampling, label samples/sweeps and define which results have equilibrated. For a scan or static calculation, label the independent parameter and avoid artificial physical-time playback.

## 6. Derive the representation and checks from the specification

Every visual channel needs a meaning: source quantity, mapping/scale, units where applicable, coordinate space and uncertainty/aggregation. Classify it as a computed quantity, imported observation or disclosed symbol. Put conventions that affect interpretation beside the view, with detail in the model panel. A palette does not become a literal physical property because it resembles a photograph.

Choose checks for the stated claims using [verification](verification.md). Conservation is required only where the equations conserve that quantity; otherwise check the correct source/loss balance. Refinement tests discretization, sampling tests estimates, and measurement comparisons test physical applicability. Do not let one stand in for the others.

A verification result attaches to the model, parameters, conditions and tested regime. If those change, withdraw unsupported claims or recompute the relevant evidence. If a signal cannot be resolved, report the missing evidence and withhold the claim instead of forcing the expected shape.

Use the [compose contract](compose.md) to implement the calculation, observables, tests and page. Reuse a matching solver or suitable UI components independently; a new model does not require a new topic-specific skill guide.

## Optional domain notes

Read a note only when it helps the model already selected. If none fits, use the same workflow with researched sources.

- [Gravity](gravity.md)
- [Electromagnetism](electromagnetism.md)
- [Radiation and relativity](radiation-relativity.md)
- [Particle physics](particle-physics.md)
- [Scalar/Higgs examples](higgs.md)
- [Molecular dynamics](molecular.md)
- [Entropy examples](entropy.md)
- [Quantum mechanics](quantum.md)
