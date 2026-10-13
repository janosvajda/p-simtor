# Build p-simtor — An Installable Scientific Physics Simulation Skill

You are an expert AI software engineer, computational physicist, numerical mathematician, and scientific visualization developer.

Your task is to create **p-simtor**, a fully functional, independently designed, installable AI agent skill for generating mathematically accurate, scientifically defensible, and visually exceptional interactive physics simulations.

This is a real software development task. Build the actual skill, not a proposal, mockup, or standalone demonstration website.

## 1. Mandatory reference: Anatomy

**Reference repository:**
https://github.com/wheresryan22/anatomy

Before designing or implementing p-simtor, thoroughly examine this repository.

Understand:

- What makes Anatomy an AI agent skill rather than a conventional application.
- How users install it directly from GitHub.
- How it is invoked through Claude Code, Codex, and OpenCode.
- How its SKILL.md guides an agent through creating working artifacts.
- How its instructions, reusable capabilities, and supporting tools work together.
- How it generates complete, standalone interactive outputs.
- How it achieves high visual quality and verifies generated work.
- How it packages its functionality for reuse across different projects.

**Use Anatomy as a reference for the type of product, installation experience, agent workflow, craftsmanship, and quality standard.**

Do not copy Anatomy's directory structure, implementation, internal abstractions, visualization methods, or design decisions without independent technical justification.

Anatomy focuses on interactive illustrations of machines and mechanisms.

**p-simtor focuses on computational physics, particle dynamics, mathematical models, and scientific simulations.**

These are fundamentally different technical domains and require independently designed systems.

The goal is to achieve a comparable level of skill usability and refinement while creating entirely original scientific capabilities.

## 2. Non-negotiable requirement: An installable agent skill

p-simtor MUST be an actual installable Agent Skill, following the open Agent Skills specification.

It must NOT be:

- Merely a physics simulation website.
- A standalone graphical application advertised as a skill.
- A library without usable agent instructions.
- A collection of example animations.
- A generic prompt with no supporting implementation.

The final GitHub repository must itself be directly installable as a skill, just as Anatomy is.

After publishing to GitHub, a user should be able to run a command in this form:

```bash
npx skills add <github-owner>/p-simtor -g -a claude-code
```

And similarly install it in compatible versions of Codex and OpenCode.

The installed skill should be discoverable and invocable as:

**Claude Code / OpenCode**
```
/p-simtor simulate a gravitational three-body system
```

**Codex**
```
$p-simtor simulate a hydrogen atom's electron probability density
```

Follow each agent's actual skill invocation conventions.

Create a properly specified SKILL.md containing the required metadata, a useful description for automatic skill discovery, instructions for invoking the skill, and a complete workflow for creating scientific simulations.

Determine the remaining packaging and installation requirements by researching the Agent Skills standard and the current supported agent implementations.

The skill must carry or correctly provision whatever resources it needs to generate working simulations.

**Installation must be tested, not merely documented.**

Verify that the skill can be installed from its repository, discovered by a compatible agent, invoked using natural language, and used to produce a working simulation in a fresh project.

If specific agent environments are unavailable, distinguish actual verification from documented but untested compatibility.

## 3. Core vision

p-simtor transforms descriptions of physical phenomena into computationally grounded interactive visualizations.

A user describes a physical system.

The AI agent uses p-simtor to:

1. Identify the physical phenomenon.
2. Research the correct governing theory.
3. Select an appropriate mathematical representation.
4. Define the model's assumptions and applicable scales.
5. Choose reliable numerical or analytical methods.
6. Calculate the physical evolution or scientific observables.
7. Verify the results against mathematical or scientific benchmarks.
8. Produce an exceptional interactive visualization.
9. Provide meaningful controls and measurements.
10. Deliver a complete, reusable, runnable simulation.

The agent should handle different requests intelligently rather than relying on a finite collection of predefined animations.

**The mathematics determines the physical behavior. The visuals reveal that behavior.**

## 4. Scientific scope

Design p-simtor to support scientific phenomena across scales.

### Classical and gravitational dynamics

- Two-body and three-body problems.
- N-body gravitational interactions.
- Orbital motion and stability.
- Chaotic dynamical systems.
- Particle trajectories.
- Collisions and close encounters.
- Energy and momentum conservation.
- Phase-space and trajectory analysis.

### Atomic and molecular physics

- Atomic interactions.
- Molecular dynamics.
- Intermolecular potentials.
- Thermal motion.
- Vibrational behavior.
- Molecule rotation.
- Particle collisions.
- Statistical mechanics.
- Physically meaningful molecular configurations.

Use validated force fields where quantitative molecular accuracy is claimed.

### Quantum and atomic systems

- Atomic orbitals.
- Electron probability densities.
- Wavefunctions.
- Quantum states.
- Energy levels.
- Quantum superposition.
- Wavepacket propagation.
- Quantum tunneling.
- Probability distributions and time evolution.

Do not visualize electrons as classical orbiting planets while presenting the result as correct quantum mechanics.

### Electromagnetism and nuclear dynamics

- Charged-particle motion.
- Electric and magnetic fields.
- Lorentz-force trajectories.
- Proton and ion dynamics.
- Nuclear interactions using appropriately validated models.
- Relativistic corrections where required.

### Quarks, gluons, and particle physics

- Proton structure.
- Quark and gluon interactions.
- Color charge.
- Confinement.
- Quantum chromodynamics.
- Particle interactions and scattering.
- Appropriate field and probability visualizations.

Recognize that quantitatively solving QCD and visualizing quark/gluon dynamics are substantially different from classical particle integration.

Do not fake quark trajectories, invent interaction physics, or claim that simplified classical models reproduce quantum field theory.

When exact simulation is computationally infeasible, use defensible approximations, established theoretical results, or validated research datasets, and identify the visualization's scientific status.

**Each domain must use mathematics appropriate to its physics. Do not force all scales into a single classical particle model.**

## 5. Mathematical accuracy

Scientific correctness is p-simtor's highest priority.

All simulations must be based on recognized physical theories.

Use appropriate mathematical methods depending on the problem, potentially including:

- Newton's equations of motion.
- Hamiltonian and Lagrangian mechanics.
- Numerical solutions of differential equations.
- Symplectic and adaptive integrators.
- Molecular force-field calculations.
- Schrödinger equation solvers.
- Analytical quantum wavefunctions.
- Electromagnetic field equations.
- Relativistic dynamics.
- Validated approximations to quantum field models.

These are examples, not instructions to adopt a predetermined numerical architecture.

Research and select methods based on their accuracy, computational cost, stability, and applicability.

Correctly handle:

- Physical constants.
- Units and dimensional consistency.
- Initial conditions.
- Boundary conditions.
- Numerical precision.
- Timestep selection.
- Error accumulation.
- Conservation laws.
- Solver stability.
- Reproducibility.

Numerical simulation speed must be independent of the display's frame rate.

Do not animate fabricated paths and describe them as computed physics.

Do not allow rendering effects to modify scientific quantities.

## 6. Scientific verification

Every generated simulation must communicate how trustworthy its results are.

Distinguish:

1. Validated quantitative simulations.
2. Approximate numerical simulations.
3. Analytical scientific visualizations.
4. Conceptual illustrations.

For numerical simulations, implement applicable checks such as:

- Energy conservation.
- Momentum conservation.
- Angular momentum conservation.
- Convergence under timestep refinement.
- Analytical reference comparisons.
- Known solutions and published benchmarks.
- Probability normalization.
- Model-specific error diagnostics.

There is no single verification method appropriate for all physical phenomena.

Select checks according to the underlying mathematics.

The generated simulation should make its scientific assumptions and limitations accessible to the user.

Never claim numerical accuracy without evidence.

## 7. Exceptional visual quality

Visual excellence is a fundamental part of p-simtor.

Every output should combine the precision of research-grade scientific visualization with the beauty of sophisticated digital art.

Aim for a distinctive visual identity:

- Beautiful interactive 2D and 3D scenes.
- Precise and elegant orbital trajectories.
- Refined particle representations.
- Sophisticated lighting and shading.
- Smooth, physically meaningful motion.
- Beautiful vector fields and field-line rendering.
- Wavefunction and probability-density visualization.
- Carefully composed scientific graphics.
- Responsive layouts.
- Minimal, highly refined controls.
- Clear mathematical annotations.
- Professional light and dark presentation.

Use shaders, volumetric techniques, high-quality rendering, and GPU acceleration where they meaningfully improve the scientific output.

Do not use a technique simply because it looks impressive.

The visualization method must suit the physical phenomenon.

A gravitational orbit, a molecular trajectory, a quantum orbital, and a quantum field should not be forced into the same generic particle animation style.

The user should immediately understand what physical quantities are being represented.

If an important visual element uses exaggerated scale, symbolic color, or nonphysical motion, disclose that convention.

The visual standard should be comparable in ambition and refinement to Anatomy, but p-simtor must develop its own visual language.

## 8. Interactive behavior

A generated simulation should be an interactive scientific instrument, not a looping video.

Depending on the physical model, users should be able to:

- Play, pause, restart, and step through time.
- Adjust physically meaningful parameters.
- Change initial conditions.
- Modify mass, velocity, charge, energy, and other relevant quantities.
- Inspect particle positions and trajectories.
- Explore different viewing perspectives.
- Display field quantities and probability distributions.
- Examine equations and numerical values.
- Track scientific observables.
- Compare different initial configurations.
- Export results and simulation data when appropriate.

Changing a physical parameter must update the actual underlying calculation.

For example, increasing one body's mass in a three-body gravitational simulation must produce the corresponding mathematical change to the system's evolution.

Do not merely change the animation speed or visual appearance.

## 9. Generated output

p-simtor must generate complete, usable artifacts.

Depending on the user's request, these may include:

- Self-contained interactive HTML simulations.
- Interactive 3D scientific visualizations.
- Reusable web components.
- Scientific animations.
- High-resolution static visualizations.
- Mathematical plots.
- Source code for the actual physical calculation.
- Exportable numerical data.

The default experience for an interactive request should be an easily opened, usable simulation with an excellent visual presentation.

The agent must be able to generate these artifacts in the user's current project.

Do not require the user to manually assemble scattered code snippets.

Do not require an unnecessarily complicated development environment for ordinary simulations.

Select the necessary implementation technologies through independent research and technical evaluation.

## 10. Skill intelligence

The value of p-simtor is not simply that it contains simulation code.

It should teach the AI agent **how to reason about, implement, evaluate, and visualize a physical system**.

The skill must recognize when a user's request:

- Has a known analytical solution.
- Requires numerical integration.
- Requires statistical or stochastic treatment.
- Involves quantum probability rather than classical trajectories.
- Requires a specialized simulation package.
- Cannot be calculated accurately within reasonable resources.
- Would be misleading if represented as ordinary particle motion.

When possible, select scientifically defensible methods autonomously.

When the requested fidelity is not feasible, explain the limitations and offer the best scientifically defensible alternative rather than generating false physics.

Keep the skill general enough to support new scientific problems beyond its original examples.

## 11. Original architecture

Design p-simtor's implementation from its requirements.

Do not reproduce Anatomy's folders, file hierarchy, drawing tools, camera system, or architectural conventions.

Investigate the best implementation approach for:

- An installable AI skill.
- Reusable computational physics functionality.
- Mathematical solvers and scientific verification.
- Browser-based scientific visualization.
- Advanced quantum or field-based representations.
- Efficient rendering.
- Scientific data formats.
- Automated tests.
- Portable generated artifacts.

Determine whether reusable libraries, generated code, numerical backends, shaders, reference implementations, or combinations of these provide the strongest foundation.

Choose programming languages and dependencies based on evidence rather than preference.

Explain important choices and tradeoffs.

Keep the finished skill usable and maintainable.

Do not introduce unnecessary frameworks, services, or infrastructure.

## 12. Build a real working first version

Do not stop after writing documentation or designing the architecture.

Implement the skill and demonstrate that it works.

Start with a gravitational three-body simulation as the first substantial end-to-end test.

It must include:

- Actual numerical integration of gravitational equations.
- Configurable physical parameters.
- Meaningful initial conditions.
- Accurate trajectories.
- Smooth interactive visualization.
- A suitable integrator.
- Energy and momentum diagnostics.
- Reproducible simulation results.
- Numerical validation.
- A polished visual interface.
- A runnable output artifact generated through the installed skill.

This first implementation is a verification target, not a restriction that p-simtor may only generate gravitational simulations.

Design the skill's instructions and supporting capabilities so they can grow into molecular, quantum, atomic, nuclear, and particle physics without scientifically invalid generalizations.

## 13. Installation, invocation, and regression testing

Test the complete user experience.

The final skill should support this workflow:

**Installation**

The user installs p-simtor directly from GitHub using a supported agent skill installer.

**Invocation**

The user invokes p-simtor through their coding agent.

**Generation**

The agent follows p-simtor's scientific and visual workflow to generate a simulation.

**Execution**

The resulting simulation runs and responds to user input.

**Validation**

The physical system behaves according to its documented mathematical model, within measured numerical tolerances.

**Delivery**

The user receives the generated files, usage instructions, and scientific information.

Include agent-skill evaluations covering requests such as:

- "Simulate a chaotic three-body gravitational system."
- "Create an interactive Lennard-Jones particle simulation."
- "Visualize the probability density of a hydrogen 2p orbital."
- "Simulate a proton moving through a uniform magnetic field."
- "Explain quark confinement with an interactive, scientifically qualified visualization."

Where functionality is not yet implemented, document the limitation explicitly. Do not substitute fake scientific demonstrations or silently advertise unsupported capabilities.

## 14. Success criteria

The project is successful only when:

- p-simtor is a real, installable Agent Skill.
- The repository can be installed directly from GitHub.
- The skill is correctly discovered and invoked by a compatible coding agent.
- It generates working artifacts rather than only textual explanations.
- Its numerical output comes from authentic mathematical models.
- Simulations can be interacted with and modified.
- It includes scientific verification and clear model limitations.
- The visual output is exceptionally polished.
- Its architecture is original and appropriate for computational physics.
- It supports extension into additional physics domains.
- Installation and at least one full generation workflow have been tested.

## 15. Execute the task

Begin by analyzing the Anatomy repository specifically as an example of a successful installable AI agent skill.

Identify what makes its user experience, installation, instructions, generated artifacts, and quality assurance effective.

Then independently research the scientific and computational requirements of p-simtor.

Design an original implementation that delivers an equivalent standard of skill usability while specializing in physics simulation.

**Create the repository, write the skill, implement the necessary functionality, produce the first working simulation, and test the installation and invocation process.**

Work iteratively. Prioritize complete, verified functionality over ambitious but nonfunctional scaffolding.

Do not ask the user to choose a directory structure, programming language, or framework unless a decision genuinely cannot be made from technical evidence.

Do not finish with a design proposal when implementation is possible.

**Project name: p-simtor**

**Product type: Installable AI Agent Skill**

**Reference product: https://github.com/wheresryan22/anatomy**

**Primary purpose: Scientifically rigorous, mathematically driven, visually extraordinary physics simulations.**

**Deliverable: A functional, tested, installable skill — not merely a simulation application.**