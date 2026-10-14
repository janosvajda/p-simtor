# Visual language

p-simtor pages look like **precision instruments**: a dark observatory field (or warm paper in light mode), luminous computed paths, quiet hairline structure and numbers you can trust at a glance. Beauty comes from fidelity: smooth curves because the integrator resolved them, glow because a body is a point source of attention, not decoration.

## The shell (ui/instrument.css)

- Layout: a large stage (the physics) plus a 372 px control column. On phones it stacks: stage at 64vh, controls below. No horizontal scroll (the `check` command fails on overflow).
- Header: title, one-line subtitle naming the model and method, **trust-tier badge** (click opens the model panel), theme toggle.
- HUD (top-left of stage): the labelled progress coordinate (physical time, sampling or scan), the primary diagnostic and resolution. Legend (top-right): units and visual conventions in one or two lines. Scale bar (bottom-left). Transport (bottom-centre).
- Colours are CSS tokens (`--ink`, `--muted`, `--accent`, `--body-0…5`, `--grid`, `--glow`). Both themes are first-class: always read colours through `view.color('--token')` so they switch live.
- Numbers: tabular figures, scientific notation as 2.4×10⁻¹¹ (`sci()`), units always shown.

## Derive the view from the state and observables

Select a representation after the mathematical specification is established. Similar-looking pictures do not imply similar physics.

| Object being represented | Suitable view | Required interpretation |
|---|---|---|
| Physical coordinates or trajectories | Computed positions, paths, phase-space plots | Units, coordinate system, projection, marker-size and time scaling |
| Scalar/vector/tensor field | Labelled scalar maps, arrows, streamlines, components or justified invariants | Quantity, units, component/basis, averaging and scale; arrow geometry must have a defined meaning |
| Amplitude or probability | Components, phase/magnitude, density maps, samples or level sets | Amplitude versus probability, normalization, sampling meaning and transfer function |
| Populations or statistical estimates | Distributions, histories, correlations and estimates with error bars | State space, ensemble, estimator, aggregation and sampling uncertainty |
| Relation or abstract coordinates | Parameter plots, surfaces, diagrams or eigenstate/constraint views | Physical space versus internal/parameter space; evaluated relation versus physical evolution |
| Observed/fitted data | Data points, uncertainty, fitted/interpolated/extrapolated values labelled separately | Provenance, fitted parameters, validity domain and limits of extrapolation |

For every channel—position, height, colour, size, arrows, opacity and playback—record:

- the state variable, computed observable, dataset or symbolic label it comes from;
- the mapping, normalization, units and uncertainty/averaging where applicable;
- whether the geometry depicts physical space, internal coordinates or a diagram;
- which part is a physical prediction and which part is a representation convention.

Use this mapping to write the adjacent legend and `model.conventions`. Colour might encode a numerical value, phase, categorical identity or a purely symbolic label. Claim optical colour only when it is actually derived from appropriate optical data/calculation. A mathematical charge/category label does not imply a visible pigment; an internal coordinate does not imply a direction in physical space. Disclose basis-dependent components and projected quantities instead of presenting representation choices as measured properties.

Label physical time, sample/sweep number, iteration, parameter scans and static views separately. A computed configuration is not necessarily an instantaneous physical snapshot. Label any conceptual overlay so it cannot borrow the computation's credibility.

## When to use the GPU (`ui/gl.mjs`)

Choose the renderer from the size and shape of the state, as a cost estimate (see method selection §5):

- **Points:** up to ~10⁴ moving bodies draw well in 2D canvas; beyond that, or when depth and glow carry meaning (a galaxy, a gas, a sampled probability cloud), use `createPointRenderer` with an `orbitCamera`.
- **2D fields:** above ~200 × 200 cells updated every frame, use `createFieldRenderer` instead of `putImageData`.
- **3D densities:** a scalar on a 3D grid (orbital |ψ|², gas density, a temperature field) is shown honestly as a ray-marched volume (`createVolumeRenderer`), not as a cloud of invented particles. Sampled points are a different object: label which one the view shows.
- Additive glow suits emission on dark backgrounds; in the light theme use `blend: 'alpha'` (read `view.color('--glow')`).
- Brightness, opacity, glow and colour-map range are representation choices: put them in the legend and `model.conventions`. A log colour scale must say so.

## Craft rules

- **Trails**: draw in ~24–28 chunks with alpha ∝ (1 − age)^1.6 and slightly thicker heads. Round joins. Never smooth: the polyline vertices are the integrator's states.
- **Glow** (dark theme only, `--glow = 1`): additive radial gradient (≈6× marker radius) plus a solid core plus a small specular highlight. In light theme, use a solid core with a soft halo ring and a hairline outline.
- **Reference overlays** (analytic solutions) are dashed and muted; computed results are solid and coloured. The eye should see the agreement at once.
- **Grids**: hairline, two weights (minor / every 5th major), and they fade with zoom. A grid is a measuring device: it must use nice steps (1-2-5).
- **Motion**: choose readable playback for the computed dynamics or sampling and disclose its meaning and compression. Respect `prefers-reduced-motion` (start paused).
- **Typography**: sans for prose, mono for numbers, MathML for equations. No emoji.
- **Disclosure**: important conventions and reductions appear beside the view; the model panel contains their detail. Include exaggerated marker size, log scaling, time compression and symbolic colour where used.
- Restraint: no lens flares, bloom on everything, random particles, or starfields behind physics that has no stars. Every pixel either shows computed physics, a measurement aid, or a disclosed convention.

## Before delivering, look at the screenshots and ask

1. Within three seconds, can a scientist tell what quantity each visual element represents?
2. Do computed results visibly agree with any analytic overlay?
3. Is anything clipped, overlapping, unreadable, or empty in either theme or on the phone?
4. Is every number labelled with units, and every convention disclosed?
