# Visual language

p-simtor pages look like **precision instruments**: a dark observatory field (or warm paper in light mode), luminous computed paths, quiet hairline structure and numbers you can trust at a glance. Beauty comes from fidelity: smooth curves because the integrator resolved them, glow because a body is a point source of attention, not decoration.

## The shell (ui/instrument.css)

- Layout: a large stage (the physics) plus a 372 px control column. On phones it stacks: stage at 64vh, controls below. No horizontal scroll (the `check` command fails on overflow).
- Header: title, one-line subtitle naming the model and method, **trust-tier badge** (click opens the model panel), theme toggle.
- HUD (top-left of stage): physical time, the primary error metric, step size. Legend (top-right): units and visual conventions in one or two lines. Scale bar (bottom-left). Transport (bottom-centre).
- Colours are CSS tokens (`--ink`, `--muted`, `--accent`, `--body-0…5`, `--grid`, `--glow`). Both themes are first-class: always read colours through `view.color('--token')` so they switch live.
- Numbers: tabular figures, scientific notation as 2.4×10⁻¹¹ (`sci()`), units always shown.

## Choose the representation from the physics

| Phenomenon | Primary visual | Supporting visuals |
|---|---|---|
| Trajectories (gravity, charges) | Fading polylines through computed states; luminous point markers | Energy-error plot, phase-space or separation plot |
| Fields | Field lines (density ∝ strength) or normalised arrows with magnitude colour | Equipotentials, |F| colour map with a colour bar |
| Probability densities | Volumetric glow or 90% isosurface, phase as hue; point cloud for sampling intuition | Radial distribution P(r) with nodes, energy level |
| Wavefunctions in 1D | Re/Im or |ψ| with phase hue, \|ψ\|² filled; potential drawn behind | Norm, ⟨x⟩, T/R over time |
| Many-particle systems | Discs or spheres coloured by speed or species | Live histograms against theory (MB), g(r), T(t) |
| Potentials and coupling (QCD) | Curves with shaded validity regions; parameter sweep slider | Annotated thresholds with citations |

## Craft rules

- **Trails**: draw in ~24–28 chunks with alpha ∝ (1 − age)^1.6 and slightly thicker heads. Round joins. Never smooth: the polyline vertices are the integrator's states.
- **Glow** (dark theme only, `--glow = 1`): additive radial gradient (≈6× marker radius) plus a solid core plus a small specular highlight. In light theme, use a solid core with a soft halo ring and a hairline outline.
- **Reference overlays** (analytic solutions) are dashed and muted; computed results are solid and coloured. The eye should see the agreement at once.
- **Grids**: hairline, two weights (minor / every 5th major), and they fade with zoom. A grid is a measuring device: it must use nice steps (1-2-5).
- **Motion**: playback speed is set so the main timescale (orbit, gyration, oscillation) takes 2–10 s. Respect `prefers-reduced-motion` (start paused).
- **Typography**: sans for prose, mono for numbers, MathML for equations. No emoji.
- **Disclosure**: any exaggeration (marker size, log scaling, time compression, symbolic colour) appears in the legend or model panel.
- Restraint: no lens flares, bloom on everything, random particles, or starfields behind physics that has no stars. Every pixel either shows computed physics, a measurement aid, or a disclosed convention.

## Before delivering, look at the screenshots and ask

1. Within three seconds, can a scientist tell what quantity each visual element represents?
2. Do computed results visibly agree with any analytic overlay?
3. Is anything clipped, overlapping, unreadable, or empty in either theme or on the phone?
4. Is every number labelled with units, and every convention disclosed?
