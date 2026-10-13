# Gravity: two-, three- and N-body

## Model

Newtonian point masses: r̈ᵢ = Σⱼ G mⱼ (rⱼ − rᵢ)/|rⱼ − rᵢ|³ (optionally Plummer-softened: |r|² → |r|² + ε²). Conserved: E, P, L. Direct O(N²) summation, up to 64 bodies in the template.

Not modelled (say so where relevant): general relativity (perihelion precession of Mercury, compact binaries, anything with v/c ≳ 1e-3), tides, finite body size, collisions, mass loss, radiation pressure.

## The built-in instrument

```bash
node "$PSIM/scripts/psim.mjs" gravity --preset pythagorean --out simulations/three-body.html
node "$PSIM/scripts/psim.mjs" gravity --scenario sim.json --out simulations/x.html [--with-presets]
```

Presets (Hénon units, G = 1): `figure-eight` (Chenciner–Montgomery), `pythagorean` (Burrau; Szebehely & Peters 1967), `lagrange` (equilateral, unstable for equal masses), `hierarchical` (stable triple), `broucke-a2` (Broucke 1975). For a three-body request without specific initial conditions, use the default (`pythagorean`, the classic chaotic problem). The page includes all presets so the user can compare them.

The page gives: play, pause, step, restart and speed; per-body mass, position and velocity editing (any edit re-integrates from t = 0, and bodies can be dragged when paused at t = 0); integrator choice and accuracy; live |ΔE/E₀|, ΔP and ΔL with a plot; a chaos twin with a separation plot and growth rate; the build-time verification report plus an in-page convergence study; MathML equations; CSV, JSON and PNG export; pan, zoom and 3D tilt.

## Scenario JSON

```json
{
  "key": "sun-jupiter-saturn",
  "title": "Sun, Jupiter and Saturn",
  "pageTitle": "Jupiter–Saturn great inequality",
  "subtitle": "Newtonian three-body, AU–M☉–yr units",
  "note": "One-line description shown under the scenario picker",
  "source": "Citation for the initial conditions",
  "units": { "G": 39.47841760435743, "system": "AU, M☉, yr (G = 4π²)", "length": "AU", "time": "yr", "mass": "M☉" },
  "integrator": "logh4",
  "stepsPerCrossing": 2000,
  "softening": 0,
  "duration": 300,
  "speed": 10,
  "trailTime": 60,
  "bodies": [
    { "name": "Sun", "m": 1, "x": [0, 0, 0], "v": [0, 0, 0] },
    { "name": "Jupiter", "m": 9.547919e-4, "x": [5.2026, 0, 0], "v": [0, 2.7546, 0] }
  ],
  "notes": ["Extra assumptions shown in the model panel"],
  "includePresets": false
}
```

- `duration`: physical time verified at build time, and the playback scale (default speed is duration/30 per second).
- `stepsPerCrossing`: accuracy knob. Raise it (8000–32000) for longer predictability horizons in chaotic systems.
- `softening > 0` sets the tier to `approximate` automatically. Explain it in `notes`.
- Units: the template is unit-agnostic. Choose G and labels consistently (G = 4π² for AU, M☉, yr; G = 1 for Hénon units). Circular speed: v = √(G M / r).
- Planetary data: take values from JPL Horizons or a cited ephemeris if the user wants real positions. Never invent "realistic" orbital elements. If idealised (circular, coplanar), say so in `note`.

## Integrator guidance

| Situation | Choose |
|---|---|
| 2–10 bodies, close encounters possible (default) | `logh4` |
| Hierarchical, well-separated (planetary systems) | `logh4`; `yoshida4` also fine |
| Softened clusters, N > 10 | `yoshida4` with softening ≈ 0.01–0.05 × system size |
| Teaching the difference | Show `leapfrog` vs `logh4` energy plots on `pythagorean` |

LogH needs a non-degenerate system: B = −E₀ must keep T + B > 0, which holds automatically since T + B = U > 0.

## Reading the verification report

- Energy error ≲ 1e-10 with `logh4` at defaults; the Pythagorean run gives ~2e-11 through all its encounters.
- Convergence order ≈ 4 (`logh4`, `yoshida4`) or 2 (`logh`, `leapfrog`).
- **Predictability horizon**: chaotic presets (Pythagorean ≈ 59, Lagrange ≈ 22 at defaults) diverge between step sizes after this time. Say this in the reply. The qualitative outcome (ejection, binary formation) is the prediction, not the exact path after the horizon.
- Determinism and frame-rate independence must both be `true`.

## Visual notes

Bodies are point masses. Marker size ∝ m^(1/3) and is disclosed. Trails are integrator states only. The view is orthographic, so in-plane lengths match the scale bar. For real systems with huge dynamic range (Sun–planet distances), don't draw bodies to scale. Disclose the marker convention, and consider a log-radius view only if you label it.
