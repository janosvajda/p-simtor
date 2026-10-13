# Changelog

## 0.1.0 — 2026-10-13

First working version.

- Agent Skill packaging (`SKILL.md`, Agent Skills spec-compliant frontmatter), installed with its own `scripts/install.mjs` (Node built-ins only; no dependencies).
- Gravity instrument: Newtonian N-body (≤ 64 bodies) with Leapfrog, Yoshida-4, LogH and LogH+Yoshida-4 integrators, Hermite dense output, frame-rate-independent stepping, five cited three-body presets, custom scenarios, chaos twin, in-page convergence study, CSV/JSON/PNG export.
- Build-time verification: conservation, convergence order, predictability horizon, determinism, frame independence, Kepler / figure-eight / Pythagorean benchmarks.
- Compose path for any domain: `ui/instrument-app.mjs` shell, simulation contract, `verify-module`, bundler with local-module support, Node `psim:` resolver.
- Worked example: relativistic proton in uniform E and B fields (Boris pusher) with six exact checks.
- Headless-Chrome checks over the DevTools protocol: runtime errors, liveness, overflow, screenshots, UI→physics interaction tests.
- References: method selection, verification, gravity, electromagnetism, molecular dynamics, quantum, particle physics, compose contract, visual language.
