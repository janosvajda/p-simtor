// Numerical verification for gravitational N-body simulations.
// Pure functions: used by scripts/psim.mjs (baked reports), tests, and in-page checks.

import {
  createSystem, createStepper, toBarycentric, integrateTo, stateDistance, drift, INTEGRATORS,
} from './nbody.mjs';
import { keplerBodies, keplerRelative } from './kepler.mjs';
import { THREE_BODY_PRESETS } from './presets.mjs';

function systemSize(bodies) {
  let r = 0;
  for (const b of bodies) r = Math.max(r, Math.hypot(b.x[0], b.x[1], b.x[2] ?? 0));
  return r || 1;
}

/**
 * Run the same initial value problem at h, h/2, h/4 and compare trajectories at
 * `samples` times. Returns the observed order of accuracy (from the earliest
 * samples, before chaos amplifies differences) and the predictability horizon:
 * the last time at which the h/2 and h/4 runs still agree to `tol` × system size.
 */
export function convergenceStudy(bodies, opts, tEnd, { samples = 200, tol = 1e-3 } = {}) {
  let base;
  const runs = [1, 0.5, 0.25].map(f => {
    const sys = createSystem(bodies, opts);
    toBarycentric(sys);
    const st = createStepper(sys, opts);
    base ??= st.h;
    st.h = base * f;
    return { st, x: new Float64Array(3 * sys.n) };
  });
  const size = systemSize(bodies);
  const rows = [];
  let horizon = tEnd, horizonFound = false;
  for (let k = 1; k <= samples; k++) {
    const t = tEnd * k / samples;
    for (const r of runs) { r.st.advanceTo(t); r.st.sampleAt(t, r.x); }
    const e1 = stateDistance(runs[0].x, runs[1].x);
    const e2 = stateDistance(runs[1].x, runs[2].x);
    rows.push({ t, e1, e2 });
    if (!horizonFound && e2 > tol * size) { horizon = rows.length > 1 ? rows[rows.length - 2].t : 0; horizonFound = true; }
  }
  // Observed order from samples where errors are resolved above round-off and still small.
  const orders = rows
    .filter(r => r.e2 > 1e-13 * size && r.e1 < 1e-2 * size)
    .slice(0, 20)
    .map(r => Math.log2(r.e1 / r.e2));
  const order = orders.length ? median(orders) : null;
  return {
    integrator: opts.integrator,
    expectedOrder: INTEGRATORS[opts.integrator].order,
    observedOrder: order,
    horizon,
    converged: !horizonFound,
    finalDifference: rows[rows.length - 1].e2 / size,
    h: runs[0].st.h,
    steps: runs.map(r => r.st.sys.steps),
    rows,
  };
}

function median(a) {
  const s = [...a].sort((p, q) => p - q);
  return s[Math.floor(s.length / 2)];
}

/** Two independent runs must be bit-identical. */
export function determinismCheck(bodies, opts, tEnd) {
  const a = integrateTo(bodies, opts, tEnd), b = integrateTo(bodies, opts, tEnd);
  let identical = a.steps === b.steps;
  for (let i = 0; i < a.x.length; i++) if (a.x[i] !== b.x[i] || a.v[i] !== b.v[i]) identical = false;
  return { identical, steps: a.steps };
}

/**
 * Advance the same system with two very different "frame" schedules (60 Hz and
 * irregular 7–41 ms frames at different playback speeds). The physical state at
 * a common time must be bit-identical: rendering never alters the step sequence.
 */
export function frameIndependenceCheck(bodies, opts, tEnd) {
  const final = schedule => {
    const sys = createSystem(bodies, opts);
    toBarycentric(sys);
    const st = createStepper(sys, opts);
    let t = 0, k = 0;
    while (t < tEnd) {
      t = Math.min(tEnd, t + schedule(k++));
      st.advanceTo(t);
    }
    const x = new Float64Array(3 * sys.n);
    st.sampleAt(tEnd, x);
    return { x, steps: sys.steps };
  };
  const a = final(() => tEnd / 600);
  const b = final(k => tEnd / 600 * (0.2 + 2.3 * ((k * 7919) % 97) / 97));
  let identical = a.steps === b.steps;
  for (let i = 0; i < a.x.length; i++) if (a.x[i] !== b.x[i]) identical = false;
  return { identical, steps: a.steps };
}

/** Eccentric Kepler orbit vs the analytic solution after `orbits` periods. */
export function keplerBenchmark({ e = 0.9, orbits = 10, integrator = 'logh4', stepsPerCrossing = 2000 } = {}) {
  const m1 = 1, m2 = 0.5, a = 1;
  const bodies = keplerBodies(m1, m2, a, e);
  const { period } = keplerRelative(a, e, m1 + m2, 0);
  const tEnd = period * orbits;
  const sys = createSystem(bodies);
  const st = createStepper(sys, { integrator, h: undefined });
  if (stepsPerCrossing !== 2000) st.h *= 2000 / stepsPerCrossing;
  const x = new Float64Array(6);
  st.advanceTo(tEnd);
  st.sampleAt(tEnd, x);
  const ref = keplerRelative(a, e, m1 + m2, tEnd);
  const rel = [x[3] - x[0], x[4] - x[1]];
  const err = Math.hypot(rel[0] - ref.x, rel[1] - ref.y) / a;
  return { name: `Kepler e=${e}, ${orbits} orbits`, integrator, positionError: err, dE: drift(sys, st.ref).dE, steps: sys.steps };
}

/** Figure-eight: return to initial state after one published period. */
export function figureEightBenchmark({ integrator = 'logh4' } = {}) {
  const p = THREE_BODY_PRESETS['figure-eight'];
  const r = integrateTo(p.bodies, { integrator }, p.period);
  const x0 = Float64Array.from(p.bodies.flatMap(b => b.x));
  return { name: 'Figure-eight periodicity (1 period)', integrator, returnError: stateDistance(r.x, x0), dE: r.dE, steps: r.steps };
}

/** Pythagorean problem: energy through close encounters, and the published outcome. */
export function pythagoreanBenchmark({ integrator = 'logh4' } = {}) {
  const p = THREE_BODY_PRESETS.pythagorean;
  const r = integrateTo(p.bodies, { integrator }, 80);
  const d = (i, j) => Math.hypot(r.x[3 * i] - r.x[3 * j], r.x[3 * i + 1] - r.x[3 * j + 1], r.x[3 * i + 2] - r.x[3 * j + 2]);
  const ejected = d(0, 1) > 20 && d(0, 2) > 20 && d(1, 2) < 2;
  return {
    name: 'Pythagorean problem to t=80',
    integrator, dE: r.dE, steps: r.steps,
    outcome: ejected ? 'm=3 ejected, (4,5) binary' : 'unexpected final configuration',
    matchesLiterature: ejected,
  };
}

/** Everything a scenario report needs. `scenario` = { bodies, integrator, duration, G?, softening? }. */
export function verifyScenario(scenario) {
  const opts = {
    integrator: scenario.integrator ?? 'logh4', G: scenario.G ?? 1, softening: scenario.softening ?? 0,
    stepsPerCrossing: scenario.stepsPerCrossing ?? 2000,
  };
  const tEnd = scenario.duration;
  const t0 = Date.now();
  const main = integrateTo(scenario.bodies, opts, tEnd);
  const conv = convergenceStudy(scenario.bodies, opts, tEnd);
  delete conv.rows;
  return {
    generated: new Date().toISOString(),
    integrator: opts.integrator,
    integratorLabel: INTEGRATORS[opts.integrator].label,
    duration: tEnd,
    conservation: { dE: main.dE, dP: main.dP, dL: main.dL, steps: main.steps, h: main.h },
    convergence: conv,
    determinism: determinismCheck(scenario.bodies, opts, tEnd),
    frameIndependence: frameIndependenceCheck(scenario.bodies, opts, tEnd),
    benchmarks: [keplerBenchmark({ integrator: opts.integrator }), figureEightBenchmark({ integrator: opts.integrator }), pythagoreanBenchmark({ integrator: opts.integrator })],
    runtimeMs: Date.now() - t0,
  };
}
