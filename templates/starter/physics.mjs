// STARTER (psim new --starter): a minimal, verified model with every part a page needs.
// Replace the model, keep the contract. Delete these STARTER comments when done.
//
// ONE IDEA (SKILL.md step 5; write these three lines first, then design everything around them)
//   Learn:   a radioactive sample halves at a fixed rhythm, the half-life, whatever its size.
//   Picture: the atoms themselves (a random sample on a grid), each going dark when it decays;
//            the count falling against the exact law is evidence, in a chart below the view.
//   Action:  "Halve the half-life": the same fall, twice as fast.
//
// Model: dN/dt = −λN, λ = ln 2 / T½, integrated with classical RK4 at a fixed step.
// Units: time in s, N in atoms. Exact solution N(t) = N₀ 2^(−t/T½).
// Picture: SAMPLE = 400 atoms with seeded random decay times t_i = −ln(u_i)/λ (the same law, one
// realisation; seed 29 is a typical one, within 1 % of the law at 1–3 half-lives, said in the legend):
// the scene shows the thing itself, the smooth N(t) is the averaged evidence.

export const PARAMS = [
  { key: 'halfLife', label: 'Half-life T½', min: 0.5, max: 100, value: 10, log: true, unit: 's', hint: 'Time for half of the atoms to decay.' },
];

const rk4 = (f, y, h) => {
  const k1 = f(y), k2 = f(y + h / 2 * k1), k3 = f(y + h / 2 * k2), k4 = f(y + h * k3);
  return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
};

export const SAMPLE = 400;
// small seeded generator (mulberry32): the same atoms decay at the same moments on every run
const rng = seed0 => { let seed = seed0 | 0; return () => { seed = seed + 0x6D2B79F5 | 0; let x = Math.imul(seed ^ seed >>> 15, 1 | seed); x = x + Math.imul(x ^ x >>> 7, 61 | x) ^ x; return ((x ^ x >>> 14) >>> 0) / 4294967296; }; };

export function createSim({ halfLife, n0 = 1000 }) {
  const lambda = Math.LN2 / halfLife, h = halfLife / 200;
  const r = rng(29), decayAt = Float64Array.from({ length: SAMPLE }, () => -Math.log(1 - r()) / lambda);
  const exact = t => n0 * 2 ** (-t / halfLife);
  let t = 0, N = n0;
  const trace = [[0, n0]];
  return {
    halfLife, n0, exact, decayAt,
    /** Atoms of the drawn sample still intact at time tt. */
    intact: (tt = t) => decayAt.reduce((k, d) => k + (d > tt), 0),
    get t() { return t; },
    get N() { return N; },
    trace,
    step() { N = rk4(y => -lambda * y, N, h); t += h; trace.push([t, N]); },
    diagnostics: () => [
      { key: 'err', label: 'Relative error vs exact N(t)', value: Math.abs(N - exact(t)) / exact(t), kind: 'error', good: 1e-9, warn: 1e-6 },
      { key: 'frac', label: 'Atoms left', value: N / n0, kind: 'value', text: `${(100 * N / n0).toFixed(1)} %` },
    ],
    state: () => ({ t, N }),
    record: () => ({ header: ['t (s)', 'N'], rows: trace.map(([a, b]) => [a, b]) }),
  };
}

export function verify() {
  // against the exact solution, and the measured RK4 order (expected 4) by step halving
  const T = 10, run = h => { let N = 1; for (let t = 0; t < 5 * T - 1e-12; t += h) N = rk4(y => -Math.LN2 / T * y, N, h); return N; };
  const exact = 2 ** -5, e1 = Math.abs(run(T / 20) - exact), e2 = Math.abs(run(T / 40) - exact);
  const order = Math.log2(e1 / e2);
  const sim = createSim({ halfLife: T, n0: 1000 });
  while (sim.t < T - 1e-9) sim.step();
  const half = Math.abs(sim.N / 500 - 1);
  // the drawn sample follows the same law: intact fraction after one half-life within 3σ of ½
  const left = sim.intact(T) / SAMPLE, sigma = Math.sqrt(0.25 / SAMPLE), z = Math.abs(left - 0.5) / sigma;
  return {
    summary: `RK4 decay vs the exact half-life law: order ${order.toFixed(2)}, N(T½) within ${half.toExponential(1)} of N₀/2`,
    checks: [
      { label: 'N(T½) = N₀/2 (exact half-life law)', value: half, pass: half < 1e-9 },
      { label: 'Measured RK4 order (expected 4)', value: order, pass: order > 3.7 },
      { label: `Drawn sample of ${SAMPLE} atoms after one half-life: intact fraction vs ½ (in standard deviations)`, value: z, pass: z < 3 },
    ],
    notes: ['Exact solution N₀ 2^(−t/T½); order from step halving over five half-lives.', 'The drawn atoms are one seeded random sample of the same law, not the RK4 count.'],
  };
}
