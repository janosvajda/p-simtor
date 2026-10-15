// Stochastic dynamics and equilibrium sampling, independent of subject. All randomness comes
// from a seeded generator (core/numerics.mjs rng), so every run is reproducible.
//
//   gillespie        exact stochastic simulation of discrete events (reactions, decays, births/deaths)
//   eulerMaruyama    dX = f(X,t) dt + g(X,t) dW    (Itô SDE, strong order ½)
//   langevinBAOAB    m dv = F dt − γ m v dt + √(2γ m kT) dW   (canonical sampling; Leimkuhler–Matthews)
//   metropolis       samples any target density known up to a constant
//   metropolisSites  single-site updates of a discrete state (spins, occupations, lattice models)
//   autocorrelationTime / effectiveSampleSize   honest error bars for correlated samples

import { rng as makeRng } from './numerics.mjs';

const needRand = (rand, who) => { if (typeof rand !== 'function') throw new Error(`${who}: pass rand = rng(seed) from core/numerics.mjs`); };

/**
 * Gillespie's direct method. state: array of counts; reactions: [{ rate: x => propensity, change: [Δ...] }].
 * Each step fires one event at an exponentially distributed waiting time (exact, no time step).
 * Returns { t, x, step(), advanceTo(t, maxEvents), events }.
 */
export function gillespie({ x0, reactions, rand, t0 = 0 }) {
  needRand(rand, 'gillespie');
  if (!reactions?.length) throw new Error('gillespie: reactions are required');
  const x = Array.from(x0);
  const s = {
    t: t0, x, events: 0, extinct: false,
    step() {
      const a = reactions.map(r => Math.max(0, r.rate(x)));
      const total = a.reduce((p, q) => p + q, 0);
      if (!(total > 0)) { s.extinct = true; s.t = Infinity; return false; }
      s.t += -Math.log(1 - rand()) / total;
      let pick = rand() * total, k = 0;
      while (pick > a[k] && k < a.length - 1) pick -= a[k++];
      reactions[k].change.forEach((d, i) => { x[i] += d; });
      s.events++;
      return true;
    },
    /** Advance to time t: the state at t is the state after the last event before t. */
    advanceTo(target, maxEvents = Infinity) {
      let n = 0;
      for (;;) {
        if (n++ >= maxEvents) return false;
        const save = x.slice(), tSave = s.t;
        if (!s.step()) return true;
        if (s.t > target) { x.splice(0, x.length, ...save); s.events--; s.t = tSave; return true; }
      }
    },
  };
  return s;
}

/** Euler–Maruyama for the Itô SDE dX = f dt + g dW (diagonal noise). f, g return arrays. */
export function eulerMaruyama({ f, g, x0, dt, rand, t0 = 0 }) {
  needRand(rand, 'eulerMaruyama');
  if (!(dt > 0)) throw new Error('eulerMaruyama: dt must be positive');
  const x = Float64Array.from(x0), sq = Math.sqrt(dt);
  const s = {
    t: t0, x, steps: 0,
    step() {
      const drift = f(x, s.t), noise = g(x, s.t);
      for (let i = 0; i < x.length; i++) x[i] += drift[i] * dt + noise[i] * sq * rand.normal();
      s.t += dt; s.steps++;
    },
    advanceTo(target) { while (s.t < target - 1e-12 * dt) s.step(); },
  };
  return s;
}

/**
 * Langevin dynamics with the BAOAB splitting (Leimkuhler & Matthews 2013): samples the canonical
 * distribution with configurational error O(dt²) and is exact for harmonic potentials in x.
 * force(x, out) writes the force; mass, gamma (friction rate), kT.
 */
export function langevinBAOAB({ x0, v0, force, mass = 1, gamma, kT, dt, rand, t0 = 0 }) {
  needRand(rand, 'langevinBAOAB');
  if (!(gamma >= 0 && kT >= 0 && dt > 0 && mass > 0)) throw new Error('langevinBAOAB: need gamma ≥ 0, kT ≥ 0, dt > 0, mass > 0');
  const n = x0.length, x = Float64Array.from(x0), v = Float64Array.from(v0 ?? new Float64Array(n)), F = new Float64Array(n);
  const c1 = Math.exp(-gamma * dt), c2 = Math.sqrt((1 - c1 * c1) * kT / mass);
  force(x, F);
  const s = {
    t: t0, x, v, steps: 0,
    /** Energy given to the system by the bath (friction + noise) since t0: H(t) − H(t0) − heat = O(dt²). */
    heat: 0,
    step() {
      for (let i = 0; i < n; i++) v[i] += 0.5 * dt * F[i] / mass;            // B
      for (let i = 0; i < n; i++) x[i] += 0.5 * dt * v[i];                    // A
      for (let i = 0; i < n; i++) {                                           // O (exact OU step)
        const old = v[i];
        v[i] = c1 * old + c2 * rand.normal();
        s.heat += 0.5 * mass * (v[i] * v[i] - old * old);
      }
      for (let i = 0; i < n; i++) x[i] += 0.5 * dt * v[i];                    // A
      force(x, F);
      for (let i = 0; i < n; i++) v[i] += 0.5 * dt * F[i] / mass;            // B
      s.t += dt; s.steps++;
    },
    advanceTo(target) { while (s.t < target - 1e-12 * dt) s.step(); },
  };
  return s;
}

/**
 * Random-walk Metropolis for a target with log-density logP (up to a constant). propose(x, rand)
 * must be symmetric (default: Gaussian steps of size `stepSize` in every coordinate).
 * Returns { x, logp, step(), accepted, steps, acceptance() }.
 */
export function metropolis({ logP, x0, rand, stepSize = 1, propose }) {
  needRand(rand, 'metropolis');
  const x = Float64Array.from(x0);
  let lp = logP(x);
  if (!Number.isFinite(lp)) throw new Error('metropolis: the starting point must have finite log-probability');
  const y = new Float64Array(x.length);
  const s = {
    x, steps: 0, accepted: 0,
    get logp() { return lp; },
    step() {
      if (propose) y.set(propose(x, rand)); else for (let i = 0; i < x.length; i++) y[i] = x[i] + stepSize * rand.normal();
      const ly = logP(y);
      if (Math.log(1 - rand()) < ly - lp) { x.set(y); lp = ly; s.accepted++; }
      s.steps++;
    },
    acceptance: () => (s.steps ? s.accepted / s.steps : 0),
  };
  return s;
}

/**
 * Metropolis with single-site updates of a discrete state (Ising/Potts spins, lattice occupations,
 * any finite-valued sites). propose(i, state, rand) returns a new value for site i;
 * deltaLogP(i, value, state) returns ln P(new) − ln P(old) for changing only site i
 * (for a Boltzmann target: −ΔE/kT). Proposals must be symmetric. One step() is a sweep of
 * n proposals, sites visited in order ('sequential') or at random ('random').
 * Returns { state, step(), sweeps, accepted, proposals, acceptance() }.
 */
export function metropolisSites({ state, propose, deltaLogP, rand, order = 'sequential' }) {
  needRand(rand, 'metropolisSites');
  if (!state?.length || typeof propose !== 'function' || typeof deltaLogP !== 'function') throw new Error('metropolisSites: need state, propose and deltaLogP');
  if (!['sequential', 'random'].includes(order)) throw new Error('metropolisSites: order must be sequential or random');
  const n = state.length;
  const s = {
    state, sweeps: 0, accepted: 0, proposals: 0,
    step() {
      for (let k = 0; k < n; k++) {
        const i = order === 'random' ? Math.floor(rand() * n) : k;
        const value = propose(i, state, rand);
        const d = deltaLogP(i, value, state);
        if (d >= 0 || Math.log(1 - rand()) < d) { state[i] = value; s.accepted++; }
        s.proposals++;
      }
      s.sweeps++;
    },
    acceptance: () => (s.proposals ? s.accepted / s.proposals : 0),
  };
  return s;
}

/**
 * Integrated autocorrelation time τ of a scalar series (Sokal's self-consistent window, c = 5):
 * the number of steps between effectively independent samples (τ = 1 for independent draws).
 */
export function autocorrelationTime(xs, c = 5) {
  const n = xs.length;
  if (n < 10) throw new Error('autocorrelationTime: need at least 10 samples');
  const mean = xs.reduce((a, b) => a + b, 0) / n;
  const d = Float64Array.from(xs, v => v - mean);
  const c0 = d.reduce((s, v) => s + v * v, 0) / n;
  if (c0 === 0) return 1;
  let tau = 1;
  for (let lag = 1; lag < n / 2; lag++) {
    let ck = 0;
    for (let i = 0; i + lag < n; i++) ck += d[i] * d[i + lag];
    tau += 2 * (ck / (n - lag)) / c0;
    if (lag >= c * tau) break;
  }
  return Math.max(1, tau);
}

/** Effective number of independent samples n/τ, and the standard error of the mean using it. */
export function effectiveSampleSize(xs) {
  const tau = autocorrelationTime(xs), n = xs.length;
  const mean = xs.reduce((a, b) => a + b, 0) / n;
  const variance = xs.reduce((s, v) => s + (v - mean) ** 2, 0) / (n - 1);
  const ess = n / tau;
  return { tau, ess, mean, error: Math.sqrt(variance / ess) };
}

/**
 * Stationary distribution of a continuous-time Markov jump process with rates[i][j] (i → j, i ≠ j;
 * diagonals ignored) by the Grassmann–Taksar–Heyman state reduction: no subtractions, so populations
 * spanning 10⁻³⁰ … 1 (thermal tautomers, rare states) keep full relative accuracy, where solving
 * πQ = 0 directly loses them. Irreducible chains only. Returns a Float64Array summing to 1.
 */
export function stationaryDistribution(rates) {
  const n = rates.length, Q = rates.map(r => Float64Array.from(r));
  for (let i = 0; i < n; i++) { Q[i][i] = 0; for (let j = 0; j < n; j++) if (!(Q[i][j] >= 0)) throw new Error('stationaryDistribution: rates must be ≥ 0'); }
  for (let k = n - 1; k > 0; k--) {
    let s = 0; for (let j = 0; j < k; j++) s += Q[k][j];
    if (!(s > 0)) throw new Error('stationaryDistribution: the chain is not irreducible');
    for (let i = 0; i < k; i++) { const f = Q[i][k] / s; if (f) for (let j = 0; j < k; j++) if (j !== i) Q[i][j] += f * Q[k][j]; }
  }
  const pi = new Float64Array(n); pi[0] = 1;
  for (let k = 1; k < n; k++) {
    let num = 0, den = 0;
    for (let i = 0; i < k; i++) num += pi[i] * Q[i][k];
    for (let j = 0; j < k; j++) den += Q[k][j];
    pi[k] = num / den;
  }
  let sum = 0; for (const v of pi) sum += v;
  return pi.map(v => v / sum);
}

export { makeRng as rng };
