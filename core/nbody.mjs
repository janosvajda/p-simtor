// Gravitational N-body dynamics: forces, invariants and integrators.
//
// State is stored in flat Float64Arrays (x0,y0,z0,x1,...) so the same code runs
// unchanged in Node (tests, verification) and in the browser (generated pages).
// Everything here is pure computation: no DOM, no timers, no randomness.
//
// Units are whatever the caller chooses, as long as they are consistent with G.
// The three-body template uses Hénon N-body units (G = 1).

/** Create a system. bodies: [{ m, x:[x,y,z], v:[vx,vy,vz] }] */
export function createSystem(bodies, { G = 1, softening = 0 } = {}) {
  const n = bodies.length;
  const sys = {
    n, G,
    eps2: softening * softening,
    t: 0,
    m: new Float64Array(n),
    x: new Float64Array(3 * n),
    v: new Float64Array(3 * n),
    a: new Float64Array(3 * n), // scratch, valid after computeAccelerations
    steps: 0,
  };
  bodies.forEach((b, i) => {
    sys.m[i] = b.m;
    for (let k = 0; k < 3; k++) {
      sys.x[3 * i + k] = b.x[k] ?? 0;
      sys.v[3 * i + k] = b.v[k] ?? 0;
    }
  });
  return sys;
}

export function cloneSystem(s) {
  return {
    ...s,
    m: s.m.slice(), x: s.x.slice(), v: s.v.slice(), a: s.a.slice(),
  };
}

/**
 * Pairwise Newtonian (optionally Plummer-softened) accelerations into s.a.
 * Returns U = +sum G m_i m_j / r_ij (the positive "force function"; potential energy is -U).
 */
export function computeAccelerations(s) {
  const { n, x, m, a, G, eps2 } = s;
  a.fill(0);
  let U = 0;
  for (let i = 0; i < n; i++) {
    const xi = x[3 * i], yi = x[3 * i + 1], zi = x[3 * i + 2];
    for (let j = i + 1; j < n; j++) {
      const dx = x[3 * j] - xi, dy = x[3 * j + 1] - yi, dz = x[3 * j + 2] - zi;
      const r2 = dx * dx + dy * dy + dz * dz + eps2;
      const r = Math.sqrt(r2);
      const inv3 = G / (r2 * r);
      a[3 * i] += m[j] * dx * inv3; a[3 * i + 1] += m[j] * dy * inv3; a[3 * i + 2] += m[j] * dz * inv3;
      a[3 * j] -= m[i] * dx * inv3; a[3 * j + 1] -= m[i] * dy * inv3; a[3 * j + 2] -= m[i] * dz * inv3;
      U += G * m[i] * m[j] / r;
    }
  }
  return U;
}

export function kineticEnergy(s) {
  let T = 0;
  for (let i = 0; i < s.n; i++) {
    const vx = s.v[3 * i], vy = s.v[3 * i + 1], vz = s.v[3 * i + 2];
    T += 0.5 * s.m[i] * (vx * vx + vy * vy + vz * vz);
  }
  return T;
}

export function forceFunction(s) {
  const { n, x, m, G, eps2 } = s;
  let U = 0;
  for (let i = 0; i < n; i++) {for (let j = i + 1; j < n; j++) {
    const dx = x[3 * j] - x[3 * i], dy = x[3 * j + 1] - x[3 * i + 1], dz = x[3 * j + 2] - x[3 * i + 2];
    U += G * m[i] * m[j] / Math.sqrt(dx * dx + dy * dy + dz * dz + eps2);
  }}
  return U;
}

export function totalEnergy(s) { return kineticEnergy(s) - forceFunction(s); }

export function totalMomentum(s) {
  const p = [0, 0, 0];
  for (let i = 0; i < s.n; i++) for (let k = 0; k < 3; k++) p[k] += s.m[i] * s.v[3 * i + k];
  return p;
}

export function angularMomentum(s) {
  const L = [0, 0, 0];
  for (let i = 0; i < s.n; i++) {
    const [x, y, z] = [s.x[3 * i], s.x[3 * i + 1], s.x[3 * i + 2]];
    const m = s.m[i];
    const [vx, vy, vz] = [s.v[3 * i], s.v[3 * i + 1], s.v[3 * i + 2]];
    L[0] += m * (y * vz - z * vy);
    L[1] += m * (z * vx - x * vz);
    L[2] += m * (x * vy - y * vx);
  }
  return L;
}

export function centerOfMass(s) {
  const c = [0, 0, 0], u = [0, 0, 0];
  let M = 0;
  for (let i = 0; i < s.n; i++) {
    M += s.m[i];
    for (let k = 0; k < 3; k++) { c[k] += s.m[i] * s.x[3 * i + k]; u[k] += s.m[i] * s.v[3 * i + k]; }
  }
  return { M, x: c.map(q => q / M), v: u.map(q => q / M) };
}

/** Shift to the barycentric frame: centre of mass at rest at the origin. */
export function toBarycentric(s) {
  const { x: c, v: u } = centerOfMass(s);
  for (let i = 0; i < s.n; i++) {for (let k = 0; k < 3; k++) {
    s.x[3 * i + k] -= c[k];
    s.v[3 * i + k] -= u[k];
  }}
  return s;
}

/** Snapshot of conserved quantities, used as the reference for drift diagnostics. */
export function invariants(s) {
  const T = kineticEnergy(s), U = forceFunction(s);
  return { t: s.t, E: T - U, T, U, P: totalMomentum(s), L: angularMomentum(s) };
}

/**
 * Drift of conserved quantities relative to a reference snapshot.
 * Energy is relative; momentum and angular momentum are normalised by natural
 * scales of the system, because their reference values are often exactly zero.
 */
export function drift(s, ref) {
  const now = invariants(s);
  const scaleP = momentumScale(s);
  const scaleL = angularMomentumScale(s);
  const dP = Math.hypot(now.P[0] - ref.P[0], now.P[1] - ref.P[1], now.P[2] - ref.P[2]);
  const dL = Math.hypot(now.L[0] - ref.L[0], now.L[1] - ref.L[1], now.L[2] - ref.L[2]);
  return {
    // relative to |E₀|, or to T₀ + U₀ when the system is (nearly) marginally bound
    dE: Math.abs(now.E - ref.E) / (Math.abs(ref.E) > 1e-9 * (ref.T + ref.U) ? Math.abs(ref.E) : (ref.T + ref.U)),
    dP: dP / scaleP,
    dL: dL / scaleL,
    now,
  };
}

function momentumScale(s) {
  let p = 0;
  for (let i = 0; i < s.n; i++) p += s.m[i] * Math.hypot(s.v[3 * i], s.v[3 * i + 1], s.v[3 * i + 2]);
  return p || 1;
}

function angularMomentumScale(s) {
  let l = 0;
  for (let i = 0; i < s.n; i++) {
    l += s.m[i] * Math.hypot(s.x[3 * i], s.x[3 * i + 1], s.x[3 * i + 2])
      * Math.hypot(s.v[3 * i], s.v[3 * i + 1], s.v[3 * i + 2]);
  }
  return l || 1;
}

/** Crossing (dynamical) time for a bound system: G M^{5/2} / (2|E|)^{3/2}. */
export function crossingTime(s) {
  const M = s.m.reduce((p, q) => p + q, 0);
  const E = totalEnergy(s);
  return s.G * Math.pow(M, 2.5) / Math.pow(2 * Math.abs(E), 1.5);
}

// ---------------------------------------------------------------------------
// Integrators
//
// Every integrator is a function (s, h, ctx) -> void that advances s by one step.
// For fixed-step methods h is a physical time step. For the logarithmic-
// Hamiltonian methods h is a step in fictitious time s, and the physical step
// adapts automatically (dt ≈ h / U), shrinking during close encounters.
// ---------------------------------------------------------------------------

function kick(s, dt) {
  const { v, a } = s;
  for (let i = 0; i < v.length; i++) v[i] += a[i] * dt;
}

function driftPositions(s, dt) {
  const { x, v } = s;
  for (let i = 0; i < x.length; i++) x[i] += v[i] * dt;
}

/** Kick-drift-kick leapfrog (velocity Verlet). 2nd order, symplectic, time-reversible. */
export function leapfrogStep(s, h) {
  computeAccelerations(s);
  kick(s, h / 2);
  driftPositions(s, h);
  computeAccelerations(s);
  kick(s, h / 2);
  s.t += h;
  s.steps++;
}

// Yoshida (1990) 4th-order triple-jump composition coefficients.
const CBRT2 = Math.cbrt(2);
const Y_W1 = 1 / (2 - CBRT2);
const Y_W0 = -CBRT2 / (2 - CBRT2);
const YOSHIDA4 = [Y_W1, Y_W0, Y_W1];

/** Yoshida 4th-order composition of leapfrog. Symplectic, fixed step. */
export function yoshida4Step(s, h) {
  for (const w of YOSHIDA4) {
    leapfrogStep(s, w * h);
    s.steps--;
  }
  s.steps++;
}

/**
 * Logarithmic-Hamiltonian (LogH) leapfrog, Mikkola & Tanikawa (1999),
 * Preto & Tremaine (1999). Integrates Γ = ln(T + B) − ln(U) in fictitious time,
 * where B = −E₀ is the (conserved) binding energy. Drift-kick-drift form:
 *   drift: dt = h/2 / (T + B);  x += v dt;  t += dt
 *   kick:  dt = h / U(x);       v += a dt
 * Exact Kepler orbits up to a phase error, and well-behaved through close
 * encounters without explicit regularisation of coordinates.
 */
export function loghStep(s, h, ctx) {
  const B = ctx.B;
  let dt = (h / 2) / (kineticEnergy(s) + B);
  driftPositions(s, dt);
  s.t += dt;
  const U = computeAccelerations(s);
  kick(s, h / U);
  dt = (h / 2) / (kineticEnergy(s) + B);
  driftPositions(s, dt);
  s.t += dt;
  s.steps++;
}

/** 4th-order Yoshida composition of the (symmetric) LogH leapfrog. Default for few-body problems. */
export function logh4Step(s, h, ctx) {
  for (const w of YOSHIDA4) {
    loghStep(s, w * h, ctx);
    s.steps--;
  }
  s.steps++;
}

export const INTEGRATORS = {
  leapfrog: { step: leapfrogStep, order: 2, adaptive: false, symplectic: true,
    label: 'Leapfrog (KDK)', note: 'Fixed step, 2nd order, symplectic. Fails in close encounters unless dt is tiny.' },
  yoshida4: { step: yoshida4Step, order: 4, adaptive: false, symplectic: true,
    label: 'Yoshida 4', note: 'Fixed step, 4th order, symplectic. Good for softened many-body systems.' },
  logh: { step: loghStep, order: 2, adaptive: true, symplectic: true,
    label: 'LogH leapfrog', note: 'Time-transformed, 2nd order. Physical step shrinks automatically near close approaches.' },
  logh4: { step: logh4Step, order: 4, adaptive: true, symplectic: true,
    label: 'LogH + Yoshida 4', note: 'Time-transformed 4th-order composition. Recommended for unsoftened few-body dynamics.' },
};

/**
 * A stepper wraps a system and an integrator and keeps the previous step's
 * state so positions can be evaluated at any physical time between steps
 * (cubic Hermite interpolation, error O(dt^4)). This is what makes the
 * simulation independent of display frame rate: the step sequence depends only
 * on the initial value problem and h, never on when frames are drawn.
 */
export function createStepper(sys, { integrator = 'logh4', h: hGiven, stepsPerCrossing = 2000 } = {}) {
  const spec = INTEGRATORS[integrator];
  if (!spec) throw new Error(`Unknown integrator "${integrator}"`);
  const ctx = { B: -totalEnergy(sys) };
  const h = hGiven ?? defaultStep(sys, integrator, stepsPerCrossing);
  if (!(Number.isFinite(h) && h !== 0)) throw new Error(`createStepper: step ${h} is not usable (coincident bodies or zero energy scale?)`);
  const prev = { t: sys.t, x: sys.x.slice(), v: sys.v.slice() };
  const ref = invariants(sys);
  return {
    sys, integrator, spec, h, ctx, prev, ref,
    step() {
      prev.t = sys.t; prev.x.set(sys.x); prev.v.set(sys.v);
      const steps0 = sys.steps;
      spec.step(sys, this.h, ctx);
      if (!Number.isFinite(sys.t) || (sys.t - prev.t) * Math.sign(this.h) <= 0
        || !sys.x.every(Number.isFinite) || !sys.v.every(Number.isFinite) || !Number.isFinite(totalEnergy(sys))) {
        sys.t = prev.t; sys.x.set(prev.x); sys.v.set(prev.v); sys.steps = steps0;
        computeAccelerations(sys);
        throw new Error('createStepper: numerical step failed; reduce the step size or resolve the collision');
      }
    },
    /** Advance until sys.t >= tTarget. Returns false if maxSteps was hit first. */
    advanceTo(tTarget, maxSteps = Infinity) {
      let k = 0;
      while (sys.t < tTarget) {
        if (k++ >= maxSteps) return false;
        this.step();
      }
      return true;
    },
    /** Positions (and velocities) at physical time t, with prev.t <= t <= sys.t. */
    sampleAt(t, outX, outV) {
      hermite(prev.t, prev.x, prev.v, sys.t, sys.x, sys.v, t, outX, outV);
    },
  };
}

/**
 * Shortest pairwise dynamical time min √(r_ij³ / G(m_i + m_j)). Unlike the crossing
 * time, it stays meaningful for extreme mass ratios (a planet around a star), where
 * the total energy, and hence the crossing time, is dominated by the light body.
 */
export function pairTime(s) {
  let t = Infinity;
  for (let i = 0; i < s.n; i++) {for (let j = i + 1; j < s.n; j++) {
    const r = Math.hypot(s.x[3 * j] - s.x[3 * i], s.x[3 * j + 1] - s.x[3 * i + 1], s.x[3 * j + 2] - s.x[3 * i + 2]);
    t = Math.min(t, Math.sqrt(r * r * r / (s.G * (s.m[i] + s.m[j]))));
  }}
  return t;
}

/** Shortest pair crossing time r_ij / |v_i − v_j|; resolves fast unbound flybys. */
export function pairCrossingTime(s) {
  let t = Infinity;
  for (let i = 0; i < s.n; i++) {for (let j = i + 1; j < s.n; j++) {
    const speed = Math.hypot(s.v[3 * j] - s.v[3 * i], s.v[3 * j + 1] - s.v[3 * i + 1], s.v[3 * j + 2] - s.v[3 * i + 2]);
    if (speed > 0) t = Math.min(t, Math.hypot(s.x[3 * j] - s.x[3 * i], s.x[3 * j + 1] - s.x[3 * i + 1], s.x[3 * j + 2] - s.x[3 * i + 2]) / speed);
  }}
  return t;
}

/**
 * A sensible default step: dt = t_char / stepsPerCrossing with t_char = min(crossing
 * time, 2π × shortest pair orbital or crossing time). LogH methods: h chosen so that dt ≈ that value at the
 * initial configuration; because dt = h/U adapts, close encounters get shorter steps.
 */
export function defaultStep(sys, integrator, stepsPerCrossing = 2000) {
  const tc = Math.min(crossingTime(sys), 2 * Math.PI * pairTime(sys), 2 * Math.PI * pairCrossingTime(sys));
  const dt = tc / stepsPerCrossing;
  if (!INTEGRATORS[integrator].adaptive) return dt;
  return dt * forceFunction(sys);
}

/** Cubic Hermite interpolation of positions/velocities between two step states. */
export function hermite(t0, x0, v0, t1, x1, v1, t, outX, outV) {
  const dt = t1 - t0;
  if (!(dt > 0)) {
    outX.set(x1);
    if (outV) outV.set(v1);
    return;
  }
  const u = Math.min(1, Math.max(0, (t - t0) / dt));
  const u2 = u * u, u3 = u2 * u;
  const h00 = 2 * u3 - 3 * u2 + 1, h10 = u3 - 2 * u2 + u, h01 = -2 * u3 + 3 * u2, h11 = u3 - u2;
  for (let i = 0; i < x0.length; i++) {
    outX[i] = h00 * x0[i] + h10 * dt * v0[i] + h01 * x1[i] + h11 * dt * v1[i];
  }
  if (outV) {
    const d00 = (6 * u2 - 6 * u) / dt, d10 = 3 * u2 - 4 * u + 1, d01 = (-6 * u2 + 6 * u) / dt, d11 = 3 * u2 - 2 * u;
    for (let i = 0; i < x0.length; i++) {
      outV[i] = d00 * x0[i] + d10 * v0[i] + d01 * x1[i] + d11 * v1[i];
    }
  }
}

/**
 * Integrate a copy of the system to time tEnd and return final state and drift.
 * Used by verification scripts and by in-page convergence checks. With
 * trackConservation, maxDrift records the worst error at every computed step.
 */
export function integrateTo(bodies, opts, tEnd) {
  const sys = createSystem(bodies, opts);
  if (opts.barycentric !== false) toBarycentric(sys);
  const st = createStepper(sys, opts);
  const x = new Float64Array(3 * sys.n), v = new Float64Array(3 * sys.n);
  let maxDrift;
  if (opts.trackConservation) {
    maxDrift = { dE: 0, dP: 0, dL: 0 };
    while (sys.t < tEnd) {
      st.step();
      const errors = drift(sys, st.ref);
      for (const key of ['dE', 'dP', 'dL']) maxDrift[key] = Math.max(maxDrift[key], errors[key]);
    }
  } else st.advanceTo(tEnd);
  st.sampleAt(tEnd, x, v);
  const d = drift(sys, st.ref);
  return { x, v, steps: sys.steps, h: st.h, dE: d.dE, dP: d.dP, dL: d.dL, maxDrift, sys };
}

/** Max position difference between two flat state vectors. */
export function stateDistance(xa, xb) {
  let d = 0;
  for (let i = 0; i < xa.length; i += 3) {
    d = Math.max(d, Math.hypot(xa[i] - xb[i], xa[i + 1] - xb[i + 1], xa[i + 2] - xb[i + 2]));
  }
  return d;
}
