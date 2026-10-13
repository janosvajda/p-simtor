// Analytic two-body (Kepler) solution, used as a reference for numerical integrators.

/** Solve Kepler's equation E − e sin E = M for elliptic orbits (Newton iteration). */
/**
 * Eccentricities computed from state vectors carry round-off: a circular orbit can come out as
 * e = −1e−16. Treat |e| below 1e−12 as 0; reject only physically impossible values.
 */
export function checkEccentricity(e, who) {
  if (Math.abs(e) < 1e-12) return 0;
  if (!(e > 0 && e < 1)) throw new Error(`${who}: elliptic orbits need 0 ≤ e < 1 (got ${e})`);
  return e;
}

export function solveKepler(meanAnomaly, eIn, tol = 1e-15) {
  const e = checkEccentricity(eIn, 'solveKepler');
  // reduce to [0, 2π): with E₀ = π for high e, Newton converges monotonically only there
  const M = ((meanAnomaly % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  let E = e < 0.8 ? M : Math.PI;
  for (let k = 0; k < 100; k++) {
    const f = E - e * Math.sin(E) - M;
    const dE = f / (1 - e * Math.cos(E));
    E -= dE;
    if (Math.abs(dE) < tol) break;
  }
  return E;
}

/**
 * Relative position of body 2 with respect to body 1 at time t, for an orbit
 * that starts at pericentre on the +x axis and moves counter-clockwise.
 * mu = G (m1 + m2).
 */
export function keplerRelative(a, e, mu, t) {
  const n = Math.sqrt(mu / (a * a * a));
  const E = solveKepler(n * t, e);
  const x = a * (Math.cos(E) - e);
  const y = a * Math.sqrt(1 - e * e) * Math.sin(E);
  const r = a * (1 - e * Math.cos(E));
  const vx = -Math.sqrt(mu * a) / r * Math.sin(E);
  const vy = Math.sqrt(mu * a) / r * Math.sqrt(1 - e * e) * Math.cos(E);
  return { x, y, vx, vy, period: 2 * Math.PI / n };
}

/** Bodies for a barycentric two-body orbit starting at pericentre. */
export function keplerBodies(m1, m2, a, eIn, G = 1) {
  const e = checkEccentricity(eIn, 'keplerBodies');
  const mu = G * (m1 + m2);
  const r = a * (1 - e);
  const v = Math.sqrt(mu * (1 + e) / r);
  const f1 = m2 / (m1 + m2), f2 = m1 / (m1 + m2);
  return [
    { m: m1, x: [-f1 * r, 0, 0], v: [0, -f1 * v, 0] },
    { m: m2, x: [f2 * r, 0, 0], v: [0, f2 * v, 0] },
  ];
}
