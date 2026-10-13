// Analytic two-body (Kepler) solution, used as a reference for numerical integrators.

/** Solve Kepler's equation E − e sin E = M for elliptic orbits (Newton iteration). */
export function solveKepler(M, e, tol = 1e-15) {
  M = M % (2 * Math.PI);
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
export function keplerBodies(m1, m2, a, e, G = 1) {
  const mu = G * (m1 + m2);
  const r = a * (1 - e);
  const v = Math.sqrt(mu * (1 + e) / r);
  const f1 = m2 / (m1 + m2), f2 = m1 / (m1 + m2);
  return [
    { m: m1, x: [-f1 * r, 0, 0], v: [0, -f1 * v, 0] },
    { m: m2, x: [f2 * r, 0, 0], v: [0, f2 * v, 0] },
  ];
}
