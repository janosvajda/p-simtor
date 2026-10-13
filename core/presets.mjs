// Initial conditions for gravitational few-body systems, with their sources.
// All presets use Hénon N-body units: G = 1, lengths and masses dimensionless.

export const THREE_BODY_PRESETS = {
  'figure-eight': {
    title: 'Figure-eight choreography',
    source: 'Chenciner & Montgomery (2000), Ann. Math. 152, 881. Initial conditions from Simó (2002).',
    note: 'Three equal masses chase each other along one figure-eight curve. Periodic, and linearly stable.',
    period: 6.32591398,
    duration: 6.32591398 * 4,
    bodies: [
      { m: 1, x: [0.97000436, -0.24308753, 0], v: [0.466203685, 0.43236573, 0] },
      { m: 1, x: [-0.97000436, 0.24308753, 0], v: [0.466203685, 0.43236573, 0] },
      { m: 1, x: [0, 0, 0], v: [-0.93240737, -0.86473146, 0] },
    ],
  },
  'pythagorean': {
    title: 'Pythagorean problem (Burrau)',
    source: 'Burrau (1913); Szebehely & Peters (1967), Astron. J. 72, 876.',
    note: 'Masses 3, 4, 5 start at rest on a 3-4-5 triangle. A long chaotic dance of close encounters ends near t ≈ 60–70 with body 3 ejected and bodies 4 and 5 bound as a binary.',
    duration: 80,
    stepsPerCrossing: 16000, // resolves the energy-error spike during the closest encounter
    bodies: [
      { m: 3, x: [1, 3, 0], v: [0, 0, 0] },
      { m: 4, x: [-2, -1, 0], v: [0, 0, 0] },
      { m: 5, x: [1, -1, 0], v: [0, 0, 0] },
    ],
  },
  'lagrange': {
    title: 'Lagrange equilateral triangle',
    source: 'Lagrange (1772). Linear stability requires (m1+m2+m3)² > 27(m1m2+m2m3+m3m1) (Routh 1875).',
    note: 'Equal masses on a rotating equilateral triangle. An exact solution, but unstable for equal masses: round-off grows until the configuration breaks up.',
    duration: 60,
    stepsPerCrossing: 4000,
    bodies: lagrangeTriangle([1, 1, 1], 1),
  },
  'hierarchical': {
    title: 'Hierarchical triple',
    source: 'Constructed: inner Kepler binary (a = 0.1, e = 0) orbited by a light outer body (a = 1, e = 0.2).',
    note: 'A tight inner binary with a distant third body. Long-term stable; the outer orbit sees the binary almost as a point mass.',
    duration: 40,
    bodies: hierarchicalTriple(),
  },
  'broucke-a2': {
    title: 'Broucke–Hénon orbit',
    source: 'Broucke (1975), Celest. Mech. 12, 439; family A, orbit A2 (as tabulated by Šuvakov & Dmitrašinović).',
    note: 'A retrograde periodic orbit of three equal masses with two bodies on nested loops.',
    duration: 20,
    bodies: [
      { m: 1, x: [0.3361300950, 0, 0], v: [0, 1.5324315370, 0] },
      { m: 1, x: [0.7699893804, 0, 0], v: [0, -0.6287350978, 0] },
      { m: 1, x: [-1.1061194753, 0, 0], v: [0, -0.9036964391, 0] },
    ],
  },
};

/** Equilateral triangle in the plane rotating rigidly about the barycentre. */
export function lagrangeTriangle(m, side, G = 1) {
  const M = m[0] + m[1] + m[2];
  const pts = [0, 1, 2].map(k => {
    const th = Math.PI / 2 + k * 2 * Math.PI / 3;
    return [side / Math.sqrt(3) * Math.cos(th), side / Math.sqrt(3) * Math.sin(th), 0];
  });
  const c = [0, 1, 2].map(i => pts.reduce((s, p, j) => s + m[j] * p[i], 0) / M);
  // Angular velocity of the rigid rotation: Ω² = G M / side³.
  const omega = Math.sqrt(G * M / (side * side * side));
  return m.map((mi, j) => {
    const r = [pts[j][0] - c[0], pts[j][1] - c[1], 0];
    return { m: mi, x: r, v: [-omega * r[1], omega * r[0], 0] };
  });
}

function hierarchicalTriple() {
  const m1 = 1, m2 = 1, m3 = 0.1;
  const aIn = 0.1, aOut = 1, eOut = 0.2;
  const vIn = Math.sqrt((m1 + m2) / aIn);
  const Mtot = m1 + m2 + m3;
  // Outer orbit starts at pericentre; relative velocity from vis-viva.
  const rOut = aOut * (1 - eOut);
  const vOut = Math.sqrt(Mtot * (1 + eOut) / rOut);
  const fIn = m3 / Mtot, fOut = (m1 + m2) / Mtot;
  return [
    { m: m1, x: [-fIn * rOut - aIn / 2, 0, 0], v: [0, -fIn * vOut - vIn / 2, 0] },
    { m: m2, x: [-fIn * rOut + aIn / 2, 0, 0], v: [0, -fIn * vOut + vIn / 2, 0] },
    { m: m3, x: [fOut * rOut, 0, 0], v: [0, fOut * vOut, 0] },
  ];
}
