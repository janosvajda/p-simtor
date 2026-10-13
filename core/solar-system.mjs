// Solar System initial conditions from published Keplerian elements.
//
// Elements: E. M. Standish, "Keplerian Elements for Approximate Positions of the Major
// Planets", JPL Solar System Dynamics, Table 1 (valid 1800–2050 AD), J2000 ecliptic and
// equinox. Planet/Sun mass ratios: IAU 2009 / JPL DE430.
// Units: AU, Julian years, solar masses, G = GM☉·yr²/AU³ = 39.476926 (≈ 4π²; see constants.mjs).

import { G_AU_MSUN_YR } from './constants.mjs';
import { solveKepler } from './kepler.mjs';

const DEG = Math.PI / 180;

// a (AU), e, I (deg), L mean longitude (deg), ϖ longitude of perihelion (deg), Ω node (deg),
// plus rates per Julian century, and the Sun/planet mass ratio.
export const PLANETS = [
  { name: 'Mercury', a: 0.38709927, e: 0.20563593, I: 7.00497902, L: 252.25032350, w: 77.45779628, O: 48.33076593, dL: 149472.67411175, sunOverM: 6023600 },
  { name: 'Venus', a: 0.72333566, e: 0.00677672, I: 3.39467605, L: 181.97909950, w: 131.60246718, O: 76.67984255, dL: 58517.81538729, sunOverM: 408523.719 },
  { name: 'Earth', a: 1.00000261, e: 0.01671123, I: -0.00001531, L: 100.46457166, w: 102.93768193, O: 0.0, dL: 35999.37244981, sunOverM: 328900.56, note: 'Earth–Moon barycentre' },
  { name: 'Mars', a: 1.52371034, e: 0.09339410, I: 1.84969142, L: -4.55343205, w: -23.94362959, O: 49.55953891, dL: 19140.30268499, sunOverM: 3098703.59 },
  { name: 'Jupiter', a: 5.20288700, e: 0.04838624, I: 1.30439695, L: 34.39644051, w: 14.72847983, O: 100.47390909, dL: 3034.74612775, sunOverM: 1047.348644 },
  { name: 'Saturn', a: 9.53667594, e: 0.05386179, I: 2.48599187, L: 49.95424423, w: 92.59887831, O: 113.66242448, dL: 1222.49362201, sunOverM: 3497.9018 },
  { name: 'Uranus', a: 19.18916464, e: 0.04725744, I: 0.77263783, L: 313.23810451, w: 170.95427630, O: 74.01692503, dL: 428.48202785, sunOverM: 22902.98 },
  { name: 'Neptune', a: 30.06992276, e: 0.00859048, I: 1.77004347, L: -55.12002969, w: 44.96476227, O: 131.78422574, dL: 218.45945325, sunOverM: 19412.26 },
];

/**
 * Known natural satellites per planet, as of February 2024 (NASA Solar System Exploration /
 * IAU Minor Planet Center). Counts of small outer-planet moons keep growing: quote the date.
 */
export const MOONS = {
  asOf: '2024-02', source: 'NASA Solar System Exploration',
  count: { Mercury: 0, Venus: 0, Earth: 1, Mars: 2, Jupiter: 95, Saturn: 146, Uranus: 28, Neptune: 16 },
  largest: { Earth: 'Moon', Mars: 'Phobos', Jupiter: 'Ganymede', Saturn: 'Titan', Uranus: 'Titania', Neptune: 'Triton' },
};

/**
 * Heliocentric ecliptic J2000 position (AU) and velocity (AU/yr) of a planet at T Julian
 * centuries from J2000 (only the mean longitude is advanced; the other elements are J2000).
 */
export function planetState(p, T = 0) {
  const m = 1 / p.sunOverM;
  const mu = G_AU_MSUN_YR * (1 + m);
  const L = (p.L + p.dL * T) * DEG, w = p.w * DEG, O = p.O * DEG, I = p.I * DEG;
  const omega = w - O;                                  // argument of perihelion
  const M = L - w;                                      // mean anomaly
  const E = solveKepler(((M % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI), p.e);
  const n = Math.sqrt(mu / (p.a ** 3));
  const Edot = n / (1 - p.e * Math.cos(E));
  // orbital plane (x toward perihelion)
  const xp = p.a * (Math.cos(E) - p.e), yp = p.a * Math.sqrt(1 - p.e * p.e) * Math.sin(E);
  const vxp = -p.a * Math.sin(E) * Edot, vyp = p.a * Math.sqrt(1 - p.e * p.e) * Math.cos(E) * Edot;
  // rotate by ω, I, Ω into the ecliptic frame
  const co = Math.cos(omega), so = Math.sin(omega), cO = Math.cos(O), sO = Math.sin(O), cI = Math.cos(I), sI = Math.sin(I);
  const rot = (x, y) => [
    (co * cO - so * sO * cI) * x + (-so * cO - co * sO * cI) * y,
    (co * sO + so * cO * cI) * x + (-so * sO + co * cO * cI) * y,
    (so * sI) * x + (co * sI) * y,
  ];
  return { x: rot(xp, yp), v: rot(vxp, vyp), m, mu };
}

/** Sidereal period (yr) from Kepler's third law with the planet's mass included. */
export const keplerPeriod = p => 2 * Math.PI * Math.sqrt(p.a ** 3 / (G_AU_MSUN_YR * (1 + 1 / p.sunOverM)));

/** A gravity-template scenario (psim gravity --scenario) for the Sun and chosen planets at J2000. */
export function solarSystemScenario({ planets = PLANETS.map(p => p.name), years = 30, stepsPerCrossing = 400 } = {}) {
  const chosen = PLANETS.filter(p => planets.includes(p.name));
  return {
    key: 'solar-system',
    title: 'Solar System at J2000',
    pageTitle: 'The Solar System',
    subtitle: `Sun and ${chosen.length} planets from JPL elements (J2000) · full N-body · AU, years`,
    note: 'Every planet pulls on every other: watch Jupiter and Saturn perturb each other. Positions start at 2000-01-01 12:00 TT.',
    source: 'Standish, JPL Keplerian elements Table 1 (1800–2050); masses IAU 2009 / DE430',
    units: { G: G_AU_MSUN_YR, system: 'AU, solar masses, Julian years (G = GM☉ yr²/AU³ ≈ 4π²)', length: 'AU', time: 'yr', mass: 'M☉' },
    integrator: 'logh4',
    stepsPerCrossing,
    duration: years,
    speed: 1,
    trailTime: Math.max(2, years / 2),
    bodies: [
      { name: 'Sun', m: 1, x: [0, 0, 0], v: [0, 0, 0] },
      ...chosen.map(p => { const s = planetState(p); return { name: p.name, m: s.m, x: s.x, v: s.v }; }),
    ],
    notes: [
      'Newtonian point masses; general relativity adds 43″/century to Mercury’s perihelion (not modelled here; see tests/real-world.test.mjs).',
      'Earth is the Earth–Moon barycentre. Moons, asteroids and Pluto are omitted.',
      'Initial conditions use JPL fitted approximate Keplerian elements, converted to two-body state vectors. These are not osculating ephemeris elements; numerical convergence does not measure their observational accuracy.',
    ],
  };
}

/** Three views of the same J2000 initial conditions, for one page with a scenario picker. */
export function solarSystemScenarios() {
  const inner = ['Mercury', 'Venus', 'Earth', 'Mars'], outer = ['Jupiter', 'Saturn', 'Uranus', 'Neptune'];
  return [
    { ...solarSystemScenario({ planets: inner, years: 4, stepsPerCrossing: 400 }), key: 'inner', title: 'Inner planets (Mercury to Mars)', speed: 0.25, trailTime: 2,
      note: 'Four rocky planets over four years. Mercury laps the Sun every 88 days; Mars takes 1.88 years.' },
    { ...solarSystemScenario({ planets: outer, years: 170, stepsPerCrossing: 400 }), key: 'outer', title: 'Outer planets (Jupiter to Neptune)', speed: 6, trailTime: 170,
      note: 'The giants over one Neptune year (165 yr). Jupiter and Saturn are near their 5:2 resonance: their mutual tug is visible in the energy exchange.' },
    { ...solarSystemScenario({ years: 30 }), key: 'all', title: 'Whole Solar System', speed: 1, trailTime: 15 },
  ];
}
