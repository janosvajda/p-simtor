// Clocks in weak gravitational fields (special + general relativity, post-Newtonian order).
//
// A clock moving at speed v at gravitational potential U = GM/r (positive) ticks, relative to
// Geocentric Coordinate Time t (TCG; Terrestrial Time TT = TCG·(1 − L_G)), at
//     dτ/dt = 1 − (U + v²/2) / c²              (accurate to O(1/c⁴) ≈ 10⁻¹⁹ near Earth)
// A clock at rest on the rotating geoid ticks at dτ/dt = 1 − W₀/c², where W₀ is the geoid
// potential (gravity + centrifugal). Constants: IERS Conventions 2010, IAU 2000, GPS ICD.

import { c } from './constants.mjs';

export const GM_EARTH = 3.986004418e14;            // m³/s² (IERS 2010 / WGS 84)
export const W0 = 62636856.0;                      // m²/s², geoid potential (IERS 2010)
export const L_G = 6.969290134e-10;                // IAU 2000 defining constant, = W₀/c²
export const OMEGA_EARTH = 7.2921151467e-5;        // rad/s (WGS 84)
export const R_EARTH = 6378137.0;                  // m, equatorial radius (WGS 84)
export const GPS_A = 26561.75e3;                   // m, nominal GPS semi-major axis
export const GPS_F0 = 10.23e6;                     // Hz, GPS fundamental frequency

/** Rate of a ground clock on the geoid relative to coordinate time: 1 − W₀/c². */
export const groundRate = () => 1 - W0 / (c * c);

/** Fractional rate of an orbiting clock relative to coordinate time, split into its terms. */
export function orbitRateTerms(r, v2) {
  const grav = GM_EARTH / r / (c * c);             // gravitational slowing (smaller far from Earth)
  const kin = v2 / 2 / (c * c);                    // time dilation from speed
  return { grav, kin, rate: 1 - grav - kin };
}

/**
 * Satellite clock minus ground clock, fractional rate: (W₀ − GM/r − v²/2)/c².
 * Positive: the satellite clock runs fast.
 */
export const fractionalOffset = (r, v2) => (W0 - GM_EARTH / r - v2 / 2) / (c * c);

/** Same for a circular orbit of radius a (v² = GM/a): (W₀ − 3GM/2a)/c². */
export const circularOffset = a => (W0 - 1.5 * GM_EARTH / a) / (c * c);

/** GPS ICD relativistic correction constant F = −2√(GM)/c² (s/√m), IS-GPS-200 §20.3.3.3.3.1. */
export const F_ICD = -2 * Math.sqrt(GM_EARTH) / (c * c);

/** Periodic eccentricity term Δt_r = F e √a sin E (s). */
export const eccentricityTerm = (e, a, E) => F_ICD * e * Math.sqrt(a) * Math.sin(E);

export const MICRO_PER_DAY = 86400e6;              // multiply a fractional rate to get µs/day
