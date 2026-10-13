// Light and matter around a non-rotating black hole (Schwarzschild metric).
// Geometric units G = c = 1; lengths in units of the mass M (so r_s = 2).
//
// Photon orbits obey the exact orbit equation (u = 1/r, φ = orbital angle):
//   (du/dφ)² = 1/b² − u²(1 − 2Mu)      ⇔      d²u/dφ² = −u + 3Mu²
// with b the impact parameter. Everything here is exact GR for test particles;
// the ray tracer in the black-hole page runs the same equation on the GPU.

import { brent, gaussLegendre, ellipticK, dopri5 } from './numerics.mjs';

export const M = 1;
export const R_S = 2 * M;                          // event horizon
export const R_PHOTON = 3 * M;                     // unstable circular photon orbit
export const B_CRIT = 3 * Math.sqrt(3) * M;        // critical impact parameter (shadow edge)
export const R_ISCO = 6 * M;                       // innermost stable circular orbit

/**
 * Trace a photon coming in from infinity with impact parameter b, integrating
 * d²u/dφ² = −u + 3Mu² with fixed-step RK4 in φ (the same scheme as the GPU shader).
 * Returns { captured, deflection, rMin, steps }; deflection = Δφ − π for escaping rays.
 */
export function tracePhoton(b, { dphi = 1e-3, maxPhi = 40 } = {}) {
  let u = 0, w = 1 / b, phi = 0, rMin = Infinity, steps = 0;
  const f = q => -q + 3 * M * q * q;          // d²u/dφ² as a function of u
  while (phi < maxPhi) {
    const k1u = w, k1w = f(u);
    const k2u = w + dphi / 2 * k1w, k2w = f(u + dphi / 2 * k1u);
    const k3u = w + dphi / 2 * k2w, k3w = f(u + dphi / 2 * k2u);
    const k4u = w + dphi * k3w, k4w = f(u + dphi * k3u);
    const un = u + dphi / 6 * (k1u + 2 * k2u + 2 * k3u + k4u);
    const wn = w + dphi / 6 * (k1w + 2 * k2w + 2 * k3w + k4w);
    steps++;
    if (un >= 1 / R_S) return { captured: true, rMin: R_S, steps, phi };
    if (un <= 0 && phi > 0) {
      // linear interpolation to u = 0 for the exit angle
      const t = u / (u - un);
      return { captured: false, deflection: phi + t * dphi - Math.PI, rMin, steps };
    }
    u = un; w = wn; phi += dphi;
    if (u > 0) rMin = Math.min(rMin, 1 / u);
  }
  return { captured: false, deflection: NaN, rMin, steps, trapped: true };
}

/** Closest approach for b > B_CRIT: largest root u₀ < 1/(3M) of 1/b² − u² + 2Mu³ = 0. */
export function turningPoint(b) {
  if (!(b > B_CRIT)) return NaN;                    // captured: no turning point
  // Solve for z = b u, which remains O(1) even in very weak fields. A fixed
  // absolute tolerance on u would otherwise swamp the tiny deflection.
  const z = brent(q => 1 - q * q + 2 * (M / b) * q * q * q, 1, Math.sqrt(3));
  return z / b;
}

/**
 * Exact deflection by quadrature (independent of the ODE tracer):
 * α(b) = 2 ∫₀^{u₀} du / √(1/b² − u² + 2Mu³) − π, with u = u₀(1 − s²) removing the
 * endpoint singularity. Factor the polynomial at its turning point and subtract
 * the flat-space integrand algebraically, avoiding subtraction of two angles ≈ π.
 */
export function deflectionExact(b, n = 64) {
  if (!(b > B_CRIT)) return Infinity;               // the photon is captured
  const u0 = turningPoint(b);
  return 4 * gaussLegendre(s => {
    const x = 1 - s * s;
    const flat = 1 + x, delta = 2 * M * u0 * (1 + x + x * x);
    const a = Math.sqrt(flat), curved = Math.sqrt(flat - delta);
    return delta / (a * curved * (a + curved));
  }, 0, 1, n);
}

/** Weak-field series α = 4M/b + (15π/4)(M/b)² + (128/3)(M/b)³ + … (Keeton & Petters 2005). */
export const deflectionWeak = (b, order = 3) => {
  const x = M / b;
  return 4 * x + (order >= 2 ? 15 * Math.PI / 4 * x * x : 0) + (order >= 3 ? 128 / 3 * x * x * x : 0);
};

/** Angular radius of the shadow seen by a static observer at r_o: sin α = (B_CRIT/r_o)√(1 − 2M/r_o). */
export const shadowAngle = ro => {
  if (!(ro > R_S)) return NaN;                      // no static observers inside the horizon
  const a = Math.asin(Math.min(1, B_CRIT / ro * Math.sqrt(1 - 2 * M / ro)));   // round-off can exceed 1 at r₀ = 3M
  return ro < R_PHOTON ? Math.PI - a : a;           // inside the photon sphere the shadow exceeds 90°
};

/** Impact parameter of a ray leaving a static observer at r_o at angle α from the radial direction. */
export const impactFromAngle = (ro, alpha) => ro * Math.sin(alpha) / Math.sqrt(1 - 2 * M / ro);

// ---------------------------------------------------------------- matter: circular orbits

/** Coordinate angular velocity of a circular geodesic: Ω = √(M/r³) (exact in Schwarzschild). */
export const orbitalOmega = r => Math.sqrt(M / (r * r * r));
/** u^t of a circular orbit: 1/√(1 − 3M/r). */
export const orbitalUt = r => 1 / Math.sqrt(1 - 3 * M / r);
/** Specific energy and angular momentum of circular orbits. */
export const orbitalEnergy = r => (1 - 2 * M / r) / Math.sqrt(1 - 3 * M / r);
export const orbitalL = r => Math.sqrt(M * r) / Math.sqrt(1 - 3 * M / r);

/**
 * Redshift factor g = E_observed / E_emitted for light from a circular-orbit emitter at r,
 * received by a static observer at r_o, for a photon with L_z/E = λ (along the orbit's axis):
 *   g = 1 / ( u^t (1 − Ω λ) √(1 − 2M/r_o) )
 * Combines gravitational redshift, transverse and longitudinal Doppler shifts exactly.
 */
export const diskRedshift = (r, lambda, ro = Infinity) =>
  1 / (orbitalUt(r) * (1 - orbitalOmega(r) * lambda) * (Number.isFinite(ro) ? Math.sqrt(1 - 2 * M / ro) : 1));

/**
 * Thin-disk temperature profile with zero torque at the ISCO (Shakura & Sunyaev 1973,
 * Newtonian form): T(r) = T* (r/r_isco)^{−3/4} (1 − √(r_isco/r))^{1/4}, normalised so its
 * maximum is T_max. A modelling approximation (the full relativistic Page–Thorne profile
 * differs by tens of per cent near the ISCO): disclose it.
 */
export function diskTemperature(r, Tmax) {
  if (r < R_ISCO) return 0;
  const shape = x => Math.pow(x, -0.75) * Math.pow(1 - Math.sqrt(1 / x), 0.25);
  const xPeak = 49 / 36;                           // maximum of the shape function
  return Tmax * shape(r / R_ISCO) / shape(xPeak);
}

/** Effective potential for massive particles, per unit mass: V² = (1 − 2M/r)(1 + L²/r²). */
export const effectivePotential = (r, L) => (1 - 2 * M / r) * (1 + L * L / (r * r));

// ---------------------------------------------------------------- massive bodies: relativistic orbits
//
// Bound orbits of a test mass, described by semi-latus rectum p and eccentricity e (in units of
// M = GM/c²): turning points r = p/(1 ± e). Exact relations from Cutler, Kennefick & Poisson,
// Phys. Rev. D 50, 3816 (1994). Valid for p > 6 + 2e (outside the separatrix).

/** Specific angular momentum² and energy² of the bound orbit (p, e). */
export function orbitConstants(p, e) {
  const d = p - 3 - e * e;
  return { L2: p * p / d, E2: ((p - 2) ** 2 - 4 * e * e) / (p * d) };
}

/**
 * Exact periapsis advance per radial period (radians):
 * Δφ = 4 √(p / (p − 6 + 2e)) K(4e / (p − 6 + 2e)) − 2π.
 */
export function precessionPerOrbit(p, e) {
  // separatrix p = 6 + 2e (where 4e/(p − 6 + 2e) reaches 1 and K diverges): no bound orbit there
  if (!(p > 6 + 2 * e)) return Infinity;
  const q = p - 6 + 2 * e;
  return 4 * Math.sqrt(p / q) * ellipticK(Math.min(4 * e / q, 1)) - 2 * Math.PI;
}

/** Weak-field series 6π/p + (3π/2)(18 + e²)/p² (Damour & Schäfer 1988), p in units of M. */
export const precessionWeak = (p, e) => 6 * Math.PI / p + 1.5 * Math.PI * (18 + e * e) / (p * p);

/**
 * In SI-style inputs: precession per orbit (rad) for semi-major axis a, eccentricity e around a
 * body with gravitational parameter GM. Weak or strong field alike (uses the exact formula).
 */
export function precessionSI({ GM, a, e, c = 299792458 }) {
  const Mgeo = GM / (c * c);
  return precessionPerOrbit(a * (1 - e * e) / Mgeo, e);
}

/**
 * Integrate the orbit equation d²u/dφ² = M/L² − u + 3Mu² (u = 1/r, M = 1) from periapsis,
 * together with coordinate time: dt/dφ = (E/L) r² / (1 − 2M/r) (Newtonian: r²/L).
 * Returns sampled points [φ, r, t], the angles and times of successive periapsis passages, and
 * at(t) giving the position at any time (for animation). precessionPerOrbit gives the exact rate.
 * Units: lengths in M = GM/c², times in GM/c³ (multiply by GM/c³ seconds for SI).
 */
export function relativisticOrbit(p, e, { orbits = 3, samplesPerOrbit = 400, newtonian = false } = {}) {
  if (!(Number.isFinite(e) && e >= 0 && e < 1 && Number.isFinite(p) && p > (newtonian ? 0 : 6 + 2 * e))) throw new Error('relativisticOrbit: need 0 ≤ e < 1 and p outside the bound-orbit separatrix');
  if (!(Number.isSafeInteger(orbits) && orbits > 0 && Number.isSafeInteger(samplesPerOrbit) && samplesPerOrbit >= 8)) throw new Error('relativisticOrbit: positive integer orbits and at least 8 samples per orbit are required');
  const { L2, E2 } = newtonian ? { L2: p, E2: 1 } : orbitConstants(p, e);
  const L = Math.sqrt(L2), EoverL = Math.sqrt(E2) / L;
  const k = newtonian ? 0 : 3;
  const points = [[0, p / (1 + e), 0]], periapsides = [0], periapsisTimes = [0];
  let last = null;
  // Strong-field radial periods can cover many azimuthal revolutions. Use the
  // actual radial period, with room to bracket the last periapsis crossing.
  const radialAngle = newtonian ? 2 * Math.PI : 2 * Math.PI + precessionPerOrbit(p, e);
  const hmax = 2 * Math.PI / samplesPerOrbit;
  const phiEnd = orbits * radialAngle + 2 * hmax;
  const dtdphi = u => (newtonian ? 1 / (L * u * u) : EoverL / (u * u * (1 - 2 * u)));
  if (e === 0) {
    // A circle has no radial periapsis; count complete azimuthal revolutions.
    const dt = dtdphi(1 / p);
    for (let i = 1; i <= orbits * samplesPerOrbit; i++) {
      const phi = 2 * Math.PI * i / samplesPerOrbit;
      points.push([phi, p, phi * dt]);
    }
    for (let i = 1; i <= orbits; i++) { periapsides.push(2 * Math.PI * i); periapsisTimes.push(2 * Math.PI * i * dt); }
  } else {
    dopri5((_phi, [u, w]) => [w, 1 / L2 - u + k * u * u, dtdphi(u)], 0, [(1 + e) / p, 0, 0], phiEnd, {
      rtol: 1e-12, atol: 1e-16, hmax,
      onStep(phi, y) {
        points.push([phi, 1 / y[0], y[2]]);
        if (last && last.w > 0 && y[1] <= 0) {
          // Locate w = 0 within the accepted step using cubic Hermite output.
          // A single Newton correction leaves an O(h³) event-angle error.
          const h = phi - last.phi;
          const interp = (a, da, b, db, s) => (2 * s ** 3 - 3 * s * s + 1) * a
            + (s ** 3 - 2 * s * s + s) * h * da + (-2 * s ** 3 + 3 * s * s) * b + (s ** 3 - s * s) * h * db;
          const f0 = 1 / L2 - last.u + k * last.u * last.u, f1 = 1 / L2 - y[0] + k * y[0] * y[0];
          const s = brent(q => interp(last.w, f0, y[1], f1, q), 0, 1);
          periapsides.push(last.phi + s * h);
          periapsisTimes.push(interp(last.t, dtdphi(last.u), y[2], dtdphi(y[0]), s));
          if (periapsides.length > orbits) return false;
        }
        last = { phi, u: y[0], w: y[1], t: y[2] };
        return true;
      },
    });
  }
  if (periapsides.length !== orbits + 1) throw new Error('relativisticOrbit: failed to resolve every requested periapsis');
  /** Position at coordinate time t (linear interpolation between samples: ~1e−5 relative in r, sub-pixel). */
  const at = t => {
    let lo = 0, hi = points.length - 1;
    if (t <= 0) return { phi: 0, r: points[0][1], x: points[0][1], y: 0 };
    if (t >= points[hi][2]) { const [phi, r] = points[hi]; return { phi, r, x: r * Math.cos(phi), y: r * Math.sin(phi) }; }
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; points[mid][2] <= t ? (lo = mid) : (hi = mid); }
    const [p0, r0, t0] = points[lo], [p1, r1, t1] = points[hi], f = (t - t0) / (t1 - t0);
    const phi = p0 + f * (p1 - p0), r = r0 + f * (r1 - r0);
    return { phi, r, x: r * Math.cos(phi), y: r * Math.sin(phi) };
  };
  return { points, periapsides, periapsisTimes, at, duration: points[points.length - 1][2] };
}
