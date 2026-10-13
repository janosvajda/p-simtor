// Thermal (black-body) radiation and the colour of light.
//
// Planck's law, its integrals (Stefan–Boltzmann, photon statistics), Wien's
// displacement law from the exact root, and the perceived colour of a spectrum
// via CIE 1931 colour matching → sRGB. SI units throughout.

import { h, c, kB } from './constants.mjs';
import { brent, gaussLegendre } from './numerics.mjs';

export const c1L = 2 * h * c * c;                 // 2hc²  (W m² sr⁻¹), spectral radiance prefactor
export const c2 = h * c / kB;                     // hc/k  (m K), second radiation constant

/** Spectral radiance per unit wavelength B_λ(λ, T), W sr⁻¹ m⁻³ (λ in m). */
export function planckLambda(lambda, T) {
  return c1L / Math.pow(lambda, 5) / Math.expm1(c2 / (lambda * T));
}

/** Spectral radiance per unit frequency B_ν(ν, T), W sr⁻¹ m⁻² Hz⁻¹. */
export function planckNu(nu, T) {
  return 2 * h * nu * nu * nu / (c * c) / Math.expm1(h * nu / (kB * T));
}

/** Rayleigh–Jeans limit (classical; diverges at short λ: the ultraviolet catastrophe). */
export const rayleighJeans = (lambda, T) => 2 * c * kB * T / Math.pow(lambda, 4);
/** Wien approximation (valid for hc/λkT ≫ 1). */
export const wienApprox = (lambda, T) => c1L / Math.pow(lambda, 5) * Math.exp(-c2 / (lambda * T));

/** x such that x = 5(1 − e^{−x}): the peak of B_λ is at λ_max = hc/(x k T). */
export const WIEN_X = brent(x => x - 5 * (1 - Math.exp(-x)), 1, 10);
/** Wien displacement constant b = hc/(x k), m K. */
export const WIEN_B = c2 / WIEN_X;
/** Peak of B_ν: x = 3(1 − e^{−x}); ν_max = x k T / h. */
export const WIEN_X_NU = brent(x => x - 3 * (1 - Math.exp(-x)), 1, 10);

/** Stefan–Boltzmann constant from first principles: σ = 2π⁵k⁴ / (15 h³ c²). */
export const SIGMA_SB = 2 * Math.pow(Math.PI, 5) * Math.pow(kB, 4) / (15 * h * h * h * c * c);

export const peakWavelength = T => WIEN_B / T;
/** Radiant exitance (power per area leaving a black surface), W m⁻². */
export const exitance = T => SIGMA_SB * T * T * T * T;
/** Photon number density in a cavity, m⁻³: n = 16π ζ(3) (kT/hc)³. */
export const ZETA3 = 1.2020569031595942;
export const photonDensity = T => 16 * Math.PI * ZETA3 * Math.pow(kB * T / (h * c), 3);
/** Mean photon energy in units of kT: π⁴/(30 ζ(3)) ≈ 2.701. */
export const MEAN_PHOTON_ENERGY_KT = Math.pow(Math.PI, 4) / (30 * ZETA3);

/** ∫ B_λ dλ over [λ0, λ1] (m), by composite Gauss–Legendre in ln λ. W sr⁻¹ m⁻². */
export function bandRadiance(T, l0, l1, panels = 24) {
  const a = Math.log(l0), b = Math.log(l1);
  let s = 0;
  for (let i = 0; i < panels; i++) {
    s += gaussLegendre(u => { const l = Math.exp(u); return planckLambda(l, T) * l; }, a + (b - a) * i / panels, a + (b - a) * (i + 1) / panels, 24);
  }
  return s;
}

// ---------------------------------------------------------------- colour

// CIE 1931 2° colour matching functions, multi-lobe Gaussian fit
// (Wyman, Sloan & Shirley, J. Computer Graphics Techniques 2(2), 2013, eqs. 2–4). λ in nm.
const g = (x, mu, s1, s2) => { const t = (x - mu) / (x < mu ? s1 : s2); return Math.exp(-0.5 * t * t); };
export function cieXYZbar(nm) {
  // Beyond 700 nm the CIE 1931 spectral locus is a single point, x = 0.7347, y = 0.2653
  // (x̄/ȳ = 2.769, z̄ = 0 in the tabulated data). The Gaussian fit's tails do not keep that
  // ratio, which turns cold black bodies green; impose it, keeping the fit's ȳ fall-off.
  if (nm > 700) {
    const y = 0.821 * g(nm, 568.8, 46.9, 40.5) + 0.286 * g(nm, 530.9, 16.3, 31.1);
    return [y * (0.734690 / 0.265310), y, 0];
  }
  return [
    1.056 * g(nm, 599.8, 37.9, 31.0) + 0.362 * g(nm, 442.0, 16.0, 26.7) - 0.065 * g(nm, 501.1, 20.4, 26.2),
    0.821 * g(nm, 568.8, 46.9, 40.5) + 0.286 * g(nm, 530.9, 16.3, 31.1),
    1.217 * g(nm, 437.0, 11.8, 36.0) + 0.681 * g(nm, 459.0, 26.0, 13.8),
  ];
}

/** CIE XYZ tristimulus of a spectral power distribution S(λ in nm), 360–830 nm. */
export function spectrumToXYZ(S, step = 1) {
  let X = 0, Y = 0, Z = 0;
  for (let nm = 360; nm <= 830; nm += step) {
    const [x, y, z] = cieXYZbar(nm), s = S(nm);
    X += s * x; Y += s * y; Z += s * z;
  }
  return [X * step, Y * step, Z * step];
}

export function chromaticity([X, Y, Z]) {
  const s = X + Y + Z;
  return [X / s, Y / s];
}

/** XYZ → linear sRGB, D65 white (7-digit matrix derived from the sRGB primaries; IEC 61966-2-1 lists it to 4 digits). */
export function xyzToLinearSRGB([X, Y, Z]) {
  return [
    3.2404542 * X - 1.5371385 * Y - 0.4985314 * Z,
    -0.9692660 * X + 1.8760108 * Y + 0.0415560 * Z,
    0.0556434 * X - 0.2040259 * Y + 1.0572252 * Z,
  ];
}

export const srgbGamma = v => v <= 0.0031308 ? 12.92 * v : 1.055 * Math.pow(v, 1 / 2.4) - 0.055;

/**
 * Displayable colour of a black body at temperature T: chromaticity is exact (within the
 * CMF fit), brightness is normalised so the largest channel is 1 (disclose this: real
 * brightness spans ~16 orders of magnitude over 1000–40000 K). Out-of-gamut colours (below
 * ~1900 K) have negative channels clipped, which keeps their hue. Chromaticity accuracy vs the
 * CIE Planckian locus: ≤ 0.001 above 2500 K, ≤ 0.008 at 1000 K (fit near 700 nm).
 * Returns { rgb: [0–255]×3, xy, inGamut }.
 */
export function blackbodyColor(T) {
  // Chromaticity only depends on the spectrum's shape, so rescale it to avoid underflow for
  // cold bodies (at 2.7 K, B_λ in the visible is ~10⁻²³⁰⁰: exactly 0 in floating point).
  // B_λ ∝ λ⁻⁵ e^{−x}/(1 − e^{−x}) with x = c₂/(λT); multiply by e^{x₀}, x₀ at 830 nm.
  const x0 = c2 / (830e-9 * T);
  const xyz = spectrumToXYZ(nm => {
    const x = c2 / (nm * 1e-9 * T);
    return Math.pow(nm / 830, -5) * Math.exp(x0 - x) / -Math.expm1(-x);
  });
  let lin = xyzToLinearSRGB(xyz);
  const inGamut = Math.min(...lin) >= 0;
  lin = lin.map(v => Math.max(0, v));                // clip out-of-gamut channels: keeps the hue
  const mx = Math.max(...lin);
  const rgb = lin.map(v => Math.round(255 * srgbGamma(v / mx)));
  return { rgb, xy: chromaticity(xyz), inGamut };
}

/** Colour of monochromatic light at λ (nm), for drawing a spectrum strip (normalised, disclosed). */
export function wavelengthColor(nm) {
  const lin = xyzToLinearSRGB(cieXYZbar(nm)).map(v => Math.max(0, v));   // clip: keeps the hue
  const mx = Math.max(...lin) || 1;
  return lin.map(v => Math.round(255 * srgbGamma(v / mx)));
}
