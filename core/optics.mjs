// Waves and optics. SI units throughout (lengths in metres, angles in radians) unless stated.
//   fresnelCS                 Fresnel integrals C(x), S(x) (series + continued fraction, ~1e-15)
//   singleSlit, doubleSlit, grating, circularAperture, fresnelEdge   closed-form diffraction patterns
//   fraunhofer                far field of a sampled 1D/2D aperture by FFT (exact DFT of the samples)
//   angularSpectrum           exact scalar propagation of a sampled complex field (periodic window)
//   gaussianBeam, beamFromQ, propagateQ, abcd, abcdSystem, cavityMode   Gaussian beams and resonators
//   reflect, refract, fresnelCoefficients, brewsterAngle, criticalAngle  rays at interfaces
//   traceRay, paraxialSystem  sequential tracing through spherical surfaces; paraxial ABCD of the same
//   sellmeier, GLASS, abbeNumber, rainbow    dispersion and the Descartes rainbow
//   thinFilm                  transfer-matrix method for multilayers (s/p, oblique, absorbing layers)
//   jones                     Jones calculus and Stokes parameters
//   fdtd2D, yeeDispersion     2D Yee FDTD (TM or TE) with PEC / periodic / Mur walls
//
// Complex numbers are plain { re, im } objects. Time convention e^{i(kz − ωt)}: an absorbing
// medium has N = n + iκ with κ ≥ 0, and evanescent fields decay as e^{−κz}.

import { fft, besselJ } from './numerics.mjs';
import { fftND } from './fields.mjs';
import { c as C0, eps0, mu0 } from './constants.mjs';

// ---------------------------------------------------------------- small helpers

const cx = (re, im = 0) => ({ re, im });
const cadd = (a, b) => cx(a.re + b.re, a.im + b.im);
const csub = (a, b) => cx(a.re - b.re, a.im - b.im);
const cmul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
const cdiv = (a, b) => {
  const scale = Math.max(Math.abs(b.re), Math.abs(b.im));
  if (!(scale > 0 && Number.isFinite(scale))) throw new Error('optics: division by a zero or non-finite complex number');
  const br = b.re / scale, bi = b.im / scale, ar = a.re / scale, ai = a.im / scale, d = br * br + bi * bi;
  return cx((ar * br + ai * bi) / d, (ai * br - ar * bi) / d);
};
const cscale = (a, s) => cx(a.re * s, a.im * s);
const cconj = a => cx(a.re, -a.im);
const cabs2 = a => a.re * a.re + a.im * a.im;
const csqrt = a => {
  const r = Math.hypot(a.re, a.im);
  if (r === 0) return cx(0, 0);
  // Recover the smaller component from 2 Re(sqrt z) Im(sqrt z) = Im z.
  if (a.re >= 0) { const re = Math.sqrt(r / 2 + a.re / 2); return cx(re, a.im / (2 * re)); }
  const im = (a.im < 0 ? -1 : 1) * Math.sqrt(r / 2 - a.re / 2);
  return cx(Math.abs(a.im / (2 * im)), im);
};
const cexp = a => { const m = Math.exp(a.re); return cx(m * Math.cos(a.im), m * Math.sin(a.im)); };
const toC = (v, who) => {
  if (typeof v === 'number') { if (!Number.isFinite(v)) throw new Error(`${who}: index must be finite`); return cx(v, 0); }
  if (v && Number.isFinite(v.re) && Number.isFinite(v.im)) return cx(v.re, v.im);
  throw new Error(`${who}: expected a number or a complex { re, im }`);
};
const pos = (v, name, who) => { if (!(Number.isFinite(v) && v > 0)) throw new Error(`${who}: ${name} must be positive and finite`); };
const fin = (v, name, who) => { if (!Number.isFinite(v)) throw new Error(`${who}: ${name} must be finite`); };

// Exact two-product decomposition x² = high + low, then reduction modulo four before multiplying
// by π/2. This retains the Fresnel phase even when x² has lost its fractional bits in a double.
function fresnelPhase(x) {
  const split = 134217729 * x, hi = split - (split - x), lo = x - hi, square = x * x;
  const error = ((hi * hi - square) + 2 * hi * lo) + lo * lo;
  return Math.PI / 2 * ((square % 4 + error % 4) % 4);
}

// ---------------------------------------------------------------- Fresnel integrals

/**
 * Fresnel integrals C(x) = ∫₀ˣ cos(πt²/2) dt, S(x) = ∫₀ˣ sin(πt²/2) dt (Abramowitz & Stegun 7.3.1–2
 * normalisation), any finite real x; odd in x, C, S → ½ as x → ∞. Power series for |x| ≤ 1.5, complex
 * continued fraction for erfc beyond (Numerical Recipes §6.8); absolute accuracy ~1e-15.
 * ±Infinity gives ±½; NaN throws.
 */
export function fresnelCS(x) {
  if (Number.isNaN(x)) throw new Error('fresnelCS: x must be a number');
  if (!Number.isFinite(x)) return { C: Math.sign(x) / 2, S: Math.sign(x) / 2 };
  const ax = Math.abs(x);
  let Cv, Sv;
  if (ax >= 2 ** 54) { Cv = 0.5; Sv = 0.5; }  // tail < half an ulp of 1/2; avoid overflowing x²
  else if (ax < 1e-150) { Cv = ax; Sv = 0; }
  else if (ax <= 1.5) {
    const fact = Math.PI / 2 * ax * ax;
    let sum = 0, sums = 0, sumc = ax, sign = 1, odd = true, term = ax, n = 3;
    for (let k = 1; k < 200; k++) {
      term *= fact / k;
      sum += sign * term / n;
      const test = Math.abs(sum) * 1e-17;
      if (odd) { sign = -sign; sums = sum; sum = sumc; } else { sumc = sum; sum = sums; }
      if (term < test) break;
      odd = !odd; n += 2;
    }
    Cv = sumc; Sv = sums;
  } else {
    const pix2 = Math.PI * ax * ax;
    let b = cx(1, -pix2), cc = cx(1e300, 0), d = cdiv(cx(1), b), h = d, n = -1;
    for (let k = 2; k < 1000; k++) {
      n += 2;
      const a = -n * (n + 1);
      b = cx(b.re + 4, b.im);
      d = cdiv(cx(1), cadd(cscale(d, a), b));
      cc = cadd(b, cdiv(cx(a), cc));
      const del = cmul(cc, d);
      h = cmul(h, del);
      if (Math.abs(del.re - 1) + Math.abs(del.im) < 1e-16) break;
    }
    h = cmul(h, cx(ax, -ax));
    const phase = fresnelPhase(ax);
    const cs = cmul(cx(0.5, 0.5), csub(cx(1), cmul(cx(Math.cos(phase), Math.sin(phase)), h)));
    Cv = cs.re; Sv = cs.im;
  }
  return x < 0 ? { C: -Cv, S: -Sv } : { C: Cv, S: Sv };
}

// ---------------------------------------------------------------- closed-form diffraction

const sinc = b => (Math.abs(b) < 1e-8 ? 1 - b * b / 6 : Math.sin(b) / b);

/**
 * Fraunhofer single slit of width a (m) at wavelength λ (m): I/I₀ = sinc²β, β = πa sinθ/λ.
 * theta: observation angle (rad) from the normal, normal incidence. Zeros at a sinθ = mλ, m ≠ 0.
 */
export function singleSlit(theta, { a, lambda }) {
  fin(theta, 'theta', 'singleSlit'); pos(a, 'a', 'singleSlit'); pos(lambda, 'lambda', 'singleSlit');
  return sinc(Math.PI * a * Math.sin(theta) / lambda) ** 2;
}

/**
 * Fraunhofer double slit (width a, centre spacing d > a, metres): I/I₀ = sinc²β · cos²γ,
 * γ = πd sinθ/λ. Normalised to 1 at θ = 0.
 */
export function doubleSlit(theta, { a, d, lambda }) {
  return grating(theta, { a, d, N: 2, lambda });
}

/**
 * Fraunhofer grating of N slits of width a, period d (metres): I/I₀ = sinc²β · [sin(Nγ)/(N sinγ)]²,
 * γ = πd sinθ/λ; principal maxima at d sinθ = mλ (limit taken exactly), N − 1 zeros between them.
 * Normalised to 1 at θ = 0. Requires integer N ≥ 1 and d ≥ a.
 */
export function grating(theta, { a, d, N, lambda }) {
  const who = 'grating';
  fin(theta, 'theta', who); pos(a, 'a', who); pos(d, 'd', who); pos(lambda, 'lambda', who);
  if (!(Number.isInteger(N) && N >= 1)) throw new Error('grating: N must be an integer ≥ 1');
  if (N > 1 && d < a) throw new Error('grating: slits overlap (d < a)');
  const s = Math.sin(theta), g = Math.PI * d * s / lambda;
  const m = Math.round(g / Math.PI), e = g - m * Math.PI;      // distance from the nearest principal maximum
  // Reduce γ modulo π, then use sinc. A Taylor expansion based on e alone fails when N e is large.
  const array = (sinc(N * e) / sinc(e)) ** 2;
  return sinc(Math.PI * a * s / lambda) ** 2 * array;
}

/**
 * Fraunhofer circular aperture of diameter D (Airy pattern): I/I₀ = [2J₁(v)/v]², v = πD sinθ/λ.
 * First dark ring at sinθ = (j₁,₁/π) λ/D = 1.219670 λ/D; 83.8 % of the power lies inside it.
 */
export function circularAperture(theta, { D, lambda }) {
  fin(theta, 'theta', 'circularAperture'); pos(D, 'D', 'circularAperture'); pos(lambda, 'lambda', 'circularAperture');
  const v = Math.PI * D * Math.sin(theta) / lambda;
  if (Math.abs(v) < 1e-6) return 1 - v * v / 4;
  return (2 * besselJ(1, v) / v) ** 2;
}

/**
 * Fresnel diffraction at a straight edge (half-plane x < 0 opaque), plane wave λ, screen at z:
 * I/I₀ = ½[(½ + C(v))² + (½ + S(v))²], v = x √(2/(λz)), x (m) measured from the geometric shadow
 * edge (x > 0 illuminated). I/I₀ = ¼ at the edge, peak 1.3704 at v = 1.2172, → 1 far out.
 */
export function fresnelEdge(x, { z, lambda }) {
  fin(x, 'x', 'fresnelEdge'); pos(z, 'z', 'fresnelEdge'); pos(lambda, 'lambda', 'fresnelEdge');
  const { C, S } = fresnelCS(x * Math.sqrt(2 / (lambda * z)));
  return 0.5 * ((0.5 + C) ** 2 + (0.5 + S) ** 2);
}

// ---------------------------------------------------------------- sampled fields: FFT optics

function sampledField({ re, im, shape, dx }, who) {
  if (!Array.isArray(shape) || !(shape.length === 1 || shape.length === 2)) throw new Error(`${who}: shape must be [n] or [nx, ny]`);
  const size = shape.reduce((p, n) => {
    if (!(Number.isInteger(n) && n >= 2 && (n & (n - 1)) === 0)) throw new Error(`${who}: every extent must be a power of two ≥ 2`);
    return p * n;
  }, 1);
  if (!re || re.length !== size) throw new Error(`${who}: re must have product(shape) entries`);
  if (im && im.length !== size) throw new Error(`${who}: im must match re`);
  pos(dx, 'dx', who);
  const R = Float64Array.from(re), I = im ? Float64Array.from(im) : new Float64Array(size);
  for (let k = 0; k < size; k++) if (!Number.isFinite(R[k]) || !Number.isFinite(I[k])) throw new Error(`${who}: field values must be finite`);
  return { R, I, size };
}
function transform(R, I, shape, inverse) {
  if (shape.length === 1) fft(R, I, inverse); else fftND(R, I, shape, inverse);
}

/**
 * Fraunhofer (far-field) pattern of a sampled aperture t(x[,y]) (complex: re, optional im), grid
 * shape [n] or [nx, ny] (axis 0 fastest, powers of two), spacing dx (m), wavelength λ (m).
 * Returns { re, im, intensity, normalized, fx, sinTheta } in centred order (zero frequency at index
 * n/2): U(f) = Σ t e^{−2πi f·x} dx^d (the exact DFT of the samples, aperture centred on index n/2),
 * f = sinθ/λ on the FFT grid f_k = (k − n/2)/(n dx). The samples' DFT approximates the continuous
 * transform only for |f| ≪ 1/(2dx); sinθ beyond ±1 is evanescent (reported, not physical).
 */
export function fraunhofer({ re, im, shape, dx, lambda }) {
  const who = 'fraunhofer';
  pos(lambda, 'lambda', who);
  const { R, I, size } = sampledField({ re, im, shape, dx }, who);
  // Centre the aperture: multiply by (−1)^{index sum} so that output index n/2 is zero frequency,
  // and phase back so that the aperture origin is at index n/2.
  const nx = shape[0], ny = shape[1] ?? 1;
  for (let j = 0; j < ny; j++) { for (let i = 0; i < nx; i++) {
    const k = i + nx * j, s = (i + j) % 2 ? -1 : 1;
    R[k] *= s; I[k] *= s;
  } }
  transform(R, I, shape, false);
  const area = dx ** shape.length, intensity = new Float64Array(size);
  let max = 0;
  for (let j = 0; j < ny; j++) { for (let i = 0; i < nx; i++) {
    const k = i + nx * j;
    // origin shift by n/2 samples in each axis: factor e^{+iπ(k − n/2)} = (−1)^{k − n/2}
    const s = ((i - nx / 2) + (shape.length === 2 ? j - ny / 2 : 0)) % 2 ? -1 : 1;
    R[k] *= s * area; I[k] *= s * area;
    intensity[k] = R[k] * R[k] + I[k] * I[k];
    if (intensity[k] > max) max = intensity[k];
  } }
  const normalized = intensity.map(v => (max > 0 ? v / max : 0));
  const fx = shape.map(n => Float64Array.from({ length: n }, (_, k) => (k - n / 2) / (n * dx)));
  const sinTheta = fx.map(f => f.map(v => v * lambda));
  return { re: R, im: I, intensity, normalized, fx: shape.length === 1 ? fx[0] : fx, sinTheta: shape.length === 1 ? sinTheta[0] : sinTheta };
}

/**
 * Exact scalar (Rayleigh–Sommerfeld, angular-spectrum) propagation of a monochromatic field over a
 * distance z ≥ 0 (m) in a medium of index n: U(k⊥, z) = U(k⊥, 0) e^{i z √((nk₀)² − k⊥²)}, evanescent
 * components decay as e^{−z √(k⊥² − (nk₀)²)}. Grid shape [n] (2D x–z problems) or [nx, ny], powers
 * of two, spacing dx (m). The window is periodic: light leaving one side re-enters the other, so keep
 * the beam well inside (pad) and compare against exact results. Returns { re, im } (new arrays);
 * power Σ|U|² is conserved exactly when no evanescent content is present.
 */
export function angularSpectrum({ re, im, shape, dx, lambda, z, n = 1 }) {
  const who = 'angularSpectrum';
  pos(lambda, 'lambda', who); pos(n, 'n', who);
  if (!(Number.isFinite(z) && z >= 0)) throw new Error('angularSpectrum: z must be finite and ≥ 0');
  const { R, I } = sampledField({ re, im, shape, dx }, who);
  transform(R, I, shape, false);
  const k = 2 * Math.PI * n / lambda, nx = shape[0], ny = shape[1] ?? 1;
  const freq = (i, m) => 2 * Math.PI * (i < m / 2 ? i : i - m) / (m * dx);
  for (let j = 0; j < ny; j++) {
    const ky = shape.length === 2 ? freq(j, ny) : 0;
    for (let i = 0; i < nx; i++) {
      const kx = freq(i, nx), q = k * k - kx * kx - ky * ky, idx = i + nx * j;
      let hr, hi;
      if (q >= 0) { const ph = z * Math.sqrt(q); hr = Math.cos(ph); hi = Math.sin(ph); }
      else { hr = Math.exp(-z * Math.sqrt(-q)); hi = 0; }
      const a = R[idx], b = I[idx];
      R[idx] = a * hr - b * hi; I[idx] = a * hi + b * hr;
    }
  }
  transform(R, I, shape, true);
  return { re: R, im: I };
}

// ---------------------------------------------------------------- Gaussian beams and ABCD

/**
 * Fundamental Gaussian beam with waist w0 (1/e² intensity radius, m) at z = 0, vacuum wavelength λ,
 * medium index n. Returns { w0, zR = πnw0²/λ, divergence = λ/(πnw0) (half-angle, paraxial),
 * w(z), R(z) (Infinity at the waist), gouy(z) = atan(z/zR), q(z) = z + i zR }.
 * Paraxial: valid for w0 ≫ λ/n (w0 ≳ 2λ keeps errors below ~1 %).
 */
export function gaussianBeam({ w0, lambda, n = 1 }) {
  pos(w0, 'w0', 'gaussianBeam'); pos(lambda, 'lambda', 'gaussianBeam'); pos(n, 'n', 'gaussianBeam');
  const zR = Math.PI * n * w0 * w0 / lambda;
  return {
    w0, zR, divergence: lambda / (Math.PI * n * w0),
    w: z => w0 * Math.sqrt(1 + (z / zR) ** 2),
    R: z => (z === 0 ? Infinity : z * (1 + (zR / z) ** 2)),
    gouy: z => Math.atan(z / zR),
    q: z => cx(z, zR),
  };
}

/**
 * Beam parameters from the complex q = z + i zR (m) in a medium of index n at vacuum wavelength λ:
 * 1/q = 1/R − iλ/(πnw²). Returns { w, R, z (distance past the waist), zR, w0 }. Throws if Im q ≤ 0
 * (not a physical, confined beam).
 */
export function beamFromQ(q, lambda, n = 1) {
  pos(lambda, 'lambda', 'beamFromQ'); pos(n, 'n', 'beamFromQ');
  if (!(Number.isFinite(q?.re) && Number.isFinite(q?.im) && q.im > 0)) throw new Error('beamFromQ: need q = { re, im } with Im q > 0');
  const inv = cdiv(cx(1), q);
  const w = Math.sqrt(-lambda / (Math.PI * n * inv.im));
  return { w, R: inv.re === 0 ? Infinity : 1 / inv.re, z: q.re, zR: q.im, w0: Math.sqrt(lambda * q.im / (Math.PI * n)) };
}

/** q' = (Aq + B)/(Cq + D) for a ray-transfer matrix M = [[A, B], [C, D]] (real-angle convention). */
export function propagateQ(q, M) {
  const [[A, B], [Cm, D]] = M;
  return cdiv(cadd(cscale(q, A), cx(B)), cadd(cscale(q, Cm), cx(D)));
}

/**
 * Ray-transfer (ABCD) matrices acting on (y, θ) with θ the real ray angle (not n·θ). R > 0 means the
 * centre of curvature lies downstream (convex towards incoming light for an interface; a concave
 * mirror is R > 0 and focuses with f = R/2). Infinity is a flat surface.
 */
export const abcd = Object.freeze({
  freeSpace: d => { fin(d, 'd', 'abcd.freeSpace'); return [[1, d], [0, 1]]; },
  thinLens: f => { if (f === 0 || Number.isNaN(f)) throw new Error('abcd.thinLens: f must be non-zero'); return [[1, 0], [-1 / f, 1]]; },
  curvedMirror: R => { if (R === 0 || Number.isNaN(R)) throw new Error('abcd.curvedMirror: R must be non-zero'); return [[1, 0], [-2 / R, 1]]; },
  flatInterface: (n1, n2) => { pos(n1, 'n1', 'abcd'); pos(n2, 'n2', 'abcd'); return [[1, 0], [0, n1 / n2]]; },
  curvedInterface: (n1, n2, R) => {
    pos(n1, 'n1', 'abcd'); pos(n2, 'n2', 'abcd');
    if (R === 0 || Number.isNaN(R)) throw new Error('abcd.curvedInterface: R must be non-zero');
    return [[1, 0], [(n1 - n2) / (n2 * R), n1 / n2]];
  },
});

/** Product of ABCD matrices listed in the order light meets them (returns M_last · … · M_first). */
export function abcdSystem(...elements) {
  if (!elements.length) throw new Error('abcdSystem: need at least one matrix');
  let M = [[1, 0], [0, 1]];
  for (const E of elements) {
    M = [[E[0][0] * M[0][0] + E[0][1] * M[1][0], E[0][0] * M[0][1] + E[0][1] * M[1][1]],
      [E[1][0] * M[0][0] + E[1][1] * M[1][0], E[1][0] * M[0][1] + E[1][1] * M[1][1]]];
  }
  return M;
}

/**
 * Two-mirror resonator of length L (m), mirror radii R1, R2 (concave > 0, Infinity flat), vacuum
 * wavelength λ, index n. g_i = 1 − L/R_i. Stable when 0 < g1g2 < 1, or g1 = g2 = 0 (symmetric
 * confocal). Returns { g1, g2, zR, w0, waistFromMirror1, w1, w2, gouy (one-way Gouy phase),
 * transverseModeSpacing (Hz: c·gouy/(2πnL)) }. Throws for unstable or marginal (g1g2 = 0 or 1
 * otherwise) cavities, where the mode size is zero or infinite.
 */
export function cavityMode({ L, R1, R2, lambda, n = 1 }) {
  const who = 'cavityMode';
  pos(L, 'L', who); pos(lambda, 'lambda', who); pos(n, 'n', who);
  for (const R of [R1, R2]) if (R === 0 || Number.isNaN(R)) throw new Error('cavityMode: mirror radii must be non-zero (Infinity = flat)');
  const g1 = 1 - L / R1, g2 = 1 - L / R2, g = g1 * g2;
  let zR, z1;
  if (g1 === g2 && g1 > -1 && g1 < 1) { zR = L / 2 * Math.sqrt((1 + g1) / (1 - g1)); z1 = L / 2; }
  else {
    if (!(g > 0 && g < 1)) throw new Error(`cavityMode: resonator not stable (g1·g2 = ${g}; need 0 < g1g2 < 1)`);
    const den = g1 + g2 - 2 * g;
    zR = L * Math.sqrt(g * (1 - g)) / Math.abs(den);
    z1 = L * g2 * (1 - g1) / den;
  }
  const lam = lambda / n, w0 = Math.sqrt(lam * zR / Math.PI), w = z => w0 * Math.sqrt(1 + (z / zR) ** 2);
  const gouy = Math.atan((L - z1) / zR) + Math.atan(z1 / zR);
  return { g1, g2, zR, w0, waistFromMirror1: z1, w1: w(z1), w2: w(L - z1), gouy, transverseModeSpacing: C0 * gouy / (2 * Math.PI * n * L) };
}

// ---------------------------------------------------------------- rays at interfaces

const dot = (a, b) => a.reduce((s, v, i) => s + v * b[i], 0);
const unit = (v, who) => {
  const m = Math.sqrt(dot(v, v));
  if (!(m > 0 && Number.isFinite(m))) throw new Error(`${who}: vectors must be finite and non-zero`);
  return v.map(x => x / m);
};
function vecs(d, nrm, who) {
  if (!Array.isArray(d) || !Array.isArray(nrm) || d.length !== nrm.length || !(d.length === 2 || d.length === 3)) throw new Error(`${who}: d and normal must be 2D or 3D arrays of equal length`);
  return [unit(d, who), unit(nrm, who)];
}

/** Specular reflection of direction d (2D or 3D array) at a surface with normal nrm (either side). Unit result. */
export function reflect(d, nrm) {
  const [u, m] = vecs(d, nrm, 'reflect'), k = 2 * dot(u, m);
  return u.map((v, i) => v - k * m[i]);
}

/**
 * Snell refraction of direction d from index n1 into n2 at a surface with normal nrm (either
 * orientation; 2D or 3D). Returns the unit transmitted direction, or null on total internal
 * reflection (then use reflect).
 */
export function refract(d, nrm, n1, n2) {
  pos(n1, 'n1', 'refract'); pos(n2, 'n2', 'refract');
  const [u, m0] = vecs(d, nrm, 'refract');
  let m = m0, ci = -dot(u, m);
  if (ci < 0) { m = m.map(v => -v); ci = -ci; }
  const eta = n1 / n2, k = 1 - eta * eta * (1 - ci * ci);
  if (k < 0) return null;
  const f = eta * ci - Math.sqrt(k);
  return unit(u.map((v, i) => eta * v + f * m[i]), 'refract');
}

/**
 * Fresnel amplitude and power coefficients at a planar interface between real indices n1 → n2,
 * incidence angle θ1 ∈ [0, π/2] (rad). Hecht's sign convention (r_p = (n2 cosθ1 − n1 cosθ2)/(n2 cosθ1 + n1 cosθ2),
 * so r_p = −r_s at normal incidence). Amplitudes are complex { re, im }: beyond the critical angle
 * cosθ2 = +i√(…) (time factor e^{−iωt}), |r| = 1 and the phases are the Goos–Hänchen-free TIR phases.
 * Returns { rs, rp, ts, tp, Rs, Rp, Ts, Tp, theta2 (null under TIR), tir }; R + T = 1 exactly.
 */
export function fresnelCoefficients(theta1, n1, n2) {
  pos(n1, 'n1', 'fresnelCoefficients'); pos(n2, 'n2', 'fresnelCoefficients');
  if (!(theta1 >= 0 && theta1 <= Math.PI / 2)) throw new Error('fresnelCoefficients: theta1 must lie in [0, π/2]');
  if (n1 === n2) return { rs: cx(0), rp: cx(0), ts: cx(1), tp: cx(1), Rs: 0, Rp: 0, Ts: 1, Tp: 1, theta2: theta1, tir: false };
  const s1 = Math.sin(theta1), c1 = Math.cos(theta1), s2 = n1 * s1 / n2, tir = s2 > 1;
  const c2 = tir ? cx(0, Math.sqrt(s2 * s2 - 1)) : cx(Math.sqrt(1 - s2 * s2), 0);
  const n1c1 = cx(n1 * c1), n2c2 = cscale(c2, n2), n2c1 = cx(n2 * c1), n1c2 = cscale(c2, n1);
  const rs = cdiv(csub(n1c1, n2c2), cadd(n1c1, n2c2)), ts = cdiv(cx(2 * n1 * c1), cadd(n1c1, n2c2));
  const rp = cdiv(csub(n2c1, n1c2), cadd(n2c1, n1c2)), tp = cdiv(cx(2 * n1 * c1), cadd(n2c1, n1c2));
  const f = tir ? 0 : (n2 * c2.re) / (n1 * c1 || 1);
  const Rs = cabs2(rs), Rp = cabs2(rp);
  return { rs, rp, ts, tp, Rs, Rp, Ts: c1 === 0 ? 0 : f * cabs2(ts), Tp: c1 === 0 ? 0 : f * cabs2(tp), theta2: tir ? null : Math.asin(s2), tir };
}

/** Brewster angle atan(n2/n1) (rad): R_p = 0 there. */
export function brewsterAngle(n1, n2) { pos(n1, 'n1', 'brewsterAngle'); pos(n2, 'n2', 'brewsterAngle'); return Math.atan2(n2, n1); }

/** Critical angle asin(n2/n1) (rad) for n1 > n2; throws when n1 ≤ n2 (no total internal reflection). */
export function criticalAngle(n1, n2) {
  pos(n1, 'n1', 'criticalAngle'); pos(n2, 'n2', 'criticalAngle');
  if (!(n1 > n2)) throw new Error('criticalAngle: total internal reflection needs n1 > n2');
  return Math.asin(n2 / n1);
}

// ---------------------------------------------------------------- sequential ray tracing

function checkSurfaces(surfaces, n0, who) {
  pos(n0, 'n0', who);
  if (!Array.isArray(surfaces) || !surfaces.length) throw new Error(`${who}: need a non-empty surface list`);
  let last = -Infinity;
  for (const s of surfaces) {
    fin(s.z, 'surface z', who); pos(s.n, 'surface n', who);
    if (s.R === 0 || Number.isNaN(s.R) || s.R === undefined) throw new Error(`${who}: surface R must be non-zero (Infinity = flat)`);
    if (!(s.z >= last)) throw new Error(`${who}: surfaces must be ordered by vertex z`);
    if (s.aperture !== undefined) pos(s.aperture, 'aperture', who);
    last = s.z;
  }
}

/**
 * Exact (real) ray trace through coaxial spherical surfaces along +z. ray = { p: [x, y, z], d: [dx, dy, dz] }
 * (dz > 0; 3D, use x = dx = 0 for meridional 2D rays); surfaces = [{ z (vertex, m), R (centre at z + R,
 * Infinity flat), n (index after the surface), aperture? (clear semi-diameter) }] ordered by z; n0 =
 * index before the first surface. Returns { points (start and each hit), p, d (unit), n, status }:
 * status 'ok', 'missed' (no intersection with the sphere), 'clipped' (outside an aperture) or 'tir';
 * on failure p, d are the last valid ray. Non-'ok' rays must not be drawn as transmitted.
 */
export function traceRay(ray, surfaces, n0 = 1) {
  checkSurfaces(surfaces, n0, 'traceRay');
  let p = Array.from(ray?.p ?? []), d = Array.from(ray?.d ?? []);
  if (p.length !== 3 || d.length !== 3 || !p.every(Number.isFinite)) throw new Error('traceRay: ray needs 3D p and d');
  d = unit(d, 'traceRay');
  if (!(d[2] > 0)) throw new Error('traceRay: ray must travel towards +z');
  const points = [p.slice()];
  let n = n0;
  for (const s of surfaces) {
    const c = Number.isFinite(s.R) ? 1 / s.R : 0, P = [p[0], p[1], p[2] - s.z];
    const B = d[2] - c * dot(d, P), Cq = c * dot(P, P) - 2 * P[2], disc = B * B - c * Cq;
    if (disc < 0) return { points, p, d, n, status: 'missed' };
    const den = B + Math.sign(B || 1) * Math.sqrt(disc);
    if (den === 0) return { points, p, d, n, status: 'missed' };
    const t = Cq / den, X = [P[0] + t * d[0], P[1] + t * d[1], P[2] + t * d[2]];
    const hit = [X[0], X[1], X[2] + s.z];
    if (s.aperture !== undefined && Math.hypot(X[0], X[1]) > s.aperture) return { points, p, d, n, status: 'clipped' };
    const N = [-c * X[0], -c * X[1], 1 - c * X[2]];
    const dn = refract(d, N, n, s.n);
    points.push(hit);
    if (!dn) return { points, p: hit, d, n, status: 'tir' };
    p = hit; d = dn; n = s.n;
  }
  return { points, p, d, n, status: 'ok' };
}

/**
 * Paraxial ABCD matrix (real-angle convention) of the same surface list from the first to the last
 * vertex, with { efl = −1/C (rear effective focal length), bfd = −A/C (back focal distance from the last
 * vertex) }; C = 0 (afocal) gives efl = bfd = Infinity.
 */
export function paraxialSystem(surfaces, n0 = 1) {
  checkSurfaces(surfaces, n0, 'paraxialSystem');
  const els = [];
  let n = n0, z = surfaces[0].z;
  for (const s of surfaces) {
    if (s.z > z) els.push(abcd.freeSpace(s.z - z));
    els.push(Number.isFinite(s.R) ? abcd.curvedInterface(n, s.n, s.R) : abcd.flatInterface(n, s.n));
    n = s.n; z = s.z;
  }
  const M = abcdSystem(...els), [[A], [Cm]] = M;
  return { matrix: M, efl: Cm === 0 ? Infinity : -1 / Cm, bfd: Cm === 0 ? Infinity : -A / Cm };
}

// ---------------------------------------------------------------- dispersion and rainbows

/**
 * Sellmeier glass data: n²(λ) − 1 = Σ B_i λ²/(λ² − C_i), λ in µm, C_i in µm².
 * BK7: Schott N-BK7 datasheet (n_d = 1.5168, ν_d = 64.17). fusedSilica: Malitson, JOSA 55, 1205 (1965)
 * (n_d = 1.4585), valid 0.21–3.71 µm.
 */
export const GLASS = Object.freeze({
  BK7: Object.freeze({ B: [1.03961212, 0.231792344, 1.01046945], C: [0.00600069867, 0.0200179144, 103.560653], range: [0.3e-6, 2.5e-6], source: 'Schott N-BK7 datasheet' }),
  fusedSilica: Object.freeze({ B: [0.6961663, 0.4079426, 0.8974794], C: [0.0684043 ** 2, 0.1162414 ** 2, 9.896161 ** 2], range: [0.21e-6, 3.71e-6], source: 'Malitson, JOSA 55, 1205 (1965)' }),
});

/**
 * Refractive index n(λ) from Sellmeier coefficients { B, C } (C in µm², see GLASS), λ the vacuum
 * wavelength in metres. Throws outside the coefficient set's stated range (if given) or at a pole.
 */
export function sellmeier(lambda, glass) {
  pos(lambda, 'lambda', 'sellmeier');
  if (!glass?.B?.length || glass.B.length !== glass.C?.length) throw new Error('sellmeier: need { B: [...], C: [...] } of equal length');
  if (glass.range && (lambda < glass.range[0] || lambda > glass.range[1])) throw new Error(`sellmeier: λ = ${lambda} m outside the fitted range ${glass.range.join('–')} m`);
  const l2 = (lambda * 1e6) ** 2;
  let n2 = 1;
  for (let i = 0; i < glass.B.length; i++) n2 += glass.B[i] * l2 / (l2 - glass.C[i]);
  if (!(n2 > 0 && Number.isFinite(n2))) throw new Error('sellmeier: λ at or past an absorption pole');
  return Math.sqrt(n2);
}

/** Abbe number ν_d = (n_d − 1)/(n_F − n_C) with d = 587.5618 nm, F = 486.1327 nm, C = 656.2725 nm. */
export function abbeNumber(glass) {
  const nd = sellmeier(587.5618e-9, glass), nF = sellmeier(486.1327e-9, glass), nC = sellmeier(656.2725e-9, glass);
  return (nd - 1) / (nF - nC);
}

/**
 * Descartes rainbow of order k (k internal reflections, k ≥ 1) for a sphere of relative index n > 1.
 * Deviation D(i) = 2(i − r) + k(π − 2r), sin i = n sin r; minimum at cos²i = (n² − 1)/(k(k + 2)).
 * Returns { incidence, refraction, deviation (rad, minimum), angle = |D − π| reduced to [0, π]
 * (angle from the antisolar point: ≈ 42.0° primary, 50.9° secondary for n = 1.333) }. Geometric
 * optics only: supernumerary bows and the Airy shift need wave theory.
 */
export function rainbow(n, k = 1) {
  pos(n, 'n', 'rainbow');
  if (!(n > 1)) throw new Error('rainbow: n must exceed 1');
  if (!(Number.isInteger(k) && k >= 1)) throw new Error('rainbow: order k must be an integer ≥ 1');
  const c2 = (n * n - 1) / (k * (k + 2));
  if (c2 > 1) throw new Error('rainbow: no stationary deviation for this n and order');
  const i = Math.acos(Math.sqrt(c2)), r = Math.asin(Math.sin(i) / n), D = 2 * (i - r) + k * (Math.PI - 2 * r);
  let a = (D - Math.PI) % (2 * Math.PI);
  if (a < 0) a += 2 * Math.PI;
  if (a > Math.PI) a = 2 * Math.PI - a;
  return { incidence: i, refraction: r, deviation: D, angle: a };
}

// ---------------------------------------------------------------- multilayers (transfer matrix)

/**
 * Thin-film stack by the equivalent scattering recursion (Born & Wolf §1.6; Macleod,
 * Thin-Film Optical Filters ch. 2), using decaying propagation factors to avoid overflow in opaque
 * films. Light from a semi-infinite real medium n0 at angle theta0 (rad)
 * through layers [{ n, d }] (n real or complex { re, im } with im = κ ≥ 0 absorbing, d physical thickness
 * in m) into the substrate ns (real or complex), vacuum wavelength λ (m). polarization 's' | 'p' |
 * 'unpolarized'. Returns { r, t (complex, Hecht conventions as fresnelCoefficients; omitted for
 * 'unpolarized'), R, T, A = 1 − R − T }. T is 0 when the substrate wave is evanescent; lossless stacks
 * give R + T = 1 to round-off.
 */
export function thinFilm({ lambda, n0 = 1, layers = [], ns, theta0 = 0, polarization = 's' }) {
  const who = 'thinFilm';
  pos(lambda, 'lambda', who); pos(n0, 'n0', who);
  if (!(theta0 >= 0 && theta0 < Math.PI / 2)) throw new Error('thinFilm: theta0 must lie in [0, π/2)');
  if (polarization === 'unpolarized') {
    const s = thinFilm({ lambda, n0, layers, ns, theta0, polarization: 's' }), p = thinFilm({ lambda, n0, layers, ns, theta0, polarization: 'p' });
    return { R: (s.R + p.R) / 2, T: (s.T + p.T) / 2, A: (s.A + p.A) / 2 };
  }
  if (polarization !== 's' && polarization !== 'p') throw new Error("thinFilm: polarization must be 's', 'p' or 'unpolarized'");
  if (!Array.isArray(layers)) throw new Error('thinFilm: layers must be an array');
  const s0 = n0 * Math.sin(theta0), c0 = Math.cos(theta0);
  // Work in Macleod's convention (e^{i(ωt − kz)}, N = n − iκ) and conjugate the amplitudes at the end.
  const medium = (v) => {
    const N = cconj(toC(v, who));
    if (N.re <= 0 || N.im > 0) throw new Error('thinFilm: indices need Re n > 0 and κ ≥ 0');
    let Nc = csqrt(csub(cmul(N, N), cx(s0 * s0)));                 // N cosθ
    if (Math.abs(Nc.re) < 1e-15 * Math.hypot(Nc.re, Nc.im) && Nc.im > 0) Nc = cx(-Nc.re, -Nc.im);
    const eta = polarization === 's' ? Nc : cdiv(cmul(N, N), Nc);
    return { N, Nc, eta };
  };
  const eta0 = polarization === 's' ? cx(n0 * c0) : cx(n0 / c0);
  const media = [{ eta: eta0 }], phases = [cx(1)];
  for (const layer of layers) {
    if (!(Number.isFinite(layer?.d) && layer.d >= 0)) throw new Error('thinFilm: layer thickness d must be finite and ≥ 0');
    const { Nc, eta } = medium(layer.n), delta = cscale(Nc, 2 * Math.PI * layer.d / lambda);
    media.push({ eta });
    phases.push(cexp(cx(delta.im, -delta.re)));  // exp(−iδ), with Im δ ≤ 0: magnitude ≤ 1
  }
  const sub = medium(ns);
  media.push(sub); phases.push(cx(1));
  let rM = cx(0), tM = cx(1);
  for (let i = media.length - 2; i >= 0; i--) {
    const e = media[i].eta, en = media[i + 1].eta, sum = cadd(e, en), phase = phases[i + 1];
    const ri = cdiv(csub(e, en), sum), ti = cdiv(cscale(e, 2), sum);
    const back = cmul(rM, cmul(phase, phase)), den = cadd(cx(1), cmul(ri, back));
    tM = cdiv(cmul(cmul(ti, phase), tM), den);
    rM = cdiv(cadd(ri, back), den);
  }
  const R = cabs2(rM), T = sub.eta.re / eta0.re * cabs2(tM);
  let r = cconj(rM), t = cconj(tM);
  if (polarization === 'p') {
    r = cscale(r, -1);
    const cs = cconj(cdiv(sub.Nc, sub.N));                            // cosθ_s (physics convention)
    t = cmul(t, cdiv(cx(c0), cs));
  }
  return { r, t, R, T, A: 1 - R - T };
}

// ---------------------------------------------------------------- polarisation

/**
 * Jones calculus. Vectors [Ex, Ey] and 2×2 matrices hold complex { re, im }. Angles in rad from x.
 * Convention: e^{i(kz − ωt)}; [1, i]/√2 rotates counter-clockwise for an observer facing the source,
 * S3 = 2 Im(Ex* Ey) = +1 ("left-circular" in Born & Wolf / Hecht, "right-hand" in IEEE usage).
 * Angles: polarizer θ = transmission axis; waveplate / halfWave / quarterWave θ = fast axis, the slow
 * axis (θ + 90°) is delayed by the retardance Γ (factor e^{+iΓ}); rotator φ rotates counter-clockwise.
 */
export const jones = Object.freeze({
  linear: theta => [cx(Math.cos(theta)), cx(Math.sin(theta))],
  circular: (sign = 1) => { if (sign !== 1 && sign !== -1) throw new Error('jones.circular: sign must be ±1'); return [cx(Math.SQRT1_2), cx(0, sign * Math.SQRT1_2)]; },
  polarizer: theta => {
    const c = Math.cos(theta), s = Math.sin(theta);
    return [[cx(c * c), cx(c * s)], [cx(c * s), cx(s * s)]];
  },
  rotator: phi => [[cx(Math.cos(phi)), cx(-Math.sin(phi))], [cx(Math.sin(phi)), cx(Math.cos(phi))]],
  waveplate: (retardance, theta) => {
    const c = Math.cos(theta), s = Math.sin(theta), e = cexp(cx(0, retardance));   // fast axis phase 1, slow axis e^{iΓ}
    const one = cx(1);
    const a = cadd(cscale(one, c * c), cscale(e, s * s)), d = cadd(cscale(one, s * s), cscale(e, c * c));
    const b = cscale(csub(one, e), c * s);
    return [[a, b], [b, d]];
  },
  halfWave: theta => jones.waveplate(Math.PI, theta),
  quarterWave: theta => jones.waveplate(Math.PI / 2, theta),
  /** Product of matrices in the order light meets them. */
  system: (...Ms) => Ms.reduce((acc, M) => [[cadd(cmul(M[0][0], acc[0][0]), cmul(M[0][1], acc[1][0])), cadd(cmul(M[0][0], acc[0][1]), cmul(M[0][1], acc[1][1]))],
    [cadd(cmul(M[1][0], acc[0][0]), cmul(M[1][1], acc[1][0])), cadd(cmul(M[1][0], acc[0][1]), cmul(M[1][1], acc[1][1]))]], [[cx(1), cx(0)], [cx(0), cx(1)]]),
  apply: (M, v) => [cadd(cmul(M[0][0], v[0]), cmul(M[0][1], v[1])), cadd(cmul(M[1][0], v[0]), cmul(M[1][1], v[1]))],
  intensity: v => cabs2(v[0]) + cabs2(v[1]),
  /** Stokes [S0, S1, S2, S3] of a Jones vector (fully polarised: S0² = S1² + S2² + S3²). */
  stokes: v => {
    const [a, b] = v, ab = cmul(cconj(a), b);
    return [cabs2(a) + cabs2(b), cabs2(a) - cabs2(b), 2 * ab.re, 2 * ab.im];
  },
});

// ---------------------------------------------------------------- FDTD

/**
 * Numerical dispersion of the 2D Yee scheme: ω with (n/(cΔt))² sin²(ωΔt/2) = sin²(kxΔx/2)/Δx² + sin²(kyΔy/2)/Δy²
 * (Taflove & Hagness, Computational Electrodynamics, eq. 4.17). Units SI (rad/m, m, s → rad/s).
 * Throws when no real ω exists (Courant limit exceeded).
 */
export function yeeDispersion({ kx, ky = 0, dx, dy = dx, dt, n = 1 }) {
  pos(dx, 'dx', 'yeeDispersion'); pos(dy, 'dy', 'yeeDispersion'); pos(dt, 'dt', 'yeeDispersion'); pos(n, 'n', 'yeeDispersion');
  fin(kx, 'kx', 'yeeDispersion'); fin(ky, 'ky', 'yeeDispersion');
  const s = (C0 * dt / n) * Math.sqrt((Math.sin(kx * dx / 2) / dx) ** 2 + (Math.sin(ky * dy / 2) / dy) ** 2);
  if (s > 1) throw new Error('yeeDispersion: no real frequency (Courant limit exceeded)');
  return 2 * Math.asin(s) / dt;
}

const WALLS = ['pec', 'periodic', 'mur'];

/**
 * 2D finite-difference time-domain solver (Yee 1966) on an nx × ny grid (index i + nx·j), spacing
 * dx, dy (m), SI fields. mode 'TM': Ez (V/m) at nodes (i, j), Hx at (i, j+½), Hy at (i+½, j).
 * mode 'TE': Hz (A/m) at (i+½, j+½), Ex at (i+½, j), Ey at (i, j+½). epsR: relative permittivity per
 * cell (Float64Array nx·ny, default 1; TE uses the average of the two adjacent cells). μ = μ0.
 * boundary: 'pec' | 'periodic' | 'mur' or { x, y } per axis. Mur is first order (Mur 1981): exact
 * only at normal incidence with cΔt = Δx, reflecting a few % at oblique incidence; size domains so
 * that this is negligible or say so. dt (s) defaults to courant × the stability limit
 * Δt_max = √(min εr)/(c√(1/Δx² + 1/Δy²)); dt > Δt_max throws.
 * sources: [{ i, j, f: t => value, hard? }] add f(t) (V/m in TM to Ez, A/m in TE to Hz) at node
 * (i, j) after each update (hard: overwrite); i or j may be 'all' for a line source.
 * Returns { Ez, Hx, Hy | Hz, Ex, Ey (live Float64Arrays, writable for initial conditions), dt, dtMax,
 * t, steps, step(k = 1), advanceTo(t), energy() }. energy() is the discrete invariant
 * ½Σ(ε E^{n}·E^{n} + μ H^{n−½}·H^{n+½}) dx dy (J/m): exactly constant (round-off) for PEC or periodic
 * walls without sources.
 */
export function fdtd2D({ nx, ny, dx, dy = dx, dt, courant = 0.5, mode = 'TM', epsR, boundary = 'mur', sources = [] }) {
  const who = 'fdtd2D';
  if (!(Number.isInteger(nx) && nx >= 2 && Number.isInteger(ny) && ny >= 2)) throw new Error('fdtd2D: nx, ny must be integers ≥ 2');
  pos(dx, 'dx', who); pos(dy, 'dy', who);
  if (mode !== 'TM' && mode !== 'TE') throw new Error("fdtd2D: mode must be 'TM' or 'TE'");
  const size = nx * ny, eps = epsR ? Float64Array.from(epsR) : new Float64Array(size).fill(1);
  if (eps.length !== size || !eps.every(v => Number.isFinite(v) && v > 0)) throw new Error('fdtd2D: epsR must have nx·ny positive finite entries');
  const bx = typeof boundary === 'string' ? boundary : boundary?.x, by = typeof boundary === 'string' ? boundary : boundary?.y;
  if (!WALLS.includes(bx) || !WALLS.includes(by)) throw new Error("fdtd2D: boundary must be 'pec', 'periodic' or 'mur' (or { x, y })");
  let minEps = Infinity;
  for (const v of eps) if (v < minEps) minEps = v;
  const dtMax = Math.sqrt(minEps) / (C0 * Math.sqrt(1 / dx ** 2 + 1 / dy ** 2));
  const step0 = dt ?? courant * dtMax;
  if (!(Number.isFinite(step0) && step0 > 0)) throw new Error('fdtd2D: dt (or courant) must be positive');
  if (step0 > dtMax) throw new Error(`fdtd2D: Courant condition violated (dt = ${step0} s > dt_max = ${dtMax} s)`);
  const srcs = sources.map(s => {
    if (typeof s?.f !== 'function') throw new Error('fdtd2D: each source needs f(t)');
    const is = s.i === 'all' ? [...Array(nx).keys()] : [s.i], js = s.j === 'all' ? [...Array(ny).keys()] : [s.j];
    if (![...is].every(i => Number.isInteger(i) && i >= 0 && i < nx) || ![...js].every(j => Number.isInteger(j) && j >= 0 && j < ny)) throw new Error('fdtd2D: source position outside the grid');
    const cells = [];
    for (const j of js) for (const i of is) cells.push(i + nx * j);
    return { cells, f: s.f, hard: !!s.hard };
  });
  const T = step0, px = bx === 'periodic', py = by === 'periodic';
  const ip = i => (i + 1 < nx ? i + 1 : (px ? 0 : -1)), im = i => (i > 0 ? i - 1 : (px ? nx - 1 : -1));
  const jp = j => (j + 1 < ny ? j + 1 : (py ? 0 : -1)), jm = j => (j > 0 ? j - 1 : (py ? ny - 1 : -1));
  const dA = dx * dy;
  const sim = { dt: T, dtMax, t: 0, steps: 0 };

  if (mode === 'TM') {
    const Ez = new Float64Array(size), Hx = new Float64Array(size), Hy = new Float64Array(size);
    const ce = Float64Array.from(eps, e => T / (eps0 * e)), ch = T / mu0;
    const hUpdate = (hx, hy) => {
      for (let j = 0; j < ny; j++) {
        const j1 = jp(j);
        for (let i = 0; i < nx; i++) {
          const k = i + nx * j, i1 = ip(i);
          if (j1 >= 0) hx[k] -= ch / dy * (Ez[i + nx * j1] - Ez[k]);
          if (i1 >= 0) hy[k] += ch / dx * (Ez[i1 + nx * j] - Ez[k]);
        }
      }
    };
    const murCoef = k => { const v = C0 / Math.sqrt(eps[k]); return [(v * T - dx) / (v * T + dx), (v * T - dy) / (v * T + dy)]; };
    sim.Ez = Ez; sim.Hx = Hx; sim.Hy = Hy;
    sim.step = (count = 1) => {
      for (let s = 0; s < count; s++) {
        hUpdate(Hx, Hy);
        const old = Ez.slice();
        for (let j = 0; j < ny; j++) {
          const j0 = jm(j);
          for (let i = 0; i < nx; i++) {
            const i0 = im(i), k = i + nx * j;
            if (i0 < 0 || j0 < 0 || ip(i) < 0 || jp(j) < 0) continue;      // wall node
            Ez[k] += ce[k] * ((Hy[k] - Hy[i0 + nx * j]) / dx - (Hx[k] - Hx[i + nx * j0]) / dy);
          }
        }
        wallsE(Ez, old, 'x', 'y', murCoef);
        sim.steps++; sim.t = sim.steps * T;
        for (const src of srcs) { const v = src.f(sim.t); for (const k of src.cells) Ez[k] = src.hard ? v : Ez[k] + v; }
      }
      return sim;
    };
    sim.energy = () => {
      const hx = Hx.slice(), hy = Hy.slice();
      hUpdate(hx, hy);
      let w = 0;
      for (let k = 0; k < size; k++) w += eps0 * eps[k] * Ez[k] * Ez[k] + mu0 * (Hx[k] * hx[k] + Hy[k] * hy[k]);
      return w * dA / 2;
    };
  } else {
    const Hz = new Float64Array(size), Ex = new Float64Array(size), Ey = new Float64Array(size);
    const ex = new Float64Array(size), ey = new Float64Array(size);
    for (let j = 0; j < ny; j++) { for (let i = 0; i < nx; i++) {
      const k = i + nx * j, i1 = ip(i), j1 = jp(j);
      ex[k] = T / (eps0 * (i1 >= 0 ? (eps[k] + eps[i1 + nx * j]) / 2 : eps[k]));
      ey[k] = T / (eps0 * (j1 >= 0 ? (eps[k] + eps[i + nx * j1]) / 2 : eps[k]));
    } }
    const ch = T / mu0;
    const hUpdate = hz => {
      for (let j = 0; j < ny; j++) {
        const j1 = jp(j);
        if (j1 < 0) continue;
        for (let i = 0; i < nx; i++) {
          const i1 = ip(i);
          if (i1 < 0) continue;
          const k = i + nx * j;
          hz[k] -= ch * ((Ey[i1 + nx * j] - Ey[k]) / dx - (Ex[i + nx * j1] - Ex[k]) / dy);
        }
      }
    };
    const murCoef = k => { const v = C0 / Math.sqrt(eps[k]); return [(v * T - dx) / (v * T + dx), (v * T - dy) / (v * T + dy)]; };
    sim.Hz = Hz; sim.Ex = Ex; sim.Ey = Ey;
    sim.step = (count = 1) => {
      for (let s = 0; s < count; s++) {
        hUpdate(Hz);
        const oldX = Ex.slice(), oldY = Ey.slice();
        for (let j = 0; j < ny; j++) {
          const j0 = jm(j), j1 = jp(j);
          for (let i = 0; i < nx; i++) {
            const i0 = im(i), i1 = ip(i), k = i + nx * j;
            if (i1 >= 0 && j0 >= 0 && j1 >= 0) Ex[k] += ex[k] * (Hz[k] - Hz[i + nx * j0]) / dy;   // Ex at (i+½, j): y-walls hold it
            if (j1 >= 0 && i0 >= 0 && i1 >= 0) Ey[k] -= ey[k] * (Hz[k] - Hz[i0 + nx * j]) / dx;   // Ey at (i, j+½): x-walls hold it
          }
        }
        wallsE(Ey, oldY, 'x', null, murCoef);
        wallsE(Ex, oldX, null, 'y', murCoef);
        sim.steps++; sim.t = sim.steps * T;
        for (const src of srcs) { const v = src.f(sim.t); for (const k of src.cells) Hz[k] = src.hard ? v : Hz[k] + v; }
      }
      return sim;
    };
    sim.energy = () => {
      const hz = Hz.slice();
      hUpdate(hz);
      let w = 0;
      for (let k = 0; k < size; k++) w += T * (Ex[k] * Ex[k] / ex[k] + Ey[k] * Ey[k] / ey[k]) + mu0 * Hz[k] * hz[k];
      return w * dA / 2;
    };
  }

  /** Wall treatment for a tangential E component F (new values) with previous values old. */
  function wallsE(F, old, xAxis, yAxis, murCoef) {
    if (xAxis && !px) { for (let j = 0; j < ny; j++) {
      for (const [b, nb] of [[0, 1], [nx - 1, nx - 2]]) {
        const k = b + nx * j, kn = nb + nx * j;
        F[k] = bx === 'pec' ? 0 : old[kn] + murCoef(k)[0] * (F[kn] - old[k]);
      }
    } }
    if (yAxis && !py) { for (let i = 0; i < nx; i++) {
      for (const [b, nb] of [[0, 1], [ny - 1, ny - 2]]) {
        const k = i + nx * b, kn = i + nx * nb;
        F[k] = by === 'pec' ? 0 : old[kn] + murCoef(k)[1] * (F[kn] - old[k]);
      }
    } }
  }

  sim.advanceTo = (time) => {
    fin(time, 'time', 'fdtd2D.advanceTo');
    const target = Math.round(time / T);
    if (Math.abs(target * T - time) > 1e-9 * Math.max(T, Math.abs(time))) throw new Error('fdtd2D.advanceTo: time must be a whole number of steps');
    if (target < sim.steps) throw new Error('fdtd2D.advanceTo: cannot go back in time');
    return sim.step(target - sim.steps);
  };
  return sim;
}
