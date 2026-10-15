// Multidimensional regular-grid mathematics. Axis 0 is fastest:
// shape [nx, ny, nz] means index = i + nx * (j + ny * k).
// Fields are flat arrays. Spacing is a positive scalar or one value per axis.
// Irregular domains: grid.mask (1 = inside, length = grid size) restricts the derivative operators
// to the masked cells; maskBc 'dirichlet' (zero outside, e.g. a clamped edge) or 'neumann' (zero
// flux). The boundary is a staircase of cells: first-order in dx for curved boundaries.

import { fft } from './numerics.mjs';

const finite = (value, who) => {
  if (!Number.isFinite(value)) throw new Error(`${who}: numerical range exceeded; use finite inputs and suitable units`);
  return value;
};
const positive = (value, who) => {
  if (!(Number.isFinite(value) && value > 0)) throw new Error(`${who}: expected a positive finite parameter`);
};
function grid(shape, spacing = 1, bc = 'periodic', mask = null, maskBc = 'dirichlet') {
  if (!shape?.length) throw new Error('fields: shape must contain positive integer extents');
  const dimensions = Array.from(shape);
  let size = 1;
  const strides = dimensions.map(length => {
    if (!(Number.isSafeInteger(length) && length > 0)) throw new Error('fields: shape must contain positive integer extents');
    const stride = size; size *= length;
    if (!Number.isSafeInteger(size)) throw new Error('fields: grid size exceeds the integer range');
    return stride;
  });
  const dx = typeof spacing === 'number' ? dimensions.map(() => spacing) : Array.from(spacing);
  if (dx.length !== dimensions.length || !dx.every(value => Number.isFinite(value) && value > 0)) throw new Error('fields: spacing must be positive and match the shape');
  const boundaries = typeof bc === 'string' ? dimensions.map(() => bc) : Array.from(bc);
  if (boundaries.length !== dimensions.length || !boundaries.every(value => ['periodic', 'dirichlet', 'neumann'].includes(value))) throw new Error('fields: boundary must be periodic, dirichlet or neumann for each axis');
  if (mask && mask.length !== size) throw new Error('fields: mask length must equal the grid size');
  if (!['dirichlet', 'neumann'].includes(maskBc)) throw new Error('fields: maskBc must be dirichlet or neumann');
  return { shape: dimensions, spacing: dx, strides, boundaries, size, mask, maskBc };
}
function field(values, size, who) {
  if (!values || values.length !== size) throw new Error(`${who}: field length must equal the grid size`);
  for (const value of values) finite(value, who);
}
function overlaps(a, b) {
  return a === b || (ArrayBuffer.isView(a) && ArrayBuffer.isView(b) && a.buffer === b.buffer && a.byteOffset < b.byteOffset + b.byteLength && b.byteOffset < a.byteOffset + a.byteLength);
}
function output(out, input, who) {
  if (!out || out.length !== input.length || overlaps(out, input)) throw new Error(`${who}: output must have matching length and not overlap the input`);
}
function neighbor(values, index, axis, side, spec) {
  const stride = spec.strides[axis], length = spec.shape[axis], coordinate = Math.floor(index / stride) % length;
  let other;
  if (coordinate + side >= 0 && coordinate + side < length) other = index + side * stride;
  else if (spec.boundaries[axis] === 'periodic') other = index - side * (length - 1) * stride;
  else return spec.boundaries[axis] === 'neumann' ? values[index] : 0;
  if (spec.mask && !spec.mask[other]) return spec.maskBc === 'neumann' ? values[index] : 0;
  return values[other];
}

/**
 * Compact second-order Laplacian on any-dimensional grid. Dirichlet: zero ghost values;
 * Neumann: copied edge-cell ghosts (zero face flux), conserving the unweighted cell sum.
 * The output must not alias input. Mixed per-axis boundaries are allowed.
 */
export function laplacianND(values, { shape, spacing = 1, bc = 'periodic', mask, maskBc }, out = new Float64Array(values.length)) {
  const spec = grid(shape, spacing, bc, mask, maskBc);
  field(values, spec.size, 'laplacianND'); output(out, values, 'laplacianND');
  for (let index = 0; index < spec.size; index++) {
    if (mask && !mask[index]) { out[index] = 0; continue; }
    let result = 0;
    for (let axis = 0; axis < spec.shape.length; axis++) {
      const dx = spec.spacing[axis];
      result += ((neighbor(values, index, axis, -1, spec) - values[index]) / dx + (neighbor(values, index, axis, 1, spec) - values[index]) / dx) / dx;
    }
    out[index] = finite(result, 'laplacianND');
  }
  return out;
}

/** Central second-order gradient, one flat component per axis, with the same ghost rules. */
export function gradientND(values, { shape, spacing = 1, bc = 'periodic', mask, maskBc }) {
  const spec = grid(shape, spacing, bc, mask, maskBc);
  field(values, spec.size, 'gradientND');
  return spec.shape.map((_length, axis) => Float64Array.from({ length: spec.size }, (_, index) => (mask && !mask[index]) ? 0 :
    finite((neighbor(values, index, axis, 1, spec) - neighbor(values, index, axis, -1, spec)) / spec.spacing[axis] / 2, 'gradientND')));
}

/** Central divergence. div(grad) uses a wider stencil than the compact laplacianND. */
export function divergenceND(components, { shape, spacing = 1, bc = 'periodic', mask, maskBc }) {
  const spec = grid(shape, spacing, bc, mask, maskBc);
  if (components?.length !== spec.shape.length) throw new Error('divergenceND: one component per axis required');
  for (const values of components) field(values, spec.size, 'divergenceND');
  return Float64Array.from({ length: spec.size }, (_, index) => {
    if (mask && !mask[index]) return 0;
    let result = 0;
    for (let axis = 0; axis < spec.shape.length; axis++) result += (neighbor(components[axis], index, axis, 1, spec) - neighbor(components[axis], index, axis, -1, spec)) / spec.spacing[axis] / 2;
    return finite(result, 'divergenceND');
  });
}

function fftGrid(shape) {
  const spec = grid(shape);
  if (!spec.shape.every(length => length <= 2 ** 29 && Number.isInteger(Math.log2(length)))) throw new Error('fftND: every extent must be a supported power of two');
  return spec;
}
const rows = new Map();
/** Separable FFT over every axis; callers check finiteness once afterwards (not per element). */
function transform(re, im, spec, inverse, axes = spec.shape.map((_, k) => k)) {
  for (const axis of axes) {
    const length = spec.shape[axis], stride = spec.strides[axis], block = length * stride;
    if (!rows.has(length)) rows.set(length, [new Float64Array(length), new Float64Array(length)]);
    const [rowRe, rowIm] = rows.get(length);
    if (stride === 1) {
      for (let base = 0; base < spec.size; base += block) {
        const r = re.subarray(base, base + length), m = im.subarray(base, base + length);
        fft(r, m, inverse);
      }
      continue;
    }
    for (let base = 0; base < spec.size; base += block) { for (let offset = 0; offset < stride; offset++) {
      for (let i = 0; i < length; i++) { const index = base + offset + i * stride; rowRe[i] = re[index]; rowIm[i] = im[index]; }
      fft(rowRe, rowIm, inverse);
      for (let i = 0; i < length; i++) { const index = base + offset + i * stride; re[index] = rowRe[i]; im[index] = rowIm[i]; }
    } }
  }
}

/** In-place separable radix-2 complex FFT; inverse includes 1/product(shape). */
export function fftND(re, im, shape, inverse = false) {
  const spec = fftGrid(shape);
  field(re, spec.size, 'fftND'); field(im, spec.size, 'fftND');
  if (overlaps(re, im)) throw new Error('fftND: real and imaginary storage must not overlap');
  transform(re, im, spec, inverse);
  field(re, spec.size, 'fftND'); field(im, spec.size, 'fftND');
  return { re, im };
}

function eigenvalues(spec, discretization) {
  if (!['spectral', 'finite-difference'].includes(discretization)) throw new Error('fields: discretization must be spectral or finite-difference');
  return Float64Array.from({ length: spec.size }, (_, index) => {
    let value = 0;
    for (let axis = 0; axis < spec.shape.length; axis++) {
      const n = spec.shape[axis], i = Math.floor(index / spec.strides[axis]) % n, mode = i < n / 2 ? i : i - n;
      const k = discretization === 'spectral' ? 2 * Math.PI * mode / n / spec.spacing[axis] : 2 * Math.sin(Math.PI * mode / n) / spec.spacing[axis];
      value += k * k;
    }
    return finite(value, 'fields Fourier operator');
  });
}
/** Kinetic energy ħ²Σ k_a²/(2 m_a) of every Fourier mode; mass is a number or one value per axis. */
function kineticEnergies(spec, mass, hbar) {
  const m = typeof mass === 'number' ? spec.shape.map(() => mass) : Array.from(mass);
  if (m.length !== spec.shape.length || !m.every(v => v > 0)) throw new Error('fields: mass must be positive (a number or one per axis)');
  return Float64Array.from({ length: spec.size }, (_, index) => {
    let e = 0;
    for (let axis = 0; axis < spec.shape.length; axis++) {
      const n = spec.shape[axis], i = Math.floor(index / spec.strides[axis]) % n, mode = i < n / 2 ? i : i - n;
      const k = 2 * Math.PI * mode / n / spec.spacing[axis];
      e += hbar * hbar * k * k / (2 * m[axis]);
    }
    return e;
  });
}
function mean(values) {
  let result = 0, correction = 0;
  for (const value of values) {
    const term = value / values.length - correction, next = result + term;
    correction = (next - result) - term; result = next;
  }
  return result;
}
function volume(spec) {
  const value = spec.spacing.reduce((product, dx) => product * dx, 1);
  positive(value, 'fields cell volume'); return value;
}

/**
 * Solve (-Laplacian + screening^2) u = source on a periodic power-of-two grid.
 * Spectral (default) or compact finite-difference symbol. Unscreened Poisson requires a
 * zero-mean source (rtol relative to maximum source magnitude); fixes the solution mean to zero.
 * Returns { u, sourceMean, discretization }. No incompatible mean is silently removed.
 */
export function solvePoissonPeriodic(source, { shape, spacing = 1, screening = 0, discretization = 'spectral', rtol = 1e-12, atol = 0 }) {
  const spec = grid(shape, spacing); fftGrid(shape); field(source, spec.size, 'solvePoissonPeriodic');
  if (!(Number.isFinite(screening) && screening >= 0 && Number.isFinite(rtol) && rtol >= 0 && Number.isFinite(atol) && atol >= 0)) throw new Error('solvePoissonPeriodic: invalid screening or tolerance');
  const screen2 = screening * screening;
  if (!Number.isFinite(screen2) || (screening > 0 && screen2 === 0)) throw new Error('solvePoissonPeriodic: screening scale outside numerical range');
  let scale = 0;
  for (const value of source) scale = Math.max(scale, Math.abs(value));
  const divisor = scale || 1, re = Float64Array.from(source, value => value / divisor), im = new Float64Array(spec.size);
  const scaledMean = mean(re), sourceMean = finite(scaledMean * divisor, 'solvePoissonPeriodic');
  if (screening === 0 && Math.abs(scaledMean) > rtol + atol / divisor) throw new Error('solvePoissonPeriodic: unscreened periodic source must have zero mean');
  const lambda = eigenvalues(spec, discretization);
  transform(re, im, spec, false);
  for (let i = 0; i < spec.size; i++) {
    const denominator = finite(lambda[i] + screen2, 'solvePoissonPeriodic');
    if (denominator === 0) re[i] = im[i] = 0;
    else { re[i] = finite(re[i] / denominator, 'solvePoissonPeriodic'); im[i] = finite(im[i] / denominator, 'solvePoissonPeriodic'); }
  }
  transform(re, im, spec, true);
  return { u: Float64Array.from(re, value => finite(value * divisor, 'solvePoissonPeriodic')), sourceMean, discretization };
}

/**
 * Linear periodic diffusion du/dt = D Laplacian u by exact Fourier-mode propagation.
 * Spectral or compact finite-difference spatial operator, no source. advanceTo evaluates
 * from the original coefficients, so results at a given time do not depend on frame schedules.
 */
export function diffusionPeriodic({ u0, shape, spacing = 1, D, dt, t0 = 0, discretization = 'spectral' }) {
  const spec = grid(shape, spacing); fftGrid(shape); field(u0, spec.size, 'diffusionPeriodic');
  if (!(Number.isFinite(D) && D >= 0 && Number.isFinite(t0))) throw new Error('diffusionPeriodic: D must be finite/non-negative and t0 finite');
  positive(dt, 'diffusionPeriodic dt');
  let scale = 0;
  for (const value of u0) scale = Math.max(scale, Math.abs(value));
  const divisor = scale || 1, initialRe = Float64Array.from(u0, value => value / divisor), initialIm = new Float64Array(spec.size);
  transform(initialRe, initialIm, spec, false);
  const lambda = eigenvalues(spec, discretization), u = Float64Array.from(u0);
  const re = new Float64Array(spec.size), im = new Float64Array(spec.size);
  const s = {
    t: t0, u, steps: 0, shape: spec.shape.slice(), discretization,
    advanceTo(target) {
      if (!(Number.isFinite(target) && target >= s.t)) throw new Error('diffusionPeriodic: target must be finite and not earlier than current time');
      const elapsed = finite(target - t0, 'diffusionPeriodic');
      for (let i = 0; i < spec.size; i++) {
        const factor = lambda[i] === 0 || D === 0 ? 1 : Math.exp(-D * lambda[i] * elapsed);
        re[i] = initialRe[i] * factor; im[i] = initialIm[i] * factor;
      }
      transform(re, im, spec, true);
      for (let i = 0; i < spec.size; i++) re[i] = finite(re[i] * divisor, 'diffusionPeriodic');
      u.set(re); s.t = target; s.steps++; return true;
    },
    step() { const target = finite(s.t + dt, 'diffusionPeriodic'); if (target === s.t) throw new Error('diffusionPeriodic: step underflow'); return s.advanceTo(target); },
    total() { return finite(mean(u) * spec.size * volume(spec), 'diffusionPeriodic total'); },
  };
  return s;
}

/**
 * Periodic time-dependent Schrodinger equation in any dimension, static real V and isotropic
 * mass. Strang split FFT, second order in dt, unitary up to round-off. Inputs are copied.
 * Returns t/re/im, step/advanceTo, probability() and norm(); no automatic normalization.
 *
 * absorb: { width, strength } adds a complex absorbing potential −iW near both ends of every axis
 * (width: grid cells per axis or one number; W = strength·s², s = depth into the layer ∈ [0, 1]),
 * so outgoing waves leave instead of wrapping round. The norm then decreases: absorbed() is the
 * probability removed. Choose strength comparable to the packet's kinetic energy and a layer of
 * several wavelengths; check the reflected remainder by changing the layer.
 */
export function schrodingerND({ re0, im0, V, shape, spacing = 1, dt, mass = 1, hbar = 1, t0 = 0, absorb }) {
  const spec = grid(shape, spacing); fftGrid(shape); field(re0, spec.size, 'schrodingerND');
  if (im0) field(im0, spec.size, 'schrodingerND');
  if (V) field(V, spec.size, 'schrodingerND potential');
  positive(dt, 'schrodingerND dt'); positive(hbar, 'schrodingerND hbar'); finite(t0, 'schrodingerND t0');
  const Tk = kineticEnergies(spec, mass, hbar);                // per-axis masses allowed (anisotropic kinetic energy)
  const halfV = Float64Array.from({ length: spec.size }, (_, i) => finite(-(V?.[i] ?? 0) / hbar * (dt / 2), 'schrodingerND potential phase'));
  const kinetic = Float64Array.from(Tk, value => finite(-value / hbar * dt, 'schrodingerND kinetic phase'));
  const re = Float64Array.from(re0), im = im0 ? Float64Array.from(im0) : new Float64Array(spec.size);
  const workRe = new Float64Array(spec.size), workIm = new Float64Array(spec.size), cellVolume = volume(spec);
  let halfDamp = null;                                    // e^{−W dt/2ħ} per cell, applied with each potential half-step
  if (absorb) {
    const widths = typeof absorb.width === 'number' ? spec.shape.map(() => absorb.width) : Array.from(absorb.width ?? []);
    if (widths.length !== spec.shape.length || !widths.every(w => Number.isInteger(w) && w >= 0) || widths.some((w, k) => 2 * w >= spec.shape[k])) throw new Error('schrodingerND: absorb.width must be whole cells per axis, less than half the axis');
    if (!(Number.isFinite(absorb.strength) && absorb.strength > 0)) throw new Error('schrodingerND: absorb.strength must be positive');
    halfDamp = Float64Array.from({ length: spec.size }, (_, index) => {
      let W = 0;
      for (let axis = 0; axis < spec.shape.length; axis++) {
        const w = widths[axis];
        if (!w) continue;
        const n = spec.shape[axis], i = Math.floor(index / spec.strides[axis]) % n, depth = Math.max(w - i, i - (n - 1 - w), 0) / w;
        W += absorb.strength * depth * depth;
      }
      return Math.exp(-W * dt / (2 * hbar));
    });
  }
  // precomputed factors: potential half-step (with the absorbing damping folded in) and kinetic step
  const vC = Float64Array.from(halfV, (p, i) => Math.cos(p) * (halfDamp ? halfDamp[i] : 1));
  const vS = Float64Array.from(halfV, (p, i) => Math.sin(p) * (halfDamp ? halfDamp[i] : 1));
  const kC = Float64Array.from(kinetic, Math.cos), kS = Float64Array.from(kinetic, Math.sin);
  const multiply = (C, S) => {
    for (let i = 0; i < spec.size; i++) {
      const c = C[i], d = S[i], a = workRe[i], b = workIm[i];
      workRe[i] = a * c - b * d; workIm[i] = a * d + b * c;
    }
  };
  const s = {
    t: t0, re, im, steps: 0, shape: spec.shape.slice(),
    step() {
      const nextTime = finite(s.t + dt, 'schrodingerND');
      if (nextTime === s.t) throw new Error('schrodingerND: step underflow');
      workRe.set(re); workIm.set(im);
      multiply(vC, vS); transform(workRe, workIm, spec, false); multiply(kC, kS); transform(workRe, workIm, spec, true); multiply(vC, vS);
      let check = 0;
      for (let i = 0; i < spec.size; i++) check += workRe[i] + workIm[i];
      if (!Number.isFinite(check)) throw new Error(`schrodingerND: the wavefunction became non-finite at t = ${s.t}`);
      re.set(workRe); im.set(workIm); s.t = nextTime; s.steps++;
    },
    advanceTo(target, maxSteps = 1000000) {
      if (!(Number.isFinite(target) && target >= s.t && Number.isSafeInteger(maxSteps) && maxSteps > 0)) throw new Error('schrodingerND: finite future target and positive integer budget required');
      let steps = 0;
      while (s.t < target) { if (steps++ === maxSteps) return false; s.step(); }
      return true;
    },
    probability() { return Float64Array.from(re, (value, i) => finite(value * value + im[i] * im[i], 'schrodingerND probability')); },
    norm() { return finite(mean(s.probability()) * spec.size * cellVolume, 'schrodingerND norm'); },
    /** Probability removed by the absorbing layers so far (0 without absorb). */
    absorbed() { return halfDamp ? Math.max(0, initialNorm - s.norm()) : 0; },
    /**
     * ⟨H⟩ = ⟨T⟩ + ⟨V⟩ per unit norm (kinetic part spectral, exact for the grid). The split-operator
     * scheme keeps it bounded with an O(Δt²) oscillation, not exactly constant: judge it by order.
     */
    energy() {
      workRe.set(re); workIm.set(im); transform(workRe, workIm, spec, false);
      let t = 0, w = 0, v = 0, nrm = 0;
      for (let i = 0; i < spec.size; i++) { const a = workRe[i] * workRe[i] + workIm[i] * workIm[i]; t += a * Tk[i]; w += a; }
      for (let i = 0; i < spec.size; i++) { const a = re[i] * re[i] + im[i] * im[i]; v += a * (V?.[i] ?? 0); nrm += a; }
      return t / w + v / nrm;
    },
  };
  const initialNorm = s.norm();
  return s;
}

/** In-place FFT along one axis only (fields.mjs storage); inverse includes 1/n. For alternating-direction splittings. */
export function fftAxis(re, im, shape, axis, inverse = false) {
  const spec = fftGrid(shape);
  if (!(Number.isInteger(axis) && axis >= 0 && axis < spec.shape.length)) throw new Error('fftAxis: axis out of range');
  field(re, spec.size, 'fftAxis'); field(im, spec.size, 'fftAxis');
  transform(re, im, spec, inverse, [axis]);
  field(re, spec.size, 'fftAxis'); field(im, spec.size, 'fftAxis');
  return { re, im };
}

/**
 * Topological charge of a complex 2D field on every grid plaquette: the phase winding
 * (Σ of wrapped phase differences around the cell)/2π, an integer (±1 at a vortex core, 0
 * elsewhere). Returns { winding: Int8Array (size (nx−1)(ny−1), index i + (nx−1)j), vortices:
 * [{ x, y, charge }] at plaquette centres in grid units }. Cells where |ψ|² < threshold·max are
 * skipped (phase undefined in empty regions).
 */
export function phaseWinding(re, im, { shape, threshold = 1e-3 }) {
  const spec = grid(shape);
  if (spec.shape.length !== 2) throw new Error('phaseWinding: a 2D shape [nx, ny] is required');
  field(re, spec.size, 'phaseWinding'); field(im, spec.size, 'phaseWinding');
  const [nx, ny] = spec.shape, ph = Float64Array.from(re, (r, i) => Math.atan2(im[i], r));
  let maxRho = 0; for (let i = 0; i < spec.size; i++) maxRho = Math.max(maxRho, re[i] * re[i] + im[i] * im[i]);
  const wrap = d => d - 2 * Math.PI * Math.round(d / (2 * Math.PI));
  const winding = new Int8Array((nx - 1) * (ny - 1)), vortices = [];
  for (let j = 0; j < ny - 1; j++) {for (let i = 0; i < nx - 1; i++) {
    const a = i + nx * j, b = a + 1, c = a + 1 + nx, d = a + nx;
    if (Math.min(...[a, b, c, d].map(k => re[k] * re[k] + im[k] * im[k])) < threshold * maxRho) continue;
    const w = Math.round((wrap(ph[b] - ph[a]) + wrap(ph[c] - ph[b]) + wrap(ph[d] - ph[c]) + wrap(ph[a] - ph[d])) / (2 * Math.PI));
    if (w) { winding[i + (nx - 1) * j] = w; vortices.push({ x: i + 0.5, y: j + 0.5, charge: w }); }
  }}
  return { winding, vortices };
}

/**
 * Probability current j = (ħ/m) Im(ψ* ∇ψ) on a periodic power-of-two grid, with the derivative taken
 * spectrally (consistent with schrodingerND's kinetic operator; Nyquist mode dropped). Returns one
 * Float64Array per axis. The flux through a line or surface, integrated over time, measures what
 * leaves a region independently of the norm: use it to check absorbing layers and transmission.
 */
export function probabilityCurrent(re, im, { shape, spacing = 1, mass = 1, hbar = 1 }) {
  const spec = grid(shape, spacing); fftGrid(shape);
  field(re, spec.size, 'probabilityCurrent'); field(im, spec.size, 'probabilityCurrent');
  const out = [];
  for (let axis = 0; axis < spec.shape.length; axis++) {
    const n = spec.shape[axis], stride = spec.strides[axis], k0 = 2 * Math.PI / (n * spec.spacing[axis]);
    const dRe = Float64Array.from(re), dIm = Float64Array.from(im);
    transform(dRe, dIm, spec, false);
    for (let i = 0; i < spec.size; i++) {
      const m = Math.floor(i / stride) % n, k = m === n / 2 ? 0 : (m < n / 2 ? m : m - n) * k0;
      const a = dRe[i], b = dIm[i];
      dRe[i] = -k * b; dIm[i] = k * a;                         // multiply by i k
    }
    transform(dRe, dIm, spec, true);
    out.push(Float64Array.from(re, (r, i) => hbar / mass * (r * dIm[i] - im[i] * dRe[i])));
  }
  return out;
}

/**
 * Level-set contours of a 2D field (storage as above: index i + nx·j) by marching squares with
 * linear interpolation; saddles are resolved by the cell average. Returns polylines in physical
 * coordinates x = origin[0] + i·dx, y = origin[1] + j·dy: [{ points: [[x, y], …], closed }].
 * Accuracy O(dx²) for smooth fields; a contour is a drawing aid, so state the level's meaning.
 */
export function contourLines(values, { shape, spacing = 1, origin = [0, 0] }, level) {
  const spec = grid(shape, spacing);
  if (spec.shape.length !== 2) throw new Error('contourLines: a 2D shape [nx, ny] is required');
  field(values, spec.size, 'contourLines'); finite(level, 'contourLines level');
  const [nx, ny] = spec.shape, [dx, dy] = spec.spacing, v = (i, j) => values[i + nx * j];
  const point = new Map();                                      // edge id → crossing point
  const cross = (i0, j0, i1, j1) => {
    const id = i1 > i0 ? 2 * (i0 + nx * j0) : 2 * (i0 + nx * j0) + 1;
    if (!point.has(id)) {
      const a = v(i0, j0), b = v(i1, j1), f = a === b ? 0.5 : (level - a) / (b - a);
      point.set(id, [origin[0] + (i0 + f * (i1 - i0)) * dx, origin[1] + (j0 + f * (j1 - j0)) * dy]);
    }
    return id;
  };
  const links = new Map(), link = (p, q) => { for (const [a, b] of [[p, q], [q, p]]) { if (!links.has(a)) links.set(a, []); links.get(a).push(b); } };
  for (let j = 0; j < ny - 1; j++) {for (let i = 0; i < nx - 1; i++) {
    const c = [v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1)];
    const code = c.reduce((s, x, k) => s | (x > level ? 1 << k : 0), 0);
    if (code === 0 || code === 15) continue;
    const e = [() => cross(i, j, i + 1, j), () => cross(i + 1, j, i + 1, j + 1), () => cross(i, j + 1, i + 1, j + 1), () => cross(i, j, i, j + 1)];
    const centreHigh = (c[0] + c[1] + c[2] + c[3]) / 4 > level;
    const table = {
      1: [[3, 0]], 2: [[0, 1]], 3: [[3, 1]], 4: [[1, 2]], 6: [[0, 2]], 7: [[3, 2]], 8: [[2, 3]],
      9: [[0, 2]], 11: [[1, 2]], 12: [[1, 3]], 13: [[0, 1]], 14: [[3, 0]],
      5: centreHigh ? [[3, 2], [1, 0]] : [[3, 0], [1, 2]], 10: centreHigh ? [[0, 3], [2, 1]] : [[0, 1], [2, 3]],
    };
    for (const [p, q] of table[code]) link(e[p](), e[q]());
  }}
  const used = new Set(), lines = [];
  const walk = start => {
    const ids = [start]; used.add(start);
    for (let cur = start; ;) {
      const nextId = (links.get(cur) ?? []).find(n => !used.has(n));
      if (nextId === undefined) break;
      ids.push(nextId); used.add(nextId); cur = nextId;
    }
    return ids;
  };
  for (const [id, nbrs] of links) if (nbrs.length === 1 && !used.has(id)) lines.push({ points: walk(id).map(k => point.get(k)), closed: false });
  for (const id of links.keys()) {
    if (used.has(id)) continue;
    const ids = walk(id);
    lines.push({ points: [...ids, ids[0]].map(k => point.get(k)), closed: true });
  }
  return lines;
}

/**
 * The density level whose super-level set holds `fraction` of the total (e.g. 0.9 for "the
 * surface enclosing 90 % of the probability"), for non-negative cell values on a uniform grid.
 * Cell-resolution estimate: refine the grid to check it.
 */
export function densityLevel(values, fraction) {
  if (!(fraction > 0 && fraction < 1)) throw new Error('densityLevel: fraction must be in (0, 1)');
  const sorted = Float64Array.from(values, x => { if (!(x >= 0)) throw new Error('densityLevel: values must be non-negative and finite'); return x; }).sort().reverse();
  const total = sorted.reduce((s, x) => s + x, 0);
  if (!(total > 0)) throw new Error('densityLevel: total must be positive');
  let acc = 0;
  for (const x of sorted) { acc += x; if (acc >= fraction * total) return x; }
  return 0;
}

/**
 * Two coupled surfaces (nonadiabatic dynamics: conical intersections, avoided crossings, curve
 * crossing): i ħ ∂ψ_k/∂t = T ψ_k + Σ_l V_kl ψ_l with a real symmetric 2×2 diabatic potential matrix
 * (V11, V22, V12 arrays on the grid). Strang splitting: exact exp(−iVΔt/2ħ) of the 2×2 matrix at
 * every point, shared spectral kinetic step (mass a number or per axis). Second order in Δt,
 * unitary to round-off. Returns { t, psi: [{ re, im }, { re, im }], step(), advanceTo(t),
 * populations() → [P1, P2] (diabatic), norm(), energy() }.
 */
export function coupledSchrodinger2({ psi0, V11, V22, V12, shape, spacing = 1, dt, mass = 1, hbar = 1, t0 = 0 }) {
  const spec = grid(shape, spacing); fftGrid(shape);
  if (!Array.isArray(psi0) || psi0.length !== 2) throw new Error('coupledSchrodinger2: psi0 must be [{ re, im }, { re, im }]');
  for (const a of [V11, V22, V12]) field(a, spec.size, 'coupledSchrodinger2 potential');
  positive(dt, 'coupledSchrodinger2 dt');
  const n = spec.size, vol = volume(spec), tau = dt / (2 * hbar);
  const psi = psi0.map(p => ({ re: Float64Array.from(p.re), im: p.im ? Float64Array.from(p.im) : new Float64Array(n) }));
  const Tk = kineticEnergies(spec, mass, hbar), kC = Float64Array.from(Tk, e => Math.cos(-e * dt / hbar)), kS = Float64Array.from(Tk, e => Math.sin(-e * dt / hbar));
  // exp(−iHτ) = e^{−i m τ} [cos(rτ) − i sin(rτ)(d σ_z + c σ_x)/r], H = [[m + d, c], [c, m − d]]
  const P = { pr: new Float64Array(n), pi: new Float64Array(n), a: new Float64Array(n), b: new Float64Array(n), c: new Float64Array(n) };
  for (let i = 0; i < n; i++) {
    const m = (V11[i] + V22[i]) / 2, d = (V11[i] - V22[i]) / 2, c = V12[i], r = Math.hypot(d, c);
    P.pr[i] = Math.cos(-m * tau); P.pi[i] = Math.sin(-m * tau);
    P.a[i] = Math.cos(r * tau); const sr = r > 0 ? Math.sin(r * tau) / r : tau;
    P.b[i] = sr * d; P.c[i] = sr * c;                            // −i·(b σ_z + c σ_x)
  }
  const potential = () => {
    const [u, w] = psi;
    for (let i = 0; i < n; i++) {
      // [x1; x2] = (a − i b σz − i c σx)[u; w]
      const ur = u.re[i], ui = u.im[i], wr = w.re[i], wi = w.im[i], a = P.a[i], b = P.b[i], c = P.c[i];
      const x1r = a * ur + b * ui + c * wi, x1i = a * ui - b * ur - c * wr;
      const x2r = a * wr - b * wi + c * ui, x2i = a * wi + b * wr - c * ur;
      const pr = P.pr[i], pi = P.pi[i];
      u.re[i] = pr * x1r - pi * x1i; u.im[i] = pr * x1i + pi * x1r;
      w.re[i] = pr * x2r - pi * x2i; w.im[i] = pr * x2i + pi * x2r;
    }
  };
  const kinetic = () => {
    for (const p of psi) {
      transform(p.re, p.im, spec, false);
      for (let i = 0; i < n; i++) { const a = p.re[i], b = p.im[i]; p.re[i] = a * kC[i] - b * kS[i]; p.im[i] = a * kS[i] + b * kC[i]; }
      transform(p.re, p.im, spec, true);
    }
  };
  const pop = p => { let s = 0; for (let i = 0; i < n; i++) s += p.re[i] * p.re[i] + p.im[i] * p.im[i]; return s * vol; };
  const s = {
    t: t0, psi, steps: 0,
    step() {
      potential(); kinetic(); potential();
      if (!Number.isFinite(psi[0].re[0] + psi[1].re[0])) throw new Error(`coupledSchrodinger2: non-finite state at t = ${s.t}`);
      s.t += dt; s.steps++;
    },
    advanceTo(target) { while (s.t < target - 1e-12 * dt) s.step(); return s; },
    populations: () => psi.map(pop),
    norm: () => pop(psi[0]) + pop(psi[1]),
    energy() {
      let e = 0;
      for (const p of psi) {
        const r = Float64Array.from(p.re), im = Float64Array.from(p.im); transform(r, im, spec, false);
        let t = 0; for (let i = 0; i < n; i++) t += (r[i] * r[i] + im[i] * im[i]) * Tk[i];
        e += t / n * vol;                                          // Parseval: Σ|ψ̃|² = n Σ|ψ|²
      }
      const [u, w] = psi;
      for (let i = 0; i < n; i++) {
        e += vol * (V11[i] * (u.re[i] ** 2 + u.im[i] ** 2) + V22[i] * (w.re[i] ** 2 + w.im[i] ** 2) + 2 * V12[i] * (u.re[i] * w.re[i] + u.im[i] * w.im[i]));
      }
      return e / s.norm();
    },
  };
  return s;
}

