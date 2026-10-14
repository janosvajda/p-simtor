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
function transform(re, im, spec, inverse) {
  for (let axis = 0; axis < spec.shape.length; axis++) {
    const length = spec.shape[axis], stride = spec.strides[axis], block = length * stride;
    const rowRe = new Float64Array(length), rowIm = new Float64Array(length);
    for (let base = 0; base < spec.size; base += block) { for (let offset = 0; offset < stride; offset++) {
      for (let i = 0; i < length; i++) { const index = base + offset + i * stride; rowRe[i] = re[index]; rowIm[i] = im[index]; }
      fft(rowRe, rowIm, inverse);
      for (let i = 0; i < length; i++) { const index = base + offset + i * stride; re[index] = finite(rowRe[i], 'fftND'); im[index] = finite(rowIm[i], 'fftND'); }
    } }
  }
}

/** In-place separable radix-2 complex FFT; inverse includes 1/product(shape). */
export function fftND(re, im, shape, inverse = false) {
  const spec = fftGrid(shape);
  field(re, spec.size, 'fftND'); field(im, spec.size, 'fftND');
  if (overlaps(re, im)) throw new Error('fftND: real and imaginary storage must not overlap');
  transform(re, im, spec, inverse);
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
  positive(dt, 'schrodingerND dt'); positive(mass, 'schrodingerND mass'); positive(hbar, 'schrodingerND hbar'); finite(t0, 'schrodingerND t0');
  const lambda = eigenvalues(spec, 'spectral');
  const halfV = Float64Array.from({ length: spec.size }, (_, i) => finite(-(V?.[i] ?? 0) / hbar * (dt / 2), 'schrodingerND potential phase'));
  const kinetic = Float64Array.from(lambda, value => finite(-hbar / mass * (value / 2) * dt, 'schrodingerND kinetic phase'));
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
  const damp = () => { if (halfDamp) for (let i = 0; i < spec.size; i++) { workRe[i] *= halfDamp[i]; workIm[i] *= halfDamp[i]; } };
  const rotate = phases => {
    for (let i = 0; i < spec.size; i++) {
      const c = Math.cos(phases[i]), d = Math.sin(phases[i]), a = workRe[i], b = workIm[i];
      workRe[i] = a * c - b * d; workIm[i] = a * d + b * c;
    }
  };
  const s = {
    t: t0, re, im, steps: 0, shape: spec.shape.slice(),
    step() {
      const nextTime = finite(s.t + dt, 'schrodingerND');
      if (nextTime === s.t) throw new Error('schrodingerND: step underflow');
      workRe.set(re); workIm.set(im);
      rotate(halfV); damp(); transform(workRe, workIm, spec, false); rotate(kinetic); transform(workRe, workIm, spec, true); rotate(halfV); damp();
      field(workRe, spec.size, 'schrodingerND'); field(workIm, spec.size, 'schrodingerND');
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
  };
  const initialNorm = s.norm();
  return s;
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
