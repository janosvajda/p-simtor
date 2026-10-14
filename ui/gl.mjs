// GPU rendering (WebGL2) for any subject, organised by what is drawn, not by topic:
//
//   createGLLayer(canvas2d)       a WebGL2 canvas behind the shell's 2D canvas (labels stay in 2D)
//   orbitCamera({...})            3D camera with drag-to-rotate and wheel zoom; project() for 2D labels
//   createPointRenderer(gl)       glowing point sprites: 10⁵–10⁶ bodies, stars, samples, molecules
//   createFieldRenderer(gl)       2D scalar field as a colour-mapped image (heat maps, |ψ|², density)
//   createVolumeRenderer(gl)      3D scalar field by ray marching (orbital clouds, densities, gas)
//   colormap(name, t)             scientific colour maps (perceptually uniform), also for 2D canvas
//
// The renderers only READ arrays the simulation computed; they never change physics. Every
// mapping (range, log scale, colour map, opacity, glow) is a representation choice: state it
// in the legend. Matrices are column-major Float32Array(16), as WebGL expects.

// ---------------------------------------------------------------- matrix maths (pure; tested in Node)

export function perspective(fovy, aspect, near, far) {
  const f = 1 / Math.tan(fovy / 2), nf = 1 / (near - far);
  return new Float32Array([f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) * nf, -1, 0, 0, 2 * far * near * nf, 0]);
}

export function lookAt(eye, target, up) {
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const unit = a => { const l = Math.hypot(...a) || 1; return a.map(v => v / l); };
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const z = unit(sub(eye, target)), x = unit(cross(up, z)), y = cross(z, x);
  return new Float32Array([x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0, -dot(x, eye), -dot(y, eye), -dot(z, eye), 1]);
}

export function multiply(a, b) {
  const out = new Float32Array(16);
  for (let c = 0; c < 4; c++) {for (let r = 0; r < 4; r++) {
    let s = 0;
    for (let k = 0; k < 4; k++) s += a[k * 4 + r] * b[c * 4 + k];
    out[c * 4 + r] = s;
  }}
  return out;
}

/** Inverse of a 4×4 matrix (cofactor expansion, in float64). */
export function invert(m) {
  const a = Array.from(m), inv = new Array(16);
  inv[0] = a[5] * a[10] * a[15] - a[5] * a[11] * a[14] - a[9] * a[6] * a[15] + a[9] * a[7] * a[14] + a[13] * a[6] * a[11] - a[13] * a[7] * a[10];
  inv[4] = -a[4] * a[10] * a[15] + a[4] * a[11] * a[14] + a[8] * a[6] * a[15] - a[8] * a[7] * a[14] - a[12] * a[6] * a[11] + a[12] * a[7] * a[10];
  inv[8] = a[4] * a[9] * a[15] - a[4] * a[11] * a[13] - a[8] * a[5] * a[15] + a[8] * a[7] * a[13] + a[12] * a[5] * a[11] - a[12] * a[7] * a[9];
  inv[12] = -a[4] * a[9] * a[14] + a[4] * a[10] * a[13] + a[8] * a[5] * a[14] - a[8] * a[6] * a[13] - a[12] * a[5] * a[10] + a[12] * a[6] * a[9];
  inv[1] = -a[1] * a[10] * a[15] + a[1] * a[11] * a[14] + a[9] * a[2] * a[15] - a[9] * a[3] * a[14] - a[13] * a[2] * a[11] + a[13] * a[3] * a[10];
  inv[5] = a[0] * a[10] * a[15] - a[0] * a[11] * a[14] - a[8] * a[2] * a[15] + a[8] * a[3] * a[14] + a[12] * a[2] * a[11] - a[12] * a[3] * a[10];
  inv[9] = -a[0] * a[9] * a[15] + a[0] * a[11] * a[13] + a[8] * a[1] * a[15] - a[8] * a[3] * a[13] - a[12] * a[1] * a[11] + a[12] * a[3] * a[9];
  inv[13] = a[0] * a[9] * a[14] - a[0] * a[10] * a[13] - a[8] * a[1] * a[14] + a[8] * a[2] * a[13] + a[12] * a[1] * a[10] - a[12] * a[2] * a[9];
  inv[2] = a[1] * a[6] * a[15] - a[1] * a[7] * a[14] - a[5] * a[2] * a[15] + a[5] * a[3] * a[14] + a[13] * a[2] * a[7] - a[13] * a[3] * a[6];
  inv[6] = -a[0] * a[6] * a[15] + a[0] * a[7] * a[14] + a[4] * a[2] * a[15] - a[4] * a[3] * a[14] - a[12] * a[2] * a[7] + a[12] * a[3] * a[6];
  inv[10] = a[0] * a[5] * a[15] - a[0] * a[7] * a[13] - a[4] * a[1] * a[15] + a[4] * a[3] * a[13] + a[12] * a[1] * a[7] - a[12] * a[3] * a[5];
  inv[14] = -a[0] * a[5] * a[14] + a[0] * a[6] * a[13] + a[4] * a[1] * a[14] - a[4] * a[2] * a[13] - a[12] * a[1] * a[6] + a[12] * a[2] * a[5];
  inv[3] = -a[1] * a[6] * a[11] + a[1] * a[7] * a[10] + a[5] * a[2] * a[11] - a[5] * a[3] * a[10] - a[9] * a[2] * a[7] + a[9] * a[3] * a[6];
  inv[7] = a[0] * a[6] * a[11] - a[0] * a[7] * a[10] - a[4] * a[2] * a[11] + a[4] * a[3] * a[10] + a[8] * a[2] * a[7] - a[8] * a[3] * a[6];
  inv[11] = -a[0] * a[5] * a[11] + a[0] * a[7] * a[9] + a[4] * a[1] * a[11] - a[4] * a[3] * a[9] - a[8] * a[1] * a[7] + a[8] * a[3] * a[5];
  inv[15] = a[0] * a[5] * a[10] - a[0] * a[6] * a[9] - a[4] * a[1] * a[10] + a[4] * a[2] * a[9] + a[8] * a[1] * a[6] - a[8] * a[2] * a[5];
  const det = a[0] * inv[0] + a[1] * inv[4] + a[2] * inv[8] + a[3] * inv[12];
  if (!det) throw new Error('invert: singular matrix');
  return new Float32Array(inv.map(v => v / det));
}

/**
 * Orbit camera around `target`, z up (physics convention). yaw/pitch in radians.
 * matrices(aspect) → { view, proj, viewProj, eye }; project([x,y,z], W, H) → CSS-pixel position
 * (and depth) for 2D labels drawn over the GL scene; attach(element) adds drag-rotate and
 * wheel-zoom and returns a detach function.
 */
export function orbitCamera({ target = [0, 0, 0], distance = 5, yaw = 0.6, pitch = 0.45, fov = 0.8, near, far } = {}) {
  const cam = {
    target: [...target], distance, yaw, pitch, fov,
    eye() {
      const cp = Math.cos(cam.pitch);
      return [cam.target[0] + cam.distance * cp * Math.cos(cam.yaw), cam.target[1] + cam.distance * cp * Math.sin(cam.yaw), cam.target[2] + cam.distance * Math.sin(cam.pitch)];
    },
    matrices(aspect) {
      const eye = cam.eye();
      const view = lookAt(eye, cam.target, [0, 0, 1]);
      const proj = perspective(cam.fov, aspect, near ?? cam.distance * 1e-3, far ?? cam.distance * 1e3);
      return { view, proj, viewProj: multiply(proj, view), eye };
    },
    project(p, W, H) {
      const m = cam.matrices(W / H).viewProj;
      const x = m[0] * p[0] + m[4] * p[1] + m[8] * p[2] + m[12], y = m[1] * p[0] + m[5] * p[1] + m[9] * p[2] + m[13];
      const w = m[3] * p[0] + m[7] * p[1] + m[11] * p[2] + m[15];
      return { x: (x / w * 0.5 + 0.5) * W, y: (1 - (y / w * 0.5 + 0.5)) * H, depth: w, visible: w > 0 };
    },
    attach(el) {
      let last = null;
      const down = e => { last = [e.clientX, e.clientY]; el.setPointerCapture?.(e.pointerId); };
      const move = e => {
        if (!last) return;
        cam.yaw -= (e.clientX - last[0]) * 0.006;
        cam.pitch = Math.max(-1.55, Math.min(1.55, cam.pitch + (e.clientY - last[1]) * 0.006));
        last = [e.clientX, e.clientY];
      };
      const up = () => { last = null; };
      const wheel = e => { e.preventDefault(); cam.distance *= Math.exp(e.deltaY * 0.001); };
      el.addEventListener('pointerdown', down); el.addEventListener('pointermove', move);
      el.addEventListener('pointerup', up); el.addEventListener('pointercancel', up);
      el.addEventListener('wheel', wheel, { passive: false });
      return () => { el.removeEventListener('pointerdown', down); el.removeEventListener('pointermove', move); el.removeEventListener('pointerup', up); el.removeEventListener('pointercancel', up); el.removeEventListener('wheel', wheel); };
    },
  };
  return cam;
}

// ---------------------------------------------------------------- colour maps

// Anchors sampled from matplotlib's perceptually uniform maps (van der Walt & Smith, CC0) and
// Moreland's diverging "cool-warm". Linear interpolation between anchors.
const MAPS = {
  viridis: ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'],
  inferno: ['#000004', '#1b0c41', '#4a0c6b', '#781c6d', '#a52c60', '#cf4446', '#ed6925', '#fb9b06', '#f7d13d', '#fcffa4'],
  magma: ['#000004', '#180f3d', '#440f76', '#721f81', '#9e2f7f', '#cd4071', '#f1605d', '#fd9668', '#feca8d', '#fcfdbf'],
  coolwarm: ['#3b4cc0', '#6282ea', '#8db0fe', '#b8d0f9', '#dddddd', '#f5c4ad', '#f49a7b', '#de604d', '#b40426'],
};
const hex = h => [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16) / 255);
export const COLORMAPS = Object.keys(MAPS);

/** Colour of map `name` at t ∈ [0, 1] (clamped): [r, g, b] in [0, 1]. */
export function colormap(name, t) {
  const a = MAPS[name];
  if (!a) throw new Error(`colormap: unknown map "${name}" (${COLORMAPS.join(', ')})`);
  const x = Math.min(1, Math.max(0, Number.isFinite(t) ? t : 0)) * (a.length - 1), i = Math.min(a.length - 2, Math.floor(x)), f = x - i;
  const p = hex(a[i]), q = hex(a[i + 1]);
  return [p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1]), p[2] + f * (q[2] - p[2])];
}

/** Map a value to [0, 1] for a colour map: linear or log10 between range[0] and range[1]. */
export function normalize(v, [lo, hi], log = false) {
  if (log) { if (!(v > 0 && lo > 0)) return 0; return (Math.log10(v) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo)); }
  return (v - lo) / (hi - lo);
}

// ---------------------------------------------------------------- GL plumbing

function compile(gl, vs, fs) {
  const shader = (type, src) => {
    const s = gl.createShader(type);
    gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(`GL shader: ${gl.getShaderInfoLog(s)}`);
    return s;
  };
  const p = gl.createProgram();
  gl.attachShader(p, shader(gl.VERTEX_SHADER, vs)); gl.attachShader(p, shader(gl.FRAGMENT_SHADER, fs));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(`GL program: ${gl.getProgramInfoLog(p)}`);
  const u = {};
  for (let i = 0; i < gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS); i++) { const n = gl.getActiveUniform(p, i).name; u[n] = gl.getUniformLocation(p, n); }
  return { p, u };
}

const cmapTextures = new WeakMap();
function cmapTexture(gl, name) {
  let cache = cmapTextures.get(gl);
  if (!cache) cmapTextures.set(gl, cache = {});
  if (cache[name]) return cache[name];
  const data = new Uint8Array(256 * 4);
  for (let i = 0; i < 256; i++) { const c = colormap(name, i / 255); data.set([c[0] * 255, c[1] * 255, c[2] * 255, 255].map(Math.round), 4 * i); }
  const tex = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, 256, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, data);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  return (cache[name] = tex);
}
const rgb = c => (typeof c === 'string' ? hex(c.length === 4 ? `#${c[1]}${c[1]}${c[2]}${c[2]}${c[3]}${c[3]}` : c) : c);

/**
 * A WebGL2 canvas stacked under `canvas2d` (same box, same device-pixel ratio). The 2D canvas
 * stays on top for labels, axes and pointer events. Returns { gl, canvas, ok, reason, sync() };
 * when WebGL2 is unavailable, ok is false and reason says why (show it; draw a 2D fallback).
 */
export function createGLLayer(canvas2d) {
  const canvas = document.createElement('canvas');
  canvas.className = 'psim-gl';
  canvas.setAttribute('aria-hidden', 'true');
  canvas2d.parentNode.insertBefore(canvas, canvas2d);
  const gl = canvas.getContext('webgl2', { alpha: true, premultipliedAlpha: true, antialias: true, preserveDrawingBuffer: true });
  const layer = {
    gl, canvas, ok: !!gl, reason: gl ? '' : 'WebGL2 is not available in this browser',
    /** Match the 2D canvas size; set the viewport and clear to transparent. */
    sync(dpr = devicePixelRatio || 1) {
      const w = Math.round(canvas2d.clientWidth * dpr), h = Math.round(canvas2d.clientHeight * dpr);
      if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
      if (!gl) return;
      gl.viewport(0, 0, w, h); gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    },
  };
  return layer;
}

// ---------------------------------------------------------------- points

/**
 * Point sprites. update(positions [x,y,z,…], { values, colors }) uploads up to any count;
 * draw(camera, { aspect, size, sizeUnit: 'px' | 'world', color, colormap, range, log, blend,
 * opacity }) with blend 'additive' (glow on dark backgrounds) or 'alpha' (light backgrounds).
 * values + colormap colour each point by a computed quantity; otherwise colors (r,g,b per point)
 * or a single color. Sizes in 'world' units shrink with distance (true perspective).
 */
export function createPointRenderer(gl) {
  const { p, u } = compile(gl, `#version 300 es
in vec3 aPos; in float aVal; in vec3 aCol;
uniform mat4 uVP; uniform float uSize, uWorld, uScale, uLog; uniform vec2 uRange;
out float vVal; out vec3 vCol;
void main() {
  gl_Position = uVP * vec4(aPos, 1.0);
  float px = uWorld > 0.5 ? uSize * uScale / max(gl_Position.w, 1e-6) : uSize;
  gl_PointSize = clamp(px, 1.0, 256.0);
  float v = uLog > 0.5 ? (log(max(aVal, 1e-30)) / log(10.0) - uRange.x) / (uRange.y - uRange.x) : (aVal - uRange.x) / (uRange.y - uRange.x);
  vVal = clamp(v, 0.0, 1.0); vCol = aCol;
}`, `#version 300 es
precision highp float;
in float vVal; in vec3 vCol;
uniform sampler2D uMap; uniform float uMode, uGlow, uOpacity; uniform vec3 uColor;
out vec4 frag;
void main() {
  vec2 d = gl_PointCoord * 2.0 - 1.0; float r2 = dot(d, d);
  if (r2 > 1.0) discard;
  vec3 c = uMode < 0.5 ? uColor : uMode < 1.5 ? texture(uMap, vec2(vVal, 0.5)).rgb : vCol;
  float a = uGlow > 0.5 ? exp(-r2 * 5.0) + 0.6 * smoothstep(0.25, 0.0, r2) : smoothstep(1.0, 0.75, r2);
  a *= uOpacity;
  frag = vec4(c * a, a);
}`);
  const vao = gl.createVertexArray(), bufs = {};
  gl.bindVertexArray(vao);
  for (const [name, size] of [['aPos', 3], ['aVal', 1], ['aCol', 3]]) {
    const loc = gl.getAttribLocation(p, name);
    bufs[name] = { buf: gl.createBuffer(), loc, size };
    if (loc < 0) continue;
    gl.bindBuffer(gl.ARRAY_BUFFER, bufs[name].buf);
    gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 0, 0);
  }
  gl.bindVertexArray(null);
  let count = 0, has = { values: false, colors: false };
  const upload = (name, data) => {
    gl.bindBuffer(gl.ARRAY_BUFFER, bufs[name].buf);
    gl.bufferData(gl.ARRAY_BUFFER, data instanceof Float32Array ? data : Float32Array.from(data), gl.DYNAMIC_DRAW);
  };
  return {
    get count() { return count; },
    update(positions, { values, colors } = {}) {
      count = Math.floor(positions.length / 3);
      upload('aPos', positions);
      has = { values: !!values, colors: !!colors };
      upload('aVal', values ?? new Float32Array(count));
      upload('aCol', colors ?? new Float32Array(3 * count));
    },
    draw(camera, { aspect, size = 3, sizeUnit = 'px', color = '#ffd27a', colormap: map = 'inferno', range = [0, 1], log = false, blend = 'additive', opacity = 1 } = {}) {
      if (!count) return;
      const { viewProj, proj } = camera.matrices(aspect ?? gl.drawingBufferWidth / gl.drawingBufferHeight);
      gl.useProgram(p); gl.bindVertexArray(vao);
      gl.uniformMatrix4fv(u.uVP, false, viewProj);
      gl.uniform1f(u.uSize, size * (sizeUnit === 'px' ? (devicePixelRatio || 1) : 1));
      gl.uniform1f(u.uWorld, sizeUnit === 'world' ? 1 : 0);
      gl.uniform1f(u.uScale, proj[5] * gl.drawingBufferHeight / 2);
      gl.uniform1f(u.uLog, log ? 1 : 0);
      gl.uniform2f(u.uRange, ...(log ? range.map(Math.log10) : range));
      gl.uniform1f(u.uMode, has.values ? 1 : has.colors ? 2 : 0);
      gl.uniform3f(u.uColor, ...rgb(color));
      gl.uniform1f(u.uGlow, blend === 'additive' ? 1 : 0);
      gl.uniform1f(u.uOpacity, opacity);
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, cmapTexture(gl, map)); gl.uniform1i(u.uMap, 0);
      gl.enable(gl.BLEND); gl.disable(gl.DEPTH_TEST);
      if (blend === 'additive') gl.blendFunc(gl.ONE, gl.ONE); else gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
      gl.drawArrays(gl.POINTS, 0, count);
      gl.bindVertexArray(null);
    },
  };
}

// ---------------------------------------------------------------- 2D fields

/**
 * A 2D scalar field (index i + nx·j, as in core/fields.mjs) drawn as a colour-mapped image.
 * update(values, nx, ny, { range, log }) normalises on the CPU (NaN or masked cells are
 * transparent); draw({ rect: [x, y, w, h] in CSS px (default: whole canvas), colormap, smooth,
 * opacity }) with smooth false showing the grid cells honestly.
 */
export function createFieldRenderer(gl) {
  const { p, u } = compile(gl, `#version 300 es
in vec2 aXY; uniform vec4 uRect; out vec2 vUV;
void main() { vUV = aXY; vec2 xy = uRect.xy + aXY * uRect.zw; gl_Position = vec4(xy, 0.0, 1.0); }`, `#version 300 es
precision highp float;
in vec2 vUV; uniform sampler2D uField, uMap; uniform float uOpacity; out vec4 frag;
void main() { float v = texture(uField, vUV).r; if (v < 0.0) discard; vec3 c = texture(uMap, vec2(v, 0.5)).rgb; frag = vec4(c * uOpacity, uOpacity); }`);
  const vao = gl.createVertexArray(), buf = gl.createBuffer(), tex = gl.createTexture();
  gl.bindVertexArray(vao); gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]), gl.STATIC_DRAW);
  const loc = gl.getAttribLocation(p, 'aXY'); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
  gl.bindVertexArray(null);
  let ready = false;
  return {
    update(values, nx, ny, { range, log = false } = {}) {
      if (values.length !== nx * ny) throw new Error('fieldRenderer: values must have nx·ny entries');
      let lo = range?.[0], hi = range?.[1];
      if (lo == null || hi == null) {
        lo = Infinity; hi = -Infinity;
        for (const v of values) if (Number.isFinite(v) && (!log || v > 0)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
        if (!(hi > lo)) hi = lo + 1;
      }
      const data = Float32Array.from(values, v => (Number.isFinite(v) ? Math.min(1, Math.max(0, normalize(v, [lo, hi], log))) : -1));
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.R16F, nx, ny, 0, gl.RED, gl.FLOAT, data);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      ready = true;
      return { range: [lo, hi] };
    },
    draw({ rect, colormap: map = 'viridis', smooth = true, opacity = 1 } = {}) {
      if (!ready) return;
      const W = gl.drawingBufferWidth / (devicePixelRatio || 1), H = gl.drawingBufferHeight / (devicePixelRatio || 1);
      const [x, y, w, h] = rect ?? [0, 0, W, H];
      gl.useProgram(p); gl.bindVertexArray(vao);
      // CSS px (y down) → clip space; the texture's row j = 0 is drawn at the bottom (y up, physics)
      gl.uniform4f(u.uRect, 2 * x / W - 1, 1 - 2 * (y + h) / H, 2 * w / W, 2 * h / H);
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, tex);
      const f = smooth ? gl.LINEAR : gl.NEAREST;
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, f); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, f);
      gl.uniform1i(u.uField, 0);
      gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, cmapTexture(gl, map)); gl.uniform1i(u.uMap, 1);
      gl.uniform1f(u.uOpacity, opacity);
      gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.disable(gl.DEPTH_TEST);
      gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
      gl.bindVertexArray(null);
    },
  };
}

// ---------------------------------------------------------------- 3D volumes

/**
 * A 3D scalar field (index i + nx·(j + ny·k)) rendered by ray marching through its box.
 * update(values, nx, ny, nz, { range, log }); draw(camera, { box: [[x0,y0,z0],[x1,y1,z1]],
 * colormap, opacity, steps, mode }) with mode 'emission' (glowing gas: colour ∝ density, adds up)
 * or 'absorption' (emission–absorption: nearer matter hides farther). Opacity is a display
 * transfer function, not a physical optical depth unless you derive it from one: say which.
 */
export function createVolumeRenderer(gl) {
  const { p, u } = compile(gl, `#version 300 es
in vec2 aXY; out vec2 vNDC;
void main() { vNDC = aXY * 2.0 - 1.0; gl_Position = vec4(vNDC, 0.0, 1.0); }`, `#version 300 es
precision highp float; precision highp sampler3D;
in vec2 vNDC;
uniform mat4 uInvVP; uniform vec3 uMin, uMax; uniform sampler3D uVol; uniform sampler2D uMap;
uniform float uOpacity, uSteps, uMode;
out vec4 frag;
void main() {
  vec4 a = uInvVP * vec4(vNDC, -1.0, 1.0), b = uInvVP * vec4(vNDC, 1.0, 1.0);
  vec3 ro = a.xyz / a.w, rd = normalize(b.xyz / b.w - ro);
  vec3 t0 = (uMin - ro) / rd, t1 = (uMax - ro) / rd;
  vec3 tn = min(t0, t1), tf = max(t0, t1);
  float near = max(max(tn.x, tn.y), max(tn.z, 0.0)), far = min(min(tf.x, tf.y), tf.z);
  if (far <= near) discard;
  float dt = (far - near) / uSteps, len = length(uMax - uMin);
  vec3 col = vec3(0.0); float alpha = 0.0;
  for (float s = 0.5; s < 1024.0; s += 1.0) {
    if (s > uSteps) break;
    vec3 pos = ro + rd * (near + s * dt);
    float v = texture(uVol, (pos - uMin) / (uMax - uMin)).r;
    if (v <= 0.0) continue;
    vec3 c = texture(uMap, vec2(v, 0.5)).rgb;
    float k = uOpacity * v * dt / len * 8.0;
    if (uMode < 0.5) { col += c * k; alpha += k; }
    else { float ak = 1.0 - exp(-k); col += (1.0 - alpha) * c * ak; alpha += (1.0 - alpha) * ak; if (alpha > 0.99) break; }
  }
  alpha = clamp(alpha, 0.0, 1.0);
  frag = vec4(min(col, vec3(1.0)), alpha);
}`);
  const vao = gl.createVertexArray(), buf = gl.createBuffer(), tex = gl.createTexture();
  gl.bindVertexArray(vao); gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]), gl.STATIC_DRAW);
  const loc = gl.getAttribLocation(p, 'aXY'); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
  gl.bindVertexArray(null);
  let ready = false;
  return {
    update(values, nx, ny, nz, { range, log = false } = {}) {
      if (values.length !== nx * ny * nz) throw new Error('volumeRenderer: values must have nx·ny·nz entries');
      let lo = range?.[0], hi = range?.[1];
      if (lo == null || hi == null) {
        lo = Infinity; hi = -Infinity;
        for (const v of values) if (Number.isFinite(v) && (!log || v > 0)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
        if (!(hi > lo)) hi = lo + 1;
      }
      const data = Float32Array.from(values, v => (Number.isFinite(v) ? Math.min(1, Math.max(0, normalize(v, [lo, hi], log))) : 0));
      gl.bindTexture(gl.TEXTURE_3D, tex);
      gl.texImage3D(gl.TEXTURE_3D, 0, gl.R16F, nx, ny, nz, 0, gl.RED, gl.FLOAT, data);
      for (const w of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T, gl.TEXTURE_WRAP_R]) gl.texParameteri(gl.TEXTURE_3D, w, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR); gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      ready = true;
      return { range: [lo, hi] };
    },
    draw(camera, { aspect, box = [[-1, -1, -1], [1, 1, 1]], colormap: map = 'inferno', opacity = 1, steps = 160, mode = 'emission' } = {}) {
      if (!ready) return;
      const { viewProj } = camera.matrices(aspect ?? gl.drawingBufferWidth / gl.drawingBufferHeight);
      gl.useProgram(p); gl.bindVertexArray(vao);
      gl.uniformMatrix4fv(u.uInvVP, false, invert(viewProj));
      gl.uniform3f(u.uMin, ...box[0]); gl.uniform3f(u.uMax, ...box[1]);
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_3D, tex); gl.uniform1i(u.uVol, 0);
      gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, cmapTexture(gl, map)); gl.uniform1i(u.uMap, 1);
      gl.uniform1f(u.uOpacity, opacity); gl.uniform1f(u.uSteps, Math.min(1000, steps)); gl.uniform1f(u.uMode, mode === 'emission' ? 0 : 1);
      gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.disable(gl.DEPTH_TEST);
      gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
      gl.bindVertexArray(null);
    },
  };
}
