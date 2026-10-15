// GPU rendering (WebGL2) for any subject, organised by what is drawn, not by topic:
//
//   createGLLayer(canvas2d)       a WebGL2 canvas behind the shell's 2D canvas (labels stay in 2D)
//   orbitCamera({...})            3D camera with drag-to-rotate and wheel zoom; project() for 2D labels
//   createPointRenderer(gl)       glowing point sprites: 10⁵–10⁶ bodies, stars, samples, molecules
//   createFieldRenderer(gl)       2D scalar field as a colour-mapped image (heat maps, |ψ|², density)
//   createVolumeRenderer(gl)      3D scalar field by ray marching (orbital clouds, densities, gas)
//   createSphereRenderer(gl)      lit spheres (planets, moons, balls, atoms) with exact silhouettes
//   createLineRenderer(gl)        thick 3D polylines (orbits, field lines, trails, axes), solid or dashed
//   createMeshRenderer(gl)        lit solid objects (geometry.*: box, cylinder, cone, sphere, disc, plane, arrow)
//   colormap(name, t)             scientific colour maps (perceptually uniform), also for 2D canvas
//   drawColorBar(ctx, {...})      the matching colour-bar legend on the 2D canvas (every mapped value needs one)
//   polylines.*                   line shapes for createLineRenderer: circle, boxEdges, grid, arrow
//   alignZ(dir)                   rotation taking +z to a direction (orient arrows, cones, cylinders)
//
// Large coordinates (planetary distances): spheres and lines are drawn relative to the camera
// target, subtracted in float64 before upload, so precision does not degrade far from the origin.
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
export function orbitCamera({ target = [0, 0, 0], distance = 5, yaw = 0.6, pitch = 0.45, fov = 0.8, near, far, minDistance = 0, maxDistance = Infinity, minPitch = -1.55, maxPitch = 1.55 } = {}) {
  const clampD = d => Math.min(maxDistance, Math.max(minDistance, d));
  const cam = {
    target: [...target], distance, yaw, pitch, fov,
    /**
     * Frame a bounding box [[x0,y0,z0],[x1,y1,z1]] so it fits the view at this aspect ratio
     * (both the width and the height of the screen), keeping the current yaw and pitch.
     */
    fit(box, aspect = 1, margin = 1.1) {
      const c = [0, 1, 2].map(k => (box[0][k] + box[1][k]) / 2);
      const r = Math.hypot(...[0, 1, 2].map(k => box[1][k] - box[0][k])) / 2 * margin;
      const half = Math.min(cam.fov / 2, Math.atan(Math.tan(cam.fov / 2) * aspect));
      cam.target = c; cam.distance = clampD(r / Math.sin(half));
      return cam;
    },
    eye() {
      const cp = Math.cos(cam.pitch);
      return [cam.target[0] + cam.distance * cp * Math.cos(cam.yaw), cam.target[1] + cam.distance * cp * Math.sin(cam.yaw), cam.target[2] + cam.distance * Math.sin(cam.pitch)];
    },
    /** relative: true builds the view with the target at the origin (for target-relative data). */
    matrices(aspect, { relative = false } = {}) {
      const eye = cam.eye(), off = relative ? cam.target : [0, 0, 0];
      const view = lookAt([eye[0] - off[0], eye[1] - off[1], eye[2] - off[2]], [cam.target[0] - off[0], cam.target[1] - off[1], cam.target[2] - off[2]], [0, 0, 1]);
      const proj = perspective(cam.fov, aspect, near ?? cam.distance * 1e-3, far ?? cam.distance * 1e3);
      return { view, proj, viewProj: multiply(proj, view), eye };
    },
    project(p, W, H) {
      const m = cam.matrices(W / H).viewProj;
      const x = m[0] * p[0] + m[4] * p[1] + m[8] * p[2] + m[12], y = m[1] * p[0] + m[5] * p[1] + m[9] * p[2] + m[13];
      const w = m[3] * p[0] + m[7] * p[1] + m[11] * p[2] + m[15];
      return { x: (x / w * 0.5 + 0.5) * W, y: (1 - (y / w * 0.5 + 0.5)) * H, depth: w, visible: w > 0 };
    },
    /**
     * Drag to rotate, wheel or pinch to zoom. wheel: 'auto' (default) zooms with the plain wheel
     * only when the page around the canvas does not scroll; on a scrolling page (workbench layout,
     * phones) plain wheel scrolls the page and Ctrl/⌘ + wheel or a trackpad pinch zooms.
     * 'always' or 'modifier' force one behaviour.
     */
    attach(el, { wheel: wheelMode = 'auto' } = {}) {
      let last = null;
      const down = e => { last = [e.clientX, e.clientY]; el.setPointerCapture?.(e.pointerId); };
      const move = e => {
        if (!last) return;
        cam.yaw -= (e.clientX - last[0]) * 0.006;
        cam.pitch = Math.max(minPitch, Math.min(maxPitch, cam.pitch + (e.clientY - last[1]) * 0.006));
        last = [e.clientX, e.clientY];
      };
      const up = () => { last = null; };
      const scrolls = () => { for (let n = el.parentElement; n; n = n.parentElement) { const o = getComputedStyle(n).overflowY; if ((o === 'auto' || o === 'scroll') && n.scrollHeight > n.clientHeight + 1) return true; } return document.scrollingElement.scrollHeight > innerHeight + 1; };
      const wheel = e => {
        const zoom = wheelMode === 'always' || e.ctrlKey || e.metaKey || (wheelMode === 'auto' && !scrolls());
        if (!zoom) return;                                    // let the page scroll
        e.preventDefault(); cam.distance = clampD(cam.distance * Math.exp(e.deltaY * 0.001));
      };
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
export const COLORMAPS = [...Object.keys(MAPS), 'phase'];

/** OKLab → sRGB (Ottosson 2020), clipped to the displayable range. */
function oklabToSrgb(L, a, b) {
  const l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3, m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3, s = (L - 0.0894841775 * a - 1.2914855480 * b) ** 3;
  const lin = [4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s, -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s, -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s];
  return lin.map(x => { const c = Math.min(1, Math.max(0, x)); return c <= 0.0031308 ? 12.92 * c : 1.055 * c ** (1 / 2.4) - 0.055; });
}

/**
 * Colour of map `name` at t ∈ [0, 1]: [r, g, b] in [0, 1]. Sequential and diverging maps clamp t;
 * 'phase' is cyclic (t and t + 1 give the same colour): a constant-lightness hue circle in OKLab,
 * for angles such as the phase of a wavefunction, where no direction may look "brighter".
 */
export function colormap(name, t) {
  if (name === 'phase') { const h = 2 * Math.PI * ((Number.isFinite(t) ? t : 0) % 1); return oklabToSrgb(0.72, 0.11 * Math.cos(h), 0.11 * Math.sin(h)); }
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

export function compile(gl, vs, fs) {
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
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, name === 'phase' ? gl.REPEAT : gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  return (cache[name] = tex);
}
/**
 * Any colour → [r, g, b, a] in 0–1: '#rgb', '#rrggbb', 'rgb()/rgba()', CSS names, a theme token
 * '--name' (read from the page, so view.color('--x') and '--x' both work) or an [r, g, b(, a)] array.
 */
export function cssColor(c) {
  if (Array.isArray(c) || ArrayBuffer.isView(c)) return [c[0], c[1], c[2], c[3] ?? 1];
  let str = String(c).trim();
  if (str.startsWith('--') && typeof getComputedStyle === 'function') str = getComputedStyle(document.documentElement).getPropertyValue(str).trim();
  if (/^#[0-9a-f]{3}$/i.test(str)) str = `#${str[1]}${str[1]}${str[2]}${str[2]}${str[3]}${str[3]}`;
  if (/^#[0-9a-f]{6}$/i.test(str)) return [...hex(str), 1];
  if (/^#[0-9a-f]{8}$/i.test(str)) return [...hex(str.slice(0, 7)), parseInt(str.slice(7), 16) / 255];
  const m = str.match(/^rgba?\(\s*([\d.]+)[\s,]+([\d.]+)[\s,]+([\d.]+)(?:[\s,/]+([\d.]+%?))?\s*\)$/i);
  if (m) return [m[1] / 255, m[2] / 255, m[3] / 255, m[4] == null ? 1 : m[4].endsWith('%') ? parseFloat(m[4]) / 100 : +m[4]];
  if (typeof document !== 'undefined') {                      // CSS names and other syntaxes via the browser
    const ctx = document.createElement('canvas').getContext('2d');
    ctx.fillStyle = '#010203'; ctx.fillStyle = str;
    if (ctx.fillStyle !== '#010203' || /^#010203$/i.test(str)) return cssColor(ctx.fillStyle);
  }
  throw new Error(`gl: cannot read the colour "${c}"`);
}
const rgb = c => cssColor(c).slice(0, 3);

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
  const gl = canvas.getContext('webgl2', { alpha: true, premultipliedAlpha: true, antialias: true, preserveDrawingBuffer: true, stencil: true });
  const layer = {
    gl, canvas, ok: !!gl, reason: gl ? '' : 'WebGL2 is not available in this browser',
    /** Match the 2D canvas size; set the viewport and clear to transparent. */
    sync(dpr = devicePixelRatio || 1) {
      const w = Math.round(canvas2d.clientWidth * dpr), h = Math.round(canvas2d.clientHeight * dpr);
      if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
      if (!gl) return;
      gl.viewport(0, 0, w, h); gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT | gl.STENCIL_BUFFER_BIT);
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
 * soft: true gives 'alpha' points the same Gaussian profile as 'additive' ones, so a point looks
 * the same size in the light and dark themes. count draws only the first count points (a cloud
 * that grows from precomputed events: upload once with update(), then raise count each frame);
 * append(positions, { values, colors }) adds new points as they happen. colors are Float32 RGB
 * in 0–1, three per point.
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
uniform sampler2D uMap; uniform float uMode, uGlow, uOpacity, uSoft; uniform vec3 uColor;
out vec4 frag;
void main() {
  vec2 d = gl_PointCoord * 2.0 - 1.0; float r2 = dot(d, d);
  if (r2 > 1.0) discard;
  vec3 c = uMode < 0.5 ? uColor : uMode < 1.5 ? texture(uMap, vec2(vVal, 0.5)).rgb : vCol;
  float a = uGlow > 0.5 ? exp(-r2 * 5.0) + 0.6 * smoothstep(0.25, 0.0, r2)
          : uSoft > 0.5 ? exp(-r2 * 3.0) * smoothstep(1.0, 0.8, r2) : smoothstep(1.0, 0.75, r2);
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
  let count = 0, has = { values: false, colors: false }, cap = 0;
  const cpu = { aPos: new Float32Array(0), aVal: new Float32Array(0), aCol: new Float32Array(0) };
  const upload = (name, data) => {
    gl.bindBuffer(gl.ARRAY_BUFFER, bufs[name].buf);
    gl.bufferData(gl.ARRAY_BUFFER, data instanceof Float32Array ? data : Float32Array.from(data), gl.DYNAMIC_DRAW);
  };
  return {
    get count() { return count; },
    update(positions, { values, colors } = {}) {
      count = Math.floor(positions.length / 3); cap = 0;
      upload('aPos', positions);
      has = { values: !!values, colors: !!colors };
      upload('aVal', values ?? new Float32Array(count));
      upload('aCol', colors ?? new Float32Array(3 * count));
    },
    /**
     * Add points to the end (events arriving one by one: detector hits, decays, samples) without
     * re-uploading the old ones; capacity doubles when needed. values/colors as in update().
     */
    append(positions, { values, colors } = {}) {
      const k = Math.floor(positions.length / 3);
      if (!k) return;
      if (count === 0) has = { values: !!values, colors: !!colors };
      if (count + k > cap) {
        if (cap === 0 && count > 0) throw new Error('pointRenderer: append() after update() needs the points again: call append() only, or update() with all points');
        cap = Math.max(1024, 2 * (count + k));
        for (const [name, size] of [['aPos', 3], ['aVal', 1], ['aCol', 3]]) {
          const grown = new Float32Array(cap * size); grown.set(cpu[name].subarray(0, count * size)); cpu[name] = grown;
          upload(name, grown);
        }
      }
      const put = (name, size, data) => {
        const arr = data ? Float32Array.from(data) : new Float32Array(k * size);
        cpu[name].set(arr, count * size);
        gl.bindBuffer(gl.ARRAY_BUFFER, bufs[name].buf); gl.bufferSubData(gl.ARRAY_BUFFER, count * size * 4, arr);
      };
      put('aPos', 3, positions); put('aVal', 1, values); put('aCol', 3, colors);
      count += k;
    },
    /** Remove all points (keeps the allocated capacity for append). */
    clear() { count = 0; },
    draw(camera, { aspect, size = 3, sizeUnit = 'px', color = '#ffd27a', colormap: map = 'inferno', range = [0, 1], log = false, blend = 'additive', opacity = 1, soft = false, count: shown = count } = {}) {
      const n = Math.max(0, Math.min(count, Math.floor(shown)));
      if (!n) return;
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
      gl.uniform1f(u.uSoft, soft ? 1 : 0);
      gl.uniform1f(u.uOpacity, opacity);
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, cmapTexture(gl, map)); gl.uniform1i(u.uMap, 0);
      gl.enable(gl.BLEND); gl.disable(gl.DEPTH_TEST);
      if (blend === 'additive') gl.blendFunc(gl.ONE, gl.ONE); else gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
      gl.drawArrays(gl.POINTS, 0, n);
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
 * smooth: false shows each cell as a uniform block (honest for coarse grids, e.g. a lattice).
 * update(…, { signed: true }) shows a field of either sign (charge density, vorticity, Re ψ): the
 * colour map's middle is zero, colour gives the sign and opacity |value|; range is symmetric,
 * [−M, M] with M = max|value| unless given. Use a diverging map ('coolwarm').
 * Values ≤ 0 are transparent otherwise (opacity ∝ normalised value).
 */
export function createVolumeRenderer(gl) {
  const { p, u } = compile(gl, `#version 300 es
in vec2 aXY; out vec2 vNDC;
void main() { vNDC = aXY * 2.0 - 1.0; gl_Position = vec4(vNDC, 0.0, 1.0); }`, `#version 300 es
precision highp float; precision highp sampler3D;
in vec2 vNDC;
uniform mat4 uInvVP; uniform vec3 uMin, uMax; uniform sampler3D uVol; uniform sampler2D uMap;
uniform float uOpacity, uSteps, uMode, uSigned;
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
    float raw = texture(uVol, (pos - uMin) / (uMax - uMin)).r;
    float v = uSigned > 0.5 ? abs(raw) : raw;
    if (v <= 0.0) continue;
    vec3 c = texture(uMap, vec2(uSigned > 0.5 ? 0.5 + 0.5 * raw : v, 0.5)).rgb;
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
  let ready = false, signed = false;
  return {
    update(values, nx, ny, nz, { range, log = false, signed: isSigned = false } = {}) {
      if (values.length !== nx * ny * nz) throw new Error('volumeRenderer: values must have nx·ny·nz entries');
      signed = isSigned;
      if (signed) {
        let M = range ? Math.max(Math.abs(range[0]), Math.abs(range[1])) : 0;
        if (!range) for (const v of values) if (Number.isFinite(v)) M = Math.max(M, Math.abs(v));
        if (!(M > 0)) M = 1;
        const data = Float32Array.from(values, v => (Number.isFinite(v) ? Math.max(-1, Math.min(1, v / M)) : 0));
        gl.bindTexture(gl.TEXTURE_3D, tex);
        gl.texImage3D(gl.TEXTURE_3D, 0, gl.R16F, nx, ny, nz, 0, gl.RED, gl.FLOAT, data);
        for (const w of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T, gl.TEXTURE_WRAP_R]) gl.texParameteri(gl.TEXTURE_3D, w, gl.CLAMP_TO_EDGE);
        ready = true;
        return { range: [-M, M] };
      }
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
    draw(camera, { aspect, box = [[-1, -1, -1], [1, 1, 1]], colormap: map = 'inferno', opacity = 1, steps = 160, mode = 'emission', smooth = true } = {}) {
      if (!ready) return;
      const f = smooth ? gl.LINEAR : gl.NEAREST;
      gl.bindTexture(gl.TEXTURE_3D, tex); gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, f); gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, f);
      const { viewProj } = camera.matrices(aspect ?? gl.drawingBufferWidth / gl.drawingBufferHeight);
      gl.useProgram(p); gl.bindVertexArray(vao);
      gl.uniformMatrix4fv(u.uInvVP, false, invert(viewProj));
      gl.uniform3f(u.uMin, ...box[0]); gl.uniform3f(u.uMax, ...box[1]);
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_3D, tex); gl.uniform1i(u.uVol, 0);
      gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, cmapTexture(gl, map)); gl.uniform1i(u.uMap, 1);
      gl.uniform1f(u.uOpacity, opacity); gl.uniform1f(u.uSteps, Math.min(1000, steps)); gl.uniform1f(u.uMode, mode === 'emission' ? 0 : 1);
      gl.uniform1f(u.uSigned, signed ? 1 : 0);
      gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.disable(gl.DEPTH_TEST);
      gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
      gl.bindVertexArray(null);
    },
  };
}

// ---------------------------------------------------------------- spheres

/**
 * Lit spheres drawn as exact ray-cast impostors: true silhouettes, per-pixel shading and correct
 * overlap (depth) between spheres. draw(camera, spheres, { aspect, light, ambient }) where each
 * sphere is { center: [x, y, z], radius, color, emissive? } and light is a point [x, y, z]
 * (e.g. the Sun's centre: each sphere is lit from its own direction, giving the right
 * day/night line) or { direction: [x, y, z] } for parallel light. emissive spheres glow
 * unshaded (the Sun). Brightness ∝ max(0, n·l): the geometric share of light, not a photograph.
 */
export function createSphereRenderer(gl) {
  const { p, u } = compile(gl, `#version 300 es
in vec2 aCorner; in vec3 aCenter; in float aRadius; in vec3 aColor; in float aEmit;
uniform mat4 uView, uProj;
out vec3 vPos; out vec3 vCenter; out float vRadius; out vec3 vColor; out float vEmit;
void main() {
  vec3 c = (uView * vec4(aCenter, 1.0)).xyz;
  float d = length(c), r = aRadius;
  float s = r * d / sqrt(max(d * d - r * r, 1e-12)) * 1.05;
  vec3 q = c + vec3(aCorner * s, 0.0);
  vPos = q; vCenter = c; vRadius = r; vColor = aColor; vEmit = aEmit;
  gl_Position = uProj * vec4(q, 1.0);
}`, `#version 300 es
precision highp float;
in vec3 vPos; in vec3 vCenter; in float vRadius; in vec3 vColor; in float vEmit;
uniform mat4 uProj; uniform vec3 uLight; uniform float uPoint, uAmbient;
out vec4 frag;
void main() {
  vec3 dir = normalize(vPos);
  float b = dot(dir, vCenter), cc = dot(vCenter, vCenter) - vRadius * vRadius, disc = b * b - cc;
  if (disc < 0.0) discard;
  vec3 h = dir * (b - sqrt(disc));
  vec3 n = (h - vCenter) / vRadius;
  vec3 l = uPoint > 0.5 ? normalize(uLight - vCenter) : normalize(uLight);
  float k = vEmit > 0.5 ? 1.0 : uAmbient + (1.0 - uAmbient) * max(dot(n, l), 0.0);
  vec4 clip = uProj * vec4(h, 1.0);
  gl_FragDepth = clip.z / clip.w * 0.5 + 0.5;
  frag = vec4(vColor * k, 1.0);
}`);
  const vao = gl.createVertexArray(), corner = gl.createBuffer(), inst = gl.createBuffer();
  gl.bindVertexArray(vao);
  gl.bindBuffer(gl.ARRAY_BUFFER, corner); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), gl.STATIC_DRAW);
  const cl = gl.getAttribLocation(p, 'aCorner'); gl.enableVertexAttribArray(cl); gl.vertexAttribPointer(cl, 2, gl.FLOAT, false, 0, 0);
  gl.bindBuffer(gl.ARRAY_BUFFER, inst);
  const stride = 8 * 4;
  for (const [name, size, off] of [['aCenter', 3, 0], ['aRadius', 1, 3], ['aColor', 3, 4], ['aEmit', 1, 7]]) {
    const loc = gl.getAttribLocation(p, name);
    gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, size, gl.FLOAT, false, stride, off * 4); gl.vertexAttribDivisor(loc, 1);
  }
  gl.bindVertexArray(null);
  return {
    draw(camera, spheres, { aspect, light = { direction: [1, 0, 0] }, ambient = 0.08 } = {}) {
      if (!spheres.length) return;
      const t = camera.target, data = new Float32Array(8 * spheres.length);
      spheres.forEach((sp, i) => {
        const c = rgb(sp.color ?? '#cfd8e3');
        data.set([sp.center[0] - t[0], sp.center[1] - t[1], sp.center[2] - t[2], sp.radius, c[0], c[1], c[2], sp.emissive ? 1 : 0], 8 * i);
      });
      const { view, proj } = camera.matrices(aspect ?? gl.drawingBufferWidth / gl.drawingBufferHeight, { relative: true });
      gl.useProgram(p); gl.bindVertexArray(vao);
      gl.bindBuffer(gl.ARRAY_BUFFER, inst); gl.bufferData(gl.ARRAY_BUFFER, data, gl.DYNAMIC_DRAW);
      gl.uniformMatrix4fv(u.uView, false, view); gl.uniformMatrix4fv(u.uProj, false, proj);
      const point = Array.isArray(light);
      const L = point ? light.map((v, k) => v - t[k]) : light.direction;
      // light position (point) or direction (parallel) into view space
      const lv = point
        ? [0, 1, 2].map(r => view[r] * L[0] + view[4 + r] * L[1] + view[8 + r] * L[2] + view[12 + r])
        : [0, 1, 2].map(r => view[r] * L[0] + view[4 + r] * L[1] + view[8 + r] * L[2]);
      gl.uniform3f(u.uLight, ...lv); gl.uniform1f(u.uPoint, point ? 1 : 0); gl.uniform1f(u.uAmbient, ambient);
      gl.disable(gl.BLEND); gl.enable(gl.DEPTH_TEST); gl.depthFunc(gl.LEQUAL); gl.depthMask(true);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, spheres.length);
      gl.bindVertexArray(null);
    },
  };
}

// ---------------------------------------------------------------- lines

/**
 * Thick polylines in 3D (orbits, trajectories, field lines, axes), constant width in screen
 * pixels. update(lines) with lines = [[x,y,z, x,y,z, …], …] (or one flat array), where a line may
 * also be { points, color, opacity } to give it its own colour (e.g. age-faded tracks in one draw);
 * draw(camera, { aspect, width, color, opacity, dash: [on, off] in world units, depthTest }).
 * A line's own opacity multiplies the draw opacity; lines without a colour use the draw colour.
 * Dashed lines are for reference overlays (analytic curves), solid for computed paths.
 */
export function createLineRenderer(gl) {
  const { p, u } = compile(gl, `#version 300 es
in vec2 aCorner; in vec3 aA; in vec3 aB; in vec2 aS; in vec4 aC;
uniform mat4 uVP; uniform vec2 uViewport; uniform float uWidth;
out float vS; out vec4 vC;
void main() {
  vec4 a = uVP * vec4(aA, 1.0), b = uVP * vec4(aB, 1.0);
  vec4 me = aCorner.x < 0.5 ? a : b;
  vec2 sa = a.xy / max(a.w, 1e-6) * uViewport, sb = b.xy / max(b.w, 1e-6) * uViewport;
  vec2 d = sb - sa; float len = length(d);
  vec2 nrm = len > 1e-6 ? vec2(-d.y, d.x) / len : vec2(0.0, 1.0);
  vec2 off = nrm * aCorner.y * uWidth / uViewport;
  gl_Position = vec4(me.xy + off * me.w, me.z, me.w);
  vS = aCorner.x < 0.5 ? aS.x : aS.y; vC = aC;
}`, `#version 300 es
precision highp float;
in float vS; in vec4 vC; uniform float uOpacity; uniform vec2 uDash;
out vec4 frag;
void main() {
  if (uDash.x > 0.0 && mod(vS, uDash.x + uDash.y) > uDash.x) discard;
  float a = vC.a * uOpacity;
  frag = vec4(vC.rgb * a, a);
}`);
  const vao = gl.createVertexArray(), corner = gl.createBuffer(), inst = gl.createBuffer();
  gl.bindVertexArray(vao);
  gl.bindBuffer(gl.ARRAY_BUFFER, corner); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, -1, 0, 1, 1, -1, 1, 1]), gl.STATIC_DRAW);
  const cl = gl.getAttribLocation(p, 'aCorner'); gl.enableVertexAttribArray(cl); gl.vertexAttribPointer(cl, 2, gl.FLOAT, false, 0, 0);
  gl.bindBuffer(gl.ARRAY_BUFFER, inst);
  for (const [name, size, off] of [['aA', 3, 0], ['aB', 3, 3], ['aS', 2, 6], ['aC', 4, 8]]) {
    const loc = gl.getAttribLocation(p, name);
    gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 48, off * 4); gl.vertexAttribDivisor(loc, 1);
  }
  gl.bindVertexArray(null);
  let lines = [];   // rebuilt each draw: the camera target (the float64 origin) may move
  return {
    update(data) { lines = Array.isArray(data) && data.length && typeof data[0] !== 'number' ? data : [data]; },
    draw(camera, { aspect, width = 2, color = '#9fb3c8', opacity = 1, dash = null, depthTest = true } = {}) {
      const t = camera.target;
      // segments relative to the camera target (float64 subtraction), with cumulative length for dashes
      const pts = l => (Array.isArray(l) || ArrayBuffer.isView(l) ? l : l.points);
      let n = 0;
      for (const l of lines) n += Math.max(0, Math.floor(pts(l).length / 3) - 1);
      if (!n) return;
      {
        const buf = new Float32Array(12 * n), base = rgb(color);
        let k = 0;
        for (const line of lines) {
          const l = pts(line), c = line.color ? rgb(line.color) : base, o = line.opacity ?? 1;
          let s = 0;
          for (let i = 0; i + 5 < l.length; i += 3) {
            const a = [l[i] - t[0], l[i + 1] - t[1], l[i + 2] - t[2]], b = [l[i + 3] - t[0], l[i + 4] - t[1], l[i + 5] - t[2]];
            const len = Math.hypot(l[i + 3] - l[i], l[i + 4] - l[i + 1], l[i + 5] - l[i + 2]);
            buf.set([...a, ...b, s, s + len, c[0], c[1], c[2], o], 12 * k++); s += len;
          }
        }
        gl.bindBuffer(gl.ARRAY_BUFFER, inst); gl.bufferData(gl.ARRAY_BUFFER, buf, gl.DYNAMIC_DRAW);
      }
      const dpr = devicePixelRatio || 1;
      const { viewProj } = camera.matrices(aspect ?? gl.drawingBufferWidth / gl.drawingBufferHeight, { relative: true });
      gl.useProgram(p); gl.bindVertexArray(vao);
      gl.uniformMatrix4fv(u.uVP, false, viewProj);
      gl.uniform2f(u.uViewport, gl.drawingBufferWidth / 2, gl.drawingBufferHeight / 2);
      gl.uniform1f(u.uWidth, width * dpr / 2);
      gl.uniform1f(u.uOpacity, opacity);
      gl.uniform2f(u.uDash, dash?.[0] ?? 0, dash?.[1] ?? 0);
      gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
      if (depthTest) { gl.enable(gl.DEPTH_TEST); gl.depthMask(false); } else gl.disable(gl.DEPTH_TEST);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, n);
      gl.depthMask(true);
      gl.bindVertexArray(null);
    },
  };
}

// ---------------------------------------------------------------- solid objects (meshes)

/**
 * Geometry for solid objects, z up, centred on the origin unless stated. Each returns
 * { positions, normals, indices } (Float32Array/Uint32Array) for createMeshRenderer.add().
 *   box(sx, sy, sz) · cylinder(r, h, n) (axis z, from z = 0 to h) · cone(r, h, n) (apex at z = h)
 *   sphere(r, n) · disc(r, n) (in the xy-plane, normal +z) · plane(sx, sy) (xy-plane, normal +z)
 *   arrow(length, { shaft, head, headLength, n }) along +z from the origin (orient with alignZ)
 */
export const geometry = (() => {
  const make = (pos, nrm, idx) => ({ positions: Float32Array.from(pos), normals: Float32Array.from(nrm), indices: Uint32Array.from(idx) });
  const quad = (P, N, I, a, b, c, d, n) => { const o = P.length / 3; P.push(...a, ...b, ...c, ...d); for (let k = 0; k < 4; k++) N.push(...n); I.push(o, o + 1, o + 2, o, o + 2, o + 3); };
  return {
    box(sx = 1, sy = 1, sz = 1) {
      const x = sx / 2, y = sy / 2, z = sz / 2, P = [], N = [], I = [];
      quad(P, N, I, [-x, -y, z], [x, -y, z], [x, y, z], [-x, y, z], [0, 0, 1]);
      quad(P, N, I, [-x, y, -z], [x, y, -z], [x, -y, -z], [-x, -y, -z], [0, 0, -1]);
      quad(P, N, I, [x, -y, -z], [x, y, -z], [x, y, z], [x, -y, z], [1, 0, 0]);
      quad(P, N, I, [-x, y, -z], [-x, -y, -z], [-x, -y, z], [-x, y, z], [-1, 0, 0]);
      quad(P, N, I, [x, y, -z], [-x, y, -z], [-x, y, z], [x, y, z], [0, 1, 0]);
      quad(P, N, I, [-x, -y, -z], [x, -y, -z], [x, -y, z], [-x, -y, z], [0, -1, 0]);
      return make(P, N, I);
    },
    cylinder(r = 0.5, h = 1, n = 32) {
      const P = [], N = [], I = [];
      for (let i = 0; i < n; i++) {
        const a0 = 2 * Math.PI * i / n, a1 = 2 * Math.PI * (i + 1) / n, c0 = Math.cos(a0), s0 = Math.sin(a0), c1 = Math.cos(a1), s1 = Math.sin(a1);
        const o = P.length / 3;
        P.push(r * c0, r * s0, 0, r * c1, r * s1, 0, r * c1, r * s1, h, r * c0, r * s0, h); N.push(c0, s0, 0, c1, s1, 0, c1, s1, 0, c0, s0, 0);
        I.push(o, o + 1, o + 2, o, o + 2, o + 3);
        for (const [zz, nz] of [[h, 1], [0, -1]]) { const q = P.length / 3; P.push(0, 0, zz, r * c0, r * s0, zz, r * c1, r * s1, zz); N.push(0, 0, nz, 0, 0, nz, 0, 0, nz); I.push(...(nz > 0 ? [q, q + 1, q + 2] : [q, q + 2, q + 1])); }
      }
      return make(P, N, I);
    },
    cone(r = 0.5, h = 1, n = 32) {
      const P = [], N = [], I = [], k = r / h;
      for (let i = 0; i < n; i++) {
        const a0 = 2 * Math.PI * i / n, a1 = 2 * Math.PI * (i + 1) / n, am = (a0 + a1) / 2, o = P.length / 3;
        const nn = a => { const l = Math.hypot(1, k); return [Math.cos(a) / l, Math.sin(a) / l, k / l]; };
        P.push(r * Math.cos(a0), r * Math.sin(a0), 0, r * Math.cos(a1), r * Math.sin(a1), 0, 0, 0, h); N.push(...nn(a0), ...nn(a1), ...nn(am));
        I.push(o, o + 1, o + 2);
        const q = P.length / 3; P.push(0, 0, 0, r * Math.cos(a1), r * Math.sin(a1), 0, r * Math.cos(a0), r * Math.sin(a0), 0); N.push(0, 0, -1, 0, 0, -1, 0, 0, -1); I.push(q, q + 1, q + 2);
      }
      return make(P, N, I);
    },
    sphere(r = 0.5, n = 24) {
      const P = [], N = [], I = [];
      for (let j = 0; j <= n; j++) {for (let i = 0; i <= 2 * n; i++) {
        const th = Math.PI * j / n, ph = Math.PI * i / n, v = [Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th)];
        P.push(r * v[0], r * v[1], r * v[2]); N.push(...v);
      }}
      const w = 2 * n + 1;
      for (let j = 0; j < n; j++) for (let i = 0; i < 2 * n; i++) { const a = j * w + i; I.push(a, a + w, a + 1, a + 1, a + w, a + w + 1); }
      return make(P, N, I);
    },
    disc(r = 0.5, n = 48) {
      const P = [0, 0, 0], N = [0, 0, 1], I = [];
      for (let i = 0; i <= n; i++) { const a = 2 * Math.PI * i / n; P.push(r * Math.cos(a), r * Math.sin(a), 0); N.push(0, 0, 1); if (i) I.push(0, i, i + 1); }
      return make(P, N, I);
    },
    plane(sx = 1, sy = 1) { const P = [], N = [], I = []; quad(P, N, I, [-sx / 2, -sy / 2, 0], [sx / 2, -sy / 2, 0], [sx / 2, sy / 2, 0], [-sx / 2, sy / 2, 0], [0, 0, 1]); return make(P, N, I); },
    arrow(length = 1, { shaft = 0.03 * length, head = 0.08 * length, headLength = 0.22 * length, n = 24 } = {}) {
      const hl = Math.min(headLength, length), a = this.cylinder(shaft, length - hl, n), b = this.cone(head, hl, n);
      const off = a.positions.length / 3;
      return {
        positions: Float32Array.from([...a.positions, ...Array.from(b.positions, (v, i) => (i % 3 === 2 ? v + length - hl : v))]),
        normals: Float32Array.from([...a.normals, ...b.normals]),
        indices: Uint32Array.from([...a.indices, ...Array.from(b.indices, i => i + off)]),
      };
    },
  };
})();

/** Rotation matrix (3×3, row-major) taking +z to the direction dir; pass it as transform({ rotate }). */
export function alignZ(dir) {
  const l = Math.hypot(...dir);
  if (!(l > 0)) throw new Error('alignZ: direction must be non-zero');
  const [x, y, z] = dir.map(v => v / l);
  if (z < -1 + 1e-12) return [1, 0, 0, 0, -1, 0, 0, 0, -1];
  // Rodrigues rotation about ẑ × d by the angle between them
  const k = 1 / (1 + z);
  return [1 - x * x * k, -x * y * k, x, -x * y * k, 1 - y * y * k, y, -x, -y, z];
}

/**
 * Point lists for createLineRenderer.update(): each returns one polyline [x,y,z, …] or an array of them.
 *   circle(center, normal, r, n)        closed circle (great circles, cone rims, orbits)
 *   boxEdges([min, max])                the 12 edges of a box, as separate segments
 *   grid(dims, { spacing, origin })     all lattice lines of a 1–3D grid with dims points per axis
 *   arrow(from, to, { head, normal })   shaft and a two-stroke head in the plane containing normal
 */
export const polylines = {
  circle(center = [0, 0, 0], normal = [0, 0, 1], r = 1, n = 96) {
    const R = alignZ(normal), out = [];
    for (let i = 0; i <= n; i++) {
      const a = 2 * Math.PI * i / n, p = [r * Math.cos(a), r * Math.sin(a), 0];
      for (let k = 0; k < 3; k++) out.push(center[k] + R[3 * k] * p[0] + R[3 * k + 1] * p[1] + R[3 * k + 2] * p[2]);
    }
    return out;
  },
  boxEdges([lo, hi]) {
    const c = i => [i & 1 ? hi[0] : lo[0], i & 2 ? hi[1] : lo[1], i & 4 ? hi[2] : lo[2]];
    const out = [];
    for (let i = 0; i < 8; i++) for (const b of [1, 2, 4]) if (!(i & b)) out.push([...c(i), ...c(i | b)]);
    return out;
  },
  grid(dims, { spacing = 1, origin = [0, 0, 0] } = {}) {
    const d = [...dims, 1, 1].slice(0, 3), sp = typeof spacing === 'number' ? [spacing, spacing, spacing] : spacing, out = [];
    for (let ax = 0; ax < dims.length; ax++) {
      const [b, c] = [0, 1, 2].filter(k => k !== ax);
      for (let i = 0; i < d[b]; i++) {for (let j = 0; j < d[c]; j++) {
        const p = [0, 0, 0], q = [0, 0, 0];
        p[b] = q[b] = i; p[c] = q[c] = j; p[ax] = 0; q[ax] = d[ax] - 1;
        out.push([0, 1, 2].map(k => origin[k] + p[k] * sp[k]).concat([0, 1, 2].map(k => origin[k] + q[k] * sp[k])));
      }}
    }
    return out;
  },
  arrow(from, to, { head = 0.2, normal = [0, 0, 1] } = {}) {
    const d = [0, 1, 2].map(k => to[k] - from[k]), L = Math.hypot(...d);
    if (!(L > 0)) return [[...from, ...to]];
    let side = [d[1] * normal[2] - d[2] * normal[1], d[2] * normal[0] - d[0] * normal[2], d[0] * normal[1] - d[1] * normal[0]];
    if (Math.hypot(...side) < 1e-9 * L) side = Math.abs(d[0]) < 0.9 * L ? [0, d[2], -d[1]] : [-d[2], 0, d[0]];
    const sl = Math.hypot(...side), h = head * L;
    const wing = sgn => [0, 1, 2].map(k => to[k] - d[k] / L * h + sgn * side[k] / sl * h * 0.5);
    return [[...from, ...to], [...wing(1), ...to, ...wing(-1)]];
  },
};

/**
 * A colour-bar legend on the 2D canvas for a colormap-mapped quantity, matching the GPU colours.
 * ctx in CSS pixels (the shell's draw ctx). Vertical by default (low at the bottom), or
 * horizontal. Ticks at nice values (decades when log). Returns the bar's rectangle.
 */
export function drawColorBar(ctx, { x, y, width = 12, height = 120, colormap: map = 'viridis', range = [0, 1], log = false, label = '', ticks = 5, horizontal = false, color, format } = {}) {
  const css = name => (typeof getComputedStyle === 'function' ? getComputedStyle(document.documentElement).getPropertyValue(name).trim() : '');
  const text = color ?? (css('--muted') || '#8a94a6');
  const [lo, hi] = log ? range.map(Math.log10) : range;
  const steps = Math.max(2, Math.round(horizontal ? width : height));
  for (let i = 0; i < steps; i++) {
    const c = colormap(map, (i + 0.5) / steps).map(v => Math.round(255 * v));
    ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`;
    if (horizontal) ctx.fillRect(x + i * width / steps, y, width / steps + 0.5, height);
    else ctx.fillRect(x, y + height - (i + 1) * height / steps, width, height / steps + 0.5);
  }
  ctx.strokeStyle = text; ctx.lineWidth = 1; ctx.strokeRect(x + 0.5, y + 0.5, width - 1, height - 1);
  const vals = [];
  if (log) { for (let e = Math.ceil(lo - 1e-9); e <= hi + 1e-9; e++) vals.push(e); }
  else {
    const raw = (hi - lo) / Math.max(1, ticks - 1), mag = 10 ** Math.floor(Math.log10(raw)), step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= raw * 0.999) ?? raw;
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9 * step; v += step) vals.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  }
  const fmt = format ?? (v => (log ? `10${String(v).replace('-', '⁻').replace(/\d/g, d => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d])}` : Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-2 && v !== 0) ? v.toExponential(1) : String(+v.toPrecision(4))));
  ctx.fillStyle = text; ctx.font = '11px system-ui, sans-serif';
  for (const v of vals) {
    const f = (v - lo) / (hi - lo);
    if (horizontal) { const px = x + f * width; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillRect(px, y + height, 1, 3); ctx.fillText(fmt(v), px, y + height + 4); }
    else { const py = y + height - f * height; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillRect(x + width, py, 3, 1); ctx.fillText(fmt(v), x + width + 5, py); }
  }
  if (label) {
    ctx.textAlign = horizontal ? 'center' : 'left'; ctx.textBaseline = 'bottom';
    ctx.fillText(label, horizontal ? x + width / 2 : x, y - 4);
  }
  return { x, y, width, height };
}

/**
 * Model matrix (column-major, float64 array of 16) from translate [x,y,z], rotate (radians about
 * x, then y, then z, or a row-major 3×3 rotation matrix of 9 numbers, e.g. alignZ(dir)) and scale (number or [sx,sy,sz]). compose(a, b) = a·b
 * builds hierarchies (e.g. a seat hanging from a rotating rope).
 */
export function transform({ translate = [0, 0, 0], rotate = [0, 0, 0], scale = 1 } = {}) {
  const s = typeof scale === 'number' ? [scale, scale, scale] : scale;
  let R;
  if (rotate.length === 9) {
    // a 3×3 rotation matrix, row-major (e.g. from alignZ), stored here column by column
    R = [rotate[0], rotate[3], rotate[6], rotate[1], rotate[4], rotate[7], rotate[2], rotate[5], rotate[8]];
  } else {
    const [ax, ay, az] = rotate;
    const cx = Math.cos(ax), sx = Math.sin(ax), cy = Math.cos(ay), sy = Math.sin(ay), cz = Math.cos(az), sz = Math.sin(az);
    // R = Rz · Ry · Rx (rotate about x first), column by column
    R = [cz * cy, sz * cy, -sy, cz * sy * sx - sz * cx, sz * sy * sx + cz * cx, cy * sx, cz * sy * cx + sz * sx, sz * sy * cx - cz * sx, cy * cx];
  }
  return [R[0] * s[0], R[1] * s[0], R[2] * s[0], 0, R[3] * s[1], R[4] * s[1], R[5] * s[1], 0, R[6] * s[2], R[7] * s[2], R[8] * s[2], 0, translate[0], translate[1], translate[2], 1];
}
export function compose(a, b) {
  const out = new Array(16).fill(0);
  for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) for (let k = 0; k < 4; k++) out[c * 4 + r] += a[k * 4 + r] * b[c * 4 + k];
  return out;
}

/**
 * Lit solid objects. add(geometry) → mesh; draw(camera, items, options) with items
 * [{ mesh, matrix (from transform/compose), color, opacity }] and options:
 *   light: { direction: [x,y,z] } (towards the light, e.g. sunDirection() from core/sun.mjs) or a point
 *   ambient (default 0.25) · shadows: { z, color, opacity }: exact hard shadows of every item on the
 *   horizontal plane z = shadows.z for the light (no double-darkening where shadows overlap).
 * Shading is Lambert n·l plus ambient: a readable picture of the lighting geometry, not a
 * photometric rendering. Drawn relative to camera.target in float64 (large coordinates stay sharp).
 */
export function createMeshRenderer(gl) {
  const { p, u } = compile(gl, `#version 300 es
in vec3 aPos; in vec3 aNrm;
uniform mat4 uVP, uModel; uniform mat3 uNrm; uniform float uFlat;
out vec3 vN; out vec3 vW;
void main() { vec4 w = uModel * vec4(aPos, 1.0); vW = w.xyz; vN = uNrm * aNrm; gl_Position = uVP * w; }`, `#version 300 es
precision highp float;
in vec3 vN; in vec3 vW;
uniform vec3 uColor, uLight; uniform float uPoint, uAmbient, uOpacity, uFlat;
out vec4 frag;
void main() {
  if (uFlat > 0.5) { frag = vec4(uColor * uOpacity, uOpacity); return; }
  vec3 n = normalize(vN), l = uPoint > 0.5 ? normalize(uLight - vW) : normalize(uLight);
  if (!gl_FrontFacing) n = -n;
  float k = uAmbient + (1.0 - uAmbient) * max(dot(n, l), 0.0);
  frag = vec4(uColor * k * uOpacity, uOpacity);
}`);
  const meshes = [];
  const inv3T = m => {                         // inverse-transpose of the upper 3×3 (normal matrix)
    const a = m[0], b = m[4], c = m[8], d = m[1], e = m[5], f = m[9], g = m[2], h = m[6], i = m[10];
    const A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g, det = a * A + b * B + c * C;
    if (!det) return new Float32Array([1, 0, 0, 0, 1, 0, 0, 0, 1]);
    const D = -(b * i - c * h), E = a * i - c * g, F = -(a * h - b * g), G = b * f - c * e, H = -(a * f - c * d), I = a * e - b * d;
    return new Float32Array([A / det, B / det, C / det, D / det, E / det, F / det, G / det, H / det, I / det]);
  };
  const relative = (m, t) => { const r = m.slice(); r[12] -= t[0]; r[13] -= t[1]; r[14] -= t[2]; return new Float32Array(r); };
  return {
    add({ positions, normals, indices }) {
      const vao = gl.createVertexArray(); gl.bindVertexArray(vao);
      for (const [name, data] of [['aPos', positions], ['aNrm', normals]]) {
        const b = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
        const loc = gl.getAttribLocation(p, name); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0);
      }
      const ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, indices, gl.STATIC_DRAW);
      gl.bindVertexArray(null);
      const mesh = { vao, count: indices.length };
      meshes.push(mesh);
      return mesh;
    },
    draw(camera, items, { aspect, light = { direction: [0.4, -0.3, 1] }, ambient = 0.25, shadows = null } = {}) {
      const t = camera.target, { viewProj } = camera.matrices(aspect ?? gl.drawingBufferWidth / gl.drawingBufferHeight, { relative: true });
      const point = Array.isArray(light), L = point ? light.map((v, k) => v - t[k]) : light.direction;
      gl.useProgram(p);
      gl.uniformMatrix4fv(u.uVP, false, viewProj); gl.uniform3f(u.uLight, ...L); gl.uniform1f(u.uPoint, point ? 1 : 0); gl.uniform1f(u.uAmbient, ambient);
      gl.enable(gl.DEPTH_TEST); gl.depthFunc(gl.LEQUAL);
      const opaque = items.filter(it => (it.opacity ?? 1) >= 1), glass = items.filter(it => (it.opacity ?? 1) < 1);
      const one = (it, model, color, opacity, flat) => {
        gl.uniformMatrix4fv(u.uModel, false, model); gl.uniformMatrix3fv(u.uNrm, false, inv3T(it.matrix));
        gl.uniform3f(u.uColor, ...rgb(color)); gl.uniform1f(u.uOpacity, opacity); gl.uniform1f(u.uFlat, flat ? 1 : 0);
        gl.bindVertexArray(it.mesh.vao); gl.drawElements(gl.TRIANGLES, it.mesh.count, gl.UNSIGNED_INT, 0);
      };
      gl.disable(gl.BLEND); gl.depthMask(true);
      for (const it of opaque) one(it, relative(it.matrix, t), it.color ?? '#cfd8e3', 1, false);
      if (shadows && !point && L[2] > 1e-6) {
        // project every vertex along the light onto z = z₀: x' = x − (z − z₀)·l/l_z (exact for parallel light)
        const z0 = shadows.z ?? 0, S = [1, 0, 0, 0, 0, 1, 0, 0, -L[0] / L[2], -L[1] / L[2], 0, 0, z0 * L[0] / L[2], z0 * L[1] / L[2], z0, 1];
        gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.depthMask(false);
        gl.enable(gl.STENCIL_TEST); gl.clear(gl.STENCIL_BUFFER_BIT); gl.stencilFunc(gl.EQUAL, 0, 0xff); gl.stencilOp(gl.KEEP, gl.KEEP, gl.INCR);
        gl.enable(gl.POLYGON_OFFSET_FILL); gl.polygonOffset(-1, -1);
        const op = shadows.opacity ?? 0.35;
        for (const it of items) if (!it.noShadow) one(it, relative(compose(S, it.matrix), t), shadows.color ?? '#000000', op, true);
        gl.disable(gl.POLYGON_OFFSET_FILL); gl.disable(gl.STENCIL_TEST); gl.depthMask(true);
      }
      if (glass.length) {
        gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.depthMask(false);
        for (const it of glass) one(it, relative(it.matrix, t), it.color ?? '#9fd0ff', it.opacity, false);
        gl.depthMask(true);
      }
      gl.bindVertexArray(null);
    },
  };
}
