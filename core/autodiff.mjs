// Automatic differentiation: exact derivatives (to rounding) of code, not of a formula.
// Plain ES module (browser + Node), no dependencies.
//
//   Dual, dual                forward mode: a + b·ε with ε² = 0; nests (Dual of Duals) for
//                             second derivatives (hyper-dual numbers)
//   Var (via gradient)        reverse mode: a tape of elementary operations swept backwards
//   add sub mul div neg pow sqrt cbrt exp expm1 log log1p sin cos tan asin acos atan atan2
//   sinh cosh tanh abs hypot min max square sum dot value
//                             arithmetic on numbers, Duals and tape Vars alike
//   derivative, secondDerivative, gradientForward, jacobian, jvp, gradient, hessian,
//   hessianVectorProduct     drivers
//
// Write the function once with these operations (JavaScript has no operator overloading):
//   const energy = r => sub(pow(r, -12), pow(r, -6));
// then gradient(energy, …) gives value and gradient by the reverse sweep (cost ~ 3–5 function
// evaluations whatever the number of variables); derivative/jvp give one directional
// derivative per forward pass. Accuracy: rounding error of the derivative expressions only,
// no truncation error. Branches (if, Math.max on value(x)) are differentiated piecewise: the
// result is the derivative of the branch taken; at a kink it is one-sided.
//
// Nesting rule: a tape Var is always the outermost layer and Duals live inside it (forward-over-
// reverse, used by hessianVectorProduct). Duals created by the drivers below are tagged only
// by nesting depth; do not close over a Dual from an outer derivative() call inside an inner
// one ("perturbation confusion") — use secondDerivative/hessian, which build the nesting
// themselves.

// ---------------------------------------------------------------- types

/** Dual number re + du·ε (ε² = 0). re and du may themselves be Duals (nested/hyper-dual). */
export class Dual {
  constructor(re, du = 0) { this.re = re; this.du = du; }
  toString() { return `${this.re} + ${this.du}ε`; }
}

/** Shorthand: new Dual(x, dx). */
export const dual = (x, dx = 1) => new Dual(x, dx);

class Tape { constructor() { this.nodes = []; } }

/** Reverse-mode variable: a node on a tape. Created by gradient()/hessianVectorProduct(). */
export class Var {
  constructor(tape, val, parents = null, partials = null) {
    this.tape = tape; this.val = val; this.parents = parents; this.partials = partials;
    this.index = tape.nodes.length; tape.nodes.push(this);
  }
  toString() { return `Var(${this.val})`; }
}

const isVar = x => x instanceof Var;
const isDual = x => x instanceof Dual;
const re = x => (isDual(x) ? x.re : x);
const du = x => (isDual(x) ? x.du : 0);

/** Primal (plain number) value of a number, Dual or Var, however deeply nested. */
export function value(x) {
  let v = x;
  for (;;) {
    if (typeof v === 'number') return v;
    if (isVar(v)) v = v.val;
    else if (isDual(v)) v = v.re;
    else throw new Error('autodiff: expected a number, Dual or Var');
  }
}

const tapeOf = (a, b) => {
  const t = isVar(a) ? a.tape : b.tape;
  if (isVar(a) && isVar(b) && a.tape !== b.tape) throw new Error('autodiff: Vars from different tapes cannot be mixed');
  return t;
};

// ---------------------------------------------------------------- arithmetic

/** a + b */
export function add(a, b) {
  if (typeof a === 'number' && typeof b === 'number') return a + b;
  if (isVar(a) || isVar(b)) {
    const t = tapeOf(a, b), av = isVar(a) ? a.val : a, bv = isVar(b) ? b.val : b;
    return new Var(t, add(av, bv), [a, b], [1, 1]);
  }
  return new Dual(add(re(a), re(b)), add(du(a), du(b)));
}
/** a − b */
export function sub(a, b) {
  if (typeof a === 'number' && typeof b === 'number') return a - b;
  if (isVar(a) || isVar(b)) {
    const t = tapeOf(a, b), av = isVar(a) ? a.val : a, bv = isVar(b) ? b.val : b;
    return new Var(t, sub(av, bv), [a, b], [1, -1]);
  }
  return new Dual(sub(re(a), re(b)), sub(du(a), du(b)));
}
/** a · b */
export function mul(a, b) {
  if (typeof a === 'number' && typeof b === 'number') return a * b;
  if (isVar(a) || isVar(b)) {
    const t = tapeOf(a, b), av = isVar(a) ? a.val : a, bv = isVar(b) ? b.val : b;
    return new Var(t, mul(av, bv), [a, b], [bv, av]);
  }
  return new Dual(mul(re(a), re(b)), add(mul(re(a), du(b)), mul(du(a), re(b))));
}
/** a / b */
export function div(a, b) {
  if (typeof a === 'number' && typeof b === 'number') return a / b;
  if (isVar(a) || isVar(b)) {
    const t = tapeOf(a, b), av = isVar(a) ? a.val : a, bv = isVar(b) ? b.val : b, q = div(av, bv);
    return new Var(t, q, [a, b], [div(1, bv), neg(div(q, bv))]);
  }
  const q = div(re(a), re(b));
  return new Dual(q, div(sub(du(a), mul(q, du(b))), re(b)));
}
/** −a */
export function neg(a) { return typeof a === 'number' ? -a : mul(-1, a); }
/** a² */
export function square(a) { return mul(a, a); }

const unary = (f0, df) => {
  const g = x => {
    if (typeof x === 'number') return f0(x);
    if (isVar(x)) return new Var(x.tape, g(x.val), [x], [df(x.val)]);
    if (isDual(x)) return new Dual(g(x.re), mul(df(x.re), x.du));
    throw new Error('autodiff: expected a number, Dual or Var');
  };
  return g;
};

/** √x (derivative 1/(2√x); infinite at 0). */
export const sqrt = unary(Math.sqrt, x => div(0.5, sqrt(x)));
/** ∛x */
export const cbrt = unary(Math.cbrt, x => div(1 / 3, square(cbrt(x))));
/** eˣ */
export const exp = unary(Math.exp, x => exp(x));
/** eˣ − 1, accurate for small x */
export const expm1 = unary(Math.expm1, x => exp(x));
/** ln x */
export const log = unary(Math.log, x => div(1, x));
/** ln(1 + x), accurate for small x */
export const log1p = unary(Math.log1p, x => div(1, add(1, x)));
/** sin x */
export const sin = unary(Math.sin, x => cos(x));
/** cos x */
export const cos = unary(Math.cos, x => neg(sin(x)));
/** tan x */
export const tan = unary(Math.tan, x => div(1, square(cos(x))));
/** arcsin x (|x| < 1 for a finite derivative) */
export const asin = unary(Math.asin, x => div(1, sqrt(sub(1, square(x)))));
/** arccos x */
export const acos = unary(Math.acos, x => neg(div(1, sqrt(sub(1, square(x))))));
/** arctan x */
export const atan = unary(Math.atan, x => div(1, add(1, square(x))));
/** sinh x */
export const sinh = unary(Math.sinh, x => cosh(x));
/** cosh x */
export const cosh = unary(Math.cosh, x => sinh(x));
/** tanh x */
// sech² x = 4e^(−2|x|)/(1 + e^(−2|x|))² retains the derivative when tanh rounds to ±1.
export const tanh = unary(Math.tanh, x => {
  const t = exp(mul(value(x) < 0 ? 2 : -2, x));
  return div(mul(4, t), square(add(1, t)));
});
/** |x| (derivative sign(x); 0 at x = 0 is a subgradient convention, not a classical derivative) */
export const abs = unary(Math.abs, x => Math.sign(value(x)));

/**
 * aᵇ. A number exponent uses b·a^{b−1} (valid for any real a where aᵇ is defined, including
 * a ≤ 0 with integer b); a differentiated exponent uses exp(b·ln a) and needs a > 0.
 */
export function pow(a, b) {
  if (typeof a === 'number' && typeof b === 'number') return a ** b;
  if (typeof b === 'number') {
    if (b === 0) return isVar(a) ? new Var(a.tape, pow(a.val, 0), [a], [0]) : isDual(a) ? new Dual(pow(a.re, 0), 0) : 1;
    if (isVar(a)) return new Var(a.tape, pow(a.val, b), [a], [mul(b, pow(a.val, b - 1))]);
    return new Dual(pow(a.re, b), mul(mul(b, pow(a.re, b - 1)), a.du));
  }
  return exp(mul(b, log(a)));
}

/** atan2(y, x), the angle of the point (x, y). */
export function atan2(y, x) {
  if (typeof y === 'number' && typeof x === 'number') return Math.atan2(y, x);
  const scale = Math.max(Math.abs(value(x)), Math.abs(value(y))) || 1;
  if (isVar(y) || isVar(x)) {
    const t = tapeOf(y, x), yv = isVar(y) ? y.val : y, xv = isVar(x) ? x.val : x;
    const xs = div(xv, scale), ys = div(yv, scale), r2 = add(square(xs), square(ys));
    return new Var(t, atan2(yv, xv), [y, x], [div(div(xs, r2), scale), neg(div(div(ys, r2), scale))]);
  }
  const xs = div(re(x), scale), ys = div(re(y), scale), r2 = add(square(xs), square(ys));
  return new Dual(atan2(re(y), re(x)), sub(mul(div(div(xs, r2), scale), du(y)), mul(div(div(ys, r2), scale), du(x))));
}

/** √(Σ aᵢ²) without differentiating through intermediate overflow (scaled by the primal max). */
export function hypot(...args) {
  if (args.every(a => typeof a === 'number')) return Math.hypot(...args);
  const scale = Math.max(...args.map(a => Math.abs(value(a)))) || 1;
  let s = 0;
  for (const a of args) s = add(s, square(div(a, scale)));
  return mul(scale, sqrt(s));
}

/** min(a, b) / max(a, b) by primal comparison (derivative of the branch taken). */
export const min = (a, b) => (value(a) <= value(b) ? a : b);
export const max = (a, b) => (value(a) >= value(b) ? a : b);
/** Σ xᵢ */
export function sum(xs) { let s = 0; for (const x of xs) s = add(s, x); return s; }
/** Σ aᵢbᵢ */
export function dot(a, b) { let s = 0; for (let i = 0; i < a.length; i++) s = add(s, mul(a[i], b[i])); return s; }

// ---------------------------------------------------------------- drivers

const vectorIn = (x, who) => {
  if (!x || !(x.length > 0)) throw new Error(`${who}: x must be a non-empty vector`);
  const v = Array.from(x);
  if (!v.every(Number.isFinite)) throw new Error(`${who}: x must be finite`);
  return v;
};
const fnIn = (f, who) => { if (typeof f !== 'function') throw new Error(`${who}: f must be a function`); };

/** f′(x) of a scalar function by one forward pass: { value, derivative }. */
export function derivative(f, x) {
  fnIn(f, 'derivative');
  const y = f(new Dual(x, 1));
  return { value: value(re(y)), derivative: isDual(y) ? value(y.du) : 0 };
}

/** f(x), f′(x), f″(x) by one pass with nested duals (x + ε₁ + ε₂): { value, first, second }. */
export function secondDerivative(f, x) {
  fnIn(f, 'secondDerivative');
  const y = f(new Dual(new Dual(x, 1), new Dual(1, 0)));
  const outer = isDual(y) ? y : new Dual(y, 0), inner = isDual(outer.re) ? outer.re : new Dual(outer.re, 0);
  const d = isDual(outer.du) ? outer.du : new Dual(outer.du, 0);
  return { value: value(inner.re), first: value(inner.du), second: value(d.du) };
}

/** Directional derivative ∇f(x)·v (or J(x)·v for vector f) in one forward pass: { value, jvp }. */
export function jvp(f, x, v) {
  fnIn(f, 'jvp');
  const xs = vectorIn(x, 'jvp'), vs = vectorIn(v, 'jvp');
  if (vs.length !== xs.length) throw new Error('jvp: x and v must have equal length');
  const y = f(xs.map((xi, i) => new Dual(xi, vs[i])));
  if (Array.isArray(y) || ArrayBuffer.isView(y)) return { value: Array.from(y, yi => value(re(yi))), jvp: Array.from(y, yi => value(du(yi))) };
  return { value: value(re(y)), jvp: value(du(y)) };
}

/** Gradient by n forward passes: { value, gradient } (prefer gradient() for large n). */
export function gradientForward(f, x) {
  fnIn(f, 'gradientForward');
  const xs = vectorIn(x, 'gradientForward'), g = new Float64Array(xs.length);
  let v = NaN;
  for (let i = 0; i < xs.length; i++) {
    const y = f(xs.map((xj, j) => new Dual(xj, i === j ? 1 : 0)));
    v = value(re(y)); g[i] = value(du(y));
  }
  return { value: v, gradient: g };
}

/** Jacobian of a vector function F: Rⁿ → Rᵐ by n forward passes: { value, jacobian (m rows) }. */
export function jacobian(F, x) {
  fnIn(F, 'jacobian');
  const xs = vectorIn(x, 'jacobian'), n = xs.length;
  let val = null, J = null;
  for (let j = 0; j < n; j++) {
    const y = Array.from(F(xs.map((xk, k) => new Dual(xk, k === j ? 1 : 0))));
    if (!J) { val = y.map(yi => value(re(yi))); J = y.map(() => new Float64Array(n)); }
    y.forEach((yi, i) => { J[i][j] = value(du(yi)); });
  }
  return { value: val, jacobian: J };
}

const backward = (tape, out, seed) => {
  const adj = new Array(tape.nodes.length).fill(0);
  adj[out.index] = seed;
  for (let k = tape.nodes.length - 1; k >= 0; k--) {
    const node = tape.nodes[k], a = adj[k];
    if (!node.parents || a === 0) continue;
    for (let p = 0; p < node.parents.length; p++) {
      const parent = node.parents[p];
      if (isVar(parent)) adj[parent.index] = add(adj[parent.index], mul(a, node.partials[p]));
    }
  }
  return adj;
};

/**
 * Value and gradient of a scalar f: Rⁿ → R by reverse mode (one taped forward pass and one
 * backward sweep; memory ∝ number of operations). f receives an array of Vars and must combine
 * them only with this module's operations (Math.* on a Var throws or returns NaN).
 * Returns { value, gradient: Float64Array }.
 */
export function gradient(f, x) {
  fnIn(f, 'gradient');
  const xs = vectorIn(x, 'gradient'), tape = new Tape(), vars = xs.map(v => new Var(tape, v));
  const y = f(vars);
  if (!isVar(y)) {
    if (typeof y !== 'number') throw new Error('gradient: f must return a scalar');
    return { value: y, gradient: new Float64Array(xs.length) };
  }
  const adj = backward(tape, y, 1);
  return { value: value(y), gradient: Float64Array.from(vars, v => value(adj[v.index])) };
}

/**
 * Hessian-vector product H(x)·v by forward-over-reverse (one taped pass on Duals x + vε):
 * { value, gradient, hvp }. Cost ~ a few gradient evaluations, independent of n.
 */
export function hessianVectorProduct(f, x, v) {
  fnIn(f, 'hessianVectorProduct');
  const xs = vectorIn(x, 'hessianVectorProduct'), vs = vectorIn(v, 'hessianVectorProduct');
  if (vs.length !== xs.length) throw new Error('hessianVectorProduct: x and v must have equal length');
  const tape = new Tape(), vars = xs.map((xi, i) => new Var(tape, new Dual(xi, vs[i])));
  const y = f(vars), n = xs.length;
  if (!isVar(y)) return { value: value(y), gradient: new Float64Array(n), hvp: new Float64Array(n) };
  const adj = backward(tape, y, 1);
  return { value: value(y), gradient: Float64Array.from(vars, w => value(re(adj[w.index]))), hvp: Float64Array.from(vars, w => value(du(adj[w.index]))) };
}

/**
 * Full Hessian of f: Rⁿ → R by n forward-over-reverse passes (rows H·e_j): { value, gradient,
 * hessian (rows) }. Exactly symmetric up to rounding; cost ~ n gradients.
 */
export function hessian(f, x) {
  const xs = vectorIn(x, 'hessian'), n = xs.length, H = [];
  let val = NaN, grad = null;
  for (let j = 0; j < n; j++) {
    const r = hessianVectorProduct(f, xs, xs.map((_, k) => (k === j ? 1 : 0)));
    val = r.value; grad = r.gradient; H.push(r.hvp);
  }
  return { value: val, gradient: grad, hessian: H };
}
