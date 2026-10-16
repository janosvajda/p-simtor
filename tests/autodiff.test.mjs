// Automatic differentiation checked against analytic derivatives (every elementary function,
// first and second order), analytic Lennard-Jones forces, exact Hessians and finite differences.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import * as AD from '../core/autodiff.mjs';
import { finiteGradient, finiteHessian } from '../core/optimize.mjs';
import { rng } from '../core/numerics.mjs';

const { add, sub, mul, div, pow, sqrt, exp, log, sin, cos, square, sum, dot, value } = AD;
const close = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${a} vs ${b} (tolerance ${tol})`);

// [name, ad function, exact f′, exact f″, sample points]
const CASES = [
  ['sqrt', AD.sqrt, x => 0.5 / Math.sqrt(x), x => -0.25 * x ** -1.5, [0.3, 2, 17]],
  ['cbrt', AD.cbrt, x => x ** (-2 / 3) / 3, x => -2 / 9 * x ** (-5 / 3), [0.3, 2, 17]],
  ['exp', AD.exp, Math.exp, Math.exp, [-3, 0, 1.7]],
  ['expm1', AD.expm1, Math.exp, Math.exp, [-1e-8, 0.5]],
  ['log', AD.log, x => 1 / x, x => -1 / (x * x), [0.1, 1, 50]],
  ['log1p', AD.log1p, x => 1 / (1 + x), x => -1 / (1 + x) ** 2, [-0.5, 1e-9, 3]],
  ['sin', AD.sin, Math.cos, x => -Math.sin(x), [-2, 0.4, 3]],
  ['cos', AD.cos, x => -Math.sin(x), x => -Math.cos(x), [-2, 0.4, 3]],
  ['tan', AD.tan, x => 1 / Math.cos(x) ** 2, x => 2 * Math.tan(x) / Math.cos(x) ** 2, [-1.2, 0.3, 1]],
  ['asin', AD.asin, x => 1 / Math.sqrt(1 - x * x), x => x / (1 - x * x) ** 1.5, [-0.9, 0, 0.6]],
  ['acos', AD.acos, x => -1 / Math.sqrt(1 - x * x), x => -x / (1 - x * x) ** 1.5, [-0.9, 0, 0.6]],
  ['atan', AD.atan, x => 1 / (1 + x * x), x => -2 * x / (1 + x * x) ** 2, [-4, 0.2, 9]],
  ['sinh', AD.sinh, Math.cosh, Math.sinh, [-2, 0.5]],
  ['cosh', AD.cosh, Math.sinh, Math.cosh, [-2, 0.5]],
  ['tanh', AD.tanh, x => 1 / Math.cosh(x) ** 2, x => -2 * Math.tanh(x) / Math.cosh(x) ** 2, [-2, 0.5]],
  ['abs', AD.abs, Math.sign, () => 0, [-2, 0.5]],
  ['pow 3.5', x => pow(x, 3.5), x => 3.5 * x ** 2.5, x => 8.75 * x ** 1.5, [0.2, 3]],
  ['pow −3', x => pow(x, -3), x => -3 * x ** -4, x => 12 * x ** -5, [-1.5, 0.7]],
  ['x^x', x => pow(x, x), x => x ** x * (Math.log(x) + 1), x => x ** x * ((Math.log(x) + 1) ** 2 + 1 / x), [0.5, 2]],
  ['quotient', x => div(sin(x), add(x, 2)), x => (Math.cos(x) * (x + 2) - Math.sin(x)) / (x + 2) ** 2,
    x => (-Math.sin(x) * (x + 2) ** 2 - 2 * (Math.cos(x) * (x + 2) - Math.sin(x))) / (x + 2) ** 3, [0.1, 1.3]],
];

test('forward and reverse mode reproduce analytic first derivatives of every elementary function', () => {
  for (const [name, f, d1, , xs] of CASES) {
    for (const x of xs) {
      const fwd = AD.derivative(f, x), rev = AD.gradient(v => f(v[0]), [x]);
      close(fwd.derivative, d1(x), 1e-14);
      close(rev.gradient[0], d1(x), 1e-14);
      close(fwd.value, value(f(x)), 1e-15);
      assert.ok(Number.isFinite(fwd.derivative), name);
    }
  }
  // two-argument functions
  const g = AD.gradient(([y, x]) => AD.atan2(y, x), [0.7, -1.3]);
  close(g.gradient[0], -1.3 / (0.7 ** 2 + 1.3 ** 2), 1e-15); close(g.gradient[1], -0.7 / (0.7 ** 2 + 1.3 ** 2), 1e-15);
  const h = AD.gradient(v => AD.hypot(...v), [3e200, 4e200]);   // no overflow in x²
  close(h.value, 5e200, 1e-15); close(h.gradient[0], 0.6, 1e-15); close(h.gradient[1], 0.8, 1e-15);
});

test('nested (hyper-dual) forward mode reproduces analytic second derivatives', () => {
  for (const [name, f, d1, d2, xs] of CASES) {
    for (const x of xs) {
      const r = AD.secondDerivative(f, x);
      close(r.first, d1(x), 1e-14);
      close(r.second, d2(x), 1e-13);
      assert.ok(Number.isFinite(r.second), name);
    }
  }
});

test('reverse-mode gradient of a Lennard-Jones cluster equals minus the analytic forces', () => {
  const rand = rng(7), N = 13, x0 = [];
  // jittered cubic lattice near the LJ minimum spacing (no overlapping pairs: a fair FD check)
  for (let i = 0; i < N; i++) x0.push(...[i % 3, Math.floor(i / 3) % 3, Math.floor(i / 9)].map(c => 1.12 * c + 0.2 * (rand() - 0.5)));
  const energy = x => {
    let E = 0;
    for (let i = 0; i < N; i++) {for (let j = i + 1; j < N; j++) {
      const r2 = add(add(square(sub(x[3 * i], x[3 * j])), square(sub(x[3 * i + 1], x[3 * j + 1]))), square(sub(x[3 * i + 2], x[3 * j + 2])));
      const s6 = pow(r2, -3);
      E = add(E, mul(4, sub(square(s6), s6)));
    }}
    return E;
  };
  const { value: E, gradient: g } = AD.gradient(energy, x0);
  let Eexact = 0;
  const F = new Float64Array(3 * N);
  for (let i = 0; i < N; i++) {for (let j = i + 1; j < N; j++) {
    const d = [0, 1, 2].map(k => x0[3 * i + k] - x0[3 * j + k]), r2 = d[0] ** 2 + d[1] ** 2 + d[2] ** 2, s6 = r2 ** -3;
    Eexact += 4 * (s6 * s6 - s6);
    const c = 24 * (2 * s6 * s6 - s6) / r2;
    for (let k = 0; k < 3; k++) { F[3 * i + k] += c * d[k]; F[3 * j + k] -= c * d[k]; }
  }}
  close(E, Eexact, 1e-13);
  const scale = Math.max(...F.map(Math.abs));
  for (let k = 0; k < 3 * N; k++) assert.ok(Math.abs(g[k] + F[k]) <= 1e-13 * scale, `component ${k}`);
  // forward mode agrees to rounding; finite differences to their truncation error (steep r⁻¹² derivatives)
  const fwd = AD.gradientForward(energy, x0).gradient, fd = finiteGradient(x => value(energy(Array.from(x))), x0, { order: 4 });
  for (let k = 0; k < 3 * N; k++) { assert.ok(Math.abs(fwd[k] - g[k]) <= 1e-13 * scale); assert.ok(Math.abs(fd[k] - g[k]) <= 1e-6 * scale, `${fd[k]} vs ${g[k]}`); }
});

test('Hessian, Hessian-vector product, Jacobian and JVP match exact results', () => {
  const rosen = ([a, b]) => add(mul(100, square(sub(b, square(a)))), square(sub(1, a)));
  const x = [0.3, 0.7], H = AD.hessian(rosen, x);
  const exact = [[1200 * 0.09 - 400 * 0.7 + 2, -400 * 0.3], [-400 * 0.3, 200]];
  for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) close(H.hessian[i][j], exact[i][j], 1e-14);
  const dir = [0.6, -1.1], hv = AD.hessianVectorProduct(rosen, x, dir).hvp;
  for (let i = 0; i < 2; i++) close(hv[i], exact[i][0] * dir[0] + exact[i][1] * dir[1], 1e-14);
  // a 6-variable function against a finite-difference Hessian
  const rand = rng(2), y = Array.from({ length: 6 }, () => rand() - 0.5);
  const f6 = z => add(mul(exp(dot(z, z)), sin(sum(z))), log(add(3, mul(z[0], z[5]))));
  const H6 = AD.hessian(f6, y).hessian, Hfd = finiteHessian(z => value(f6(Array.from(z))), y);
  for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) { close(H6[i][j], H6[j][i], 1e-14); close(H6[i][j], Hfd[i][j], 1e-6); }
  // polar → Cartesian Jacobian and a directional derivative
  const polar = ([r, t]) => [mul(r, cos(t)), mul(r, sin(t))];
  const J = AD.jacobian(polar, [2, 0.4]).jacobian;
  close(J[0][0], Math.cos(0.4), 1e-15); close(J[0][1], -2 * Math.sin(0.4), 1e-15);
  close(J[1][0], Math.sin(0.4), 1e-15); close(J[1][1], 2 * Math.cos(0.4), 1e-15);
  const d = AD.jvp(polar, [2, 0.4], [1, 0.5]).jvp;
  close(d[0], J[0][0] + 0.5 * J[0][1], 1e-15); close(d[1], J[1][0] + 0.5 * J[1][1], 1e-15);
  // a constant function has zero gradient; bad input is rejected
  assert.deepEqual(Array.from(AD.gradient(() => 5, [1, 2]).gradient), [0, 0]);
  assert.throws(() => AD.gradient(() => 'x', [1]), /scalar/);
  assert.throws(() => AD.gradient(v => v[0], [NaN]), /finite/);
  close(value(sqrt(new AD.Dual(4, 1))), 2, 0);
});
