// Shared real linear algebra and nonlinear equations. Dense matrices are arrays of rows;
// vectors may be ordinary arrays or Float64Array. Inputs are not modified.

const EPS = Number.EPSILON;
const finiteVector = (v, n, who) => {
  if (!v || v.length !== n || !Array.from(v).every(Number.isFinite)) throw new Error(`${who}: expected ${n} finite values`);
};
const matrixShape = (A, who) => {
  const m = A?.length, n = A?.[0]?.length;
  if (!(m > 0 && n > 0)) throw new Error(`${who}: matrix must be non-empty`);
  for (const row of A) finiteVector(row, n, who);
  return [m, n];
};
const norm = v => {
  let result = 0;
  for (const value of v) result = Math.hypot(result, value);
  return result;
};
const dot = (a, b) => {
  let sum = 0, correction = 0;
  for (let i = 0; i < a.length; i++) {
    const term = a[i] * b[i], next = sum + term;
    correction += Math.abs(sum) >= Math.abs(term) ? (sum - next) + term : (term - next) + sum;
    sum = next;
  }
  return sum + correction;
};
const tolerance = (value, who) => {
  if (!(Number.isFinite(value) && value >= 0)) throw new Error(`${who}: tolerance must be finite and non-negative`);
};
const budget = (value, who) => {
  if (!(Number.isSafeInteger(value) && value > 0)) throw new Error(`${who}: budget must be a positive integer`);
};
const finiteResult = (value, who) => {
  if (!Number.isFinite(value)) throw new Error(`${who}: numerical range exceeded; rescale the problem`);
  return value;
};

/** Rectangular matrix-vector product, with compensated summation in each row. */
export function matrixVector(A, x) {
  const [m, n] = matrixShape(A, 'matrixVector');
  finiteVector(x, n, 'matrixVector');
  return Float64Array.from({ length: m }, (_, i) => finiteResult(dot(A[i], x), 'matrixVector'));
}

/**
 * Solve square A x = b by row equilibration and Gaussian elimination with partial pivoting.
 * pivotTolerance is relative to each equilibrated row (default 64 epsilon).
 * Throws for singular/numerically rank-deficient input; no inverse is formed.
 */
export function solveLinear(A, b, { pivotTolerance = 64 * EPS } = {}) {
  const [m, n] = matrixShape(A, 'solveLinear');
  if (m !== n) throw new Error('solveLinear: matrix must be square');
  finiteVector(b, n, 'solveLinear'); tolerance(pivotTolerance, 'solveLinear');
  const M = [], rhs = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let scale = 0;
    for (const value of A[i]) scale = Math.max(scale, Math.abs(value));
    if (scale === 0) throw new Error('solveLinear: singular matrix (zero row)');
    M.push(Float64Array.from(A[i], value => value / scale));
    rhs[i] = finiteResult(b[i] / scale, 'solveLinear');
  }
  for (let k = 0; k < n; k++) {
    let pivot = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[pivot][k])) pivot = i;
    if (!(Math.abs(M[pivot][k]) > pivotTolerance)) throw new Error('solveLinear: singular or numerically rank-deficient matrix');
    if (pivot !== k) { [M[k], M[pivot]] = [M[pivot], M[k]]; [rhs[k], rhs[pivot]] = [rhs[pivot], rhs[k]]; }
    for (let i = k + 1; i < n; i++) {
      const factor = M[i][k] / M[k][k];
      M[i][k] = 0;
      for (let j = k + 1; j < n; j++) M[i][j] = finiteResult(M[i][j] - factor * M[k][j], 'solveLinear');
      rhs[i] = finiteResult(rhs[i] - factor * rhs[k], 'solveLinear');
    }
  }
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let value = rhs[i];
    for (let j = i + 1; j < n; j++) value -= M[i][j] * x[j];
    x[i] = finiteResult(value / M[i][i], 'solveLinear');
  }
  return x;
}

/**
 * Full-column-rank least squares (m >= n), using Householder QR with column pivoting.
 * Columns and rhs are scaled independently. Returns { x, residualNorm, rank, permutation }.
 * Rank-deficient systems throw: this routine is not a pseudoinverse/minimum-norm solver.
 */
export function leastSquares(A, b, { rankTolerance = 64 * EPS } = {}) {
  const [m, n] = matrixShape(A, 'leastSquares');
  if (m < n) throw new Error('leastSquares: need at least as many rows as columns');
  finiteVector(b, m, 'leastSquares'); tolerance(rankTolerance, 'leastSquares');
  const scales = new Float64Array(n), permutation = Array.from({ length: n }, (_, i) => i);
  let bScale = 0;
  for (const value of b) bScale = Math.max(bScale, Math.abs(value));
  if (bScale === 0) bScale = 1;
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < m; i++) scales[j] = Math.max(scales[j], Math.abs(A[i][j]));
    if (scales[j] === 0) throw new Error('leastSquares: numerically rank-deficient matrix (zero column)');
  }
  const R = A.map(row => Float64Array.from(row, (value, j) => value / scales[j]));
  const y = Float64Array.from(b, value => value / bScale);
  const columnNorm = (j, start) => { let value = 0; for (let i = start; i < m; i++) value = Math.hypot(value, R[i][j]); return value; };
  let reference = 0;
  for (let j = 0; j < n; j++) reference = Math.max(reference, columnNorm(j, 0));
  for (let k = 0; k < n; k++) {
    let pivot = k;
    for (let j = k + 1; j < n; j++) if (columnNorm(j, k) > columnNorm(pivot, k)) pivot = j;
    if (!(columnNorm(pivot, k) > rankTolerance * reference)) throw new Error('leastSquares: numerically rank-deficient matrix');
    if (pivot !== k) {
      for (const row of R) [row[k], row[pivot]] = [row[pivot], row[k]];
      [permutation[k], permutation[pivot]] = [permutation[pivot], permutation[k]];
    }
    const v = Float64Array.from({ length: m - k }, (_, i) => R[k + i][k]);
    const alpha = (v[0] >= 0 ? -1 : 1) * norm(v);
    v[0] -= alpha;
    const vNorm = norm(v);
    for (let i = 0; i < v.length; i++) v[i] /= vNorm;
    for (let j = k; j < n; j++) {
      let projection = 0;
      for (let i = 0; i < v.length; i++) projection += v[i] * R[k + i][j];
      for (let i = 0; i < v.length; i++) R[k + i][j] -= 2 * v[i] * projection;
    }
    let projection = 0;
    for (let i = 0; i < v.length; i++) projection += v[i] * y[k + i];
    for (let i = 0; i < v.length; i++) y[k + i] -= 2 * v[i] * projection;
    R[k][k] = alpha;
    for (let i = k + 1; i < m; i++) R[i][k] = 0;
  }
  const z = new Float64Array(n), x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let value = y[i];
    for (let j = i + 1; j < n; j++) value -= R[i][j] * z[j];
    z[i] = value / R[i][i];
  }
  for (let i = 0; i < n; i++) x[permutation[i]] = finiteResult(z[i] * bScale / scales[permutation[i]], 'leastSquares');
  const residual = Float64Array.from(matrixVector(A, x), (value, i) => value - b[i]);
  return { x, residualNorm: finiteResult(norm(residual), 'leastSquares'), rank: n, permutation };
}

/**
 * Real symmetric dense eigensystem by cyclic Jacobi rotations. Values sorted ascending;
 * vectors[i] is the unit eigenvector of values[i]. O(n^3) per sweep, for modest dense systems.
 * rtol controls the off-diagonal Frobenius norm relative to the original Frobenius norm.
 */
export function eigenSymmetric(A, { rtol = 1e-12, symmetryTolerance = 64 * EPS, maxSweeps = 100 } = {}) {
  const [m, n] = matrixShape(A, 'eigenSymmetric');
  if (m !== n) throw new Error('eigenSymmetric: matrix must be square');
  tolerance(rtol, 'eigenSymmetric'); tolerance(symmetryTolerance, 'eigenSymmetric'); budget(maxSweeps, 'eigenSymmetric');
  let scale = 0;
  for (const row of A) for (const value of row) scale = Math.max(scale, Math.abs(value));
  const divisor = scale || 1;
  const M = A.map(row => Float64Array.from(row, value => value / divisor));
  const V = Array.from({ length: n }, (_row, i) => Float64Array.from({ length: n }, (_column, j) => i === j ? 1 : 0));
  for (let i = 0; i < n; i++) { for (let j = i + 1; j < n; j++) {
    if (Math.abs(M[i][j] - M[j][i]) > symmetryTolerance) throw new Error('eigenSymmetric: matrix must be symmetric');
    M[i][j] = M[j][i] = (M[i][j] + M[j][i]) / 2;
  } }
  let reference = 0;
  for (const row of M) for (const value of row) reference = Math.hypot(reference, value);
  const target = rtol * reference;
  const offNorm = () => { let value = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) value = Math.hypot(value, Math.SQRT2 * M[i][j]); return value; };
  let sweeps = 0;
  while (offNorm() > target) {
    if (sweeps === maxSweeps) throw new Error('eigenSymmetric: sweep budget exhausted before convergence');
    for (let p = 0; p < n; p++) { for (let q = p + 1; q < n; q++) {
      const a = M[p][q];
      if (Math.abs(a) <= target / (n * n)) continue;
      const delta = (M[q][q] - M[p][p]) / 2;
      const t = delta === 0 ? 1 : a / (delta + Math.sign(delta) * Math.hypot(delta, a));
      const c = 1 / Math.hypot(1, t), s = t * c;
      M[p][p] -= t * a; M[q][q] += t * a; M[p][q] = M[q][p] = 0;
      for (let k = 0; k < n; k++) {
        if (k !== p && k !== q) {
          const u = M[k][p], v = M[k][q];
          M[k][p] = M[p][k] = c * u - s * v;
          M[k][q] = M[q][k] = s * u + c * v;
        }
        const u = V[k][p], v = V[k][q];
        V[k][p] = c * u - s * v; V[k][q] = s * u + c * v;
      }
    } }
    sweeps++;
  }
  const order = Array.from({ length: n }, (_, i) => i).sort((i, j) => M[i][i] - M[j][j]);
  return {
    values: Float64Array.from(order, i => finiteResult(M[i][i] * divisor, 'eigenSymmetric')),
    vectors: order.map(i => Float64Array.from(V, row => row[i])),
    sweeps, offDiagonalNorm: finiteResult(offNorm() * divisor, 'eigenSymmetric'),
  };
}

/**
 * Cholesky factor of a symmetric positive-definite matrix: A = L Lᵀ, L lower triangular (arrays of
 * rows). Throws when A is not positive definite (a pivot ≤ rtol·max diagonal), e.g. a noisy
 * correlation matrix at large time separation: drop that operator or that time.
 */
export function cholesky(A, { symmetryTolerance = 64 * EPS, rtol = 64 * EPS } = {}) {
  const [m, n] = matrixShape(A, 'cholesky');
  if (m !== n) throw new Error('cholesky: matrix must be square');
  let scale = 0;
  for (let i = 0; i < n; i++) {
    scale = Math.max(scale, Math.abs(A[i][i]));
    for (let j = 0; j < i; j++) if (Math.abs(A[i][j] - A[j][i]) > symmetryTolerance * Math.max(1, Math.abs(A[i][j]))) throw new Error('cholesky: matrix must be symmetric');
  }
  const L = Array.from({ length: n }, () => new Float64Array(n));
  for (let j = 0; j < n; j++) {
    let d = A[j][j];
    for (let k = 0; k < j; k++) d -= L[j][k] * L[j][k];
    if (!(d > rtol * scale)) throw new Error(`cholesky: matrix is not positive definite (pivot ${j}: ${d})`);
    L[j][j] = Math.sqrt(d);
    for (let i = j + 1; i < n; i++) {
      let v = (A[i][j] + A[j][i]) / 2;
      for (let k = 0; k < j; k++) v -= L[i][k] * L[j][k];
      L[i][j] = v / L[j][j];
    }
  }
  return L;
}

/**
 * Generalised symmetric-definite eigenproblem A v = λ B v (B positive definite), by Cholesky
 * B = LLᵀ and the ordinary problem L⁻¹AL⁻ᵀ w = λ w, v = L⁻ᵀw. Values ascending; vectors[i] is
 * B-normalised (vᵀBv = 1). The variational (GEVP) method for correlation matrices is
 * A = C(t), B = C(t₀): the largest λ = e^{−E₀(t−t₀)} gives the ground state.
 */
export function eigenGeneralized(A, B, options = {}) {
  const [n] = matrixShape(A, 'eigenGeneralized');
  if (matrixShape(B, 'eigenGeneralized')[0] !== n || A[0].length !== n) throw new Error('eigenGeneralized: A and B must be square and the same size');
  const L = cholesky(B, options);
  // X = L⁻¹ A (forward substitution column by column), then C = X L⁻ᵀ = L⁻¹ (L⁻¹A)ᵀ
  const lower = M => {
    const X = Array.from({ length: n }, () => new Float64Array(n));
    for (let c = 0; c < n; c++) {
      for (let i = 0; i < n; i++) {
        let v = M[i][c];
        for (let k = 0; k < i; k++) v -= L[i][k] * X[k][c];
        X[i][c] = v / L[i][i];
      }
    }
    return X;
  };
  const X = lower(A);
  const C = lower(Array.from({ length: n }, (_a, i) => Float64Array.from({ length: n }, (_b, j) => X[j][i])));
  for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) C[i][j] = C[j][i] = (C[i][j] + C[j][i]) / 2;
  const { values, vectors } = eigenSymmetric(C, { rtol: options.rtol ?? 1e-12 });
  // v = L⁻ᵀ w (back substitution)
  return {
    values,
    vectors: vectors.map(w => {
      const v = new Float64Array(n);
      for (let i = n - 1; i >= 0; i--) { let x = w[i]; for (let k = i + 1; k < n; k++) x -= L[k][i] * v[k]; v[i] = x / L[i][i]; }
      return v;
    }),
  };
}

/** Dense matrix or a linear callback (x, out) -> out/array/undefined. Callback must not mutate x. */
function linearOperator(A, n, who) {
  if (typeof A !== 'function') {
    const [m, columns] = matrixShape(A, who);
    if (m !== n || columns !== n) throw new Error(`${who}: operator must be ${n} by ${n}`);
    return (x, out) => { for (let i = 0; i < n; i++) out[i] = dot(A[i], x); finiteVector(out, n, who); return out; };
  }
  return (x, out) => {
    out.fill(NaN);
    const result = A(x, out) ?? out;
    finiteVector(result, n, who);
    if (result !== out) out.set(result);
    return out;
  };
}

/**
 * Preconditioned conjugate gradient for symmetric positive-definite A. Dense or matrix-free.
 * preconditioner applies M ~= A^-1, also SPD. True residual is checked before success.
 * Returns { x, converged, iterations, residualNorm }; exhausted budgets return converged:false.
 */
export function conjugateGradient(A, b, { x0, preconditioner, rtol = 1e-10, atol = 0, maxIterations = Math.max(20, 4 * b.length) } = {}) {
  const n = b?.length;
  if (!(n > 0)) throw new Error('conjugateGradient: rhs must be non-empty');
  finiteVector(b, n, 'conjugateGradient');
  if (x0) finiteVector(x0, n, 'conjugateGradient');
  tolerance(rtol, 'conjugateGradient'); tolerance(atol, 'conjugateGradient'); budget(maxIterations, 'conjugateGradient');
  const apply = linearOperator(A, n, 'conjugateGradient');
  const precondition = preconditioner ? linearOperator(preconditioner, n, 'conjugateGradient preconditioner') : (r, z) => { z.set(r); return z; };
  let bScale = 0;
  for (const value of b) bScale = Math.max(bScale, Math.abs(value));
  const bNorm = finiteResult(norm(b), 'conjugateGradient');
  const target = Math.max(rtol * bNorm, atol);
  if (bScale === 0) return { x: new Float64Array(n), converged: true, iterations: 0, residualNorm: 0 };
  const y = x0 ? Float64Array.from(x0, value => value / bScale) : new Float64Array(n);
  finiteVector(y, n, 'conjugateGradient');
  const rhs = Float64Array.from(b, value => value / bScale), Ay = new Float64Array(n), r = new Float64Array(n), z = new Float64Array(n), p = new Float64Array(n), Ap = new Float64Array(n);
  apply(y, Ay);
  for (let i = 0; i < n; i++) r[i] = rhs[i] - Ay[i];
  const report = iterations => {
    const x = Float64Array.from(y, value => finiteResult(value * bScale, 'conjugateGradient'));
    apply(x, Ay);
    const residual = Float64Array.from(b, (value, i) => value - Ay[i]);
    const residualNorm = finiteResult(norm(residual), 'conjugateGradient');
    return { x, converged: residualNorm <= target, iterations, residualNorm, residual };
  };
  if (norm(r) <= target / bScale) { const result = report(0); if (result.converged) return result; }
  precondition(r, z); p.set(z);
  let rz = finiteResult(dot(r, z), 'conjugateGradient');
  if (!(rz > 0)) throw new Error('conjugateGradient: non-positive preconditioned residual (need SPD preconditioner)');
  for (let iteration = 1; iteration <= maxIterations; iteration++) {
    apply(p, Ap);
    const curvature = finiteResult(dot(p, Ap), 'conjugateGradient');
    if (!(curvature > 0)) throw new Error('conjugateGradient: non-positive curvature (need SPD operator)');
    const alpha = rz / curvature;
    for (let i = 0; i < n; i++) { y[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
    finiteVector(y, n, 'conjugateGradient'); finiteVector(r, n, 'conjugateGradient');
    let restart = false;
    if (norm(r) <= target / bScale || iteration === maxIterations) {
      const result = report(iteration);
      if (result.converged || iteration === maxIterations) return result;
      for (let i = 0; i < n; i++) r[i] = result.residual[i] / bScale;
      restart = true;
    }
    precondition(r, z);
    const next = finiteResult(dot(r, z), 'conjugateGradient');
    if (!(next > 0)) throw new Error('conjugateGradient: non-positive preconditioned residual (need SPD preconditioner)');
    const beta = restart ? 0 : next / rz;
    for (let i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
    rz = next;
  }
}

/**
 * Damped Newton solve f(x)=0, with an analytic Jacobian or central finite differences.
 * Use dimensionless/scaled equations; xScale sets coordinate scales for differencing.
 * Success means ||f|| <= max(atol, rtol*||f(x0)||), not a bound on distance to a root.
 */
export function newtonSolve(f, x0, { jacobian, xScale, rtol = 1e-10, atol = 1e-12, maxIterations = 50, maxBacktracks = 25 } = {}) {
  const n = x0?.length;
  if (typeof f !== 'function' || !(n > 0)) throw new Error('newtonSolve: function and non-empty initial state required');
  finiteVector(x0, n, 'newtonSolve');
  if (jacobian !== undefined && typeof jacobian !== 'function') throw new Error('newtonSolve: jacobian must be a function');
  tolerance(rtol, 'newtonSolve'); tolerance(atol, 'newtonSolve'); budget(maxIterations, 'newtonSolve'); budget(maxBacktracks, 'newtonSolve');
  const scales = xScale ? Float64Array.from(xScale) : new Float64Array(n).fill(1);
  finiteVector(scales, n, 'newtonSolve');
  if (!scales.every(value => value > 0)) throw new Error('newtonSolve: xScale must be positive');
  let x = Float64Array.from(x0), residual = Float64Array.from(f(x));
  finiteVector(residual, n, 'newtonSolve');
  let residualNorm = finiteResult(norm(residual), 'newtonSolve');
  const target = Math.max(atol, rtol * residualNorm);
  const result = (iterations, reason) => ({ x, residual, residualNorm, iterations, converged: residualNorm <= target, reason });
  for (let iteration = 0; iteration < maxIterations; iteration++) {
    if (residualNorm <= target) return result(iteration, 'residual tolerance');
    let J;
    if (jacobian) J = jacobian(x);
    else {
      J = Array.from({ length: n }, () => new Float64Array(n));
      for (let j = 0; j < n; j++) {
        const h = Math.cbrt(EPS) * Math.max(Math.abs(x[j]), scales[j]);
        const plus = x.slice(), minus = x.slice(); plus[j] += h; minus[j] -= h;
        if (!Number.isFinite(h) || plus[j] === minus[j]) throw new Error('newtonSolve: derivative step cannot be resolved');
        const fp = Float64Array.from(f(plus)), fm = f(minus);
        finiteVector(fp, n, 'newtonSolve derivative'); finiteVector(fm, n, 'newtonSolve derivative');
        for (let i = 0; i < n; i++) J[i][j] = (fp[i] - fm[i]) / (plus[j] - minus[j]);
      }
    }
    const delta = solveLinear(J, Float64Array.from(residual, value => -value));
    let accepted = false, alpha = 1;
    for (let attempt = 0; attempt < maxBacktracks; attempt++, alpha /= 2) {
      const trial = Float64Array.from(x, (value, i) => value + alpha * delta[i]);
      if (!trial.every(Number.isFinite)) continue;
      const next = f(trial);
      if (!next || next.length !== n) throw new Error('newtonSolve: function must return one residual per coordinate');
      if (!Array.from(next).every(Number.isFinite)) continue;
      const nextNorm = norm(next);
      if (nextNorm <= (1 - 1e-4 * alpha) * residualNorm || nextNorm <= target) {
        x = trial; residual = Float64Array.from(next); residualNorm = nextNorm; accepted = true; break;
      }
    }
    if (!accepted) return result(iteration, 'line search stalled');
  }
  return result(maxIterations, residualNorm <= target ? 'residual tolerance' : 'iteration budget exhausted');
}

/**
 * Matrix exponential e^A of a real square matrix (Higham 2005: scaling and squaring with the
 * degree-13 Padé approximant). x(t) = e^{At}x₀ solves any linear system x' = Ax: decay chains,
 * first-order kinetics, compartment models, continuous-time Markov chains. Accuracy is normwise
 * (relative to ‖e^A‖); check tiny components against a structured solution when they matter.
 * Returns the matrix as an array of Float64Array rows.
 */
export function expm(A) {
  const [m, n] = matrixShape(A, 'expm');
  if (m !== n) throw new Error('expm: matrix must be square');
  const b = [64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000, 10559470521600, 670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1];
  const mul = (X, Y) => X.map(row => { const out = new Float64Array(n); for (let k = 0; k < n; k++) { const r = row[k]; if (r) for (let j = 0; j < n; j++) out[j] += r * Y[k][j]; } return out; });
  const comb = terms => Array.from({ length: n }, (_, i) => { const out = new Float64Array(n); for (const [c, X] of terms) for (let j = 0; j < n; j++) out[j] += c * (X === null ? (i === j ? 1 : 0) : X[i][j]); return out; });
  const norm1 = Math.max(...Array.from({ length: n }, (_, j) => A.reduce((s, row) => s + Math.abs(row[j]), 0)));
  const s = Math.max(0, Math.ceil(Math.log2(norm1 / 5.371920351148152)));
  const S = A.map(row => Float64Array.from(row, v => v / 2 ** s));
  const A2 = mul(S, S), A4 = mul(A2, A2), A6 = mul(A4, A2);
  const U = mul(S, comb([[1, mul(A6, comb([[b[13], A6], [b[11], A4], [b[9], A2]]))], [b[7], A6], [b[5], A4], [b[3], A2], [b[1], null]]));
  const V = comb([[1, mul(A6, comb([[b[12], A6], [b[10], A4], [b[8], A2]]))], [b[6], A6], [b[4], A4], [b[2], A2], [b[0], null]]);
  const P = comb([[1, V], [-1, U]]), Q = comb([[1, V], [1, U]]);
  // R = P⁻¹Q column by column, then undo the scaling by repeated squaring
  const cols = Array.from({ length: n }, (_, j) => solveLinear(P, Q.map(row => row[j])));
  let R = Array.from({ length: n }, (_, i) => Float64Array.from(cols, col => col[i]));
  for (let k = 0; k < s; k++) R = mul(R, R);
  for (const row of R) for (const v of row) finiteResult(v, 'expm');
  return R;
}
