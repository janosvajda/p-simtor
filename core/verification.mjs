// Shared verification gates for the CLI and browser. Reports are evidence only
// when every check is explicit and the tested settings still match the simulation.

export function validateModuleReport(report) {
  if (!report || !Array.isArray(report.checks) || !report.checks.length) throw new Error('verification: need a non-empty checks array');
  for (const [i, check] of report.checks.entries()) {
    const where = `verification check ${i + 1}`;
    if (!check || typeof check.label !== 'string' || !check.label.trim()) throw new Error(`${where}: label is required`);
    if (typeof check.pass !== 'boolean') throw new Error(`${where}: pass must be explicitly true or false`);
    if (typeof check.value === 'number') {
      if (!Number.isFinite(check.value)) throw new Error(`${where}: measured value must be finite`);
    } else if (typeof check.value !== 'string' || !check.value.trim() || /NaN|Infinity|∞|^\s*—\s*$/.test(check.value)) {
      throw new Error(`${where}: measured value must be finite and present`);
    }
    if (check.state != null && check.state !== 'warn') throw new Error(`${where}: optional state must be 'warn'`);
  }
  return report;
}

/** Reasons a gravity report cannot support a successful verification claim. */
export function gravityReportProblems(report) {
  const problems = [];
  const finiteError = (value, limit, label) => {
    if (!(Number.isFinite(value) && value >= 0)) problems.push(`${label} is missing or non-finite`);
    else if (value > limit) problems.push(`${label} ${value.toExponential(2)} exceeds ${limit.toExponential(0)}`);
  };
  const c = report?.conservation, cv = report?.convergence;
  finiteError(c?.dE, report?.settings?.softening > 0 ? 1e-4 : 1e-8, 'energy error');
  finiteError(c?.dP, 1e-10, 'momentum drift');
  finiteError(c?.dL, 1e-10, 'angular momentum drift');
  if (!(Number.isFinite(c?.h) && c.h > 0)) problems.push('step parameter is missing or non-finite');
  if (!(Number.isFinite(report?.duration) && report.duration > 0)) problems.push('verification duration is invalid');
  if (!(Number.isFinite(cv?.observedOrder) && Number.isFinite(cv?.expectedOrder) && cv.expectedOrder > 0)) problems.push('convergence order is unresolved or non-finite');
  else if (Math.abs(cv.observedOrder - cv.expectedOrder) > 0.6) problems.push('observed convergence order differs from theory');
  if (!(Number.isFinite(cv?.horizon) && cv.horizon >= 0 && cv.horizon <= report.duration && typeof cv.converged === 'boolean')) problems.push('convergence horizon is invalid');
  finiteError(cv?.finalDifference, Infinity, 'final refinement difference');
  if (report?.determinism?.identical !== true || report?.frameIndependence?.identical !== true) problems.push('reproducibility check failed');
  if (!Array.isArray(report?.benchmarks) || !report.benchmarks.length) problems.push('benchmarks are missing');
  else {for (const b of report.benchmarks) {
    const name = b?.name ?? 'unnamed benchmark';
    finiteError(b?.dE, 1e-9, `${name}: energy error`);
    if (b && Object.hasOwn(b, 'positionError')) finiteError(b.positionError, 1e-3, `${name}: position error`);
    else if (b && Object.hasOwn(b, 'returnError')) finiteError(b.returnError, 1e-6, `${name}: return error`);
    else if (b?.matchesLiterature !== true) problems.push(`benchmark failed: ${name}`);
  }}
  return problems;
}

/** Match solver settings, allowing only round-off in the computed step parameter. */
export function gravityReportMatches(report, settings) {
  if (!report?.settings || report.integrator !== settings.integrator) return false;
  for (const key of ['integrator', 'G', 'softening', 'stepsPerCrossing']) {
    if (report.settings[key] !== settings[key]) return false;
  }
  const h = report.conservation?.h;
  return Number.isFinite(h) && h > 0 && Number.isFinite(settings.h) && Math.abs(h - settings.h) <= 1e-12 * Math.abs(h);
}
