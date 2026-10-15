// QCD running coupling: β coefficients against the SU(3) closed forms and the published β3,
// the integrator against the 1-loop closed form and the exact implicit 2-loop solution,
// reversibility through flavour thresholds, physical sanity at m_τ, and the PDG input table.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { suN, betaCoefficients, runAlphaS, alphaS1Loop, lambdaQCD1Loop, exactTwoLoop } from '../core/qcd.mjs';
import { PDG } from '../core/constants.mjs';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const ZETA3 = 1.2020569031595942;

test('β0–β2 match the SU(3) closed forms for nf = 0..6', () => {
  for (let nf = 0; nf <= 6; nf++) {
    const [b0, b1, b2] = betaCoefficients(nf);
    assert.ok(rel(b0, 11 - 2 * nf / 3) < 1e-14);
    assert.ok(rel(b1, 102 - 38 * nf / 3) < 1e-14);
    assert.ok(rel(b2, 2857 / 2 - 5033 * nf / 18 + 325 * nf * nf / 54) < 1e-13);
  }
});

test('Casimir-form β3 reproduces the SU(3) polynomial (van Ritbergen–Vermaseren–Larin 1997)', () => {
  const su3 = nf => (149753 / 6 + 3564 * ZETA3) - (1078361 / 162 + 6508 / 27 * ZETA3) * nf
    + (50065 / 162 + 6472 / 81 * ZETA3) * nf * nf + 1093 / 729 * nf ** 3;
  for (let nf = 0; nf <= 6; nf++) assert.ok(rel(betaCoefficients(nf)[3], su3(nf)) < 1e-13, `nf = ${nf}`);
  assert.ok(Math.abs(betaCoefficients(5)[3] - 4826.156) < 1e-3, 'β3(nf = 5) ≈ 4826.156');
  assert.ok(Math.abs(betaCoefficients(3)[3] - 12090.378) < 1e-3, 'β3(nf = 3) ≈ 12090.378');
  // the explicit SU(3) group object gives the same numbers; Casimirs alone leave β3 unknown
  assert.deepEqual(betaCoefficients(4, { group: suN(3) }), betaCoefficients(4));
  assert.ok(Number.isNaN(betaCoefficients(4, { group: { CA: 3, CF: 4 / 3, TF: 1 / 2 } })[3]));
  // SU(2) pure gauge: β0 = 22/3, β1 = 136/3
  const [s0, s1] = betaCoefficients(0, { group: suN(2) });
  assert.ok(rel(s0, 22 / 3) < 1e-15 && rel(s1, 136 / 3) < 1e-15);
});

test('1-loop running equals the closed form at fixed nf and piecewise through thresholds', () => {
  const run = runAlphaS({ loops: 1, nf: 5 });
  const L = lambdaQCD1Loop(0.118, 91.1876, 5);
  for (const Q of [5, 20, 91.1876, 500, 1e4]) assert.ok(rel(run.at(Q), alphaS1Loop(Q, L, 5)) < 1e-10, `Q = ${Q}`);
  assert.ok(rel(alphaS1Loop(91.1876, L, 5), 0.118) < 1e-14, 'Λ inverse round trip');
  // with c/b thresholds the 1-loop solution is a chain of closed forms with continuous matching
  const thr = runAlphaS({ loops: 1, thresholds: { c: 1.3, b: 4.2 } });
  const L4 = lambdaQCD1Loop(alphaS1Loop(4.2, L, 5), 4.2, 4);
  const L3 = lambdaQCD1Loop(alphaS1Loop(1.3, L4, 4), 1.3, 3);
  assert.ok(rel(thr.at(2), alphaS1Loop(2, L4, 4)) < 1e-10);
  assert.ok(rel(thr.at(1.1), alphaS1Loop(1.1, L3, 3)) < 1e-10);
});

test('2-loop running equals the exact implicit solution', () => {
  for (const nf of [3, 5]) {
    const run = runAlphaS({ loops: 2, nf, alpha0: 0.2, mu0: 10 });
    for (const Q of [1.5, 3, 10, 50, 1e3, 1e5]) assert.ok(rel(run.at(Q), exactTwoLoop(0.2, 10, Q, nf)) < 1e-10, `nf = ${nf}, Q = ${Q}`);
  }
});

test('round trip M_Z → low scale → M_Z through thresholds closes to 1e-10', () => {
  for (const loops of [2, 4]) {
    for (const Qlow of [1.5, 1.0]) {
      const down = runAlphaS({ loops }).at(Qlow);
      const up = runAlphaS({ loops, alpha0: down, mu0: Qlow }).at(PDG.mZ.value);
      assert.ok(rel(up, PDG.alphaS_MZ.value) < 1e-10, `${loops} loops via ${Qlow} GeV: ${up}`);
    }
  }
});

test('α_s is continuous at thresholds, switches nf there, and the cache is self-consistent', () => {
  const run = runAlphaS();
  const mb = PDG.mb_mb.value;
  assert.equal(run.nf(mb * (1 - 1e-9)), 4);
  assert.equal(run.nf(mb), 5);
  assert.equal(run.nf(1), 3);
  assert.equal(run.nf(200), 6);
  assert.ok(rel(run.at(mb * (1 - 1e-12)), run.at(mb * (1 + 1e-12))) < 1e-10);
  const a1 = run.at(3.3), fresh = runAlphaS().at(3.3);   // after many cached queries vs a cold start
  for (const Q of [2, 7, 1.4, 50, 3.31]) run.at(Q);
  assert.ok(rel(run.at(3.3), a1) < 1e-11 && rel(a1, fresh) < 1e-11);
});

test('α_s(m_τ) at 4 loops from 0.1180 lies in the τ-determination band', () => {
  const a = runAlphaS({ alpha0: 0.1180 }).at(1.77686);
  assert.ok(a > 0.30 && a < 0.34, `α_s(m_τ) = ${a}`);
});

test('α_s decreases monotonically with Q (asymptotic freedom)', () => {
  const run = runAlphaS();
  let prev = Infinity;
  for (let Q = 0.9; Q < 5000; Q *= 1.15) { const a = run.at(Q); assert.ok(a < prev, `Q = ${Q}`); prev = a; }
});

test('errors: Landau pole, non-positive Q, bad options', () => {
  const run = runAlphaS();
  assert.throws(() => run.at(0.3), /Landau pole/);
  assert.throws(() => run.at(0), /Q must be/);
  assert.throws(() => run.at(-1), /Q must be/);
  assert.throws(() => runAlphaS({ loops: 1, nf: 3 }).at(0.05), /Landau pole/);
  assert.throws(() => alphaS1Loop(0.1, 0.2, 3), /Landau pole/);
  assert.throws(() => exactTwoLoop(0.118, 91.1876, 0.01, 5), /Landau pole/);
  assert.throws(() => runAlphaS({ loops: 5 }), /loops/);
  assert.throws(() => runAlphaS({ group: { CA: 3, CF: 4 / 3, TF: 1 / 2 } }), /β3 unknown/);
  assert.ok(runAlphaS({ loops: 3, group: { CA: 3, CF: 4 / 3, TF: 1 / 2 } }).at(10) > 0);
});

test('PDG entries all have finite value and error, a unit string and a source', () => {
  assert.ok(Object.isFrozen(PDG));
  for (const [k, v] of Object.entries(PDG)) {
    assert.ok(Number.isFinite(v.value) && Number.isFinite(v.error) && v.error >= 0, k);
    assert.equal(typeof v.unit, 'string', k);
    assert.match(v.source, /^(PDG 2024|CODATA 2022)/, k);
  }
  assert.ok(rel(PDG.mZ.value, 91.1876) < 1e-12 && rel(PDG.alphaS_MZ.value, 0.118) < 1e-12);
});
