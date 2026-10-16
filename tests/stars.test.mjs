import { test } from 'node:test';
import assert from 'node:assert/strict';
import { G, Msun, mp, kB } from '../core/constants.mjs';
import { brent } from '../core/numerics.mjs';
import {
  laneEmden, polytrope, polytropeMass, chandrasekharMass, whiteDwarf, whiteDwarfOfMass, semfBinding, mostStableZ,
  bindingFromMass, AME2020_BA, sahaRatio, sahaHydrogen, halfIonisationTemperature, epsPP, epsCNO,
  stefanBoltzmannLuminosity, effectiveTemperature, eddingtonLuminosity, kelvinHelmholtzTime, nuclearTime,
  mainSequenceLuminosity, mainSequenceRadius, mainSequenceLifetime, R_SUN, L_SUN, TEFF_SUN, SIGMA_T, SSM_CENTRE, YEAR,
} from '../core/stars.mjs';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const close = (a, b, tol, msg = '') => assert.ok(rel(a, b) <= tol, `${msg} ${a} vs ${b} (rel tol ${tol})`);

test('Lane–Emden exact solutions: n = 0 (1 − ξ²/6), n = 1 (sin ξ/ξ), n = 5 ((1 + ξ²/3)^{-1/2})', () => {
  const l0 = laneEmden(0), l1 = laneEmden(1), l5 = laneEmden(5);
  close(l0.xi1, Math.sqrt(6), 1e-12); close(l0.omega, 2 * Math.sqrt(6), 1e-12); close(l0.densityRatio, 1, 1e-12);
  close(l1.xi1, Math.PI, 1e-12); close(l1.omega, Math.PI, 1e-10); close(l1.densityRatio, Math.PI ** 2 / 3, 1e-10);
  for (const xi of [0.5, 1.7, 2.4]) {
    assert.ok(Math.abs(l0.theta(xi) - (1 - xi * xi / 6)) < 1e-9);
    assert.ok(Math.abs(l1.theta(xi) - Math.sin(xi) / xi) < 1e-9, `n=1 ξ=${xi}`);
    assert.ok(Math.abs(l1.dtheta(xi) - (xi * Math.cos(xi) - Math.sin(xi)) / (xi * xi)) < 1e-9);
  }
  assert.ok(Number.isNaN(l0.theta(2.5)), 'outside the star');
  assert.equal(l5.xi1, Infinity);
  for (const xi of [1, 10, 40]) assert.ok(Math.abs(l5.theta(xi) - 1 / Math.sqrt(1 + xi * xi / 3)) < 1e-9, `n=5 ξ=${xi}`);
  close(l5.omega, 50 ** 3 / 3 * (1 + 2500 / 3) ** -1.5, 1e-9, 'n=5: −ξ²θ′ → √3 as ξ → ∞');
  assert.throws(() => laneEmden(-1), /n must/);
});

test('Lane–Emden tabulated constants (Chandrasekhar 1939 Table 4; Horedt 2004 Table 2.4.1): n = 1.5, 2, 3, 4', () => {
  const tab = { 1.5: [3.65375, 2.71406, 5.99071], 2: [4.35287, 2.41105, 11.40254], 3: [6.89685, 2.01824, 54.1825], 4: [14.97155, 1.79723, 622.408] };
  for (const [n, [xi1, om, dr]] of Object.entries(tab)) {
    const l = laneEmden(Number(n));
    close(l.xi1, xi1, 2e-6, `ξ1 n=${n}`); close(l.omega, om, 3e-6, `ω n=${n}`); close(l.densityRatio, dr, 2e-6, `ρc/ρ̄ n=${n}`);
  }
  close(laneEmden(3).Wn, 11.05, 1e-3, 'Pc = 11.05 GM²/R⁴ for n = 3');
});

test('Lane–Emden convergence: ξ1 and ω for n = 3 settle under tolerance refinement', () => {
  const a = laneEmden(3, { rtol: 1e-9 }), b = laneEmden(3, { rtol: 1e-11 }), ref = laneEmden(3, { rtol: 1e-13 });
  assert.ok(Math.abs(b.omega - ref.omega) < Math.abs(a.omega - ref.omega) / 10, 'error falls with rtol');
  assert.ok(Math.abs(b.xi1 - ref.xi1) < 1e-9);
});

test('polytrope: mass integral of the profile = M, hydrostatic equilibrium dP/dr = −Gmρ/r², n = 1 analytic profile', () => {
  const s = polytrope({ n: 3, M: Msun, R: R_SUN });
  close(polytropeMass(s), Msun, 1e-8);
  close(s.rhoC / s.rhoMean, 54.18, 1e-4, 'ρc/ρ̄ = 54.18 for n = 3');
  for (const x of [0.1, 0.3, 0.7]) {
    const r = x * R_SUN, dr = 1e-5 * R_SUN, p = s.profile(r);
    const dPdr = (s.profile(r + dr).P - s.profile(r - dr).P) / (2 * dr);
    close(dPdr, -G * p.m * p.rho / (r * r), 1e-6, `r = ${x} R`);
  }
  const s1 = polytrope({ n: 1, M: Msun, R: R_SUN }), r = 0.4 * R_SUN, k = Math.PI * r / R_SUN;
  close(s1.profile(r).rho, s1.rhoC * Math.sin(k) / k, 1e-9);
  assert.throws(() => polytrope({ n: 5, M: Msun, R: R_SUN }), /n must/);
});

test('APPROXIMATE n = 3 Sun vs Standard Solar Model BP2000: centre within a factor 2.1', () => {
  const s = polytrope({ n: 3, M: Msun, R: R_SUN, mu: 0.62 });
  for (const k of ['Tc', 'rhoC', 'Pc']) {
    const r = s[k] / SSM_CENTRE[k];
    assert.ok(r > 1 / 2.1 && r < 2.1, `${k}: ${s[k]} vs ${SSM_CENTRE[k]}`);
  }
  assert.ok(s.Tc > 1e7 && s.Tc < 1.6e7, 'central temperature ~1e7 K: hot enough for the pp chain');
});

test('Chandrasekhar mass: 5.825/μe² M☉ with m_u (1.4563 M☉ for μe = 2); 1.4353 M☉ with m_p', () => {
  close(chandrasekharMass(1) / Msun, 5.83, 1e-3, 'textbook 5.83/μe²');
  close(chandrasekharMass(2) / Msun, 1.4563, 1e-4);
  close(chandrasekharMass(2, { mUnit: mp }) / Msun, 1.4353, 1e-4);
  close(chandrasekharMass(2) / chandrasekharMass(2.15), (2.15 / 2) ** 2, 1e-14);
});

test('white dwarf (exact degenerate EOS): xc → ∞ tends to M_Ch with R → 0; xc ≪ 1 is the n = 3/2 polytrope', () => {
  const Mch = chandrasekharMass(2);
  const m = [10, 100, 1e3, 1e4].map(xc => whiteDwarf({ xc }).mass);
  for (let i = 1; i < m.length; i++) assert.ok(m[i] > m[i - 1] && m[i] < Mch, 'mass rises monotonically below M_Ch');
  close(m[3], Mch, 1e-6);
  assert.ok(whiteDwarf({ xc: 1e4 }).radius < whiteDwarf({ xc: 1e3 }).radius / 9);
  // non-relativistic limit: P = K ρ^{5/3}, K = (8A/5)/B^{5/3}
  const w = whiteDwarf({ xc: 1e-3 });
  const A = Math.PI * 9.1093837015e-31 ** 4 * 299792458 ** 5 / (3 * 6.62607015e-34 ** 3);
  const B = 2 * 1.66053906660e-27 * 8 * Math.PI / 3 * (9.1093837015e-31 * 299792458 / 6.62607015e-34) ** 3;
  const K = 8 * A / 5 / B ** (5 / 3), le = laneEmden(1.5), alpha = Math.sqrt(2.5 * K * w.rhoC ** (-1 / 3) / (4 * Math.PI * G));
  close(w.radius, alpha * le.xi1, 2e-6); close(w.mass, 4 * Math.PI * alpha ** 3 * w.rhoC * le.omega, 2e-6);
  close(whiteDwarf({ rhoC: w.rhoC }).mass, w.mass, 1e-9);
});

test('white dwarf vs Sirius B (Bond et al., ApJ 840, 70 (2017): M = 1.018 M☉, R = 0.0081 R☉), labelled zero-T C/O model', () => {
  const w = whiteDwarfOfMass(1.018 * Msun);
  close(w.mass, 1.018 * Msun, 1e-8);
  close(w.radius / R_SUN, 0.0081, 0.05);
  assert.throws(() => whiteDwarfOfMass(1.5 * Msun), /M_Ch/);
});

test('AME2020 table is consistent with NUCLIDE_MASS atomic masses (H-2, C-12, O-16 to 0.01 keV); Ni-62 > Fe-58 > Fe-56', () => {
  close(bindingFromMass(1, 2, 2.01410177812) / 2 * 1000, AME2020_BA.H2.BA, 1e-7);
  close(bindingFromMass(6, 12, 12) / 12 * 1000, AME2020_BA.C12.BA, 1e-7);
  close(bindingFromMass(8, 16, 15.99491461957) / 16 * 1000, AME2020_BA.O16.BA, 1e-7);
  const { Ni62, Fe58, Fe56 } = AME2020_BA;
  assert.ok(Ni62.BA > Fe58.BA && Fe58.BA > Fe56.BA);
  for (const v of Object.values(AME2020_BA)) assert.ok(v.BA <= Ni62.BA);
});

test('Bethe–Weizsäcker (Rohlf coefficients) vs AME2020: within 1 % for A ≥ 40, peak of B/A at A ≈ 56–62, light nuclei poor', () => {
  for (const k of ['Ca40', 'Fe56', 'Fe58', 'Ni62', 'Pb208', 'U238']) {
    const { Z, A, BA } = AME2020_BA[k];
    close(semfBinding(Z, A) / A * 1000, BA, 0.01, k);
  }
  assert.ok(rel(semfBinding(2, 4) / 4 * 1000, AME2020_BA.He4.BA) > 0.1, 'He-4 is badly described (disclose)');
  let best = { b: 0 };
  for (let A = 10; A <= 250; A++) for (let Z = 1; Z < A; Z++) { const b = semfBinding(Z, A) / A; if (b > best.b) best = { b, Z, A }; }
  assert.ok(best.A >= 54 && best.A <= 64, `SEMF peak at A = ${best.A}`);
  // valley of stability: within 1.5 of the stable isobar (pairing and shell effects not modelled)
  for (const [A, Zs] of [[56, 26], [62, 28], [127, 53], [208, 82], [238, 92]]) assert.ok(Math.abs(mostStableZ(A) - Zs) < 1.5, `A=${A}: ${mostStableZ(A)}`);
  // pairing term: +aP/√A even–even, −aP/√A odd–odd, 0 odd A (Rohlf aP = 11.18 MeV)
  const k = { aV: 15.75, aS: 17.8, aC: 0.711, aA: 23.7, aP: 11.18 }, smooth = (Z, A) => semfBinding(Z, A, { ...k, aP: 0 });
  close(semfBinding(26, 56) - smooth(26, 56), 11.18 / Math.sqrt(56), 1e-12);
  close(semfBinding(29, 58) - smooth(29, 58), -11.18 / Math.sqrt(58), 1e-12);
  assert.equal(semfBinding(6, 13), smooth(6, 13));
  assert.throws(() => semfBinding(7, 6), /Z/);
});

test('Saha: x = 0.5 at the half-ionisation temperature; x → 1 hot, → 0 cold; Carroll & Ostlie (P_e = 20 Pa → 50 % at ≈ 9600 K)', () => {
  for (const n of [1e10, 1e20, 1e26]) {
    const T = halfIonisationTemperature(n);
    assert.ok(Math.abs(sahaHydrogen(T, n).x - 0.5) < 1e-10);
  }
  assert.ok(halfIonisationTemperature(1e26) > halfIonisationTemperature(1e20), 'denser gas ionises at higher T');
  assert.ok(sahaHydrogen(3e4, 1e20).x > 0.999 && sahaHydrogen(3000, 1e20).x < 1e-6);
  const s = sahaHydrogen(9000, 1e20); close(s.x * s.x / (1 - s.x), sahaRatio(9000, 1) / 1e20, 1e-12);
  const T50 = brent(T => sahaRatio(T, 20 / (kB * T)) - 1, 5000, 20000);
  close(T50, 9600, 0.02, 'C&O §8.1 example');
});

test('energy generation power laws (APPROXIMATE): exponents 4 and 19.9, pp/CNO crossover near 1.8e7 K for X = 0.7, X_CNO = 0.01', () => {
  close(epsPP(1e5, 0.7, 2e7) / epsPP(1e5, 0.7, 1e7), 16, 1e-14);
  close(epsCNO(1e5, 0.7, 0.01, 2e7) / epsCNO(1e5, 0.7, 0.01, 1e7), 2 ** 19.9, 1e-12);
  const Tx = brent(T => Math.log(epsPP(1e5, 0.7, T) / epsCNO(1e5, 0.7, 0.01, T)), 1e7, 3e7);
  assert.ok(Tx > 1.6e7 && Tx < 2e7, `${Tx}`);
  const e = epsPP(SSM_CENTRE.rhoC, 0.34, SSM_CENTRE.Tc);
  assert.ok(e > 3e-4 && e < 5e-3, 'solar centre ~1e-3 W/kg (order of magnitude)');
});

test('Stefan–Boltzmann: IAU 2015 B3 nominal L☉ and R☉ give T_eff = 5772 K; inverse functions agree', () => {
  assert.ok(Math.abs(effectiveTemperature(L_SUN, R_SUN) - TEFF_SUN) < 0.5);
  close(stefanBoltzmannLuminosity(R_SUN, effectiveTemperature(L_SUN, R_SUN)), L_SUN, 1e-14);
  assert.throws(() => effectiveTemperature(-1, R_SUN), /> 0/);
});

test('Eddington luminosity 1.257e31 W per M☉ (pure H); σ_T = 6.6524587321e-29 m² (CODATA 2018)', () => {
  close(SIGMA_T, 6.6524587321e-29, 1e-9);
  close(eddingtonLuminosity(Msun), 1.2571e31, 1e-3);
  close(eddingtonLuminosity(Msun, { X: 0.7 }), eddingtonLuminosity(Msun) * 2 / 1.7, 1e-14);
  close(eddingtonLuminosity(Msun) / L_SUN, 3.28e4, 2e-3);
});

test('timescales: Kelvin–Helmholtz GM²/RL ≈ 31.4 Myr and nuclear 0.1·0.007·Mc²/L ≈ 10.4 Gyr for the Sun', () => {
  close(kelvinHelmholtzTime(Msun, R_SUN, L_SUN) / YEAR / 1e6, 31.4, 2e-3);
  close(kelvinHelmholtzTime(Msun, R_SUN, L_SUN, { q: 1.5 }), 1.5 * kelvinHelmholtzTime(Msun, R_SUN, L_SUN), 1e-15);
  close(nuclearTime(Msun, L_SUN) / YEAR / 1e9, 10.36, 1e-3);
});

test('main-sequence scalings (APPROXIMATE): Sun at 1, continuous-ish and monotone, massive stars short-lived', () => {
  assert.equal(mainSequenceLuminosity(1), 1); assert.equal(mainSequenceRadius(1), 1);
  close(mainSequenceLifetime(1), 10, 1e-15);
  let prev = 0;
  for (let m = 0.1; m < 100; m *= 1.1) { const L = mainSequenceLuminosity(m); assert.ok(L > prev); prev = L; }
  assert.ok(mainSequenceLifetime(10) < 0.05 && mainSequenceLifetime(0.5) > 50);
  assert.throws(() => mainSequenceLuminosity(0.01), /M/);
});
