import { test } from 'node:test';
import assert from 'node:assert/strict';
import { dopri5, brent } from '../core/numerics.mjs';
import {
  cosmology, PLANCK18, PLANCK18_DERIVED, omegaGammaH2, E, hubbleParameter, hubbleTime, hubbleDistance, MPC, GYR,
  HUBBLE_TIME_UNIT, age, lookbackTime, comovingDistance, transverseComovingDistance, angularDiameterDistance,
  luminosityDistance, distanceModulus, particleHorizon, zMatterRadiation, zMatterLambda, decelerationParameter,
  zAcceleration, cmbTemperature, growthFactor, growthRate, growthFactorCPT, soundHorizon, soundHorizonDragFit,
  friedmannTrajectory,
} from '../core/cosmology.mjs';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const close = (a, b, tol, msg = '') => assert.ok(rel(a, b) <= tol, `${msg} ${a} vs ${b} (rel tol ${tol})`);

const EdS = cosmology({ H0: 70, Om: 1, OL: 0, Or: 0 });
const milne = cosmology({ H0: 70, Om: 0, OL: 0, Or: 0 });
const deSitter = cosmology({ H0: 70, Om: 0, OL: 1, Or: 0 });
const radOnly = cosmology({ H0: 70, Om: 0, OL: 0, Or: 1 });

test('units: Mpc (IAU 2015 B2), Gyr, Hubble time unit 977.79 Gyr; Ωγh² = 2.473e-5 at 2.7255 K (PDG)', () => {
  close(MPC, 3.08567758149137e22, 1e-14);
  close(GYR, 3.15576e16, 1e-15);
  close(HUBBLE_TIME_UNIT, 977.792221, 1e-8);
  close(omegaGammaH2(2.7255), 2.473e-5, 2e-4);
  close(hubbleDistance(EdS), 299792.458 / 70, 1e-15);
});

test('model construction: flat closes the budget, curvature derived, invalid input throws', () => {
  assert.equal(PLANCK18.Ok, 0);
  close(PLANCK18.Om + PLANCK18.Or + PLANCK18.OL, 1, 1e-15);
  close(PLANCK18.OL, 0.6889, 2e-4, 'Planck Table 2 ΩΛ (rounded)');
  assert.equal(cosmology({ H0: 70, Om: 2, OL: 0, Or: 0 }).Ok, -1);
  assert.throws(() => cosmology({ H0: -1, Om: 0.3, flat: true }), /H0/);
  assert.throws(() => cosmology({ H0: 70, Om: 0.3 }), /OL/);
  assert.throws(() => E(-1.5, PLANCK18), /z/);
  assert.equal(E(0, PLANCK18), 1);
  close(hubbleParameter(1, EdS), 70 * 2 ** 1.5, 1e-15);
});

test('exact ages: EdS 2/(3H0), Milne 1/H0, radiation-only 1/(2H0), flat ΛCDM closed form; de Sitter throws', () => {
  close(age(EdS), 2 / 3 * hubbleTime(EdS), 1e-13);
  close(age(EdS, 3), 2 / 3 * hubbleTime(EdS) / 8, 1e-13);
  close(age(milne), hubbleTime(milne), 1e-13);
  close(age(radOnly), hubbleTime(radOnly) / 2, 1e-13);
  const L = cosmology({ H0: 70, Om: 0.3, flat: true, Or: 0 });
  close(age(L), 2 / (3 * Math.sqrt(0.7)) * Math.asinh(Math.sqrt(0.7 / 0.3)) * hubbleTime(L), 1e-13);
  assert.throws(() => age(deSitter), /no Big Bang/);
});

test('closed matter-only universe: age on the cycloid H0 t = Ωm/(2(Ωm−1)^{3/2})(η − sin η)', () => {
  const Om = 3, p = cosmology({ H0: 70, Om, OL: 0, Or: 0 });
  const A = Om / (2 * (Om - 1)), B = Om / (2 * (Om - 1) ** 1.5);
  for (const z of [0, 1, 5]) {
    const eta = Math.acos(1 - 1 / ((1 + z) * A));
    close(age(p, z), B * (eta - Math.sin(eta)) * hubbleTime(p), 1e-12, `z=${z}`);
  }
});

test('lookback time: Milne z/(1+z), de Sitter ln(1+z), EdS (2/3)(1 − (1+z)^{-3/2}); t0 − t(z) consistency', () => {
  for (const z of [0.1, 1, 10]) {
    close(lookbackTime(z, milne), hubbleTime(milne) * z / (1 + z), 1e-13);
    close(lookbackTime(z, deSitter), hubbleTime(deSitter) * Math.log1p(z), 1e-12);
    close(lookbackTime(z, EdS), 2 / 3 * hubbleTime(EdS) * (1 - (1 + z) ** -1.5), 1e-13);
    close(lookbackTime(z, PLANCK18), age(PLANCK18) - age(PLANCK18, z), 1e-12);
  }
});

test('distances, exact: EdS 2D_H(1 − 1/√(1+z)), de Sitter D_H z, Milne D_L = D_H z(1 + z/2) (open curvature)', () => {
  for (const z of [0.5, 2, 10]) {
    close(comovingDistance(z, EdS), 2 * hubbleDistance(EdS) * (1 - 1 / Math.sqrt(1 + z)), 1e-13);
    close(comovingDistance(z, deSitter), hubbleDistance(deSitter) * z, 1e-13);
    close(luminosityDistance(z, milne), hubbleDistance(milne) * z * (1 + z / 2), 1e-12);
    close(comovingDistance(z, milne), hubbleDistance(milne) * Math.log1p(z), 1e-12);
  }
});

test('distances with curvature: Mattig formula for ΩΛ = 0, open and closed (Hogg 1999 astro-ph/9905116, eq. 16 note)', () => {
  for (const Om of [0.3, 2]) {
    const p = cosmology({ H0: 70, Om, OL: 0, Or: 0 });
    for (const z of [0.5, 1, 3]) {
      const dm = hubbleDistance(p) * 2 * (2 - Om * (1 - z) - (2 - Om) * Math.sqrt(1 + Om * z)) / (Om * Om * (1 + z));
      close(transverseComovingDistance(z, p), dm, 1e-12, `Om=${Om} z=${z}`);
    }
  }
});

test('Etherington reciprocity D_L = (1+z)² D_A; distance modulus 5 log₁₀(D_L/10 pc)', () => {
  const z = 1.3;
  close(luminosityDistance(z, PLANCK18), (1 + z) ** 2 * angularDiameterDistance(z, PLANCK18), 1e-14);
  close(distanceModulus(z, PLANCK18), 5 * Math.log10(luminosityDistance(z, PLANCK18) * 1e6 / 10), 1e-14);
  assert.throws(() => distanceModulus(0, PLANCK18), /z/);
});

test('astropy FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0) documented values (astropy docs, cosmology): D_C and age', () => {
  const p = cosmology({ H0: 70, Om: 0.3, flat: true, Or: 0 });
  close(comovingDistance(0.5, p), 1888.62539593, 1e-9);
  close(comovingDistance(1.0, p), 3303.82880589, 1e-9);
  close(comovingDistance(1.5, p), 4363.85576624, 1e-9);
  close(age(p), 13.466983947, 1e-9);
});

test('Planck 2018 (arXiv:1807.06209 Table 2): age 13.787 ± 0.020 Gyr and z_eq 3387 ± 21 within 1σ', () => {
  const { age: A, zEq } = PLANCK18_DERIVED;
  assert.ok(Math.abs(age(PLANCK18) - A.value) < A.error, `${age(PLANCK18)}`);
  // massive ν counted in Ωm (Planck excludes it at z_eq): +16 shift, inside 1σ
  assert.ok(Math.abs(zMatterRadiation(PLANCK18) - zEq.value) < zEq.error, `${zMatterRadiation(PLANCK18)}`);
});

test('horizons: EdS 2c/H0, radiation-only c/H0, none in Milne (infinite)', () => {
  close(particleHorizon(EdS), 2 * hubbleDistance(EdS), 1e-13);
  close(particleHorizon(radOnly), hubbleDistance(radOnly), 1e-13);
  assert.equal(particleHorizon(milne), Infinity);
  assert.ok(particleHorizon(PLANCK18) > 14000 && particleHorizon(PLANCK18) < 14300, 'Planck: ≈ 14.2 Gpc comoving');
});

test('equality and acceleration redshifts: exact without radiation, q(z_acc) = 0 with radiation', () => {
  const p = cosmology({ H0: 70, Om: 0.3, flat: true, Or: 0 });
  close(1 + zAcceleration(p), Math.cbrt(2 * 0.7 / 0.3), 1e-12);
  close(1 + zMatterLambda(p), Math.cbrt(0.7 / 0.3), 1e-15);
  close(decelerationParameter(0, p), 0.15 - 0.7, 1e-14);
  assert.ok(Math.abs(decelerationParameter(zAcceleration(PLANCK18), PLANCK18)) < 1e-12);
  assert.throws(() => zAcceleration(EdS), /never accelerates/);
  close(zMatterRadiation(cosmology({ H0: 70, Om: 1, OL: 0, Or: 1e-3 })), 999, 1e-12);
});

test('CMB temperature scales as T0(1+z): 2973 K at Planck z* = 1089.80', () => {
  close(cmbTemperature(1089.80), 2.7255 * 1090.80, 1e-15);
  assert.throws(() => cmbTemperature(-2), /z/);
});

test('growth factor: EdS D = a; open ΩΛ = 0 closed form (Heath 1977 / Peebles 1980 §11)', () => {
  for (const a of [0.01, 0.3, 1]) close(growthFactor(a, EdS), a, 1e-13);
  close(growthRate(0.5, EdS), 1, 1e-12);
  const Om = 0.3, p = cosmology({ H0: 70, Om, OL: 0, Or: 0 });
  const exact = a => { const x = (1 / Om - 1) * a; return 1 + 3 / x + 3 * Math.sqrt(1 + x) / x ** 1.5 * Math.log(Math.sqrt(1 + x) - Math.sqrt(x)); };
  close(growthFactor(1, p) / growthFactor(0.25, p), exact(1) / exact(0.25), 1e-9);
});

test('growth factor for ΛCDM vs direct ODE D″ + (2 + dlnE/dlna)D′ − (3/2)Ωm(a)D = 0; f(a) vs numerical derivative', () => {
  const p = cosmology({ H0: 67.66, Om: 0.3111, flat: true, Or: 0 });
  const e2 = a => p.Om / a ** 3 + p.OL;
  const f = (x, y) => { const a = Math.exp(x), dl = -1.5 * p.Om / a ** 3 / e2(a); return [y[1], -(2 + dl) * y[1] + 1.5 * p.Om / a ** 3 / e2(a) * y[0]]; };
  const a0 = 1e-3, { y } = dopri5(f, Math.log(a0), [a0, a0], 0, { rtol: 1e-12, atol: 1e-16 });
  close(growthFactor(1, p), y[0], 1e-7);
  const h = 1e-4, num = (Math.log(growthFactor(Math.exp(h), p)) - Math.log(growthFactor(Math.exp(-h), p))) / (2 * h);
  close(growthRate(1, p), num, 1e-7);
  close(growthRate(1, p), Math.pow(p.Om, 0.55), 0.01, 'Linder γ = 0.55 approximation');
  close(growthFactorCPT(1, p), growthFactor(1, p), 0.01, 'Carroll–Press–Turner fit (~1 %)');
  assert.throws(() => growthFactor(1, milne), /Om > 0/);
});

test('sound horizon: integral at Planck z_drag and z* within 0.2 %; Aubourg et al. 2015 fit within 0.05 %', () => {
  const { zDrag, rDrag, zStar, rStar, Och2 } = PLANCK18_DERIVED;
  close(soundHorizon(zDrag.value, PLANCK18), rDrag.value, 2e-3);
  close(soundHorizon(zStar.value, PLANCK18), rStar.value, 2e-3);
  close(soundHorizonDragFit({ obh2: 0.02242, ocbh2: 0.02242 + Och2.value }), rDrag.value, 5e-4);
  assert.throws(() => soundHorizon(1000, EdS), /Obh2/);
});

test('Friedmann ODE: closed cycloid a = 1 − cos η through recollapse to the crunch at 2π/H0 (Ωm = 2)', () => {
  const p = cosmology({ H0: 70, Om: 2, OL: 0, Or: 0 }), tH = hubbleTime(p);
  const err = rtol => {
    const tr = friedmannTrajectory(p, { tMax: 100, n: 40, rtol });
    let m = 0;
    for (let i = 0; i < tr.a.length - 2; i++) {
      const T = tr.t[i] / tH, eta = brent(e => e - Math.sin(e) - T, 0, 2 * Math.PI);
      m = Math.max(m, Math.abs(tr.a[i] - (1 - Math.cos(eta))));
    }
    return { m, crunch: tr.crunch };
  };
  const coarse = err(1e-8), fine = err(1e-12);
  assert.ok(fine.m < 3e-7, `${fine.m}`);
  assert.ok(coarse.m / fine.m > 1000, 'converges with tolerance');
  close(fine.crunch / tH, 2 * Math.PI, 1e-5);
});

test('Friedmann ODE: EdS a = (t/t0)^{2/3}, Milne a = H0 t, de Sitter a = e^{H0 t} (from today), Planck a(t0) = 1', () => {
  const t0 = age(EdS), tr = friedmannTrajectory(EdS, { tMax: 2 * t0, n: 11, rtol: 1e-12 });
  tr.t.forEach((t, i) => close(tr.a[i], (t / t0) ** (2 / 3), 1e-8));
  const tm = friedmannTrajectory(milne, { tMax: 30, n: 7, rtol: 1e-12 });
  tm.t.forEach((t, i) => close(tm.a[i], t / hubbleTime(milne), 1e-8));
  const td = friedmannTrajectory(deSitter, { tMax: 20, n: 6, origin: 'today', rtol: 1e-12 });
  td.t.forEach((t, i) => close(td.a[i], Math.exp(t / hubbleTime(deSitter)), 1e-8));
  const tp = friedmannTrajectory(PLANCK18, { tMax: age(PLANCK18), n: 3, rtol: 1e-12 });
  close(tp.a[2], 1, 1e-7); close(tp.H[2], 67.66, 1e-6);
  assert.equal(tp.crunch, null);
  assert.throws(() => friedmannTrajectory(PLANCK18, {}), /tMax/);
});
