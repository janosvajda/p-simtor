import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  createParticles, gravity, coulomb, uniformField, centralBody, drag, springs, externalPotential,
} from '../core/particles.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);

test('two-body gravity: Kepler period, energy and momentum conserved (symplectic path)', () => {
  // equal masses on a circular orbit of separation 1: ω² = G(m1+m2)/r³ = 2
  const v = Math.sqrt(2) / 2;
  const sys = createParticles({ bodies: [{ m: 1, x: [-0.5, 0, 0], v: [0, -v, 0] }, { m: 1, x: [0.5, 0, 0], v: [0, v, 0] }], forces: [gravity()], h: 1e-3 });
  assert.equal(sys.method, 'yoshida4');
  const P = 2 * Math.PI / Math.sqrt(2);
  sys.advanceTo(P);
  const x = sys.sampleAt(P);
  close(x[3], 0.5, 1e-9, 'returns after one period');
  assert.ok(sys.energyBalanceError() < 1e-12 && sys.momentumDrift() < 1e-13);
});

test('Coulomb repulsion: energy conserved through a head-on approach', () => {
  const sys = createParticles({
    bodies: [{ m: 1, charge: 1, x: [-5, 0, 0], v: [1, 0, 0] }, { m: 1, charge: 1, x: [5, 0, 0], v: [-1, 0, 0] }],
    forces: [coulomb({ k: 1 })], h: 1e-3,
  });
  sys.advanceTo(10);
  // closest approach where kinetic energy 1 + U(10) = 1/r_min: r_min = 1/(1 + 0.1)
  assert.ok(sys.energyBalanceError() < 1e-10);
  assert.ok(sys.x[3] - sys.x[0] > 0, 'they bounce back');
});

test('uniform B: gyration radius m v / (qB) and period 2πm/(qB); magnetic force does no work', () => {
  const sys = createParticles({ bodies: [{ m: 2, charge: 1, x: [0, 0, 0], v: [3, 0, 0.5] }], forces: [uniformField({ B: [0, 0, 1.5] })], h: 0.01, rtol: 1e-12 });
  assert.equal(sys.method, 'dopri5');
  const r = 2 * 3 / 1.5, T = 2 * Math.PI * 2 / 1.5;
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  while (sys.t < T) { sys.step(); xmin = Math.min(xmin, sys.x[0]); xmax = Math.max(xmax, sys.x[0]); ymin = Math.min(ymin, sys.x[1]); ymax = Math.max(ymax, sys.x[1]); }
  close(Math.max(xmax - xmin, ymax - ymin) / 2, r, 2e-3, 'radius');
  const back = sys.sampleAt(T);
  close(Math.hypot(back[0], back[1]), 0, 1e-8, 'closed orbit after one period');
  close(back[2], 0.5 * T, 1e-8, 'drift along B');
  assert.ok(sys.energyBalanceError() < 1e-10);
  assert.equal(sys.momentumDrift(), null, 'external field: momentum not conserved, so not reported');
});

test('E × B drift velocity = E/B', () => {
  const sys = createParticles({ bodies: [{ m: 1, charge: 1, x: [0, 0, 0], v: [0, 0, 0] }], forces: [uniformField({ E: [0, 0.2, 0], B: [0, 0, 1] })], h: 0.01, rtol: 1e-11 });
  const T = 2 * Math.PI;                                   // starting at rest, the guiding centre drifts at E/B
  sys.advanceTo(10 * T);
  const x = sys.sampleAt(10 * T);
  close(x[0] / (10 * T), 0.2, 1e-8, 'drift speed along E × B (+x)');
});

test('drag: work–energy balance holds and velocity decays as e^(−γt)', () => {
  const sys = createParticles({ bodies: [{ m: 3, x: [0, 0, 0], v: [2, 0, 0] }], forces: [drag({ gamma: 0.4 })], h: 0.01 });
  sys.advanceTo(5);
  const v = sys.v[0];
  close(v, 2 * Math.exp(-0.4 * sys.t), 1e-9);
  assert.ok(sys.work < 0 && sys.energyBalanceError() < 1e-9, `${sys.energyBalanceError()}`);
});

test('springs: two masses oscillate at ω = √(2k/m)', () => {
  const k = 5, m = 2, w = Math.sqrt(2 * k / m);
  const sys = createParticles({ bodies: [{ m, x: [-1.1, 0, 0], v: [0, 0, 0] }, { m, x: [1.1, 0, 0], v: [0, 0, 0] }], forces: [springs([[0, 1, k, 2]])], h: 1e-3 });
  sys.advanceTo(2 * Math.PI / w);
  close(sys.sampleAt(2 * Math.PI / w)[3], 1.1, 1e-8);
  assert.ok(sys.energyBalanceError() < 1e-12);
});

test('external potential: harmonic trap behaves like the analytic oscillator', () => {
  const trap = externalPotential({ U: (p, b) => 0.5 * b.m * (p[0] ** 2 + p[1] ** 2 + p[2] ** 2), gradU: (p, b, g) => { g[0] = b.m * p[0]; g[1] = b.m * p[1]; g[2] = b.m * p[2]; } });
  const sys = createParticles({ bodies: [{ m: 1.5, x: [1, 0, 0], v: [0, 2, 0] }], forces: [trap], h: 1e-3 });
  sys.advanceTo(1.3);
  const x = sys.sampleAt(1.3);
  close(x[0], Math.cos(1.3), 1e-9); close(x[1], 2 * Math.sin(1.3), 1e-9);
});

test('central body with relativity: perihelion advance 6πGM/(c²a(1−e²)) and 1PN energy conserved', () => {
  // GM = 1 units; c chosen so the effect is resolvable: GM/(c²p) = 1e-4
  const GM = 1, a = 1, e = 0.3, p = a * (1 - e * e), c = Math.sqrt(GM / (1e-4 * p));
  const rp = a * (1 - e), vp = Math.sqrt(GM * (1 + e) / rp);
  // angle of the next periapsis (minimum of r: r·v changes sign − → +), located by bisection
  // on the dense output, so the result does not depend on where steps happen to fall
  const periapsisAngle = relativity => {
    const sys = createParticles({ bodies: [{ m: 1, x: [rp, 0, 0], v: [0, vp, 0] }], forces: [centralBody({ GM, relativity, c })], h: 1e-3, rtol: 1e-12, atol: 1e-14 });
    const rv = t => { const { x, v } = sys.sampleStateAt(t); return x[0] * v[0] + x[1] * v[1]; };
    let tPrev = 0, prev = 0;
    for (;;) {
      tPrev = sys.t; sys.step();
      const now = sys.x[0] * sys.v[0] + sys.x[1] * sys.v[1];
      if (sys.t > 3 && prev < 0 && now >= 0) {
        let lo = tPrev, hi = sys.t;
        for (let i = 0; i < 60; i++) { const mid = (lo + hi) / 2; rv(mid) < 0 ? (lo = mid) : (hi = mid); }
        const { x } = sys.sampleStateAt((lo + hi) / 2);
        return { angle: Math.atan2(x[1], x[0]), balance: sys.energyBalanceError() };
      }
      prev = now;
    }
  };
  const gr = periapsisAngle(true), newton = periapsisAngle(false);
  const measured = gr.angle - newton.angle;                 // differential: shared integration error cancels
  const expected = 6 * Math.PI * GM / (c * c * p);
  close(measured, expected, 0.01 * expected, 'perihelion advance (first order; next order ~ GM/c²p)');
  assert.ok(gr.balance < 1e-8, `1PN energy ${gr.balance}`);
});

test('bad input is rejected with clear messages', () => {
  assert.throws(() => createParticles({ bodies: [{ m: 0, x: [0, 0], v: [0, 0] }], forces: [], h: 0.1 }), /positive mass/);
  assert.throws(() => createParticles({ bodies: [{ m: 1, x: [NaN, 0], v: [0, 0] }], forces: [], h: 0.1 }), /finite/);
  assert.throws(() => createParticles({ bodies: [{ m: 1, x: [0, 0], v: [0, 0] }], forces: [drag({ gamma: 1 })], method: 'yoshida4', h: 0.1 }), /velocity-dependent/);
});
