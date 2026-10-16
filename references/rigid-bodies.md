# Rigid bodies, rotations and tops

Status: generic tools (`core/quaternion.mjs` for orientations, `core/rigid.mjs` for inertia, torque-free motion and the heavy top). No ready-made recipe. Read [method selection](method-selection.md) for the general integrator choice and [math toolkit](math-toolkit.md) for the shared numerics.

A spinning body is two coupled problems: the **angular momentum in the body frame** (Euler's equations) and the **orientation** that carries it into space. Both have exact invariants, and an instrument that lets them drift (a top that slowly speeds up, a book that stops flipping) is visibly wrong. Use the structure-preserving stepper, not a generic RK on Euler angles: Euler angles are singular at θ = 0 (gimbal lock), and RK drifts in energy and |L|.

## Candidate models (pick by the question)

| Question | Model | Observables and benchmark |
|---|---|---|
| Why does a book / phone / T-handle flip when spun about its middle axis? (Dzhanibekov effect, tennis-racket theorem) | Torque-free asymmetric body, I1 < I2 < I3, L0 near axis 2 | Flip interval = half the elliptic period, 2K(m)/λ (`torqueFreeExact`); rotation about axes 1 and 3 stays stable; E, \|L\| and spatial L constant |
| Free precession, wobbling plates, Chandler wobble | Torque-free symmetric body | Body-frame precession of L about the symmetry axis at (I3 − I1)/I1 · ω3 exactly; the axis cones about L at \|L\|/I1 |
| Why does a spinning top not fall? Gyroscope precession | Heavy symmetric (Lagrange) top about a fixed pivot | Steady precession Ω from I1 cosθ Ω² − I3ω3 Ω + mgl = 0, → mgl/(I3ω3) for a fast top; nutation band [θmin, θmax] and period 4K(m)/√(β(u3 − u1)) (`lagrangeTopNutation`); Casimirs L·Γ, \|Γ\| = 1, L3 |
| Moments of inertia of a shape or an assembly | `inertiaBox`, `inertiaCylinder`, `inertiaSphere`, `inertiaCone`, `compositeBody` (parallel-axis theorem), `principalAxes` | Symmetry limits (cube, isotropic cylinder h² = 3r²), slice integration, assembly from halves |
| Animate an orientation smoothly | Quaternion keyframes with `qSlerp`; body-rate integration with `qIntegrate` | Constant angular speed between keys; exact for constant ω |

Rolling contact (rolling coin, rattleback, tippe top) is non-holonomic and not provided: say so if asked, and do not fake it with the free-body stepper.

## Method

- **Orientation:** unit quaternion q = [w, x, y, z] mapping body to space; never integrate Euler angles. `qIntegrate(q, ω, dt)` is exact for constant body ω. `qToMatrix` for drawing; `qRotate` to move body points.
- **Torque-free motion:** `freeRigidBody` (Lie–Poisson splitting: each piece of H = Σ L_i²/2I_i is an exact rotation). Symplectic and reversible; |L| and spatial L conserved to round-off; energy error bounded, O(h²) (`'strang'`) or O(h⁴) (`'yoshida4'`). Exact for symmetric bodies.
- **Heavy top:** `heavyTop` with moments about the pivot and mgl; set the initial state from Euler angles and rates with `lagrangeTopInitial`. Gravity is −z.
- **Exact reference:** `torqueFreeExact` (Jacobi elliptic functions, any initial L) and `jacobiElliptic`/`ellipticF` for your own closed forms.
- Moments of inertia about the **pivot** for the heavy top: `parallelAxis(Icm, m, d)`.

## Verification

- Integrator against `torqueFreeExact` over several periods, with the observed order (`refinementOrder`): 2 for `'strang'`, 4 for `'yoshida4'`.
- Conserved to round-off: |L| and `spatialMomentum()` (free body); |Γ|, L·Γ and, for a symmetric top, L3 (heavy top). Energy error bounded without drift over 10⁴–10⁵ steps.
- Time reversal: set `s.h = −h` and retrace to round-off.
- Dzhanibekov flip interval against 2K(m)/λ (`period / 2`).
- Lagrange top: integrated θ extremes and nutation period against `lagrangeTopNutation`; steady precession at the slow root keeps θ fixed; a dropped fast top precesses on average at mgl/(I3ω3).
- Inertia: symmetry limits and an independent route (slice integration or assembly from parts).

## Representation

- Draw the body with its principal axes and the angular-momentum vector L (fixed in space for a free body) and, optionally, ω (it moves). Make the conservation of L visible: it is the reason the motion looks surprising.
- Show the body-frame trajectory of L on the energy ellipsoid ∩ momentum sphere (polhodes): the separatrix through the intermediate axis explains the flip.
- For tops, draw the trace of the symmetry axis tip on a sphere: loops, cusps or waves depending on the initial φ̇, and label the steady-precession rate.
- State the frames: "body frame" plots rotate with the object; "space frame" plots do not.

**Visible scope:** "Ideal rigid body: no friction at the pivot, no air drag, no deformation; rolling contact not modelled." Say which body (moments, mgl) the numbers refer to.

## API

`core/quaternion.mjs` (q = [w, x, y, z], active rotation, body → space):

- `qmul(a, b, out?)` → a ⊗ b. `qconj(q)`, `qnorm(q)`, `qnormalize(q, out?)` (throws for zero), `qinverse(q)`.
- `qFromAxisAngle(axis, angle)` → unit q (axis any non-zero length). `qToAxisAngle(q)` → `{ axis, angle ∈ [0, π] }`.
- `qToMatrix(q)` → 3×3 rows, orthogonal with det +1 (normalises q). `qFromMatrix(R, { tolerance })` → q with w ≥ 0 (Shepperd; throws if R is not a proper rotation).
- `qRotate(q, v, out?)` → q v q* for unit q. `qSlerp(a, b, t)` → shorter-arc interpolation, constant angular speed.
- `qExp(q)`, `qLog(q)` → principal exp/log; `qExp(qLog(q)) = q`.
- `qIntegrate(q, ω, dt, { frame: 'body' | 'space' })` → q ⊗ exp(½ω dt) or exp(½ω dt) ⊗ q; exact for constant ω, Lie–Euler otherwise (use the midpoint ω for second order).

`core/rigid.mjs` (consistent units; body frame = principal axes):

- `inertiaBox(m, [a, b, c])`, `inertiaSphere(m, r, { shell })`, `inertiaCylinder(m, r, h, { innerRadius })` → 3×3 about the centre of mass. `inertiaCone(m, r, h)` → `{ inertia, comFromBase: h/4 }`.
- `parallelAxis(Icm, m, d)` → inertia about the point at −d from the centre of mass. `compositeBody([{ mass, inertia, position, orientation? }])` → `{ mass, com, inertia }` (negative mass for cavities).
- `principalAxes(I)` → `{ moments (ascending), axes (right-handed), rotation, quaternion }`; throws unless positive definite.
- `jacobiElliptic(u, m)` → `{ sn, cn, dn, am }`, 0 ≤ m ≤ 1. `ellipticF(φ, m)` → F(φ | m) for any real φ.
- `torqueFreeExact({ inertia: [I1 ≤ I2 ≤ I3], L0 })` → `{ at(t), period, m, lambda, axis, energy, L2 }`; throws on the separatrix.
- `freeRigidBody({ inertia, L0, q0, h, method: 'strang' | 'yoshida4' })` and `heavyTop({ inertia, mgl, chi, L0, q0, h, method })` → stepper `{ t, q, L, h, step(), advanceTo(t, maxSteps), omega(), energy(), spatialMomentum(), rotation(), gamma() }`.
- `lagrangeTopInitial({ I1, I3, theta, phi, psi, thetaDot, phiDot, psiDot })` → `{ q, L, omega }` (z-x-z Euler angles).
- `lagrangeTopNutation({ I1, I3, mgl, theta, thetaDot, phiDot, psiDot })` → `{ thetaMin, thetaMax, period, roots, pPhi, pPsi, energy }`; needs mgl > 0, 0 < θ < π.
- `steadyPrecession({ I1, I3, mgl, omega3, theta })` → `{ slow, fast, approx: mgl/(I3ω3) }`; throws if the spin is too slow.
