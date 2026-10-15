// Spin ½ and other two-level systems (qubits): states, rotations, Born probabilities and
// projective measurement. Pure functions; a state is a spinor [aRe, aIm, bRe, bIm] = a|↑⟩ + b|↓⟩
// in the σ_z basis, normalised to |a|² + |b|² = 1.
//
//   ketAlong(n, s)            eigenstate of σ·n̂ with eigenvalue s = ±1 (spin up/down along n̂)
//   rotate(ψ, n, φ)           exp(−i φ σ·n̂/2) ψ: precession by φ about n̂ (Larmor: φ = ω t, ω = gμ_B B/ħ)
//   bloch(ψ)                  Bloch vector ⟨σ⟩ (a unit vector for a pure state)
//   probability(ψ, n, s)      Born probability of the outcome s = ±1 along n̂: (1 + s n̂·⟨σ⟩)/2
//   measure(ψ, n, rand)       sample an outcome; returns { s, state } with the collapsed state
//
// For a moving particle with a spin (Stern–Gerlach, neutron interferometry, NMR) the usual
// semiclassical model is a classical centre-of-mass path per branch plus this internal state,
// with the branch chosen by the Born rule where the beams separate. Say so beside the view.

const unit = n => { const l = Math.hypot(n[0], n[1], n[2]); if (!(l > 0)) throw new Error('spin: direction must be non-zero'); return [n[0] / l, n[1] / l, n[2] / l]; };

export function ketAlong(n, s = 1) {
  const [x, y, z] = unit(n);
  const th = Math.acos(Math.max(-1, Math.min(1, z))), ph = Math.atan2(y, x);
  // |+n⟩ = (cos θ/2, e^{iφ} sin θ/2), |−n⟩ = (sin θ/2, −e^{iφ} cos θ/2)
  const c = Math.cos(th / 2), sn = Math.sin(th / 2);
  return s > 0 ? Float64Array.from([c, 0, sn * Math.cos(ph), sn * Math.sin(ph)])
    : Float64Array.from([sn, 0, -c * Math.cos(ph), -c * Math.sin(ph)]);
}

export function rotate(psi, n, phi) {
  const [x, y, z] = unit(n), c = Math.cos(phi / 2), s = Math.sin(phi / 2);
  // U = c·1 − i s (σ·n) = [[c − i s z, −i s (x − i y)], [−i s (x + i y), c + i s z]]
  const [ar, ai, br, bi] = psi;
  const u11r = c, u11i = -s * z, u12r = -s * y, u12i = -s * x, u21r = s * y, u21i = -s * x, u22r = c, u22i = s * z;
  return Float64Array.from([
    u11r * ar - u11i * ai + u12r * br - u12i * bi, u11r * ai + u11i * ar + u12r * bi + u12i * br,
    u21r * ar - u21i * ai + u22r * br - u22i * bi, u21r * ai + u21i * ar + u22r * bi + u22i * br,
  ]);
}

export function bloch(psi) {
  const [ar, ai, br, bi] = psi;
  return [2 * (ar * br + ai * bi), 2 * (ar * bi - ai * br), ar * ar + ai * ai - br * br - bi * bi];
}

export function probability(psi, n, s = 1) {
  const u = unit(n), b = bloch(psi);
  return Math.min(1, Math.max(0, (1 + s * (u[0] * b[0] + u[1] * b[1] + u[2] * b[2])) / 2));
}

export function measure(psi, n, rand) {
  const s = rand() < probability(psi, n, 1) ? 1 : -1;
  return { s, state: ketAlong(n, s) };
}
