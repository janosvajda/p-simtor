# Gravitational waves and colliding black holes

Status: verified tools in `core/gravitational-waves.mjs` (tested in `tests/gravitational-waves.test.mjs` against the Newtonian limit, the vacuum field equations, the quadrupole formula, the Hulse–Taylor pulsar and GW150914's chirp mass). Units G = c = 1 for binaries: masses, lengths and times in the total mass M (1 M = GM/c² ≈ 1.48 km per M☉, or GM/c³ ≈ 4.93 µs per M☉).

## What can be computed in a browser, and what cannot

| Stage | Physics | Here? |
|---|---|---|
| Inspiral, separation ≫ M | Leading-order (quadrupole) radiation, Peters inspiral; post-Newtonian corrections of order v² = M/a | **Yes**, `circularInspiral`; disclose that v² ≈ 0.07 at 14 M and 0.17 at the ISCO |
| The field around the pair | Linearised tidal field E_ij (tendexes) in near, transition and wave zones | **Yes**, `binaryTidalField`; frame-drag field B_ij (vortexes) is not provided |
| Late inspiral, merger, ring-down | Full nonlinear Einstein equations (numerical relativity, supercomputers) | **No**. Stop the run at the ISCO (6 M) or earlier, with a marker naming the merger as not computed. Show the merger only as labelled published data (for example the GW150914 waveform from GWOSC, Abbott et al. 2016) if at all |
| Ring-down of the final hole | Quasinormal modes of a Kerr hole (needs spin-weighted spheroidal harmonics and Leaver's method) | Not provided; a single damped sinusoid with tabulated frequencies is a labelled fit, not a computation |
| Near a singularity (BKL/Mixmaster) | Homogeneous cosmological models, an ODE system | Computable as a separate topic (a set of ODEs), not part of a black-hole collision |

A request about "the nonlinear dynamics of curved spacetime" (Thorne's colloquium; Owen et al. PRL 106, 151101 (2011)) is about the merger, which no page here can compute. Build the inspiral and its field, say in the `Fact check:` line that the nonlinear stage is not computed, and end the run where it starts (SKILL.md step 7, *most representative*).

## API (`core/gravitational-waves.mjs`)

- `circularInspiral({ m1, m2, a0, aEnd = 6M })` → `{ a(t), phi(t), omega(t), tOf(a), tEnd, tCoalesce, luminosity(t), orbitsLeft(t), strain(t, r, ι), mu, eta, M }`. Exact closed-form solution of the leading-order equations: drive the page's state from it (it is computed, not keyframed; there is no step and no convergence order: verify it against its ODE, see verification.md).
- `binaryTidalField(ins, t, x, y, z)` → E_ij (row-major, 9 numbers) at the retarded time, adiabatic in the slowly changing orbit; the static M/r part is excluded (add it if the near zone matters). `quadrupoleTidalAmplitude(x, y, z, ω, Q)` is the general oscillating-quadrupole field.
- `radiatedPowerFromField(ω, Q)`: the energy flux of the computed field through a far sphere. Use it in `verify()`: it must equal `quadrupoleLuminosity(a, m1, m2)` (agrees to 5e-9).
- `chirpMass`, `symmetricMassRatio`, `quadrupoleLuminosity`; SI: `gwFrequencySI(a, M☉)`, `iscoFrequencySI(M☉)`, `chirpRateSI(f, Mc☉)`, `timeToCoalescenceSI(f, Mc☉)`.
- `periodDecay({ Pb, e, mp, mc })` and `eccentricityEnhancement(e)`: the binary-pulsar test (Peters & Mathews 1963).

## Real-world checks for `verify()`

- **Hulse–Taylor pulsar** (PSR B1913+16): P_b = 0.322997448918 d, e = 0.6171340, m_p = 1.4398, m_c = 1.3886 M☉ → dP_b/dt = −2.4026×10⁻¹²; observed (Galaxy-corrected) ratio to GR 0.9983 ± 0.0016 (Weisberg & Huang 2016, ApJ 829, 55).
- **GW150914**: 36 + 29 M☉ (source frame), chirp mass 28 (+1.8 −1.5) M☉; over 0.2 s the signal swept upwards from 35 to 250 Hz, its amplitude peaking near 150 Hz (Abbott et al. 2016, PRL 116, 061102). The leading-order ISCO frequency of 65 M☉ is 68 Hz (detector-frame masses are ~10 % higher, lowering it further): everything above is the late inspiral and merger, which leading-order formulas do not reach.
- **Equal masses, 14 M → 6 M**: 2900 M of time and 12.84 orbits.

## Pitfalls

- **Signs and phases of h₊, h×** differ between texts. The module's strain is h_ij = (2/r)Ï_ij^TT with φ from the x axis; E_ij = −½ ḧ_ij far away. Only matters when you compare the two on one page.
- **The field falls as 1/r** in the wave zone: a picture of E_ij over a wide region is dominated by the centre. Scale by r (and say so in the legend) or use a symmetric log colour scale; never clip silently.
- **Retardation is essential**: the spiral arms of the field are the retarded phase 2φ(t − r). Evaluating at t instead of t − r draws a rigidly rotating bar and loses the waves.
- **Adiabatic approximation**: the field uses the orbit's a and ω at the retarded time as if constant over one wave period; valid while the orbit shrinks little per orbit (true until close to the ISCO). Disclose it.
- **Sizes**: holes drawn at their true radius (2m each) are tiny at 14 M separation; enlarge and label "(enlarged)".
- **Exaggerated ring distortion**: real strains near the holes are ~0.1, far away ~10⁻²¹ on Earth; any visible distortion of test rings is enlarged and must say so.
