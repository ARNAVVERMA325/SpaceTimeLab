# Spacetime Lab — Milestone Status

**Companion to:** `CLAUDE.md` (project rules) and `ROADMAP.md` (implementation roadmap).

**Purpose.** A fresh Claude Code session starting mid-milestone has no memory of prior
sessions. This file is the authoritative record of *validated* state per milestone and
sub-gate, per `ROADMAP.md` § "Per-milestone tracking".

**Update policy.** Update at the **close of each gate**, not continuously. This file
records validated state, not work-in-progress. A gate moves to `gate passed` only when
every validation test listed under it is passing — not when the code merely runs or the
render looks correct (`CLAUDE.md` §16: "A result is not considered validated merely
because it 'looks right.'").

**Last updated:** 2026-09-24 — M3 closed: observers in motion, frequency shift, Novikov-Thorne thin disk, colorimetry, and parallel CPU rendering.

---

## State legend

| State | Meaning |
| --- | --- |
| `not started` | No implementation work begun. |
| `in progress` | Implementation begun; validation gate **not** closed. |
| `gate passed` | Every validation test for this gate is passing. |

---

## Summary

| Milestone / sub-gate | State |
| --- | --- |
| M1 — Minkowski baseline & engine plumbing | **`gate passed`** |
| M2A — Exterior Schwarzschild, CPU reference | **`gate passed`** |
| M2B — WebGPU/WGSL parallelization | `not started` |
| M3 — Physical observer, tetrad frames & accretion disk | **`gate passed`** |
| M4A — Kerr, Boyer–Lindquist exterior | `not started` |
| M4B — Kerr–Schild horizon-penetrating integration | `not started` |
| M5A — 3+1 embedding & GWOSC strain (committed scope) | `not started` |
| M5B — SXS / EHT (stretch, exploratory) | `not started` |

Milestones 1, 2A and 3 are closed: all of their validation tests pass in CI, 261 tests
across 23 files. Everything below M3 remains unimplemented, and those rows read
`not run` because the corresponding physics code does not exist yet.

M3 was taken out of order, ahead of M2B. The roadmap sequences 2B before 3, but 2B is a
WebGPU port whose only purpose is speed, and porting a pipeline that did not yet have
observers in motion, frequency shift or an emission model would have meant porting it
twice. M2B's prerequisite is unchanged and still met: the CPU reference is the ground
truth it will be checked against, and it is now a larger reference.

**M2B is unblocked.** The roadmap is explicit that no shader work starts before the
2A gate closes, and it has: the CPU reference is the ground truth every GPU result will
be checked against.

---

## Declared conventions

These must be fixed **before** M1.1 lands, since `CLAUDE.md` §2 forbids silently mixing
conventions and §19 warns that renaming later touches every downstream module.

All conventions are declared in one place, `src/physics/conventions.ts`, and read from
there by every module rather than assumed locally.

| Convention | Value | Status |
| --- | --- | --- |
| Metric signature | `(-, +, +, +)` — fixed by `CLAUDE.md` §3 (`g_mu_nu u^mu u^nu = -1`) | in code, asserted by tests |
| Unit system | Geometric units, `G = c = 1`; horizon at `r = 2M` (`CLAUDE.md` §6.1) | in code |
| Coordinate ordering | Carried per chart on the model (`CoordinateChart.coordinateNames`); `(t, x, y, z)` for Minkowski Cartesian, `(t, r, theta, phi)` for spherical-type charts from M2A | in code |
| Integration parameter | Affine parameter `lambda` for null geodesics; proper time `tau` for timelike. Enforced by `parameterName(kind)` | in code, asserted by tests |
| Index / naming convention | `CLAUDE.md` §19 names used verbatim (`g_mu_nu`, `g_inv_mu_nu`, `christoffel`, `four_velocity_u`, `null_wavevector_k`, `energy_E`, `angular_momentum_Lz`) | in code |
| Index position (variance) | Tracked on `FourVector`; raising an already-raised index throws | in code, asserted by tests |
| Floating point | IEEE-754 binary64 on the CPU reference path — higher precision than f32, still finite (`CLAUDE.md` §8) | in code |
| Carter constant convention | Not yet chosen — required by `CLAUDE.md` §16 before M4A | **open** |
| Chart kind | Recorded per chart (`cartesian` / `spherical`), since some routines are valid for only one | in code, enforced |

---

## M1 — Minkowski baseline & engine plumbing

**State:** `gate passed` — 84 tests across 8 files, green in CI.

| Validation test (`CLAUDE.md` §16) | Status | Measured |
| --- | --- | --- |
| Christoffel symbols identically zero in Minkowski (analytical check) | `pass` | exactly `0` |
| Null normalization `g_mu_nu k^mu k^nu = 0` | `pass` | ≤ `4.7e-16` |
| Timelike normalization `g_mu_nu u^mu u^nu = -1` | `pass` | within `1e-14` relative |
| Straight-line propagation in flat space across 10^6 steps | `pass` | `1.1e-14` relative over coordinate distance `1.4e3` |
| Flat-space limit: connection vanishes, curvature vanishes | `pass` | exactly `0` |
| `energy_E` and `angular_momentum_Lz` conservation | `pass` | exactly `0` drift over 10^6 steps |
| Convergence: RK4 global error under timestep halving | `pass` | ratio in `[14, 18]`, i.e. 4th order |
| Tetrad orthonormality `g_mu_nu e^mu_(a) e^nu_(b) = eta_ab` | `pass` | exactly `0` |
| Backward-traced image vs. direct analytic projection | `pass` | identical, pixel for pixel |
| Wavevector preserved bit-for-bit along a flat-space ray | `pass` | exact equality |

**Definition of Done:** unit test suite passes in CI; canvas renders an undistorted grid
via backward ray tracing. **Both met.**

**What landed:**

- `MetricTensor`, `ChristoffelSymbols`, `FourVector`, `PhaseSpaceState`, with `CLAUDE.md`
  §19 naming from the first commit rather than renamed later.
- Minkowski model with analytically zero Christoffel symbols, plus a general numerical
  connection from central-differenced metric derivatives for models that will not have
  closed-form symbols.
- Fixed-step RK4 and adaptive RKF45 behind an integrator interface (`CLAUDE.md` §7.3),
  so a symplectic method can be added for long-lived bound orbits without touching call
  sites.
- An integration driver owning termination policy, explicit NaN/infinity detection, and
  `CLAUDE.md` §17 failure diagnostics (integrator, model, chart, parameter, step size,
  failed quantity).
- General 4x4 metric inversion with partial pivoting, exercised on off-diagonal metrics
  now so that Kerr in M4A does not meet it for the first time inside a render.
- Observer layer with a validated orthonormal tetrad and past-directed null ray
  generation.
- Backward null-geodesic ray tracer and Canvas2D output, with a provenance panel
  covering every item `CLAUDE.md` §22 requires.

**Decisions closed:**

- Flat-space propagation tolerance: `1e-12` relative, set from a measured `1.1e-14` with
  about two orders of magnitude of headroom. Justification is recorded alongside the
  value in `src/physics/validation/tolerances.ts`.
- RKF45 propagates the **fifth-order** solution (local extrapolation), so the embedded
  difference estimates the error of the fourth-order solution that is *not* propagated
  and acts as a conservative proxy. Stated explicitly because it changes what the
  reported error norm means; Fehlberg's original formulation propagates the fourth-order
  solution instead.
- Backward-traced rays use a **past-directed** wavevector, `k^(a) = (-1, n)` in the
  observer frame. Geometrically equivalent to the future-directed choice, since time
  reversal maps null geodesics to null geodesics, but the sense matters for the
  emitter-to-observer frequency shift in M3 and is fixed now to avoid a sign error there.

**Known limitations, deliberately not papered over:**

- The observer tetrad is the identity frame, which is correct and exactly orthonormal in
  Cartesian Minkowski but is *not* the general static-observer construction. That is a
  M3 deliverable (`ROADMAP.md` 3.1) and has not been faked.
- No frequency shift, emission model or radiative transfer. M1 has no emitting matter and
  no relative motion, so there is nothing to shift — steps 5 and 6 of `CLAUDE.md` §9 are
  absent rather than approximated.
- Render throughput is roughly 390k rays in about 7.5 s on the reference machine. That is
  a CPU reference path, not an interactive target; 60 FPS is a M2B goal and `CLAUDE.md`
  §21 is explicit that it is a rendering goal, not a scientific-validity requirement.
- Curvature visible in the rendered grid lines is rectilinear projection of a sphere, not
  light deflection. The UI says so, and the projection property it rests on — that a
  pinhole camera maps great circles, and only great circles, to straight image lines — is
  asserted by tests.

---

## M2A — Exterior Schwarzschild, CPU reference

**State:** `gate passed`
**Depends on:** M1 (passed).

| Validation test (`ROADMAP.md` 2A.4) | Status | Measured |
| --- | --- | --- |
| Photon sphere locks at `r = 3M` | `pass` | exact fixed point in the Hamiltonian form, held over ~27 orbits |
| Photon-sphere Lyapunov exponent vs. `1/(3 sqrt(3) M)` | `pass` | `1.2e-5` relative; growth per half orbit `e^pi` |
| Capture threshold vs. `b_c = 3 sqrt(3) M` | `pass` | `1.8e-10` relative, by bisection |
| `energy_E = -p_t` relative drift `< 1e-6` | `pass` | **exactly zero** — `p_t` bit-identical after integration |
| `angular_momentum_Lz = p_phi` relative drift `< 1e-6` | `pass` | **exactly zero** — `p_phi` bit-identical after integration |
| Weak-field deflection vs. analytical `alpha ~ 4GM/(c^2 b)` | `pass` | ratio → 1 as `1/b`: 1.0302, 1.0030, 1.00029, 1.000029 |
| Traced deflection vs. exact orbit-equation quadrature | `pass` | `2e-13` to `6.4e-9` relative, `r_0` from 3.2M to 10^4 M |
| Flat-space limit recovered as `M -> 0` (`CLAUDE.md` §16) | `pass` | mass-sourced symbols vanish linearly in M |
| Analytic vs. numerically differenced Christoffels | `pass` | `<= 1.6e-8` relative |
| Null normalization along traced rays | `pass` | `<= 5.1e-10` over a full image at render tolerance `1e-10` |
| Background direction independent of stopping radius | `pass` | `< 1e-10` from R = 100M to 5000M; uncorrected bias falls as `1/R^2` |
| Shadow angular radius vs. `sin(psi) = b_c sqrt(f)/r` | `pass` | `< 1e-6` relative |
| Orbital-plane reduction vs. full 3D integration | `pass` | exit directions agree to `< 1e-8`; ray stays in-plane to `< 1e-9` |

**Extended validation** — items the roadmap or `CLAUDE.md` §16 name that the first pass
through 2A did not actually deliver. Every reference below was computed independently in
mpmath at 40–60 digits, sharing no code with the TypeScript under test.

| Validation | Status | Measured |
| --- | --- | --- |
| Einstein ring, on-axis source, exact finite-distance lens equation (2A.3) | `pass` | observed angle to `< 1e-12` rad for rings n = 0–3 at three observer/source distances |
| Relativistic ring spacing `(b_{n+1} - b_c)/(b_n - b_c)` → `e^{-2 pi}` | `pass` | `1e-6` vs. mpmath at n = 1→2; `1e-4` vs. `e^{-2 pi}` at n = 2→3 |
| Primary ring vs. weak-field `theta_E = sqrt(4M D_LS/(D_OL D_OS))` | `pass` | ratio 1.107 → 1.033 → < 1.02 at 100M, 1000M, 10^4 M |
| Strong-deflection limit, Bozza 2002: `alpha + ln(b/b_c - 1) -> b_bar = -0.40023` | `pass` | residual shrinks monotonically to `< 1e-5` at `b/b_c - 1 = 1e-6` |
| Near-critical deflection accuracy | `pass` | bounded by a fixed effective impact-parameter error `< 2e-12`; see below |
| Circular orbits: Kepler's law `dphi/dt = sqrt(M/r^3)`, exact in these coordinates | `pass` | `< 1e-11` at r = 7, 10, 20, 100 M |
| ISCO at 6M, `E = sqrt(8/9)`, `L = 2 sqrt(3) M` | `pass` | closed form, and E, L minimal there along the circular family |
| Radial epicyclic frequency `Omega sqrt(1 - 6M/r)` (stable side) | `pass` | `< 1e-6` at p = 8, 10, 20 M |
| Instability growth rate `sqrt(M(6M - r))/r^2` (unstable side) | `pass` | `< 1e-3` at r = 5M |
| Periapsis precession vs. exact elliptic-integral formula | `pass` | `2e-14`–`4e-13` relative, p from 10M to 1000M |
| Precession → Einstein's `6 pi M/p` in the weak field | `pass` | ratio `< 1 + 1e-3` at p = 10^4 M |
| Radial free fall: proper time and coordinate time | `pass` | `< 1e-11` and `< 1e-10` vs. closed forms |
| Symplectic Gauss–Legendre vs. RK4 over 300 orbits (`CLAUDE.md` §7.2–7.3) | `pass` | mass-shell error flat at `1.0e-10` vs. RK4 growing linearly to `2.1e-8` |

A result worth recording from this work: **near the critical curve the deflection is
ill-conditioned, and the error follows a law rather than drifting.** With
`alpha = -ln(b/b_c - 1) + ...`, `d alpha/db = -1/(b - b_c)`, so any integration error acts
as a fixed effective impact-parameter error `delta b`, amplified by `1/(b - b_c)`. At
tolerance `1e-12`, `|delta alpha| (b - b_c)` is `6.6e-13` at every `b/b_c - 1` from `1e-4`
to `1e-10`, and the same `delta b` appears independently as the absolute error of every
ring solution. The relative error in `alpha` therefore grows as the ray approaches `b_c`
(to `5.6e-5` at `b/b_c - 1 = 1e-10`), and that is the problem's sensitivity, not a defect.
The bound is placed on `delta b`, the quantity the integrator controls.

**Definition of Done:** all of 2A.4 passes on the CPU reference. **Met.** This is the
ground truth every later GPU or Kerr result is checked against.

**What landed:**

- Schwarzschild metric, inverse and analytical Christoffel symbols in `(t, r, theta, phi)`.
  The symbols were derived symbolically from `CLAUDE.md` §2's definition and cross-checked
  by confirming `R_mu_nu = 0` (vacuum) and `K = 48 M^2 / r^6`.
- An independent analytic reference: the exact deflection quadrature, the finite-radius
  tail integral, and the closed-form photon-sphere constants. Validated against values
  computed separately before anything was checked against it.
- Exterior domain checking, a static-observer tetrad, an inward-facing pinhole screen,
  the orbital-plane reduction, and a chart-agnostic raytracer producing the shadow and
  the lensed background.

**Engine, rebuilt after the gate first closed** (2026-09-23). A second pass, looking
for what a referee would object to, found the gate passing on an engine with avoidable
weaknesses. All M2A rows above were re-measured on the new engine.

| Change | Why | Measured effect |
| --- | --- | --- |
| Hamiltonian formulation `(x^mu, p_mu)`, `H = 1/2 g^{mu nu} p_mu p_nu`, now the default | `CLAUDE.md` §2 prefers it; momenta conjugate to ignorable coordinates are conserved by construction | `E`, `L_z` exact; 12x smaller deflection error at equal tolerance, 30% fewer steps. The null constraint becomes the one drifting invariant, and is reported |
| Dormand–Prince 5(4) with FSAL and Shampine's continuous extension, replacing Fehlberg 4(5) as the default | Fehlberg's error constants suit the fourth-order solution, not the propagated fifth | 5th order confirmed; dense output `O(h^5)`; at every tolerance tried, smaller error for fewer steps than Fehlberg |
| Exact event location: Brent on the dense output, then a Newton-refined landing step | Terminations previously overshot by up to a step — the cause of the 150% deflection error found during 2A | Events land on their surface to `< 1e-11`, independent of tolerance |
| Exact asymptotic background direction via the orbit-equation tail integral | Sampling the sky at finite `R` biased every ray by the deflection still to come | Bias removed; result independent of `R` to `1e-10` |
| Stratified supersampling in linear light | `CLAUDE.md` §9 step 7 | Converges at the stratified-sampling rate `O(N^-3/4)`, ratios 2.5–3.0 per doubling |
| Curved-space convergence test | `CLAUDE.md` §16 asked for convergence; it existed only for a harmonic oscillator | RK4 4th order on a Schwarzschild deflection: ratios 15.8, 15.9 |

Combined, a full Schwarzschild render at tolerance `1e-10` is 2.8x faster with a 33x
smaller worst null residual than before, which lets the interactive render meet the
reference `1e-9` gate. The separate "preview" tolerance introduced when the gate first
closed has been retired.

Two tests were found to be checking numerical artifacts rather than physics and were
replaced: the photon-sphere "instability" test relied on Christoffel rounding to seed the
departure (the Hamiltonian form holds the orbit exactly, which is the correct solution),
and is now a measurement of the Lyapunov exponent.

**Decisions closed:**

- **Horizon detection** reads `f = 1 - 2M/r` from the metric components in use, never the
  determinant. `CLAUDE.md` §6.1 forbids the determinant test, and the suite asserts why:
  `det g = -r^4 sin^2(theta)` is about `-16 M^4` at `r = 2M`, nowhere near zero, while `f`
  has vanished.
- **Ray capture** uses the exact condition `r < 3M` with `k^r < 0`, not a tuned radius.
  The null effective potential `f/r^2` increases inward of `3M`, so an inward-moving
  photon there can never turn around. Rays stop while the chart is still well behaved
  instead of being integrated toward `r = 2M`.
- **Polar-axis handling.** Spherical charts degenerate on the axis. The renderer works
  in each ray's own orbital plane — exact under spherical symmetry — so `sin(theta) = 1`
  throughout and the axis is never approached. Cross-validated against full 3D
  integration, and refused outright for a chart or a symmetry that does not support it.
- **One tolerance for render and tests.** The interactive render integrates at `1e-10` and
  is judged against the same `null-normalization-traced` gate (`1e-9`) as the test suite,
  measured `5.1e-10`. The earlier separate preview tolerance is retired.

**Known limitations, deliberately not papered over:**

- Throughput is about 3.5k rays/s single-threaded on the reference machine (30k rays in
  8 s). This is the CPU reference path; 60 FPS is M2B's goal and `CLAUDE.md` §21 is
  explicit it is not a validity requirement.
- The finest bands near the shadow edge alias at any sample count: each successive band
  is `e^pi ~ 23` times thinner than the last, so no finite ray count resolves them all.
  Supersampling averages them correctly in linear light; it does not and cannot resolve
  them.
- `equatorialNullRay` builds `k^r` from `E^2 - f L^2/r^2`, which cancels catastrophically
  for a ray both near-critical and near `r = 3M` (the difference falls below rounding at
  `delta r ~ 1e-9 M`). It throws rather than returning a wrong ray; the Lyapunov test
  builds its initial data from the factored form instead.
- The exact-deflection quadrature is refused for `r_0 < 3.05M`, where the integrand
  approaches a double root, rather than returning a quietly inaccurate value. It also
  loses *relative* accuracy in the far weak field, where `alpha = 4*Integral - pi`
  cancels almost completely; absolute accuracy stays near `1e-15`.
- Throughput above is single-threaded. The parallel CPU path landed with M3 and gives
  roughly the core count; it is still not M2B.
- The static-observer tetrad at the close of 2A was only the shell observer. Freely-falling
  observers, frequency shift, beaming and the accretion disk landed with M3.

---

## M2B — WebGPU/WGSL parallelization

**State:** `not started` — unblocked; M2A's gate is closed.
**Depends on:** M2A (passed). The CPU reference is now the cross-validation target.

| Validation test (`ROADMAP.md` 2B.2) | Status |
| --- | --- |
| GPU photon sphere vs. CPU reference | `not run` |
| GPU deflection angle vs. CPU reference | `not run` |
| GPU conserved-quantity drift vs. CPU reference | `not run` |
| Required shader numeric features confirmed on target devices | `not run` |

**Definition of Done:** interactive shadow + Einstein-ring rendering at 60 FPS **and**
GPU results matching the CPU reference within a documented, justified tolerance.

**Open decisions:**

- **GPU/CPU tolerance gap (f32 shader vs. f64 CPU reference): not yet chosen.** This is
  the tolerance decision `ROADMAP.md` calls out by name. It must be documented and
  justified here when M2B opens, not smoothed over (`CLAUDE.md` §8, §17, §24).

---

## M3 — Physical observer, tetrad frames & accretion disk

**State:** `gate passed` — 2026-09-24.

Every reference value below is either a closed form or an independent second route through
the code; nothing is checked against a stored output of the same function.

| Validation test | Status | Measured |
| --- | --- | --- |
| Tetrad orthonormality `g_mu_nu e^mu_(a) e^nu_(b) = eta_ab` (static, boosted, free-fall) | `pass` | `< 1e-14` at r = 3.5M to 1000M |
| Static-observer normalization `u^mu = (1/sqrt(f), 0, 0, 0)` | `pass` | to 14 digits; spatial components exactly `0` |
| Free-faller's `u^mu = (1/f, -sqrt(2M/r), 0, 0)`, and it is a geodesic | `pass` | agrees with the `E = 1` timelike geodesic to `1e-13` |
| Redshift factor vs. the analytical static-observer relation `sqrt(f_emit/f_obs)` | `pass` | `< 1e-12` on a traced radial ray, r = 8M to 20M |
| Frequency shift: static-frame route vs. covariant `u^t (p_t + Omega p_phi)` | `pass` | `< 1e-12` relative, on rays traced in the full world chart |
| Frequency shift on meridional rays vs. `sqrt(1 - 3M/r) / sqrt(f_obs)` | `pass` | `< 1e-10` |
| Infaller vs. hoverer at one event: local Lorentz factor `gamma (1 - beta . n)` | `pass` | `< 1e-14` for five directions |
| Infaller's view of the sky behind: `1 / (1 + sqrt(2M/r))` | `pass` | `< 1e-14` at r = 5M, 20M, 100M |
| Aberration of the shadow edge: `cos psi' = (cos psi + v)/(1 + v cos psi)` | `pass` | `< 1e-8` relative, by bisection on capture in the infaller's frame |
| Flat-space limit: no shift for co-stationary observer and emitter, any direction | `pass` | exactly `1` in Minkowski; `< 1e-15` at equal radius in Schwarzschild |
| Flat-space limit: orbital Doppler shift vanishes as `sqrt(M/r)` | `pass` | monotone, `< 1e-3` at M/r = 1e-7 |
| Novikov–Thorne flux vs. an independent mpmath evaluation of the closed form | `pass` | `< 1e-12` relative, r = 6M to 1000M |
| Disk energy balance: `L_inf / Mdot = 1 - sqrt(8/9)` | `pass` | `< 1e-9`; efficiency `0.0571909584179366` |
| Flux peak radius `r = 9.550928 M`, zero torque at the ISCO | `pass` | `< 2e-6`; `F(6M) = 0`, `F >= 0` everywhere |
| Far-field flux `-> (3 M Mdot / 8 pi r^3)(1 - C sqrt(M/r))`, `C = sqrt6 + sqrt3 ln(1 + sqrt2)` | `pass` | monotone convergence to `< 1e-3` |
| Disk temperature in physical units vs. an independent CODATA 2018 calculation | `pass` | `< 1e-9` relative; scales as `(f_Edd / M)^(1/4)` to 12 digits |
| Planck spectrum vs. the Stefan–Boltzmann and Wien laws (neither used in the code) | `pass` | `< 1e-9` and `< 1e-6` |
| Blackbody chromaticity vs. `colour-science`'s independent 1 nm integration | `pass` | `< 1e-4` in CIE 1931 (x, y) |
| Fast colour table vs. direct 5 nm integration | `pass` | `< 5e-6` relative on X, Y, Z over 1000–60,000 K |
| Rendered disk: both redshift and blueshift at 80 degrees | `pass` | `g` from `0.49` to `1.44` |
| Rendered disk: approaching side brighter (Doppler beaming) | `pass` | `4.3x` in display-linear luminance at 8,830 K |
| Rendered disk: the same beaming is weak for a hot disk (Rayleigh–Jeans band) | `pass` | `1.8x` at 39,500 K, identical kinematics |
| Rendered disk: far side lensed over the shadow (Luminet 1979) | `pass` | lit pixels above the shadow centre |
| Rendered disk: face-on view symmetric and wholly redshifted | `pass` | `g_max < 1`; left/right within `10%` |
| Parallel render bit-identical to serial | `pass` | every byte and every linear-radiance value, all four scene kinds |
| Null normalization over a whole disk image | `pass` | `<= 6.0e-10` against the `1e-9` gate |

**Definition of Done:** toggling static vs. free-falling observer visibly and correctly
changes beaming/Doppler shift; disk shows the classical asymmetric intensity pattern.
**Met**, and both halves are checked numerically rather than by eye: the aberration and
Doppler rows above are the "correctly", and the app exposes the toggle.

**What landed:**

- **Observers in motion.** `boostTetrad` Lorentz-boosts an orthonormal frame in place;
  the freely-falling observer is the static frame boosted by `-sqrt(2M/r) r-hat`. Both
  the aberration of the shadow and the shift of the sky are consequences, not special
  cases in the renderer.
- **Frequency shift** `g = (k . u_obs) / (k . u_emit)`, computed once and carrying
  gravitational redshift, transverse Doppler and line-of-sight Doppler together. Rays are
  launched with `k . u_obs = 1` by construction, so `g = 1 / (k . u_emit)`, evaluated in
  the static orthonormal frame at the emission event — which keeps it chart-independent
  and indifferent to the orbital-plane reduction.
- **Novikov–Thorne / Page–Thorne thin disk** in closed form between the ISCO and 20M,
  with the physical temperature scale set by the black-hole mass and an Eddington
  fraction, from CODATA 2018 and IAU 2015 values.
- **Colorimetry rather than a colour ramp.** The observed spectrum of a blackbody at `T`
  seen with ratio `g` is a blackbody at `g T`; that spectrum is integrated against the
  CIE 1931 2-degree colour-matching functions and converted to sRGB per IEC 61966-2-1.
  The CMF table is generated from `colour-science` by `scripts/generate_cie_cmf.py` and
  carries its provenance in the file.
- **A disclosed display mapping.** Exposure is relative to the luminance of the hottest
  ring seen at rest; the Reinhard curve is applied at encode only. The linear radiance
  behind every pixel is kept unmodified, so the tone curve cannot leak into a measurement.
- **Parallel CPU rendering.** Rows are dealt round-robin to a pool of Web Workers, which
  rebuild the scene from a plain-data `SceneDescription` — the same function the main
  thread uses. Interleaved rather than banded, because rays near the shadow edge cost many
  times more than rays that miss the hole.

**Measured performance** (4 cores in this container, 200x150 at 1 sample/px, tolerance
`1e-10`): disk scene 4.1 s, Schwarzschild sky 3.6 s, Minkowski 1.1 s — about 7.3k rays/s
against 3.5k single-threaded. A 320x240 disk render at 4 samples/px takes 30 s.

**Known limitations, deliberately not papered over:**

- The disk is a **Novikov–Thorne model, not an accretion simulation**: no disk
  atmosphere, no electron scattering or spectral hardening, no limb darkening, no
  self-irradiation or returning radiation, no emission inside the ISCO, no plasma and no
  radiative transfer. It does not affect the spacetime. Thin disks are a reasonable model
  only at moderate accretion rates. The UI says all of this.
- **Colour is not a measurement.** The sRGB values are a colorimetric rendering of a
  computed spectrum, subject to gamut clamping, an exposure choice and a tone curve.
  The frequency ratio and the observed temperature are the physical outputs; the picture
  is a display of them.
- The background grid still has **no spectrum**, so it is drawn unshifted. A ray that
  ends on the sky carries no radiometry — only disk light does.
- What a colour image shows depends on where the visible band sits on the Planck curve:
  the bolometric boost is always `g^4`, but in-band luminance goes as `g` deep in the
  Rayleigh–Jeans tail. A 40,000 K disk therefore looks almost uniform however fast its
  gas moves. This is measured (the two beaming rows above), and stated in the UI, because
  it is the kind of thing that otherwise looks like a rendering bug.
- Workers give the CPU path roughly the core count. **This is not M2B**: WebGPU is
  untouched, and 60 FPS remains out of reach on this path.

---

## M4A — Kerr, Boyer–Lindquist exterior (`r > r_+`)

**State:** `not started` — unblocked; M2A's gate is closed.
**Note:** Kerr is axisymmetric but *not* spherically symmetric, so its geodesics do not
lie in planes through the centre and the orbital-plane reduction does not apply. The
reduction refuses a model that does not declare spherical symmetry, so this cannot be
reached for by accident. Kerr in Boyer–Lindquist also has a `g_t_phi` cross term, which
the diagonal static-tetrad construction refuses; M4A needs a genuine orthonormalization.

| Validation test (`ROADMAP.md` 4A) | Status |
| --- | --- |
| `carter_Q` conservation (documented convention) | `not run` |
| `energy_E`, `angular_momentum_Lz` conservation outside `r_+` | `not run` |
| Frame-dragging precession vs. analytical limit | `not run` |
| Schwarzschild limit recovered as `a -> 0` | `not run` |
| Asymmetric shadow vs. Thorne/DNGR reference values | `not run` |

**Definition of Done:** exterior Kerr results validated and stable; shadow benchmarked
against published reference values.

**Open decisions:**

- Carter constant convention: not yet chosen. `CLAUDE.md` §16 requires exactly one
  clearly documented convention.

---

## M4B — Kerr–Schild horizon-penetrating integration

**State:** `not started`
**Blocked on:** M4A gate. **Do not open this gate until M4A's conservation tests are
green** — that is the whole reason this file exists (`ROADMAP.md` § "Per-milestone
tracking").

| Validation test (`ROADMAP.md` 4B) | Status |
| --- | --- |
| Traversal through `r_+` with no NaN / overflow / shader crash | `not run` |
| Kerr–Schild vs. Boyer–Lindquist agreement outside `r_+` | `not run` |

**Definition of Done:** smooth traversal through `r_+` **and** confirmed exterior
agreement with 4A. Treated as research-adjacent, not a routine port; budget time
disproportionate to its bullet count.

---

## M5A — 3+1 embedding & GWOSC strain (committed scope)

**State:** `not started`

| Validation / disclosure requirement | Status |
| --- | --- |
| `spatial_metric_gamma`, `lapse_alpha`, `shift_beta` extraction | `not run` |
| Flamm-paraboloid embedding labeled as an embedding of a spatial slice (`CLAUDE.md` §5.2) | `not done` |
| UI states coordinate system, foliation, time parameter, geometry type | `not done` |
| GWOSC strain labeled "linearized-wave visualization driven by measured detector strain" (`CLAUDE.md` §12) | `not done` |
| GW model assumptions stated: propagation direction, polarization basis, linearized regime, gauge, plane-wave approximation | `not done` |

**Definition of Done:** user can switch between Optical, 3+1 Embedding, and Curvature
heatmap views; GWOSC strain drives real-time grid oscillations under clearly stated
assumptions.

---

## M5B — SXS / EHT (stretch, exploratory)

**State:** `not started`

Treated as exploratory rather than a hard commitment until data access and licensing are
confirmed in practice (`ROADMAP.md` 5B).

| Disclosure requirement | Status |
| --- | --- |
| SXS data labeled "Numerical Relativity Data", never live browser calculation (`CLAUDE.md` §13) | `not done` |
| SXS provenance recorded: simulation ID, parameters, data product type, units, resolution, citation | `not done` |
| EHT framed as model-to-data comparison with omitted physics disclosed (`CLAUDE.md` §14) | `not done` |

---

## Open numerical-tolerance decisions (consolidated)

Carried here so a session picking up mid-project can see them in one place.

| # | Decision | Milestone | Status |
| --- | --- | --- | --- |
| 1 | Flat-space 10^6-step propagation tolerance | M1 | **closed** — `1e-12` relative, from a measured `1.1e-14` |
| 2 | Conserved-quantity drift budget | M2A | **closed** — `1e-7` relative, 100x tighter than the roadmap's `1e-6`; measured `<= 1.4e-9` |
| 3 | Weak-field deflection benchmark tolerance | M2A | **closed** — stated as convergence, not equality; traced vs. exact quadrature at `1e-7` relative |
| 4 | **GPU f32 vs. CPU f64 cross-validation tolerance** | M2B | open |
| 5 | Carter constant convention (not a tolerance, but must be fixed before M4A) | M4A | open |
| 6 | Kerr–Schild vs. Boyer–Lindquist exterior agreement tolerance | M4B | open |
| 7 | Preview vs. reference null-residual tolerance | M2A | **closed** — preview tolerance retired; render meets the `1e-9` reference gate at `5.1e-10` |
| 8 | Fast blackbody-colour table vs. direct spectral integration | M3 | **closed** — `5e-6` relative on X, Y, Z, from a measured worst case below `3e-6`; quadratic log-log interpolation at 256 nodes per decade |
| 9 | Disk energy-balance check (quadrature over an infinite domain) | M3 | **closed** — `1e-9` relative, after substituting `r = 1000/w^2` to remove the `sqrt(u)` tail behaviour |

`CLAUDE.md` §17: there is no universal numerical-error threshold. Each entry above must
be justified by the relevant numerical method and quantity when it is closed.
