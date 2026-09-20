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

**Last updated:** 2026-09-20 — Milestone 1 gate closed.

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
| M2A — Exterior Schwarzschild, CPU reference | `not started` |
| M2B — WebGPU/WGSL parallelization | `not started` |
| M3 — Physical observer, tetrad frames & accretion disk | `not started` |
| M4A — Kerr, Boyer–Lindquist exterior | `not started` |
| M4B — Kerr–Schild horizon-penetrating integration | `not started` |
| M5A — 3+1 embedding & GWOSC strain (committed scope) | `not started` |
| M5B — SXS / EHT (stretch, exploratory) | `not started` |

Milestone 1 is closed: all of its validation tests pass in CI. Everything below M1
remains unimplemented, and those rows read `not run` because the corresponding physics
code does not exist yet.

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

**State:** `not started` — unblocked; M1's gate is closed.
**Depends on:** M1 (passed).

| Validation test (`ROADMAP.md` 2A.4) | Status |
| --- | --- |
| Photon sphere locks at `r = 3M` | `not run` |
| `energy_E = -p_t` relative drift `< 1e-6` | `not run` |
| `angular_momentum_Lz = p_phi` relative drift `< 1e-6` | `not run` |
| Weak-field deflection vs. analytical `alpha ~ 4GM/(c^2 b)` | `not run` |
| Flat-space limit recovered as `M -> 0` (`CLAUDE.md` §16) | `not run` |

**Definition of Done:** all of the above pass on the CPU reference. This is the ground
truth every later GPU or Kerr result is checked against. **No shader work starts before
this gate closes.**

**Open decisions:**

- Horizon-detection predicate: must be derived from the metric / inverse-metric
  components actually in use. `CLAUDE.md` §6.1 explicitly forbids detecting the horizon
  via a vanishing metric determinant (the Schwarzschild determinant goes as
  `-r^4 sin^2(theta)` and does not vanish at `r = 2M`).

---

## M2B — WebGPU/WGSL parallelization

**State:** `not started`
**Blocked on:** M2A gate (hard gate — CPU reference must be green first).

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

**State:** `not started`
**Blocked on:** M2A gate.

| Validation test | Status |
| --- | --- |
| Tetrad orthonormality `g_mu_nu e^mu_(a) e^nu_(b) = eta_ab` | `not run` |
| Static-observer normalization `u^mu = (1/sqrt(-g_00), 0, 0, 0)` | `not run` |
| Redshift factor against analytical static-observer relation | `not run` |
| Flat-space limit: no beaming / no Doppler shift for a static observer | `not run` |

**Definition of Done:** toggling static vs. free-falling observer visibly and correctly
changes beaming/Doppler shift; disk shows the classical asymmetric intensity pattern.

---

## M4A — Kerr, Boyer–Lindquist exterior (`r > r_+`)

**State:** `not started`
**Blocked on:** M2A gate.

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
| 2 | Conserved-quantity drift budget (roadmap states `< 1e-6` relative — confirm per quantity, `CLAUDE.md` §17) | M2A | open |
| 3 | Weak-field deflection benchmark tolerance | M2A | open |
| 4 | **GPU f32 vs. CPU f64 cross-validation tolerance** | M2B | open |
| 5 | Carter constant convention (not a tolerance, but must be fixed before M4A) | M4A | open |
| 6 | Kerr–Schild vs. Boyer–Lindquist exterior agreement tolerance | M4B | open |

`CLAUDE.md` §17: there is no universal numerical-error threshold. Each entry above must
be justified by the relevant numerical method and quantity when it is closed.
