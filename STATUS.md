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

**Last updated:** 2026-09-20 — repository bootstrap; governance documents added.

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
| M1 — Minkowski baseline & engine plumbing | `not started` |
| M2A — Exterior Schwarzschild, CPU reference | `not started` |
| M2B — WebGPU/WGSL parallelization | `not started` |
| M3 — Physical observer, tetrad frames & accretion disk | `not started` |
| M4A — Kerr, Boyer–Lindquist exterior | `not started` |
| M4B — Kerr–Schild horizon-penetrating integration | `not started` |
| M5A — 3+1 embedding & GWOSC strain (committed scope) | `not started` |
| M5B — SXS / EHT (stretch, exploratory) | `not started` |

Nothing below has been implemented or validated yet. Every validation row reads
`not run` because no physics code exists in this repository.

---

## Declared conventions

These must be fixed **before** M1.1 lands, since `CLAUDE.md` §2 forbids silently mixing
conventions and §19 warns that renaming later touches every downstream module.

| Convention | Value | Status |
| --- | --- | --- |
| Metric signature | `(-, +, +, +)` — implied by `CLAUDE.md` §3 (`g_mu_nu u^mu u^nu = -1`) | declared, not yet enforced in code |
| Unit system | Geometric units, `G = c = 1`; horizon at `r = 2M` (`CLAUDE.md` §6.1) | declared, not yet enforced in code |
| Coordinate ordering | `(t, r, theta, phi)` for spherical-type charts (`ROADMAP.md` 2A.1) | declared, not yet enforced in code |
| Integration parameter | Affine parameter `lambda` for null geodesics; proper time `tau` for timelike | declared, not yet enforced in code |
| Index / naming convention | `CLAUDE.md` §19 names (`g_mu_nu`, `g_inv_mu_nu`, `christoffel`, `four_velocity_u`, `null_wavevector_k`, ...) | declared, not yet enforced in code |
| Carter constant convention | Not yet chosen — required by `CLAUDE.md` §16 before M4A | **open** |

---

## M1 — Minkowski baseline & engine plumbing

**State:** `not started`

| Validation test (`CLAUDE.md` §16) | Status |
| --- | --- |
| Null normalization `g_mu_nu k^mu k^nu = 0` | `not run` |
| Timelike normalization `g_mu_nu u^mu u^nu = -1` | `not run` |
| Straight-line propagation in flat space across 10^6 steps | `not run` |
| Christoffel symbols identically zero in Minkowski (analytical check) | `not run` |

**Definition of Done:** unit test suite passes in CI; canvas renders an undistorted grid
via backward ray tracing.

**Open decisions:**

- Tolerance for the 10^6-step flat-space propagation test: not yet chosen. Per
  `CLAUDE.md` §17, this must be justified per quantity — not a blanket `1e-5`.

---

## M2A — Exterior Schwarzschild, CPU reference

**State:** `not started`
**Blocked on:** M1 gate.

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
| 1 | Flat-space 10^6-step propagation tolerance | M1 | open |
| 2 | Conserved-quantity drift budget (roadmap states `< 1e-6` relative — confirm per quantity, `CLAUDE.md` §17) | M2A | open |
| 3 | Weak-field deflection benchmark tolerance | M2A | open |
| 4 | **GPU f32 vs. CPU f64 cross-validation tolerance** | M2B | open |
| 5 | Carter constant convention (not a tolerance, but must be fixed before M4A) | M4A | open |
| 6 | Kerr–Schild vs. Boyer–Lindquist exterior agreement tolerance | M4B | open |

`CLAUDE.md` §17: there is no universal numerical-error threshold. Each entry above must
be justified by the relevant numerical method and quantity when it is closed.
