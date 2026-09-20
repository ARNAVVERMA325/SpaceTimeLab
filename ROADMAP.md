# Spacetime Lab — Implementation Roadmap v2

**Companion Document to:** `CLAUDE.md` **Execution Strategy:** Strict Vertical Slices (Pipeline → Physics → Validation → Render for every milestone).

2026-09-20 · @Someone

This revises the original roadmap in three places, all aimed at the same failure mode CLAUDE.md warns against: hiding a validation gap behind a good-looking render.

1. Milestone 2's GPU port is now its own milestone (2B), gated on matching the CPU reference within tolerance, not bundled as a sub-bullet of the Schwarzschild milestone.
2. Milestone 4's Kerr-Schild horizon crossing is split out as its own hard gate (4B), separate from the Boyer-Lindquist exterior work (4A), since it is the single hardest numerical problem in the whole roadmap.
3. Milestone 5 is split into committed scope (3+1 embedding, GWOSC strain) and stretch/exploratory scope (SXS, EHT overlays), since data availability and licensing for the latter aren't fully knowable in advance.

Each milestone still follows the vertical-slice discipline: pipeline, then physics, then validation, then render — never render before validation.

## Milestone 1 — Minkowski Baseline & Engine Plumbing

**Goal:** build the minimal end-to-end pipeline in flat spacetime, establishing architectural patterns and the automated test suite everything later depends on.

**1.1 Data structures & naming.** `MetricTensor`, `ChristoffelSymbols`, `FourVector`, `PhaseSpaceState`. Enforce the standard naming from CLAUDE.md §19 (`g_mu_nu`, `g_inv_mu_nu`, `christoffel`, `four_velocity_u`, `null_wavevector_k`) from day one — renaming later touches every downstream module.

**1.2 Metric.** Minkowski spacetime, g\_mu\_nu = diag(-1, 1, 1, 1). Christoffel symbols identically zero (analytical check, not numerical).

**1.3 Integrator core (CPU).** Fixed-step RK4 and adaptive RKF45 for the first-order geodesic equations. Both are needed later (RK45 for ray tracing, a symplectic option deferred to when long-lived orbits actually require it — CLAUDE.md §7.2 explicitly warns against mandating one integrator everywhere).

**1.4 Validation suite.**

- Photon wavevector normalization: g\_mu\_nu k^mu k^nu = 0
- Massive-particle normalization: g\_mu\_nu u^mu u^nu = -1
- Straight-line propagation in flat space held across 10^6 integration steps

**1.5 Minimal visualizer.** Canvas2D/WebGL canvas rendering backward-traced rays in flat space onto a 2D grid texture — deliberately undramatic; the point is proving the pipeline, not the picture.

**Definition of Done:** unit test suite passes in CI; canvas renders an undistorted grid via backward ray tracing.

## Milestone 2A — Exterior Schwarzschild, CPU Reference

**Goal:** curved spacetime outside the horizon (r > 2M), validated on CPU before any shader work begins. This gate is CPU-only on purpose — see Milestone 2B below.

**2A.1 Schwarzschild metric module.** Analytical g\_mu\_nu and inverse g^mu\_nu in standard spherical coordinates (t, r, theta, phi); analytical Christoffel symbols.

**2A.2 Integrator adaptations.** Termination conditions at r → 2M+ (approaching the horizon) and r > r\_max (asymptotic infinity). Per CLAUDE.md §6.1, the horizon is never detected by a vanishing metric determinant — detect the coordinate breakdown from the actual metric/inverse-metric components in use.

**2A.3 CPU reference raytracer.** Backward-trace null geodesics from a pinhole camera to a background celestial sphere/starfield grid. Compute Einstein-ring deflection angles for axial rays.

**2A.4 Physics validation.**

- Photon-sphere check: unstable circular photon orbit numerically locks at r = 3M
- Conserved quantities: track E = -p\_t and L\_z = p\_phi along trajectories; relative drift under 10^-6
- Deflection-angle benchmark: weak-field deflection against the analytical limit alpha \~ 4GM/(c^2 b)

**Definition of Done:** all 2A.4 tests pass on the CPU reference implementation; this is the ground truth every later GPU or Kerr result gets checked against. No shader work starts before this gate closes.

## Milestone 2B — WebGPU/WGSL Parallelization (separate milestone)

The original plan folded this into M2 as a sub-bullet ("2.5"). It's promoted to its own milestone here: porting a metric evaluation plus an adaptive RKF45 integrator into WGSL compute shaders is a different skillset from the CPU work — no easy NaN inspection, and precision behavior that can differ across GPUs (CLAUDE.md §8 flags this directly: don't assume shader numeric features are portable, and don't claim WASM/GPU precision is automatically "high precision").

**2B.1 Port.** Schwarzschild metric evaluation and the RKF45 integrator into WGSL compute/fragment shaders, targeting 60 FPS at 1080p.

**2B.2 Cross-validation gate.** GPU output checked against the CPU reference from 2A on the same test cases (photon sphere, deflection angle, conserved-quantity drift), not just "looks right" — CLAUDE.md §16 is explicit that visual plausibility isn't validation. Document any tolerance gap between f32 shader precision and the f64 CPU reference rather than silently smoothing over it.

**2B.3 Device checks.** Confirm required shader numeric features on the actual target devices being supported; don't assume portability.

**Definition of Done:** interactive black-hole shadow and Einstein-ring rendering at 60 FPS, AND GPU results matching the CPU reference within a documented, justified tolerance. Shipping the GPU version without this cross-check is exactly the "hide numerical instability behind visual smoothing" failure CLAUDE.md §17/§24 warns against.

## Milestone 3 — Physical Observer, Tetrad Frames & Accretion Disk

**Goal:** replace the naive coordinate camera with a physically rigorous observer tetrad frame (CLAUDE.md §4's "no God camera" for observer-based rendering).

**3.1 Tetrad frame engine.** Orthonormal tetrads e^mu\_(a) for: static observers (u^mu = (1/sqrt(-g\_00), 0, 0, 0)) and freely-falling observers (radial drop geodesics).

**3.2 Relativistic ray generation.** Transform pixel screen coordinates in the observer's local rest frame into global four-momentum k^mu.

**3.3 Frequency shift & intensity.** Emitter-to-observer redshift factor (1+z) = (u^mu k\_mu)\_emitter / (u^mu k\_mu)\_observer; Doppler boosting and relativistic beaming, intensity scaling as (1+z)^-4.

**3.4 Thin accretion disk.** Equatorial geometrically-thin disk between ISCO (r = 6M) and r = 20M; render the asymmetric Doppler bright/dark sides.

**Definition of Done:** toggling static vs. free-falling observer visibly and correctly changes beaming/Doppler shift; disk shows the classical asymmetric intensity pattern.

## Milestone 4 — Kerr Spacetime & Horizon Crossing

**Goal:** frame-dragging around spinning black holes and horizon-penetrating coordinates. Split into 4A and 4B because they are not comparable in difficulty — 4B is likely the hardest single task in the entire roadmap and deserves its own explicit budget and gate rather than sharing a bullet with the metric implementation.

### 4A — Boyer-Lindquist exterior (r > r\_+)

- Implement the Kerr metric in Boyer-Lindquist coordinates for the exterior region.
- Track the Carter constant Q alongside E and L\_z.
- Implement the Kretschmann scalar K = R^abcd R\_abcd for curvature-field visualization.
- Render the asymmetric shadow shape from spin parameter a, and the ergosphere boundary r\_E(theta) = M + sqrt(M^2 - a^2 cos^2 theta).
- Validate: frame-dragging precession rate against the analytical limit; Q, E, L\_z conservation, all outside r\_+.

**Definition of Done (4A):** exterior Kerr results validated and stable; asymmetric shadow benchmarked against Thorne/DNGR reference values.

### 4B — Kerr-Schild horizon-penetrating integration (its own hard gate)

- Implement Kerr-Schild horizon-penetrating coordinates so integration continues through the event horizon without NaN/overflow.
- Cross-validate every trajectory against the 4A Boyer-Lindquist results in the region outside r\_+, where both formulations must agree — this is the check that catches a horizon-penetrating implementation that's subtly wrong but happens to not blow up.
- Budget disproportionate time here relative to its bullet count; treat it as research-adjacent, not a routine port.

**Definition of Done (4B):** smooth camera/particle traversal through r\_+ with no shader crashes or NaNs, AND exterior agreement with 4A confirmed before this is considered done.

## Milestone 5 — Observational Data & Alternative Views

**Goal:** ingest external scientific datasets and expose non-optical geometry views. Split into committed and stretch scope because SXS/EHT data availability and licensing for a browser app can't be fully confirmed until in progress.

### 5A — Committed scope

- **3+1 spatial foliation & embedding view.** Extract the spatial 3-metric gamma\_ij, lapse alpha, shift beta^i. Render the Flamm's-paraboloid 2D isometric embedding into flat 3D space, with strict UI labeling per CLAUDE.md §5.2 (state coordinate system, foliation, time parameter, and that this is an embedding — never imply it's literally what spacetime looks like).
- **GWOSC LIGO strain ingestion.** Parse 1D time-series strain data h(t); drive a linearized transverse-traceless weak-field metric perturbation model. UI must state the model assumptions explicitly (propagation direction, polarization basis, linearized regime, gauge, plane-wave approximation) per CLAUDE.md §12 — label it "linearized-wave visualization driven by measured detector strain," never "the measured spacetime."

**Definition of Done (5A):** user can switch between Optical, 3+1 Embedding, and Curvature heatmap views; GWOSC strain drives real-time grid oscillations under clearly stated assumptions.

### 5B — Stretch / exploratory scope

- **SXS numerical-relativity data.** Ingest precomputed SXS waveform datasets, labeled explicitly as "Numerical Relativity Data," never as a live browser calculation (CLAUDE.md §13).
- **EHT comparison overlay.** Side-by-side modal comparing vacuum Kerr ray-marching against EHT observational reconstructions, framed as a model-to-data comparison with omitted physics (emission, plasma, radiative transfer, instrumental effects) disclosed (CLAUDE.md §14).

Treat 5B as exploratory rather than a hard commitment until data access and licensing are confirmed in practice.

## Per-milestone tracking

A fresh Claude Code session starting mid-milestone has no memory of prior sessions. Keep a `STATUS.md` (or equivalent) in the repo root that records, per milestone/sub-gate: current state (not started / in progress / gate passed), which validation tests from CLAUDE.md §16 are passing, and any open numerical-tolerance decisions (e.g. the 2B GPU/CPU tolerance, once chosen). Update it at the close of each gate, not continuously — it should reflect validated state, not work-in-progress.

This keeps the numerical-validation layer (CLAUDE.md §20) visibly authoritative: a session picking up M4B, for instance, should be able to read that 4A's conservation tests are green before touching horizon-penetrating code, rather than re-deriving that from git history.
