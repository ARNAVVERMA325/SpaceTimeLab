# Spacetime Lab

A small interactive computational-physics laboratory for exploring General Relativity
through mathematically defined spacetime models, numerical computation, and interactive
visualization.

This is a scientific visualization laboratory, not a decorative gravity animation. The
pipeline is deliberately one-directional:

**physical model → metric / field equations → numerical computation → validated result → visualization**

Nothing is rendered before it is validated.

## Project documents

| File | Purpose |
| --- | --- |
| [`CLAUDE.md`](CLAUDE.md) | Project rules: scientific principles, conventions, coordinate and integrator policy, validation requirements, UI honesty rules. |
| [`ROADMAP.md`](ROADMAP.md) | Milestone plan in strict vertical slices, with the GPU port (2B) and Kerr–Schild horizon crossing (4B) as their own hard gates. |
| [`STATUS.md`](STATUS.md) | Per-milestone validated state: which gates are closed, which validation tests pass, and the open numerical-tolerance decisions. |

`STATUS.md` is the file to read first when picking up work mid-project. It records
validated state, not work in progress.

## Current state

**Milestone 1 — Minkowski baseline & engine plumbing: gate passed.**
**Milestone 2A — Exterior Schwarzschild, CPU reference: gate passed.**
**Milestone 3 — Physical observers, frequency shift & thin accretion disk: gate passed.**
**Milestone 4A — Kerr, Boyer–Lindquist exterior: gate passed.**

Four scenes through one pipeline. Flat spacetime bends nothing, so the sky grid arrives
exactly as a pinhole camera projects it. Schwarzschild bends light into a black hole
shadow whose angular radius matches the closed-form prediction
`sin(psi) = b_c sqrt(f) / r` to better than a part in a million. A Novikov–Thorne thin
disk adds an emitting surface: its colour and brightness come from the frequency shift
`g = (k . u_obs) / (k . u_emit)` of the rays that reach it, and the camera can be a
hovering observer or one falling freely from rest at infinity. Kerr adds spin: the shadow
is displaced and flattened on the prograde side, matching Bardeen's analytic critical
curve to within a pixel, and the camera is a zero-angular-momentum observer because a
rotating spacetime admits no diagonal static frame.

M3 and M4A were taken ahead of M2B deliberately — porting the pipeline to shaders before
it had observers in motion, an emission model and a second spacetime would have meant
porting it several times. Milestone 2B (the WebGPU port) remains unblocked and untouched:
the roadmap allows no shader work until the CPU reference is green, and it now is, with
much more in it. M4B (horizon-penetrating Kerr–Schild) is unblocked too. See `STATUS.md`.

## Running it

```bash
npm install
npm run dev        # development server
npm run validate   # typecheck + full validation suite
npm test           # validation suite only
npm run build      # typecheck + production build
```

The validation suite includes a 10^6-step flat-space propagation gate and several full
image renders, so a full run takes roughly twelve seconds.

Rendering runs in a pool of Web Workers, one per core. That is scheduling only: the image
assembled from row blocks is bit-identical to a single-threaded render, which the suite
checks byte for byte on every scene.

## Architecture

The layers of `CLAUDE.md` §20 are directories, not a suggestion. They are kept separate
so the engine does not collapse into one monolithic "physics renderer".

```
src/physics/
  conventions.ts        Signature, units, index and normalization conventions
  core/                 MetricTensor, ChristoffelSymbols, FourVector, PhaseSpaceState
  geometry/             Differential geometry: inversion, index raising, connection
  spacetimes/           Spacetime models — the geometry, and nothing else,
                        plus their independent analytic reference results
  geodesic/             Integrators and the integration driver
  observer/             Tetrad frames, observers, screens, frequency shift
  radiation/            Physical constants, Planck spectrum, CIE 1931 colorimetry
  validation/           Normalization, conserved quantities, tolerances, health
src/visualization/      Ray tracing, orbital-plane reduction, disk emission, display
                        mapping, scene descriptions, canvas output
src/workers/            Render workers for the parallel CPU path
src/ui/                 Provenance and validation reporting
src/data/               Generated colorimetric tables; external datasets (Milestone 5)
```

Two conventions worth knowing before reading the code:

- **Naming follows `CLAUDE.md` §19 verbatim** — `g_mu_nu`, `g_inv_mu_nu`, `christoffel`,
  `four_velocity_u`, `null_wavevector_k`, `energy_E`, `angular_momentum_Lz`.
- **A null tangent is never called a four-velocity.** `PhaseSpaceState.tangent` is
  reached through `null_wavevector_k()` or `four_velocity_u()`, which check the
  worldline kind and throw on a mismatch (`CLAUDE.md` §3, §24).

## Validation

Every tolerance is declared in `src/physics/validation/tolerances.ts` with the quantity
it bounds, whether it is absolute or relative, and why it has the value it does. There
is no universal error threshold (`CLAUDE.md` §17), and the values are set from measured
drift rather than picked round.

Milestone 1 — flat spacetime:

| Check | Result |
| --- | --- |
| Christoffel symbols vanish in Minkowski | Exactly zero, analytically |
| Null normalization `g_mu_nu k^mu k^nu = 0` | ≤ 1e-14 |
| Timelike normalization `g_mu_nu u^mu u^nu = -1` | ≤ 1e-14 relative |
| Straight-line propagation, 10^6 steps | 1.1e-14 relative over a coordinate distance of 1.4e3 |
| `energy_E`, `angular_momentum_Lz` conservation | Exactly zero drift |
| RK4 convergence order | 4th order confirmed under timestep halving |
| Backward-traced image vs. analytic projection | Identical, pixel for pixel |

Milestone 2A — exterior Schwarzschild:

| Check | Result |
| --- | --- |
| Photon sphere locks at `r = 3M` | `\|r - 3M\| = 5.9e-13` over 4.6 orbits, then departs as an unstable orbit must |
| Capture threshold vs. `b_c = 3√3 M` | 1.8e-10 relative, by bisecting traced rays |
| `energy_E`, `angular_momentum_Lz` drift | ≤ 6.3e-12 and ≤ 1.4e-9 — the roadmap asks for 1e-6 |
| Weak-field deflection vs. `4M/b` | Ratio → 1 as 1/b: 1.0302, 1.0030, 1.00029, 1.000029 |
| Traced deflection vs. exact quadrature | 2e-13 to 6.4e-9 relative, from `r_0 = 3.2M` to `10^4 M` |
| Shadow angular radius vs. closed form | < 1e-6 relative |
| Orbital-plane reduction vs. full 3D | Exit directions agree to < 1e-8 |
| Ricci-flatness and `K = 48M²/r⁶` | Confirmed symbolically during derivation |

Milestone 3 — observers, frequency shift and the thin disk:

| Check | Result |
| --- | --- |
| Tetrad orthonormality, static / boosted / free-falling | < 1e-14 |
| Redshift vs. `sqrt(f_emit / f_obs)` on a traced ray | < 1e-12 |
| Static-frame shift vs. covariant `u^t (p_t + Omega p_phi)` | < 1e-12 relative |
| Infaller vs. hoverer at one event vs. `gamma (1 - beta · n)` | < 1e-14 |
| Aberration of the shadow edge under free fall | < 1e-8 relative |
| Flat-space limit: no shift for a co-stationary emitter | Exactly 1 |
| Novikov–Thorne flux vs. independent mpmath evaluation | < 1e-12 relative |
| Disk energy balance `L_inf / Mdot = 1 - sqrt(8/9)` | < 1e-9 |
| Disk temperature vs. independent CODATA 2018 calculation | < 1e-9 relative |
| Blackbody chromaticity vs. `colour-science` at 1 nm | < 1e-4 in CIE (x, y) |
| Parallel render vs. serial | Bit-identical, all scenes |
| Inspector's b = L / E vs. r sin(psi) / sqrt(f) | < 1e-12 |
| Inspector's outcome vs. the b_c capture threshold | Agrees on every pixel sampled |

Milestone 4A — Kerr:

| Check | Result |
| --- | --- |
| Metric, inverse, derivatives, Christoffels vs. independent SymPy | < 1e-12 relative |
| Ricci-flat; Kretschmann vs. the full R_abcd R^abcd contraction | Zero; < 1e-15 relative |
| Schwarzschild limit, component for component | < 1e-14 relative |
| E, L_z conservation | Exact by construction |
| Carter constant drift over a 300M trace | < 1e-9 relative |
| Photon-orbit radius and ISCO vs. their defining equations | < 1e-12, < 1e-11 |
| Photon orbits held over two turns | Radial wander < 1e-8 |
| Exact capture criterion vs. brute force and vs. full integration | Agrees on every ray |
| Extremal shadow spans alpha in [-2M, +7M] edge-on | < 1e-9 |
| Shadow height exactly 2 x 3 sqrt(3) M at every spin | < 1e-6 relative |
| Traced rays inside/outside the critical curve | Captured / escape, three spin-inclination pairs |
| Rendered shadow vs. the analytic extent | Within 2 pixels |

## Four notes on the renders

**In the flat scene, the grid lines curve.** That is rectilinear projection of a sphere,
not light deflection: a pinhole camera maps great circles to straight lines, so meridians
appear straight while parallels do not. In flat spacetime the traced image is identical,
pixel for pixel, to sampling the background along each pixel's initial viewing direction
with no integration at all — which is what the validation suite asserts. Any deviation
from that reference would be a defect, not lensing.

**Click a pixel.** That one ray is traced again on its own, with the same metric,
integrator and tolerance, and reported in full: in Schwarzschild its impact parameter
against `b_c = 3 sqrt(3) M`, in Kerr the two constants `xi = L_z/E` and `eta = Q/E^2` that
separate its motion, where it ended, how far it swept around the hole — rays near the
shadow edge come back having looped more than once — and, for a ray that ends on the disk,
the emission radius, the gas velocity there, the frequency shift split into its static and
orbital factors, and the emitted and observed temperatures. A picture with a panel beside
it only half answers "what am I looking at"; a single line of sight answers it for the
pixel actually in question.

**The disk's colour is a rendering of a computed spectrum, not a measurement.** What the
engine produces is the frequency ratio `g` and the observed blackbody temperature `g T`
at every disk hit; the picture is those numbers put through CIE 1931 colorimetry, an
exposure choice and a tone curve, all disclosed in the panel. How much structure you see
depends on where the visible band falls on the Planck curve: the bolometric boost is
always `g^4`, but in-band luminance goes as `g` deep in the Rayleigh–Jeans tail, so a
40,000 K disk looks nearly uniform while an 8,800 K one shows the approaching side 4.3x
brighter. Same kinematics, different band — measured, not asserted.

**In the black hole scene, the fine bands near the shadow edge are aliased.** Approaching
the capture boundary, the lensing map compresses an unbounded sequence of images of the
whole sky into a vanishing angular width, so no finite ray count resolves it. That is a
sampling limit of the render, not an error in the trajectories, and no smoothing is
applied — anti-aliasing is a rendering operation and must not stand in for resolving the
structure. The dark disc is the *shadow*, which is larger than the horizon and is not a
picture of it: it is the set of directions whose backward-traced rays end on the hole.

## License

MIT — see [`LICENSE`](LICENSE).
