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

Two spacetimes through one pipeline. Flat spacetime bends nothing, so the sky grid
arrives exactly as a pinhole camera projects it. Schwarzschild bends light into a black
hole shadow whose angular radius matches the closed-form prediction
`sin(psi) = b_c sqrt(f) / r` to better than a part in a million.

Milestone 2B (the WebGPU port) is unblocked: the roadmap allows no shader work until the
CPU reference is green, and it now is. See `STATUS.md`.

## Running it

```bash
npm install
npm run dev        # development server
npm run validate   # typecheck + full validation suite
npm test           # validation suite only
npm run build      # typecheck + production build
```

The validation suite includes a 10^6-step flat-space propagation gate, so a full run
takes roughly ten seconds.

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
  observer/             Tetrad frames, observers, screens
  validation/           Normalization, conserved quantities, tolerances, health
src/visualization/      Ray tracing, orbital-plane reduction, canvas output
src/ui/                 Provenance and validation reporting
src/data/               External scientific datasets (Milestone 5)
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

## Two notes on the renders

**In the flat scene, the grid lines curve.** That is rectilinear projection of a sphere,
not light deflection: a pinhole camera maps great circles to straight lines, so meridians
appear straight while parallels do not. In flat spacetime the traced image is identical,
pixel for pixel, to sampling the background along each pixel's initial viewing direction
with no integration at all — which is what the validation suite asserts. Any deviation
from that reference would be a defect, not lensing.

**In the black hole scene, the fine bands near the shadow edge are aliased.** Approaching
the capture boundary, the lensing map compresses an unbounded sequence of images of the
whole sky into a vanishing angular width, so no finite ray count resolves it. That is a
sampling limit of the render, not an error in the trajectories, and no smoothing is
applied — anti-aliasing is a rendering operation and must not stand in for resolving the
structure. The dark disc is the *shadow*, which is larger than the horizon and is not a
picture of it: it is the set of directions whose backward-traced rays end on the hole.

## License

MIT — see [`LICENSE`](LICENSE).
