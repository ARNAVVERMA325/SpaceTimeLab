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

Flat spacetime, end to end: the metric and connection, a CPU integrator core, the
validation suite the later milestones are checked against, and a backward null-geodesic
ray tracer rendering an undistorted sky grid.

Milestone 2A (exterior Schwarzschild, CPU reference) is next. See `STATUS.md`.

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
  spacetimes/           Spacetime models — the geometry, and nothing else
  geodesic/             Integrators and the integration driver
  observer/             Tetrad frames, observers, screens
  validation/           Normalization, conserved quantities, tolerances, health
src/visualization/      Ray tracing and canvas output
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

The Milestone 1 gates:

| Check | Result |
| --- | --- |
| Christoffel symbols vanish in Minkowski | Exactly zero, analytically |
| Null normalization `g_mu_nu k^mu k^nu = 0` | ≤ 1e-14 |
| Timelike normalization `g_mu_nu u^mu u^nu = -1` | ≤ 1e-14 relative |
| Straight-line propagation, 10^6 steps | 1.1e-14 relative over a coordinate distance of 1.4e3 |
| `energy_E`, `angular_momentum_Lz` conservation | Exactly zero drift |
| RK4 convergence order | 4th order confirmed under timestep halving |
| Backward-traced image vs. analytic projection | Identical, pixel for pixel |

## A note on the render

The grid lines curve. That is rectilinear projection of a sphere, not light deflection:
a pinhole camera maps great circles to straight lines, so meridians appear straight
while parallels do not. In flat spacetime the traced image is identical, pixel for
pixel, to sampling the background along each pixel's initial viewing direction with no
integration at all — which is what the validation suite asserts. Any deviation from that
reference would be a defect, not lensing.

## License

MIT — see [`LICENSE`](LICENSE).
