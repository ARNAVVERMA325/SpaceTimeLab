# SPACETIME LAB — CLAUDE CODE PROJECT RULES

## 0. PURPOSE

Spacetime Lab is an AnyaLabs scientific-visualization project for exploring General Relativity through mathematically defined spacetime models, numerical computation, and interactive visualization.

The project must follow this conceptual pipeline:

**Physical model → metric / field equations → numerical computation → validated result → visualization**

This is a scientific visualization laboratory, not a decorative gravity animation and not a full numerical-relativity research code unless a future scope explicitly adds that capability.

The project must never sacrifice physical meaning merely to obtain a visually impressive result.

---

# 1. NON-NEGOTIABLE SCIENTIFIC PRINCIPLES

## 1.1 No invented physics

Never introduce arbitrary equations, deformation factors, fake gravity parameters, visual warps, or unexplained corrections merely to make an effect look better.

If a calculation is too expensive:

* reduce rendering resolution,
* reduce sampling density,
* reduce temporal resolution,
* use a validated approximation,
* restrict the computational domain,
* use precomputed scientific data,
* or disable the feature.

Do **not** change the physical equations simply to improve performance.

Performance optimizations must preserve the declared mathematical model.

---

## 1.2 Do not pretend to solve Einstein's equations when the system only evaluates a known metric

An analytical metric such as Schwarzschild or Kerr is already a specified spacetime solution.

If the application computes geodesics through that metric, it is evaluating consequences of that spacetime.

It is **not** solving the full Einstein field equations dynamically.

The project may eventually support:

1. exact analytical metrics,
2. approximate / perturbative metrics,
3. imported numerical-relativity data,
4. dynamically evolved spacetime fields.

These categories must remain clearly distinguished.

---

## 1.3 Separate mathematical objects from visual representations

The following are not interchangeable:

* spacetime geometry,
* coordinate representation,
* invariant quantities,
* observer measurements,
* visualization mappings.

Whenever a visual element depends on a coordinate choice, slicing, embedding, gauge, or observer frame, the UI should say so.

---

# 2. THE CORE MATHEMATICAL FOUNDATION

The central object is the metric tensor:

$$
ds^2 = g_{\mu\nu} dx^\mu dx^\nu
$$

The engine should conceptually support:

* \(g_{\mu\nu}\)
* \(g^{\mu\nu}\)
* metric derivatives
* Christoffel symbols
* geodesic equations
* proper time
* four-momentum
* four-velocity
* Riemann tensor
* Ricci tensor
* Ricci scalar
* Weyl tensor / Weyl quantities
* curvature invariants
* Killing symmetries
* conserved quantities
* 3+1 quantities where appropriate

Christoffel symbols are defined by

$$
\Gamma^\mu_{\alpha\beta}
=
\frac12 g^{\mu\sigma}
\left(
\partial_\alpha g_{\sigma\beta}
+
\partial_\beta g_{\sigma\alpha}
-
\partial_\sigma g_{\alpha\beta}
\right)
$$

A geodesic may be written as

$$
\frac{d^2x^\mu}{d\lambda^2}
+
\Gamma^\mu_{\alpha\beta}
\frac{dx^\alpha}{d\lambda}
\frac{dx^\beta}{d\lambda}
=0
$$

or preferably, where appropriate, in Hamiltonian first-order form.

The implementation should make the mathematical convention explicit:

* metric signature,
* unit system,
* coordinate ordering,
* parameter used for integration,
* normalization convention,
* index conventions.

Do not silently mix conventions.

---

# 3. GEODESIC REPRESENTATION

For timelike geodesics:

$$
u^\mu = \frac{dx^\mu}{d\tau}
$$

with

$$
g_{\mu\nu}u^\mu u^\nu = -1
$$

for the \((-+++)\) signature.

For null geodesics, use a null tangent / wavevector

$$
k^\mu = \frac{dx^\mu}{d\lambda}
$$

with

$$
g_{\mu\nu}k^\mu k^\nu = 0
$$

Do not call \(k^\mu\) a four-velocity.

The engine must distinguish:

* massive test-particle worldlines,
* photon / null-ray trajectories,
* observer worldlines.

---

# 4. OBSERVERS AND CAMERA MODEL

An observer-based camera should be treated as a physical observer, not a fictional Euclidean "god camera."

An observer has:

* an event \(x^\mu\),
* a timelike four-velocity \(u^\mu\),
* a local orthonormal spatial frame,
* a local screen / viewing basis.

A tetrad or equivalent orthonormal observer frame should be used when computing what an observer actually sees.

For observer-based rendering, quantities such as:

* photon direction,
* frequency,
* redshift,
* aberration,
* observed angular position,
* intensity transformations

must be defined relative to the observer frame.

However, not every application view requires a physical observer.

A purely mathematical geometry-inspection mode may use an abstract visualization camera. In that case it must be clearly labeled as a visualization coordinate system rather than a physical observer.

Therefore:

**"No God Camera" applies to observer-based physical rendering, not to every UI camera in the application.**

---

# 5. SPACETIME VISUALIZATION MODES

The project should not force all aspects of GR into one literal "spacetime mesh."

Use multiple linked representations.

## 5.1 Observer / Optical View

Backward-trace null geodesics from an observer's local screen into the spacetime.

Useful for:

* gravitational lensing,
* black-hole shadows,
* accretion-disk appearance,
* redshift / blueshift,
* multiple images,
* caustic structure,
* visual appearance near compact objects.

The null-geodesic structure is geometric, but the observed image depends on the observer's worldline and local tetrad.

Do not describe the entire rendered mesh or screen representation as "gauge invariant."

---

## 5.2 3+1 Spatial View

Represent a spacetime as spatial hypersurfaces

$$
\Sigma_t
$$

with spatial metric

$$
\gamma_{ij}
$$

and, where appropriate,

* lapse \(\alpha\),
* shift \(\beta^i\),
* extrinsic curvature \(K_{ij}\).

The spatial slice and its coordinate components depend on the chosen foliation and coordinates.

Therefore a 3D grid visualization must explicitly state:

* coordinate system,
* foliation,
* time parameter,
* whether the geometry shown is coordinate geometry, proper-distance geometry, or an embedding.

Never imply that a Euclidean 3D deformation is literally "what spacetime looks like."

For Schwarzschild, a 2D equatorial isometric embedding such as the Flamm paraboloid may be shown as an educational visualization, but it must be labeled as an embedding of a spatial slice into an auxiliary Euclidean space.

---

## 5.3 Curvature / Invariant View

Use genuine scalar quantities where the goal is coordinate-independent curvature information.

Examples include:

* Ricci scalar \(R\),
* Kretschmann scalar

  $$
  K=R_{\alpha\beta\gamma\delta}R^{\alpha\beta\gamma\delta}
  $$
* contractions of the Weyl tensor,
* other explicitly defined curvature invariants.

Newman–Penrose Weyl scalars such as \(\Psi_0,\ldots,\Psi_4\) are useful but are tetrad-dependent quantities and must not be presented as ordinary coordinate-free scalar invariants.

If Weyl scalars are displayed, the chosen tetrad / frame must be documented.

---

## 5.4 Causal / Spacetime-Diagram View

Support genuinely 4D representations through linked views such as:

* spacetime diagrams,
* worldlines,
* null trajectories,
* light cones,
* time slices,
* synchronized selections between a 3D view and a spacetime diagram.

"4D" means **3 spatial dimensions + time**.

Never claim to display four spatial dimensions literally.

---

# 6. COORDINATE-SYSTEM POLICY

Coordinates are representations, not physical objects.

The engine must keep coordinate-system assumptions explicit.

Potential coordinate systems include:

* Minkowski Cartesian coordinates,
* Schwarzschild coordinates,
* ingoing / outgoing Eddington–Finkelstein coordinates,
* Painlevé–Gullstrand coordinates,
* Kerr Boyer–Lindquist coordinates,
* horizon-penetrating Kerr coordinates,
* Kerr–Schild coordinates,
* FLRW coordinates,
* other explicitly documented systems.

## 6.1 Schwarzschild horizon

Standard Schwarzschild coordinates have a coordinate singularity at

$$
r=2M
$$

in geometric units.

This is a singularity of the coordinate representation, not a physical curvature singularity.

The Schwarzschild metric determinant in standard spherical coordinates is proportional to

$$
-r^4\sin^2\theta
$$

and does not vanish merely because \(r=2M\).

Therefore code must not detect the horizon by checking for a vanishing metric determinant.

Instead detect the coordinate breakdown from the relevant metric components / inverse metric / equations being used.

---

## 6.2 Horizon-penetrating coordinates

Horizon-penetrating coordinates are appropriate when the numerical trajectory or field evolution must continue through the event horizon.

For Schwarzschild, simpler choices such as ingoing Eddington–Finkelstein or Painlevé–Gullstrand may be preferable depending on the task.

For Kerr, horizon-penetrating Kerr / Kerr–Schild-type coordinates are strong candidates for calculations that cross the horizon.

There is no rule that every calculation must use Kerr–Schild coordinates.

Selection must depend on:

* domain,
* numerical stability,
* analytical simplicity,
* desired outputs,
* implementation complexity.

Do not introduce a coordinate transformation unless its equations are verified.

---

## 6.3 Exterior-only calculations

If a feature only traces rays outside the horizon, a horizon-penetrating system may be unnecessary.

For example, a black-hole shadow renderer can terminate rays that cross the event horizon rather than integrating through the interior.

Use the simplest validated formulation appropriate to the actual task.

---

# 7. NUMERICAL INTEGRATION POLICY

There is no universal "correct integrator."

The integrator must be selected according to the mathematical formulation and task.

## 7.1 Adaptive Runge–Kutta

Adaptive RK methods such as RK45 / RKF45 are valid candidates for:

* finite geodesic integrations,
* rapidly varying trajectories,
* calculations requiring local error control,
* ray tracing where different rays encounter very different curvature scales.

They are particularly attractive when implementation simplicity and local error estimation are important.

---

## 7.2 Symplectic / geometric integration

Symplectic methods are valuable for long-term integration of Hamiltonian systems because they control qualitative long-time error behavior and can avoid secular energy drift.

They are especially worth considering for:

* long-lived bound timelike orbits,
* Hamiltonian geodesic benchmarks,
* comparisons of numerical integrators,
* long-duration conservative dynamics.

However:

**Do not require a symplectic integrator for every timelike geodesic.**

Some symplectic methods require a suitable Hamiltonian splitting, and Kerr spacetime can require specialized constructions.

"Velocity Verlet everywhere" is not an acceptable general GR rule.

"Implicit midpoint everywhere" is also not automatically appropriate because of its computational cost and nonlinear solve requirements.

The architecture should therefore allow multiple validated integrators.

---

## 7.3 Recommended conceptual strategy

Implement an integrator abstraction so that the physics engine can support, for example:

* RK4,
* adaptive RK45 / RKF45,
* an appropriate symplectic / variational method for long-term Hamiltonian tests.

Choose the method per simulation mode.

When two methods are available, benchmark them against analytical invariants and known solutions.

---

# 8. GPU / CPU / WASM POLICY

GPU acceleration is a performance mechanism, not a replacement for scientific validation.

## CPU / WASM

Suitable for:

* reference implementations,
* high-precision calculations,
* validation,
* selected particle trajectories,
* metric / tensor evaluation,
* long-lived orbits,
* preprocessing,
* scientific comparisons.

Use floating-point precision deliberately.

Do not claim that WASM automatically provides "high precision." It normally provides IEEE-style double precision (`f64`) when the implementation and language path support it; that is higher precision than `f32`, but it is still finite precision.

---

## WebGPU / WGSL

Suitable for:

* massive numbers of independent rays,
* parallel geodesic calculations,
* image synthesis,
* field sampling,
* grid visualization,
* particle batches,
* post-processing.

Do not assume arbitrary shader numeric features are universally portable.

Shader precision and supported features must be checked on the target device.

Maintain a CPU/reference implementation for physics-sensitive validation.

The renderer may use GPU approximations that are visually efficient only when they preserve the declared physical model.

---

# 9. GEODESIC RENDERING

Use the term **ray tracing / null-geodesic integration** rather than treating everything as generic "ray marching."

A physically grounded renderer should conceptually perform:

1. observer-frame ray generation,
2. conversion into spacetime phase-space initial conditions,
3. null-geodesic integration,
4. interaction / intersection tests with emitting surfaces or fields,
5. frequency / redshift calculations where applicable,
6. radiometric or color processing,
7. image accumulation / anti-aliasing.

Temporal anti-aliasing, denoising, interpolation, and similar graphics techniques are rendering operations.

They must not alter the underlying physical trajectory.

---

# 10. SPACETIME MODELS

Models must be explicitly classified.

## Exact / analytical

Potential examples:

* Minkowski spacetime,
* Schwarzschild spacetime,
* Kerr spacetime,
* specified FLRW cosmology.

## Approximate / perturbative

Examples:

* linearized gravitational-wave spacetime,
* weak-field perturbations.

These must never be described as exact solutions of the full nonlinear Einstein equations.

## Theoretical metric models

Examples:

* Morris–Thorne-type wormhole metrics.

These must be labeled as specified theoretical metric models / ansätze, with assumptions clearly stated.

Do not imply observational confirmation.

---

# 11. GENERAL RELATIVITY DATA HIERARCHY

Every result shown by the application should conceptually fall into one of these categories:

### A. Exact analytical

Derived from a mathematically defined exact solution.

### B. Numerical approximation

Obtained by numerically integrating or evaluating equations.

### C. Imported numerical-relativity result

Obtained from an external scientific simulation such as SXS.

### D. Observational data

Measured / calibrated astronomical or detector data.

### E. Derived comparison

A result produced by comparing a theoretical model with observational or numerical data.

The UI must never silently mix these categories.

---

# 12. GRAVITATIONAL-WAVE DATA

GWOSC provides public detector strain data \(h(t)\).

Detector strain is an observational data product obtained through detector calibration; it is not itself a complete four-dimensional metric field.

Therefore:

**Do not take one detector's \(h(t)\) and present it as a uniquely reconstructed global spacetime.**

A browser demonstration may use measured strain to drive a linearized gravitational-wave model, but only under explicitly stated assumptions such as:

* chosen propagation direction,
* polarization basis,
* linearized regime,
* gauge,
* spatial approximation / plane-wave assumption.

The UI should say something like:

> "Linearized-wave visualization driven by measured detector strain under the stated model assumptions."

Not:

> "This is the measured spacetime."

---

# 13. SXS NUMERICAL-RELATIVITY DATA

SXS provides publicly available numerical-relativity simulation products, including gravitational-wave waveforms and related metadata.

Imported SXS data must be labeled as:

**Numerical Relativity Data**

not as a live calculation performed by the browser.

The project should store or expose:

* simulation identifier,
* physical parameters,
* provenance,
* data product type,
* units,
* resolution information where available,
* source citation.

Do not imply that an imported waveform is equivalent to importing the complete dynamically evolved metric unless that exact data product is actually available and supported.

---

# 14. EVENT HORIZON TELESCOPE DATA

EHT provides public observational data products, including calibrated data from observing campaigns.

Use a strict distinction between:

* observed / calibrated interferometric data,
* reconstructed images,
* model-derived quantities,
* theoretical simulations.

A vacuum Kerr rendering is not automatically a prediction of the observed EHT image because actual observations involve emission physics, plasma, radiative transfer, scattering, instrumental effects, and image-reconstruction methodology.

An EHT comparison should therefore be presented as a **model-to-data comparison**, with all omitted physical effects disclosed.

---

# 15. INTERSTELLAR / DNGR INSPIRATION

The project's Interstellar connection is methodological, not a requirement to reproduce the movie.

DNGR was developed by Oliver James, Eugénie von Tunzelmann, Paul Franklin and Kip Thorne for rendering observations near Kerr black holes and wormholes using general-relativistic light propagation.

Relevant references include:

* James, von Tunzelmann, Franklin & Thorne,
  **"Gravitational lensing by spinning black holes in astrophysics, and in the movie Interstellar"**
  Classical and Quantum Gravity 32, 065001 (2015).
  DOI: 10.1088/0264-9381/32/6/065001

* James, von Tunzelmann, Franklin & Thorne,
  **"Visualizing Interstellar's Wormhole"**
  American Journal of Physics 83, 486–499 (2015).
  DOI: 10.1119/1.4916949

DNGR used ray-bundle methods for smooth high-quality rendering. Spacetime Lab does not need to reproduce DNGR internally.

A practical browser implementation may use individual-ray geodesic tracing plus sampling / accumulation methods.

Never claim to reproduce DNGR unless the implementation actually does so.

---

# 16. VALIDATION IS PART OF THE PHYSICS ENGINE

Every important numerical subsystem must have tests.

## Required validation categories

### Flat-space limit

Recover Minkowski behavior when curvature parameters vanish or the relevant limit is taken.

### Geodesic normalization

Timelike:

$$
g_{\mu\nu}u^\mu u^\nu \approx -1
$$

Null:

$$
g_{\mu\nu}k^\mu k^\nu \approx 0
$$

### Conserved quantities

Use only when the metric possesses the corresponding symmetry.

For stationary spacetimes:

$$
E=-p_t
$$

For axisymmetric spacetimes:

$$
L_z=p_\phi
$$

In Kerr, track the Carter constant using one clearly documented convention.

### Known analytical behavior

Examples include:

* Schwarzschild photon sphere,
* Schwarzschild ISCO,
* weak-field limit,
* Kerr frame-dragging behavior,
* flat-space propagation,
* analytically known redshift relations.

### Convergence tests

Changing:

* timestep,
* spatial resolution,
* ray resolution,
* numerical tolerance

should produce predictable convergence behavior where expected.

A result is not considered validated merely because it "looks right."

---

# 17. ERROR HANDLING

Do not use one universal numerical-error threshold for every quantity.

Instead use:

* absolute tolerance,
* relative tolerance,
* model-specific tolerances,
* benchmark-specific tolerances,
* accumulated invariant drift,
* convergence behavior.

The UI may display numerical-health indicators, but thresholds must be justified by the relevant numerical method and quantity.

A warning such as

> "Numerical error exceeds validated tolerance"

is preferable to pretending that `1e-5` is universally meaningful.

NaN / infinity detection must be explicit.

When numerical failure occurs, expose useful information where possible:

* integrator,
* model,
* coordinate system,
* current parameter,
* step size,
* failed quantity.

Never silently replace a failed physical calculation with a visually plausible fake result.

---

# 18. TIME DILATION / REDSHIFT LANGUAGE

Avoid oversimplified wording such as:

> "Time slows down near a black hole."

Use the physically specific statement appropriate to the situation.

For example, in Schwarzschild coordinates for a stationary observer at fixed radius outside the horizon,

$$
d\tau = \sqrt{1-\frac{2M}{r}}\,dt
$$

with \(t\) the Schwarzschild coordinate time normalized to an observer at infinity.

This statement is specific to the chosen observer and coordinate system.

Do not turn the lapse \(\alpha\) into a universal synonym for "how fast time passes."

Lapse is a quantity associated with a chosen 3+1 slicing.

Similarly, do not describe coordinate-time behavior of a freely falling observer using the same language as a stationary observer.

---

# 19. CODE NAMING

Use standard, readable mathematical naming where practical.

Examples:

* `g_mu_nu`
* `g_inv_mu_nu`
* `christoffel`
* `riemann`
* `ricci`
* `weyl`
* `kretschmann`
* `lapse_alpha`
* `shift_beta`
* `spatial_metric_gamma`
* `four_velocity_u`
* `null_wavevector_k`
* `four_momentum_p`
* `energy_E`
* `angular_momentum_Lz`
* `carter_Q`

Avoid meaningless names such as:

* `gravityStrength`
* `spaceBend`
* `warpAmount`
* `magicFactor`

unless they are clearly part of a purely visual effect and not part of the physics engine.

---

# 20. ARCHITECTURAL PRINCIPLE

The physics engine should conceptually separate:

### Metric / spacetime model

Defines the geometry.

### Differential-geometry layer

Computes metric-related quantities.

### Geodesic / dynamical layer

Integrates particle and photon trajectories.

### Observer layer

Defines physical observers, tetrads, screens, and measurements.

### Numerical-validation layer

Checks normalization, invariants, convergence, and reference solutions.

### Visualization layer

Turns validated physical results into visual representations.

### Data layer

Handles scientific external datasets and provenance.

### UI / educational layer

Explains what the user is seeing without confusing coordinates, observables, invariants, and simulated quantities.

These concerns should not become one monolithic "physics renderer."

---

# 21. PERFORMANCE POLICY

The project may target interactive frame rates, including approximately 60 FPS where practical.

However:

**60 FPS is a rendering goal, not a scientific-validity requirement.**

When the physics calculation cannot meet the rendering budget, use techniques such as:

* cached results,
* asynchronous computation,
* lower ray count,
* temporal accumulation,
* precomputed trajectories,
* precomputed scientific data,
* reduced field resolution,
* progressive refinement.

Do not alter the underlying metric or physical equations merely to preserve frame rate.

---

# 22. SCIENTIFIC HONESTY IN THE UI

Avoid:

> "Gravity pulls the light."

Prefer:

> "The photon follows a null geodesic of the chosen spacetime."

Avoid:

> "The fabric of space bends downward."

Prefer:

> "This is a visualization / embedding of a chosen spatial slice."

Avoid:

> "This is what a black hole really looks like."

Prefer:

> "Computed appearance for the selected spacetime, observer, emission model, and rendering assumptions."

Avoid:

> "Real spacetime reconstructed from LIGO."

Prefer:

> "Linearized-wave visualization driven by detector strain under the selected assumptions."

Every major visualization should make it possible for the user to inspect:

* model,
* metric,
* coordinates,
* units,
* observer,
* numerical method,
* approximation level,
* external-data provenance,
* validation status.

---

# 23. AI AGENT BEHAVIOR

Claude Code must behave as an engineering agent working inside a scientific project.

Before making major changes:

1. inspect the repository,
2. identify existing architecture,
3. identify existing numerical conventions,
4. read relevant documentation,
5. preserve working validated components,
6. avoid unnecessary rewrites.

When implementing physics:

* derive or verify the equation first,
* identify the coordinate system,
* identify the assumptions,
* identify the integration parameter,
* identify relevant invariants,
* write a validation test,
* then optimize.

If uncertain about a physics claim, do not invent an answer.

Instead:

* check project references,
* consult authoritative literature,
* document the uncertainty,
* choose the simplest defensible implementation.

If a requested feature would require full numerical relativity, radiation hydrodynamics, GRMHD, or another major research-scale subsystem, do not silently approximate it while claiming scientific equivalence.

State the approximation explicitly.

---

# 24. WHAT CLAUDE MUST NOT DO

Do not:

* invent equations,
* invent physical constants,
* change metric signatures silently,
* mix coordinate systems without transformation,
* call coordinate artifacts physical curvature,
* call tetrad-dependent quantities invariants,
* call null vectors four-velocities,
* claim observational data is a unique spacetime reconstruction,
* claim precomputed numerical-relativity data is live browser simulation,
* hide numerical instability behind visual smoothing,
* add unexplained fudge factors,
* claim a visualization is physically meaningful when it is only illustrative,
* prioritize visual resemblance over mathematical correctness.

---

# 25. CURRENT PROJECT BOUNDARIES

The project may eventually contain:

* Minkowski spacetime,
* Schwarzschild black holes,
* Kerr black holes,
* FLRW cosmology,
* linearized gravitational waves,
* theoretical wormhole metrics,
* null geodesic rendering,
* timelike geodesics,
* observer/tetrad rendering,
* curvature visualization,
* spacetime diagrams,
* 3+1 foliation views,
* gravitational lensing,
* black-hole shadows,
* redshift / frequency-shift calculations,
* accretion-disk visualization,
* GWOSC observational data,
* SXS numerical-relativity data,
* EHT observational comparisons.

These are potential capabilities, not a requirement that all must be implemented simultaneously.

Full 3D dynamical numerical relativity is outside the default scope unless explicitly introduced as a separate major subsystem.

---

# 26. FINAL DESIGN PHILOSOPHY

Spacetime Lab should answer two different questions without confusing them:

### Mathematical question

> "What does this spacetime geometry and its geodesics do?"

### Observational question

> "What would this particular observer measure or see?"

The first uses geometry, tensors, coordinates, curvature, and numerical dynamics.

The second uses an observer, tetrad, photon propagation, emission assumptions, and measurement models.

The strongest version of Spacetime Lab connects both while keeping their assumptions visible.

The project should feel less like a "black-hole animation" and more like a **small interactive computational-physics laboratory**.

When in doubt:

**Prefer the scientifically explicit answer over the visually convenient one.**
