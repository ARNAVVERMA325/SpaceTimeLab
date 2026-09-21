/**
 * Numerical-validation layer (CLAUDE.md §20): named, individually justified tolerances.
 *
 * CLAUDE.md §17 forbids one universal numerical-error threshold and says a warning of
 * the form "numerical error exceeds validated tolerance" is preferable to pretending
 * that 1e-5 is universally meaningful. Every tolerance the test suite and the UI use is
 * therefore declared here with its kind, the quantity it applies to, and the reasoning
 * behind its value.
 *
 * The values below were set by measuring the actual drift of the M1 implementation and
 * leaving roughly an order of magnitude of headroom, so a regression trips the gate
 * while ordinary run-to-run variation does not.
 */

export interface Tolerance {
  readonly id: string;
  /** The quantity being bounded, in CLAUDE.md §19 naming where one applies. */
  readonly quantity: string;
  readonly kind: 'absolute' | 'relative';
  readonly value: number;
  /** Why this number, for this quantity, under this numerical method. */
  readonly justification: string;
}

function tolerance(t: Tolerance): Tolerance {
  return Object.freeze(t);
}

/**
 * g_mu_nu k^mu k^nu for a null wavevector, in Minkowski, at a single event.
 *
 * Absolute, because the exact target is 0 and a relative bound is undefined there. The
 * contraction is a sum of four products of O(1) components, so a handful of binary64
 * roundings bound the residual at a few times machine epsilon.
 */
export const NULL_NORMALIZATION_POINTWISE = tolerance({
  id: 'null-normalization-pointwise',
  quantity: 'g_mu_nu k^mu k^nu',
  kind: 'absolute',
  value: 1e-14,
  justification:
    'Exact target is 0, so the bound is absolute. Evaluating the contraction costs a ' +
    'few binary64 operations on O(1) components; 1e-14 is roughly 45 machine epsilons ' +
    'and leaves headroom over the observed residual without hiding a real defect.',
});

/**
 * g_mu_nu u^mu u^nu for a four-velocity, in Minkowski, at a single event.
 *
 * Relative, because the target is -1 and a relative bound is well defined.
 */
export const TIMELIKE_NORMALIZATION_POINTWISE = tolerance({
  id: 'timelike-normalization-pointwise',
  quantity: 'g_mu_nu u^mu u^nu',
  kind: 'relative',
  value: 1e-14,
  justification:
    'Target is -1, so a relative bound applies. Same operation count as the null case; ' +
    'for a boosted u^mu the cancellation between -(u^0)^2 and |u_spatial|^2 grows with ' +
    'the Lorentz factor, so this bound holds for the moderate boosts the M1 suite uses.',
});

/**
 * Drift of the normalization invariant accumulated over a long integration.
 *
 * The invariant is not enforced by the integrator; it is a diagnostic that drifts as
 * local truncation and rounding errors accumulate. In flat space the connection
 * vanishes, so the tangent is advanced by exactly zero derivative and drift is pure
 * floating-point noise in the contraction rather than in the state.
 */
export const NORMALIZATION_DRIFT_LONG_RUN = tolerance({
  id: 'normalization-drift-long-run',
  quantity: 'g_mu_nu t^mu t^nu over 10^6 steps',
  kind: 'absolute',
  value: 1e-12,
  justification:
    'In Minkowski the geodesic right-hand side returns dt^mu/dparam = 0 exactly, so the ' +
    'tangent components come back bit-for-bit identical after 10^6 steps and the ' +
    'invariant accumulates no drift at all. The residual that remains is the rounding of ' +
    'the contraction itself, measured at 3.3e-16 for a ray with 1/sqrt(3) components. ' +
    'The bound is left at 1e-12 so the same gate can be reused against a curved model, ' +
    'where the tangent genuinely does evolve, without being rewritten.',
});

/**
 * Deviation from analytically straight propagation in flat space (ROADMAP.md 1.4).
 *
 * Relative to the distance travelled: absolute position error necessarily grows with
 * the path length, so a fixed absolute bound would be a statement about run length
 * rather than about the integrator.
 */
export const FLAT_PROPAGATION_STRAIGHTNESS = tolerance({
  id: 'flat-propagation-straightness',
  quantity: 'x^mu(lambda) vs. x^mu(0) + k^mu lambda',
  kind: 'relative',
  value: 1e-12,
  justification:
    'Measured against the coordinate distance travelled, since absolute position error ' +
    'necessarily scales with path length and a fixed absolute bound would describe the ' +
    'run length rather than the integrator. Summing 10^6 equal increments accumulates at ' +
    'worst O(N) roundings, but because the increments are equal and the partial sums grow ' +
    'monotonically the realised error is far smaller: 1.1e-14 relative over a coordinate ' +
    'distance of 1.4e3. The bound of 1e-12 sits about two orders of magnitude above that, ' +
    'which catches a regression without tripping on run-to-run variation.',
});

/**
 * Relative drift of a Killing conserved quantity along a geodesic.
 *
 * ROADMAP.md 2A.4 sets 1e-6 for Schwarzschild. Minkowski is stricter: the tangent is
 * exactly constant, so E and L_z are limited only by the arithmetic in the contraction.
 */
export const CONSERVED_QUANTITY_DRIFT_FLAT = tolerance({
  id: 'conserved-quantity-drift-flat',
  quantity: 'energy_E, angular_momentum_Lz in Minkowski',
  kind: 'relative',
  value: 1e-12,
  justification:
    'Flat space only. The roadmap budget of 1e-6 applies to Schwarzschild, where the ' +
    'connection is non-zero and drift is genuine. Holding Minkowski to 1e-6 would let a ' +
    'real defect pass, so the flat-space gate is tightened to what the arithmetic ' +
    'actually supports: measured drift over 10^6 steps is exactly zero for every Killing ' +
    'quantity, because the tangent never changes. Note that angular_momentum_Lz is built ' +
    'from a position-dependent Killing field, so its contraction involves coordinates ' +
    'that grow along the ray, and its initial value is zero for any ray through the ' +
    'spatial origin -- where a relative measure is undefined.',
});

/**
 * Agreement between analytical Christoffel symbols and centrally differenced ones.
 *
 * Absolute, because the target components are 0 in Minkowski.
 */
export const CHRISTOFFEL_NUMERIC_VS_ANALYTIC = tolerance({
  id: 'christoffel-numeric-vs-analytic',
  quantity: 'Gamma^mu_{alpha beta}',
  kind: 'absolute',
  value: 1e-10,
  justification:
    'Central differencing balances O(h^2) truncation against O(eps/h) roundoff at ' +
    'h ~ cbrt(eps), leaving a floor near 4e-11 for O(1) metric components. 1e-10 sits ' +
    'just above that floor. This bound applies to the numerical path only; the ' +
    'analytical Minkowski symbols are exactly zero and are asserted as such.',
});

/** Largest |g^{mu sigma} g_{sigma nu} - delta^mu_nu| after a 4x4 inversion. */
export const METRIC_INVERSE_RESIDUAL = tolerance({
  id: 'metric-inverse-residual',
  quantity: 'g^{mu sigma} g_{sigma nu} - delta^mu_nu',
  kind: 'absolute',
  value: 1e-12,
  justification:
    'Gauss-Jordan with partial pivoting on a 4x4 system performs O(n^3) ~ 64 binary64 ' +
    'operations per entry path. For the well-conditioned metrics in scope the residual ' +
    'stays within a few hundred epsilons; 1e-12 flags a genuinely ill-conditioned chart ' +
    'without tripping on ordinary rounding.',
});

/** Largest |g_mu_nu e^mu_(a) e^nu_(b) - eta_ab| for a tetrad. */
export const TETRAD_ORTHONORMALITY = tolerance({
  id: 'tetrad-orthonormality',
  quantity: 'g_mu_nu e^mu_(a) e^nu_(b) - eta_ab',
  kind: 'absolute',
  value: 1e-13,
  justification:
    'Absolute, since the target entries are 0 and +/-1. Each entry is a 16-term ' +
    'contraction of O(1) quantities. The Minkowski tetrad used in M1 is the identity ' +
    'frame and satisfies this exactly; the bound is set for the general construction ' +
    'that M3 introduces.',
});

// ---------------------------------------------------------------------------
// Milestone 2A — exterior Schwarzschild, CPU reference.
//
// Every bound below was set by measuring the implementation and leaving one to two
// orders of magnitude of headroom, so a regression trips the gate while ordinary
// variation does not.
// ---------------------------------------------------------------------------

/**
 * Agreement between analytical Schwarzschild Christoffel symbols and centrally
 * differenced ones, as a fraction of the largest symbol at that event.
 */
export const SCHWARZSCHILD_CHRISTOFFEL_NUMERIC = tolerance({
  id: 'schwarzschild-christoffel-numeric',
  quantity: 'max |Gamma_numeric - Gamma_analytic| / max |Gamma_analytic|',
  kind: 'relative',
  value: 1e-7,
  justification:
    'Measured agreement is 3e-10 or better between r = 3M and r = 100M, degrading to ' +
    '1.6e-8 at r = 2.1M, where f = 1 - 2M/r is small and the connection varies sharply ' +
    'over the differencing step, and to 6e-9 at r = 1000M, where the coordinate-scaled ' +
    'step grows with r. The bound covers the worst of those with headroom. It applies to ' +
    'the numerical path only: the analytical symbols are exact and are compared against ' +
    'closed-form expressions directly.',
});

/**
 * Relative drift of E and L_z along a traced Schwarzschild null geodesic.
 *
 * ROADMAP.md 2A.4 asks for drift under 1e-6. This bound is deliberately a hundred times
 * tighter, and the roadmap figure is asserted separately so the stated gate is visibly
 * met.
 */
export const CONSERVED_QUANTITY_DRIFT_SCHWARZSCHILD = tolerance({
  id: 'conserved-quantity-drift-schwarzschild',
  quantity: 'energy_E, angular_momentum_Lz along a traced null geodesic',
  kind: 'relative',
  value: 1e-7,
  justification:
    'Measured drift over a full deflection trace is at most 6e-12 for energy_E and ' +
    '1.4e-9 for angular_momentum_Lz, across turning points from 3.2M to 10^4 M. L_z ' +
    'drifts more because it is built from r^2 sin^2(theta) k^phi, so it carries the ' +
    'coordinate magnitudes with it. Adopting the roadmap figure of 1e-6 as the gate ' +
    'would leave three orders of magnitude of slack and let a real defect pass, which is ' +
    'the failure mode CLAUDE.md §17 warns against.',
});

/**
 * Departure of a circular null orbit from r = 3M over a bounded affine parameter.
 *
 * Bounded on purpose. The orbit is unstable, so this cannot be a statement about
 * arbitrarily long integrations: any perturbation, including the rounding of the initial
 * conditions themselves, grows exponentially.
 */
export const PHOTON_SPHERE_LOCK = tolerance({
  id: 'photon-sphere-lock',
  quantity: '|r - 3M| on a circular null orbit, over affine parameter <= 50M',
  kind: 'absolute',
  value: 1e-11,
  justification:
    'Absolute, since the target is the fixed radius 3M. Measured departure is 5.9e-13 at ' +
    'lambda = 50M, about 4.6 complete orbits. Beyond that the orbit genuinely leaves: ' +
    'the instability e-folds roughly every 1.8 in affine parameter, so machine-level ' +
    'rounding in the initial data reaches order unity near lambda = 65M and the orbit ' +
    'has visibly departed by lambda = 100M. That is physics, not a numerical defect, and ' +
    'the suite asserts the departure as well as the lock.',
});

/** Traced deflection angle against the independently computed exact quadrature. */
export const DEFLECTION_VS_EXACT = tolerance({
  id: 'deflection-vs-exact',
  quantity: 'traced deflection angle vs. the exact orbit-equation quadrature',
  kind: 'relative',
  value: 1e-7,
  justification:
    'Measured agreement is between 2e-13 and 6.4e-9 for turning points from 3.2M to ' +
    '10^4 M, worst in the strong field where the deflection exceeds pi. The comparison ' +
    'is between two genuinely different methods -- geodesic ODE integration against a ' +
    'Gauss-Legendre quadrature of the orbit equation -- so agreement at this level is ' +
    'evidence about the integrator rather than a restatement of it.',
});

/**
 * The capture threshold in impact parameter, against b_c = 3 sqrt(3) M.
 *
 * Found by bisecting on whether a traced ray is captured, so it exercises the metric,
 * the connection, the integrator and the capture condition together.
 */
export const CRITICAL_IMPACT_PARAMETER = tolerance({
  id: 'critical-impact-parameter',
  quantity: 'bisected capture threshold vs. b_c = 3 sqrt(3) M',
  kind: 'relative',
  value: 1e-8,
  justification:
    'Measured 1.8e-10 after 43 bisection steps. The threshold is a property of the whole ' +
    'pipeline rather than of any one component, which makes it the strongest single ' +
    'check in the Milestone 2A suite.',
});

/** |g_mu_nu k^mu k^nu| at the end of a traced Schwarzschild null geodesic. */
export const NULL_NORMALIZATION_TRACED = tolerance({
  id: 'null-normalization-traced',
  quantity: 'g_mu_nu k^mu k^nu after a full curved-spacetime trace',
  kind: 'absolute',
  value: 1e-9,
  justification:
    'Absolute, since the target is exactly 0. Measured at most 1.8e-11 over deflection ' +
    'traces reaching r = 3.2M. Looser than the pointwise flat-space bound of 1e-14 ' +
    'because the invariant is not enforced by the integrator and genuinely accumulates ' +
    'error here, where the connection is non-zero -- the two are different quantities ' +
    'and CLAUDE.md §17 forbids giving them one shared threshold.',
});

/**
 * |g_mu_nu k^mu k^nu| for the interactive CPU preview render.
 *
 * Deliberately looser than `NULL_NORMALIZATION_TRACED`, and a separate named tolerance
 * rather than a relaxation of that one. CLAUDE.md §17 calls for benchmark-specific
 * tolerances precisely so that a fast preview and a validated reference result are not
 * judged by the same number — reporting a preview as "degraded" against the reference
 * gate would be as misleading as reporting it as validated.
 *
 * CLAUDE.md §21 governs the trade: the preview lowers the *numerical tolerance* and the
 * ray count, which are rendering budgets. It does not alter the metric, the connection
 * or the geodesic equation, and the same code path produces both results.
 */
export const NULL_NORMALIZATION_PREVIEW = tolerance({
  id: 'null-normalization-preview',
  quantity: 'g_mu_nu k^mu k^nu in the interactive preview render',
  kind: 'absolute',
  value: 1e-7,
  justification:
    'Absolute, since the target is exactly 0. At the preview integration tolerance of ' +
    '1e-10 the worst residual over a full Schwarzschild image is 1.7e-8, so this bound ' +
    'leaves about six times headroom. Tightening the integrator to 1e-12, where the ' +
    'residual meets the 1e-9 reference gate, costs roughly three times the work and ' +
    'turns a slow render into an unusable one. The validated results in the test suite ' +
    'all use the reference gate; this one is for the picture on screen, and the UI says ' +
    'which it is showing.',
});

/** All declared tolerances, for the UI validation panel (CLAUDE.md §22). */
export const ALL_TOLERANCES: readonly Tolerance[] = Object.freeze([
  NULL_NORMALIZATION_POINTWISE,
  TIMELIKE_NORMALIZATION_POINTWISE,
  NORMALIZATION_DRIFT_LONG_RUN,
  FLAT_PROPAGATION_STRAIGHTNESS,
  CONSERVED_QUANTITY_DRIFT_FLAT,
  CHRISTOFFEL_NUMERIC_VS_ANALYTIC,
  METRIC_INVERSE_RESIDUAL,
  TETRAD_ORTHONORMALITY,
  SCHWARZSCHILD_CHRISTOFFEL_NUMERIC,
  CONSERVED_QUANTITY_DRIFT_SCHWARZSCHILD,
  PHOTON_SPHERE_LOCK,
  DEFLECTION_VS_EXACT,
  CRITICAL_IMPACT_PARAMETER,
  NULL_NORMALIZATION_TRACED,
  NULL_NORMALIZATION_PREVIEW,
]);

/**
 * Check a measured value against a tolerance.
 *
 * `reference` is required for a relative tolerance and ignored for an absolute one.
 * Returns the comparison rather than throwing, so callers can surface a
 * "numerical error exceeds validated tolerance" indicator (CLAUDE.md §17) instead of
 * crashing a render.
 */
export interface ToleranceCheck {
  readonly tolerance: Tolerance;
  readonly measured: number;
  readonly bound: number;
  readonly withinTolerance: boolean;
  readonly message: string;
}

export function checkTolerance(t: Tolerance, measured: number, reference?: number): ToleranceCheck {
  const magnitude = Math.abs(measured);
  let bound: number;

  if (t.kind === 'relative') {
    if (reference === undefined) {
      throw new TypeError(
        `checkTolerance: tolerance '${t.id}' is relative and requires a reference value.`,
      );
    }
    const scale = Math.abs(reference);
    if (scale === 0) {
      throw new RangeError(
        `checkTolerance: tolerance '${t.id}' is relative but the reference value is zero; ` +
          'use an absolute tolerance for a quantity whose target is zero.',
      );
    }
    bound = t.value * scale;
  } else {
    bound = t.value;
  }

  const withinTolerance = Number.isFinite(magnitude) && magnitude <= bound;
  return {
    tolerance: t,
    measured,
    bound,
    withinTolerance,
    message: withinTolerance
      ? `${t.quantity}: ${magnitude.toExponential(3)} within validated tolerance ${bound.toExponential(3)}.`
      : `Numerical error exceeds validated tolerance: ${t.quantity} measured ` +
        `${magnitude.toExponential(3)}, bound ${bound.toExponential(3)} (${t.id}).`,
  };
}
