import { NORMALIZATION_TARGET } from '../conventions.js';
import { contractLower } from '../core/metric-tensor.js';
import type { PhaseSpaceState } from '../core/phase-space.js';
import type { SpacetimeModel } from '../spacetimes/spacetime-model.js';

/**
 * Numerical-validation layer (CLAUDE.md §20, §16): geodesic normalization.
 *
 * The invariant g_mu_nu t^mu t^nu is 0 for a null worldline and -1 for a timelike one
 * in the (-,+,+,+) signature. It is not enforced by the integrator, so it is a genuine
 * independent check on the trajectory rather than a restatement of the update rule.
 */

/** The raw invariant g_mu_nu t^mu t^nu at the state's event. */
export function normalizationInvariant(model: SpacetimeModel, state: PhaseSpaceState): number {
  const metric = model.metricAt(state.position_x);
  return contractLower(metric, state.tangent, state.tangent);
}

/** The target value of the invariant for this worldline kind. */
export function normalizationTarget(state: PhaseSpaceState): number {
  return NORMALIZATION_TARGET[state.kind];
}

/** Signed residual: invariant minus its target. Zero for an exactly normalized state. */
export function normalizationResidual(model: SpacetimeModel, state: PhaseSpaceState): number {
  return normalizationInvariant(model, state) - normalizationTarget(state);
}

/**
 * Rescale a timelike tangent so that g_mu_nu u^mu u^nu = -1 exactly.
 *
 * Used when constructing initial conditions, never to paper over drift during an
 * integration: silently renormalizing a drifting trajectory would hide numerical
 * instability, which CLAUDE.md §24 forbids.
 */
export function normalizeTimelike(model: SpacetimeModel, state: PhaseSpaceState): PhaseSpaceState {
  if (state.kind !== 'timelike') {
    throw new TypeError('normalizeTimelike: only a timelike worldline can be normalized to -1.');
  }
  const invariant = normalizationInvariant(model, state);
  if (!(invariant < 0)) {
    throw new RangeError(
      `normalizeTimelike: g_mu_nu u^mu u^nu = ${invariant}, which is not negative. ` +
        'This tangent is not timelike and cannot be normalized to -1.',
    );
  }
  const scale = 1 / Math.sqrt(-invariant);
  return {
    ...state,
    tangent: [
      state.tangent[0] * scale,
      state.tangent[1] * scale,
      state.tangent[2] * scale,
      state.tangent[3] * scale,
    ],
  };
}
