import { contractLower } from '../core/metric-tensor.js';
import type { PhaseSpaceState } from '../core/phase-space.js';
import type { KillingVector, SpacetimeModel } from '../spacetimes/spacetime-model.js';

/**
 * Numerical-validation layer (CLAUDE.md §20, §16): Killing conserved quantities.
 *
 * For a Killing vector field xi^mu, the scalar xi^mu p_mu = g_mu_nu xi^mu t^nu is
 * constant along a geodesic. The model supplies the sign convention so that the
 * stationary case yields E = -p_t and the axisymmetric case L_z = p_phi, matching
 * CLAUDE.md §16.
 *
 * These checks are only meaningful when the metric actually possesses the symmetry,
 * which is why `model.killingVectors` is the only source of fields used here.
 */
export function conservedQuantity(
  model: SpacetimeModel,
  killing: KillingVector,
  state: PhaseSpaceState,
): number {
  const metric = model.metricAt(state.position_x);
  const xi = killing.at(state.position_x);
  return killing.sign * contractLower(metric, xi, state.tangent);
}

export interface ConservedQuantityDrift {
  readonly id: string;
  readonly name: string;
  readonly initial: number;
  readonly final: number;
  readonly absoluteDrift: number;
  /** |final - initial| / |initial|, or NaN when the initial value is zero. */
  readonly relativeDrift: number;
}

/**
 * Drift of every conserved quantity the model supports, between two states.
 *
 * `relativeDrift` is NaN when the initial value is zero: a relative measure is
 * undefined there, and reporting NaN is more honest than substituting the absolute
 * value and letting a caller compare it against a relative bound.
 */
export function conservedQuantityDrift(
  model: SpacetimeModel,
  initial: PhaseSpaceState,
  final: PhaseSpaceState,
): readonly ConservedQuantityDrift[] {
  return model.killingVectors.map((killing) => {
    const q0 = conservedQuantity(model, killing, initial);
    const q1 = conservedQuantity(model, killing, final);
    const absoluteDrift = Math.abs(q1 - q0);
    return {
      id: killing.id,
      name: killing.conservedQuantityName,
      initial: q0,
      final: q1,
      absoluteDrift,
      relativeDrift: q0 === 0 ? Number.NaN : absoluteDrift / Math.abs(q0),
    };
  });
}
