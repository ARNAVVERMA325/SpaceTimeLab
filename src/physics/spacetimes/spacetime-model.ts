import type { ConventionSet } from '../conventions.js';
import type { ChristoffelSymbols } from '../core/christoffel.js';
import type { MetricTensor } from '../core/metric-tensor.js';
import type { Vec4 } from '../core/indices.js';

/**
 * How a spacetime model was obtained (CLAUDE.md §10 and §11).
 *
 * The UI must never silently mix these categories, so every model carries its own
 * classification and the provenance panel renders it.
 */
export type ModelClassification =
  /** A mathematically defined exact solution, e.g. Minkowski, Schwarzschild, Kerr. */
  | 'exact-analytical'
  /** A linearized or weak-field model. Never an exact solution of the full equations. */
  | 'approximate-perturbative'
  /** A specified theoretical metric ansatz, e.g. Morris-Thorne. Not observationally confirmed. */
  | 'theoretical-metric-model'
  /** Imported from an external numerical-relativity simulation, e.g. SXS. */
  | 'imported-numerical-relativity';

/**
 * A coordinate chart (CLAUDE.md §6).
 *
 * Coordinates are representations, not physical objects, so the chart travels with
 * the model and is shown in the UI rather than being assumed by downstream code.
 */
export interface CoordinateChart {
  readonly id: string;
  readonly displayName: string;
  /** Coordinate names in index order 0..3, e.g. ['t','x','y','z'] or ['t','r','theta','phi']. */
  readonly coordinateNames: readonly [string, string, string, string];
  /** Whether integration can continue through an event horizon in this chart (CLAUDE.md §6.2). */
  readonly horizonPenetrating: boolean;
  readonly notes: string;
}

export type DomainExitCode =
  | 'coordinate-breakdown'
  | 'curvature-singularity'
  | 'outside-chart';

/**
 * Whether an event lies inside the chart's valid domain.
 *
 * CLAUDE.md §6.1 is explicit that a horizon must be detected from the metric
 * components or equations actually in use, never from a vanishing metric determinant.
 * Each model implements this check for its own chart.
 */
export type DomainStatus =
  | { readonly inDomain: true }
  | { readonly inDomain: false; readonly code: DomainExitCode; readonly reason: string };

export const IN_DOMAIN: DomainStatus = Object.freeze({ inDomain: true });

/**
 * A Killing vector field and the conserved quantity it generates (CLAUDE.md §16).
 *
 * The conserved quantity along a geodesic with tangent t^mu is
 *
 *   sign * g_mu_nu xi^mu t^nu
 *
 * The explicit `sign` carries the convention difference between E = -p_t (stationarity,
 * sign -1) and L_z = p_phi (axisymmetry, sign +1), so neither is hard-coded downstream.
 *
 * CLAUDE.md §16 also warns that conserved quantities may only be used when the metric
 * actually possesses the corresponding symmetry, which is why this list belongs to the
 * model rather than to the validation layer.
 */
export interface KillingVector {
  readonly id: string;
  readonly displayName: string;
  /** Name of the conserved quantity, using CLAUDE.md §19 naming, e.g. 'energy_E'. */
  readonly conservedQuantityName: string;
  readonly sign: 1 | -1;
  /** The field xi^mu evaluated at an event. */
  at(x: Vec4): Vec4;
}

/**
 * A spacetime model: the geometry layer of CLAUDE.md §20.
 *
 * A model defines the geometry and nothing else. It does not integrate, render, or
 * decide termination policy; those belong to the geodesic, visualization and driver
 * layers respectively.
 */
export interface SpacetimeModel {
  readonly id: string;
  readonly displayName: string;
  readonly classification: ModelClassification;
  readonly chart: CoordinateChart;
  readonly conventions: ConventionSet;
  /** Model parameters in the declared unit system, e.g. { M: 1 } or { M: 1, a: 0.9 }. */
  readonly parameters: Readonly<Record<string, number>>;
  /** Killing vectors this metric actually possesses. */
  readonly killingVectors: readonly KillingVector[];
  /** Short statement of what this model is and is not, for the UI (CLAUDE.md §22). */
  readonly description: string;

  /** g_mu_nu and g_inv_mu_nu at an event. Allocates; not for the hot path. */
  metricAt(x: Vec4): MetricTensor;

  /** Gamma^mu_{alpha beta} at an event. Allocates; not for the hot path. */
  christoffelAt(x: Vec4): ChristoffelSymbols;

  /**
   * Write Gamma^mu_{alpha beta} into a caller-owned 64-entry buffer.
   *
   * This is the allocation-free path the geodesic right-hand side calls once per
   * stage, per step. It must produce the same values as `christoffelAt`.
   */
  christoffelInto(x: Vec4, out: Float64Array): void;

  /** Whether the event lies in this chart's valid domain (CLAUDE.md §6.1). */
  domainCheck(x: Vec4): DomainStatus;
}
