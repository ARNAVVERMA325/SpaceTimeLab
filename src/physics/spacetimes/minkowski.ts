import { CONVENTIONS } from '../conventions.js';
import { ChristoffelSymbols } from '../core/christoffel.js';
import { MetricTensor } from '../core/metric-tensor.js';
import type { Vec4 } from '../core/indices.js';
import {
  IN_DOMAIN,
  type CartesianVec3,
  type ChartGeometry,
  type CoordinateChart,
  type DomainStatus,
  type KillingVector,
  type SpacetimeModel,
  type Symmetries,
} from './spacetime-model.js';

/**
 * Minkowski spacetime in Cartesian coordinates (ROADMAP.md 1.2).
 *
 * Chart: (t, x, y, z), global and geodesically complete.
 * Metric: g_mu_nu = diag(-1, 1, 1, 1) in the (-,+,+,+) signature.
 *
 * The metric components are constant in this chart, so every partial derivative
 * d_alpha g_{sigma beta} vanishes and the Christoffel symbols of CLAUDE.md §2 are
 * identically zero — analytically, not to within some tolerance. ROADMAP.md 1.2 asks
 * for exactly that check, and `christoffelAt` returns exact zeros rather than a
 * numerically differentiated approximation.
 *
 * Minkowski spacetime is flat: the Riemann tensor vanishes identically, so every
 * curvature invariant is zero. It is the flat-space limit that CLAUDE.md §16 requires
 * every curved model to reduce to.
 */
export const MINKOWSKI_CARTESIAN_CHART: CoordinateChart = Object.freeze<CoordinateChart>({
  id: 'minkowski-cartesian',
  displayName: 'Minkowski Cartesian (t, x, y, z)',
  kind: 'cartesian',
  coordinateNames: ['t', 'x', 'y', 'z'],
  horizonPenetrating: false,
  notes:
    'Global inertial chart. No horizon and no coordinate singularity anywhere, so ' +
    'horizon-penetrating coordinates are not applicable to this model.',
});

/**
 * Reference components, held as a frozen plain array because a typed array with
 * elements cannot be frozen and these must not be mutable shared state.
 *
 * In an orthonormal Cartesian chart the inverse metric has the same components as the
 * metric: diag(-1,1,1,1) is its own inverse, so one constant serves for both.
 */
const G_COMPONENTS: readonly number[] = Object.freeze([
  -1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1,
]);

/**
 * Killing vectors of Minkowski spacetime used for conserved-quantity validation.
 *
 * Minkowski has the full ten-parameter Poincare symmetry. The four listed here are the
 * ones the M1 validation suite uses; CLAUDE.md §16 permits conserved quantities only
 * where the metric genuinely has the corresponding symmetry, and all four below are
 * genuine Killing fields of this metric.
 */
const KILLING_VECTORS: readonly KillingVector[] = Object.freeze([
  {
    id: 'd_dt',
    displayName: 'Time translation d/dt (stationarity)',
    conservedQuantityName: 'energy_E',
    // E = -p_t, so the conserved quantity is -(g_mu_nu xi^mu t^nu) for xi = d/dt.
    sign: -1,
    at: (): Vec4 => [1, 0, 0, 0],
  },
  {
    id: 'd_dx',
    displayName: 'Spatial translation d/dx',
    conservedQuantityName: 'momentum_px',
    sign: 1,
    at: (): Vec4 => [0, 1, 0, 0],
  },
  {
    id: 'd_dy',
    displayName: 'Spatial translation d/dy',
    conservedQuantityName: 'momentum_py',
    sign: 1,
    at: (): Vec4 => [0, 0, 1, 0],
  },
  {
    id: 'rotation_z',
    displayName: 'Rotation about the z axis (axisymmetry)',
    // In Cartesian coordinates the axial Killing field d/dphi is xi^mu = (0, -y, x, 0).
    conservedQuantityName: 'angular_momentum_Lz',
    sign: 1,
    at: (x: Vec4): Vec4 => [0, -x[2], x[1], 0],
  },
]);

/**
 * In a Cartesian chart the auxiliary visualization axes are the coordinates themselves,
 * and the metric is already orthonormal, so every conversion is the identity.
 */
const MINKOWSKI_GEOMETRY: ChartGeometry = {
  spatialRadius: (x: Vec4): number => Math.hypot(x[1], x[2], x[3]),
  toCartesianPosition: (x: Vec4): CartesianVec3 => [x[1], x[2], x[3]],
  toCartesianDirection: (_x: Vec4, tangent: Vec4): CartesianVec3 => [
    tangent[1],
    tangent[2],
    tangent[3],
  ],
};

const MINKOWSKI_SYMMETRIES: Symmetries = Object.freeze({
  stationary: true,
  axisymmetric: true,
  // Minkowski is spherically symmetric about any chosen origin. The flag is true, but
  // the Cartesian chart has no polar-axis singularity for the orbital-plane reduction
  // to avoid, so the raytracer has no reason to use it here.
  sphericallySymmetric: true,
});

class MinkowskiSpacetime implements SpacetimeModel {
  readonly id = 'minkowski';
  readonly displayName = 'Minkowski spacetime';
  readonly classification = 'exact-analytical' as const;
  readonly chart = MINKOWSKI_CARTESIAN_CHART;
  readonly conventions = CONVENTIONS;
  readonly parameters: Readonly<Record<string, number>> = Object.freeze({});
  readonly killingVectors = KILLING_VECTORS;
  readonly symmetries = MINKOWSKI_SYMMETRIES;
  readonly geometry = MINKOWSKI_GEOMETRY;
  readonly description =
    'Exact vacuum solution with zero curvature. Evaluating geodesics here is evaluating ' +
    'the consequences of a specified flat geometry, not solving the Einstein field ' +
    'equations dynamically (CLAUDE.md §1.2).';

  metricAt(_x: Vec4): MetricTensor {
    // Fresh copies: MetricTensor exposes mutable typed arrays, and callers must not be
    // able to corrupt the shared module-level constants.
    return new MetricTensor(Float64Array.from(G_COMPONENTS), Float64Array.from(G_COMPONENTS));
  }

  christoffelAt(_x: Vec4): ChristoffelSymbols {
    // Analytically zero: the metric components are constant in this chart.
    return ChristoffelSymbols.zero();
  }

  christoffelInto(_x: Vec4, out: Float64Array): void {
    if (out.length !== 64) {
      throw new RangeError('christoffelInto expects a 64-entry buffer.');
    }
    out.fill(0);
  }

  domainCheck(_x: Vec4): DomainStatus {
    // The Cartesian chart covers all of Minkowski spacetime.
    return IN_DOMAIN;
  }
}

export const minkowski: SpacetimeModel = new MinkowskiSpacetime();
