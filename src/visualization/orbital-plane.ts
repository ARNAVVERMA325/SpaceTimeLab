import type { Vec4 } from '../physics/core/indices.js';
import { nullState, type PhaseSpaceState } from '../physics/core/phase-space.js';
import type { CartesianVec3, SpacetimeModel } from '../physics/spacetimes/spacetime-model.js';

/**
 * Orbital-plane reduction for a spherically symmetric spacetime.
 *
 * In a spherically symmetric metric every geodesic lies in a plane through the centre.
 * That is a theorem, not an approximation: the rotation group acts, the total angular
 * momentum is a conserved vector, and the motion stays orthogonal to it.
 *
 * The reduction exists for a numerical reason. A spherical chart degenerates on the
 * polar axis, where g_phi_phi = r^2 sin^2(theta) vanishes and
 * Gamma^phi_{theta phi} = cot(theta) diverges; a ray whose plane happens to pass near
 * the axis will stall an error-controlled integrator or fail outright, even though
 * nothing physical is happening there. Rotating each ray into its own equatorial plane
 * keeps sin(theta) = 1 for the whole trace, so the axis is never approached at all.
 *
 * The state in the rotated chart obeys the same geodesic equation with the same model:
 * only the initial data is re-expressed. Nothing about the physics is special-cased, and
 * `tests/visualization/orbital-plane.test.ts` checks the reduced trace against a full
 * three-dimensional integration of the same ray.
 *
 * This applies to spherical symmetry alone. Kerr is axisymmetric but not spherically
 * symmetric, its geodesics do not lie in planes, and Milestone 4 must integrate them in
 * full; `reduceToOrbitalPlane` refuses a model that does not declare the symmetry.
 */

export interface OrbitalPlaneFrame {
  /** Unit vector toward the ray's starting position. */
  readonly e1: CartesianVec3;
  /** Unit vector in the plane, along the ray's initial tangential motion. */
  readonly e2: CartesianVec3;
  /** Unit normal of the plane, e1 x e2. */
  readonly normal: CartesianVec3;
}

export interface ReducedRay {
  readonly frame: OrbitalPlaneFrame;
  /** The same ray, re-expressed at theta = pi/2, phi = 0 in the rotated chart. */
  readonly initial: PhaseSpaceState;
  /** Orthonormal tangential speed, retained for diagnostics. */
  readonly tangentialSpeed: number;
}

function cross(a: CartesianVec3, b: CartesianVec3): CartesianVec3 {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

function scale(v: CartesianVec3, s: number): CartesianVec3 {
  return [v[0] * s, v[1] * s, v[2] * s];
}

/** Any unit vector orthogonal to `v`, chosen from the smallest component for stability. */
function anyPerpendicular(v: CartesianVec3): CartesianVec3 {
  const ax = Math.abs(v[0]);
  const ay = Math.abs(v[1]);
  const az = Math.abs(v[2]);
  const axis: CartesianVec3 = ax <= ay && ax <= az ? [1, 0, 0] : ay <= az ? [0, 1, 0] : [0, 0, 1];
  const perpendicular = cross(v, axis);
  const norm = Math.hypot(...perpendicular);
  return scale(perpendicular, 1 / norm);
}

/**
 * Below this tangential speed a ray counts as purely radial and its plane is degenerate.
 *
 * A radial ray has no preferred plane, so any plane containing it will do. The threshold
 * is relative to the total spatial speed, which makes it scale-free.
 */
export const RADIAL_RAY_THRESHOLD = 1e-14;

export function reduceToOrbitalPlane(model: SpacetimeModel, state: PhaseSpaceState): ReducedRay {
  if (!model.symmetries.sphericallySymmetric) {
    throw new RangeError(
      `reduceToOrbitalPlane: ${model.id} does not declare spherical symmetry. Geodesics ` +
        'of an axisymmetric-only spacetime such as Kerr do not lie in planes through the ' +
        'centre, and reducing them would be wrong rather than merely inaccurate.',
    );
  }

  if (model.chart.kind !== 'spherical') {
    throw new RangeError(
      `reduceToOrbitalPlane: the ${model.chart.id} chart is ${model.chart.kind}, not ` +
        'spherical. The reduction re-expresses the ray at (t, r, pi/2, 0), which names a ' +
        'radius and two angles; in a Cartesian chart those slots hold x, y and z, and the ' +
        'result would be silently meaningless. A Cartesian chart has no polar axis to ' +
        'avoid, so it has no use for the reduction either.',
    );
  }

  const position = model.geometry.toCartesianPosition(state.position_x);
  const direction = model.geometry.toCartesianDirection(state.position_x, state.tangent);
  const r = Math.hypot(...position);

  if (!(r > 0) || !Number.isFinite(r)) {
    throw new RangeError(`reduceToOrbitalPlane: degenerate starting position, radius ${r}.`);
  }

  const e1 = scale(position, 1 / r);
  const radialComponent = direction[0] * e1[0] + direction[1] * e1[1] + direction[2] * e1[2];
  const perpendicular: CartesianVec3 = [
    direction[0] - radialComponent * e1[0],
    direction[1] - radialComponent * e1[1],
    direction[2] - radialComponent * e1[2],
  ];
  const tangentialSpeed = Math.hypot(...perpendicular);
  const totalSpeed = Math.hypot(...direction);

  const e2 =
    tangentialSpeed > RADIAL_RAY_THRESHOLD * Math.max(totalSpeed, 1)
      ? scale(perpendicular, 1 / tangentialSpeed)
      : anyPerpendicular(e1);

  const frame: OrbitalPlaneFrame = { e1, e2, normal: cross(e1, e2) };

  // In the rotated chart the ray starts at theta = pi/2, phi = 0. The time and radial
  // components of the tangent are untouched by a spatial rotation about the centre; all
  // of the tangential motion is now azimuthal, and at theta = pi/2 the orthonormal
  // azimuthal component is r k^phi, so k^phi = v_tangential / r.
  const reducedTangent: Vec4 = [state.tangent[0], state.tangent[1], 0, tangentialSpeed / r];

  return {
    frame,
    initial: nullState([state.position_x[0], r, Math.PI / 2, 0], reducedTangent, state.parameter),
    tangentialSpeed,
  };
}

/** Map a vector from the rotated chart's Cartesian axes back to world axes. */
export function liftToWorld(frame: OrbitalPlaneFrame, local: CartesianVec3): CartesianVec3 {
  const { e1, e2, normal } = frame;
  return [
    local[0] * e1[0] + local[1] * e2[0] + local[2] * normal[0],
    local[0] * e1[1] + local[1] * e2[1] + local[2] * normal[1],
    local[0] * e1[2] + local[1] * e2[2] + local[2] * normal[2],
  ];
}

/** The world-frame propagation direction of a state expressed in the rotated chart. */
export function liftDirection(
  model: SpacetimeModel,
  frame: OrbitalPlaneFrame,
  state: PhaseSpaceState,
): CartesianVec3 {
  return liftToWorld(frame, model.geometry.toCartesianDirection(state.position_x, state.tangent));
}

/** The world-frame position of a state expressed in the rotated chart. */
export function liftPosition(
  model: SpacetimeModel,
  frame: OrbitalPlaneFrame,
  state: PhaseSpaceState,
): CartesianVec3 {
  return liftToWorld(frame, model.geometry.toCartesianPosition(state.position_x));
}
