import type { Vec4 } from '../core/indices.js';
import { nullState, type PhaseSpaceState } from '../core/phase-space.js';
import { HAMILTONIAN, type GeodesicFormulation } from '../geodesic/formulation.js';
import { integrateGeodesic, type IntegrationLimits } from '../geodesic/integrate.js';
import type { Integrator } from '../geodesic/integrators/integrator.js';
import { asymptoticSweepTail } from './schwarzschild-analytic.js';
import type { CartesianVec3 } from './spacetime-model.js';
import { photonSphereRadius, type SchwarzschildModel } from './schwarzschild.js';

/**
 * Equatorial null-ray setup and deflection measurement for Schwarzschild
 * (ROADMAP.md 2A.3).
 *
 * Equatorial without loss of generality: Schwarzschild is spherically symmetric, so
 * every geodesic lies in a plane through the centre, and that plane can always be taken
 * to be theta = pi/2.
 */

/**
 * A null geodesic in the equatorial plane with energy E and impact parameter b.
 *
 * The two Killing constants fix the tangent completely. With E = -p_t and L = p_phi,
 *
 *   k^t   = E / f
 *   k^phi = L / r^2
 *   k^r   = +/- sqrt( E^2 - f L^2 / r^2 )
 *
 * and b = L / E. The radial component follows from the null condition rather than being
 * chosen: substituting the three into g_mu_nu k^mu k^nu gives
 * -E^2/f + (k^r)^2/f + L^2/r^2 = 0, which is satisfied exactly by the expression above.
 */
export function equatorialNullRay(
  model: SchwarzschildModel,
  r: number,
  impactParameter: number,
  direction: 'ingoing' | 'outgoing',
  energy_E = 1,
  phi = 0,
): PhaseSpaceState {
  const f = model.lapseFunction(r);
  if (!(f > 0)) {
    throw new RangeError(`equatorialNullRay: r = ${r} is at or inside the horizon.`);
  }

  const angular_momentum_Lz = impactParameter * energy_E;
  const radialSquared =
    energy_E * energy_E - (f * angular_momentum_Lz * angular_momentum_Lz) / (r * r);

  if (radialSquared < 0) {
    throw new RangeError(
      `equatorialNullRay: (k^r)^2 = ${radialSquared} < 0 at r = ${r} for b = ` +
        `${impactParameter}. This radius lies inside the turning point, so no null ` +
        'geodesic with that impact parameter passes through it.',
    );
  }

  const radial = Math.sqrt(radialSquared);
  const null_wavevector_k: Vec4 = [
    energy_E / f,
    direction === 'ingoing' ? -radial : radial,
    0,
    angular_momentum_Lz / (r * r),
  ];

  return nullState([0, r, Math.PI / 2, phi], null_wavevector_k);
}

/**
 * Whether a photon at this event is already certain to be captured.
 *
 * Exact for Schwarzschild, not a tuned cutoff. The null effective potential
 * V(r) = f / r^2 has V'(r) = (6M - 2r) / r^4, which is positive for r < 3M, so V is
 * strictly increasing there. A photon moving inward below r = 3M therefore sees
 * E^2 - L^2 V(r) grow monotonically as r decreases and can never reach a turning point.
 *
 * Terminating on this condition means a shadow renderer never has to integrate down to
 * r = 2M, where the chart degenerates and the connection diverges — the ray is already
 * known to be captured while the geometry is still perfectly well behaved. CLAUDE.md
 * §6.3 endorses exactly this: terminate rays that cross the horizon rather than
 * integrating through the interior, using the simplest validated formulation.
 */
export function isCaptured(model: SchwarzschildModel, position_x: Vec4, tangent: Vec4): boolean {
  return position_x[1] < photonSphereRadius(model.M) && tangent[1] < 0;
}

export interface DeflectionMeasurement {
  /** True when the ray fell in rather than escaping. */
  readonly captured: boolean;
  /** Which formulation of the geodesic equation was integrated. */
  readonly formulation: string;
  /** The asymptotic deflection angle, in radians. Undefined for a captured ray. */
  readonly deflectionAngle?: number;
  /** Total azimuthal sweep actually integrated, between the start and end radii. */
  readonly sweptAngle: number;
  /** Radius at the located turning point; NaN for a captured ray, which has none. */
  readonly minimumRadius: number;
  /** Radius at which the trace actually ended, which overshoots the requested one. */
  readonly endRadius?: number;
  readonly steps: number;
  readonly rejectedSteps: number;
  /** Relative drift of E over the whole trace. */
  readonly energyDrift: number;
  /** Relative drift of L_z over the whole trace. */
  readonly angularMomentumDrift: number;
  /** |g_mu_nu k^mu k^nu| at the end of the trace. */
  readonly nullResidual: number;
}

export interface DeflectionOptions {
  readonly model: SchwarzschildModel;
  readonly integrator: Integrator;
  readonly impactParameter: number;
  /** Radius at which the ray starts inbound and at which it is considered escaped. */
  readonly startRadius: number;
  readonly limits?: Partial<IntegrationLimits>;
  /** Defaults to the Hamiltonian formulation. */
  readonly formulation?: GeodesicFormulation;
}

/**
 * Trace one equatorial null geodesic inward from `startRadius` and measure its
 * deflection.
 *
 * The ray is launched inbound, integrated until it either returns past `startRadius`
 * moving outward or is captured, and the swept azimuth is converted to an asymptotic
 * deflection by adding back the exact sweep beyond each endpoint:
 *
 *   alpha = Delta_phi + T(b, r_start) + T(b, r_end) - pi
 *
 * with T from `asymptoticSweepTail`. The two tails are evaluated at the radii the trace
 * actually had, which matters more than it looks: the final step overshoots the target
 * radius by an arbitrary amount, so assuming a symmetric trace and using 2 T(b, r_start)
 * leaves an error of order b (1/r_start - 1/r_end). That error is invisible in the
 * integrator's own diagnostics -- it never rejects a step, because nothing is wrong with
 * the integration -- and it is amplified by the cancellation below, which is how it came
 * to dominate the result before being tracked down.
 *
 * Note the cancellation. Delta_phi is close to pi - T(b, r_start) - T(b, r_end), so
 * alpha is a small difference of quantities of order pi, amplified by roughly pi/alpha.
 * In the weak field that factor runs to thousands, which is why an apparently tiny
 * inconsistency at the endpoints swamps the answer, and why the endpoint treatment has
 * to be exact rather than merely close.
 *
 * Because T is the exact Schwarzschild integral rather than its flat-space limit, the
 * start radius does not have to be enormous; it only has to lie outside the turning
 * point.
 */
export function measureDeflection(options: DeflectionOptions): DeflectionMeasurement {
  const { model, integrator, impactParameter, startRadius } = options;

  const initial = equatorialNullRay(model, startRadius, impactParameter, 'ingoing');
  const energy_E = 1;
  const angular_momentum_Lz = impactParameter;
  const formulation = options.formulation ?? HAMILTONIAN;

  const result = integrateGeodesic({
    model,
    integrator,
    initial,
    session: formulation.bind(model),
    limits: {
      initialStep: 1e-3,
      parameterMax: 200 * startRadius,
      maxSteps: 2_000_000,
      minStep: 1e-14,
      ...options.limits,
    },
    events: [
      {
        // The turning point: k^r passes from negative to positive. Located exactly, so
        // the closest approach is measured rather than sampled at step boundaries.
        id: 'periapsis',
        value: (_x, tangent) => tangent[1],
        direction: 1,
        terminal: false,
      },
      {
        // Return to the starting radius, moving outward. Terminal, and located exactly,
        // so the trace ends on r = r_start instead of one step beyond it.
        id: 'escape',
        value: (x) => x[1] - startRadius,
        direction: 1,
        terminal: true,
      },
    ],
    terminator: (position_x, tangent) => isCaptured(model, position_x, tangent),
  });

  const periapsis = result.events.find((event) => event.id === 'periapsis');
  const minimumRadius = periapsis ? periapsis.state.position_x[1] : Number.NaN;

  const final = result.final;
  const captured = isCaptured(model, final.position_x, final.tangent);

  // E = -p_t and L_z = p_phi. In the Hamiltonian formulation these are read straight
  // from the integrated covariant momentum, where they are conserved by construction; in
  // the Lagrangian formulation they are rebuilt from the contravariant tangent.
  let finalEnergy: number;
  let finalAngularMomentum: number;
  if (result.formulation === 'hamiltonian') {
    finalEnergy = -result.packed[4];
    finalAngularMomentum = result.packed[7];
  } else {
    const f = model.lapseFunction(final.position_x[1]);
    finalEnergy = f * final.tangent[0];
    finalAngularMomentum =
      final.position_x[1] * final.position_x[1] * Math.sin(final.position_x[2]) ** 2 * final.tangent[3];
  }

  const metric = model.metricAt(final.position_x);
  let nullResidual = 0;
  for (let mu = 0; mu < 4; mu += 1) {
    nullResidual += metric.g_mu_nu[mu * 4 + mu] * final.tangent[mu] * final.tangent[mu];
  }

  const sweptAngle = final.position_x[3] - initial.position_x[3];

  const measurement: DeflectionMeasurement = {
    captured,
    formulation: result.formulation,
    sweptAngle,
    minimumRadius,
    steps: result.steps,
    rejectedSteps: result.rejectedSteps,
    energyDrift: Math.abs(finalEnergy - energy_E) / Math.abs(energy_E),
    angularMomentumDrift:
      Math.abs(finalAngularMomentum - angular_momentum_Lz) / Math.abs(angular_momentum_Lz),
    nullResidual: Math.abs(nullResidual),
  };

  if (captured) return measurement;

  const endRadius = final.position_x[1];
  const tailIn = asymptoticSweepTail(model.M, impactParameter, startRadius);
  const tailOut = asymptoticSweepTail(model.M, impactParameter, endRadius);
  return {
    ...measurement,
    endRadius,
    deflectionAngle: sweptAngle + tailIn + tailOut - Math.PI,
  };
}

/**
 * The direction a ray is travelling *at infinity*, from its state at finite radius R.
 *
 * A ray is still being bent at any finite radius, so reading the background off the
 * local propagation direction at r = R biases the image by the deflection still to come.
 * The bias is systematic and scales as M b / R^2: 2.5e-4 rad for b = 10M at R = 200M,
 * about a third of a pixel at 1000 px across. It is removed exactly here.
 *
 * For an outgoing ray, the direction of motion at infinity is parallel to the position
 * vector at infinity. By spherical symmetry the position keeps sweeping, within the
 * orbital plane, through the exact tail angle T(b, R) = asymptoticSweepTail(M, b, R), so
 *
 *   n_infinity = cos(T) r_hat + sin(T) t_hat
 *
 * with r_hat the radial unit vector at R and t_hat the in-plane unit vector along the
 * tangential motion. The impact parameter comes from the static-frame direction at R:
 * a photon arriving at angle psi from radial has b = R sin(psi) / sqrt(f(R)).
 *
 * In flat space T = asin(b/R) = psi, and this reduces to the local direction, as it must.
 * `tests/visualization/asymptotic-direction.test.ts` checks the correction the useful
 * way: the corrected direction must not depend on which R the trace stopped at.
 */
export function asymptoticDirection(
  model: SchwarzschildModel,
  positionWorld: CartesianVec3,
  directionWorld: CartesianVec3,
): CartesianVec3 {
  const R = Math.hypot(positionWorld[0], positionWorld[1], positionWorld[2]);
  const speed = Math.hypot(directionWorld[0], directionWorld[1], directionWorld[2]);
  if (!(R > 0) || !(speed > 0)) {
    throw new RangeError('asymptoticDirection: degenerate position or direction.');
  }

  const rHat: CartesianVec3 = [positionWorld[0] / R, positionWorld[1] / R, positionWorld[2] / R];
  const radial =
    (directionWorld[0] * rHat[0] + directionWorld[1] * rHat[1] + directionWorld[2] * rHat[2]) / speed;
  if (!(radial > 0)) {
    throw new RangeError(
      'asymptoticDirection: the ray is not moving outward, so it has no asymptotic ' +
        'direction from here.',
    );
  }

  const perpendicular: CartesianVec3 = [
    directionWorld[0] / speed - radial * rHat[0],
    directionWorld[1] / speed - radial * rHat[1],
    directionWorld[2] / speed - radial * rHat[2],
  ];
  const sinPsi = Math.hypot(perpendicular[0], perpendicular[1], perpendicular[2]);
  if (sinPsi < 1e-15) return rHat; // purely radial: no further bending.

  const tHat: CartesianVec3 = [
    perpendicular[0] / sinPsi,
    perpendicular[1] / sinPsi,
    perpendicular[2] / sinPsi,
  ];
  const f = model.lapseFunction(R);
  const b = (R * sinPsi) / Math.sqrt(f);
  const T = asymptoticSweepTail(model.M, b, R);
  const c = Math.cos(T);
  const s = Math.sin(T);
  return [c * rHat[0] + s * tHat[0], c * rHat[1] + s * tHat[1], c * rHat[2] + s * tHat[2]];
}
