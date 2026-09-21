import type { Vec4 } from '../core/indices.js';
import { nullState, type PhaseSpaceState } from '../core/phase-space.js';
import { geodesicDerivative } from '../geodesic/geodesic-system.js';
import { integrateGeodesic, type IntegrationLimits } from '../geodesic/integrate.js';
import type { Integrator } from '../geodesic/integrators/integrator.js';
import { asymptoticSweepTail } from './schwarzschild-analytic.js';
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
  /** The asymptotic deflection angle, in radians. Undefined for a captured ray. */
  readonly deflectionAngle?: number;
  /** Total azimuthal sweep actually integrated, between the start and end radii. */
  readonly sweptAngle: number;
  /** Smallest radius reached. */
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

  let minimumRadius = startRadius;
  let turnedAround = false;

  const result = integrateGeodesic({
    model,
    integrator,
    initial,
    derivative: geodesicDerivative(model),
    limits: {
      initialStep: 1e-3,
      parameterMax: 200 * startRadius,
      maxSteps: 2_000_000,
      minStep: 1e-14,
      ...options.limits,
    },
    terminator: (position_x, tangent) => {
      const r = position_x[1];
      if (r < minimumRadius) minimumRadius = r;
      if (tangent[1] > 0) turnedAround = true;
      if (isCaptured(model, position_x, tangent)) return true;
      return turnedAround && r >= startRadius;
    },
  });

  const final = result.final;
  const captured = isCaptured(model, final.position_x, final.tangent);

  // E and L_z read straight off the tangent, using the same relations the setup used.
  const f = model.lapseFunction(final.position_x[1]);
  const finalEnergy = f * final.tangent[0];
  const finalAngularMomentum =
    final.position_x[1] *
    final.position_x[1] *
    Math.sin(final.position_x[2]) ** 2 *
    final.tangent[3];

  const metric = model.metricAt(final.position_x);
  let nullResidual = 0;
  for (let mu = 0; mu < 4; mu += 1) {
    nullResidual += metric.g_mu_nu[mu * 4 + mu] * final.tangent[mu] * final.tangent[mu];
  }

  const sweptAngle = final.position_x[3] - initial.position_x[3];

  const measurement: DeflectionMeasurement = {
    captured,
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
