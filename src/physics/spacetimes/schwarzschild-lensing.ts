import { HAMILTONIAN } from '../geodesic/formulation.js';
import { brentRoot, integrateGeodesic, type IntegrationLimits } from '../geodesic/integrate.js';
import type { Integrator } from '../geodesic/integrators/integrator.js';
import { criticalImpactParameter, photonSphereRadius, type SchwarzschildModel } from './schwarzschild.js';
import { equatorialNullRay, isCaptured } from './schwarzschild-rays.js';

/**
 * Gravitational lensing by a Schwarzschild black hole: the strong-deflection limit and
 * Einstein rings for a source on the optical axis (ROADMAP.md 2A.3).
 *
 * Everything here is computed by tracing null geodesics. The closed forms are included
 * as limits to compare against, never as substitutes: the weak-field Einstein angle and
 * Bozza's strong-deflection expansion are each accurate only at one end of the range.
 */

/**
 * Strong-deflection coefficients for Schwarzschild.
 *
 * Bozza, "Gravitational lensing in the strong field limit", Phys. Rev. D 66, 103001
 * (2002): as the impact parameter approaches b_c = 3 sqrt(3) M,
 *
 *   alpha(b) = -a_bar ln(b / b_c - 1) + b_bar + O((b - b_c) ln(b - b_c))
 *
 * with a_bar = 1 and b_bar = -pi + ln(216 (7 - 4 sqrt(3))) ~ -0.40023. Verified
 * independently at 60 digits with mpmath: alpha + ln(b/b_c - 1) - b_bar falls from 2.5e-2
 * at b/b_c - 1 = 1e-2 to 8.9e-12 at 1e-12.
 *
 * a_bar = 1 is the lensing face of the photon-sphere Lyapunov exponent: a ray that loops
 * once more around the black hole needs an impact parameter e^(2 pi / a_bar) = e^(2 pi)
 * times closer to b_c, the same e^pi per half orbit at which perturbations of the
 * circular photon orbit grow.
 */
export const STRONG_DEFLECTION = Object.freeze({
  aBar: 1,
  // Written as ln 216 - ln(7 + 4 sqrt 3) rather than ln(216 (7 - 4 sqrt 3)). The two are
  // equal, since (7 - 4 sqrt 3)(7 + 4 sqrt 3) = 1, but 7 - 4 sqrt 3 = 7 - 6.928... cancels
  // about two digits in binary64; the sum has no cancellation at all.
  bBar: -Math.PI + Math.log(216) - Math.log(7 + 4 * Math.sqrt(3)),
  reference:
    'V. Bozza, "Gravitational lensing in the strong field limit", Phys. Rev. D 66, 103001 (2002).',
});

/** Bozza's strong-deflection approximation to the deflection angle. Valid only as b -> b_c. */
export function strongDeflectionAngle(M: number, b: number): number {
  const bc = criticalImpactParameter(M);
  if (!(b > bc)) {
    throw new RangeError(`strongDeflectionAngle: b = ${b} is not above b_c = ${bc}; the ray is captured.`);
  }
  return -STRONG_DEFLECTION.aBar * Math.log(b / bc - 1) + STRONG_DEFLECTION.bBar;
}

/**
 * The weak-field Einstein angle for a point lens, theta_E = sqrt(4 M D_LS / (D_OL D_OS)).
 *
 * Distances are taken as the coordinate radii of observer and source, which is the
 * flat-space reading the formula assumes; it is the leading term of an expansion in
 * M / b, and at r_O = r_S = 100M it is already 11% low.
 */
export function weakFieldEinsteinAngle(M: number, rObserver: number, rSource: number): number {
  return Math.sqrt((4 * M * rSource) / (rObserver * (rObserver + rSource)));
}

/** The angle from the lens at which a static observer at r_O sees a ray of impact b. */
export function observedAngle(model: SchwarzschildModel, rObserver: number, b: number): number {
  return Math.asin((b * Math.sqrt(model.lapseFunction(rObserver))) / rObserver);
}

export interface SweepOptions {
  readonly model: SchwarzschildModel;
  readonly integrator: Integrator;
  readonly impactParameter: number;
  readonly rObserver: number;
  readonly rSource: number;
  readonly limits?: Partial<IntegrationLimits>;
}

/**
 * The azimuth swept by a null geodesic from the observer to the source radius.
 *
 * The ray leaves the observer at r_O heading inward, passes its turning point, and is
 * stopped exactly on r = r_S on the way out by an event. Returns undefined if it is
 * captured instead. By time-reversal symmetry this is also the sweep of the physical
 * photon travelling from source to observer.
 */
export function traceSweep(options: SweepOptions): number | undefined {
  const { model, integrator, impactParameter, rObserver, rSource } = options;
  const initial = equatorialNullRay(model, rObserver, impactParameter, 'ingoing');
  let captured = false;
  const result = integrateGeodesic({
    model,
    integrator,
    initial,
    session: HAMILTONIAN.bind(model),
    limits: {
      initialStep: 1e-3,
      parameterMax: 100 * (rObserver + rSource),
      maxSteps: 4_000_000,
      minStep: 1e-14,
      ...options.limits,
    },
    events: [{ id: 'source', value: (x) => x[1] - rSource, direction: 1, terminal: true }],
    terminator: (x, k) => {
      captured = isCaptured(model, x, k);
      return captured;
    },
  });
  if (captured || result.reason !== 'event') return undefined;
  return result.final.position_x[3] - initial.position_x[3];
}

export interface EinsteinRing {
  /** Ring order: the number of complete loops around the black hole. */
  readonly order: number;
  readonly impactParameter: number;
  /** Angular radius of the ring as seen by a static observer at r_O. */
  readonly observedAngle: number;
  /** The traced sweep at the solution, which should equal (2n + 1) pi. */
  readonly sweep: number;
  /** Rays traced to find it. */
  readonly evaluations: number;
}

/**
 * Solve the exact finite-distance lens equation for a source on the optical axis.
 *
 * Observer at r_O on one side of the black hole, point source at r_S directly behind it.
 * By symmetry every ray reaching the source lies in a plane containing the axis, and the
 * image of the source is a ring. The ray of order n loops n times around the hole, so its
 * sweep from observer to source is (2n + 1) pi. That condition is solved for the impact
 * parameter by Brent's method, with every evaluation a fully traced geodesic.
 *
 * Near b_c the sweep diverges logarithmically, so the relativistic rings n >= 1 are
 * searched in the variable x = ln(b / b_c - 1), where it is close to linear.
 */
export function solveEinsteinRing(
  model: SchwarzschildModel,
  integrator: Integrator,
  order: number,
  rObserver: number,
  rSource: number,
): EinsteinRing {
  if (!Number.isInteger(order) || order < 0) {
    throw new RangeError(`solveEinsteinRing: ring order must be a non-negative integer, received ${order}.`);
  }
  const M = model.M;
  const bc = criticalImpactParameter(M);
  const target = (2 * order + 1) * Math.PI;
  let evaluations = 0;

  const sweepFor = (b: number): number => {
    evaluations += 1;
    const sweep = traceSweep({ model, integrator, impactParameter: b, rObserver, rSource });
    // A captured ray sweeps "infinitely" far: treat it as overshooting the target.
    return sweep === undefined ? Number.POSITIVE_INFINITY : sweep;
  };

  let b: number;
  if (order === 0) {
    // The sweep falls monotonically from +infinity at b_c. The largest usable b is set by
    // the ray still having to pass through both r_O and r_S.
    const rMin = Math.min(rObserver, rSource);
    const hi = 0.999 * rMin * Math.sqrt(1 / (1 - (2 * M) / rMin));
    const lo = bc * (1 + 1e-3);
    const g = (value: number): number => sweepFor(value) - target;
    const gLo = g(lo);
    const gHi = g(hi);
    if (!(gLo > 0 && gHi < 0)) {
      throw new RangeError(
        `solveEinsteinRing: no primary ring bracketed for r_O = ${rObserver}, r_S = ${rSource} ` +
          `(sweep - target = ${gLo} at b = ${lo}, ${gHi} at b = ${hi}).`,
      );
    }
    b = brentRoot(g, lo, hi, gLo, gHi, 1e-15);
  } else {
    const guess = STRONG_DEFLECTION.bBar - 2 * Math.PI * order;
    const g = (x: number): number => sweepFor(bc * (1 + Math.exp(x))) - target;
    const lo = guess - 2;
    const hi = guess + 2;
    const gLo = g(lo);
    const gHi = g(hi);
    if (!(gLo > 0 && gHi < 0)) {
      throw new RangeError(
        `solveEinsteinRing: ring ${order} not bracketed around ln(b/b_c - 1) = ${guess}.`,
      );
    }
    const x = brentRoot(g, lo, hi, gLo, gHi, 1e-15);
    b = bc * (1 + Math.exp(x));
  }

  if (!(b > bc) || b > photonSphereRadius(M) * 1e9) {
    throw new RangeError(`solveEinsteinRing: implausible solution b = ${b}.`);
  }

  return {
    order,
    impactParameter: b,
    observedAngle: observedAngle(model, rObserver, b),
    sweep: sweepFor(b),
    evaluations,
  };
}
