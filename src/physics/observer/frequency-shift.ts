import type { CartesianVec3 } from '../spacetimes/spacetime-model.js';

/**
 * Frequency shift between emitter and observer (ROADMAP.md 3.3).
 *
 *   g = nu_obs / nu_emit = (k . u_obs) / (k . u_emit)
 *
 * with the same wavevector k at both ends, parallel-transported along the null geodesic by
 * the integration itself. The ratio does not depend on k's normalization or on whether it
 * is future- or past-directed, since both change the numerator and denominator alike.
 *
 * The renderer launches every ray with k^(a) = (-1, n) in the observer's own frame, so
 * k . u_obs = 1 by construction and g = 1 / (k . u_emit). The emitter's four-velocity is
 * handled in the static orthonormal frame at the emission event, where the photon is
 * described by its static-frame energy and direction and the emitter by a Lorentz factor
 * and a velocity. That keeps the calculation chart-independent: nothing below cares
 * whether the ray was integrated in a rotated orbital plane.
 *
 * The same combination 1 + z = 1 / g carries every effect at once — gravitational
 * redshift, the transverse Doppler shift of the orbital motion, and the ordinary Doppler
 * shift that makes the approaching side of a disk brighter and bluer. None is added
 * separately.
 */

export interface StaticFramePhoton {
  /** k^(0) in the static frame. Negative for the renderer's past-directed rays. */
  readonly energy: number;
  /** k^(i) in the static frame, expressed on the world Cartesian axes. */
  readonly momentum: CartesianVec3;
}

export interface StaticFrameEmitter {
  readonly gamma: number;
  /** Velocity relative to the static observer at the same event, world Cartesian axes. */
  readonly velocity: CartesianVec3;
}

/** k . u_emit in the static orthonormal frame, with eta = diag(-1, 1, 1, 1). */
export function photonDotEmitter(photon: StaticFramePhoton, emitter: StaticFrameEmitter): number {
  const { momentum } = photon;
  const { velocity, gamma } = emitter;
  const dot = momentum[0] * velocity[0] + momentum[1] * velocity[1] + momentum[2] * velocity[2];
  return gamma * (-photon.energy + dot);
}

/** g = nu_obs / nu_emit for a ray launched with k . u_obs = 1. */
export function frequencyRatio(photon: StaticFramePhoton, emitter: StaticFrameEmitter): number {
  return 1 / photonDotEmitter(photon, emitter);
}

/**
 * A prograde Keplerian emitter at the disk position `positionWorld`, as seen by the static
 * observer there.
 *
 * Coordinate angular velocity Omega = sqrt(M / r^3); the static observer measures a speed
 * v = r Omega / sqrt(1 - 2M/r) along +phi-hat, i.e. counter-clockwise about the world +z
 * axis. At the ISCO, v = 1/2.
 */
export function keplerianEmitter(M: number, positionWorld: CartesianVec3): StaticFrameEmitter {
  const [x, y] = positionWorld;
  const rho = Math.hypot(x, y);
  const r = Math.hypot(x, y, positionWorld[2]);
  const f = 1 - (2 * M) / r;
  const speed = (r * Math.sqrt(M / (r * r * r))) / Math.sqrt(f);
  if (!(speed < 1)) {
    throw new RangeError(`keplerianEmitter: orbital speed ${speed} at r = ${r} is not subluminal.`);
  }
  return {
    gamma: 1 / Math.sqrt(1 - speed * speed),
    velocity: [(-y / rho) * speed, (x / rho) * speed, 0],
  };
}

/** A static emitter: at rest relative to the static observer. */
export const STATIC_EMITTER: StaticFrameEmitter = Object.freeze({ gamma: 1, velocity: [0, 0, 0] as CartesianVec3 });
