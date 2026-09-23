import { timelikeState, type PhaseSpaceState } from '../core/phase-space.js';
import type { SchwarzschildModel } from './schwarzschild.js';

/**
 * Timelike geodesics of Schwarzschild: circular orbits, the ISCO, epicyclic motion,
 * periapsis precession and radial infall (CLAUDE.md §16, "Schwarzschild ISCO").
 *
 * All closed forms here were derived symbolically from the radial equation
 *
 *   (dr/dtau)^2 = E^2 - V(r),   V(r) = (1 - 2M/r)(1 + L^2 / r^2)
 *
 * and the precession formula was checked three independent ways — the elliptic-integral
 * closed form, a quadrature in the relativistic anomaly, and a direct quadrature in r —
 * which agreed to 22 digits in mpmath.
 */

export interface CircularOrbit {
  readonly radius: number;
  /** Specific energy E = -u_t. */
  readonly energy_E: number;
  /** Specific angular momentum L = u_phi. */
  readonly angular_momentum_Lz: number;
  /** Angular velocity in Schwarzschild coordinate time, dphi/dt. */
  readonly omega: number;
  /** u^t = dt/dtau for this orbit. */
  readonly ut: number;
}

/**
 * The circular timelike orbit at radius r > 3M.
 *
 *   L^2 = M r^2 / (r - 3M),   E^2 = (r - 2M)^2 / (r (r - 3M)),   dphi/dt = sqrt(M / r^3)
 *
 * The last is Kepler's third law, and it holds *exactly* in Schwarzschild coordinate
 * time — a coincidence of this chart that is worth knowing, and that must not be read as
 * saying Newtonian gravity is exact: the orbit is stable only for r > 6M, and in proper
 * time the frequency differs by u^t = 1 / sqrt(1 - 3M/r).
 */
export function circularOrbit(M: number, r: number): CircularOrbit {
  if (!(r > 3 * M)) {
    throw new RangeError(
      `circularOrbit: no timelike circular orbit exists at r = ${r} <= 3M; the orbit would ` +
        'have to move at or beyond the speed of light.',
    );
  }
  return {
    radius: r,
    energy_E: (r - 2 * M) / Math.sqrt(r * (r - 3 * M)),
    angular_momentum_Lz: Math.sqrt((M * r * r) / (r - 3 * M)),
    omega: Math.sqrt(M / (r * r * r)),
    ut: 1 / Math.sqrt(1 - (3 * M) / r),
  };
}

/** The innermost stable circular orbit, r = 6M, where V''(r) = 2M(r - 6M)/(r^3 (r - 3M)) vanishes. */
export function iscoRadius(M: number): number {
  return 6 * M;
}

/**
 * The radial epicyclic angular frequency in coordinate time, omega_r = Omega sqrt(1 - 6M/r).
 *
 * The frequency of small radial oscillations about a circular orbit. It goes to zero at
 * the ISCO, which is precisely what marginal stability means, and it is smaller than the
 * orbital frequency everywhere, which is why orbits precess.
 */
export function radialEpicyclicFrequency(M: number, r: number): number {
  if (!(r > 6 * M)) {
    throw new RangeError(`radialEpicyclicFrequency: circular orbits at r = ${r} <= 6M are not stable.`);
  }
  return Math.sqrt(M / (r * r * r)) * Math.sqrt(1 - (6 * M) / r);
}

export interface EccentricOrbitConstants {
  readonly energy_E: number;
  readonly angular_momentum_Lz: number;
  readonly periapsis: number;
  readonly apoapsis: number;
}

/**
 * E and L for the bound orbit with semi-latus rectum p and eccentricity e (in units of M),
 * whose turning points are r_p = p / (1 + e) and r_a = p / (1 - e):
 *
 *   E^2 = ((p - 2)^2 - 4 e^2) / (p (p - 3 - e^2)),   L^2 = p^2 M^2 / (p - 3 - e^2)
 *
 * Cutler, Kennefick & Poisson, Phys. Rev. D 50, 3816 (1994). A bound orbit requires
 * p > 6 + 2e (the separatrix).
 */
export function eccentricOrbitConstants(M: number, p: number, e: number): EccentricOrbitConstants {
  if (!(e >= 0 && e < 1)) throw new RangeError(`eccentricOrbitConstants: e = ${e} is not in [0, 1).`);
  if (!(p > 6 + 2 * e)) {
    throw new RangeError(
      `eccentricOrbitConstants: p = ${p} is inside the separatrix p = 6 + 2e = ${6 + 2 * e}; ` +
        'there is no bound orbit, only a plunge.',
    );
  }
  const q = p - 3 - e * e;
  return {
    energy_E: Math.sqrt(((p - 2) * (p - 2) - 4 * e * e) / (p * q)),
    angular_momentum_Lz: (p * M) / Math.sqrt(q),
    periapsis: (p * M) / (1 + e),
    apoapsis: (p * M) / (1 - e),
  };
}

/**
 * The complete elliptic integral of the first kind, K(m) with parameter m = k^2.
 *
 * Computed through the arithmetic-geometric mean, K(m) = pi / (2 AGM(1, sqrt(1 - m))),
 * which converges quadratically and reaches machine precision in a handful of steps.
 */
export function ellipticK(m: number): number {
  if (!(m < 1)) throw new RangeError(`ellipticK: parameter m = ${m} must be below 1.`);
  let a = 1;
  let g = Math.sqrt(1 - m);
  for (let i = 0; i < 40 && Math.abs(a - g) > 1e-16 * a; i += 1) {
    const next = 0.5 * (a + g);
    g = Math.sqrt(a * g);
    a = next;
  }
  return Math.PI / (2 * a);
}

/**
 * The exact periapsis advance per radial period, in radians.
 *
 *   Delta_phi = 4 sqrt(p / (p - 6 + 2e)) K(4e / (p - 6 + 2e)) - 2 pi
 *
 * In the weak field this tends to 6 pi M / p, Einstein's perihelion formula; at p = 1000,
 * e = 0.5 the exact value is still 0.46% above it.
 */
export function periapsisPrecessionExact(M: number, p: number, e: number): number {
  eccentricOrbitConstants(M, p, e); // validates the orbit
  const d = p - 6 + 2 * e;
  return 4 * Math.sqrt(p / d) * ellipticK((4 * e) / d) - 2 * Math.PI;
}

/** Einstein's weak-field periapsis advance, 6 pi M / p. */
export function periapsisPrecessionWeakField(M: number, p: number): number {
  return (6 * Math.PI * M) / (p * M);
}

/**
 * An equatorial timelike state with given E and L at radius r.
 *
 * u^t = E / f,  u^phi = L / r^2,  (u^r)^2 = E^2 - f (1 + L^2 / r^2).
 *
 * At a turning point the radial term is zero analytically but the subtraction leaves a
 * rounding residue of either sign, so `atTurningPoint` sets u^r = 0 exactly rather than
 * taking the square root of noise.
 */
export function equatorialTimelikeState(
  model: SchwarzschildModel,
  r: number,
  energy_E: number,
  angular_momentum_Lz: number,
  radial: 'inward' | 'outward' | 'turning-point',
  phi = 0,
): PhaseSpaceState {
  const f = model.lapseFunction(r);
  if (!(f > 0)) throw new RangeError(`equatorialTimelikeState: r = ${r} is not outside the horizon.`);
  let ur = 0;
  if (radial !== 'turning-point') {
    const ur2 = energy_E * energy_E - f * (1 + (angular_momentum_Lz * angular_momentum_Lz) / (r * r));
    if (ur2 < 0) {
      throw new RangeError(
        `equatorialTimelikeState: (u^r)^2 = ${ur2} < 0 at r = ${r}; this radius is not ` +
          'accessible with the given E and L.',
      );
    }
    ur = radial === 'inward' ? -Math.sqrt(ur2) : Math.sqrt(ur2);
  }
  return timelikeState([0, r, Math.PI / 2, phi], [energy_E / f, ur, 0, angular_momentum_Lz / (r * r)]);
}

/**
 * Proper time for radial free fall from rest at infinity, from r0 down to r:
 *
 *   tau = (2 / 3) (r0^(3/2) - r^(3/2)) / sqrt(2M)
 *
 * finite all the way through r = 2M, where Schwarzschild coordinate time diverges. The
 * two statements are about different clocks (CLAUDE.md §18): the infalling observer's own,
 * and the time coordinate normalized to a static observer at infinity.
 */
export function radialInfallProperTime(M: number, r0: number, r: number): number {
  return ((2 / 3) * (Math.pow(r0, 1.5) - Math.pow(r, 1.5))) / Math.sqrt(2 * M);
}

/**
 * Schwarzschild coordinate time for the same infall, up to an additive constant:
 *
 *   t(r) = -2M [ (2/3) x^3 + 2x + ln |(x - 1) / (x + 1)| ],   x = sqrt(r / 2M)
 *
 * The logarithm diverges as r -> 2M: in these coordinates the infalling body never
 * reaches the horizon, although its proper time to do so is finite. That is a property
 * of the time coordinate, not of the fall.
 */
export function radialInfallCoordinateTime(M: number, r: number): number {
  const x = Math.sqrt(r / (2 * M));
  return -2 * M * ((2 / 3) * x * x * x + 2 * x + Math.log(Math.abs((x - 1) / (x + 1))));
}
