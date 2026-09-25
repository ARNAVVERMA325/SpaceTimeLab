import { equatorialPhotonOrbitRadius } from './kerr.js';

/**
 * The critical curve of a Kerr black hole: the analytic shadow boundary.
 *
 * Bardeen, in *Black Holes* (Les Houches 1972), ed. DeWitt & DeWitt (1973); Cunningham &
 * Bardeen, ApJ 183, 237 (1973).
 *
 * A photon's motion in Kerr is separable, with the radial potential
 *
 *   R(r) = [E(r^2 + a^2) - a L_z]^2 - Delta [ (L_z - a E)^2 + Q ]
 *
 * A *spherical* photon orbit — constant r, generally not equatorial — needs R = R' = 0.
 * Solving those two equations for the two ratios xi = L_z/E and eta = Q/E^2 gives
 *
 *   xi(r)  = [ M (r^2 - a^2) - r Delta ] / [ a (r - M) ]
 *   eta(r) = r^3 [ 4 M Delta - r (r - M)^2 ] / [ a^2 (r - M)^2 ]
 *
 * which were derived here symbolically rather than copied. These orbits are unstable, so
 * they are exactly the rays that separate capture from escape: the shadow's edge is the
 * image of this one-parameter family, with r running between the prograde and retrograde
 * equatorial photon-orbit radii.
 *
 * For an observer far away at polar angle theta_o, the celestial coordinates are the
 * impact parameters
 *
 *   alpha = -xi / sin(theta_o)
 *   beta^2 = eta + a^2 cos^2(theta_o) - xi^2 cot^2(theta_o)
 *
 * with alpha measured perpendicular to the projected spin axis and beta parallel to it,
 * both in units of M and both defined in the limit of infinite observer distance.
 *
 * This is category A in CLAUDE.md §11: an exact analytical result. It is what the traced
 * shadow is checked against, never a replacement for tracing.
 */

export interface PhotonOrbitConstants {
  /** xi = L_z / E for a spherical photon orbit at this radius. */
  readonly xi: number;
  /** eta = Q / E^2, with Q in this project's Carter convention. */
  readonly eta: number;
}

export interface ShadowPoint {
  /** The spherical photon orbit radius this point is the image of. */
  readonly orbitRadius: number;
  /** Impact parameter perpendicular to the projected spin axis, in units of M. */
  readonly alpha: number;
  /** Impact parameter along the projected spin axis; the upper half of the curve. */
  readonly beta: number;
}

/** Below this spin the (xi, eta) parametrization is 0/0 and Schwarzschild's exact circle applies. */
const MINIMUM_SPIN = 1e-8;

function requireSpin(a: number): void {
  if (Math.abs(a) < MINIMUM_SPIN) {
    throw new RangeError(
      `kerr-shadow: the critical-curve parametrization divides by a and is singular at ` +
        `a = ${a}. The Schwarzschild shadow is the exact circle of radius 3 sqrt(3) M; use ` +
        'criticalImpactParameter from schwarzschild.ts rather than taking a limit here.',
    );
  }
}

/** xi and eta for the spherical photon orbit at radius r. */
export function sphericalPhotonOrbitConstants(M: number, a: number, r: number): PhotonOrbitConstants {
  requireSpin(a);
  const delta = r * r - 2 * M * r + a * a;
  return {
    xi: (M * (r * r - a * a) - r * delta) / (a * (r - M)),
    eta: (r * r * r * (4 * M * delta - r * (r - M) * (r - M))) / (a * a * (r - M) * (r - M)),
  };
}

/**
 * The shadow-boundary point for one spherical photon orbit, seen from polar angle
 * `observerPolarAngle` at infinite distance.
 *
 * Returns `undefined` where beta^2 < 0: those orbits are not visible from this
 * inclination, which is why the curve closes rather than running over the whole family.
 */
export function shadowBoundaryPoint(
  M: number,
  a: number,
  observerPolarAngle: number,
  orbitRadius: number,
): ShadowPoint | undefined {
  const { xi, eta } = sphericalPhotonOrbitConstants(M, a, orbitRadius);
  const s = Math.sin(observerPolarAngle);
  const c = Math.cos(observerPolarAngle);
  const cot2 = (c * c) / (s * s);
  const betaSquared = eta + a * a * c * c - xi * xi * cot2;
  if (betaSquared < 0) {
    // The curve closes where beta^2 reaches zero, and the three terms there cancel to
    // nothing. Seen exactly edge-on, cos(pi/2) is 6e-17 rather than 0 in binary64 and the
    // cancellation lands a few parts in 1e31 below zero, which would drop the very points
    // that close the curve. Anything below this scale is that cancellation, not a sign.
    const scale = Math.abs(eta) + a * a * c * c + xi * xi * cot2;
    if (betaSquared < -1e-12 * Math.max(scale, 1)) return undefined;
    return { orbitRadius, alpha: -xi / s, beta: 0 };
  }
  return { orbitRadius, alpha: -xi / s, beta: Math.sqrt(betaSquared) };
}

/**
 * The upper half of the critical curve, sampled uniformly in the orbit radius.
 *
 * The curve is symmetric under beta -> -beta — reflection through the projected spin
 * axis — so the lower half is the mirror image and is not duplicated here. It is *not*
 * symmetric in alpha: that asymmetry is the shadow's displacement, and it is the visible
 * signature of the spin.
 */
export function shadowCurve(
  M: number,
  a: number,
  observerPolarAngle: number,
  samples = 512,
): readonly ShadowPoint[] {
  const window = visibleOrbitWindow(M, a, observerPolarAngle);
  if (!window) return [];
  const points: ShadowPoint[] = [];
  for (let i = 0; i <= samples; i += 1) {
    const r = window.lo + ((window.hi - window.lo) * i) / samples;
    const point = shadowBoundaryPoint(M, a, observerPolarAngle, r);
    if (point) points.push(point);
  }
  return points;
}

/**
 * The range of spherical photon orbits actually visible from this inclination.
 *
 * Away from an equatorial view only part of the family reaches the observer: beta^2 turns
 * negative outside a window in r, and the curve closes at the two roots. Those roots are
 * found by bisection rather than by sampling, because alpha = -xi / sin(theta_o) amplifies
 * an error in r by 1/sin(theta_o), so for a nearly polar view a sampled endpoint can miss
 * the extent by a hundredth of M while the interior of the curve is fine.
 */
function visibleOrbitWindow(
  M: number,
  a: number,
  observerPolarAngle: number,
  scan = 4096,
): { readonly lo: number; readonly hi: number } | undefined {
  requireSpin(a);
  const inner = equatorialPhotonOrbitRadius(M, Math.abs(a), 'prograde');
  const outer = equatorialPhotonOrbitRadius(M, Math.abs(a), 'retrograde');
  const visible = (r: number): boolean => shadowBoundaryPoint(M, a, observerPolarAngle, r) !== undefined;

  let first = -1;
  let last = -1;
  for (let i = 0; i <= scan; i += 1) {
    const r = inner + ((outer - inner) * i) / scan;
    if (visible(r)) {
      if (first < 0) first = i;
      last = i;
    }
  }
  if (first < 0) return undefined;

  const at = (i: number): number => inner + ((outer - inner) * i) / scan;
  const refine = (outsideIndex: number, insideIndex: number): number => {
    let outside = at(outsideIndex);
    let inside = at(insideIndex);
    for (let i = 0; i < 60; i += 1) {
      const mid = 0.5 * (outside + inside);
      if (visible(mid)) inside = mid;
      else outside = mid;
    }
    return inside;
  };

  return {
    lo: first === 0 ? inner : refine(first - 1, first),
    hi: last === scan ? outer : refine(last + 1, last),
  };
}

export interface ShadowExtent {
  readonly alphaMin: number;
  readonly alphaMax: number;
  readonly betaMax: number;
  /** Displacement of the alpha-midpoint from the line of sight: the spin's signature. */
  readonly alphaCentre: number;
  readonly width: number;
  readonly height: number;
}

/** The bounding box of the critical curve, sampled from `shadowCurve`. */
export function shadowExtent(
  M: number,
  a: number,
  observerPolarAngle: number,
  samples = 8192,
): ShadowExtent {
  const points = shadowCurve(M, a, observerPolarAngle, samples);
  if (points.length === 0) {
    throw new RangeError('shadowExtent: no part of the critical curve is visible from this inclination.');
  }
  let alphaMin = Number.POSITIVE_INFINITY;
  let alphaMax = Number.NEGATIVE_INFINITY;
  let betaMax = 0;
  for (const { alpha, beta } of points) {
    alphaMin = Math.min(alphaMin, alpha);
    alphaMax = Math.max(alphaMax, alpha);
    betaMax = Math.max(betaMax, beta);
  }
  return {
    alphaMin,
    alphaMax,
    betaMax,
    alphaCentre: (alphaMin + alphaMax) / 2,
    width: alphaMax - alphaMin,
    height: 2 * betaMax,
  };
}

/**
 * The conserved quantities of the ray that arrives at celestial position (alpha, beta).
 *
 * The inverse of the relations above:
 *
 *   xi = -alpha sin(theta_o),   eta = beta^2 + cos^2(theta_o) (alpha^2 - a^2)
 *
 * with E = 1, which fixes the affine scale and nothing physical. Turning an image
 * position into initial data this way is how a traced ray is compared with the analytic
 * curve: the same (alpha, beta) goes into both.
 */
export function constantsFromImagePosition(
  a: number,
  observerPolarAngle: number,
  alpha: number,
  beta: number,
): PhotonOrbitConstants {
  const s = Math.sin(observerPolarAngle);
  const c = Math.cos(observerPolarAngle);
  return { xi: -alpha * s, eta: beta * beta + c * c * (alpha * alpha - a * a) };
}
