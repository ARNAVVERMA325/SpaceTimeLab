import { photonSphereRadius } from './schwarzschild.js';

/**
 * Closed-form and quadrature reference results for Schwarzschild null geodesics.
 *
 * These exist to check the geodesic integrator against something derived independently
 * of it. CLAUDE.md §16 requires known analytical behaviour to be tested, and a
 * reference computed by a different method is worth far more than one that reuses the
 * machinery under test.
 *
 * Everything here is equatorial (theta = pi/2), which costs no generality: Schwarzschild
 * is spherically symmetric, so every geodesic lies in a plane through the centre.
 */

/**
 * The impact parameter b = L / E for a null geodesic with turning point r_0.
 *
 * At the turning point dr/dlambda = 0, so E^2 = L^2 f(r_0) / r_0^2 and
 *
 *   b = r_0 / sqrt(1 - 2M / r_0)
 *
 * b is not simply r_0: the two agree only in the weak-field limit.
 */
export function impactParameterForTurningPoint(M: number, r0: number): number {
  const f = 1 - (2 * M) / r0;
  if (!(f > 0)) {
    throw new RangeError(
      `impactParameterForTurningPoint: r_0 = ${r0} is at or inside the horizon r = ${2 * M}.`,
    );
  }
  return r0 / Math.sqrt(f);
}

/**
 * The turning point r_0 for a null geodesic of impact parameter b, if one exists.
 *
 * Inverts the relation above. A turning point exists only for b > b_c = 3 sqrt(3) M,
 * and is then the largest root of r_0^3 = b^2 (r_0 - 2M), which lies in (3M, b].
 * Solved by bisection on that bracket, where the function is monotone.
 *
 * Returns `undefined` when no turning point exists, i.e. the ray is captured.
 */
export function turningPointForImpactParameter(
  M: number,
  b: number,
  tolerance = 1e-14,
): number | undefined {
  const rPhoton = photonSphereRadius(M);
  // g(r) = r^2 / f(r) - b^2 is zero at the turning point. On (3M, b] it increases
  // monotonically, and g(3M) = 27 M^2 - b^2 < 0 exactly when b > 3 sqrt(3) M.
  const g = (r: number): number => (r * r) / (1 - (2 * M) / r) - b * b;

  if (!(g(rPhoton) < 0)) return undefined;

  let lo = rPhoton;
  let hi = Math.max(b, rPhoton * (1 + 1e-12));
  // The turning point never exceeds b, but guard the bracket anyway.
  let guard = 0;
  while (g(hi) < 0 && guard < 200) {
    hi *= 2;
    guard += 1;
  }
  if (!(g(hi) > 0)) return undefined;

  for (let i = 0; i < 200; i += 1) {
    const mid = 0.5 * (lo + hi);
    if (g(mid) < 0) lo = mid;
    else hi = mid;
    if (hi - lo < tolerance * Math.max(1, hi)) break;
  }
  return 0.5 * (lo + hi);
}

/** Nodes and weights for n-point Gauss-Legendre quadrature on [0, 1]. */
function gaussLegendreUnitInterval(n: number): { nodes: Float64Array; weights: Float64Array } {
  const nodes = new Float64Array(n);
  const weights = new Float64Array(n);

  // Newton's method on the Legendre polynomial P_n, using the standard initial guess.
  for (let i = 0; i < n; i += 1) {
    let x = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5));
    let dp = 0;
    for (let iteration = 0; iteration < 100; iteration += 1) {
      let p0 = 1;
      let p1 = 0;
      for (let j = 0; j < n; j += 1) {
        const p2 = p1;
        p1 = p0;
        p0 = ((2 * j + 1) * x * p1 - j * p2) / (j + 1);
      }
      dp = (n * (x * p0 - p1)) / (x * x - 1);
      const dx = -p0 / dp;
      x += dx;
      if (Math.abs(dx) < 1e-15) break;
    }
    // Map the root from [-1, 1] to [0, 1].
    nodes[i] = 0.5 * (x + 1);
    weights[i] = 1 / ((1 - x * x) * dp * dp);
  }
  return { nodes, weights };
}

/**
 * The exact light-deflection angle for a null geodesic with turning point r_0.
 *
 * The orbit equation gives the total azimuthal sweep from infinity in to r_0 and back
 * out, and the deflection is that sweep less the straight-line value pi:
 *
 *   alpha = 2 * Integral_{r_0}^{inf} dr / ( r^2 sqrt( 1/b^2 - f(r)/r^2 ) ) - pi
 *
 * Substituting u = r_0 / r turns this into an integral over [0, 1] whose integrand has
 * an inverse-square-root singularity at u = 1. Factoring the cubic, which has a root
 * there, and then substituting u = 1 - s^2 removes it entirely:
 *
 *   alpha = 4 * Integral_0^1 ds / sqrt( h(1 - s^2) ) - pi
 *   h(u)  = (1 - 2w)(1 + u) - 2w u^2,   w = M / r_0
 *
 * The integrand is now smooth on [0, 1] and Gauss-Legendre converges rapidly — a
 * 20-node rule already agrees with an adaptive evaluation of the original singular form
 * to about 1e-12 for r_0 >= 3.5M.
 *
 * Accuracy degrades as r_0 approaches the photon sphere 3M, where h acquires a double
 * root at u = 1 (h(1) = 2 - 6w vanishes at w = 1/3) and the deflection diverges
 * logarithmically. `deflectionAngleExact` therefore refuses r_0 too close to 3M rather
 * than returning a quietly inaccurate number.
 *
 * There is a second, opposite limitation worth stating. The result is formed as
 * 4 * Integral - pi, a difference of two quantities of order pi. In the weak field the
 * deflection is tiny, so that subtraction cancels almost completely: the *absolute*
 * accuracy stays near 1e-15, but the *relative* accuracy degrades as roughly
 * eps * pi / alpha. At r_0 = 10^7 M, where alpha ~ 4e-7, only about eight significant
 * figures survive. That is a property of the formula, not a defect in the quadrature,
 * and it is why callers should compare against this reference with a combined absolute
 * and relative tolerance rather than a relative one alone (CLAUDE.md §17).
 */
export const DEFLECTION_QUADRATURE_NODES = 128;

/** Closest approach, as a multiple of M, below which the quadrature is not trusted. */
export const DEFLECTION_MIN_TURNING_POINT_FACTOR = 3.05;

export function deflectionAngleExact(
  M: number,
  r0: number,
  nodes = DEFLECTION_QUADRATURE_NODES,
): number {
  const rPhoton = photonSphereRadius(M);
  if (!(r0 > rPhoton)) {
    throw new RangeError(
      `deflectionAngleExact: r_0 = ${r0} is at or inside the photon sphere r = ${rPhoton}. ` +
        'No such escaping null geodesic exists; the deflection is not defined.',
    );
  }
  if (r0 < DEFLECTION_MIN_TURNING_POINT_FACTOR * M) {
    throw new RangeError(
      `deflectionAngleExact: r_0 = ${r0} is within ` +
        `${DEFLECTION_MIN_TURNING_POINT_FACTOR}M of the centre, where the integrand ` +
        'approaches a double root and this quadrature loses accuracy. Refusing rather ' +
        'than returning a silently inaccurate value.',
    );
  }

  const w = M / r0;
  const { nodes: s, weights } = gaussLegendreUnitInterval(nodes);

  let sum = 0;
  for (let i = 0; i < s.length; i += 1) {
    const u = 1 - s[i] * s[i];
    const h = (1 - 2 * w) * (1 + u) - 2 * w * u * u;
    if (!(h > 0)) {
      throw new RangeError(
        `deflectionAngleExact: the integrand became non-positive (h = ${h}) at r_0 = ${r0}. ` +
          'The turning point is too close to the photon sphere for this quadrature.',
      );
    }
    sum += weights[i] / Math.sqrt(h);
  }

  return 4 * sum - Math.PI;
}

/**
 * The azimuthal sweep of a null geodesic between radius `rStart` and infinity.
 *
 *   T(b, r_start) = Integral_{r_start}^{inf} dr / ( r^2 sqrt( 1/b^2 - f(r)/r^2 ) )
 *
 * Substituting u = r_start / r gives
 *
 *   T = Integral_0^1 du / sqrt( r_start^2/b^2 - u^2 (1 - 2 M u / r_start) )
 *
 * which is smooth whenever r_start lies outside the turning point, since the radicand
 * is then bounded away from zero on the whole interval.
 *
 * This exists so that a ray traced over a finite radial range can be compared with the
 * asymptotic deflection *exactly*. A geodesic integrated from r_start in to its turning
 * point and back out to r_start sweeps alpha + pi - 2T, so adding 2T back recovers the
 * asymptotic angle with no flat-space approximation anywhere. The obvious alternative,
 * treating the region beyond r_start as flat and adding arcsin(b / r_start), is only the
 * M -> 0 limit of this integral and leaves an O(M / r_start) error.
 */
export function asymptoticSweepTail(
  M: number,
  b: number,
  rStart: number,
  nodes = DEFLECTION_QUADRATURE_NODES,
): number {
  const { nodes: u, weights } = gaussLegendreUnitInterval(nodes);
  const scale = (rStart * rStart) / (b * b);

  let sum = 0;
  for (let i = 0; i < u.length; i += 1) {
    const ui = u[i];
    const radicand = scale - ui * ui * (1 - (2 * M * ui) / rStart);
    if (!(radicand > 0)) {
      throw new RangeError(
        `asymptoticSweepTail: the radicand became non-positive (${radicand}) at u = ${ui}. ` +
          `r_start = ${rStart} is at or inside the turning point for b = ${b}, so no ray ` +
          'reaches infinity from there along this branch.',
      );
    }
    sum += weights[i] / Math.sqrt(radicand);
  }
  return sum;
}

/**
 * The weak-field deflection alpha = 4M/b (ROADMAP.md 2A.4).
 *
 * The leading term of the expansion in M/b. In geometric units 4GM/(c^2 b) is 4M/b.
 * This is the limit the exact result must approach as b grows, not a substitute for it:
 * at b = 10M it is already about 28% low, and at the photon sphere it is meaningless.
 */
export function deflectionAngleWeakField(M: number, b: number): number {
  return (4 * M) / b;
}
