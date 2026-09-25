import type { Vec4 } from '../core/indices.js';
import { carterConstant, outerHorizonRadius } from './kerr.js';
import type { SpacetimeModel } from './spacetime-model.js';

/**
 * Deciding a Kerr ray's fate from its constants of motion, exactly.
 *
 * Carter's separation gives the radial motion of a null geodesic as
 *
 *   (Sigma dr/dlambda)^2 = R(r)
 *     = [E(r^2 + a^2) - a L_z]^2 - Delta [ (L_z - a E)^2 + Q ]
 *     = E^2 r^4 - (L_z^2 - a^2 E^2 + Q) r^2 + 2 M K r - a^2 Q,
 *
 * with K = (L_z - a E)^2 + Q. R is a quartic in r, and a turning point is a root of it.
 * So a photon moving inward at radius r reaches the horizon if and only if R has no root
 * between r_+ and r: there is nowhere for it to turn around.
 *
 * That makes capture an exact statement rather than a tuned radius. The critical points
 * of R are the roots of a cubic, solved in closed form below, so the test costs a handful
 * of arithmetic operations and no search. It is the Kerr counterpart of Schwarzschild's
 * "r < 3M with k^r < 0", which is exact for the same reason.
 *
 * Stopping a ray here rather than integrating it toward r_+ also keeps the calculation
 * away from where Boyer-Lindquist coordinates fail: as Delta -> 0 both dt/dlambda and
 * dphi/dlambda diverge while r barely moves, so an error-controlled integrator grinds to
 * a halt chasing a coordinate artifact (CLAUDE.md §6.1, §6.3).
 */

export interface ConstantsOfMotion {
  readonly energy_E: number;
  readonly angular_momentum_Lz: number;
  readonly carter_Q: number;
}

/** E = -p_t, L_z = p_phi and Q for a ray, from its position and tangent. */
export function constantsOfMotion(
  model: SpacetimeModel & { readonly M: number; readonly a: number },
  x: Vec4,
  tangent: Vec4,
): ConstantsOfMotion {
  const g = new Float64Array(16);
  model.metricInto(x, g);
  let p_t = 0;
  let p_phi = 0;
  for (let mu = 0; mu < 4; mu += 1) {
    p_t += g[mu] * tangent[mu];
    p_phi += g[12 + mu] * tangent[mu];
  }
  return {
    energy_E: -p_t,
    angular_momentum_Lz: p_phi,
    carter_Q: carterConstant(model.M, model.a, x, tangent),
  };
}

/** R(r), the radial potential of a null geodesic with these constants. */
export function radialPotential(
  M: number,
  a: number,
  constants: ConstantsOfMotion,
  r: number,
): number {
  const { energy_E: E, angular_momentum_Lz: L, carter_Q: Q } = constants;
  const K = (L - a * E) * (L - a * E) + Q;
  return (
    E * E * r * r * r * r - (L * L - a * a * E * E + Q) * r * r + 2 * M * K * r - a * a * Q
  );
}

/**
 * The real roots of x^3 + p x + q, the depressed cubic.
 *
 * Trigonometric form in the three-real-root case and Cardano otherwise, which is the
 * numerically stable pairing: Cardano's formula loses the three-root case to cancellation
 * between complex cube roots.
 */
function depressedCubicRoots(p: number, q: number): number[] {
  if (p === 0 && q === 0) return [0];
  const discriminant = 4 * p * p * p + 27 * q * q;
  if (discriminant <= 0) {
    const m = 2 * Math.sqrt(-p / 3);
    const argument = Math.max(-1, Math.min(1, (3 * q) / (p * m)));
    const phi = Math.acos(argument) / 3;
    return [0, 1, 2].map((k) => m * Math.cos(phi - (2 * Math.PI * k) / 3));
  }
  const half = -q / 2;
  const root = Math.sqrt(discriminant / 108);
  return [Math.cbrt(half + root) + Math.cbrt(half - root)];
}

/** Stationary points of R(r): the real roots of R'(r) = 4 E^2 r^3 - 2 C r + 2 M K. */
function radialPotentialStationaryPoints(M: number, a: number, constants: ConstantsOfMotion): number[] {
  const { energy_E: E, angular_momentum_Lz: L, carter_Q: Q } = constants;
  const K = (L - a * E) * (L - a * E) + Q;
  const leading = 4 * E * E;
  if (leading === 0) return [];
  // r^3 + p r + q = 0 after dividing through by 4E^2; there is no r^2 term to remove.
  const p = (-2 * (L * L - a * a * E * E + Q)) / leading;
  const q = (2 * M * K) / leading;
  return depressedCubicRoots(p, q);
}

/**
 * Whether this ray must reach the horizon.
 *
 * True when the photon is moving inward and R(r) has no root between r_+ and its present
 * radius, so no turning point stands between it and the hole.
 */
export function isCaptured(
  model: SpacetimeModel & { readonly M: number; readonly a: number },
  x: Vec4,
  tangent: Vec4,
): boolean {
  if (!(tangent[1] < 0)) return false; // moving outward, or momentarily at a turning point
  const { M, a } = model;
  const horizon = outerHorizonRadius(M, a);
  const r = x[1];
  if (!(r > horizon)) return true;

  const constants = constantsOfMotion(model, x, tangent);
  // R >= 0 wherever the photon actually is, so only an interior dip can stop it.
  for (const critical of radialPotentialStationaryPoints(M, a, constants)) {
    if (critical <= horizon || critical >= r) continue;
    if (radialPotential(M, a, constants, critical) <= 0) return false;
  }
  return radialPotential(M, a, constants, horizon) > 0;
}
