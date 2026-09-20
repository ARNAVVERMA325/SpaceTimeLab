import { ETA_AB } from '../conventions.js';
import type { Index4, Vec4 } from '../core/indices.js';
import type { MetricTensor } from '../core/metric-tensor.js';

/**
 * Observer layer (CLAUDE.md §4, §20): an orthonormal tetrad frame.
 *
 * A tetrad is four vectors e^mu_(a) satisfying
 *
 *   g_mu_nu e^mu_(a) e^nu_(b) = eta_ab
 *
 * with eta_ab = diag(-1,1,1,1). Leg (0) is the observer's four-velocity u^mu; legs
 * (1),(2),(3) span the observer's local orthonormal spatial frame.
 *
 * CLAUDE.md §4 requires quantities an observer actually measures — photon direction,
 * frequency, aberration, observed angular position — to be defined relative to this
 * frame rather than to the coordinate basis.
 *
 * Stored flat, 16 entries, index a*4 + mu.
 */
export class Tetrad {
  readonly components: Float64Array;

  constructor(components: Float64Array) {
    if (components.length !== 16) {
      throw new RangeError('Tetrad expects 16 components (4 legs x 4 coordinate components).');
    }
    this.components = components;
  }

  /** e^mu_(a) */
  e(a: Index4, mu: Index4): number {
    return this.components[a * 4 + mu];
  }

  /** Leg (a) as a contravariant four-vector. */
  leg(a: Index4): Vec4 {
    const base = a * 4;
    return [
      this.components[base],
      this.components[base + 1],
      this.components[base + 2],
      this.components[base + 3],
    ];
  }

  /** The observer's four-velocity u^mu = e^mu_(0). */
  four_velocity_u(): Vec4 {
    return this.leg(0);
  }

  /**
   * The identity frame: e^mu_(a) = delta^mu_a.
   *
   * Orthonormal precisely when the metric components equal eta_ab at the event, which
   * holds throughout Minkowski spacetime in Cartesian coordinates. This is the M1
   * observer frame; the general construction for a curved metric arrives with M3.
   */
  static identity(): Tetrad {
    const components = new Float64Array(16);
    for (let a = 0; a < 4; a += 1) components[a * 4 + a] = 1;
    return new Tetrad(components);
  }
}

/**
 * Largest |g_mu_nu e^mu_(a) e^nu_(b) - eta_ab| over all leg pairs.
 *
 * The direct test that a frame really is orthonormal, and the reason the M1 observer
 * frame can be asserted correct rather than assumed correct.
 */
export function orthonormalityResidual(metric: MetricTensor, tetrad: Tetrad): number {
  let worst = 0;
  for (let a = 0; a < 4; a += 1) {
    for (let b = 0; b < 4; b += 1) {
      let sum = 0;
      for (let mu = 0; mu < 4; mu += 1) {
        const eAmu = tetrad.components[a * 4 + mu];
        if (eAmu === 0) continue;
        for (let nu = 0; nu < 4; nu += 1) {
          const g = metric.g_mu_nu[mu * 4 + nu];
          if (g !== 0) sum += g * eAmu * tetrad.components[b * 4 + nu];
        }
      }
      const residual = Math.abs(sum - ETA_AB[a * 4 + b]);
      if (residual > worst) worst = residual;
    }
  }
  return worst;
}

/**
 * Convert frame components V^(a) into coordinate components V^mu = e^mu_(a) V^(a).
 *
 * This is step 2 of the rendering pipeline in CLAUDE.md §9: turning an observer-frame
 * ray into spacetime phase-space initial conditions.
 */
export function frameToCoordinate(tetrad: Tetrad, frameComponents: Vec4): Vec4 {
  const out: [number, number, number, number] = [0, 0, 0, 0];
  for (let mu = 0; mu < 4; mu += 1) {
    let sum = 0;
    for (let a = 0; a < 4; a += 1) {
      sum += tetrad.components[a * 4 + mu] * frameComponents[a];
    }
    out[mu] = sum;
  }
  return out;
}
