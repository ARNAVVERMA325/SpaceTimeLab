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
   * The orthonormal frame of a static observer in a diagonal metric.
   *
   * For g_mu_nu = diag(g_00, g_11, g_22, g_33) with g_00 < 0 and the spatial entries
   * positive, the coordinate basis is already orthogonal, so orthonormalizing it is
   * just a rescaling:
   *
   *   e^mu_(0) = (1 / sqrt(-g_00), 0, 0, 0)
   *   e^mu_(i) = delta^mu_i / sqrt(g_ii)
   *
   * Leg (0) is the observer's four-velocity, which for Schwarzschild is the familiar
   * u^mu = (1 / sqrt(-g_t_t), 0, 0, 0) = (1 / sqrt(f), 0, 0, 0) of a shell observer
   * hovering at fixed r.
   *
   * This is ROADMAP.md 3.1's static-observer construction, landing early because
   * Milestone 2A's camera has to be a physical observer rather than a labelled fake.
   * The freely-falling observers, and everything that depends on relative motion
   * between emitter and observer, remain Milestone 3 work.
   *
   * Throws for a metric with off-diagonal terms, rather than silently discarding them:
   * Kerr in Boyer-Lindquist has a g_t_phi cross term and needs a genuine
   * orthonormalization, not this shortcut.
   */
  static diagonalStatic(metric: MetricTensor): Tetrad {
    for (let mu = 0; mu < 4; mu += 1) {
      for (let nu = 0; nu < 4; nu += 1) {
        if (mu !== nu && metric.g_mu_nu[mu * 4 + nu] !== 0) {
          throw new RangeError(
            `Tetrad.diagonalStatic: the metric has a non-zero off-diagonal component ` +
              `g_${mu}_${nu} = ${metric.g_mu_nu[mu * 4 + nu]}. This construction assumes a ` +
              'diagonal metric; a chart with cross terms needs a full orthonormalization.',
          );
        }
      }
    }

    const g00 = metric.g_mu_nu[0];
    if (!(g00 < 0)) {
      throw new RangeError(
        `Tetrad.diagonalStatic: g_00 = ${g00} is not negative, so no static observer ` +
          'exists here. Inside a horizon the timelike Killing vector is spacelike and a ' +
          'hovering observer is impossible.',
      );
    }

    const components = new Float64Array(16);
    components[0] = 1 / Math.sqrt(-g00);
    for (let i = 1; i < 4; i += 1) {
      const gii = metric.g_mu_nu[i * 4 + i];
      if (!(gii > 0)) {
        throw new RangeError(
          `Tetrad.diagonalStatic: g_${i}_${i} = ${gii} is not positive; the spatial ` +
            'metric has degenerated in this chart at this event.',
        );
      }
      components[i * 4 + i] = 1 / Math.sqrt(gii);
    }
    return new Tetrad(components);
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
