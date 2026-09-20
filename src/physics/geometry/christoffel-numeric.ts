import { ChristoffelSymbols } from '../core/christoffel.js';
import type { Index4, Vec4 } from '../core/indices.js';
import type { SpacetimeModel } from '../spacetimes/spacetime-model.js';

/**
 * Differential-geometry layer (CLAUDE.md §20): Christoffel symbols from numerical
 * metric derivatives, for models that do not supply closed-form symbols.
 *
 * Implements CLAUDE.md §2 directly:
 *
 *   Gamma^mu_{alpha beta}
 *     = (1/2) g^{mu sigma} ( d_alpha g_{sigma beta}
 *                          + d_beta  g_{sigma alpha}
 *                          - d_sigma g_{alpha beta} )
 *
 * with each d_alpha g_{sigma beta} taken as a central difference.
 *
 * This is a cross-check and a fallback, not the preferred path. Where a model has
 * analytical Christoffel symbols it should supply them: they are exact, cheaper, and
 * free of the step-size tradeoff below.
 */

/**
 * Default per-coordinate step for the central difference.
 *
 * Central differencing has truncation error O(h^2) and roundoff error O(eps/h), which
 * balance near h ~ cbrt(eps). In binary64, cbrt(eps) ~ 6.06e-6, giving a total error
 * around 4e-11 for O(1) metric components. The step is scaled by the coordinate
 * magnitude so that it stays meaningful far from the origin, with a floor of 1 so it
 * does not collapse to zero at x^alpha = 0.
 */
export const CENTRAL_DIFFERENCE_SCALE = Math.cbrt(Number.EPSILON);

export function defaultStep(x: Vec4, alpha: Index4): number {
  return CENTRAL_DIFFERENCE_SCALE * Math.max(Math.abs(x[alpha]), 1);
}

/**
 * Compute d_alpha g_{sigma beta} at `x` by central differences.
 *
 * Returns a flat 64-entry buffer indexed alpha*16 + sigma*4 + beta.
 */
export function metricDerivatives(
  model: SpacetimeModel,
  x: Vec4,
  step: (x: Vec4, alpha: Index4) => number = defaultStep,
): Float64Array {
  const derivatives = new Float64Array(64);

  for (let alpha = 0 as Index4; alpha < 4; alpha = (alpha + 1) as Index4) {
    const h = step(x, alpha);
    if (!(h > 0) || !Number.isFinite(h)) {
      throw new RangeError(`metricDerivatives: non-positive or non-finite step for coordinate ${alpha}.`);
    }

    const forward: [number, number, number, number] = [x[0], x[1], x[2], x[3]];
    const backward: [number, number, number, number] = [x[0], x[1], x[2], x[3]];
    forward[alpha] = x[alpha] + h;
    backward[alpha] = x[alpha] - h;

    // Use the actual represented separation, not the nominal 2h: x + h may not be
    // exactly representable, and the difference of the stored endpoints is what the
    // metric was really evaluated at.
    const denominator = forward[alpha] - backward[alpha];

    const gForward = model.metricAt(forward).g_mu_nu;
    const gBackward = model.metricAt(backward).g_mu_nu;

    for (let i = 0; i < 16; i += 1) {
      derivatives[alpha * 16 + i] = (gForward[i] - gBackward[i]) / denominator;
    }
  }

  return derivatives;
}

/**
 * Christoffel symbols built from numerically differentiated metric components.
 *
 * The result is symmetric in the lower indices by construction, since the bracket is
 * symmetric under alpha <-> beta.
 */
export function christoffelFromMetricNumeric(
  model: SpacetimeModel,
  x: Vec4,
  step: (x: Vec4, alpha: Index4) => number = defaultStep,
): ChristoffelSymbols {
  const d = metricDerivatives(model, x, step);
  const gInv = model.metricAt(x).g_inv_mu_nu;
  const out = new Float64Array(64);

  for (let mu = 0; mu < 4; mu += 1) {
    for (let alpha = 0; alpha < 4; alpha += 1) {
      for (let beta = alpha; beta < 4; beta += 1) {
        let sum = 0;
        for (let sigma = 0; sigma < 4; sigma += 1) {
          const inv = gInv[mu * 4 + sigma];
          if (inv === 0) continue;
          const dAlpha_g_sigma_beta = d[alpha * 16 + sigma * 4 + beta];
          const dBeta_g_sigma_alpha = d[beta * 16 + sigma * 4 + alpha];
          const dSigma_g_alpha_beta = d[sigma * 16 + alpha * 4 + beta];
          sum += inv * (dAlpha_g_sigma_beta + dBeta_g_sigma_alpha - dSigma_g_alpha_beta);
        }
        const value = 0.5 * sum;
        out[mu * 16 + alpha * 4 + beta] = value;
        out[mu * 16 + beta * 4 + alpha] = value;
      }
    }
  }

  return new ChristoffelSymbols(out);
}

/** Largest componentwise difference between two sets of Christoffel symbols. */
export function christoffelMaxDifference(a: ChristoffelSymbols, b: ChristoffelSymbols): number {
  let worst = 0;
  for (let i = 0; i < 64; i += 1) {
    const d = Math.abs(a.components[i] - b.components[i]);
    if (d > worst) worst = d;
  }
  return worst;
}
