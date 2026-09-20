import { idx3, type Index4 } from './indices.js';

/**
 * Christoffel symbols of the second kind, Gamma^mu_{alpha beta}, at one event.
 *
 * Defined by CLAUDE.md §2:
 *
 *   Gamma^mu_{alpha beta}
 *     = (1/2) g^{mu sigma} ( d_alpha g_{sigma beta}
 *                          + d_beta  g_{sigma alpha}
 *                          - d_sigma g_{alpha beta} )
 *
 * For a Levi-Civita connection these are symmetric in the two lower indices,
 * Gamma^mu_{alpha beta} = Gamma^mu_{beta alpha}; `symmetryResidual` checks that a
 * given set actually is.
 *
 * Stored flat (64 entries, index mu*16 + alpha*4 + beta) so the geodesic right-hand
 * side can be evaluated without allocating.
 */
export class ChristoffelSymbols {
  readonly components: Float64Array;

  constructor(components: Float64Array) {
    if (components.length !== 64) {
      throw new RangeError('ChristoffelSymbols expects 64 components (4 x 4 x 4).');
    }
    this.components = components;
  }

  get(mu: Index4, alpha: Index4, beta: Index4): number {
    return this.components[idx3(mu, alpha, beta)];
  }

  /** All components identically zero, as in any chart where the connection vanishes. */
  static zero(): ChristoffelSymbols {
    return new ChristoffelSymbols(new Float64Array(64));
  }

  /** Largest |Gamma^mu_{alpha beta}| over all components. */
  maxAbs(): number {
    let m = 0;
    for (let i = 0; i < 64; i += 1) {
      const a = Math.abs(this.components[i]);
      if (a > m) m = a;
    }
    return m;
  }
}

/** Largest |Gamma^mu_{alpha beta} - Gamma^mu_{beta alpha}| over all components. */
export function symmetryResidual(christoffel: ChristoffelSymbols): number {
  let worst = 0;
  for (let mu = 0; mu < 4; mu += 1) {
    for (let alpha = 0; alpha < 4; alpha += 1) {
      for (let beta = alpha + 1; beta < 4; beta += 1) {
        const d = Math.abs(
          christoffel.components[mu * 16 + alpha * 4 + beta] -
            christoffel.components[mu * 16 + beta * 4 + alpha],
        );
        if (d > worst) worst = d;
      }
    }
  }
  return worst;
}
