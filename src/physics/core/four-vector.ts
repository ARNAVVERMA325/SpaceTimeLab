import type { Index4, Vec4 } from './indices.js';

/**
 * Index position of a four-vector's components.
 *
 * Tracked explicitly because CLAUDE.md §24 forbids mixing index positions. A
 * contravariant V^mu and a covariant V_mu are different objects even when both are
 * "four numbers", and the engine should not let them be combined by accident.
 */
export type Variance = 'contravariant' | 'covariant';

export interface FourVector {
  readonly components: Vec4;
  readonly variance: Variance;
}

/** A vector with an upper index, V^mu. */
export function contravariant(components: Vec4): FourVector {
  return { components, variance: 'contravariant' };
}

/** A vector with a lower index, V_mu. */
export function covariant(components: Vec4): FourVector {
  return { components, variance: 'covariant' };
}

export function component(v: FourVector, mu: Index4): number {
  return v.components[mu];
}

/**
 * Assert a four-vector's index position before using it.
 *
 * Throws rather than coercing: a silent index-position mismatch is exactly the class
 * of error CLAUDE.md §24 calls out, and it produces plausible-looking numbers.
 */
export function requireVariance(v: FourVector, expected: Variance, context: string): Vec4 {
  if (v.variance !== expected) {
    throw new TypeError(
      `${context}: expected a ${expected} four-vector (index ${
        expected === 'contravariant' ? 'up' : 'down'
      }), received ${v.variance}.`,
    );
  }
  return v.components;
}
