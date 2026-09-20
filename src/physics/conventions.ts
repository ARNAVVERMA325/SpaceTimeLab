/**
 * Declared mathematical conventions (CLAUDE.md §2).
 *
 * CLAUDE.md §2 requires the metric signature, unit system, coordinate ordering,
 * integration parameter, normalization convention and index convention to be made
 * explicit, and forbids silently mixing them. Every module in the physics engine
 * reads its conventions from here rather than assuming them locally, and the UI
 * provenance panel (CLAUDE.md §22) renders this object verbatim.
 */

/** Metric signature. Fixed by CLAUDE.md §3, which states g_mu_nu u^mu u^nu = -1. */
export const SIGNATURE = '(-,+,+,+)' as const;
export type Signature = typeof SIGNATURE;

/**
 * Minkowski metric eta_ab in an orthonormal frame, for the declared signature.
 * Used as the target of tetrad orthonormality checks (CLAUDE.md §4).
 *
 * Held as a frozen plain array rather than a Float64Array: a typed array with elements
 * cannot be frozen, and this constant is read-only shared state that must not be
 * mutable. It is indexed a*4 + b and only ever read.
 */
export const ETA_AB: readonly number[] = Object.freeze([
  -1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1,
]);

/**
 * Normalization targets for g_mu_nu t^mu t^nu (CLAUDE.md §3, §16).
 *
 * These are exact mathematical targets, not tolerances. Tolerances live in
 * `src/physics/validation/tolerances.ts`, where each one is justified separately
 * as CLAUDE.md §17 requires.
 */
export const NORMALIZATION_TARGET = {
  /** g_mu_nu k^mu k^nu = 0 for a null wavevector. */
  null: 0,
  /** g_mu_nu u^mu u^nu = -1 for a four-velocity in the (-,+,+,+) signature. */
  timelike: -1,
} as const;

export interface ConventionSet {
  readonly signature: Signature;
  readonly unitSystem: string;
  readonly unitSystemNote: string;
  readonly indexConvention: string;
  readonly nullNormalization: string;
  readonly timelikeNormalization: string;
  readonly floatingPoint: string;
}

/**
 * The project-wide convention set.
 *
 * `floatingPoint` is stated deliberately: CLAUDE.md §8 forbids claiming that a
 * CPU/WASM path is automatically "high precision". IEEE-754 binary64 is higher
 * precision than f32 and is still finite precision.
 */
export const CONVENTIONS: ConventionSet = Object.freeze({
  signature: SIGNATURE,
  unitSystem: 'Geometric units, G = c = 1',
  unitSystemNote:
    'Lengths, times and masses share one unit. The Schwarzschild horizon sits at r = 2M in these units.',
  indexConvention:
    'Greek indices mu, nu, alpha, beta run 0..3 over spacetime. Latin frame indices (a),(b) label tetrad legs. Repeated upper/lower pairs are summed.',
  nullNormalization: 'g_mu_nu k^mu k^nu = 0',
  timelikeNormalization: 'g_mu_nu u^mu u^nu = -1',
  floatingPoint:
    'IEEE-754 binary64 (f64) on the CPU reference path. Higher precision than f32, but still finite precision.',
});
