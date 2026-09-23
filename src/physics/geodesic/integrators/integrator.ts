import type { DerivativeFn } from '../geodesic-system.js';

/**
 * Integrator abstraction (CLAUDE.md §7.3).
 *
 * CLAUDE.md §7 is explicit that there is no universal "correct integrator" and that
 * the architecture must allow several validated ones, chosen per simulation mode. The
 * engine therefore programs against this interface rather than against RK4 or RKF45
 * directly, so that a symplectic method can be added for long-lived bound orbits
 * (§7.2) without touching call sites.
 */

export interface StepResult {
  /**
   * Whether the step was accepted. A fixed-step method always accepts. An adaptive
   * method rejects when its local error estimate exceeds the requested tolerance, and
   * leaves `y` untouched so the caller can retry with `nextStep`.
   */
  readonly accepted: boolean;
  /** The step actually attempted. */
  readonly usedStep: number;
  /** The step size to attempt next. */
  readonly nextStep: number;
  /**
   * Scaled local error norm, 1 meaning "exactly at tolerance". NaN for methods that
   * do not estimate local error.
   */
  readonly errorNorm: number;
}

export interface Integrator {
  readonly id: string;
  readonly displayName: string;
  /** Order of the solution actually propagated. */
  readonly order: number;
  readonly adaptive: boolean;
  /** Short description of the method and its error control, for the UI (CLAUDE.md §22). */
  readonly description: string;

  /**
   * Attempt to advance `y` from parameter `p` by `h`.
   *
   * On an accepted step `y` holds the advanced state. On a rejected step `y` is
   * unchanged. The caller owns the parameter bookkeeping.
   */
  step(f: DerivativeFn, p: number, y: Float64Array, h: number): StepResult;

  /**
   * Continuous extension over the most recent *accepted* step: write the state at
   * p_old + theta * h into `out`, for theta in [0, 1].
   *
   * This is what makes event location exact. Without it a termination can only be
   * detected after the fact, and the state reported is wherever the last step happened
   * to land — the overshoot that corrupted the first deflection measurement in
   * Milestone 2A. Valid only until the next call to `step`.
   */
  interpolate?(theta: number, out: Float64Array): void;

  /** Local order of the continuous extension, when one is provided. */
  readonly denseOutputOrder?: number;
}

/**
 * Local error tolerances for adaptive methods (CLAUDE.md §17).
 *
 * Separate absolute and relative components, because CLAUDE.md §17 forbids one
 * universal numerical-error threshold. The absolute term keeps control meaningful for
 * components passing through zero; the relative term keeps it meaningful for large
 * components such as a coordinate far from the origin.
 */
export interface ErrorTolerance {
  readonly absolute: number;
  readonly relative: number;
}

export const DEFAULT_ERROR_TOLERANCE: ErrorTolerance = Object.freeze({
  absolute: 1e-10,
  relative: 1e-10,
});

/**
 * RMS of the componentwise error scaled by (atol + rtol * |y|).
 *
 * A value of 1 sits exactly at tolerance; the caller accepts at <= 1. The RMS form
 * (rather than a max) keeps a single marginal component from dominating step-size
 * control while still responding to a genuinely bad step.
 */
export function scaledErrorNorm(
  error: Float64Array,
  yOld: Float64Array,
  yNew: Float64Array,
  tolerance: ErrorTolerance,
): number {
  let sum = 0;
  const n = error.length;
  for (let i = 0; i < n; i += 1) {
    const scale =
      tolerance.absolute + tolerance.relative * Math.max(Math.abs(yOld[i]), Math.abs(yNew[i]));

    // A zero scale arises legitimately: a relative-only tolerance against a component
    // passing through zero. Dividing there would give 0/0 = NaN, which the adaptive
    // controller cannot distinguish from a genuine blow-up. Treat it exactly instead —
    // with no allowance, a non-zero error is unboundedly out of tolerance and a zero
    // error is exactly within it.
    if (scale === 0) {
      if (error[i] === 0) continue;
      return Number.POSITIVE_INFINITY;
    }

    const ratio = error[i] / scale;
    sum += ratio * ratio;
  }
  return Math.sqrt(sum / n);
}
