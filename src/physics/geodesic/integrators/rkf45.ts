import type { DerivativeFn } from '../geodesic-system.js';
import {
  DEFAULT_ERROR_TOLERANCE,
  scaledErrorNorm,
  type ErrorTolerance,
  type Integrator,
  type StepResult,
} from './integrator.js';

/**
 * Adaptive Runge-Kutta-Fehlberg 4(5) (ROADMAP.md 1.3, CLAUDE.md §7.1).
 *
 * Fehlberg's embedded pair evaluates six stages and combines them with two different
 * weight vectors, giving a fourth- and a fifth-order estimate from the same work. Their
 * difference estimates the local truncation error, which drives step-size control.
 *
 * Butcher tableau (Fehlberg 1969):
 *
 *   0     |
 *   1/4   | 1/4
 *   3/8   | 3/32        9/32
 *   12/13 | 1932/2197  -7200/2197   7296/2197
 *   1     | 439/216    -8           3680/513    -845/4104
 *   1/2   | -8/27       2          -3544/2565   1859/4104   -11/40
 *   ------+-------------------------------------------------------------------
 *   b(4)  | 25/216      0           1408/2565    2197/4104   -1/5        0
 *   b(5)  | 16/135      0           6656/12825   28561/56430 -9/50       2/55
 *
 * Propagation order: this implementation advances with the fifth-order weights (local
 * extrapolation). The embedded difference b(5) - b(4) therefore estimates the error of
 * the fourth-order solution, which is the one *not* propagated, so the estimate acts as
 * a conservative proxy for the error actually committed. This choice is stated rather
 * than left implicit because it changes what the reported `errorNorm` means; Fehlberg's
 * original formulation propagates the fourth-order solution instead.
 *
 * Step-size control uses exponent 1/5 = 1/(p+1) with p = 4 the order of the error
 * estimate, with a safety factor and clamped growth/shrink so a single anomalous step
 * cannot make the step size jump wildly.
 *
 * CLAUDE.md §7.1 lists exactly this family as appropriate for ray tracing, where
 * different rays encounter very different curvature scales and need local error control.
 */

// Nodes.
const C2 = 1 / 4;
const C3 = 3 / 8;
const C4 = 12 / 13;
const C5 = 1;
const C6 = 1 / 2;

// Stage coefficients.
const A21 = 1 / 4;
const A31 = 3 / 32;
const A32 = 9 / 32;
const A41 = 1932 / 2197;
const A42 = -7200 / 2197;
const A43 = 7296 / 2197;
const A51 = 439 / 216;
const A52 = -8;
const A53 = 3680 / 513;
const A54 = -845 / 4104;
const A61 = -8 / 27;
const A62 = 2;
const A63 = -3544 / 2565;
const A64 = 1859 / 4104;
const A65 = -11 / 40;

// Fifth-order weights (propagated).
const B5_1 = 16 / 135;
const B5_3 = 6656 / 12825;
const B5_4 = 28561 / 56430;
const B5_5 = -9 / 50;
const B5_6 = 2 / 55;

// Error weights, b(5) - b(4).
const E1 = 1 / 360;
const E3 = -128 / 4275;
const E4 = -2197 / 75240;
const E5 = 1 / 50;
const E6 = 2 / 55;

export interface RKF45Options {
  readonly tolerance?: ErrorTolerance;
  /** Safety factor applied to the predicted optimal step. */
  readonly safety?: number;
  /** Largest allowed step growth in one accepted step. */
  readonly maxGrowth?: number;
  /** Smallest allowed step shrink in one rejected step. */
  readonly minShrink?: number;
}

export class RKF45Integrator implements Integrator {
  readonly id = 'rkf45';
  readonly displayName = 'Runge-Kutta-Fehlberg 4(5) (adaptive)';
  readonly order = 5;
  readonly adaptive = true;
  readonly description =
    'Embedded Fehlberg 4(5) pair. Advances with the fifth-order weights (local ' +
    'extrapolation); the embedded difference estimates the fourth-order error and is ' +
    'used as a conservative error proxy. Step control uses separate absolute and ' +
    'relative tolerances.';

  readonly tolerance: ErrorTolerance;
  private readonly safety: number;
  private readonly maxGrowth: number;
  private readonly minShrink: number;

  private readonly k1: Float64Array;
  private readonly k2: Float64Array;
  private readonly k3: Float64Array;
  private readonly k4: Float64Array;
  private readonly k5: Float64Array;
  private readonly k6: Float64Array;
  private readonly scratch: Float64Array;
  private readonly candidate: Float64Array;
  private readonly error: Float64Array;

  constructor(dimension: number, options: RKF45Options = {}) {
    this.tolerance = options.tolerance ?? DEFAULT_ERROR_TOLERANCE;
    this.safety = options.safety ?? 0.9;
    this.maxGrowth = options.maxGrowth ?? 5;
    this.minShrink = options.minShrink ?? 0.1;

    this.k1 = new Float64Array(dimension);
    this.k2 = new Float64Array(dimension);
    this.k3 = new Float64Array(dimension);
    this.k4 = new Float64Array(dimension);
    this.k5 = new Float64Array(dimension);
    this.k6 = new Float64Array(dimension);
    this.scratch = new Float64Array(dimension);
    this.candidate = new Float64Array(dimension);
    this.error = new Float64Array(dimension);
  }

  step(f: DerivativeFn, p: number, y: Float64Array, h: number): StepResult {
    const { k1, k2, k3, k4, k5, k6, scratch, candidate, error } = this;
    const n = y.length;
    if (n !== k1.length) {
      throw new RangeError(`RKF45Integrator was constructed for dimension ${k1.length}, received ${n}.`);
    }

    f(p, y, k1);

    for (let i = 0; i < n; i += 1) scratch[i] = y[i] + h * A21 * k1[i];
    f(p + C2 * h, scratch, k2);

    for (let i = 0; i < n; i += 1) scratch[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);
    f(p + C3 * h, scratch, k3);

    for (let i = 0; i < n; i += 1) {
      scratch[i] = y[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]);
    }
    f(p + C4 * h, scratch, k4);

    for (let i = 0; i < n; i += 1) {
      scratch[i] = y[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]);
    }
    f(p + C5 * h, scratch, k5);

    for (let i = 0; i < n; i += 1) {
      scratch[i] =
        y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i]);
    }
    f(p + C6 * h, scratch, k6);

    for (let i = 0; i < n; i += 1) {
      candidate[i] =
        y[i] + h * (B5_1 * k1[i] + B5_3 * k3[i] + B5_4 * k4[i] + B5_5 * k5[i] + B5_6 * k6[i]);
      error[i] = h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k6[i]);
    }

    const errorNorm = scaledErrorNorm(error, y, candidate, this.tolerance);

    // A non-finite error norm means the stages themselves blew up. Reject hard and
    // shrink, rather than letting NaN propagate into the state (CLAUDE.md §17).
    if (!Number.isFinite(errorNorm)) {
      return {
        accepted: false,
        usedStep: h,
        nextStep: h * this.minShrink,
        errorNorm,
      };
    }

    const accepted = errorNorm <= 1;

    let factor: number;
    if (errorNorm === 0) {
      factor = this.maxGrowth;
    } else {
      factor = this.safety * Math.pow(errorNorm, -0.2);
      factor = Math.min(this.maxGrowth, Math.max(this.minShrink, factor));
    }

    if (accepted) {
      for (let i = 0; i < n; i += 1) y[i] = candidate[i];
    }

    return {
      accepted,
      usedStep: h,
      nextStep: h * factor,
      errorNorm,
    };
  }
}
