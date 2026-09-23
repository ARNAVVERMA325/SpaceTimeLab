import type { DerivativeFn } from '../geodesic-system.js';
import {
  DEFAULT_ERROR_TOLERANCE,
  scaledErrorNorm,
  type ErrorTolerance,
  type Integrator,
  type StepResult,
} from './integrator.js';

/**
 * Dormand-Prince 5(4) with FSAL and a fourth-order continuous extension.
 *
 * Dormand & Prince, "A family of embedded Runge-Kutta formulae", J. Comput. Appl. Math.
 * 6 (1980) 19-26. The continuous extension is Shampine's, "Some practical Runge-Kutta
 * formulas", Math. Comp. 46 (1986) 135-150.
 *
 * Why this replaces Fehlberg 4(5) as the default. Both are six-stage embedded pairs, but
 * Dormand-Prince was constructed to minimize the error coefficients of the *propagated*
 * fifth-order solution, where Fehlberg optimized the fourth-order one and so is poorly
 * suited to local extrapolation. It is also first-same-as-last: the derivative at the
 * new point is the seventh stage of the current step, so an accepted step costs six
 * evaluations rather than seven. It is the standard non-stiff integrator for this reason
 * (MATLAB ode45, SciPy RK45, Hairer-Norsett-Wanner DOPRI5).
 *
 * The tableau and the dense-output matrix are transcribed from SciPy's implementation
 * (scipy.integrate._ivp.rk.RK45) rather than from memory, and
 * `tests/physics/dopri5.test.ts` checks the order conditions, the consistency relations
 * and the convergence order independently.
 *
 *   0    |
 *   1/5  | 1/5
 *   3/10 | 3/40        9/40
 *   4/5  | 44/45      -56/15       32/9
 *   8/9  | 19372/6561 -25360/2187  64448/6561  -212/729
 *   1    | 9017/3168  -355/33      46732/5247   49/176     -5103/18656
 *   1    | 35/384      0           500/1113     125/192    -2187/6784    11/84
 *   -----+----------------------------------------------------------------------------
 *   b    | 35/384      0           500/1113     125/192    -2187/6784    11/84      0
 *   e    | -71/57600   0           71/16695    -71/1920    17253/339200 -22/525    1/40
 *
 * where e is SciPy's error vector: the difference between the fifth- and fourth-order
 * weights, including the FSAL stage.
 */

const C2 = 1 / 5;
const C3 = 3 / 10;
const C4 = 4 / 5;
const C5 = 8 / 9;

const A21 = 1 / 5;
const A31 = 3 / 40;
const A32 = 9 / 40;
const A41 = 44 / 45;
const A42 = -56 / 15;
const A43 = 32 / 9;
const A51 = 19372 / 6561;
const A52 = -25360 / 2187;
const A53 = 64448 / 6561;
const A54 = -212 / 729;
const A61 = 9017 / 3168;
const A62 = -355 / 33;
const A63 = 46732 / 5247;
const A64 = 49 / 176;
const A65 = -5103 / 18656;

const B1 = 35 / 384;
const B3 = 500 / 1113;
const B4 = 125 / 192;
const B5 = -2187 / 6784;
const B6 = 11 / 84;

const E1 = -71 / 57600;
const E3 = 71 / 16695;
const E4 = -71 / 1920;
const E5 = 17253 / 339200;
const E6 = -22 / 525;
const E7 = 1 / 40;

/**
 * Shampine's dense-output matrix, rows = stages 1..7, columns = powers theta^1..theta^4.
 *
 * y(p_old + theta h) = y_old + h * sum_i K_i * sum_j P[i][j] theta^(j+1).
 * At theta = 1 each row sums to the corresponding weight b_i, so the extension passes
 * exactly through the propagated solution.
 */
const P: readonly (readonly [number, number, number, number])[] = [
  [1, -8048581381 / 2820520608, 8663915743 / 2820520608, -12715105075 / 11282082432],
  [0, 0, 0, 0],
  [0, 131558114200 / 32700410799, -68118460800 / 10900136933, 87487479700 / 32700410799],
  [0, -1754552775 / 470086768, 14199869525 / 1410260304, -10690763975 / 1880347072],
  [0, 127303824393 / 49829197408, -318862633887 / 49829197408, 701980252875 / 199316789632],
  [0, -282668133 / 205662961, 2019193451 / 616988883, -1453857185 / 822651844],
  [0, 40617522 / 29380423, -110615467 / 29380423, 69997945 / 29380423],
];

export interface Dopri5Options {
  readonly tolerance?: ErrorTolerance;
  readonly safety?: number;
  readonly maxGrowth?: number;
  readonly minShrink?: number;
}

export class Dopri5Integrator implements Integrator {
  readonly id = 'dopri5';
  readonly displayName = 'Dormand-Prince 5(4) (adaptive, dense output)';
  readonly order = 5;
  readonly adaptive = true;
  readonly denseOutputOrder = 4;
  readonly description =
    'Embedded Dormand-Prince 5(4) pair, advancing with the fifth-order solution. FSAL: ' +
    'six derivative evaluations per accepted step. Fourth-order continuous extension ' +
    "(Shampine) used for exact event location. Step control uses separate absolute and " +
    'relative tolerances.';

  readonly tolerance: ErrorTolerance;
  private readonly safety: number;
  private readonly maxGrowth: number;
  private readonly minShrink: number;

  private readonly k: Float64Array[];
  private readonly scratch: Float64Array;
  private readonly candidate: Float64Array;
  private readonly error: Float64Array;

  // Derivative cache: the value of f at (cachedP, cachedY), reused when the next step
  // starts from exactly that point. Covers both FSAL after an acceptance and the retry
  // after a rejection, and is invalidated automatically if the caller changes y.
  private cachedF: DerivativeFn | undefined;
  private cachedP = Number.NaN;
  private readonly cachedY: Float64Array;
  private readonly cachedK: Float64Array;

  // Dense output of the last accepted step.
  private readonly yOld: Float64Array;
  private readonly denseWeights = new Float64Array(7);
  private hOld = Number.NaN;
  private denseValid = false;

  constructor(dimension: number, options: Dopri5Options = {}) {
    this.tolerance = options.tolerance ?? DEFAULT_ERROR_TOLERANCE;
    this.safety = options.safety ?? 0.9;
    this.maxGrowth = options.maxGrowth ?? 10;
    this.minShrink = options.minShrink ?? 0.2;

    this.k = Array.from({ length: 7 }, () => new Float64Array(dimension));
    this.scratch = new Float64Array(dimension);
    this.candidate = new Float64Array(dimension);
    this.error = new Float64Array(dimension);
    this.cachedY = new Float64Array(dimension);
    this.cachedK = new Float64Array(dimension);
    this.yOld = new Float64Array(dimension);
  }

  private cacheMatches(f: DerivativeFn, p: number, y: Float64Array): boolean {
    if (f !== this.cachedF || p !== this.cachedP) return false;
    for (let i = 0; i < y.length; i += 1) {
      if (y[i] !== this.cachedY[i]) return false;
    }
    return true;
  }

  step(f: DerivativeFn, p: number, y: Float64Array, h: number): StepResult {
    const [k1, k2, k3, k4, k5, k6, k7] = this.k;
    const { scratch, candidate, error } = this;
    const n = y.length;
    if (n !== k1.length) {
      throw new RangeError(`Dopri5Integrator was constructed for dimension ${k1.length}, received ${n}.`);
    }

    // A new step invalidates the previous step's dense output.
    this.denseValid = false;

    if (this.cacheMatches(f, p, y)) {
      k1.set(this.cachedK);
    } else {
      f(p, y, k1);
      this.cachedF = f;
      this.cachedP = p;
      this.cachedY.set(y);
      this.cachedK.set(k1);
    }

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
    f(p + h, scratch, k6);

    for (let i = 0; i < n; i += 1) {
      candidate[i] =
        y[i] + h * (B1 * k1[i] + B3 * k3[i] + B4 * k4[i] + B5 * k5[i] + B6 * k6[i]);
    }

    // FSAL stage: the derivative at the candidate point. Needed for the error estimate,
    // and becomes k1 of the next step if this one is accepted.
    f(p + h, candidate, k7);

    for (let i = 0; i < n; i += 1) {
      error[i] =
        h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k6[i] + E7 * k7[i]);
    }

    const errorNorm = scaledErrorNorm(error, y, candidate, this.tolerance);

    if (!Number.isFinite(errorNorm)) {
      return { accepted: false, usedStep: h, nextStep: h * this.minShrink, errorNorm };
    }

    const accepted = errorNorm <= 1;
    // Exponent -1/5 = -1/(q+1) with q = 4 the order of the embedded error estimate.
    let factor =
      errorNorm === 0 ? this.maxGrowth : this.safety * Math.pow(errorNorm, -0.2);
    factor = Math.min(this.maxGrowth, Math.max(this.minShrink, factor));

    if (accepted) {
      this.yOld.set(y);
      this.hOld = h;
      this.denseValid = true;
      y.set(candidate);

      this.cachedF = f;
      this.cachedP = p + h;
      this.cachedY.set(candidate);
      this.cachedK.set(k7);
    }

    return { accepted, usedStep: h, nextStep: h * factor, errorNorm };
  }

  interpolate(theta: number, out: Float64Array): void {
    if (!this.denseValid) {
      throw new Error(
        'Dopri5Integrator.interpolate: there is no accepted step to interpolate over.',
      );
    }
    const n = out.length;
    const t1 = theta;
    const t2 = t1 * theta;
    const t3 = t2 * theta;
    const t4 = t3 * theta;
    const h = this.hOld;

    // Stage weights b_i(theta) = sum_j P[i][j] theta^(j+1).
    const weights = this.denseWeights;
    for (let i = 0; i < 7; i += 1) {
      const row = P[i];
      weights[i] = row[0] * t1 + row[1] * t2 + row[2] * t3 + row[3] * t4;
    }

    for (let j = 0; j < n; j += 1) {
      let sum = 0;
      for (let i = 0; i < 7; i += 1) {
        const w = weights[i];
        if (w !== 0) sum += w * this.k[i][j];
      }
      out[j] = this.yOld[j] + h * sum;
    }
  }
}
