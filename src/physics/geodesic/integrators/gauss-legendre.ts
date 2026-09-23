import type { DerivativeFn } from '../geodesic-system.js';
import type { Integrator, StepResult } from './integrator.js';

/**
 * Two-stage Gauss-Legendre collocation: implicit, fourth order, and symplectic
 * (CLAUDE.md §7.2).
 *
 *   1/2 - sqrt(3)/6 | 1/4               1/4 - sqrt(3)/6
 *   1/2 + sqrt(3)/6 | 1/4 + sqrt(3)/6   1/4
 *   ----------------+-----------------------------------
 *                   | 1/2               1/2
 *
 * Why it exists alongside Dormand-Prince. Explicit Runge-Kutta methods are not
 * symplectic, so over a long-lived bound orbit their energy error accumulates secularly.
 * Gauss-Legendre methods are symplectic for any Hamiltonian, separable or not (Sanz-Serna
 * 1988), which is what the geodesic Hamiltonian H = 1/2 g^{mu nu} p_mu p_nu requires. By
 * backward error analysis the numerical flow exactly conserves a modified Hamiltonian
 * H + O(h^4) for exponentially long times, so the mass-shell error stays bounded instead
 * of growing. CLAUDE.md §7.2 names long-lived bound timelike orbits as the case for this.
 *
 * CLAUDE.md §7.2 also warns about the cost: the stages are defined implicitly. They are
 * solved here by fixed-point iteration to near machine precision — anything looser would
 * leave a non-symplectic residual that accumulates like any other error. If the
 * iteration does not converge the step is rejected and the driver shrinks it, rather than
 * accepting an unconverged, non-symplectic step.
 */

const SQRT3_6 = Math.sqrt(3) / 6;
const C1 = 0.5 - SQRT3_6;
const C2 = 0.5 + SQRT3_6;
const A11 = 0.25;
const A12 = 0.25 - SQRT3_6;
const A21 = 0.25 + SQRT3_6;
const A22 = 0.25;

export interface GaussLegendreOptions {
  /** Relative convergence threshold on the stage derivatives. */
  readonly iterationTolerance?: number;
  readonly maxIterations?: number;
}

export class GaussLegendre4Integrator implements Integrator {
  readonly id = 'gauss-legendre-4';
  readonly displayName = 'Gauss-Legendre 2-stage (implicit, symplectic, fixed step)';
  readonly order = 4;
  readonly adaptive = false;
  readonly symplectic = true;
  readonly description =
    'Implicit fourth-order Gauss-Legendre collocation. Symplectic for any Hamiltonian, so ' +
    'energy-like errors stay bounded over long integrations instead of drifting. Stages ' +
    'solved by fixed-point iteration to near machine precision; unconverged steps are ' +
    'rejected, not accepted.';

  private readonly iterationTolerance: number;
  private readonly maxIterations: number;
  private readonly k1: Float64Array;
  private readonly k2: Float64Array;
  private readonly k1Next: Float64Array;
  private readonly k2Next: Float64Array;
  private readonly stage: Float64Array;

  /** Fixed-point iterations used by the last step, for diagnostics. */
  lastIterations = 0;

  constructor(dimension: number, options: GaussLegendreOptions = {}) {
    this.iterationTolerance = options.iterationTolerance ?? 4 * Number.EPSILON;
    this.maxIterations = options.maxIterations ?? 100;
    this.k1 = new Float64Array(dimension);
    this.k2 = new Float64Array(dimension);
    this.k1Next = new Float64Array(dimension);
    this.k2Next = new Float64Array(dimension);
    this.stage = new Float64Array(dimension);
  }

  step(f: DerivativeFn, p: number, y: Float64Array, h: number): StepResult {
    const { k1, k2, k1Next, k2Next, stage } = this;
    const n = y.length;
    if (n !== k1.length) {
      throw new RangeError(`GaussLegendre4Integrator was constructed for dimension ${k1.length}, received ${n}.`);
    }

    // Initial guess: both stages at the explicit derivative.
    f(p, y, k1);
    k2.set(k1);

    let converged = false;
    for (let iteration = 1; iteration <= this.maxIterations; iteration += 1) {
      for (let i = 0; i < n; i += 1) stage[i] = y[i] + h * (A11 * k1[i] + A12 * k2[i]);
      f(p + C1 * h, stage, k1Next);
      for (let i = 0; i < n; i += 1) stage[i] = y[i] + h * (A21 * k1[i] + A22 * k2[i]);
      f(p + C2 * h, stage, k2Next);

      let change = 0;
      let scale = 0;
      for (let i = 0; i < n; i += 1) {
        change = Math.max(change, Math.abs(k1Next[i] - k1[i]), Math.abs(k2Next[i] - k2[i]));
        scale = Math.max(scale, Math.abs(k1Next[i]), Math.abs(k2Next[i]));
      }
      k1.set(k1Next);
      k2.set(k2Next);
      this.lastIterations = iteration;

      if (!Number.isFinite(change)) break;
      if (change <= this.iterationTolerance * Math.max(scale, 1e-300)) {
        converged = true;
        break;
      }
    }

    if (!converged) {
      return { accepted: false, usedStep: h, nextStep: h / 2, errorNorm: Number.POSITIVE_INFINITY };
    }

    for (let i = 0; i < n; i += 1) y[i] += h * 0.5 * (k1[i] + k2[i]);
    return { accepted: true, usedStep: h, nextStep: h, errorNorm: Number.NaN };
  }
}
