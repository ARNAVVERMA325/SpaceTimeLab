import type { DerivativeFn } from '../geodesic-system.js';
import type { Integrator, StepResult } from './integrator.js';

/**
 * Classical fixed-step fourth-order Runge-Kutta (ROADMAP.md 1.3).
 *
 * Butcher tableau:
 *
 *   0   |
 *   1/2 | 1/2
 *   1/2 | 0    1/2
 *   1   | 0    0    1
 *   ----+---------------------
 *       | 1/6  1/3  1/3  1/6
 *
 * Fourth-order accurate locally to O(h^5) and globally to O(h^4). It carries no error
 * estimate, so `errorNorm` is NaN and every step is accepted: with RK4 the step size
 * is the caller's responsibility.
 *
 * RK4 is not symplectic. CLAUDE.md §7.2 notes that long-lived bound timelike orbits
 * are where that matters, and that a symplectic method should be considered for them
 * rather than mandated everywhere. For the finite ray-tracing integrations of M1 and
 * M2, the relevant control is local error, not long-time phase-space structure.
 */
export class RK4Integrator implements Integrator {
  readonly id = 'rk4';
  readonly displayName = 'Runge-Kutta 4 (fixed step)';
  readonly order = 4;
  readonly adaptive = false;
  readonly description =
    'Classical fixed-step RK4. Global error O(h^4). No local error estimate; step size ' +
    'is chosen by the caller. Not symplectic.';

  private readonly k1: Float64Array;
  private readonly k2: Float64Array;
  private readonly k3: Float64Array;
  private readonly k4: Float64Array;
  private readonly scratch: Float64Array;

  constructor(dimension: number) {
    this.k1 = new Float64Array(dimension);
    this.k2 = new Float64Array(dimension);
    this.k3 = new Float64Array(dimension);
    this.k4 = new Float64Array(dimension);
    this.scratch = new Float64Array(dimension);
  }

  step(f: DerivativeFn, p: number, y: Float64Array, h: number): StepResult {
    const { k1, k2, k3, k4, scratch } = this;
    const n = y.length;
    if (n !== k1.length) {
      throw new RangeError(`RK4Integrator was constructed for dimension ${k1.length}, received ${n}.`);
    }

    f(p, y, k1);

    for (let i = 0; i < n; i += 1) scratch[i] = y[i] + 0.5 * h * k1[i];
    f(p + 0.5 * h, scratch, k2);

    for (let i = 0; i < n; i += 1) scratch[i] = y[i] + 0.5 * h * k2[i];
    f(p + 0.5 * h, scratch, k3);

    for (let i = 0; i < n; i += 1) scratch[i] = y[i] + h * k3[i];
    f(p + h, scratch, k4);

    for (let i = 0; i < n; i += 1) {
      y[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }

    return { accepted: true, usedStep: h, nextStep: h, errorNorm: Number.NaN };
  }
}
