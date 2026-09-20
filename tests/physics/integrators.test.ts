import { describe, expect, it } from 'vitest';
import type { DerivativeFn } from '../../src/physics/geodesic/geodesic-system.js';
import { RK4Integrator } from '../../src/physics/geodesic/integrators/rk4.js';
import { RKF45Integrator } from '../../src/physics/geodesic/integrators/rkf45.js';
import { scaledErrorNorm } from '../../src/physics/geodesic/integrators/integrator.js';

/**
 * ROADMAP.md 1.3 — the integrator core, tested as a numerical method in its own right.
 *
 * These are unit tests of the integrators, not of general relativity. Minkowski has a
 * vanishing connection, so a geodesic there exercises none of the integrator's
 * machinery: the right-hand side returns a constant and any consistent method is exact.
 * A harmonic oscillator with a known closed-form solution is used instead, so that the
 * order-of-accuracy and error-control claims are actually tested before the integrators
 * are trusted with curved spacetimes in Milestone 2A.
 *
 * This is the convergence-test category CLAUDE.md §16 requires: changing the timestep
 * must produce predictable convergence behaviour.
 */

/** y'' = -y, written first-order as y = [q, p] with q' = p, p' = -q. */
const harmonicOscillator: DerivativeFn = (_t, y, dydt) => {
  dydt[0] = y[1];
  dydt[1] = -y[0];
};

/** Exact solution with q(0) = 1, p(0) = 0: q(t) = cos t, p(t) = -sin t. */
function exactHarmonic(t: number): [number, number] {
  return [Math.cos(t), -Math.sin(t)];
}

function integrateFixed(steps: number, tEnd: number): Float64Array {
  const integrator = new RK4Integrator(2);
  const y = new Float64Array([1, 0]);
  const h = tEnd / steps;
  for (let n = 0; n < steps; n += 1) {
    integrator.step(harmonicOscillator, n * h, y, h);
  }
  return y;
}

function maxError(y: Float64Array, tEnd: number): number {
  const [q, p] = exactHarmonic(tEnd);
  return Math.max(Math.abs(y[0] - q), Math.abs(y[1] - p));
}

describe('RK4 (ROADMAP.md 1.3)', () => {
  it('reports itself as fixed-step with no error estimate', () => {
    const integrator = new RK4Integrator(2);
    expect(integrator.adaptive).toBe(false);
    expect(integrator.order).toBe(4);
    const result = integrator.step(harmonicOscillator, 0, new Float64Array([1, 0]), 0.1);
    expect(result.accepted).toBe(true);
    expect(Number.isNaN(result.errorNorm)).toBe(true);
  });

  it('converges at fourth order when the timestep is halved', () => {
    const tEnd = 2 * Math.PI;
    // Step counts chosen so the error is well above the f64 roundoff floor at the
    // coarse end and still above it at the fine end, where the O(h^4) signal would
    // otherwise be swamped by accumulated rounding.
    const coarse = maxError(integrateFixed(40, tEnd), tEnd);
    const medium = maxError(integrateFixed(80, tEnd), tEnd);
    const fine = maxError(integrateFixed(160, tEnd), tEnd);

    const firstRatio = coarse / medium;
    const secondRatio = medium / fine;

    // Halving h must reduce the global error by ~2^4 = 16 for a fourth-order method.
    expect(firstRatio).toBeGreaterThan(14);
    expect(firstRatio).toBeLessThan(18);
    expect(secondRatio).toBeGreaterThan(14);
    expect(secondRatio).toBeLessThan(18);
  });

  it('rejects a state whose dimension does not match its scratch buffers', () => {
    const integrator = new RK4Integrator(2);
    expect(() => integrator.step(harmonicOscillator, 0, new Float64Array(8), 0.1)).toThrow(RangeError);
  });
});

describe('RKF45 (ROADMAP.md 1.3, CLAUDE.md §7.1)', () => {
  it('reports itself as adaptive with a fifth-order propagated solution', () => {
    const integrator = new RKF45Integrator(2);
    expect(integrator.adaptive).toBe(true);
    expect(integrator.order).toBe(5);
  });

  it('keeps the solution within the requested tolerance over a full period', () => {
    const tolerance = { absolute: 1e-10, relative: 1e-10 };
    const integrator = new RKF45Integrator(2, { tolerance });
    const y = new Float64Array([1, 0]);
    const tEnd = 2 * Math.PI;

    let t = 0;
    let h = 0.1;
    let accepted = 0;
    let rejected = 0;

    while (t < tEnd && accepted + rejected < 100_000) {
      const attempt = Math.min(h, tEnd - t);
      const result = integrator.step(harmonicOscillator, t, y, attempt);
      if (result.accepted) {
        t += attempt;
        accepted += 1;
      } else {
        rejected += 1;
      }
      h = result.nextStep;
    }

    expect(t).toBeCloseTo(tEnd, 12);
    // Local error control does not bound global error directly, but over one period of
    // a non-stiff oscillator the accumulated error stays within a small multiple of the
    // per-step tolerance. A loose factor is used deliberately rather than pretending
    // local control implies a global guarantee.
    expect(maxError(y, tEnd)).toBeLessThan(1e-8);
    expect(accepted).toBeGreaterThan(0);
  });

  it('takes larger steps at a looser tolerance', () => {
    const countSteps = (absolute: number, relative: number): number => {
      const integrator = new RKF45Integrator(2, { tolerance: { absolute, relative } });
      const y = new Float64Array([1, 0]);
      const tEnd = 2 * Math.PI;
      let t = 0;
      let h = 0.1;
      let accepted = 0;
      while (t < tEnd && accepted < 100_000) {
        const attempt = Math.min(h, tEnd - t);
        const result = integrator.step(harmonicOscillator, t, y, attempt);
        if (result.accepted) {
          t += attempt;
          accepted += 1;
        }
        h = result.nextStep;
      }
      return accepted;
    };

    expect(countSteps(1e-6, 1e-6)).toBeLessThan(countSteps(1e-12, 1e-12));
  });

  it('leaves the state untouched when it rejects a step', () => {
    // A tolerance far tighter than a large step can satisfy forces a rejection.
    const integrator = new RKF45Integrator(2, { tolerance: { absolute: 1e-18, relative: 1e-18 } });
    const y = new Float64Array([1, 0]);
    const before = Array.from(y);
    const result = integrator.step(harmonicOscillator, 0, y, 1.0);
    expect(result.accepted).toBe(false);
    expect(Array.from(y)).toEqual(before);
    expect(result.nextStep).toBeLessThan(1.0);
  });

  it('rejects and shrinks rather than propagating a non-finite stage', () => {
    const blowUp: DerivativeFn = (_t, _y, dydt) => {
      dydt[0] = Number.POSITIVE_INFINITY;
      dydt[1] = Number.NaN;
    };
    const integrator = new RKF45Integrator(2);
    const y = new Float64Array([1, 0]);
    const result = integrator.step(blowUp, 0, y, 0.1);
    expect(result.accepted).toBe(false);
    expect(Number.isFinite(result.errorNorm)).toBe(false);
    expect(Array.from(y)).toEqual([1, 0]);
    expect(result.nextStep).toBeLessThan(0.1);
  });
});

describe('scaled error norm (CLAUDE.md §17)', () => {
  it('returns exactly 1 when every component sits on its tolerance', () => {
    const tolerance = { absolute: 1e-6, relative: 0 };
    const error = new Float64Array([1e-6, 1e-6]);
    const y = new Float64Array([1, 1]);
    expect(scaledErrorNorm(error, y, y, tolerance)).toBeCloseTo(1, 15);
  });

  it('uses the relative term to scale with the magnitude of the solution', () => {
    const tolerance = { absolute: 0, relative: 1e-6 };
    const small = scaledErrorNorm(
      new Float64Array([1e-6, 0]),
      new Float64Array([1, 0]),
      new Float64Array([1, 0]),
      tolerance,
    );
    const large = scaledErrorNorm(
      new Float64Array([1e-6, 0]),
      new Float64Array([1000, 0]),
      new Float64Array([1000, 0]),
      tolerance,
    );
    // The same absolute error is a thousand times less significant against a
    // thousand-times-larger solution component.
    expect(small / large).toBeCloseTo(1000, 6);
  });
});
