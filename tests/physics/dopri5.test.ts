import { describe, expect, it } from 'vitest';
import type { DerivativeFn } from '../../src/physics/geodesic/geodesic-system.js';
import { brentRoot, integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { RK4Integrator } from '../../src/physics/geodesic/integrators/rk4.js';
import { RKF45Integrator } from '../../src/physics/geodesic/integrators/rkf45.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { schwarzschild } from '../../src/physics/spacetimes/schwarzschild.js';
import {
  deflectionAngleExact,
  impactParameterForTurningPoint,
} from '../../src/physics/spacetimes/schwarzschild-analytic.js';
import {
  equatorialNullRay,
  measureDeflection,
} from '../../src/physics/spacetimes/schwarzschild-rays.js';

/**
 * Dormand-Prince 5(4), its continuous extension, and exact event location.
 *
 * Every claim in the integrator's documentation is measured here: the order of the
 * propagated solution, the order of the dense output, the FSAL evaluation count, and
 * the efficiency advantage over Fehlberg that justified replacing it.
 */

/** y'' = -y as y = [q, p]; exact solution q = cos t, p = -sin t from (1, 0). */
const oscillator: DerivativeFn = (_t, y, d) => {
  d[0] = y[1];
  d[1] = -y[0];
};

/**
 * A DOPRI5 whose error control can never reject, so a step size is imposed exactly.
 * Used only to measure order; never for real integrations.
 */
function forcedFixedStep(): Dopri5Integrator {
  return new Dopri5Integrator(2, { tolerance: { absolute: 1e6, relative: 1e6 } });
}

describe('Dormand-Prince 5(4)', () => {
  it('converges at fifth order', () => {
    const T = 2 * Math.PI;
    const errorFor = (n: number): number => {
      const integrator = forcedFixedStep();
      const y = new Float64Array([1, 0]);
      const h = T / n;
      for (let i = 0; i < n; i += 1) {
        expect(integrator.step(oscillator, i * h, y, h).accepted).toBe(true);
      }
      return Math.max(Math.abs(y[0] - 1), Math.abs(y[1]));
    };
    const errors = [10, 20, 40, 80].map(errorFor);
    // Halving h must cut the global error by about 2^5 = 32. Measured 34.9, 29.3, 31.3.
    for (let i = 1; i < errors.length; i += 1) {
      const ratio = errors[i - 1] / errors[i];
      expect(ratio).toBeGreaterThan(26);
      expect(ratio).toBeLessThan(40);
    }
  });

  it('has a continuous extension with O(h^5) local error', () => {
    for (const theta of [0.3, 0.5, 0.8]) {
      const errors = [0.4, 0.2, 0.1, 0.05].map((h) => {
        const integrator = forcedFixedStep();
        const y = new Float64Array([1, 0]);
        integrator.step(oscillator, 0, y, h);
        const out = new Float64Array(2);
        integrator.interpolate(theta, out);
        return Math.max(Math.abs(out[0] - Math.cos(theta * h)), Math.abs(out[1] + Math.sin(theta * h)));
      });
      // Fourth-order interpolant, so local error O(h^5): ratio -> 32 as h shrinks.
      const finest = errors[2] / errors[3];
      expect(finest, `theta = ${theta}`).toBeGreaterThan(28);
      expect(finest, `theta = ${theta}`).toBeLessThan(36);
    }
  });

  it('interpolates exactly through both ends of the step', () => {
    const integrator = forcedFixedStep();
    const y = new Float64Array([1, 0]);
    integrator.step(oscillator, 0, y, 0.3);
    const start = new Float64Array(2);
    const end = new Float64Array(2);
    integrator.interpolate(0, start);
    integrator.interpolate(1, end);
    expect(Array.from(start)).toEqual([1, 0]);
    // At theta = 1 each dense-output row sums to its weight b_i.
    expect(Math.abs(end[0] - y[0])).toBeLessThan(1e-15);
    expect(Math.abs(end[1] - y[1])).toBeLessThan(1e-15);
  });

  it('refuses to interpolate before any step has been accepted', () => {
    expect(() => forcedFixedStep().interpolate(0.5, new Float64Array(2))).toThrow(/no accepted step/);
  });

  it('spends six derivative evaluations per accepted step (FSAL)', () => {
    let evaluations = 0;
    const counted: DerivativeFn = (t, y, d) => {
      evaluations += 1;
      oscillator(t, y, d);
    };
    const integrator = forcedFixedStep();
    const y = new Float64Array([1, 0]);
    const steps = 50;
    // Advance the parameter the way a driver does, by accumulation. The FSAL cache is
    // keyed on the exact parameter value, and i * h differs from the accumulated sum by
    // an ulp often enough to defeat it — correctly, since it cannot know f ignores t.
    let t = 0;
    for (let i = 0; i < steps; i += 1) {
      integrator.step(counted, t, y, 0.05);
      t += 0.05;
    }
    // Seven on the first step, six thereafter: the seventh stage is f at the new point
    // and becomes the first stage of the next step.
    expect(evaluations).toBe(7 + 6 * (steps - 1));
  });

  it('recomputes the first stage when the caller changes the state between steps', () => {
    // The FSAL cache must never be reused for a point it was not computed at.
    const integrator = forcedFixedStep();
    const y = new Float64Array([1, 0]);
    integrator.step(oscillator, 0, y, 0.1);
    y[0] += 1e-3;
    const withCache = Float64Array.from(y);
    integrator.step(oscillator, 0.1, withCache, 0.1);

    const fresh = forcedFixedStep();
    const reference = Float64Array.from(y);
    fresh.step(oscillator, 0.1, reference, 0.1);
    expect(Array.from(withCache)).toEqual(Array.from(reference));
  });

  it('beats Fehlberg 4(5) on both accuracy and cost for a Schwarzschild deflection', () => {
    // The justification for making it the default, measured rather than asserted. At
    // every tolerance tried, DOPRI5 reached a smaller error for fewer steps.
    const model = schwarzschild(1);
    const b = impactParameterForTurningPoint(1, 10);
    const exact = deflectionAngleExact(1, 10);
    for (const tolerance of [1e-8, 1e-10, 1e-12]) {
      const run = (integrator: Dopri5Integrator | RKF45Integrator) =>
        measureDeflection({ model, integrator, impactParameter: b, startRadius: 1000 });
      const dp = run(new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: tolerance, relative: tolerance } }));
      const rkf = run(new RKF45Integrator(STATE_DIM, { tolerance: { absolute: tolerance, relative: tolerance } }));
      const errorDp = Math.abs(dp.deflectionAngle! - exact) / exact;
      const errorRkf = Math.abs(rkf.deflectionAngle! - exact) / exact;
      expect(errorDp, `tolerance ${tolerance}`).toBeLessThan(errorRkf);
      expect(dp.steps + dp.rejectedSteps).toBeLessThanOrEqual(rkf.steps + rkf.rejectedSteps);
    }
  });
});

describe('convergence in curved spacetime (CLAUDE.md §16)', () => {
  it('shows fourth-order convergence of RK4 on a Schwarzschild deflection', () => {
    // The harmonic-oscillator order test checks the integrator in isolation. This checks
    // it on the actual problem: a strongly bent null geodesic, with the error measured
    // against the exact quadrature. Measured ratios 15.8 and 15.9 per halving.
    const model = schwarzschild(1);
    const b = impactParameterForTurningPoint(1, 10);
    const exact = deflectionAngleExact(1, 10);
    const errors = [0.4, 0.2, 0.1].map((h) => {
      const m = measureDeflection({
        model,
        integrator: new RK4Integrator(STATE_DIM),
        impactParameter: b,
        startRadius: 200,
        limits: { initialStep: h },
      });
      return Math.abs(m.deflectionAngle! - exact);
    });
    for (let i = 1; i < errors.length; i += 1) {
      const ratio = errors[i - 1] / errors[i];
      expect(ratio).toBeGreaterThan(14);
      expect(ratio).toBeLessThan(18);
    }
  });

  it('shows the adaptive error falling with the requested tolerance', () => {
    const model = schwarzschild(1);
    const b = impactParameterForTurningPoint(1, 10);
    const exact = deflectionAngleExact(1, 10);
    let previous = Infinity;
    for (const tolerance of [1e-6, 1e-8, 1e-10, 1e-12]) {
      const m = measureDeflection({
        model,
        integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: tolerance, relative: tolerance } }),
        impactParameter: b,
        startRadius: 1000,
      });
      const error = Math.abs(m.deflectionAngle! - exact) / exact;
      expect(error, `tolerance ${tolerance}`).toBeLessThan(previous);
      previous = error;
    }
  });
});

describe('event location', () => {
  it("finds roots with Brent's method to near machine precision", () => {
    const root = brentRoot((x) => Math.cos(x) - x, 0, 1, 1, Math.cos(1) - 1);
    expect(Math.abs(Math.cos(root) - root)).toBeLessThan(1e-15);
    expect(brentRoot((x) => x * x * x - 2, 0, 2, -2, 6)).toBeCloseTo(Math.cbrt(2), 14);
    expect(() => brentRoot((x) => x, 1, 2, 1, 2)).toThrow(/does not straddle/);
  });

  it('lands a terminal event exactly on its surface instead of overshooting', () => {
    const model = schwarzschild(1);
    for (const integrator of [
      new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-10, relative: 1e-10 } }),
      new RK4Integrator(STATE_DIM),
    ]) {
      const result = integrateGeodesic({
        model,
        integrator,
        initial: equatorialNullRay(model, 30, 12, 'outgoing'),
        limits: { initialStep: 0.5, parameterMax: 1e4, maxSteps: 100_000, maxStep: 50 },
        events: [{ id: 'sphere', value: (x) => x[1] - 500, direction: 1, terminal: true }],
      });
      expect(result.reason, integrator.id).toBe('event');
      // Steps of up to 50 would overshoot by as much; the located state sits on r = 500
      // to rounding, independent of the integration tolerance.
      expect(Math.abs(result.final.position_x[1] - 500), integrator.id).toBeLessThan(1e-11);
      expect(result.events.at(-1)?.id).toBe('sphere');
    }
  });

  it('records non-terminal crossings and carries on', () => {
    // A ray that passes periapsis and then escapes: the turning point is recorded, and
    // its radius equals the analytic closest approach.
    const model = schwarzschild(1);
    const r0 = 8;
    const b = impactParameterForTurningPoint(1, r0);
    const result = integrateGeodesic({
      model,
      integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-12, relative: 1e-12 } }),
      initial: equatorialNullRay(model, 100, b, 'ingoing'),
      limits: { initialStep: 1e-3, parameterMax: 1e4, maxSteps: 100_000 },
      events: [
        { id: 'periapsis', value: (_x, k) => k[1], direction: 1, terminal: false },
        { id: 'out', value: (x) => x[1] - 100, direction: 1, terminal: true },
      ],
    });
    expect(result.events.map((e) => e.id)).toEqual(['periapsis', 'out']);
    const periapsis = result.events[0].state.position_x[1];
    expect(Math.abs(periapsis - r0) / r0).toBeLessThan(1e-9);
  });

  it('respects event direction', () => {
    const model = schwarzschild(1);
    const result = integrateGeodesic({
      model,
      integrator: new Dopri5Integrator(STATE_DIM),
      initial: equatorialNullRay(model, 30, 12, 'outgoing'),
      limits: { initialStep: 0.5, parameterMax: 200, maxSteps: 100_000 },
      // An outgoing ray crosses r = 60 upward only; a downward-only event never fires.
      events: [{ id: 'down-only', value: (x) => x[1] - 60, direction: -1, terminal: true }],
    });
    expect(result.reason).toBe('parameter-limit');
    expect(result.events).toHaveLength(0);
  });
});
