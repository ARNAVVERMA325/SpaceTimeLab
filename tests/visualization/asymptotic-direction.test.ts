import { describe, expect, it } from 'vitest';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { generateNullRay, staticObserver } from '../../src/physics/observer/observer.js';
import { schwarzschild } from '../../src/physics/spacetimes/schwarzschild.js';
import { asymptoticDirection, isCaptured } from '../../src/physics/spacetimes/schwarzschild-rays.js';
import type { CartesianVec3 } from '../../src/physics/spacetimes/spacetime-model.js';
import { DEFAULT_CELESTIAL_GRID } from '../../src/visualization/celestial-grid.js';
import { traceRay, type TraceConfig } from '../../src/visualization/raytracer.js';

/**
 * Removing the finite-radius background bias.
 *
 * A ray reaching the background sphere at radius R is still being bent, so reading the
 * sky off its local direction there samples the wrong point by the deflection still to
 * come — measured at 2.5e-4 rad for b = 10M at R = 200M. The correction replaces that
 * with the exact direction at infinity.
 *
 * The decisive test is invariance: the corrected direction must not depend on which R
 * the trace happened to stop at, while the uncorrected one visibly does.
 */

const M = 1;
const model = schwarzschild(M);
const observer = staticObserver(model, [0, 20, Math.PI / 2, 0]);

function unit(v: readonly number[]): CartesianVec3 {
  const n = Math.hypot(v[0], v[1], v[2]);
  return [v[0] / n, v[1] / n, v[2] / n];
}

function chord(a: CartesianVec3, b: CartesianVec3): number {
  const ua = unit(a);
  const ub = unit(b);
  return Math.hypot(ua[0] - ub[0], ua[1] - ub[1], ua[2] - ub[2]);
}

function configFor(radius: number, corrected: boolean): TraceConfig {
  return {
    model,
    integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-13, relative: 1e-13 } }),
    observer,
    grid: { ...DEFAULT_CELESTIAL_GRID, radius },
    limits: { initialStep: 1e-3, parameterMax: 1e7, maxSteps: 2_000_000 },
    captureTest: (x, k) => isCaptured(model, x, k),
    orbitalPlaneReduction: true,
    ...(corrected ? { asymptoticDirection: (p: CartesianVec3, d: CartesianVec3) => asymptoticDirection(model, p, d) } : {}),
  };
}

// Local directions whose rays escape (their impact parameter clears b_c from r = 20M).
const DIRECTIONS: readonly CartesianVec3[] = [
  unit([-1, 0.3, 0.45]),
  unit([-0.8, -0.5, 0.3]),
  unit([-0.4, 0.2, -0.9]),
];

describe('asymptotic direction correction', () => {
  it('makes the sampled direction independent of the background radius', () => {
    for (const local of DIRECTIONS) {
      const ray = generateNullRay(observer, local);
      const corrected = [100, 200, 1000, 5000].map((R) => traceRay(configFor(R, true), ray));
      const uncorrected = [100, 200, 1000, 5000].map((R) => traceRay(configFor(R, false), ray));
      for (const result of [...corrected, ...uncorrected]) expect(result.outcome).toBe('background');

      // Corrected: all four agree to integration accuracy.
      for (let i = 1; i < corrected.length; i += 1) {
        expect(chord(corrected[0].exitDirection, corrected[i].exitDirection), `${local}`).toBeLessThan(1e-10);
      }
      // Uncorrected: R = 100M and R = 5000M differ by the residual bending, which is
      // orders of magnitude larger than the corrected spread.
      expect(chord(uncorrected[0].exitDirection, uncorrected[3].exitDirection)).toBeGreaterThan(1e-5);
    }
  });

  it('agrees with tracing essentially to infinity and not correcting at all', () => {
    // An independent route to the same answer: at R = 10^7 M the residual bending is
    // about M b / R^2 ~ 1e-13 and the local direction already is the asymptotic one.
    for (const local of DIRECTIONS) {
      const ray = generateNullRay(observer, local);
      const near = traceRay(configFor(200, true), ray);
      const far = traceRay(configFor(1e7, false), ray);
      expect(chord(near.exitDirection, far.exitDirection), `${local}`).toBeLessThan(1e-9);
    }
  });

  it('shrinks the uncorrected bias as M b / R^2, as the analysis predicts', () => {
    const ray = generateNullRay(observer, DIRECTIONS[0]);
    const truth = traceRay(configFor(200, true), ray).exitDirection;
    const bias = [100, 200, 400].map((R) => chord(traceRay(configFor(R, false), ray).exitDirection, truth));
    // Doubling R quarters the bias.
    expect(bias[0] / bias[1]).toBeGreaterThan(3.6);
    expect(bias[0] / bias[1]).toBeLessThan(4.4);
    expect(bias[1] / bias[2]).toBeGreaterThan(3.6);
    expect(bias[1] / bias[2]).toBeLessThan(4.4);
  });

  it('reduces to the local direction in flat space', () => {
    const weak = schwarzschild(1e-12);
    const position: CartesianVec3 = [300, 40, -20];
    const direction: CartesianVec3 = unit([0.9, 0.3, 0.1]);
    const corrected = asymptoticDirection(weak, position, direction);
    expect(chord(corrected, direction)).toBeLessThan(1e-12);
  });

  it('refuses a ray that is not moving outward', () => {
    expect(() => asymptoticDirection(model, [100, 0, 0], [-1, 0.1, 0])).toThrow(/not moving outward/);
  });
});
