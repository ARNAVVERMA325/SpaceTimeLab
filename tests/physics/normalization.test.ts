import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import {
  four_velocity_u,
  null_wavevector_k,
  nullState,
  parameterName,
  timelikeState,
} from '../../src/physics/core/phase-space.js';
import { RK4Integrator } from '../../src/physics/geodesic/integrators/rk4.js';
import { RKF45Integrator } from '../../src/physics/geodesic/integrators/rkf45.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { minkowski } from '../../src/physics/spacetimes/minkowski.js';
import {
  normalizationInvariant,
  normalizationResidual,
  normalizeTimelike,
} from '../../src/physics/validation/normalization.js';
import {
  NORMALIZATION_DRIFT_LONG_RUN,
  NULL_NORMALIZATION_POINTWISE,
  TIMELIKE_NORMALIZATION_POINTWISE,
  checkTolerance,
} from '../../src/physics/validation/tolerances.js';

/**
 * ROADMAP.md 1.4 — geodesic normalization (CLAUDE.md §16).
 *
 *   Null:     g_mu_nu k^mu k^nu = 0
 *   Timelike: g_mu_nu u^mu u^nu = -1
 */

/** A null wavevector in Minkowski: |k_spatial| = k^0. */
function nullRay(direction: readonly [number, number, number]): Vec4 {
  const norm = Math.hypot(...direction);
  return [1, direction[0] / norm, direction[1] / norm, direction[2] / norm];
}

/** A four-velocity for a boost of speed v along x. */
function boostedFourVelocity(v: number): Vec4 {
  const gamma = 1 / Math.sqrt(1 - v * v);
  return [gamma, gamma * v, 0, 0];
}

describe('geodesic normalization (ROADMAP.md 1.4, CLAUDE.md §16)', () => {
  it('reports g_mu_nu k^mu k^nu = 0 for null wavevectors', () => {
    const directions: readonly (readonly [number, number, number])[] = [
      [1, 0, 0],
      [0, 1, 0],
      [1, 1, 1],
      [-3, 7, 0.5],
    ];
    for (const d of directions) {
      const state = nullState([0, 0, 0, 0], nullRay(d));
      const check = checkTolerance(
        NULL_NORMALIZATION_POINTWISE,
        normalizationResidual(minkowski, state),
      );
      expect(check.withinTolerance, `${check.message} for direction ${JSON.stringify(d)}`).toBe(true);
    }
  });

  it('reports g_mu_nu u^mu u^nu = -1 for four-velocities', () => {
    for (const v of [0, 0.1, 0.5, 0.9, 0.99]) {
      const state = timelikeState([0, 0, 0, 0], boostedFourVelocity(v));
      const invariant = normalizationInvariant(minkowski, state);
      const check = checkTolerance(TIMELIKE_NORMALIZATION_POINTWISE, invariant - -1, -1);
      expect(check.withinTolerance, `${check.message} at v = ${v}`).toBe(true);
    }
  });

  it('normalizes an un-normalized timelike tangent to exactly -1', () => {
    const raw = timelikeState([1, 2, 3, 4], [5, 1, 0, 0]);
    const normalized = normalizeTimelike(minkowski, raw);
    const check = checkTolerance(
      TIMELIKE_NORMALIZATION_POINTWISE,
      normalizationInvariant(minkowski, normalized) - -1,
      -1,
    );
    expect(check.withinTolerance, check.message).toBe(true);
  });

  it('refuses to normalize a null or spacelike tangent to -1', () => {
    const nullish = timelikeState([0, 0, 0, 0], [1, 1, 0, 0]);
    expect(() => normalizeTimelike(minkowski, nullish)).toThrow(RangeError);

    const spacelike = timelikeState([0, 0, 0, 0], [0, 1, 0, 0]);
    expect(() => normalizeTimelike(minkowski, spacelike)).toThrow(RangeError);
  });

  it('refuses to call a null tangent a four-velocity (CLAUDE.md §3, §24)', () => {
    const photon = nullState([0, 0, 0, 0], nullRay([1, 0, 0]));
    expect(() => four_velocity_u(photon)).toThrow(TypeError);
    expect(null_wavevector_k(photon)).toEqual([1, 1, 0, 0]);

    const particle = timelikeState([0, 0, 0, 0], boostedFourVelocity(0.5));
    expect(() => null_wavevector_k(particle)).toThrow(TypeError);
    expect(four_velocity_u(particle)).toEqual(boostedFourVelocity(0.5));
  });

  it('names the integration parameter per worldline kind (CLAUDE.md §2)', () => {
    expect(parameterName('null')).toBe('lambda');
    expect(parameterName('timelike')).toBe('tau');
  });

  it.each([
    ['RK4', () => new RK4Integrator(STATE_DIM), { initialStep: 1e-2, maxSteps: 20_000 }],
    ['RKF45', () => new RKF45Integrator(STATE_DIM), { initialStep: 1e-2, maxSteps: 20_000 }],
  ])('preserves the null invariant along an integrated geodesic (%s)', (_name, makeIntegrator, limits) => {
    const initial = nullState([0, 0, 0, 0], nullRay([1, 2, 3]));
    const result = integrateGeodesic({
      model: minkowski,
      integrator: makeIntegrator(),
      initial,
      limits: { ...limits, parameterMax: 100 },
    });

    expect(result.reason).toBe('parameter-limit');
    const check = checkTolerance(
      NORMALIZATION_DRIFT_LONG_RUN,
      normalizationResidual(minkowski, result.final),
    );
    expect(check.withinTolerance, check.message).toBe(true);
  });

  it.each([
    ['RK4', () => new RK4Integrator(STATE_DIM)],
    ['RKF45', () => new RKF45Integrator(STATE_DIM)],
  ])('preserves the timelike invariant along an integrated worldline (%s)', (_name, makeIntegrator) => {
    const initial = timelikeState([0, 0, 0, 0], boostedFourVelocity(0.8));
    const result = integrateGeodesic({
      model: minkowski,
      integrator: makeIntegrator(),
      initial,
      limits: { initialStep: 1e-2, parameterMax: 100, maxSteps: 20_000 },
    });

    const check = checkTolerance(
      NORMALIZATION_DRIFT_LONG_RUN,
      normalizationResidual(minkowski, result.final),
    );
    expect(check.withinTolerance, check.message).toBe(true);
  });
});
