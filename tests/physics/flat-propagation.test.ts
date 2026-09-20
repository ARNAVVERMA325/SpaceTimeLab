import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { nullState, timelikeState } from '../../src/physics/core/phase-space.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { RK4Integrator } from '../../src/physics/geodesic/integrators/rk4.js';
import { RKF45Integrator } from '../../src/physics/geodesic/integrators/rkf45.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { minkowski } from '../../src/physics/spacetimes/minkowski.js';
import { conservedQuantityDrift } from '../../src/physics/validation/conserved.js';
import { normalizationResidual } from '../../src/physics/validation/normalization.js';
import {
  CONSERVED_QUANTITY_DRIFT_FLAT,
  FLAT_PROPAGATION_STRAIGHTNESS,
  NORMALIZATION_DRIFT_LONG_RUN,
  checkTolerance,
} from '../../src/physics/validation/tolerances.js';

/**
 * ROADMAP.md 1.4 — straight-line propagation in flat space held across 10^6 steps.
 *
 * This is the gate the roadmap names. In Minkowski the Christoffel symbols vanish, so
 * the geodesic equation reduces to dt^mu/dparam = 0 and the exact solution is
 *
 *   x^mu(param) = x^mu(0) + t^mu param
 *
 * Any deviation is accumulated floating-point error in the integrator and the driver,
 * with no physics mixed in — which is exactly what makes this a clean plumbing test.
 */

const STEPS = 1_000_000;
const STEP_SIZE = 1e-3;

/** Analytic straight-line position after `param` of affine parameter. */
function exactPosition(x0: Vec4, tangent: Vec4, param: number): Vec4 {
  return [
    x0[0] + tangent[0] * param,
    x0[1] + tangent[1] * param,
    x0[2] + tangent[2] * param,
    x0[3] + tangent[3] * param,
  ];
}

function maxComponentDeviation(a: Vec4, b: Vec4): number {
  let worst = 0;
  for (let mu = 0; mu < 4; mu += 1) {
    const d = Math.abs(a[mu] - b[mu]);
    if (d > worst) worst = d;
  }
  return worst;
}

describe('flat-space propagation across 10^6 steps (ROADMAP.md 1.4)', () => {
  it('holds a null geodesic straight, and conserves E and L_z exactly', () => {
    // Offset from the spatial origin so that angular_momentum_Lz is non-zero and its
    // relative drift is actually defined. A ray through the origin has L_z = 0, where a
    // relative measure would be undefined.
    const x0: Vec4 = [0, 0, 5, 0];
    const null_wavevector_k: Vec4 = [1, 1, 0, 0];
    const initial = nullState(x0, null_wavevector_k);

    const result = integrateGeodesic({
      model: minkowski,
      integrator: new RK4Integrator(STATE_DIM),
      initial,
      limits: {
        initialStep: STEP_SIZE,
        parameterMax: STEPS * STEP_SIZE,
        maxSteps: STEPS + 10,
      },
    });

    expect(result.reason).toBe('parameter-limit');
    expect(result.steps).toBeGreaterThanOrEqual(STEPS);

    const lambda = result.final.parameter;
    const expected = exactPosition(x0, null_wavevector_k, lambda);
    const deviation = maxComponentDeviation(result.final.position_x, expected);
    const distance = Math.hypot(...expected.map((v, i) => v - x0[i]));

    const straightness = checkTolerance(
      FLAT_PROPAGATION_STRAIGHTNESS,
      deviation / distance,
      1,
    );
    expect(straightness.withinTolerance, straightness.message).toBe(true);

    // The tangent obeys dk^mu/dlambda = 0 exactly, so it must come back unchanged.
    expect(Array.from(result.final.tangent)).toEqual(Array.from(null_wavevector_k));

    const nullDrift = checkTolerance(
      NORMALIZATION_DRIFT_LONG_RUN,
      normalizationResidual(minkowski, result.final),
    );
    expect(nullDrift.withinTolerance, nullDrift.message).toBe(true);

    for (const drift of conservedQuantityDrift(minkowski, initial, result.final)) {
      if (drift.initial === 0) {
        // A relative measure is undefined; the absolute drift must still vanish.
        expect(drift.absoluteDrift, `${drift.name} drifted from zero`).toBe(0);
        continue;
      }
      const check = checkTolerance(
        CONSERVED_QUANTITY_DRIFT_FLAT,
        drift.absoluteDrift,
        drift.initial,
      );
      expect(check.withinTolerance, `${drift.name}: ${check.message}`).toBe(true);
    }

    // angular_momentum_Lz must be genuinely exercised, not trivially zero.
    const lz = conservedQuantityDrift(minkowski, initial, result.final).find(
      (d) => d.name === 'angular_momentum_Lz',
    );
    expect(lz).toBeDefined();
    expect(Math.abs(lz!.initial)).toBeGreaterThan(0);
  });

  it('holds a timelike worldline straight under the adaptive integrator', () => {
    const x0: Vec4 = [0, 3, -2, 1];
    // At rest in this frame: u^mu = (1,0,0,0), already normalized to -1.
    const four_velocity_u: Vec4 = [1, 0, 0, 0];
    const initial = timelikeState(x0, four_velocity_u);

    const result = integrateGeodesic({
      model: minkowski,
      integrator: new RKF45Integrator(STATE_DIM),
      initial,
      limits: {
        initialStep: STEP_SIZE,
        parameterMax: STEPS * STEP_SIZE,
        maxSteps: STEPS + 10,
        // Cap the step so the run genuinely takes 10^6 of them. Without a cap the
        // adaptive controller would correctly identify this problem as trivial and
        // stride across it in a handful of steps, which would not exercise the gate.
        maxStep: STEP_SIZE,
      },
    });

    expect(result.reason).toBe('parameter-limit');
    expect(result.steps).toBeGreaterThanOrEqual(STEPS);
    expect(result.rejectedSteps).toBe(0);

    const tau = result.final.parameter;
    const expected = exactPosition(x0, four_velocity_u, tau);
    const deviation = maxComponentDeviation(result.final.position_x, expected);

    const straightness = checkTolerance(
      FLAT_PROPAGATION_STRAIGHTNESS,
      deviation / Math.max(tau, 1),
      1,
    );
    expect(straightness.withinTolerance, straightness.message).toBe(true);

    const timelikeDrift = checkTolerance(
      NORMALIZATION_DRIFT_LONG_RUN,
      normalizationResidual(minkowski, result.final),
    );
    expect(timelikeDrift.withinTolerance, timelikeDrift.message).toBe(true);
  });

  it('spends no rejected steps on a problem with zero local error', () => {
    // With a vanishing connection the embedded pair agrees exactly, so the scaled error
    // norm is 0 and the controller should grow the step to its ceiling every time.
    const integrator = new RKF45Integrator(STATE_DIM);
    const result = integrateGeodesic({
      model: minkowski,
      integrator,
      initial: nullState([0, 0, 0, 0], [1, 1, 0, 0]),
      limits: { initialStep: 1e-3, parameterMax: 1e3, maxSteps: 100_000 },
    });
    expect(result.rejectedSteps).toBe(0);
    expect(result.maxErrorNorm).toBe(0);
    // Step growth is clamped at 5x per accepted step, so reaching lambda = 1000 from
    // h = 1e-3 needs far fewer steps than the fixed-step run above.
    expect(result.steps).toBeLessThan(100);
  });
});
