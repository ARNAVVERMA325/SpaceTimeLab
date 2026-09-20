import { describe, expect, it } from 'vitest';
import { contravariant, covariant, requireVariance } from '../../src/physics/core/four-vector.js';
import { MetricTensor } from '../../src/physics/core/metric-tensor.js';
import { inverseResidual, invertMetric4 } from '../../src/physics/geometry/invert-metric.js';
import {
  four_momentum_p,
  lowerIndex,
  raiseIndex,
} from '../../src/physics/geometry/raise-lower.js';
import { minkowski } from '../../src/physics/spacetimes/minkowski.js';
import { METRIC_INVERSE_RESIDUAL, checkTolerance } from '../../src/physics/validation/tolerances.js';

/**
 * Differential-geometry layer (CLAUDE.md §20).
 *
 * The inversion routine is exercised on off-diagonal metrics here, well before Kerr
 * needs it in Milestone 4A, so that a defect surfaces in an isolated unit test rather
 * than inside a spinning-black-hole render.
 */

describe('4x4 metric inversion', () => {
  it('inverts the Minkowski metric to itself', () => {
    const g = Float64Array.from([
      -1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
    ]);
    // Elimination leaves some off-diagonal zeros as -0. That is numerically identical
    // to +0 and harmless in every contraction, but strict equality distinguishes them,
    // so the sign of zero is normalized before comparing.
    const normalized = Array.from(invertMetric4(g), (v) => v + 0);
    expect(normalized).toEqual(Array.from(g));
  });

  it('inverts a diagonal metric with unequal entries', () => {
    // The shape of a static spherically symmetric metric in spherical coordinates.
    const g = Float64Array.from([
      -0.5, 0, 0, 0,
      0, 2, 0, 0,
      0, 0, 16, 0,
      0, 0, 0, 9,
    ]);
    const gInv = invertMetric4(g);
    const check = checkTolerance(METRIC_INVERSE_RESIDUAL, inverseResidual(g, gInv));
    expect(check.withinTolerance, check.message).toBe(true);
    expect(gInv[0]).toBeCloseTo(-2, 15);
    expect(gInv[5]).toBeCloseTo(0.5, 15);
  });

  it('inverts a metric with off-diagonal t-phi mixing', () => {
    // The structural form of Kerr in Boyer-Lindquist: a g_{t phi} cross term. Numbers
    // are illustrative, not a physical Kerr metric at any particular event.
    const g = Float64Array.from([
      -0.7, 0, 0, -1.3,
      0, 1.9, 0, 0,
      0, 0, 4.0, 0,
      -1.3, 0, 0, 12.5,
    ]);
    const gInv = invertMetric4(g);
    const check = checkTolerance(METRIC_INVERSE_RESIDUAL, inverseResidual(g, gInv));
    expect(check.withinTolerance, check.message).toBe(true);
  });

  it('inverts correctly when a leading diagonal entry is tiny, via partial pivoting', () => {
    const g = Float64Array.from([
      1e-13, 0, 0, 1,
      0, 1, 0, 0,
      0, 0, 1, 0,
      1, 0, 0, 1,
    ]);
    const gInv = invertMetric4(g);
    const check = checkTolerance(METRIC_INVERSE_RESIDUAL, inverseResidual(g, gInv));
    expect(check.withinTolerance, check.message).toBe(true);
  });

  it('refuses a singular metric rather than returning garbage', () => {
    const singular = Float64Array.from([
      0, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
    ]);
    expect(() => invertMetric4(singular)).toThrow(RangeError);
  });

  it('refuses a non-finite metric rather than propagating NaN', () => {
    const bad = Float64Array.from([
      Number.NaN, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
    ]);
    expect(() => invertMetric4(bad)).toThrow(RangeError);
  });

  it('rejects a wrongly sized input', () => {
    expect(() => invertMetric4(new Float64Array(9))).toThrow(RangeError);
  });

  it('refuses to build a degenerate diagonal metric', () => {
    expect(() => MetricTensor.diagonal(-1, 1, 0, 1)).toThrow(RangeError);
  });
});

describe('raising and lowering indices', () => {
  it('round-trips a vector through lower then raise', () => {
    const metric = minkowski.metricAt([0, 0, 0, 0]);
    const up = contravariant([2, -3, 5, 7]);
    const down = lowerIndex(metric, up);
    expect(down.variance).toBe('covariant');
    // In Minkowski the time component flips sign and the spatial ones do not.
    expect(Array.from(down.components)).toEqual([-2, -3, 5, 7]);

    const backUp = raiseIndex(metric, down);
    expect(backUp.variance).toBe('contravariant');
    expect(Array.from(backUp.components)).toEqual([2, -3, 5, 7]);
  });

  it('refuses to lower an index that is already down (CLAUDE.md §24)', () => {
    const metric = minkowski.metricAt([0, 0, 0, 0]);
    expect(() => lowerIndex(metric, covariant([1, 0, 0, 0]))).toThrow(TypeError);
    expect(() => raiseIndex(metric, contravariant([1, 0, 0, 0]))).toThrow(TypeError);
  });

  it('names the mismatch clearly when a variance check fails', () => {
    expect(() => requireVariance(covariant([0, 0, 0, 0]), 'contravariant', 'someRoutine')).toThrow(
      /someRoutine: expected a contravariant four-vector/,
    );
  });

  it('builds p_mu from a tangent, giving E = -p_t for a photon', () => {
    const metric = minkowski.metricAt([0, 0, 0, 0]);
    const null_wavevector_k = [3, 3, 0, 0] as const;
    const p = four_momentum_p(metric, null_wavevector_k);
    expect(p[0]).toBe(-3);
    // E = -p_t, matching CLAUDE.md §16 and the sign convention the model declares.
    expect(-p[0]).toBe(3);
  });
});
