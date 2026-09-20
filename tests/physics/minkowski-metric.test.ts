import { describe, expect, it } from 'vitest';
import { CONVENTIONS, SIGNATURE } from '../../src/physics/conventions.js';
import { symmetryResidual } from '../../src/physics/core/christoffel.js';
import { contractLower } from '../../src/physics/core/metric-tensor.js';
import { inverseResidual, invertMetric4 } from '../../src/physics/geometry/invert-metric.js';
import {
  christoffelFromMetricNumeric,
  christoffelMaxDifference,
} from '../../src/physics/geometry/christoffel-numeric.js';
import { minkowski } from '../../src/physics/spacetimes/minkowski.js';
import {
  CHRISTOFFEL_NUMERIC_VS_ANALYTIC,
  METRIC_INVERSE_RESIDUAL,
  checkTolerance,
} from '../../src/physics/validation/tolerances.js';
import type { Vec4 } from '../../src/physics/core/indices.js';

/**
 * ROADMAP.md 1.2 — Minkowski metric and its Christoffel symbols.
 *
 * The roadmap asks for the Christoffel check to be analytical, not numerical: the
 * symbols are identically zero because the metric components are constant, so the
 * assertion is exact equality, not a tolerance.
 */

const SAMPLE_EVENTS: readonly Vec4[] = [
  [0, 0, 0, 0],
  [1, 2, -3, 4],
  [-17.5, 100, 0.001, -250],
  [1e6, -1e6, 1e3, 1e-3],
];

describe('Minkowski metric (ROADMAP.md 1.2)', () => {
  it('declares the (-,+,+,+) signature the rest of the engine assumes', () => {
    expect(SIGNATURE).toBe('(-,+,+,+)');
    expect(minkowski.conventions).toBe(CONVENTIONS);
  });

  it('is classified as an exact analytical solution (CLAUDE.md §10, §11)', () => {
    expect(minkowski.classification).toBe('exact-analytical');
  });

  it('has g_mu_nu = diag(-1, 1, 1, 1) at every event', () => {
    for (const x of SAMPLE_EVENTS) {
      const metric = minkowski.metricAt(x);
      expect(Array.from(metric.g_mu_nu)).toEqual([
        -1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
      ]);
    }
  });

  it('has an inverse metric that is its own analytical inverse', () => {
    const metric = minkowski.metricAt([0, 0, 0, 0]);
    expect(Array.from(metric.g_inv_mu_nu)).toEqual(Array.from(metric.g_mu_nu));
    expect(inverseResidual(metric.g_mu_nu, metric.g_inv_mu_nu)).toBe(0);
  });

  it('agrees with the general 4x4 inversion routine', () => {
    const metric = minkowski.metricAt([0, 0, 0, 0]);
    const numericInverse = invertMetric4(metric.g_mu_nu);
    const check = checkTolerance(
      METRIC_INVERSE_RESIDUAL,
      inverseResidual(metric.g_mu_nu, numericInverse),
    );
    expect(check.withinTolerance, check.message).toBe(true);
  });

  it('hands out defensive copies so a caller cannot corrupt the shared metric', () => {
    const first = minkowski.metricAt([0, 0, 0, 0]);
    first.g_mu_nu[0] = 42;
    const second = minkowski.metricAt([0, 0, 0, 0]);
    expect(second.g_mu_nu[0]).toBe(-1);
  });

  it('has Christoffel symbols that are identically zero, exactly', () => {
    for (const x of SAMPLE_EVENTS) {
      const christoffel = minkowski.christoffelAt(x);
      // Exact equality, not a tolerance: ROADMAP.md 1.2 calls for an analytical check.
      expect(christoffel.maxAbs()).toBe(0);
      expect(symmetryResidual(christoffel)).toBe(0);
    }
  });

  it('writes the same zero symbols through the allocation-free hot path', () => {
    const buffer = new Float64Array(64).fill(7);
    minkowski.christoffelInto([3, 1, 4, 1], buffer);
    expect(buffer.every((v) => v === 0)).toBe(true);
  });

  it('rejects a wrongly sized Christoffel buffer rather than writing past it', () => {
    expect(() => minkowski.christoffelInto([0, 0, 0, 0], new Float64Array(16))).toThrow(RangeError);
  });

  it('recovers zero Christoffel symbols through numerical metric differentiation', () => {
    // A cross-check of the differential-geometry layer against the analytical answer.
    // This one is a tolerance check: central differencing has a roundoff floor.
    for (const x of SAMPLE_EVENTS) {
      const numeric = christoffelFromMetricNumeric(minkowski, x);
      const analytic = minkowski.christoffelAt(x);
      const check = checkTolerance(
        CHRISTOFFEL_NUMERIC_VS_ANALYTIC,
        christoffelMaxDifference(numeric, analytic),
      );
      expect(check.withinTolerance, `${check.message} at event ${JSON.stringify(x)}`).toBe(true);
    }
  });

  it('measures the flat-space interval correctly for a known separation', () => {
    const metric = minkowski.metricAt([0, 0, 0, 0]);
    // A displacement of 3 in t and 5 in x has ds^2 = -9 + 25 = 16 (spacelike).
    const separation: Vec4 = [3, 5, 0, 0];
    expect(contractLower(metric, separation, separation)).toBe(16);
    // A null separation: equal time and space intervals give ds^2 = 0.
    const nullSeparation: Vec4 = [5, 5, 0, 0];
    expect(contractLower(metric, nullSeparation, nullSeparation)).toBe(0);
  });

  it('reports every event as inside the chart domain (no horizon exists)', () => {
    for (const x of SAMPLE_EVENTS) {
      expect(minkowski.domainCheck(x).inDomain).toBe(true);
    }
    expect(minkowski.chart.horizonPenetrating).toBe(false);
  });
});
