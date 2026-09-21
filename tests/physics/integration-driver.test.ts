import { describe, expect, it } from 'vitest';
import { CONVENTIONS } from '../../src/physics/conventions.js';
import { ChristoffelSymbols } from '../../src/physics/core/christoffel.js';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { MetricTensor } from '../../src/physics/core/metric-tensor.js';
import { nullState } from '../../src/physics/core/phase-space.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { RK4Integrator } from '../../src/physics/geodesic/integrators/rk4.js';
import { RKF45Integrator } from '../../src/physics/geodesic/integrators/rkf45.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { minkowski } from '../../src/physics/spacetimes/minkowski.js';
import { MINKOWSKI_CARTESIAN_CHART } from '../../src/physics/spacetimes/minkowski.js';
import type { DomainStatus, SpacetimeModel } from '../../src/physics/spacetimes/spacetime-model.js';
import { IN_DOMAIN } from '../../src/physics/spacetimes/spacetime-model.js';
import { summarizeHealth } from '../../src/physics/validation/numeric-health.js';

/**
 * Geodesic driver behaviour (CLAUDE.md §17).
 *
 * CLAUDE.md §17 requires NaN and infinity detection to be explicit, requires a failure
 * to expose the integrator, model, coordinate system, parameter, step size and failed
 * quantity, and forbids replacing a failed calculation with a plausible substitute.
 * These tests use deliberately pathological models to drive each of those paths.
 */

/** Base for test-only models: Minkowski geometry with one behaviour overridden. */
function testModel(overrides: Partial<SpacetimeModel> & { id: string }): SpacetimeModel {
  return {
    displayName: overrides.id,
    classification: 'exact-analytical',
    chart: MINKOWSKI_CARTESIAN_CHART,
    conventions: CONVENTIONS,
    parameters: {},
    killingVectors: [],
    symmetries: minkowski.symmetries,
    geometry: minkowski.geometry,
    description: 'Test-only model. Not a physical spacetime.',
    metricAt: (x: Vec4) => minkowski.metricAt(x),
    christoffelAt: (x: Vec4) => minkowski.christoffelAt(x),
    christoffelInto: (x: Vec4, out: Float64Array) => minkowski.christoffelInto(x, out),
    domainCheck: (x: Vec4) => minkowski.domainCheck(x),
    ...overrides,
  } as SpacetimeModel;
}

describe('integration driver termination (CLAUDE.md §17)', () => {
  it('stops at the parameter limit and reports it', () => {
    const result = integrateGeodesic({
      model: minkowski,
      integrator: new RK4Integrator(STATE_DIM),
      initial: nullState([0, 0, 0, 0], [1, 1, 0, 0]),
      limits: { initialStep: 0.1, parameterMax: 10, maxSteps: 10_000 },
    });
    expect(result.reason).toBe('parameter-limit');
    expect(result.final.parameter).toBeCloseTo(10, 12);
    expect(result.diagnostics).toBeUndefined();
  });

  it('stops at the step budget and reports it', () => {
    const result = integrateGeodesic({
      model: minkowski,
      integrator: new RK4Integrator(STATE_DIM),
      initial: nullState([0, 0, 0, 0], [1, 1, 0, 0]),
      limits: { initialStep: 0.1, parameterMax: 1e9, maxSteps: 25 },
    });
    expect(result.reason).toBe('step-limit');
    expect(result.steps).toBe(25);
  });

  it('fires a caller-supplied terminator', () => {
    const result = integrateGeodesic({
      model: minkowski,
      integrator: new RK4Integrator(STATE_DIM),
      initial: nullState([0, 0, 0, 0], [1, 1, 0, 0]),
      limits: { initialStep: 0.1, parameterMax: 1e6, maxSteps: 10_000 },
      terminator: (position_x) => position_x[1] >= 5,
    });
    expect(result.reason).toBe('terminator');
    expect(result.final.position_x[1]).toBeGreaterThanOrEqual(5);
  });

  it('reports a domain exit with the model reason, not a silent stop', () => {
    const walled = testModel({
      id: 'walled',
      domainCheck: (x: Vec4): DomainStatus =>
        x[1] > 3
          ? {
              inDomain: false,
              code: 'coordinate-breakdown',
              reason: 'Test model: the chart is declared invalid beyond x = 3.',
            }
          : IN_DOMAIN,
    });

    const result = integrateGeodesic({
      model: walled,
      integrator: new RK4Integrator(STATE_DIM),
      initial: nullState([0, 0, 0, 0], [1, 1, 0, 0]),
      limits: { initialStep: 0.1, parameterMax: 1e6, maxSteps: 10_000 },
    });

    expect(result.reason).toBe('domain-exit');
    expect(result.diagnostics?.failedQuantity).toBe('position_x');
    expect(result.diagnostics?.detail).toContain('beyond x = 3');
    expect(result.diagnostics?.coordinateSystem).toBe('minkowski-cartesian');
  });

  it('refuses to start outside the chart domain', () => {
    const closed = testModel({
      id: 'closed',
      domainCheck: (): DomainStatus => ({
        inDomain: false,
        code: 'outside-chart',
        reason: 'Test model: no event is in the domain.',
      }),
    });

    const result = integrateGeodesic({
      model: closed,
      integrator: new RK4Integrator(STATE_DIM),
      initial: nullState([0, 0, 0, 0], [1, 1, 0, 0]),
    });

    expect(result.reason).toBe('domain-exit');
    expect(result.steps).toBe(0);
  });

  it('detects a NaN, returns the last finite state, and reports full diagnostics', () => {
    // Christoffel symbols that turn non-finite past x = 1, so the tangent update
    // produces NaN rather than the state simply growing large.
    const poisoned = testModel({
      id: 'poisoned',
      christoffelInto: (x: Vec4, out: Float64Array) => {
        out.fill(0);
        if (x[1] > 1) out[1 * 16 + 1 * 4 + 1] = Number.NaN;
      },
      christoffelAt: (x: Vec4) => {
        const c = new Float64Array(64);
        if (x[1] > 1) c[1 * 16 + 1 * 4 + 1] = Number.NaN;
        return new ChristoffelSymbols(c);
      },
    });

    const result = integrateGeodesic({
      model: poisoned,
      integrator: new RK4Integrator(STATE_DIM),
      initial: nullState([0, 0, 0, 0], [1, 1, 0, 0]),
      limits: { initialStep: 0.1, parameterMax: 100, maxSteps: 10_000 },
    });

    expect(result.reason).toBe('numerical-failure');

    // The returned state is the last finite one; nothing was invented to replace it.
    for (const v of result.final.position_x) expect(Number.isFinite(v)).toBe(true);
    for (const v of result.final.tangent) expect(Number.isFinite(v)).toBe(true);

    // CLAUDE.md §17's required diagnostic fields.
    const d = result.diagnostics;
    expect(d).toBeDefined();
    expect(d!.integrator).toBe('rk4');
    expect(d!.model).toBe('poisoned');
    expect(d!.coordinateSystem).toBe('minkowski-cartesian');
    expect(d!.parameterName).toBe('lambda');
    expect(Number.isFinite(d!.parameterValue)).toBe(true);
    expect(d!.stepSize).toBeGreaterThan(0);
    expect(d!.failedQuantity).toMatch(/^(x|tangent)\^[0-3]$/);
    expect(d!.detail).toContain('no substitute value was invented');
  });

  it('reports step underflow when the requested tolerance cannot be met', () => {
    // A right-hand side that is violently stiff in the tangent sector, with a tolerance
    // no finite step can satisfy, drives the controller into the step floor.
    const stiff = testModel({
      id: 'stiff',
      christoffelInto: (x: Vec4, out: Float64Array) => {
        out.fill(0);
        out[1 * 16 + 1 * 4 + 1] = 1e12 * Math.sin(1e9 * x[1]);
      },
    });

    const result = integrateGeodesic({
      model: stiff,
      integrator: new RKF45Integrator(STATE_DIM, {
        tolerance: { absolute: 1e-18, relative: 1e-18 },
      }),
      initial: nullState([0, 0, 0, 0], [1, 1, 0, 0]),
      limits: { initialStep: 1, parameterMax: 100, maxSteps: 10_000, minStep: 1e-12 },
    });

    expect(result.reason).toBe('step-underflow');
    expect(result.diagnostics?.failedQuantity).toBe('step size');
    expect(result.diagnostics?.detail).toContain('cannot be met');
    expect(result.rejectedSteps).toBeGreaterThan(0);
  });

  it('records the path only when asked', () => {
    const withoutPath = integrateGeodesic({
      model: minkowski,
      integrator: new RK4Integrator(STATE_DIM),
      initial: nullState([0, 0, 0, 0], [1, 1, 0, 0]),
      limits: { initialStep: 0.5, parameterMax: 5, maxSteps: 100 },
    });
    expect(withoutPath.path).toBeUndefined();

    const withPath = integrateGeodesic({
      model: minkowski,
      integrator: new RK4Integrator(STATE_DIM),
      initial: nullState([0, 0, 0, 0], [1, 1, 0, 0]),
      limits: { initialStep: 0.5, parameterMax: 5, maxSteps: 100 },
      recordPath: true,
    });
    expect(withPath.path?.length).toBe(withPath.steps + 1);
  });

  it('carries the metric determinant sign without using it as a domain test', () => {
    // CLAUDE.md §6.1: the horizon is never detected from a vanishing determinant. The
    // Minkowski determinant is -1 everywhere and the domain check ignores it entirely.
    const metric: MetricTensor = minkowski.metricAt([0, 0, 0, 0]);
    expect(metric.g(0, 0)).toBe(-1);
    expect(minkowski.domainCheck([0, 1e9, 0, 0]).inDomain).toBe(true);
  });
});

describe('numerical health summary (CLAUDE.md §17)', () => {
  it('reports ok when every check passes', () => {
    expect(summarizeHealth([{ withinTolerance: true, message: 'fine' }], false).level).toBe('ok');
  });

  it('reports degraded, with the message, when a tolerance is exceeded', () => {
    const health = summarizeHealth(
      [{ withinTolerance: false, message: 'Numerical error exceeds validated tolerance: x' }],
      false,
    );
    expect(health.level).toBe('degraded');
    expect(health.messages[0]).toContain('exceeds validated tolerance');
  });

  it('reports failed when a non-finite value appeared', () => {
    const health = summarizeHealth([{ withinTolerance: true, message: 'fine' }], true);
    expect(health.level).toBe('failed');
    expect(health.messages[0]).toContain('NaN or infinity');
  });
});
