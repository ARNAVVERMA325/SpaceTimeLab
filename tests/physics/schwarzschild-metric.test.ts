import { describe, expect, it } from 'vitest';
import { symmetryResidual } from '../../src/physics/core/christoffel.js';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { contractLower } from '../../src/physics/core/metric-tensor.js';
import {
  christoffelFromMetricNumeric,
  christoffelMaxDifference,
} from '../../src/physics/geometry/christoffel-numeric.js';
import { inverseResidual, invertMetric4 } from '../../src/physics/geometry/invert-metric.js';
import {
  horizonRadius,
  kretschmann,
  photonSphereRadius,
  POLAR_AXIS_SIN_THETA_FLOOR,
  schwarzschild,
} from '../../src/physics/spacetimes/schwarzschild.js';
import {
  METRIC_INVERSE_RESIDUAL,
  SCHWARZSCHILD_CHRISTOFFEL_NUMERIC,
  checkTolerance,
} from '../../src/physics/validation/tolerances.js';

/**
 * ROADMAP.md 2A.1 — the Schwarzschild metric module.
 *
 * The analytical Christoffel symbols are checked against closed-form expressions and,
 * independently, against symbols built from numerically differentiated metric
 * components. Both were derived from the definition in CLAUDE.md §2; the symbolic
 * derivation additionally confirmed R_mu_nu = 0 and K = 48 M^2 / r^6.
 */

const M = 1;
const model = schwarzschild(M);

const EXTERIOR_EVENTS: readonly Vec4[] = [
  [0, 3, Math.PI / 2, 0],
  [0, 6, Math.PI / 2, 1],
  [0, 10, Math.PI / 3, 2],
  [0, 100, Math.PI / 4, 3],
  [0, 2.1, Math.PI / 2, 0],
  [0, 1000, 0.3, 5],
];

describe('Schwarzschild metric components (ROADMAP.md 2A.1)', () => {
  it('is classified as an exact analytical vacuum solution', () => {
    expect(model.classification).toBe('exact-analytical');
    expect(model.parameters.M).toBe(M);
    expect(model.chart.id).toBe('schwarzschild-spherical');
    expect(model.chart.coordinateNames).toEqual(['t', 'r', 'theta', 'phi']);
  });

  it('declares the symmetries it actually has', () => {
    expect(model.symmetries.stationary).toBe(true);
    expect(model.symmetries.axisymmetric).toBe(true);
    expect(model.symmetries.sphericallySymmetric).toBe(true);
  });

  it('has the diagonal components of the line element', () => {
    for (const x of EXTERIOR_EVENTS) {
      const r = x[1];
      const theta = x[2];
      const f = 1 - (2 * M) / r;
      const g = model.metricAt(x);

      expect(g.g(0, 0)).toBeCloseTo(-f, 15);
      expect(g.g(1, 1)).toBeCloseTo(1 / f, 12);
      expect(g.g(2, 2)).toBeCloseTo(r * r, 9);
      expect(g.g(3, 3)).toBeCloseTo(r * r * Math.sin(theta) ** 2, 9);

      // Strictly diagonal: no off-diagonal terms in this chart.
      for (let mu = 0; mu < 4; mu += 1) {
        for (let nu = 0; nu < 4; nu += 1) {
          if (mu !== nu) expect(g.g(mu as 0, nu as 0)).toBe(0);
        }
      }
    }
  });

  it('has an inverse metric that really inverts it', () => {
    for (const x of EXTERIOR_EVENTS) {
      const g = model.metricAt(x);
      const check = checkTolerance(
        METRIC_INVERSE_RESIDUAL,
        inverseResidual(g.g_mu_nu, g.g_inv_mu_nu),
      );
      expect(check.withinTolerance, `${check.message} at r = ${x[1]}`).toBe(true);

      // And agrees with the general Gauss-Jordan routine, not just with itself.
      const numeric = invertMetric4(g.g_mu_nu);
      const cross = checkTolerance(METRIC_INVERSE_RESIDUAL, inverseResidual(g.g_mu_nu, numeric));
      expect(cross.withinTolerance, `${cross.message} at r = ${x[1]}`).toBe(true);
    }
  });

  it('has a metric determinant of -r^4 sin^2(theta), non-zero at the horizon', () => {
    // CLAUDE.md §6.1 in the most direct form available: the determinant does not vanish
    // at r = 2M, so a determinant test cannot detect the horizon. Asserted here so the
    // rule is enforced by the suite rather than only stated in prose.
    for (const x of EXTERIOR_EVENTS) {
      const g = model.metricAt(x);
      const determinant = g.g(0, 0) * g.g(1, 1) * g.g(2, 2) * g.g(3, 3);
      const expected = -Math.pow(x[1], 4) * Math.sin(x[2]) ** 2;
      expect(Math.abs(determinant - expected) / Math.abs(expected)).toBeLessThan(1e-12);
    }

    // Just outside the horizon the determinant is close to -(2M)^4 and nowhere near 0,
    // while f itself has almost vanished. The two carry completely different information.
    const nearHorizon: Vec4 = [0, 2 * M * (1 + 1e-9), Math.PI / 2, 0];
    const g = model.metricAt(nearHorizon);
    const determinant = g.g(0, 0) * g.g(1, 1) * g.g(2, 2) * g.g(3, 3);
    expect(Math.abs(determinant)).toBeCloseTo(16 * Math.pow(M, 4), 6);
    expect(model.lapseFunction(nearHorizon[1])).toBeLessThan(1e-8);
  });

  it('has a Kretschmann scalar that is finite at the horizon', () => {
    // The horizon is a coordinate singularity, not a curvature singularity.
    expect(Number.isFinite(kretschmann(M, horizonRadius(M)))).toBe(true);
    expect(kretschmann(M, horizonRadius(M))).toBeCloseTo(0.75, 12);
  });

  it('matches the closed-form Christoffel symbols exactly', () => {
    for (const x of EXTERIOR_EVENTS) {
      const r = x[1];
      const theta = x[2];
      const f = 1 - (2 * M) / r;
      const c = model.christoffelAt(x);

      expect(c.get(0, 0, 1)).toBeCloseTo(M / (r * r * f), 12);
      expect(c.get(1, 0, 0)).toBeCloseTo((M * f) / (r * r), 12);
      expect(c.get(1, 1, 1)).toBeCloseTo(-M / (r * r * f), 12);
      expect(c.get(1, 2, 2)).toBeCloseTo(-r * f, 10);
      expect(c.get(1, 3, 3)).toBeCloseTo(-r * f * Math.sin(theta) ** 2, 10);
      expect(c.get(2, 1, 2)).toBeCloseTo(1 / r, 15);
      expect(c.get(2, 3, 3)).toBeCloseTo(-Math.sin(theta) * Math.cos(theta), 15);
      expect(c.get(3, 1, 3)).toBeCloseTo(1 / r, 15);
      expect(c.get(3, 2, 3)).toBeCloseTo(Math.cos(theta) / Math.sin(theta), 13);
    }
  });

  it('has Christoffel symbols symmetric in the lower index pair', () => {
    for (const x of EXTERIOR_EVENTS) {
      expect(symmetryResidual(model.christoffelAt(x))).toBe(0);
    }
  });

  it('agrees with Christoffel symbols built from numerical metric derivatives', () => {
    for (const x of EXTERIOR_EVENTS) {
      const analytic = model.christoffelAt(x);
      const numeric = christoffelFromMetricNumeric(model, x);
      const scale = analytic.maxAbs();
      const check = checkTolerance(
        SCHWARZSCHILD_CHRISTOFFEL_NUMERIC,
        christoffelMaxDifference(numeric, analytic),
        scale,
      );
      expect(check.withinTolerance, `${check.message} at r = ${x[1]}`).toBe(true);
    }
  });

  it('writes the same symbols through the allocation-free hot path', () => {
    const buffer = new Float64Array(64).fill(7);
    const x: Vec4 = [0, 8, Math.PI / 3, 1];
    model.christoffelInto(x, buffer);
    const reference = model.christoffelAt(x);
    for (let i = 0; i < 64; i += 1) expect(buffer[i]).toBe(reference.components[i]);
  });

  it('reduces to flat space as M tends to zero (CLAUDE.md §16)', () => {
    // With M = 0 the metric is Minkowski written in spherical coordinates. The
    // Christoffel symbols do not vanish — spherical coordinates are curvilinear — but
    // the mass-dependent ones do, and the curvature is identically zero.
    const x: Vec4 = [0, 10, Math.PI / 3, 1];
    const r = x[1];
    const theta = x[2];

    let previous = Infinity;
    for (const tiny of [1e-3, 1e-5, 1e-7]) {
      const weak = schwarzschild(tiny);
      const c = weak.christoffelAt(x);
      // The mass-sourced symbol Gamma^r_{t t} = M f / r^2 vanishes linearly in M.
      const massSourced = Math.abs(c.get(1, 0, 0));
      expect(massSourced).toBeLessThan(previous);
      previous = massSourced;

      // The purely geometric ones survive, because they come from the coordinates
      // rather than from the mass. In flat space in spherical coordinates
      // Gamma^theta_{r theta} = Gamma^phi_{r phi} = 1/r and Gamma^phi_{theta phi} = cot(theta).
      expect(c.get(2, 1, 2)).toBeCloseTo(1 / r, 15);
      expect(c.get(3, 1, 3)).toBeCloseTo(1 / r, 15);
      expect(c.get(3, 2, 3)).toBeCloseTo(Math.cos(theta) / Math.sin(theta), 13);

      // Gamma^r_{theta theta} = -r f = -(r - 2M) approaches its flat value -r linearly
      // in M, with a deviation of exactly 2M. Asserting the rate is a stronger statement
      // than asserting closeness at any one mass.
      //
      // The bound is the cancellation floor, not a round number: recovering 2M by adding
      // r to -(r - 2M) subtracts two quantities of order r, so the absolute error is
      // about r * eps. At r = 10 that is 2.2e-15, and demanding better would be
      // demanding precision the arithmetic cannot deliver.
      const cancellationFloor = 4 * Number.EPSILON * r;
      expect(Math.abs(c.get(1, 2, 2) + r - 2 * tiny)).toBeLessThan(cancellationFloor);

      // And the curvature invariant goes to zero with M^2.
      expect(kretschmann(tiny, r)).toBeLessThan(48 * tiny * tiny);
    }

    // A null vector stays null under the M -> 0 metric, with the flat spherical form.
    const weak = schwarzschild(1e-12);
    const g = weak.metricAt(x);
    const radialNull: Vec4 = [1, 1, 0, 0];
    expect(Math.abs(contractLower(g, radialNull, radialNull))).toBeLessThan(1e-11);
  });
});

describe('Schwarzschild chart domain (ROADMAP.md 2A.2, CLAUDE.md §6.1)', () => {
  it('accepts the exterior', () => {
    for (const x of EXTERIOR_EVENTS) {
      expect(model.domainCheck(x).inDomain, `r = ${x[1]}`).toBe(true);
    }
  });

  it('rejects the horizon and the interior, naming f rather than the determinant', () => {
    for (const r of [2 * M, 2 * M * (1 - 1e-12), M, 0.1 * M]) {
      const status = model.domainCheck([0, r, Math.PI / 2, 0]);
      expect(status.inDomain, `r = ${r} should be outside the exterior chart`).toBe(false);
      if (!status.inDomain) {
        expect(status.code).toBe('coordinate-breakdown');
        expect(status.reason).toContain('f = 1 - 2M/r');
        // The message must not claim a curvature singularity at the horizon.
        expect(status.reason).toContain('not a curvature singularity');
      }
    }
  });

  it('rejects the polar axis as a coordinate artifact, not a physical boundary', () => {
    const status = model.domainCheck([0, 10, POLAR_AXIS_SIN_THETA_FLOOR / 2, 0]);
    expect(status.inDomain).toBe(false);
    if (!status.inDomain) {
      expect(status.reason).toContain('artifact of spherical coordinates');
      expect(status.reason).toContain('not a physical boundary');
    }
    // Well away from the axis the chart is fine.
    expect(model.domainCheck([0, 10, 0.1, 0]).inDomain).toBe(true);
  });

  it('rejects non-finite coordinates', () => {
    expect(model.domainCheck([0, Number.NaN, Math.PI / 2, 0]).inDomain).toBe(false);
    expect(model.domainCheck([0, Number.POSITIVE_INFINITY, Math.PI / 2, 0]).inDomain).toBe(false);
  });

  it('refuses to evaluate the metric inside the horizon rather than returning nonsense', () => {
    expect(() => model.metricAt([0, M, Math.PI / 2, 0])).toThrow(RangeError);
    expect(() => model.metricAt([0, 2 * M, Math.PI / 2, 0])).toThrow(RangeError);
  });

  it('refuses a non-positive or non-finite mass', () => {
    expect(() => schwarzschild(0)).toThrow(RangeError);
    expect(() => schwarzschild(-1)).toThrow(RangeError);
    expect(() => schwarzschild(Number.NaN)).toThrow(RangeError);
  });

  it('places the horizon and photon sphere where they belong', () => {
    expect(horizonRadius(M)).toBe(2 * M);
    expect(photonSphereRadius(M)).toBe(3 * M);
    expect(horizonRadius(7)).toBe(14);
  });
});
