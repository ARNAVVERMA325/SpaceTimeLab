import { describe, expect, it } from 'vitest';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { RKF45Integrator } from '../../src/physics/geodesic/integrators/rkf45.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { generateNullRay, staticObserver } from '../../src/physics/observer/observer.js';
import { minkowski } from '../../src/physics/spacetimes/minkowski.js';
import { schwarzschild } from '../../src/physics/spacetimes/schwarzschild.js';
import { isCaptured } from '../../src/physics/spacetimes/schwarzschild-rays.js';
import { normalizationResidual } from '../../src/physics/validation/normalization.js';
import {
  liftDirection,
  liftPosition,
  reduceToOrbitalPlane,
} from '../../src/visualization/orbital-plane.js';
import type { CartesianVec3 } from '../../src/physics/spacetimes/spacetime-model.js';
import { overrideModel } from '../helpers/model-override.js';

/**
 * The orbital-plane reduction, checked against full three-dimensional integration.
 *
 * The reduction is exact under spherical symmetry, but "exact in principle" is not a
 * test. These cases run the same ray both ways — reduced into its own equatorial plane,
 * and integrated in the original chart with all four coordinates live — and require the
 * results to agree. Without that, a subtly wrong rotation would show up only as a
 * plausible-looking image, which is exactly what CLAUDE.md §16 says does not count.
 */

const M = 1;
const model = schwarzschild(M);
const observer = staticObserver(model, [0, 20, Math.PI / 2, 0]);

function unit(v: CartesianVec3): CartesianVec3 {
  const n = Math.hypot(...v);
  return [v[0] / n, v[1] / n, v[2] / n];
}

function angleBetween(a: CartesianVec3, b: CartesianVec3): number {
  const ua = unit(a);
  const ub = unit(b);
  // Chord rather than acos: near-parallel vectors make acos ill-conditioned, and these
  // are meant to be parallel.
  return Math.hypot(ua[0] - ub[0], ua[1] - ub[1], ua[2] - ub[2]);
}

describe('orbital-plane reduction', () => {
  it('round-trips position and direction through the rotated frame', () => {
    for (const local of [
      [-1, 0.1, 0.2],
      [-0.8, -0.3, 0.5],
      [-1, 0, 0.05],
      [-0.6, 0.6, -0.4],
    ] as const) {
      const ray = generateNullRay(observer, unit(local));
      const reduced = reduceToOrbitalPlane(model, ray);

      const originalPosition = model.geometry.toCartesianPosition(ray.position_x);
      const originalDirection = model.geometry.toCartesianDirection(ray.position_x, ray.tangent);

      const liftedPosition = liftPosition(model, reduced.frame, reduced.initial);
      const liftedDirection = liftDirection(model, reduced.frame, reduced.initial);

      for (let i = 0; i < 3; i += 1) {
        expect(Math.abs(liftedPosition[i] - originalPosition[i])).toBeLessThan(1e-11);
        expect(Math.abs(liftedDirection[i] - originalDirection[i])).toBeLessThan(1e-11);
      }
      expect(angleBetween(liftedDirection, originalDirection)).toBeLessThan(1e-13);
    }
  });

  it('places the reduced ray exactly in the equatorial plane', () => {
    const ray = generateNullRay(observer, unit([-0.7, 0.4, 0.6]));
    const reduced = reduceToOrbitalPlane(model, ray);
    expect(reduced.initial.position_x[2]).toBe(Math.PI / 2);
    expect(reduced.initial.position_x[3]).toBe(0);
    expect(reduced.initial.tangent[2]).toBe(0);
    // Radius and the time and radial tangent components are untouched by the rotation.
    expect(reduced.initial.position_x[1]).toBeCloseTo(ray.position_x[1], 12);
    expect(reduced.initial.tangent[0]).toBe(ray.tangent[0]);
    expect(reduced.initial.tangent[1]).toBe(ray.tangent[1]);
  });

  it('preserves the null condition', () => {
    for (const local of [
      [-1, 0.2, 0.3],
      [-0.5, 0.8, 0.2],
      [-1, 0, 0],
    ] as const) {
      const ray = generateNullRay(observer, unit(local));
      const reduced = reduceToOrbitalPlane(model, ray);
      expect(Math.abs(normalizationResidual(model, ray))).toBeLessThan(1e-14);
      expect(Math.abs(normalizationResidual(model, reduced.initial))).toBeLessThan(1e-14);
    }
  });

  it('agrees with full three-dimensional integration of the same ray', () => {
    // Two constraints on these directions. They must clear the critical impact
    // parameter, or the rays fall in and their "exit direction" is just wherever the
    // trace happened to stop while spiralling — not a quantity worth comparing. And they
    // must keep the full-3D trace away from the polar axis, since that is where the
    // unreduced chart fails, which is the whole reason the reduction exists.
    //
    // From a static observer at r = 20M, b ~ r sin(psi) / sqrt(f), so escaping needs
    // sin(psi) > b_c sqrt(f) / r ~ 0.247: the transverse part must be a good fraction of
    // the whole.
    const directions = [
      [-1, 0.5, 0.35],
      [-0.9, -0.6, 0.3],
      [-1, 0.4, -0.45],
      [-0.8, 0.7, 0.1],
    ] as const;

    for (const local of directions) {
      const ray = generateNullRay(observer, unit(local));
      const limits = {
        initialStep: 1e-3,
        parameterMax: 100_000,
        maxSteps: 500_000,
        maxStep: 5,
      };
      const tolerance = { absolute: 1e-12, relative: 1e-12 };

      // The comparison radius is far out for a specific reason. The two traces take
      // different step sequences, so their final steps overshoot the target radius by
      // different amounts — several M apart. A ray is still being bent at any finite
      // radius, at a rate falling off as 4 M b / r^3, so two traces ending at different
      // radii legitimately point in slightly different directions. That residual is
      // physics, not disagreement, and at r = 500 it dominated: the measured difference
      // tracked the radius mismatch at 1.7e-7 per unit r, independent of the integration
      // tolerance. At r = 5000 the same mismatch contributes under 1e-9.
      const backgroundRadius = 5000;
      const stopAtBackground = (position_x: readonly [number, number, number, number]): boolean =>
        position_x[1] >= backgroundRadius;

      // Full three-dimensional integration, in the original chart.
      const full = integrateGeodesic({
        model,
        integrator: new RKF45Integrator(STATE_DIM, { tolerance }),
        initial: ray,
        limits,
        terminator: (position_x, tangent) =>
          isCaptured(model, position_x, tangent) || stopAtBackground(position_x),
      });

      // The same ray, reduced into its own orbital plane.
      const reduced = reduceToOrbitalPlane(model, ray);
      const planar = integrateGeodesic({
        model,
        integrator: new RKF45Integrator(STATE_DIM, { tolerance }),
        initial: reduced.initial,
        limits,
        terminator: (position_x, tangent) =>
          isCaptured(model, position_x, tangent) || stopAtBackground(position_x),
      });

      expect(full.reason, `full 3D trace for ${local} ended as ${full.reason}`).toBe('terminator');
      expect(planar.reason).toBe('terminator');

      // Both must have escaped to the background. A captured ray stops wherever the
      // spiral happened to be when the capture test fired, so comparing exit directions
      // of infalling rays would compare arbitrary numbers.
      expect(
        full.final.position_x[1],
        `full 3D ray for ${local} was captured`,
      ).toBeGreaterThanOrEqual(backgroundRadius);
      expect(
        planar.final.position_x[1],
        `planar ray for ${local} was captured`,
      ).toBeGreaterThanOrEqual(backgroundRadius);

      const fullDirection = model.geometry.toCartesianDirection(
        full.final.position_x,
        full.final.tangent,
      );
      const planarDirection = liftDirection(model, reduced.frame, planar.final);

      // The exit directions are what the background lookup uses, so this is the quantity
      // that actually has to agree.
      expect(
        angleBetween(fullDirection, planarDirection),
        `exit directions diverged for local direction ${local}`,
      ).toBeLessThan(1e-8);

      // The rays also end in the same direction from the centre. Compared as angular
      // positions rather than as raw points: the traces stop at radii several M apart,
      // so their positions differ by that overshoot no matter how exact the reduction
      // is. The angle still drifts slightly with radius, at dphi/dr ~ b / r^2.
      const fullPosition = model.geometry.toCartesianPosition(full.final.position_x);
      const planarPosition = liftPosition(model, reduced.frame, planar.final);
      expect(
        angleBetween(fullPosition, planarPosition),
        `exit positions diverged in angle for ${local}`,
      ).toBeLessThan(1e-5);

      // The structural claim the whole reduction rests on: the full three-dimensional
      // ray never leaves the plane the reduction picked out. Measured as the component
      // of its final position along the plane normal, which must vanish.
      const { normal } = reduced.frame;
      const outOfPlane =
        Math.abs(
          fullPosition[0] * normal[0] + fullPosition[1] * normal[1] + fullPosition[2] * normal[2],
        ) / Math.hypot(...fullPosition);
      expect(outOfPlane, `full 3D ray left its orbital plane for ${local}`).toBeLessThan(1e-9);
    }
  });

  it('keeps sin(theta) at unity, so the polar axis is never approached', () => {
    // The reason the reduction exists. A reduced trace stays exactly equatorial, so
    // cot(theta) never grows and the chart never degenerates.
    const ray = generateNullRay(observer, unit([-0.8, 0.3, 0.5]));
    const reduced = reduceToOrbitalPlane(model, ray);
    const result = integrateGeodesic({
      model,
      integrator: new RKF45Integrator(STATE_DIM, { tolerance: { absolute: 1e-12, relative: 1e-12 } }),
      initial: reduced.initial,
      limits: { initialStep: 1e-3, parameterMax: 5000, maxSteps: 500_000, maxStep: 5 },
      terminator: (position_x) => position_x[1] >= 500,
    });
    // theta drifts only by the rounding of cos(pi/2), which is 6.1e-17 rather than 0.
    expect(Math.abs(result.final.position_x[2] - Math.PI / 2)).toBeLessThan(1e-10);
  });

  it('handles a purely radial ray, which has no preferred plane', () => {
    const ray = generateNullRay(observer, [-1, 0, 0]);
    const reduced = reduceToOrbitalPlane(model, ray);
    expect(reduced.tangentialSpeed).toBeLessThan(1e-13);
    expect(reduced.initial.tangent[3]).toBeLessThan(1e-14);

    // Whatever plane was chosen, it must still contain the ray.
    const lifted = liftDirection(model, reduced.frame, reduced.initial);
    const original = model.geometry.toCartesianDirection(ray.position_x, ray.tangent);
    expect(angleBetween(lifted, original)).toBeLessThan(1e-13);
  });

  it('refuses a model that does not declare spherical symmetry', () => {
    const axisymmetricOnly = overrideModel(model, {
      id: 'pretend-kerr',
      symmetries: { stationary: true, axisymmetric: true, sphericallySymmetric: false },
    });
    const ray = generateNullRay(observer, unit([-1, 0.2, 0.1]));
    expect(() => reduceToOrbitalPlane(axisymmetricOnly, ray)).toThrow(
      /does not declare spherical symmetry/,
    );
  });

  it('refuses a Cartesian chart, where the reduced coordinates would be meaningless', () => {
    // Minkowski is spherically symmetric about the origin, so the symmetry precondition
    // holds — but the reduction writes the ray at (t, r, pi/2, 0), and in a Cartesian
    // chart those slots are x, y and z. Producing a silently wrong state would be far
    // worse than refusing, and a Cartesian chart has no polar axis to avoid anyway.
    const flatObserver = staticObserver(minkowski, [0, 10, 0, 0], 'flat');
    const ray = generateNullRay(flatObserver, unit([-1, 0.3, 0.2]));
    expect(minkowski.symmetries.sphericallySymmetric).toBe(true);
    expect(() => reduceToOrbitalPlane(minkowski, ray)).toThrow(/is cartesian, not/);
  });
});
