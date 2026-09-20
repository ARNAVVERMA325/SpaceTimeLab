import { describe, expect, it } from 'vitest';
import { contractLower } from '../../src/physics/core/metric-tensor.js';
import {
  defaultScreen,
  generateNullRay,
  localRayDirection,
  staticMinkowskiObserver,
} from '../../src/physics/observer/observer.js';
import { frameToCoordinate, orthonormalityResidual, Tetrad } from '../../src/physics/observer/tetrad.js';
import { minkowski } from '../../src/physics/spacetimes/minkowski.js';
import {
  NULL_NORMALIZATION_POINTWISE,
  TETRAD_ORTHONORMALITY,
  checkTolerance,
} from '../../src/physics/validation/tolerances.js';

/**
 * Observer layer (CLAUDE.md §4, ROADMAP.md 1.5 groundwork).
 *
 * The tetrad is validated against g_mu_nu e^mu_(a) e^nu_(b) = eta_ab rather than
 * assumed correct, and every generated ray is checked to be genuinely null.
 */

describe('tetrad frames', () => {
  it('is orthonormal against the Minkowski metric', () => {
    const metric = minkowski.metricAt([0, 0, 0, 0]);
    const check = checkTolerance(TETRAD_ORTHONORMALITY, orthonormalityResidual(metric, Tetrad.identity()));
    expect(check.withinTolerance, check.message).toBe(true);
    // In Cartesian Minkowski the identity frame is orthonormal exactly, not approximately.
    expect(orthonormalityResidual(metric, Tetrad.identity())).toBe(0);
  });

  it('detects a frame that is not orthonormal', () => {
    const metric = minkowski.metricAt([0, 0, 0, 0]);
    const skewed = Tetrad.identity();
    skewed.components[4 * 1 + 2] = 0.5; // leg (1) picks up a component along y
    expect(orthonormalityResidual(metric, skewed)).toBeGreaterThan(0.1);
  });

  it('exposes leg (0) as the observer four-velocity', () => {
    const tetrad = Tetrad.identity();
    expect(Array.from(tetrad.four_velocity_u())).toEqual([1, 0, 0, 0]);
    const metric = minkowski.metricAt([0, 0, 0, 0]);
    // u^mu is timelike and normalized to -1.
    expect(contractLower(metric, tetrad.four_velocity_u(), tetrad.four_velocity_u())).toBe(-1);
  });

  it('maps frame components to coordinate components through e^mu_(a)', () => {
    // Identity frame: coordinate components equal frame components.
    expect(Array.from(frameToCoordinate(Tetrad.identity(), [1, 2, 3, 4]))).toEqual([1, 2, 3, 4]);
  });

  it('rejects a wrongly sized tetrad buffer', () => {
    expect(() => new Tetrad(new Float64Array(12))).toThrow(RangeError);
  });
});

describe('pinhole screen ray generation (CLAUDE.md §9 steps 1-2)', () => {
  const screen = defaultScreen(64, 48, Math.PI / 3);

  it('produces unit directions for every pixel', () => {
    for (let j = 0; j < screen.heightPx; j += 7) {
      for (let i = 0; i < screen.widthPx; i += 7) {
        const d = localRayDirection(screen, i, j);
        expect(Math.hypot(...d)).toBeCloseTo(1, 14);
      }
    }
  });

  it('points the centre of the image along the screen forward axis', () => {
    const even = defaultScreen(64, 48, Math.PI / 3);
    // With an even pixel count no pixel sits exactly at the centre, so the two central
    // columns must straddle it symmetrically.
    const left = localRayDirection(even, 31, 23);
    const right = localRayDirection(even, 32, 23);
    expect(left[1]).toBeCloseTo(-right[1], 14);
    expect(left[0]).toBeCloseTo(right[0], 14);
  });

  it('spans the requested horizontal field of view', () => {
    const first = localRayDirection(screen, 0, 24);
    const last = localRayDirection(screen, screen.widthPx - 1, 24);
    const angle = Math.acos(Math.min(1, first[0] * last[0] + first[1] * last[1] + first[2] * last[2]));
    // Pixel centres sit half a pixel inside each edge, so the spanned angle is slightly
    // under the nominal field of view.
    const pixelAngle = screen.horizontalFovRad / screen.widthPx;
    expect(angle).toBeLessThan(screen.horizontalFovRad);
    expect(angle).toBeGreaterThan(screen.horizontalFovRad - 2 * pixelAngle);
  });

  it('generates wavevectors that are genuinely null', () => {
    const observer = staticMinkowskiObserver([0, 0, 0, 0]);
    const metric = minkowski.metricAt(observer.position_x);
    for (let j = 0; j < screen.heightPx; j += 5) {
      for (let i = 0; i < screen.widthPx; i += 5) {
        const ray = generateNullRay(observer, localRayDirection(screen, i, j));
        expect(ray.kind).toBe('null');
        const check = checkTolerance(
          NULL_NORMALIZATION_POINTWISE,
          contractLower(metric, ray.tangent, ray.tangent),
        );
        expect(check.withinTolerance, `${check.message} at pixel (${i}, ${j})`).toBe(true);
      }
    }
  });

  it('generates past-directed rays, as backward tracing requires', () => {
    const observer = staticMinkowskiObserver([0, 0, 0, 0]);
    const ray = generateNullRay(observer, [1, 0, 0]);
    // k^0 < 0: as the affine parameter increases, coordinate time decreases and the ray
    // runs back towards the photon's source.
    expect(ray.tangent[0]).toBeLessThan(0);
    expect(ray.tangent[1]).toBeGreaterThan(0);
  });

  it('labels a physical observer as such (CLAUDE.md §4)', () => {
    expect(staticMinkowskiObserver([0, 0, 0, 0]).kind).toBe('physical-observer');
  });

  it('rejects a degenerate screen basis', () => {
    const degenerate = { ...screen, forward: [0, 0, 0] as const, right: [0, 0, 0] as const, up: [0, 0, 0] as const };
    expect(() => localRayDirection(degenerate, 32, 24)).toThrow(RangeError);
  });
});
