import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { generateNullRay, localRayDirectionAt } from '../../src/physics/observer/observer.js';
import { kerr, outerHorizonRadius } from '../../src/physics/spacetimes/kerr.js';
import { shadowExtent } from '../../src/physics/spacetimes/kerr-shadow.js';
import { renderImage } from '../../src/visualization/raytracer.js';
import { buildScene, type SceneDescription } from '../../src/visualization/scene.js';

/**
 * The rendered Kerr image (ROADMAP.md 4A).
 *
 * The analytic critical curve says where the shadow's edge should be; these tests check
 * that the picture agrees. Nothing here is drawn from the curve — the renderer only ever
 * integrates geodesics and asks the capture criterion what happened — so the comparison
 * is between two independent calculations.
 */

const M = 1;

function scene(spin: number, inclinationDeg: number, widthPx = 121): SceneDescription {
  return {
    kind: 'kerr-sky',
    widthPx,
    heightPx: Math.round((widthPx * 3) / 4),
    samplesPerAxis: 1,
    seed: 1,
    cameraRadius: 60,
    inclinationDeg,
    spin,
  };
}

/** Columns of the shadow, found by looking for captured (black) pixels in the image. */
function shadowColumns(pixels: Uint8ClampedArray, width: number, row: number): number[] {
  const columns: number[] = [];
  for (let i = 0; i < width; i += 1) {
    const offset = (row * width + i) * 4;
    if (pixels[offset] === 0 && pixels[offset + 1] === 0 && pixels[offset + 2] === 0) columns.push(i);
  }
  return columns;
}

describe('a rendered Kerr shadow', () => {
  const built = buildScene(scene(0.9, 90));
  const image = renderImage(built.config, built.screen);

  it('traces without a single numerical failure', () => {
    expect(image.diagnostics.raysFailed).toBe(0);
    expect(image.diagnostics.health.level).toBe('ok');
    expect(image.diagnostics.raysCaptured).toBeGreaterThan(100);
    expect(image.diagnostics.raysReachingBackground).toBeGreaterThan(1000);
  });

  it('holds the null constraint over the whole image', () => {
    expect(image.diagnostics.maxNullResidual).toBeLessThan(image.diagnostics.residualTolerance.value);
  });

  it('puts the shadow where the analytic critical curve says, to within a pixel', () => {
    // The camera looks inward along -r with 'right' along -phi-hat, so the analytic alpha
    // axis runs right to left across the image: a shadow displaced to +alpha appears left
    // of centre. The comparison is of the edges in impact parameter, converted from pixels
    // through the field of view at the camera's radius.
    const { widthPx, heightPx } = built.description;
    const middle = Math.floor(heightPx / 2);
    const columns = shadowColumns(image.pixels, widthPx, middle);
    expect(columns.length).toBeGreaterThan(10);

    const first = columns[0];
    const last = columns[columns.length - 1];
    // Angle of a pixel centre from the image centre, signed along the screen's right axis.
    const angleAt = (i: number): number => {
      const direction = localRayDirectionAt(built.screen, i + 0.5, middle + 0.5);
      // Component along 'right' against the forward component gives the signed angle.
      const right = built.screen.right;
      const forward = built.screen.forward;
      const alongRight = direction[0] * right[0] + direction[1] * right[1] + direction[2] * right[2];
      const alongForward = direction[0] * forward[0] + direction[1] * forward[1] + direction[2] * forward[2];
      return Math.atan2(alongRight, alongForward);
    };
    // Impact parameter b = r sin(psi) / sqrt(-g_tt) is the Schwarzschild relation; for a
    // ZAMO in Kerr at large r the same small-angle relation b ~ r psi holds to O(M/r), so
    // this is compared with a tolerance that accounts for it.
    const r = built.description.kind === 'kerr-sky' ? built.description.cameraRadius : 0;
    const measuredLeft = -angleAt(first) * r;
    const measuredRight = -angleAt(last) * r;

    const analytic = shadowExtent(M, 0.9, Math.PI / 2, 8192);
    const pixelInB = (Math.abs(angleAt(1) - angleAt(0)) * r);
    expect(Math.abs(measuredLeft - analytic.alphaMax), 'left edge of the image').toBeLessThan(2 * pixelInB);
    expect(Math.abs(measuredRight - analytic.alphaMin), 'right edge of the image').toBeLessThan(2 * pixelInB);
  });

  it('agrees with full integration on which pixels are dark', () => {
    // The renderer decides capture from the radial potential's roots, which is exact but
    // is not the same calculation as following the ray in. Here a band of pixels across
    // the shadow edge is integrated with nothing but a horizon terminator, and the two
    // verdicts must match pixel for pixel.
    const { widthPx, heightPx } = built.description;
    const middle = Math.floor(heightPx / 2);
    const model = kerr(M, 0.9);
    const horizon = outerHorizonRadius(M, 0.9);
    const columns = shadowColumns(image.pixels, widthPx, middle);
    const edge = columns[0];

    let checked = 0;
    for (let i = edge - 3; i <= edge + 3; i += 1) {
      if (i < 0 || i >= widthPx) continue;
      const local = localRayDirectionAt(built.screen, i + 0.5, middle + 0.5);
      let reachedHorizon = false;
      integrateGeodesic({
        model,
        integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-11, relative: 1e-11 } }),
        initial: generateNullRay(built.config.observer, local),
        limits: { initialStep: 1e-3, parameterMax: 100_000, maxSteps: 400_000, maxStep: 10 },
        events: [{ id: 'sky', value: (x: Vec4) => x[1] - 400, direction: 1, terminal: true }],
        terminator: (x: Vec4) => (reachedHorizon = x[1] < horizon * 1.0005),
      });
      const offset = (middle * widthPx + i) * 4;
      const painted = image.pixels[offset] === 0 && image.pixels[offset + 1] === 0 && image.pixels[offset + 2] === 0;
      expect(painted, `pixel ${i} of row ${middle}`).toBe(reachedHorizon);
      checked += 1;
    }
    expect(checked).toBeGreaterThanOrEqual(5);
  });

  it('is asymmetric left to right, and mirrors when the spin reverses', () => {
    const { widthPx, heightPx } = built.description;
    const middle = Math.floor(heightPx / 2);
    const centreOf = (pixels: Uint8ClampedArray): number => {
      const columns = shadowColumns(pixels, widthPx, middle);
      return (columns[0] + columns[columns.length - 1]) / 2 - (widthPx - 1) / 2;
    };

    const prograde = centreOf(image.pixels);
    expect(Math.abs(prograde), 'displacement in pixels').toBeGreaterThan(1.5);

    const reversed = buildScene(scene(-0.9, 90));
    const mirrored = renderImage(reversed.config, reversed.screen);
    expect(centreOf(mirrored.pixels)).toBeCloseTo(-prograde, 6);
  });

  it('loses the displacement when the hole is viewed down its spin axis', () => {
    const axial = buildScene(scene(0.9, 2, 81));
    const rendered = renderImage(axial.config, axial.screen);
    const { widthPx, heightPx } = axial.description;
    const middle = Math.floor(heightPx / 2);
    const columns = shadowColumns(rendered.pixels, widthPx, middle);
    expect(columns.length).toBeGreaterThan(5);
    const centre = (columns[0] + columns[columns.length - 1]) / 2 - (widthPx - 1) / 2;
    expect(Math.abs(centre)).toBeLessThan(1);
  });
});

describe('the background radius the Kerr scene samples at', () => {
  it('leaves a bias well under one pixel', () => {
    // Kerr has no closed-form tail integral, so the sky is sampled along the local
    // direction at 400M. This measures what that costs, by integrating the same rays on to
    // 4000M and comparing the directions: the Schwarzschild bias falls as 1/R^2, and the
    // same behaviour here puts it far below the angular size of a pixel.
    const built = buildScene(scene(0.9, 90));
    const model = kerr(M, 0.9);
    const integrator = (): Dopri5Integrator =>
      new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-11, relative: 1e-11 } });

    const directionAt = (i: number, j: number, radius: number): readonly [number, number, number] => {
      const local = localRayDirectionAt(built.screen, i + 0.5, j + 0.5);
      const result = integrateGeodesic({
        model,
        integrator: integrator(),
        initial: generateNullRay(built.config.observer, local),
        limits: { initialStep: 1e-3, parameterMax: 200_000, maxSteps: 400_000, maxStep: 20 },
        events: [{ id: 'sky', value: (x: Vec4) => x[1] - radius, direction: 1, terminal: true }],
      });
      const world = model.geometry.toCartesianDirection(result.final.position_x, result.final.tangent);
      const norm = Math.hypot(world[0], world[1], world[2]);
      return [world[0] / norm, world[1] / norm, world[2] / norm];
    };

    const pixelAngle = built.screen.horizontalFovRad / built.screen.widthPx;
    let worst = 0;
    for (const [i, j] of [[20, 45], [95, 45], [60, 12], [30, 70]] as const) {
      const near = directionAt(i, j, 400);
      const far = directionAt(i, j, 4000);
      const dot = near[0] * far[0] + near[1] * far[1] + near[2] * far[2];
      worst = Math.max(worst, Math.acos(Math.min(1, dot)));
    }
    expect(worst / pixelAngle, 'bias in pixels').toBeLessThan(0.1);
  });
});
