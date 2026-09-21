import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { RKF45Integrator } from '../../src/physics/geodesic/integrators/rkf45.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import {
  generateNullRay,
  inwardFacingScreen,
  staticObserver,
} from '../../src/physics/observer/observer.js';
import { orthonormalityResidual } from '../../src/physics/observer/tetrad.js';
import {
  criticalImpactParameter,
  schwarzschild,
} from '../../src/physics/spacetimes/schwarzschild.js';
import { isCaptured } from '../../src/physics/spacetimes/schwarzschild-rays.js';
import {
  NULL_NORMALIZATION_TRACED,
  TETRAD_ORTHONORMALITY,
  checkTolerance,
} from '../../src/physics/validation/tolerances.js';
import { DEFAULT_CELESTIAL_GRID } from '../../src/visualization/celestial-grid.js';
import { renderImage, SHADOW_COLOR, traceRay, type TraceConfig } from '../../src/visualization/raytracer.js';

/**
 * ROADMAP.md 2A.3 — the CPU reference raytracer, validated against analysis.
 *
 * The headline result is the black-hole shadow, and its angular size is known in closed
 * form. For a static observer at radius r, a photon arriving at angle psi from the
 * inward radial direction has impact parameter b = r sin(psi) / sqrt(f), so the shadow
 * boundary — the locus b = b_c — sits at
 *
 *   sin(psi_shadow) = b_c sqrt(f(r)) / r,   b_c = 3 sqrt(3) M
 *
 * That is a prediction the renderer cannot fudge: it depends on the tetrad, the ray
 * generation, the orbital-plane reduction, the connection, the integrator and the
 * capture test all being right together.
 */

const M = 1;
const model = schwarzschild(M);
const CAMERA_RADIUS = 20;
const camera: Vec4 = [0, CAMERA_RADIUS, Math.PI / 2, 0];
const observer = staticObserver(model, camera);

/** The analytic angular radius of the shadow for a static observer at `r`. */
function shadowAngularRadius(r: number): number {
  const f = 1 - (2 * M) / r;
  return Math.asin((criticalImpactParameter(M) * Math.sqrt(f)) / r);
}

const config: TraceConfig = {
  model,
  // 1e-12 rather than something looser, for a measured reason. At 1e-11 the worst ray
  // in a full render — one passing just outside the capture boundary, which spirals
  // several times near the photon sphere before escaping — accumulates a null residual
  // of 1.5e-9, above the declared NULL_NORMALIZATION_TRACED gate of 1e-9. The honest
  // response is to integrate more tightly, not to widen the gate to fit the result
  // (CLAUDE.md §17). At 1e-12 the worst residual over the same image is 1.4e-10.
  integrator: new RKF45Integrator(STATE_DIM, {
    tolerance: { absolute: 1e-12, relative: 1e-12 },
  }),
  observer,
  grid: { ...DEFAULT_CELESTIAL_GRID, radius: 400 },
  limits: { initialStep: 1e-3, parameterMax: 20_000, maxSteps: 200_000, maxStep: 5 },
  captureTest: (position_x, tangent) => isCaptured(model, position_x, tangent),
  orbitalPlaneReduction: true,
};

/** A local direction at angle `psi` from the inward radial direction. */
function directionAtAngle(psi: number): [number, number, number] {
  // Frame legs are (1) radial, (2) polar, (3) azimuthal. Inward is -leg(1).
  return [-Math.cos(psi), 0, -Math.sin(psi)];
}

function isRayCaptured(psi: number): boolean {
  return traceRay(config, generateNullRay(observer, directionAtAngle(psi))).outcome === 'captured';
}

describe('Schwarzschild observer and ray generation', () => {
  it('builds an orthonormal static tetrad in the curved metric', () => {
    for (const r of [3.5, 6, 20, 1000]) {
      const metric = model.metricAt([0, r, Math.PI / 2, 0]);
      const tetrad = staticObserver(model, [0, r, Math.PI / 2, 0]).tetrad;
      const check = checkTolerance(TETRAD_ORTHONORMALITY, orthonormalityResidual(metric, tetrad));
      expect(check.withinTolerance, `${check.message} at r = ${r}M`).toBe(true);
    }
  });

  it('gives the static observer the shell four-velocity u^mu = (1/sqrt(f), 0, 0, 0)', () => {
    const f = model.lapseFunction(CAMERA_RADIUS);
    const u = observer.tetrad.four_velocity_u();
    expect(u[0]).toBeCloseTo(1 / Math.sqrt(f), 14);
    expect(u[1]).toBe(0);
    expect(u[2]).toBe(0);
    expect(u[3]).toBe(0);
    expect(observer.kind).toBe('physical-observer');
    // It is not a freely-falling frame, and the description must not imply otherwise.
    expect(observer.description).toContain('not a freely-falling frame');
  });

  it('refuses to place a static observer inside the horizon', () => {
    expect(() => staticObserver(model, [0, 1.5 * M, Math.PI / 2, 0])).toThrow(RangeError);
  });

  it('generates rays whose impact parameter matches b = r sin(psi) / sqrt(f)', () => {
    const f = model.lapseFunction(CAMERA_RADIUS);
    for (const psi of [0.05, 0.2, 0.4, 0.8]) {
      const ray = generateNullRay(observer, directionAtAngle(psi));
      // E = -p_t = f k^t, and L_z = r^2 sin^2(theta) k^phi. The generated ray is
      // past-directed, so E is negative; the impact parameter is the ratio's magnitude.
      const energy = f * ray.tangent[0];
      const angular = CAMERA_RADIUS * CAMERA_RADIUS * ray.tangent[3];
      const b = Math.abs(angular / energy);
      const expected = (CAMERA_RADIUS * Math.sin(psi)) / Math.sqrt(f);
      expect(Math.abs(b - expected) / expected, `psi = ${psi}`).toBeLessThan(1e-12);
    }
  });
});

describe('black-hole shadow (ROADMAP.md 2A.3)', () => {
  it('puts the shadow boundary at the analytic angular radius', () => {
    const expected = shadowAngularRadius(CAMERA_RADIUS);
    // Sanity: asin(3 sqrt(3) sqrt(0.9) / 20), about 14.27 degrees from a camera at r = 20M.
    expect(expected).toBeCloseTo(0.2490415079, 9);

    let lo = expected * 0.5;
    let hi = expected * 1.5;
    expect(isRayCaptured(lo), 'a ray well inside the shadow should be captured').toBe(true);
    expect(isRayCaptured(hi), 'a ray well outside the shadow should escape').toBe(false);

    for (let i = 0; i < 40 && hi - lo > 1e-12; i += 1) {
      const mid = 0.5 * (lo + hi);
      if (isRayCaptured(mid)) lo = mid;
      else hi = mid;
    }

    const measured = 0.5 * (lo + hi);
    expect(
      Math.abs(measured - expected) / expected,
      `shadow boundary: traced ${measured}, analytic ${expected}`,
    ).toBeLessThan(1e-6);
  });

  it('shrinks the shadow, in angle, as the observer moves away', () => {
    // The shadow subtends a smaller angle from further out, while its impact-parameter
    // radius b_c is fixed. Checked against the closed form at each radius.
    let previous = Infinity;
    for (const r of [10, 20, 50, 200]) {
      const angle = shadowAngularRadius(r);
      expect(angle).toBeLessThan(previous);
      previous = angle;
    }
    // Far away the angular radius approaches b_c / r.
    expect(shadowAngularRadius(1e6)).toBeCloseTo(criticalImpactParameter(M) / 1e6, 10);
  });

  it('traces rays that stay null all the way out', () => {
    for (const psi of [0.3, 0.5, 0.9, 1.2]) {
      const ray = traceRay(config, generateNullRay(observer, directionAtAngle(psi)));
      expect(ray.outcome, `psi = ${psi}`).toBe('background');
      const check = checkTolerance(NULL_NORMALIZATION_TRACED, ray.nullResidual);
      expect(check.withinTolerance, `psi = ${psi}: ${check.message}`).toBe(true);
    }
  });

  it('bends light: a ray outside the shadow does not arrive from where it points', () => {
    // The point of the whole exercise. A ray just outside the capture boundary is
    // deflected by a large angle, so the background it samples is nowhere near its
    // launch direction. In flat space the two would coincide exactly.
    const psi = shadowAngularRadius(CAMERA_RADIUS) * 1.02;
    const launch = directionAtAngle(psi);
    const ray = traceRay(config, generateNullRay(observer, launch));
    expect(ray.outcome).toBe('background');

    const exit = ray.exitDirection;
    const norm = Math.hypot(...exit);
    // The launch direction is expressed in the observer's frame and the exit direction
    // in world axes, but at this camera position the radial leg is the +x axis, so the
    // inward launch direction maps to the same axes up to the frame's sign conventions.
    // What matters is that the ray has been turned through a large angle overall.
    const deflectionProxy = (exit[0] / norm) * -launch[0] + (exit[2] / norm) * -launch[2];
    expect(Math.abs(deflectionProxy)).toBeLessThan(0.99);
  });
});

describe('rendered Schwarzschild image', () => {
  it('renders a centred, roughly circular shadow of the analytic size', () => {
    const fov = 4 * shadowAngularRadius(CAMERA_RADIUS);
    const screen = inwardFacingScreen(41, 41, fov);
    const rendered = renderImage(config, screen);

    expect(rendered.diagnostics.raysFailed).toBe(0);
    expect(rendered.diagnostics.raysCaptured).toBeGreaterThan(0);
    expect(rendered.diagnostics.health.level).toBe('ok');

    const isShadow = (i: number, j: number): boolean => {
      const o = (j * screen.widthPx + i) * 4;
      return (
        rendered.pixels[o] === SHADOW_COLOR.r &&
        rendered.pixels[o + 1] === SHADOW_COLOR.g &&
        rendered.pixels[o + 2] === SHADOW_COLOR.b
      );
    };

    const centre = (screen.widthPx - 1) / 2;
    expect(isShadow(centre, centre), 'the image centre should be inside the shadow').toBe(true);
    expect(isShadow(0, 0), 'the image corner should be outside the shadow').toBe(false);

    // Measure the shadow radius along four rays from the centre. A spherically symmetric
    // black hole seen by a static observer gives a circular shadow, so all four agree.
    const radii: number[] = [];
    for (const [di, dj] of [
      [1, 0],
      [-1, 0],
      [0, 1],
      [0, -1],
    ] as const) {
      let n = 0;
      while (isShadow(centre + di * (n + 1), centre + dj * (n + 1))) n += 1;
      radii.push(n);
    }
    for (const r of radii) expect(r).toBe(radii[0]);

    // And that radius matches the analytic angular size, to within the pixel grid.
    const tanHalfFov = Math.tan(fov / 2);
    const measuredAngle = Math.atan((radii[0] / centre) * tanHalfFov);
    const pixelAngle = fov / screen.widthPx;
    expect(Math.abs(measuredAngle - shadowAngularRadius(CAMERA_RADIUS))).toBeLessThan(pixelAngle);
  });

  it('shows the lensed background outside the shadow, not a blank field', () => {
    const fov = 4 * shadowAngularRadius(CAMERA_RADIUS);
    const screen = inwardFacingScreen(33, 33, fov);
    const rendered = renderImage(config, screen);

    const colours = new Set<string>();
    for (let p = 0; p < rendered.pixels.length; p += 4) {
      colours.add(`${rendered.pixels[p]},${rendered.pixels[p + 1]},${rendered.pixels[p + 2]}`);
    }
    // Shadow, grid lines and at least one background cell shade.
    expect(colours.size).toBeGreaterThanOrEqual(3);
    expect(colours.has('0,0,0')).toBe(true);
    expect(rendered.diagnostics.raysReachingBackground).toBeGreaterThan(0);
    expect(
      rendered.diagnostics.raysCaptured + rendered.diagnostics.raysReachingBackground,
    ).toBe(rendered.diagnostics.raysTraced);
  });
});
