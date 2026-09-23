import { describe, expect, it } from 'vitest';
import { HAMILTONIAN } from '../../src/physics/geodesic/formulation.js';
import { RK4Integrator } from '../../src/physics/geodesic/integrators/rk4.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import {
  defaultScreen,
  generateNullRay,
  localRayDirection,
  staticMinkowskiObserver,
} from '../../src/physics/observer/observer.js';
import { minkowski } from '../../src/physics/spacetimes/minkowski.js';
import {
  DEFAULT_CELESTIAL_GRID,
  directionToSkyAngles,
  sampleCelestialGrid,
} from '../../src/visualization/celestial-grid.js';
import { renderImage, traceRay, type TraceConfig } from '../../src/visualization/raytracer.js';
import {
  NULL_NORMALIZATION_POINTWISE,
  checkTolerance,
} from '../../src/physics/validation/tolerances.js';

/**
 * ROADMAP.md 1.5 and the Milestone 1 Definition of Done: the canvas renders an
 * undistorted grid via backward ray tracing.
 *
 * "Undistorted" is asserted numerically, not judged by eye. In flat spacetime a null
 * geodesic is a straight line, so every backward-traced ray must arrive at the
 * background along exactly the direction it left the observer, and the traced image
 * must agree pixel-for-pixel with a direct analytic projection of the same grid.
 *
 * CLAUDE.md §16 is explicit that a result is not validated because it looks right.
 * This test is the reason the M1 render can be trusted as a working pipeline rather
 * than a picture that happens to be plausible.
 */

const observer = staticMinkowskiObserver([0, 0, 0, 0]);
const grid = DEFAULT_CELESTIAL_GRID;

const config: TraceConfig = {
  model: minkowski,
  integrator: new RK4Integrator(STATE_DIM),
  observer,
  grid,
  limits: { initialStep: 0.5, parameterMax: 1000, maxSteps: 10_000 },
};

describe('backward ray tracing in flat space (ROADMAP.md 1.5)', () => {
  it('returns every ray along its original direction, with no deflection', () => {
    const screen = defaultScreen(33, 25, Math.PI / 2);
    const session = HAMILTONIAN.bind(minkowski);
    let worstChord = 0;

    for (let j = 0; j < screen.heightPx; j += 1) {
      for (let i = 0; i < screen.widthPx; i += 1) {
        const direction = localRayDirection(screen, i, j);
        const initial = generateNullRay(observer, direction);
        const ray = traceRay(config, initial, session);

        expect(ray.outcome, `pixel (${i}, ${j}) did not reach the background`).toBe('background');

        // In flat space dk^mu/dlambda = 0 exactly, so the wavevector must come back
        // bit-for-bit unchanged. This is the strongest available statement, stronger
        // than any tolerance.
        expect(
          Array.from(ray.final.tangent),
          `pixel (${i}, ${j}) wavevector changed`,
        ).toEqual([-1, ...direction]);

        // Deflection is measured as the chord between the unit directions rather than
        // as acos of their dot product. Near-parallel vectors make acos severely
        // ill-conditioned: a dot product of 1 - 1e-16 gives acos ~ sqrt(2e-16) ~ 1.5e-8,
        // which reads as a deflection of 1e-8 radians when the true deflection is zero.
        // The chord is well conditioned and equals the angle to leading order.
        const k = ray.final.tangent;
        const norm = Math.hypot(k[1], k[2], k[3]);
        const chord = Math.hypot(
          k[1] / norm - direction[0],
          k[2] / norm - direction[1],
          k[3] / norm - direction[2],
        );
        if (chord > worstChord) worstChord = chord;
      }
    }

    // A straight line in flat spacetime: the deflection is zero to the limit of binary64
    // arithmetic. Anything larger means the pipeline is bending light that the geometry
    // does not bend.
    expect(worstChord).toBeLessThan(1e-15);
  });

  it('keeps every traced ray null all the way to the background', () => {
    const screen = defaultScreen(21, 17, Math.PI / 2);
    const session = HAMILTONIAN.bind(minkowski);

    for (let j = 0; j < screen.heightPx; j += 1) {
      for (let i = 0; i < screen.widthPx; i += 1) {
        const initial = generateNullRay(observer, localRayDirection(screen, i, j));
        const ray = traceRay(config, initial, session);
        const check = checkTolerance(NULL_NORMALIZATION_POINTWISE, ray.nullResidual);
        expect(check.withinTolerance, `${check.message} at pixel (${i}, ${j})`).toBe(true);
      }
    }
  });

  it('reproduces a direct analytic projection of the grid, pixel for pixel', () => {
    const screen = defaultScreen(48, 36, Math.PI / 2);
    const rendered = renderImage(config, screen);

    // The reference: sample the background along each pixel's initial viewing direction,
    // with no integration at all. In flat space the traced image must equal this exactly.
    const reference = new Uint8ClampedArray(screen.widthPx * screen.heightPx * 4);
    for (let j = 0; j < screen.heightPx; j += 1) {
      for (let i = 0; i < screen.widthPx; i += 1) {
        const d = localRayDirection(screen, i, j);
        const color = sampleCelestialGrid(grid, d[0], d[1], d[2]);
        const offset = (j * screen.widthPx + i) * 4;
        reference[offset] = color.r;
        reference[offset + 1] = color.g;
        reference[offset + 2] = color.b;
        reference[offset + 3] = 255;
      }
    }

    let differingPixels = 0;
    for (let p = 0; p < reference.length; p += 4) {
      if (
        rendered.pixels[p] !== reference[p] ||
        rendered.pixels[p + 1] !== reference[p + 1] ||
        rendered.pixels[p + 2] !== reference[p + 2]
      ) {
        differingPixels += 1;
      }
    }
    expect(differingPixels).toBe(0);
  });

  it('renders a grid that actually contains structure', () => {
    // Guards the test above: a uniform image would match a uniform reference trivially
    // and prove nothing about the pipeline.
    const screen = defaultScreen(48, 36, Math.PI / 2);
    const rendered = renderImage(config, screen);

    const distinctColors = new Set<string>();
    for (let p = 0; p < rendered.pixels.length; p += 4) {
      distinctColors.add(`${rendered.pixels[p]},${rendered.pixels[p + 1]},${rendered.pixels[p + 2]}`);
    }
    expect(distinctColors.size).toBeGreaterThanOrEqual(3);

    const lineColor = `${grid.lineColor.r},${grid.lineColor.g},${grid.lineColor.b}`;
    expect(distinctColors.has(lineColor), 'no grid lines appear in the render').toBe(true);
  });

  it('reports healthy diagnostics with no failed rays', () => {
    const screen = defaultScreen(32, 24, Math.PI / 2);
    const { diagnostics } = renderImage(config, screen);

    expect(diagnostics.raysTraced).toBe(32 * 24);
    expect(diagnostics.raysReachingBackground).toBe(diagnostics.raysTraced);
    expect(diagnostics.raysFailed).toBe(0);
    expect(diagnostics.totalSteps).toBeGreaterThan(0);
    expect(diagnostics.health.level).toBe('ok');

    const check = checkTolerance(NULL_NORMALIZATION_POINTWISE, diagnostics.maxNullResidual);
    expect(check.withinTolerance, check.message).toBe(true);
  });

  it('marks rays that never reach the background as failures, not as sky', () => {
    // A budget too small to reach the background sphere. The renderer must not colour
    // these as though they had landed somewhere (CLAUDE.md §17, §24).
    const starved: TraceConfig = {
      ...config,
      limits: { initialStep: 0.5, parameterMax: 1, maxSteps: 10 },
    };
    const ray = traceRay(starved, generateNullRay(observer, [1, 0, 0]));
    expect(ray.outcome).toBe('escaped-budget');
    expect(ray.color).toEqual({ r: 255, g: 0, b: 220 });
  });
});

describe('rectilinear projection (why grid lines curve without any deflection)', () => {
  // The rendered grid shows curved parallels. That is projection geometry, not lensing,
  // and the UI says so. These tests hold the projection to the property that claim rests
  // on: a rectilinear camera maps great circles, and only great circles, to straight
  // image lines.
  const screen = defaultScreen(65, 49, Math.PI / 2);

  /** Unit normal of the plane through the origin best fitting a set of directions. */
  function planarityResidual(directions: readonly (readonly [number, number, number])[]): number {
    const [a, b] = [directions[0], directions[directions.length - 1]];
    const normal: [number, number, number] = [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
    ];
    const norm = Math.hypot(...normal);
    let worst = 0;
    for (const d of directions) {
      const out = Math.abs((d[0] * normal[0] + d[1] * normal[1] + d[2] * normal[2]) / norm);
      if (out > worst) worst = out;
    }
    return worst;
  }

  it('maps every image row to a great circle on the sky', () => {
    for (const j of [0, 12, 24, 36, 48]) {
      const directions = Array.from({ length: screen.widthPx }, (_, i) =>
        localRayDirection(screen, i, j),
      );
      // Coplanar through the origin means the directions trace a great circle, which is
      // exactly why a straight image row is not a line of constant polar angle.
      expect(planarityResidual(directions), `row ${j} is not a great circle`).toBeLessThan(1e-14);
    }
  });

  it('maps every image column to a great circle on the sky', () => {
    for (const i of [0, 16, 32, 48, 64]) {
      const directions = Array.from({ length: screen.heightPx }, (_, j) =>
        localRayDirection(screen, i, j),
      );
      expect(planarityResidual(directions), `column ${i} is not a great circle`).toBeLessThan(1e-14);
    }
  });

  it('shows a parallel that is not a great circle, so it must project curved', () => {
    // Directions at a fixed polar angle away from the equator: a small circle. Its
    // failure to be coplanar with the origin is precisely why it renders as a curve.
    const theta = Math.PI / 4;
    const smallCircle = Array.from({ length: 16 }, (_, n) => {
      const phi = (n / 16) * Math.PI - Math.PI / 2;
      return [
        Math.sin(theta) * Math.cos(phi),
        Math.sin(theta) * Math.sin(phi),
        Math.cos(theta),
      ] as [number, number, number];
    });
    expect(planarityResidual(smallCircle)).toBeGreaterThan(0.1);
  });
});

describe('celestial grid sampling', () => {
  it('maps the +z axis to the north pole and +x to zero azimuth', () => {
    expect(directionToSkyAngles(0, 0, 1).theta).toBeCloseTo(0, 14);
    expect(directionToSkyAngles(0, 0, -1).theta).toBeCloseTo(Math.PI, 14);
    const equator = directionToSkyAngles(1, 0, 0);
    expect(equator.theta).toBeCloseTo(Math.PI / 2, 14);
    expect(equator.phi).toBeCloseTo(0, 14);
  });

  it('returns azimuth in [0, 2pi)', () => {
    const behind = directionToSkyAngles(-1, -1e-9, 0);
    expect(behind.phi).toBeGreaterThan(Math.PI);
    expect(behind.phi).toBeLessThan(2 * Math.PI);
  });

  it('is scale invariant: only the direction matters', () => {
    const near = sampleCelestialGrid(grid, 1, 2, 3);
    const far = sampleCelestialGrid(grid, 1000, 2000, 3000);
    expect(near).toEqual(far);
  });

  it('refuses a zero-length direction rather than returning an arbitrary colour', () => {
    expect(() => directionToSkyAngles(0, 0, 0)).toThrow(RangeError);
    expect(() => sampleCelestialGrid(grid, 0, 0, 0)).toThrow(RangeError);
  });
});
