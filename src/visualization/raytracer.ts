import type { Vec4 } from '../physics/core/indices.js';
import type { PhaseSpaceState } from '../physics/core/phase-space.js';
import { geodesicDerivative } from '../physics/geodesic/geodesic-system.js';
import {
  integrateGeodesic,
  type IntegrationLimits,
  type TerminationReason,
} from '../physics/geodesic/integrate.js';
import type { Integrator } from '../physics/geodesic/integrators/integrator.js';
import {
  generateNullRay,
  localRayDirection,
  type Observer,
  type PinholeScreen,
} from '../physics/observer/observer.js';
import type { CartesianVec3, SpacetimeModel } from '../physics/spacetimes/spacetime-model.js';
import { normalizationResidual } from '../physics/validation/normalization.js';
import {
  checkTolerance,
  NULL_NORMALIZATION_TRACED,
  type Tolerance,
} from '../physics/validation/tolerances.js';
import { summarizeHealth, type NumericalHealth } from '../physics/validation/numeric-health.js';
import { sampleCelestialGrid, type CelestialGrid, type RGB } from './celestial-grid.js';
import {
  liftDirection,
  reduceToOrbitalPlane,
  type OrbitalPlaneFrame,
} from './orbital-plane.js';

/**
 * Visualization layer (CLAUDE.md §20): backward null-geodesic ray tracing.
 *
 * This follows the pipeline of CLAUDE.md §9 in order — observer-frame ray generation,
 * conversion into phase-space initial conditions, null-geodesic integration, then
 * intersection with the background. It is ray tracing by null-geodesic integration, not
 * generic ray marching, and the term is used deliberately.
 *
 * The tracer is chart-agnostic: it asks the model for the radial coordinate and for the
 * Cartesian direction of a tangent, so the same code renders flat spacetime in Cartesian
 * coordinates and Schwarzschild in spherical coordinates.
 *
 * Frequency shift, emission and radiometric processing (steps 5 and 6 of §9) are not
 * implemented. Milestones 1 and 2A have no emitting matter and no relative motion
 * between emitter and observer, so there is nothing to shift. They arrive with Milestone
 * 3 rather than being approximated now.
 */

export type RayOutcome =
  /** Reached the background sphere; the colour is a genuine sample of the grid. */
  | 'background'
  /** Fell into the black hole. No background light reaches the observer along this ray. */
  | 'captured'
  /** Ran out of parameter or steps before reaching the background. */
  | 'escaped-budget'
  /** Left the chart's valid domain. */
  | 'domain-exit'
  /** Integration failed numerically. */
  | 'numerical-failure';

export interface RayResult {
  readonly outcome: RayOutcome;
  readonly color: RGB;
  readonly steps: number;
  /** g_mu_nu k^mu k^nu at the end of the trace. Exactly 0 for an ideal null geodesic. */
  readonly nullResidual: number;
  readonly terminationReason: TerminationReason;
  readonly final: PhaseSpaceState;
  /** World-frame propagation direction at the end of the trace. */
  readonly exitDirection: CartesianVec3;
  /** The orbital plane used, when the reduction was applied. */
  readonly frame?: OrbitalPlaneFrame;
}

/**
 * Colour used where no physical result exists.
 *
 * Deliberately a flat, obviously non-physical magenta rather than something that blends
 * into the image: CLAUDE.md §17 and §24 forbid replacing a failed calculation with a
 * visually plausible fake. A failed ray should look wrong.
 */
export const FAILED_RAY_COLOR: RGB = Object.freeze({ r: 255, g: 0, b: 220 });

/**
 * Colour for a captured ray: the black-hole shadow.
 *
 * Unlike the failure colour this one is a real physical result. Tracing backward from
 * the observer, the ray ends on the horizon rather than on the background, so no light
 * from the celestial sphere arrives along that line of sight. The region is dark because
 * nothing in this model emits, not because the calculation gave up — and it is the
 * shadow, which is larger than the horizon, not a picture of the horizon itself.
 */
export const SHADOW_COLOR: RGB = Object.freeze({ r: 0, g: 0, b: 0 });

/**
 * A model-supplied test for a ray that can no longer escape.
 *
 * Kept out of the tracer because it is physics, not rendering policy: for Schwarzschild
 * it is the exact statement that a photon inside r = 3M moving inward cannot turn
 * around. CLAUDE.md §6.3 endorses terminating such rays rather than integrating through
 * the horizon in a chart that does not cover it.
 */
export type CaptureTest = (position_x: Vec4, tangent: Vec4) => boolean;

export interface TraceConfig {
  readonly model: SpacetimeModel;
  readonly integrator: Integrator;
  readonly observer: Observer;
  readonly grid: CelestialGrid;
  readonly limits?: Partial<IntegrationLimits>;
  readonly captureTest?: CaptureTest;
  /**
   * Integrate each ray in its own orbital plane.
   *
   * Valid only for a spherically symmetric model, and checked against the model's
   * declared symmetries. Needed for a spherical chart, whose polar axis would otherwise
   * stall rays that pass near it; pointless for a Cartesian chart, which has no axis.
   */
  readonly orbitalPlaneReduction?: boolean;
  /**
   * Which declared tolerance the image-wide null residual is judged against.
   *
   * Defaults to the reference gate. A preview render that deliberately integrates more
   * loosely passes its own tolerance here, so the reported health describes what was
   * actually asked for rather than failing against a standard it was never meant to meet.
   */
  readonly residualTolerance?: Tolerance;
}

/**
 * Trace one backward null geodesic from the observer to the background.
 *
 * The background is sampled along the ray's outgoing spatial direction rather than along
 * its position vector: the direction is what the observer's line of sight actually maps
 * to, and in a curved spacetime the two differ — which is precisely the lensing the
 * image is meant to show.
 */
export function traceRay(
  config: TraceConfig,
  initial: PhaseSpaceState,
  derivative = geodesicDerivative(config.model),
): RayResult {
  const { model, integrator, grid, captureTest } = config;

  let frame: OrbitalPlaneFrame | undefined;
  let working = initial;
  if (config.orbitalPlaneReduction) {
    const reduced = reduceToOrbitalPlane(model, initial);
    frame = reduced.frame;
    working = reduced.initial;
  }

  let captured = false;
  const result = integrateGeodesic({
    model,
    integrator,
    initial: working,
    derivative,
    limits: config.limits,
    terminator: (position_x, tangent) => {
      if (captureTest?.(position_x, tangent)) {
        captured = true;
        return true;
      }
      return model.geometry.spatialRadius(position_x) >= grid.radius;
    },
  });

  const nullResidual = normalizationResidual(model, result.final);
  const exitDirection = frame
    ? liftDirection(model, frame, result.final)
    : model.geometry.toCartesianDirection(result.final.position_x, result.final.tangent);

  const base = {
    steps: result.steps,
    nullResidual,
    terminationReason: result.reason,
    final: result.final,
    exitDirection,
    ...(frame ? { frame } : {}),
  };

  if (captured) {
    return { ...base, outcome: 'captured', color: SHADOW_COLOR };
  }

  let outcome: RayOutcome;
  switch (result.reason) {
    case 'terminator':
      outcome = 'background';
      break;
    case 'domain-exit':
      outcome = 'domain-exit';
      break;
    case 'numerical-failure':
    case 'step-underflow':
      outcome = 'numerical-failure';
      break;
    case 'parameter-limit':
    case 'step-limit':
      outcome = 'escaped-budget';
      break;
  }

  if (outcome !== 'background') {
    return { ...base, outcome, color: FAILED_RAY_COLOR };
  }

  return {
    ...base,
    outcome,
    color: sampleCelestialGrid(grid, exitDirection[0], exitDirection[1], exitDirection[2]),
  };
}

export interface RenderDiagnostics {
  readonly raysTraced: number;
  readonly raysReachingBackground: number;
  readonly raysFailed: number;
  /** Rays that ended on the horizon: the black-hole shadow. */
  readonly raysCaptured: number;
  readonly totalSteps: number;
  /** Largest |g_mu_nu k^mu k^nu| over every traced ray. */
  readonly maxNullResidual: number;
  readonly health: NumericalHealth;
  /** The tolerance the residual was judged against. */
  readonly residualTolerance: Tolerance;
}

export interface RenderResult {
  readonly widthPx: number;
  readonly heightPx: number;
  /**
   * RGBA8, row-major from the top-left pixel.
   *
   * Pinned to a plain ArrayBuffer rather than the default ArrayBufferLike so the buffer
   * can be handed straight to an ImageData, which does not accept a SharedArrayBuffer.
   */
  readonly pixels: Uint8ClampedArray<ArrayBuffer>;
  readonly diagnostics: RenderDiagnostics;
}

/**
 * Mutable state shared across the bands of one image.
 *
 * Rendering is split into row bands so a browser can yield between them and stay
 * responsive. Splitting changes nothing physical: each ray is traced independently, and
 * the band boundaries are purely a scheduling detail.
 */
export interface RenderAccumulator {
  readonly widthPx: number;
  readonly heightPx: number;
  readonly pixels: Uint8ClampedArray<ArrayBuffer>;
  readonly derivative: ReturnType<typeof geodesicDerivative>;
  readonly residualTolerance: Tolerance;
  raysReachingBackground: number;
  raysFailed: number;
  raysCaptured: number;
  totalSteps: number;
  maxNullResidual: number;
}

export function createAccumulator(config: TraceConfig, screen: PinholeScreen): RenderAccumulator {
  return {
    widthPx: screen.widthPx,
    heightPx: screen.heightPx,
    pixels: new Uint8ClampedArray(screen.widthPx * screen.heightPx * 4),
    derivative: geodesicDerivative(config.model),
    residualTolerance: config.residualTolerance ?? NULL_NORMALIZATION_TRACED,
    raysReachingBackground: 0,
    raysFailed: 0,
    raysCaptured: 0,
    totalSteps: 0,
    maxNullResidual: 0,
  };
}

/**
 * Trace the rows [rowStart, rowEnd) of an image into the accumulator.
 *
 * One backward null geodesic per pixel, following CLAUDE.md §9 in order.
 */
export function renderBand(
  config: TraceConfig,
  screen: PinholeScreen,
  accumulator: RenderAccumulator,
  rowStart: number,
  rowEnd: number,
): void {
  const end = Math.min(rowEnd, screen.heightPx);
  for (let j = Math.max(0, rowStart); j < end; j += 1) {
    for (let i = 0; i < screen.widthPx; i += 1) {
      const direction = localRayDirection(screen, i, j);
      const initial = generateNullRay(config.observer, direction);
      const ray = traceRay(config, initial, accumulator.derivative);

      accumulator.totalSteps += ray.steps;
      if (ray.outcome === 'background') accumulator.raysReachingBackground += 1;
      if (ray.outcome === 'captured') accumulator.raysCaptured += 1;
      if (ray.outcome === 'numerical-failure') accumulator.raysFailed += 1;

      const residual = Math.abs(ray.nullResidual);
      if (Number.isFinite(residual) && residual > accumulator.maxNullResidual) {
        accumulator.maxNullResidual = residual;
      }

      const offset = (j * screen.widthPx + i) * 4;
      accumulator.pixels[offset] = ray.color.r;
      accumulator.pixels[offset + 1] = ray.color.g;
      accumulator.pixels[offset + 2] = ray.color.b;
      accumulator.pixels[offset + 3] = 255;
    }
  }
}

/**
 * Close out an accumulator into a result, checking the image-wide null residual.
 *
 * CLAUDE.md §16 is explicit that a result is not validated merely because it looks
 * right, so the render reports its invariant alongside the picture rather than leaving
 * the picture to speak for itself.
 */
export function finalizeRender(accumulator: RenderAccumulator): RenderResult {
  const nullCheck = checkTolerance(accumulator.residualTolerance, accumulator.maxNullResidual);
  return {
    widthPx: accumulator.widthPx,
    heightPx: accumulator.heightPx,
    pixels: accumulator.pixels,
    diagnostics: {
      raysTraced: accumulator.widthPx * accumulator.heightPx,
      raysReachingBackground: accumulator.raysReachingBackground,
      raysFailed: accumulator.raysFailed,
      raysCaptured: accumulator.raysCaptured,
      totalSteps: accumulator.totalSteps,
      maxNullResidual: accumulator.maxNullResidual,
      health: summarizeHealth([nullCheck], accumulator.raysFailed > 0),
      residualTolerance: accumulator.residualTolerance,
    },
  };
}

/** Render a whole image in one pass. */
export function renderImage(config: TraceConfig, screen: PinholeScreen): RenderResult {
  const accumulator = createAccumulator(config, screen);
  renderBand(config, screen, accumulator, 0, screen.heightPx);
  return finalizeRender(accumulator);
}
