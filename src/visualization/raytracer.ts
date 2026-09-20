import type { Vec4 } from '../physics/core/indices.js';
import type { PhaseSpaceState } from '../physics/core/phase-space.js';
import { geodesicDerivative } from '../physics/geodesic/geodesic-system.js';
import {
  integrateGeodesic,
  type IntegrationLimits,
  type TerminationReason,
} from '../physics/geodesic/integrate.js';
import type { Integrator } from '../physics/geodesic/integrators/integrator.js';
import { generateNullRay, localRayDirection, type Observer, type PinholeScreen } from '../physics/observer/observer.js';
import type { SpacetimeModel } from '../physics/spacetimes/spacetime-model.js';
import { normalizationResidual } from '../physics/validation/normalization.js';
import { checkTolerance, NULL_NORMALIZATION_POINTWISE } from '../physics/validation/tolerances.js';
import { summarizeHealth, type NumericalHealth } from '../physics/validation/numeric-health.js';
import { sampleCelestialGrid, type CelestialGrid, type RGB } from './celestial-grid.js';

/**
 * Visualization layer (CLAUDE.md §20): backward null-geodesic ray tracing.
 *
 * This follows the pipeline of CLAUDE.md §9 in order — observer-frame ray generation,
 * conversion into phase-space initial conditions, null-geodesic integration, then
 * intersection with the background. It is ray tracing by null-geodesic integration, not
 * generic ray marching, and the term is used deliberately.
 *
 * Frequency shift, emission and radiometric processing (steps 5 and 6 of §9) are not
 * implemented here: Milestone 1 has no emitting matter and no relative motion, so there
 * is nothing to shift. They arrive with Milestone 3 rather than being approximated now.
 */

export type RayOutcome =
  /** Reached the background sphere; the colour is a genuine sample of the grid. */
  | 'background'
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
}

/**
 * Colour used where no physical result exists.
 *
 * Deliberately a flat, obviously non-physical magenta rather than something that blends
 * into the image: CLAUDE.md §17 and §24 forbid replacing a failed calculation with a
 * visually plausible fake. A failed ray should look wrong.
 */
export const FAILED_RAY_COLOR: RGB = Object.freeze({ r: 255, g: 0, b: 220 });

export interface TraceConfig {
  readonly model: SpacetimeModel;
  readonly integrator: Integrator;
  readonly observer: Observer;
  readonly grid: CelestialGrid;
  readonly limits?: Partial<IntegrationLimits>;
}

/** Squared coordinate distance from the spatial origin, in the active chart. */
function spatialRadiusSquared(position_x: Vec4): number {
  return position_x[1] * position_x[1] + position_x[2] * position_x[2] + position_x[3] * position_x[3];
}

/**
 * Trace one backward null geodesic from the observer to the background.
 *
 * The terminator fires when the ray reaches the background sphere's coordinate radius.
 * The background is then sampled along the ray's outgoing spatial direction rather than
 * along its position vector: the direction is what an observer's line of sight actually
 * maps to, and in a curved spacetime the two differ.
 */
export function traceRay(config: TraceConfig, initial: PhaseSpaceState, derivative = geodesicDerivative(config.model)): RayResult {
  const { model, integrator, grid } = config;
  const radiusSquared = grid.radius * grid.radius;

  const result = integrateGeodesic({
    model,
    integrator,
    initial,
    derivative,
    limits: config.limits,
    terminator: (position_x) => spatialRadiusSquared(position_x) >= radiusSquared,
  });

  const nullResidual = normalizationResidual(model, result.final);

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
    return {
      outcome,
      color: FAILED_RAY_COLOR,
      steps: result.steps,
      nullResidual,
      terminationReason: result.reason,
      final: result.final,
    };
  }

  const k = result.final.tangent;
  const color = sampleCelestialGrid(grid, k[1], k[2], k[3]);

  return {
    outcome,
    color,
    steps: result.steps,
    nullResidual,
    terminationReason: result.reason,
    final: result.final,
  };
}

export interface RenderDiagnostics {
  readonly raysTraced: number;
  readonly raysReachingBackground: number;
  readonly raysFailed: number;
  readonly totalSteps: number;
  /** Largest |g_mu_nu k^mu k^nu| over every traced ray. */
  readonly maxNullResidual: number;
  readonly health: NumericalHealth;
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
  raysReachingBackground: number;
  raysFailed: number;
  totalSteps: number;
  maxNullResidual: number;
}

export function createAccumulator(config: TraceConfig, screen: PinholeScreen): RenderAccumulator {
  return {
    widthPx: screen.widthPx,
    heightPx: screen.heightPx,
    pixels: new Uint8ClampedArray(screen.widthPx * screen.heightPx * 4),
    derivative: geodesicDerivative(config.model),
    raysReachingBackground: 0,
    raysFailed: 0,
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
  const nullCheck = checkTolerance(NULL_NORMALIZATION_POINTWISE, accumulator.maxNullResidual);
  return {
    widthPx: accumulator.widthPx,
    heightPx: accumulator.heightPx,
    pixels: accumulator.pixels,
    diagnostics: {
      raysTraced: accumulator.widthPx * accumulator.heightPx,
      raysReachingBackground: accumulator.raysReachingBackground,
      raysFailed: accumulator.raysFailed,
      totalSteps: accumulator.totalSteps,
      maxNullResidual: accumulator.maxNullResidual,
      health: summarizeHealth([nullCheck], accumulator.raysFailed > 0),
    },
  };
}

/** Render a whole image in one pass. */
export function renderImage(config: TraceConfig, screen: PinholeScreen): RenderResult {
  const accumulator = createAccumulator(config, screen);
  renderBand(config, screen, accumulator, 0, screen.heightPx);
  return finalizeRender(accumulator);
}
