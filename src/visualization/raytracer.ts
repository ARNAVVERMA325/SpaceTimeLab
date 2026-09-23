import type { Vec4 } from '../physics/core/indices.js';
import type { PhaseSpaceState } from '../physics/core/phase-space.js';
import {
  HAMILTONIAN,
  type FormulationSession,
  type GeodesicFormulation,
} from '../physics/geodesic/formulation.js';
import {
  integrateGeodesic,
  type IntegrationLimits,
  type TerminationReason,
} from '../physics/geodesic/integrate.js';
import type { Integrator } from '../physics/geodesic/integrators/integrator.js';
import {
  generateNullRay,
  localRayDirectionAt,
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
  decodeSrgb,
  hashToUnit,
  LINEAR_BLACK,
  linearToSrgb8,
  type LinearRGB,
} from './color.js';
import {
  liftDirection,
  liftPosition,
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
  /**
   * Direction used for the background lookup, in world axes: the direction at infinity
   * when the scene supplies an asymptotic correction, otherwise the local direction at
   * the background radius.
   */
  readonly exitDirection: CartesianVec3;
  /** Whether `exitDirection` was corrected to infinity. */
  readonly asymptoticallyCorrected: boolean;
  /** This sample's contribution in linear light. */
  readonly radiance: LinearRGB;
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
  /** How the geodesic equation is written. Defaults to Hamiltonian. */
  readonly formulation?: GeodesicFormulation;
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
  /**
   * Map a ray's position and local direction at the background radius to its direction
   * at infinity. Without it the background is sampled along the local direction, which
   * carries the deflection still to come beyond that radius as a systematic bias.
   */
  readonly asymptoticDirection?: (positionWorld: CartesianVec3, directionWorld: CartesianVec3) => CartesianVec3;
  /** Rays per pixel. Defaults to one ray through each pixel centre. */
  readonly sampling?: SamplingPlan;
}

/**
 * Stratified supersampling (CLAUDE.md §9, step 7: image accumulation / anti-aliasing).
 *
 * Each pixel is divided into n x n cells and one ray is traced through a jittered point
 * in each. Every sample is an independent, fully integrated geodesic; accumulation only
 * averages their results, in linear light. No trajectory is altered, interpolated or
 * smoothed — which is the line CLAUDE.md §9 draws between rendering operations and
 * physics.
 *
 * Jitter is a deterministic hash of (pixel, sample, seed), so renders are reproducible
 * bit for bit.
 */
export interface SamplingPlan {
  readonly samplesPerAxis: number;
  readonly seed?: number;
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
  session: FormulationSession = (config.formulation ?? HAMILTONIAN).bind(config.model),
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
    session,
    limits: config.limits,
    events: [
      {
        // Arrival at the background sphere, located exactly rather than detected one
        // step late: the ray is placed on r = R, so the direction read off it belongs to
        // that radius and not to wherever the last step overshot to.
        id: 'background',
        value: (position_x) => model.geometry.spatialRadius(position_x) - grid.radius,
        direction: 1,
        terminal: true,
      },
    ],
    terminator: captureTest
      ? (position_x, tangent) => {
          if (captureTest(position_x, tangent)) {
            captured = true;
            return true;
          }
          return false;
        }
      : undefined,
  });

  const nullResidual = normalizationResidual(model, result.final);
  const localDirection = frame
    ? liftDirection(model, frame, result.final)
    : model.geometry.toCartesianDirection(result.final.position_x, result.final.tangent);

  const reachedBackground = !captured && result.reason === 'event';
  let exitDirection = localDirection;
  let asymptoticallyCorrected = false;
  if (reachedBackground && config.asymptoticDirection) {
    const positionWorld = frame
      ? liftPosition(model, frame, result.final)
      : model.geometry.toCartesianPosition(result.final.position_x);
    exitDirection = config.asymptoticDirection(positionWorld, localDirection);
    asymptoticallyCorrected = true;
  }

  const base = {
    steps: result.steps,
    nullResidual,
    terminationReason: result.reason,
    final: result.final,
    exitDirection,
    asymptoticallyCorrected,
    ...(frame ? { frame } : {}),
  };

  if (captured) {
    return { ...base, outcome: 'captured', color: SHADOW_COLOR, radiance: LINEAR_BLACK };
  }

  let outcome: RayOutcome;
  switch (result.reason) {
    case 'event':
      outcome = 'background';
      break;
    case 'terminator':
      // Only the capture test is a terminator here, and it is handled above.
      outcome = 'numerical-failure';
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
    return { ...base, outcome, color: FAILED_RAY_COLOR, radiance: decodeSrgb(FAILED_RAY_COLOR) };
  }

  const color = sampleCelestialGrid(grid, exitDirection[0], exitDirection[1], exitDirection[2]);
  return { ...base, outcome, color, radiance: decodeSrgb(color) };
}

export interface RenderDiagnostics {
  /** Total geodesics integrated: pixels times samples per pixel. */
  readonly raysTraced: number;
  readonly samplesPerPixel: number;
  readonly raysReachingBackground: number;
  readonly raysFailed: number;
  /** Rays that ended on the horizon: the black-hole shadow. */
  readonly raysCaptured: number;
  /** Pixels containing at least one failed sample; drawn in the failure colour. */
  readonly pixelsWithFailures: number;
  readonly totalSteps: number;
  /** Largest |g_mu_nu k^mu k^nu| over every traced ray. */
  readonly maxNullResidual: number;
  readonly health: NumericalHealth;
  /** The tolerance the residual was judged against. */
  readonly residualTolerance: Tolerance;
  /** Whether background lookups used directions at infinity. */
  readonly asymptoticallyCorrected: boolean;
}

export interface RenderResult {
  readonly widthPx: number;
  readonly heightPx: number;
  /**
   * RGBA8 sRGB, row-major from the top-left pixel.
   *
   * Pinned to a plain ArrayBuffer rather than the default ArrayBufferLike so the buffer
   * can be handed straight to an ImageData, which does not accept a SharedArrayBuffer.
   */
  readonly pixels: Uint8ClampedArray<ArrayBuffer>;
  /** Mean linear-light RGB per pixel, three entries per pixel, before encoding. */
  readonly linear: Float64Array;
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
  readonly samplesPerAxis: number;
  readonly seed: number;
  readonly pixels: Uint8ClampedArray<ArrayBuffer>;
  readonly linear: Float64Array;
  readonly session: FormulationSession;
  readonly residualTolerance: Tolerance;
  raysReachingBackground: number;
  raysFailed: number;
  raysCaptured: number;
  pixelsWithFailures: number;
  totalSteps: number;
  maxNullResidual: number;
}

export function createAccumulator(config: TraceConfig, screen: PinholeScreen): RenderAccumulator {
  const samplesPerAxis = Math.max(1, Math.floor(config.sampling?.samplesPerAxis ?? 1));
  return {
    widthPx: screen.widthPx,
    heightPx: screen.heightPx,
    samplesPerAxis,
    seed: config.sampling?.seed ?? 0,
    pixels: new Uint8ClampedArray(screen.widthPx * screen.heightPx * 4),
    linear: new Float64Array(screen.widthPx * screen.heightPx * 3),
    session: (config.formulation ?? HAMILTONIAN).bind(config.model),
    residualTolerance: config.residualTolerance ?? NULL_NORMALIZATION_TRACED,
    raysReachingBackground: 0,
    raysFailed: 0,
    raysCaptured: 0,
    pixelsWithFailures: 0,
    totalSteps: 0,
    maxNullResidual: 0,
  };
}

/**
 * Trace the rows [rowStart, rowEnd) of an image into the accumulator.
 *
 * Each pixel's samples are all traced before moving on, so every pixel in a finished
 * band is final and can be displayed immediately.
 */
export function renderBand(
  config: TraceConfig,
  screen: PinholeScreen,
  accumulator: RenderAccumulator,
  rowStart: number,
  rowEnd: number,
): void {
  const n = accumulator.samplesPerAxis;
  const samples = n * n;
  const end = Math.min(rowEnd, screen.heightPx);

  for (let j = Math.max(0, rowStart); j < end; j += 1) {
    for (let i = 0; i < screen.widthPx; i += 1) {
      let r = 0;
      let g = 0;
      let b = 0;
      let failed = false;

      for (let sy = 0; sy < n; sy += 1) {
        for (let sx = 0; sx < n; sx += 1) {
          // One ray through the pixel centre when n = 1, so a single-sample render is
          // unchanged; otherwise one jittered point per stratum.
          const index = sy * n + sx;
          const jx = n === 1 ? 0.5 : hashToUnit(i, j, 2 * index, accumulator.seed);
          const jy = n === 1 ? 0.5 : hashToUnit(i, j, 2 * index + 1, accumulator.seed);
          const direction = localRayDirectionAt(screen, i + (sx + jx) / n, j + (sy + jy) / n);
          const ray = traceRay(config, generateNullRay(config.observer, direction), accumulator.session);

          accumulator.totalSteps += ray.steps;
          if (ray.outcome === 'background') accumulator.raysReachingBackground += 1;
          if (ray.outcome === 'captured') accumulator.raysCaptured += 1;
          if (ray.outcome !== 'background' && ray.outcome !== 'captured') {
            accumulator.raysFailed += 1;
            failed = true;
          }

          const residual = Math.abs(ray.nullResidual);
          if (Number.isFinite(residual) && residual > accumulator.maxNullResidual) {
            accumulator.maxNullResidual = residual;
          }

          r += ray.radiance.r;
          g += ray.radiance.g;
          b += ray.radiance.b;
        }
      }

      const pixel = j * screen.widthPx + i;
      accumulator.linear[pixel * 3] = r / samples;
      accumulator.linear[pixel * 3 + 1] = g / samples;
      accumulator.linear[pixel * 3 + 2] = b / samples;

      const offset = pixel * 4;
      if (failed) {
        // A pixel with any failed sample is drawn entirely in the failure colour rather
        // than averaged: blending would dilute a failure into a plausible tint, which is
        // precisely what CLAUDE.md §17 and §24 forbid.
        accumulator.pixelsWithFailures += 1;
        accumulator.pixels[offset] = FAILED_RAY_COLOR.r;
        accumulator.pixels[offset + 1] = FAILED_RAY_COLOR.g;
        accumulator.pixels[offset + 2] = FAILED_RAY_COLOR.b;
      } else {
        accumulator.pixels[offset] = linearToSrgb8(r / samples);
        accumulator.pixels[offset + 1] = linearToSrgb8(g / samples);
        accumulator.pixels[offset + 2] = linearToSrgb8(b / samples);
      }
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
export function finalizeRender(accumulator: RenderAccumulator, config?: TraceConfig): RenderResult {
  const nullCheck = checkTolerance(accumulator.residualTolerance, accumulator.maxNullResidual);
  const samplesPerPixel = accumulator.samplesPerAxis * accumulator.samplesPerAxis;
  return {
    widthPx: accumulator.widthPx,
    heightPx: accumulator.heightPx,
    pixels: accumulator.pixels,
    linear: accumulator.linear,
    diagnostics: {
      raysTraced: accumulator.widthPx * accumulator.heightPx * samplesPerPixel,
      samplesPerPixel,
      raysReachingBackground: accumulator.raysReachingBackground,
      raysFailed: accumulator.raysFailed,
      raysCaptured: accumulator.raysCaptured,
      pixelsWithFailures: accumulator.pixelsWithFailures,
      totalSteps: accumulator.totalSteps,
      maxNullResidual: accumulator.maxNullResidual,
      health: summarizeHealth([nullCheck], accumulator.raysFailed > 0),
      residualTolerance: accumulator.residualTolerance,
      asymptoticallyCorrected: config?.asymptoticDirection !== undefined,
    },
  };
}

/** Render a whole image in one pass. */
export function renderImage(config: TraceConfig, screen: PinholeScreen): RenderResult {
  const accumulator = createAccumulator(config, screen);
  renderBand(config, screen, accumulator, 0, screen.heightPx);
  return finalizeRender(accumulator, config);
}
