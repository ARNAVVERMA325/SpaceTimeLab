import type { Vec4 } from '../core/indices.js';
import { nullState, type PhaseSpaceState } from '../core/phase-space.js';
import type { SpacetimeModel } from '../spacetimes/spacetime-model.js';
import { boostTetrad, frameToCoordinate, Tetrad } from './tetrad.js';

/**
 * Observer layer (CLAUDE.md §4, §20).
 *
 * An observer is an event x^mu, a timelike four-velocity u^mu, and a local orthonormal
 * spatial frame — here carried together as a tetrad whose leg (0) is u^mu.
 *
 * CLAUDE.md §4 draws a line this module respects: a physical observer is what
 * observer-based rendering must use, while a purely mathematical geometry-inspection
 * view may use an abstract visualization camera provided it is labelled as a
 * visualization coordinate system. `kind` records which one an instance is, so the UI
 * can label it correctly rather than implying every camera is a physical observer.
 */
export type ObserverKind = 'physical-observer' | 'visualization-camera';

export interface Observer {
  readonly id: string;
  readonly displayName: string;
  readonly kind: ObserverKind;
  readonly position_x: Vec4;
  readonly tetrad: Tetrad;
  /** What this observer is, in words, for the provenance panel (CLAUDE.md §22). */
  readonly description: string;
}

/**
 * A static observer at rest in Minkowski Cartesian coordinates.
 *
 * In this chart the metric components are already eta_ab everywhere, so the coordinate
 * basis is orthonormal and the tetrad is the identity frame: u^mu = (1,0,0,0) and the
 * spatial legs are the coordinate axes. The general static-observer construction,
 * u^mu = (1/sqrt(-g_00), 0, 0, 0) with a Gram-Schmidt spatial triad, is a Milestone 3
 * deliverable (ROADMAP.md 3.1) and is deliberately not faked here.
 */
export function staticMinkowskiObserver(position_x: Vec4): Observer {
  return {
    id: 'static-minkowski',
    displayName: 'Static observer at rest in the global inertial frame',
    kind: 'physical-observer',
    position_x,
    tetrad: Tetrad.identity(),
    description:
      'At rest in Minkowski Cartesian coordinates, u^mu = (1,0,0,0). The coordinate ' +
      'basis is already orthonormal in this chart, so the tetrad is the identity frame.',
  };
}

/**
 * A static observer hovering at fixed spatial coordinates in a diagonal static metric.
 *
 * Its four-velocity is u^mu = (1 / sqrt(-g_00), 0, 0, 0), and its spatial frame is the
 * orthonormalized coordinate basis. For Schwarzschild this is the shell observer: to
 * stay at fixed r outside a black hole it must accelerate, which is why its frame is
 * not the frame of a freely-falling observer and why Milestone 3 adds the latter rather
 * than reusing this one.
 *
 * The frame's legs are ordered as the chart's coordinates: for Schwarzschild, leg (1)
 * points along increasing r, leg (2) along increasing theta, leg (3) along increasing phi.
 */
export function staticObserver(model: SpacetimeModel, position_x: Vec4, id = 'static'): Observer {
  const domain = model.domainCheck(position_x);
  if (!domain.inDomain) {
    throw new RangeError(
      `staticObserver: the requested event is outside the chart domain. ${domain.reason}`,
    );
  }

  const tetrad = Tetrad.diagonalStatic(model.metricAt(position_x));
  const names = model.chart.coordinateNames;

  return {
    id,
    displayName: `Static observer at fixed ${names[1]} = ${position_x[1]}`,
    kind: 'physical-observer',
    position_x,
    tetrad,
    description:
      `At rest in the ${model.chart.displayName} chart, with four-velocity ` +
      'u^mu = (1/sqrt(-g_00), 0, 0, 0). Its spatial frame is the orthonormalized ' +
      `coordinate basis: leg (1) along ${names[1]}, leg (2) along ${names[2]}, leg (3) ` +
      `along ${names[3]}. Holding this position requires proper acceleration; this is ` +
      'not a freely-falling frame.',
  };
}

/**
 * An observer falling radially from rest at infinity (ROADMAP.md 3.1).
 *
 * Its four-velocity in Schwarzschild coordinates is u^mu = (1/f, -sqrt(2M/r), 0, 0): energy
 * E = 1, no angular momentum. Relative to the static observer at the same event it moves
 * radially inward at v = sqrt(2M/r), so its frame is the static frame boosted by
 * beta = -v along leg (1). At r = 20M that is 0.32c; it approaches c as r -> 2M, where
 * the static frame itself ceases to exist.
 *
 * This is a geodesic observer — it feels no acceleration — in contrast to the static one,
 * which must accelerate to hover. The two see different images from the same event:
 * aberration concentrates the view ahead of the motion, so the infalling observer sees the
 * shadow smaller, and the Doppler shift brightens and blueshifts what lies ahead.
 */
export function freeFallingObserver(
  model: SpacetimeModel,
  position_x: Vec4,
  mass: number,
  id = 'free-fall',
): Observer {
  const hover = staticObserver(model, position_x, id);
  const v = Math.sqrt((2 * mass) / position_x[1]);
  const tetrad = boostTetrad(hover.tetrad, [-v, 0, 0]);
  return {
    id,
    displayName: `Observer falling radially from rest at infinity, now at r = ${position_x[1]}`,
    kind: 'physical-observer',
    position_x,
    tetrad,
    description:
      'Freely falling, with energy E = 1 and no angular momentum: u^mu = (1/f, -sqrt(2M/r), 0, 0). ' +
      `Moving inward at ${v.toFixed(4)}c relative to a static observer at the same event, so ` +
      'its frame is the static frame Lorentz-boosted along the radial leg. A geodesic ' +
      'observer: it feels no acceleration.',
  };
}

/**
 * A pinhole screen in the observer's local rest frame (ROADMAP.md 3.2 groundwork).
 *
 * `forward`, `right` and `up` are unit 3-vectors in the observer's local orthonormal
 * spatial frame — components against tetrad legs (1),(2),(3), not coordinate components.
 */
export interface PinholeScreen {
  readonly widthPx: number;
  readonly heightPx: number;
  /** Horizontal field of view, radians. */
  readonly horizontalFovRad: number;
  readonly forward: readonly [number, number, number];
  readonly right: readonly [number, number, number];
  readonly up: readonly [number, number, number];
}

export function defaultScreen(widthPx: number, heightPx: number, horizontalFovRad: number): PinholeScreen {
  return {
    widthPx,
    heightPx,
    horizontalFovRad,
    forward: [1, 0, 0],
    right: [0, 1, 0],
    up: [0, 0, 1],
  };
}

/**
 * A screen for an observer in a spherical chart, looking inward at the centre.
 *
 * Frame legs are (1) = radial, (2) = polar, (3) = azimuthal, so looking inward is
 * forward = -leg(1). Up is taken as -leg(2), which points toward decreasing theta, i.e.
 * north. Right is then -leg(3) so that right x up = forward, matching the handedness
 * `defaultScreen` uses and keeping the rendered image unmirrored.
 */
export function inwardFacingScreen(
  widthPx: number,
  heightPx: number,
  horizontalFovRad: number,
): PinholeScreen {
  return {
    widthPx,
    heightPx,
    horizontalFovRad,
    forward: [-1, 0, 0],
    right: [0, 0, -1],
    up: [0, -1, 0],
  };
}

/**
 * The unit viewing direction for a pixel, in the observer's local spatial frame.
 *
 * Pixel centres are sampled at (i + 0.5, j + 0.5). Row 0 is the top of the image, so
 * the vertical normalized coordinate runs from +1 down to -1.
 */
export function localRayDirection(
  screen: PinholeScreen,
  i: number,
  j: number,
): [number, number, number] {
  return localRayDirectionAt(screen, i + 0.5, j + 0.5);
}

/**
 * The unit viewing direction through a continuous screen position (u, v).
 *
 * u runs over [0, widthPx) left to right and v over [0, heightPx) top to bottom, so the
 * centre of pixel (i, j) is (i + 0.5, j + 0.5). Subpixel sampling uses this directly.
 */
export function localRayDirectionAt(
  screen: PinholeScreen,
  u: number,
  v: number,
): [number, number, number] {
  const { widthPx, heightPx } = screen;
  const i = u - 0.5;
  const j = v - 0.5;
  const tanHalfFov = Math.tan(screen.horizontalFovRad / 2);
  const aspect = widthPx / heightPx;

  const xNdc = (u / widthPx) * 2 - 1;
  const yNdc = 1 - (v / heightPx) * 2;

  const sx = xNdc * tanHalfFov;
  const sy = (yNdc * tanHalfFov) / aspect;

  const d: [number, number, number] = [
    screen.forward[0] + sx * screen.right[0] + sy * screen.up[0],
    screen.forward[1] + sx * screen.right[1] + sy * screen.up[1],
    screen.forward[2] + sx * screen.right[2] + sy * screen.up[2],
  ];

  const norm = Math.hypot(d[0], d[1], d[2]);
  if (!(norm > 0) || !Number.isFinite(norm)) {
    throw new RangeError(`localRayDirection: degenerate screen basis at pixel (${i}, ${j}).`);
  }
  return [d[0] / norm, d[1] / norm, d[2] / norm];
}

/**
 * Turn a local viewing direction into a null geodesic initial condition.
 *
 * Steps 1 and 2 of the rendering pipeline in CLAUDE.md §9. In the observer's
 * orthonormal frame the wavevector is
 *
 *   k^(a) = (-1, n^1, n^2, n^3),   |n| = 1
 *
 * which is null, since eta_ab k^(a) k^(b) = -(-1)^2 + |n|^2 = 0.
 *
 * The sign of k^(0) is deliberate. Backward ray tracing follows the photon from the
 * observer back towards its source, so the wavevector is *past*-directed: as the affine
 * parameter increases the coordinate time decreases and the ray moves outward along n.
 * The traced curve is the photon's own worldline — time reversal maps null geodesics to
 * null geodesics — so the geometry is unaffected by this choice, but the sense matters
 * for the emitter-to-observer frequency shift that Milestone 3 introduces, and encoding
 * it now avoids a sign error there.
 */
export function generateNullRay(
  observer: Observer,
  direction: readonly [number, number, number],
): PhaseSpaceState {
  const frameComponents: Vec4 = [-1, direction[0], direction[1], direction[2]];
  const null_wavevector_k = frameToCoordinate(observer.tetrad, frameComponents);
  return nullState(observer.position_x, null_wavevector_k, 0);
}
