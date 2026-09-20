/**
 * Visualization layer (CLAUDE.md §20): the background the backward-traced rays land on.
 *
 * A celestial sphere of coordinate radius `radius`, carrying a latitude/longitude grid.
 * ROADMAP.md 1.5 asks for exactly this and calls it deliberately undramatic: in flat
 * spacetime the rendered grid must come back undistorted, which is what proves the
 * pipeline rather than the picture.
 *
 * This is a visualization mapping, not physics. The grid is a texture on a sphere at a
 * chosen coordinate radius, not a physical emitting surface, and it carries no emission
 * model or redshift — those arrive with Milestone 3.
 */

export interface RGB {
  readonly r: number;
  readonly g: number;
  readonly b: number;
}

export interface CelestialGrid {
  /** Coordinate radius at which rays are considered to have reached the background. */
  readonly radius: number;
  /** Number of latitude bands. */
  readonly latitudeDivisions: number;
  /** Number of longitude sectors. */
  readonly longitudeDivisions: number;
  /** Half-width of a grid line, in radians of the relevant angle. */
  readonly lineHalfWidthRad: number;
  readonly lineColor: RGB;
  readonly backgroundColorA: RGB;
  readonly backgroundColorB: RGB;
}

export const DEFAULT_CELESTIAL_GRID: CelestialGrid = Object.freeze({
  radius: 100,
  latitudeDivisions: 12,
  longitudeDivisions: 24,
  lineHalfWidthRad: 0.012,
  lineColor: { r: 235, g: 240, b: 248 },
  backgroundColorA: { r: 14, g: 18, b: 30 },
  backgroundColorB: { r: 26, g: 34, b: 54 },
});

/** Angular position on the celestial sphere for a direction from the origin. */
export interface SkyAngles {
  /** Polar angle from the +z axis, in [0, pi]. */
  readonly theta: number;
  /** Azimuth in the x-y plane, in [0, 2pi). */
  readonly phi: number;
}

export function directionToSkyAngles(x: number, y: number, z: number): SkyAngles {
  const r = Math.hypot(x, y, z);
  if (!(r > 0) || !Number.isFinite(r)) {
    throw new RangeError('directionToSkyAngles: direction has zero or non-finite length.');
  }
  const theta = Math.acos(Math.min(1, Math.max(-1, z / r)));
  let phi = Math.atan2(y, x);
  if (phi < 0) phi += 2 * Math.PI;
  return { theta, phi };
}

/**
 * Distance from `angle` to the nearest multiple of `spacing`, in the same units.
 */
function distanceToNearestGridLine(angle: number, spacing: number): number {
  const offset = angle / spacing;
  return Math.abs(offset - Math.round(offset)) * spacing;
}

/**
 * Sample the grid for a direction, returning a colour.
 *
 * Longitude lines are drawn with their angular half-width scaled by sin(theta), so that
 * meridians keep a roughly constant apparent thickness instead of fanning out near the
 * poles.
 */
export function sampleCelestialGrid(grid: CelestialGrid, x: number, y: number, z: number): RGB {
  const { theta, phi } = directionToSkyAngles(x, y, z);

  const latitudeSpacing = Math.PI / grid.latitudeDivisions;
  const longitudeSpacing = (2 * Math.PI) / grid.longitudeDivisions;

  const latitudeDistance = distanceToNearestGridLine(theta, latitudeSpacing);
  const sinTheta = Math.max(Math.sin(theta), 1e-6);
  const longitudeDistance = distanceToNearestGridLine(phi, longitudeSpacing) * sinTheta;

  if (latitudeDistance < grid.lineHalfWidthRad || longitudeDistance < grid.lineHalfWidthRad) {
    return grid.lineColor;
  }

  // Checker the cells so the grid stays readable where lines are sparse.
  const latitudeCell = Math.floor(theta / latitudeSpacing);
  const longitudeCell = Math.floor(phi / longitudeSpacing);
  return (latitudeCell + longitudeCell) % 2 === 0 ? grid.backgroundColorA : grid.backgroundColorB;
}
