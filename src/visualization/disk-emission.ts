import type { PhaseSpaceState } from '../physics/core/phase-space.js';
import { frequencyRatio, keplerianEmitter } from '../physics/observer/frequency-shift.js';
import { blackbodyXYZ, blackbodyXYZFast, xyzToLinearSrgb } from '../physics/radiation/blackbody.js';
import {
  NOVIKOV_THORNE_PEAK_RADIUS,
  type ThinDiskScale,
} from '../physics/spacetimes/novikov-thorne.js';
import type { CartesianVec3, SpacetimeModel } from '../physics/spacetimes/spacetime-model.js';
import type { LinearRGB } from './color.js';
import { liftDirection, liftPosition, type OrbitalPlaneFrame } from './orbital-plane.js';

/**
 * Emission from a Novikov-Thorne disk, and the mapping from radiance to display.
 *
 * The physics ends at `diskSample`, which returns absolute radiance: the CIE XYZ of a
 * blackbody at the observed temperature g T(r), in cd m^-2, converted to linear sRGB.
 * Everything after that — exposure and tone curve — is a display mapping, chosen to fit a
 * range of radiances that spans orders of magnitude into what a screen can show. It is
 * applied once, at the end, and disclosed in the UI (CLAUDE.md §22). It changes how the
 * picture looks, never what was computed.
 */

export interface ThinDiskEmitter {
  /** Mass parameter in geometric units, matching the spacetime model. */
  readonly mass: number;
  readonly innerRadius: number;
  readonly outerRadius: number;
  /** Physical temperature law from mass and accretion rate. */
  readonly scale: ThinDiskScale;
}

export interface DiskSample {
  readonly radius: number;
  /** g = nu_obs / nu_emit. */
  readonly frequencyRatio: number;
  readonly emittedTemperatureK: number;
  readonly observedTemperatureK: number;
  /** Absolute linear sRGB, in units where 1 is 1 cd m^-2 of luminance. */
  readonly radiance: LinearRGB;
}

/**
 * The radiance a ray carries back from the point where it meets the disk.
 *
 * The photon's static-frame energy is sqrt(f) k^t — unchanged by the orbital-plane
 * rotation, since k^t is — and its static-frame direction is the lifted world direction.
 * The emitter is the prograde Keplerian gas at that point. Their contraction gives g, and
 * by the invariance of I_nu / nu^3 the observed spectrum is a blackbody at g T(r).
 */
export function diskSample(
  model: SpacetimeModel,
  disk: ThinDiskEmitter,
  state: PhaseSpaceState,
  frame: OrbitalPlaneFrame | undefined,
): DiskSample {
  const position: CartesianVec3 = frame
    ? liftPosition(model, frame, state)
    : model.geometry.toCartesianPosition(state.position_x);
  const direction: CartesianVec3 = frame
    ? liftDirection(model, frame, state)
    : model.geometry.toCartesianDirection(state.position_x, state.tangent);

  const r = Math.hypot(position[0], position[1], position[2]);
  const f = 1 - (2 * disk.mass) / r;
  const photon = { energy: Math.sqrt(f) * state.tangent[0], momentum: direction };
  const g = frequencyRatio(photon, keplerianEmitter(disk.mass, position));

  const emittedTemperatureK = disk.scale.temperatureK(r / disk.mass);
  const observedTemperatureK = g * emittedTemperatureK;
  const rgb =
    observedTemperatureK > 0
      ? xyzToLinearSrgb(blackbodyXYZFast(observedTemperatureK))
      : { r: 0, g: 0, b: 0 };

  return { radius: r, frequencyRatio: g, emittedTemperatureK, observedTemperatureK, radiance: rgb };
}

/** The luminance of the disk's hottest ring, seen at rest: the natural exposure reference. */
export function diskReferenceLuminance(disk: ThinDiskEmitter): number {
  return blackbodyXYZ(disk.scale.temperatureK(NOVIKOV_THORNE_PEAK_RADIUS)).Y;
}

/**
 * How linear radiance becomes a displayed colour.
 *
 * Display value = radiance * 2^exposureStops / referenceLuminance, so at 0 stops the
 * disk's hottest ring, seen at rest, sits at display white. Then an optional tone curve,
 * then sRGB encoding.
 */
export interface DisplayMapping {
  readonly referenceLuminance: number;
  readonly exposureStops: number;
  readonly toneMap: 'none' | 'reinhard';
}

export function displayGain(mapping: DisplayMapping): number {
  return Math.pow(2, mapping.exposureStops) / mapping.referenceLuminance;
}

/**
 * Reinhard's global operator on luminance, L' = L / (1 + L), preserving chromaticity.
 *
 * Reinhard et al., "Photographic tone reproduction for digital images", SIGGRAPH 2002.
 * Luminance-based so that it compresses brightness without shifting hue: a blueshifted
 * ring stays blue as it saturates, instead of washing out channel by channel. Negative,
 * out-of-gamut components are clamped to zero first, which is the one place this mapping
 * alters colour rather than brightness.
 */
export function toneMapReinhard(color: LinearRGB): LinearRGB {
  const r = Math.max(0, color.r);
  const g = Math.max(0, color.g);
  const b = Math.max(0, color.b);
  const luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b;
  if (!(luminance > 0)) return { r: 0, g: 0, b: 0 };
  const scale = 1 / (1 + luminance);
  return { r: r * scale, g: g * scale, b: b * scale };
}
