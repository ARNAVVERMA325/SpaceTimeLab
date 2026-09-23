/**
 * Colour encoding for the visualization layer.
 *
 * Averaging must happen in linear light. The sRGB values in an image file are
 * gamma-encoded, so averaging them directly darkens every edge — a mixed pixel of black
 * and white comes out at sRGB 128, which is about 22% of white's luminance rather than
 * 50%. Supersampled pixels are therefore accumulated in linear RGB and encoded once.
 *
 * Transfer functions are the piecewise definitions of IEC 61966-2-1 (sRGB).
 */

/** Relative linear-light RGB, where 1 is display white. May exceed 1 before encoding. */
export interface LinearRGB {
  readonly r: number;
  readonly g: number;
  readonly b: number;
}

export const LINEAR_BLACK: LinearRGB = Object.freeze({ r: 0, g: 0, b: 0 });

/** Decode one 8-bit sRGB channel to linear light in [0, 1]. */
export function srgbToLinear(channel8: number): number {
  const c = channel8 / 255;
  return c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4);
}

/** Encode linear light to one 8-bit sRGB channel, clamping to the displayable range. */
export function linearToSrgb8(value: number): number {
  const v = Math.min(1, Math.max(0, value));
  const c = v <= 0.0031308 ? 12.92 * v : 1.055 * Math.pow(v, 1 / 2.4) - 0.055;
  return Math.round(c * 255);
}

export function decodeSrgb(color: { readonly r: number; readonly g: number; readonly b: number }): LinearRGB {
  return { r: srgbToLinear(color.r), g: srgbToLinear(color.g), b: srgbToLinear(color.b) };
}

/**
 * A deterministic hash to [0, 1), for reproducible sample jitter.
 *
 * "lowbias32" (C. Wellons), an integer hash with good avalanche behaviour. Determinism
 * matters here: two renders of the same scene must be bit-identical, so that a parallel
 * render can be checked against a serial one exactly and a regression shows up as a
 * changed pixel rather than as noise.
 */
export function hashToUnit(a: number, b: number, c: number, seed: number): number {
  let x = (a * 0x9e3779b1) ^ (b * 0x85ebca77) ^ (c * 0xc2b2ae3d) ^ (seed * 0x27d4eb2f);
  x ^= x >>> 16;
  x = Math.imul(x, 0x7feb352d);
  x ^= x >>> 15;
  x = Math.imul(x, 0x846ca68b);
  x ^= x >>> 16;
  return (x >>> 0) / 4294967296;
}
