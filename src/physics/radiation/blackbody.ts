import { BLACKBODY_CHROMATICITY_REFERENCE, CIE1931_CMF_5NM } from '../../data/cie1931-cmf.js';
import { PHYSICAL_CONSTANTS } from './constants.js';

/**
 * Blackbody radiation and its colour.
 *
 * Planck's law, integrated against the CIE 1931 2-degree colour-matching functions. The
 * result is absolute: CIE XYZ scaled by K_m = 683 lm W^-1, so that Y is the photometric
 * luminance in cd m^-2 of a surface radiating as a blackbody at temperature T. Colour and
 * brightness therefore come out of the same physics rather than being chosen separately.
 *
 * The Lorentz invariance of I_nu / nu^3 makes this exact under frequency shift: a
 * blackbody at T_emit, seen with frequency ratio g = nu_obs / nu_emit, has observed
 * specific intensity g^3 B_{nu/g}(T_emit) = B_nu(g T_emit). The observed spectrum is again
 * Planckian, at g T_emit. So Doppler beaming, gravitational shift and the change of colour
 * all follow from evaluating this function at the shifted temperature, with no separate
 * beaming factor to get wrong.
 */

const { h, c, kB, Km } = PHYSICAL_CONSTANTS;

/** Planck spectral radiance B_lambda(T) in W m^-2 sr^-1 m^-1, at wavelength lambda in metres. */
export function planckSpectralRadiance(wavelengthM: number, temperatureK: number): number {
  if (!(temperatureK > 0)) return 0;
  const x = (h * c) / (wavelengthM * kB * temperatureK);
  // expm1 keeps precision in the Rayleigh-Jeans limit, where x is small and exp(x) - 1
  // would cancel.
  if (x > 700) return 0;
  return (2 * h * c * c) / (Math.pow(wavelengthM, 5) * Math.expm1(x));
}

export interface XYZ {
  readonly X: number;
  readonly Y: number;
  readonly Z: number;
}

/**
 * Absolute CIE XYZ of a blackbody: X = K_m integral B_lambda(T) x-bar(lambda) d lambda, and so on.
 *
 * Integrated with Simpson's rule over the 5 nm CIE table, 360-830 nm. Y is luminance in
 * cd m^-2. `tests/physics/blackbody.test.ts` checks the chromaticities against
 * colour-science's own 1 nm integration.
 */
export function blackbodyXYZ(temperatureK: number): XYZ {
  const table = CIE1931_CMF_5NM;
  const n = table.length; // 95 points, 94 intervals: Simpson needs an even count.
  const step = 5e-9;
  let X = 0;
  let Y = 0;
  let Z = 0;
  for (let i = 0; i < n; i += 1) {
    const [nm, xb, yb, zb] = table[i];
    const weight = i === 0 || i === n - 1 ? 1 : i % 2 === 1 ? 4 : 2;
    const B = planckSpectralRadiance(nm * 1e-9, temperatureK);
    X += weight * B * xb;
    Y += weight * B * yb;
    Z += weight * B * zb;
  }
  const scale = (Km * step) / 3;
  return { X: X * scale, Y: Y * scale, Z: Z * scale };
}

export function chromaticity(xyz: XYZ): { readonly x: number; readonly y: number } {
  const sum = xyz.X + xyz.Y + xyz.Z;
  return { x: xyz.X / sum, y: xyz.Y / sum };
}

/**
 * CIE XYZ to linear sRGB, IEC 61966-2-1:1999 (D65 white).
 *
 * Out-of-gamut colours give negative components. They are clamped to zero by the display
 * mapping, not here: this function is colorimetry, and clamping is a display decision.
 */
export function xyzToLinearSrgb(xyz: XYZ): { readonly r: number; readonly g: number; readonly b: number } {
  return {
    r: 3.2406 * xyz.X - 1.5372 * xyz.Y - 0.4986 * xyz.Z,
    g: -0.9689 * xyz.X + 1.8758 * xyz.Y + 0.0415 * xyz.Z,
    b: 0.0557 * xyz.X - 0.204 * xyz.Y + 1.057 * xyz.Z,
  };
}

/**
 * Memoized blackbody XYZ on a logarithmic temperature grid, with quadratic interpolation
 * in log XYZ against log T.
 *
 * A full image evaluates the blackbody once per disk-hitting ray; 95 exponentials each is
 * affordable but wasteful. The table spans 300 K to 10^8 K at 256 points per decade.
 * Quadratic rather than linear interpolation, for a measured reason: at low temperature
 * the visible band sits deep in the Wien tail, where log XYZ is strongly curved in log T,
 * and linear interpolation at 64 points per decade was off by up to 5.9e-3 there. Low
 * temperatures do occur — a Novikov-Thorne disk's flux falls to zero at its inner edge —
 * so the table has to be good there too. Quadratic at 256 per decade, about 1400 nodes:
 * worst case under 3e-6 in any channel, the Z (blue) channel near 400 K being hardest.
 */
const TABLE_LOG_MIN = Math.log(300);
const TABLE_LOG_MAX = Math.log(1e8);
const TABLE_POINTS_PER_DECADE = 256;
const TABLE_SIZE = Math.ceil(((TABLE_LOG_MAX - TABLE_LOG_MIN) / Math.LN10) * TABLE_POINTS_PER_DECADE) + 1;
const TABLE_STEP = (TABLE_LOG_MAX - TABLE_LOG_MIN) / (TABLE_SIZE - 1);
let table: Float64Array | undefined;

function buildTable(): Float64Array {
  const t = new Float64Array(TABLE_SIZE * 3);
  for (let i = 0; i < TABLE_SIZE; i += 1) {
    const xyz = blackbodyXYZ(Math.exp(TABLE_LOG_MIN + i * TABLE_STEP));
    t[i * 3] = Math.log(xyz.X);
    t[i * 3 + 1] = Math.log(xyz.Y);
    t[i * 3 + 2] = Math.log(xyz.Z);
  }
  return t;
}

export function blackbodyXYZFast(temperatureK: number): XYZ {
  const logT = Math.log(temperatureK);
  if (!(logT >= TABLE_LOG_MIN && logT <= TABLE_LOG_MAX)) return blackbodyXYZ(temperatureK);
  table ??= buildTable();
  const position = (logT - TABLE_LOG_MIN) / TABLE_STEP;
  // Three-point (quadratic) interpolation centred on the nearest node.
  const i = Math.min(TABLE_SIZE - 2, Math.max(1, Math.round(position)));
  const u = position - i;
  const t = table;
  const at = (k: number): number => {
    const ym = t[(i - 1) * 3 + k];
    const y0 = t[i * 3 + k];
    const yp = t[(i + 1) * 3 + k];
    return Math.exp(y0 + 0.5 * u * (yp - ym) + 0.5 * u * u * (yp - 2 * y0 + ym));
  };
  return { X: at(0), Y: at(1), Z: at(2) };
}

export { BLACKBODY_CHROMATICITY_REFERENCE };
