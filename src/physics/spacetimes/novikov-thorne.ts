import { PHYSICAL_CONSTANTS, SOLAR_MASS_KG } from '../radiation/constants.js';
import { circularOrbit } from './schwarzschild-orbits.js';

/**
 * A geometrically thin, optically thick accretion disk around a Schwarzschild black hole:
 * the Novikov-Thorne model with the Page-Thorne radiative flux (ROADMAP.md 3.4).
 *
 * Model assumptions, stated because CLAUDE.md §14 and §22 require the omitted physics to
 * be disclosed:
 * - the disk lies in the equatorial plane and has negligible thickness;
 * - gas moves on prograde Keplerian circular geodesics, from the ISCO (6M) outward;
 * - there is no torque at the inner edge, so the flux vanishes there;
 * - each face radiates locally as a blackbody at the temperature set by sigma T^4 = F,
 *   isotropically in the gas rest frame — no limb darkening, no electron-scattering
 *   spectral hardening (colour correction), no atmosphere;
 * - no returning radiation, self-irradiation, or radiation from inside the ISCO;
 * - the disk is stationary and does not affect the spacetime.
 *
 * The flux, derived from Page & Thorne, ApJ 191, 499 (1974), with e^{nu+psi+mu} = r in the
 * Schwarzschild equatorial plane, is F = Mdot f(r) / (4 pi r) with
 *
 *   f(r) = -Omega_{,r} / (E - Omega L)^2  Integral_{r_isco}^{r} (E - Omega L) L_{,r} dr
 *
 * which for Schwarzschild integrates in closed form. With x = sqrt(r/M):
 *
 *   f(r) = 3 N(x) / (4 M x^3 (3 - x^2))
 *   N(x) = -2x + 2 sqrt 6 + sqrt 3 [ ln((x - sqrt 3)/(x + sqrt 3)) + 2 ln(1 + sqrt 2) ]
 *
 * Derived with sympy and verified in mpmath by the model's own energy balance: the
 * luminosity reaching infinity, Integral 4 pi r E F dr, equals Mdot (1 - E_isco) to all
 * 20 digits computed. It vanishes at the ISCO, peaks at r = 9.550928 M, and far out goes
 * as (3 M Mdot / 8 pi r^3)(1 - C sqrt(M/r)) with C = sqrt 6 + sqrt 3 ln(1 + sqrt 2).
 */

const SQRT2 = Math.SQRT2;
const SQRT3 = Math.sqrt(3);
const SQRT6 = Math.sqrt(6);
const TWO_LN_1_PLUS_SQRT2 = 2 * Math.log(1 + SQRT2);

/** Radius of peak flux, from mpmath root-finding on dF/dr, in units of M. */
export const NOVIKOV_THORNE_PEAK_RADIUS = 9.55092807794387;

/** Radiative efficiency of the disk: the binding energy at the ISCO, 1 - sqrt(8/9). */
export const NOVIKOV_THORNE_EFFICIENCY = 1 - Math.sqrt(8 / 9);

/** Large-radius coefficient: F -> (3 M Mdot / 8 pi r^3)(1 - C sqrt(M/r)). */
export const NOVIKOV_THORNE_ASYMPTOTIC_C = SQRT6 + SQRT3 * Math.log(1 + SQRT2);

/**
 * Flux per unit accretion rate, F / Mdot, in geometric units (1 / length^2), from one
 * face of the disk. Zero at and inside the ISCO.
 */
export function novikovThorneFluxPerAccretionRate(M: number, r: number): number {
  if (!(r > 6 * M)) return 0;
  const x = Math.sqrt(r / M);
  const n =
    -2 * x + 2 * SQRT6 + SQRT3 * (Math.log((x - SQRT3) / (x + SQRT3)) + TWO_LN_1_PLUS_SQRT2);
  const f = (3 * n) / (4 * M * x * x * x * (3 - x * x));
  return Math.max(0, f / (4 * Math.PI * r));
}

/** Specific energy E(r) of the circular orbit at r, as used in the energy balance. */
export function diskOrbitEnergy(M: number, r: number): number {
  return circularOrbit(M, r).energy_E;
}

export interface ThinDiskPhysicalParameters {
  /** Black-hole mass in solar masses. */
  readonly massSolar: number;
  /** Accretion rate as a fraction of the Eddington rate at this disk's own efficiency. */
  readonly eddingtonFraction: number;
}

export interface ThinDiskScale {
  readonly massKg: number;
  /** GM / c^2, the length unit of the geometric calculation, in metres. */
  readonly gravitationalRadiusM: number;
  readonly eddingtonLuminosityW: number;
  readonly accretionRateKgPerS: number;
  readonly accretionRateSolarPerYear: number;
  /** Peak rest-frame temperature, at r = 9.55M. */
  readonly peakTemperatureK: number;
  /** Rest-frame temperature at radius r (in units of M). */
  temperatureK(rOverM: number): number;
}

/**
 * Put the model in physical units.
 *
 * The geometric calculation is scale-free: F / Mdot depends only on r / M. A physical
 * temperature needs the mass and the accretion rate, and follows from sigma T^4 = F with
 *
 *   F_phys = Mdot c^2 (F M^2 / Mdot)(r/M) / r_g^2,   r_g = G M / c^2.
 *
 * The accretion rate is given as a fraction of Eddington, L_Edd = 4 pi G M m_p c / sigma_T,
 * with Mdot_Edd = L_Edd / (eta c^2) and eta this disk's own efficiency. Nothing in this is
 * adjustable for appearance: a 1e9 solar-mass hole at a tenth of Eddington has a peak
 * temperature near 39,500 K, and the colour on screen is what that temperature gives.
 */
export function thinDiskScale(params: ThinDiskPhysicalParameters): ThinDiskScale {
  const { massSolar, eddingtonFraction } = params;
  if (!(massSolar > 0) || !(eddingtonFraction > 0)) {
    throw new RangeError('thinDiskScale: mass and Eddington fraction must be positive.');
  }
  const { G, c, protonMass, thomsonCrossSection, sigmaSB } = PHYSICAL_CONSTANTS;
  const massKg = massSolar * SOLAR_MASS_KG;
  const rg = (G * massKg) / (c * c);
  const eddingtonLuminosityW = (4 * Math.PI * G * massKg * protonMass * c) / thomsonCrossSection;
  const accretionRateKgPerS = (eddingtonFraction * eddingtonLuminosityW) / (NOVIKOV_THORNE_EFFICIENCY * c * c);

  const temperatureK = (rOverM: number): number => {
    const fluxShape = novikovThorneFluxPerAccretionRate(1, rOverM); // F M^2 / Mdot, with M = 1
    const flux = (accretionRateKgPerS * c * c * fluxShape) / (rg * rg);
    return Math.pow(flux / sigmaSB, 0.25);
  };

  return {
    massKg,
    gravitationalRadiusM: rg,
    eddingtonLuminosityW,
    accretionRateKgPerS,
    accretionRateSolarPerYear: (accretionRateKgPerS * 365.25 * 86400) / SOLAR_MASS_KG,
    peakTemperatureK: temperatureK(NOVIKOV_THORNE_PEAK_RADIUS),
    temperatureK,
  };
}
