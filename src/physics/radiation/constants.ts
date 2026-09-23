/**
 * Physical constants in SI units, with sources (CLAUDE.md §24: never invent constants).
 *
 * The geodesic engine works in geometric units, G = c = 1, where none of these appear.
 * They are needed only where the model meets the laboratory: turning a black hole's mass
 * and accretion rate into a disk temperature, and a spectrum into a colour.
 */
export const PHYSICAL_CONSTANTS = Object.freeze({
  /** Speed of light, m s^-1. Exact (SI definition). */
  c: 299_792_458,
  /** Planck constant, J s. Exact (SI 2019). */
  h: 6.626_070_15e-34,
  /** Boltzmann constant, J K^-1. Exact (SI 2019). */
  kB: 1.380_649e-23,
  /** Stefan-Boltzmann constant, W m^-2 K^-4. Exact, derived from h, kB, c (CODATA 2018). */
  sigmaSB: 5.670_374_419e-8,
  /** Newtonian constant of gravitation, m^3 kg^-1 s^-2 (CODATA 2018). */
  G: 6.674_30e-11,
  /** Nominal solar mass parameter GM_sun, m^3 s^-2 (IAU 2015 Resolution B3). */
  GMsun: 1.327_124_4e20,
  /** Proton mass, kg (CODATA 2018). */
  protonMass: 1.672_621_923_69e-27,
  /** Thomson cross-section, m^2 (CODATA 2018). */
  thomsonCrossSection: 6.652_458_732_1e-29,
  /** Maximum luminous efficacy of radiation, lm W^-1 (SI definition of the candela, at 540 THz). */
  Km: 683,
  sources:
    'CODATA 2018 recommended values (Tiesinga et al., Rev. Mod. Phys. 93, 025010, 2021); ' +
    'IAU 2015 Resolution B3 for GM_sun; SI Brochure, 9th ed. (2019), for defined constants.',
});

/** Solar mass in kg, from the IAU nominal GM_sun and CODATA G. */
export const SOLAR_MASS_KG = PHYSICAL_CONSTANTS.GMsun / PHYSICAL_CONSTANTS.G;
