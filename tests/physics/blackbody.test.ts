import { describe, expect, it } from 'vitest';
import { CIE1931_CMF_5NM } from '../../src/data/cie1931-cmf.js';
import {
  BLACKBODY_CHROMATICITY_REFERENCE,
  blackbodyXYZ,
  blackbodyXYZFast,
  chromaticity,
  planckSpectralRadiance,
  xyzToLinearSrgb,
} from '../../src/physics/radiation/blackbody.js';
import { PHYSICAL_CONSTANTS } from '../../src/physics/radiation/constants.js';

/**
 * Blackbody radiation and colorimetry, checked against physics (the Stefan-Boltzmann and
 * Wien laws, which the implementation does not use) and against colour-science's
 * independent 1 nm integration.
 */

describe('Planck spectral radiance', () => {
  it('integrates to the Stefan-Boltzmann law, sigma T^4 / pi', () => {
    for (const T of [300, 3000, 10_000, 1e5]) {
      // Integrate in x = hc / (lambda k T) on a log grid, wide enough to capture the tails.
      let integral = 0;
      const n = 20_000;
      const lo = Math.log(1e-9 * (2.9e-3 / T) / 1e-2); // far into the Wien tail
      const hi = Math.log(1e4 * (2.9e-3 / T));
      const dl = (hi - lo) / n;
      for (let i = 0; i <= n; i += 1) {
        const lambda = Math.exp(lo + i * dl);
        const w = i === 0 || i === n ? 0.5 : 1;
        integral += w * planckSpectralRadiance(lambda, T) * lambda * dl;
      }
      const expected = (PHYSICAL_CONSTANTS.sigmaSB * T ** 4) / Math.PI;
      expect(Math.abs(integral / expected - 1), `T = ${T}`).toBeLessThan(1e-6);
    }
  });

  it("peaks where Wien's displacement law puts it, lambda_max T = 2.897771955e-3 m K", () => {
    for (const T of [3000, 5772, 10_000]) {
      let best = 0;
      let bestLambda = 0;
      for (let nm = 100; nm <= 3000; nm += 0.01) {
        const b = planckSpectralRadiance(nm * 1e-9, T);
        if (b > best) {
          best = b;
          bestLambda = nm * 1e-9;
        }
      }
      expect(Math.abs(bestLambda * T - 2.897771955e-3) / 2.897771955e-3, `T = ${T}`).toBeLessThan(1e-5);
    }
  });

  it('stays accurate in the Rayleigh-Jeans limit, where exp(x) - 1 would cancel', () => {
    const lambda = 1e-3;
    const T = 1e8;
    const rj = (2 * PHYSICAL_CONSTANTS.c * PHYSICAL_CONSTANTS.kB * T) / lambda ** 4;
    expect(Math.abs(planckSpectralRadiance(lambda, T) / rj - 1)).toBeLessThan(1e-6);
  });
});

describe('blackbody colour', () => {
  it('uses the CIE table at the standard 5 nm interval from 360 to 830 nm', () => {
    expect(CIE1931_CMF_5NM[0][0]).toBe(360);
    expect(CIE1931_CMF_5NM.at(-1)![0]).toBe(830);
    // y-bar is the photopic luminosity function: it peaks at 555 nm with value 1.
    const peak = CIE1931_CMF_5NM.find((row) => row[0] === 555)!;
    expect(peak[2]).toBeCloseTo(1, 6);
  });

  it("matches colour-science's independent 1 nm chromaticities to 5e-5", () => {
    for (const [T, x, y] of BLACKBODY_CHROMATICITY_REFERENCE) {
      const c = chromaticity(blackbodyXYZ(T));
      expect(Math.abs(c.x - x), `x at ${T} K`).toBeLessThan(5e-5);
      expect(Math.abs(c.y - y), `y at ${T} K`).toBeLessThan(5e-5);
    }
  });

  it('gives the photometric luminance of a blackbody in cd/m^2', () => {
    // Y is absolute: K_m times the luminous integral. A 5772 K blackbody (the Sun's
    // effective temperature) has a luminance near 2e9 cd/m^2; the Sun's measured disk-centre
    // luminance is about 1.6e9, lower mainly because of limb darkening.
    const Y = blackbodyXYZ(5772).Y;
    expect(Y).toBeGreaterThan(1.5e9);
    expect(Y).toBeLessThan(2.5e9);
  });

  it('maps a 6500 K blackbody close to display white', () => {
    const rgb = xyzToLinearSrgb(blackbodyXYZ(6500));
    const max = Math.max(rgb.r, rgb.g, rgb.b);
    // Not exactly D65 — the Planckian locus passes slightly off it — but within a few %.
    expect(Math.abs(rgb.r / max - 1)).toBeLessThan(0.08);
    expect(Math.abs(rgb.g / max - 1)).toBeLessThan(0.08);
    expect(Math.abs(rgb.b / max - 1)).toBeLessThan(0.08);
  });

  it('runs red when cool and blue when hot', () => {
    const cool = xyzToLinearSrgb(blackbodyXYZ(2000));
    const hot = xyzToLinearSrgb(blackbodyXYZ(40_000));
    expect(cool.r).toBeGreaterThan(cool.b);
    expect(hot.b).toBeGreaterThan(hot.r);
  });

  it('interpolates its lookup table to within 5e-6 of direct evaluation, in every channel', () => {
    // A dense sweep: an earlier, sparser version of this test passed a linear table whose
    // true worst case was 5.9e-3, because it happened to miss the worst points.
    for (let logT = Math.log(400); logT < Math.log(5e7); logT += 0.00737) {
      const T = Math.exp(logT);
      const direct = blackbodyXYZ(T);
      const fast = blackbodyXYZFast(T);
      expect(Math.abs(fast.Y / direct.Y - 1), `Y at ${T.toFixed(0)} K`).toBeLessThan(5e-6);
      expect(Math.abs(fast.X / direct.X - 1), `X at ${T.toFixed(0)} K`).toBeLessThan(5e-6);
      expect(Math.abs(fast.Z / direct.Z - 1), `Z at ${T.toFixed(0)} K`).toBeLessThan(5e-6);
    }
  });

  it('beams in-band as g in the Rayleigh-Jeans limit, far more steeply in the Wien limit', () => {
    // B_nu(g T) = g^3 B_{nu/g}(T): the observed spectrum is Planckian at g T. Its visible
    // luminance therefore scales very differently from the bolometric g^4 depending on
    // where the visible band sits. Deep in the Rayleigh-Jeans tail, B_nu is proportional to
    // T, so the in-band boost is just g — which is why a hot disk looks nearly symmetric in
    // visible light while its bolometric beaming is strong.
    const g = 1.4;
    const hot = blackbodyXYZ(1e7 * g).Y / blackbodyXYZ(1e7).Y;
    expect(Math.abs(hot / g - 1)).toBeLessThan(1e-3);
    const cool = blackbodyXYZ(3000 * g).Y / blackbodyXYZ(3000).Y;
    expect(cool).toBeGreaterThan(g ** 4);
    // And bolometrically, sigma (g T)^4 / sigma T^4 = g^4 exactly, whatever T is.
  });
});

