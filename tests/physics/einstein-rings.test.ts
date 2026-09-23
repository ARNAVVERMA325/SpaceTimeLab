import { describe, expect, it } from 'vitest';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { criticalImpactParameter, schwarzschild } from '../../src/physics/spacetimes/schwarzschild.js';
import { measureDeflection } from '../../src/physics/spacetimes/schwarzschild-rays.js';
import {
  observedAngle,
  solveEinsteinRing,
  strongDeflectionAngle,
  STRONG_DEFLECTION,
  weakFieldEinsteinAngle,
} from '../../src/physics/spacetimes/schwarzschild-lensing.js';

/**
 * ROADMAP.md 2A.3 — Einstein rings for axial rays, and the strong-deflection regime.
 *
 * Every reference value below was computed independently with mpmath at 60 significant
 * digits, from the orbit-equation quadrature with the turning point from the closed-form
 * root of r^3 - b^2 (r - 2M) = 0. No code is shared with the TypeScript under test.
 */

const M = 1;
const model = schwarzschild(M);
const bc = criticalImpactParameter(M);

function integrator(tolerance = 1e-12): Dopri5Integrator {
  return new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: tolerance, relative: tolerance } });
}

/** Exact deflection at b = b_c (1 + eps), mpmath. */
const EXACT_DEFLECTION: readonly (readonly [number, number])[] = [
  [1e-2, 4.2296618118375594197],
  [1e-3, 6.5106448266014707212],
  [1e-4, 8.8104863549956068302],
  [1e-6, 13.415285557794178382],
  [1e-8, 18.020450767385114685],
  [1e-10, 22.625620890944996692],
];

/** Einstein-ring impact parameters and observed angles for an on-axis source, mpmath. */
const RINGS: readonly {
  readonly rObserver: number;
  readonly rSource: number;
  readonly b: readonly number[];
  readonly theta: readonly number[];
}[] = [
  {
    rObserver: 100,
    rSource: 100,
    b: [15.747258207218576311, 5.202033046775031005, 5.196163367557764602, 5.196152443145311528],
    theta: [0.15652830653478451026, 0.051520288843910374619, 0.051462104868403572091, 0.051461996578875562398],
  },
  {
    rObserver: 1000,
    rSource: 1000,
    b: [46.233888675911726186, 5.202612595346619503, 5.196164441206025921],
    theta: [0.046204069423819072168, 0.005197430778756762231, 0.0051909889890986477966],
  },
  {
    rObserver: 1000,
    rSource: 100,
    b: [20.597665923099120489, 5.202315977069557990, 5.196163891825288479],
    theta: [0.020578510331689036662, 0.0051971344532443095957, 0.005190988440260171473],
  },
];

describe('strong-deflection limit (Bozza 2002)', () => {
  it('has the published coefficient b_bar = -pi + ln(216 (7 - 4 sqrt 3))', () => {
    expect(STRONG_DEFLECTION.aBar).toBe(1);
    expect(STRONG_DEFLECTION.bBar).toBeCloseTo(-0.40023003975526165, 14);
  });

  it('is approached by traced deflections as b -> b_c', () => {
    // alpha + ln(b/b_c - 1) -> b_bar. The mpmath residuals are 2.5e-2, 3.1e-3, 3.8e-4 and
    // 5.0e-6 at these epsilons; the traced ones must shrink the same way.
    let previous = Infinity;
    for (const epsilon of [1e-2, 1e-3, 1e-4, 1e-6]) {
      const m = measureDeflection({ model, integrator: integrator(), impactParameter: bc * (1 + epsilon), startRadius: 1000 });
      const residual = Math.abs(m.deflectionAngle! + Math.log(epsilon) - STRONG_DEFLECTION.bBar);
      expect(residual, `eps = ${epsilon}`).toBeLessThan(previous);
      previous = residual;
    }
    expect(previous).toBeLessThan(1e-5);
    // And the closed form agrees with itself as a function. The bound is set by the
    // input, not the formula: forming b/b_c - 1 from b = b_c (1 + 1e-6) cancels six digits,
    // leaving a relative error near eps / 1e-6 ~ 2e-10 in the logarithm's argument.
    expect(
      Math.abs(strongDeflectionAngle(M, bc * (1 + 1e-6)) - (-Math.log(1e-6) + STRONG_DEFLECTION.bBar)),
    ).toBeLessThan(1e-9);
  });

  it('bends light through several complete loops near the critical curve', () => {
    const m = measureDeflection({ model, integrator: integrator(), impactParameter: bc * (1 + 1e-10), startRadius: 1000 });
    expect(m.captured).toBe(false);
    expect(m.deflectionAngle! / (2 * Math.PI)).toBeGreaterThan(3.5);
  });
});

describe('near-critical deflection accuracy is limited by conditioning, not by defects', () => {
  it('matches the exact deflection with a constant effective impact-parameter error', () => {
    // Near b_c, alpha = -ln(b/b_c - 1) + ..., so d alpha / db = -1 / (b - b_c): any
    // integration error acts like an error delta_b in the impact parameter, amplified by
    // 1 / (b - b_c). Measured at tolerance 1e-12, |delta alpha| (b - b_c) is 6.6e-13 at
    // every epsilon from 1e-4 to 1e-10 — a fixed delta_b, which is exactly the absolute
    // error the ring solutions below show independently. So the growth of the relative
    // error as eps -> 0 is the problem's sensitivity, and the bound is placed on the
    // quantity the integrator actually controls.
    for (const [epsilon, exact] of EXACT_DEFLECTION) {
      const b = bc * (1 + epsilon);
      const m = measureDeflection({ model, integrator: integrator(), impactParameter: b, startRadius: 1000 });
      const effectiveImpactError = Math.abs(m.deflectionAngle! - exact) * (b - bc);
      expect(effectiveImpactError, `eps = ${epsilon}`).toBeLessThan(2e-12);
    }
  });
});

describe('Einstein rings for a source on the optical axis (ROADMAP.md 2A.3)', () => {
  for (const config of RINGS) {
    it(`solves the exact lens equation for r_O = ${config.rObserver}M, r_S = ${config.rSource}M`, () => {
      for (let n = 0; n < config.b.length; n += 1) {
        const ring = solveEinsteinRing(model, integrator(), n, config.rObserver, config.rSource);

        // The observable — the ring's angular radius — to 1e-12 rad at every order.
        expect(Math.abs(ring.observedAngle - config.theta[n]), `theta, n = ${n}`).toBeLessThan(1e-12);

        if (n === 0) {
          // The primary ring is conditioned the opposite way to the relativistic ones.
          // Far from b_c the sweep is shallow in b (dPhi/db ~ -4M/b^2 - 2/r ~ -0.004 at
          // r = 1000M), so a sweep error of order the tolerance moves b by a few 1e-10.
          // Bounded relatively; the angle above is what is observed.
          expect(Math.abs(ring.impactParameter - config.b[n]) / config.b[n], `b, n = 0`).toBeLessThan(1e-10);
        } else {
          // Near b_c the sweep is steep, which pins b to the fixed delta_b of the
          // conditioning law above.
          expect(Math.abs(ring.impactParameter - config.b[n]), `b, n = ${n}`).toBeLessThan(5e-12);
        }
        expect(ring.sweep).toBeCloseTo((2 * n + 1) * Math.PI, 8);
      }
    });
  }

  it('spaces successive relativistic rings by e^(-2 pi) towards the critical curve', () => {
    // (b_{n+1} - b_c) / (b_n - b_c) -> e^(-2 pi) ~ 0.0018674: one more loop needs an
    // impact parameter e^(2 pi) times closer to b_c. mpmath gives 0.00186117 for n = 1 -> 2
    // and 0.00186742 for n = 2 -> 3 at r_O = r_S = 100M; the leading-order formula's
    // correction terms account for the difference at n = 1.
    const [b1, b2, b3] = [1, 2, 3].map(
      (n) => solveEinsteinRing(model, integrator(), n, 100, 100).impactParameter,
    );
    const r12 = (b2 - bc) / (b1 - bc);
    const r23 = (b3 - bc) / (b2 - bc);
    expect(Math.abs(r12 - 0.00186117170650918) / r12).toBeLessThan(1e-6);
    expect(Math.abs(r23 - Math.exp(-2 * Math.PI)) / r23).toBeLessThan(1e-4);
  });

  it('converges on the weak-field Einstein angle as observer and source recede', () => {
    // Exact / weak-field: 1.107 at 100M and 1.033 at 1000M from mpmath. The correction is
    // of order M / b, so it keeps falling as the distances grow.
    const ratios = [100, 1000, 10000].map((d) => {
      const ring = solveEinsteinRing(model, integrator(), 0, d, d);
      return ring.observedAngle / weakFieldEinsteinAngle(M, d, d);
    });
    expect(ratios[0]).toBeCloseTo(0.15652830653478451026 / weakFieldEinsteinAngle(M, 100, 100), 10);
    expect(ratios[0]).toBeGreaterThan(ratios[1]);
    expect(ratios[1]).toBeGreaterThan(ratios[2]);
    expect(ratios[2]).toBeGreaterThan(1);
    expect(ratios[2]).toBeLessThan(1.02);
  });

  it('places every relativistic ring just outside the shadow edge', () => {
    // The shadow edge is where b = b_c. Rings of order >= 1 crowd towards it from outside.
    const edge = observedAngle(model, 100, bc);
    for (let n = 1; n <= 3; n += 1) {
      const ring = solveEinsteinRing(model, integrator(), n, 100, 100);
      expect(ring.observedAngle).toBeGreaterThan(edge);
      expect(ring.observedAngle - edge).toBeLessThan(1e-3);
    }
  });

  it('refuses a negative or fractional ring order', () => {
    expect(() => solveEinsteinRing(model, integrator(), -1, 100, 100)).toThrow(RangeError);
    expect(() => solveEinsteinRing(model, integrator(), 1.5, 100, 100)).toThrow(RangeError);
  });
});
