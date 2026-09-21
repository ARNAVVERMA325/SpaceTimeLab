import { describe, expect, it } from 'vitest';
import {
  DEFLECTION_MIN_TURNING_POINT_FACTOR,
  deflectionAngleExact,
  deflectionAngleWeakField,
  impactParameterForTurningPoint,
  turningPointForImpactParameter,
} from '../../src/physics/spacetimes/schwarzschild-analytic.js';
import {
  criticalImpactParameter,
  kretschmann,
  photonSphereRadius,
} from '../../src/physics/spacetimes/schwarzschild.js';

/**
 * The analytic reference itself, checked before anything is checked against it.
 *
 * A reference that is wrong is worse than no reference: it would certify a broken
 * integrator. These values were produced independently, by evaluating both the smoothed
 * quadrature and the original singular integral with an adaptive rule in a separate
 * computer-algebra session; the two agreed to about 1e-12.
 */

const M = 1;

/** Independently computed exact deflection angles, keyed by turning point r_0 / M. */
const REFERENCE_DEFLECTION: readonly (readonly [number, number])[] = [
  [1e7, 4.000000766e-7],
  [1e5, 4.000077811e-5],
  [1e4, 4.000778268e-4],
  [1000, 0.004007798117],
  [100, 0.04079561289],
  [50, 0.08325597088],
  [20, 0.2218761043],
  [10, 0.5002356566],
  [6, 1.014875432],
  [4, 2.184100188],
  [3.5, 3.206122742],
];

describe('Schwarzschild analytic reference', () => {
  it('reproduces independently computed exact deflection angles', () => {
    // Combined absolute and relative bound, not a relative one alone. The deflection is
    // formed as 4 * Integral - pi, a difference of two quantities of order pi, so in the
    // weak field almost everything cancels: at r_0 = 10^7 M the deflection is 4e-7 and
    // only about eight significant figures survive. The absolute accuracy stays near
    // 1e-15 throughout, which is what the absolute term captures.
    const absoluteTolerance = 1e-14;
    const relativeTolerance = 1e-9;

    for (const [r0, expected] of REFERENCE_DEFLECTION) {
      const alpha = deflectionAngleExact(M, r0);
      const bound = absoluteTolerance + relativeTolerance * Math.abs(expected);
      expect(
        Math.abs(alpha - expected),
        `r_0 = ${r0}M: got ${alpha}, expected ${expected}, bound ${bound}`,
      ).toBeLessThan(bound);
    }
  });

  it('keeps absolute accuracy where the weak-field cancellation destroys relative accuracy', () => {
    // Pins the behaviour documented above, so a future change that traded absolute
    // accuracy for something else would be caught.
    const alpha = deflectionAngleExact(M, 1e7);
    expect(Math.abs(alpha - 4.000000766e-7)).toBeLessThan(1e-14);
    // Relative accuracy at this radius is genuinely worse than 1e-10, and claiming
    // otherwise would be claiming precision the formula does not have.
    expect(Math.abs(alpha - 4.000000766e-7) / 4.000000766e-7).toBeGreaterThan(1e-10);
  });

  it('is converged in the number of quadrature nodes', () => {
    // Doubling the node count must not move the answer: the integrand is smooth after
    // the singularity-removing substitution, so Gauss-Legendre converges fast.
    for (const r0 of [3.5, 6, 20, 1000]) {
      const coarse = deflectionAngleExact(M, r0, 32);
      const fine = deflectionAngleExact(M, r0, 256);
      expect(Math.abs(coarse - fine) / fine, `r_0 = ${r0}M`).toBeLessThan(1e-12);
    }
  });

  it('approaches the weak-field limit 4M/b as the ray passes further out', () => {
    // ROADMAP.md 2A.4's benchmark, stated as a limit rather than an equality.
    let previousError = Infinity;
    for (const r0 of [100, 1000, 1e4, 1e5, 1e6]) {
      const b = impactParameterForTurningPoint(M, r0);
      const exact = deflectionAngleExact(M, r0);
      const weak = deflectionAngleWeakField(M, b);
      const relativeError = Math.abs(exact - weak) / exact;
      // Monotone improvement, and the error falls roughly as M/b.
      expect(relativeError, `r_0 = ${r0}M`).toBeLessThan(previousError);
      previousError = relativeError;
    }
    // By r_0 = 10^6 M the weak-field formula is good to better than a part in 10^5.
    const b = impactParameterForTurningPoint(M, 1e6);
    const exact = deflectionAngleExact(M, 1e6);
    expect(Math.abs(exact - deflectionAngleWeakField(M, b)) / exact).toBeLessThan(1e-5);
  });

  it('shows the weak-field formula failing in the strong field, as it must', () => {
    // Guards against the benchmark being trivially satisfied: if 4M/b agreed everywhere
    // the test above would prove nothing about the strong-field regime.
    const b = impactParameterForTurningPoint(M, 10);
    const exact = deflectionAngleExact(M, 10);
    const weak = deflectionAngleWeakField(M, b);
    expect(weak / exact).toBeLessThan(0.75);
  });

  it('relates impact parameter and turning point consistently in both directions', () => {
    for (const r0 of [3.2, 4, 6, 10, 100, 1e4]) {
      const b = impactParameterForTurningPoint(M, r0);
      const recovered = turningPointForImpactParameter(M, b);
      expect(recovered).toBeDefined();
      expect(Math.abs(recovered! - r0) / r0, `r_0 = ${r0}M`).toBeLessThan(1e-12);
    }
  });

  it('places the critical impact parameter at the photon sphere', () => {
    const bCritical = criticalImpactParameter(M);
    expect(bCritical).toBeCloseTo(3 * Math.sqrt(3), 15);
    expect(bCritical).toBeCloseTo(5.196152422706632, 14);

    // The turning point tends to the photon sphere as b tends to b_c from above.
    const justAbove = turningPointForImpactParameter(M, bCritical * (1 + 1e-9));
    expect(justAbove).toBeDefined();
    expect(justAbove!).toBeGreaterThan(photonSphereRadius(M));
    expect(justAbove! / photonSphereRadius(M) - 1).toBeLessThan(1e-3);
  });

  it('reports no turning point below the critical impact parameter (capture)', () => {
    const bCritical = criticalImpactParameter(M);
    for (const b of [0.1, 1, 3, 5, bCritical * (1 - 1e-9)]) {
      expect(turningPointForImpactParameter(M, b), `b = ${b}M should be captured`).toBeUndefined();
    }
  });

  it('refuses a turning point too close to the photon sphere rather than guessing', () => {
    expect(() => deflectionAngleExact(M, photonSphereRadius(M))).toThrow(RangeError);
    expect(() => deflectionAngleExact(M, 2.5)).toThrow(RangeError);
    expect(() => deflectionAngleExact(M, DEFLECTION_MIN_TURNING_POINT_FACTOR * M * 0.99)).toThrow(
      RangeError,
    );
  });

  it('diverges as the turning point approaches the photon sphere', () => {
    const near = deflectionAngleExact(M, 3.06);
    const far = deflectionAngleExact(M, 10);
    // Strong-field deflection exceeds a full half-turn and keeps growing.
    expect(near).toBeGreaterThan(Math.PI);
    expect(near).toBeGreaterThan(far * 5);
  });

  it('gives a Kretschmann scalar that is finite at the horizon and divergent at r = 0', () => {
    // CLAUDE.md §6.1: the horizon is a coordinate singularity, not a curvature one.
    const atHorizon = kretschmann(M, 2 * M);
    expect(Number.isFinite(atHorizon)).toBe(true);
    expect(atHorizon).toBeCloseTo(48 / 64, 15);
    expect(kretschmann(M, 0.01 * M)).toBeGreaterThan(1e10);
    // It falls off as r^-6.
    expect(kretschmann(M, 2 * M) / kretschmann(M, 4 * M)).toBeCloseTo(64, 10);
  });
});
