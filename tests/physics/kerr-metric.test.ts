import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import {
  carterConstant,
  equatorialPhotonOrbitRadius,
  ergosphereRadius,
  frameDraggingOmega,
  innerHorizonRadius,
  iscoRadius,
  kerr,
  kretschmann,
  outerHorizonRadius,
} from '../../src/physics/spacetimes/kerr.js';
import { schwarzschild, kretschmann as schwarzschildKretschmann } from '../../src/physics/spacetimes/schwarzschild.js';
import { KERR_REFERENCE_POINTS } from '../fixtures/kerr-reference.js';

/**
 * Kerr in Boyer-Lindquist coordinates (ROADMAP.md 4A).
 *
 * The closed forms in `kerr.ts` are checked against an independent SymPy derivation that
 * inverts the metric by matrix inversion and builds the Christoffel symbols from
 * CLAUDE.md §2's definition, sharing no algebraic route with them
 * (`scripts/generate_kerr_reference.py`). The orbit radii are checked against the
 * equations that define them rather than against remembered numbers.
 */

const SPINS = [0, 0.3, 0.5, 0.9, 0.998, 1];

function relativeMiss(measured: number, expected: number): number {
  const scale = Math.max(Math.abs(expected), 1e-8);
  return Math.abs(measured - expected) / scale;
}

describe('Kerr metric against an independent symbolic derivation', () => {
  it('reproduces g_mu_nu, g^{mu nu}, d_alpha g^{mu nu} and the Christoffel symbols', () => {
    for (const point of KERR_REFERENCE_POINTS) {
      const model = kerr(point.mass, point.spin);
      const x: Vec4 = [0, point.r, point.theta, 0];
      const label = `M = ${point.mass}, a = ${point.spin}, r = ${point.r}, theta = ${point.theta}`;

      const metric = model.metricAt(x);
      for (let i = 0; i < 16; i += 1) {
        expect(relativeMiss(metric.g_mu_nu[i], point.metric[i]), `g_mu_nu[${i}] at ${label}`).toBeLessThan(1e-13);
        expect(
          relativeMiss(metric.g_inv_mu_nu[i], point.inverseMetric[i]),
          `g^{mu nu}[${i}] at ${label}`,
        ).toBeLessThan(1e-13);
      }

      const derivatives = new Float64Array(64);
      model.inverseMetricDerivativesInto(x, derivatives);
      for (let i = 0; i < 64; i += 1) {
        expect(
          relativeMiss(derivatives[i], point.inverseMetricDerivatives[i]),
          `d g^{mu nu}[${i}] at ${label}`,
        ).toBeLessThan(1e-12);
      }

      const christoffel = model.christoffelAt(x);
      for (let i = 0; i < 64; i += 1) {
        expect(
          relativeMiss(christoffel.components[i], point.christoffel[i]),
          `Gamma[${i}] at ${label}`,
        ).toBeLessThan(1e-12);
      }
    }
  });

  it('inverts its own metric: g_mu_nu g^{nu sigma} = delta^sigma_mu', () => {
    for (const point of KERR_REFERENCE_POINTS) {
      const metric = kerr(point.mass, point.spin).metricAt([0, point.r, point.theta, 0]);
      for (let mu = 0; mu < 4; mu += 1) {
        for (let sigma = 0; sigma < 4; sigma += 1) {
          let sum = 0;
          for (let nu = 0; nu < 4; nu += 1) sum += metric.g_mu_nu[mu * 4 + nu] * metric.g_inv_mu_nu[nu * 4 + sigma];
          expect(Math.abs(sum - (mu === sigma ? 1 : 0))).toBeLessThan(1e-13);
        }
      }
    }
  });

  it('has determinant -Sigma^2 sin^2(theta), which does not vanish at the horizon', () => {
    // CLAUDE.md §6.1: the horizon must never be detected from a vanishing determinant.
    // In Kerr the determinant is -Sigma^2 sin^2(theta), entirely independent of Delta, so
    // at r_+ it is a perfectly ordinary negative number while Delta has vanished.
    const M = 1;
    const a = 0.9;
    const model = kerr(M, a);
    const horizon = outerHorizonRadius(M, a);
    const determinantAt = (r: number, theta: number): { measured: number; expected: number } => {
      const g = model.metricAt([0, r, theta, 0]);
      // The (t, phi) block's determinant times g_rr g_theta_theta.
      const block = g.g_mu_nu[0] * g.g_mu_nu[15] - g.g_mu_nu[3] * g.g_mu_nu[12];
      const sigma = r * r + a * a * Math.cos(theta) ** 2;
      return {
        measured: block * g.g_mu_nu[5] * g.g_mu_nu[10],
        expected: -(sigma * sigma) * Math.sin(theta) ** 2,
      };
    };

    for (const [r, theta] of [[horizon + 1e-3, 1.1], [3, 0.7], [20, Math.PI / 2]] as const) {
      const { measured, expected } = determinantAt(r, theta);
      expect(relativeMiss(measured, expected), `r = ${r}`).toBeLessThan(1e-12);
    }

    // A hair above the horizon the determinant is still an ordinary number of order
    // Sigma^2, nowhere near zero, while Delta has all but vanished. The relative bound is
    // looser here only because r itself resolves an offset of 1e-9 from r_+ to about seven
    // digits in binary64; that is the coordinate, not the formula.
    const nearHorizon = determinantAt(horizon + 1e-9, 1.1);
    expect(relativeMiss(nearHorizon.measured, nearHorizon.expected)).toBeLessThan(1e-5);
    expect(Math.abs(nearHorizon.measured)).toBeGreaterThan(1);
    // Delta is evaluated in factored form, (r - r_+)(r - r_-), which keeps full relative
    // accuracy where the expanded polynomial cancels to nothing.
    for (const epsilon of [1e-6, 1e-9, 1e-12]) {
      const r = horizon + epsilon;
      // The offset r actually carries, not the one asked for: adding 1e-12 to 1.4358...
      // is itself rounded, and no formula can recover what the coordinate no longer holds.
      const offset = r - horizon;
      const expected = offset * (offset + horizon - innerHorizonRadius(1, a));
      expect(relativeMiss(model.deltaFunction(r), expected), `Delta at r_+ + ${epsilon}`).toBeLessThan(1e-14);
      // The expanded polynomial, for comparison: by 1e-9 from the horizon it has already
      // lost more relative accuracy than the factored form ever does.
      const expanded = r * r - 2 * r + a * a;
      if (epsilon <= 1e-9) expect(relativeMiss(expanded, expected)).toBeGreaterThan(1e-9);
    }

    const atHorizon = model.metricAt([0, horizon + 1e-9, Math.PI / 2, 0]);
    expect(Math.abs(atHorizon.g_mu_nu[5])).toBeGreaterThan(1e6); // g_rr = Sigma/Delta diverges
    expect(model.deltaFunction(horizon)).toBeCloseTo(0, 12);
  });
});

describe('the Schwarzschild limit a -> 0', () => {
  it('matches the Schwarzschild metric and connection component for component', () => {
    const flat = kerr(1, 0);
    const reference = schwarzschild(1);
    for (const [r, theta] of [[3, 0.6], [8, Math.PI / 2], [50, 2.4]] as const) {
      const x: Vec4 = [0, r, theta, 0];
      const kerrMetric = flat.metricAt(x);
      const schwarzschildMetric = reference.metricAt(x);
      for (let i = 0; i < 16; i += 1) {
        expect(relativeMiss(kerrMetric.g_mu_nu[i], schwarzschildMetric.g_mu_nu[i]), `g[${i}] at r = ${r}`).toBeLessThan(1e-14);
        expect(relativeMiss(kerrMetric.g_inv_mu_nu[i], schwarzschildMetric.g_inv_mu_nu[i])).toBeLessThan(1e-14);
      }
      const kerrChristoffel = flat.christoffelAt(x);
      const schwarzschildChristoffel = reference.christoffelAt(x);
      for (let i = 0; i < 64; i += 1) {
        expect(
          relativeMiss(kerrChristoffel.components[i], schwarzschildChristoffel.components[i]),
          `Gamma[${i}] at r = ${r}`,
        ).toBeLessThan(1e-14);
      }
    }
  });

  it('loses the frame dragging and the cross term linearly in a', () => {
    const r = 10;
    const theta = Math.PI / 2;
    let previous = Number.POSITIVE_INFINITY;
    for (const a of [0.1, 0.01, 0.001, 0.0001]) {
      const g = kerr(1, a).metricAt([0, r, theta, 0]);
      const ratio = Math.abs(g.g_mu_nu[3]) / a;
      // g_t_phi = -2 M a r sin^2(theta) / Sigma -> -2 M a / r as a -> 0.
      expect(relativeMiss(ratio, (2 * 1) / r), `a = ${a}`).toBeLessThan(1e-3);
      const omega = frameDraggingOmega(1, a, r, theta);
      expect(omega).toBeLessThan(previous);
      previous = omega;
    }
    expect(frameDraggingOmega(1, 0, r, theta)).toBe(0);
  });

  it('reduces the Kretschmann scalar to 48 M^2 / r^6', () => {
    for (const r of [3, 10, 100]) {
      expect(relativeMiss(kretschmann(1, 0, r, 1.0), schwarzschildKretschmann(1, r))).toBeLessThan(1e-15);
    }
  });
});

describe('horizons, ergosphere and curvature', () => {
  it('places the horizons at M +- sqrt(M^2 - a^2), meeting at a = M', () => {
    for (const a of SPINS) {
      const outer = outerHorizonRadius(1, a);
      const inner = innerHorizonRadius(1, a);
      expect(outer * inner).toBeCloseTo(a * a, 12); // r_+ r_- = a^2
      expect(outer + inner).toBeCloseTo(2, 12); // r_+ + r_- = 2M
      expect(kerr(1, a).deltaFunction(outer)).toBeCloseTo(0, 12);
    }
    expect(outerHorizonRadius(1, 1)).toBeCloseTo(1, 15);
  });

  it('puts the static limit outside the horizon, touching it on the axis', () => {
    const a = 0.8;
    const horizon = outerHorizonRadius(1, a);
    expect(ergosphereRadius(1, a, 0)).toBeCloseTo(horizon, 12);
    expect(ergosphereRadius(1, a, Math.PI / 2)).toBeCloseTo(2, 12);
    for (const theta of [0.3, 0.9, 1.5, 2.6]) {
      expect(ergosphereRadius(1, a, theta)).toBeGreaterThanOrEqual(horizon - 1e-12);
    }
  });

  it('keeps every curvature invariant finite at the horizon and divergent only on the ring', () => {
    const a = 0.9;
    const horizon = outerHorizonRadius(1, a);
    expect(Number.isFinite(kretschmann(1, a, horizon, 1.0))).toBe(true);
    expect(Math.abs(kretschmann(1, a, horizon, 1.0))).toBeLessThan(100);
    // Approaching the ring r -> 0 in the equatorial plane, K must blow up.
    let previous = 0;
    for (const r of [0.1, 0.01, 0.001]) {
      const K = Math.abs(kretschmann(1, a, r, Math.PI / 2));
      expect(K).toBeGreaterThan(previous);
      previous = K;
    }
    expect(previous).toBeGreaterThan(1e12);
  });

  it('refuses a spin beyond extremal rather than describing a naked singularity', () => {
    expect(() => kerr(1, 1.0001)).toThrow(/\|a\| <= M/);
    expect(() => kerr(1, -1.5)).toThrow(/\|a\| <= M/);
    expect(() => kerr(1, 1)).not.toThrow();
  });

  it('reports the chart breakdown from Delta, not from the determinant', () => {
    const model = kerr(1, 0.9);
    const horizon = outerHorizonRadius(1, 0.9);
    const inside = model.domainCheck([0, horizon * 0.99, 1.1, 0]);
    expect(inside.inDomain).toBe(false);
    if (!inside.inDomain) {
      expect(inside.code).toBe('coordinate-breakdown');
      expect(inside.reason).toMatch(/Delta/);
      expect(inside.reason).toMatch(/curvature invariants are finite/);
    }
    expect(model.domainCheck([0, 5, 1e-9, 0]).inDomain).toBe(false);
    expect(model.domainCheck([0, 5, 1.1, 0]).inDomain).toBe(true);
  });
});

describe('closed-form orbit radii', () => {
  it('solves the defining cubic for the equatorial photon orbit', () => {
    // r^3 - 6 M r^2 + 9 M^2 r - 4 M a^2 = 0, the condition for a circular null orbit.
    for (const a of SPINS) {
      for (const sense of ['prograde', 'retrograde'] as const) {
        const r = equatorialPhotonOrbitRadius(1, a, sense);
        const residual = r ** 3 - 6 * r * r + 9 * r - 4 * a * a;
        expect(Math.abs(residual), `a = ${a}, ${sense}`).toBeLessThan(1e-12);
      }
    }
    expect(equatorialPhotonOrbitRadius(1, 0, 'prograde')).toBeCloseTo(3, 12);
    expect(equatorialPhotonOrbitRadius(1, 1, 'prograde')).toBeCloseTo(1, 9);
    expect(equatorialPhotonOrbitRadius(1, 1, 'retrograde')).toBeCloseTo(4, 12);
  });

  it('solves the marginal-stability condition for the ISCO', () => {
    // r^2 - 6 M r +- 8 a sqrt(M r) - 3 a^2 = 0, upper sign prograde.
    for (const a of SPINS) {
      for (const [sense, sign] of [['prograde', 1], ['retrograde', -1]] as const) {
        const r = iscoRadius(1, a, sense);
        const residual = r * r - 6 * r + sign * 8 * a * Math.sqrt(r) - 3 * a * a;
        expect(Math.abs(residual), `a = ${a}, ${sense}`).toBeLessThan(1e-11);
      }
    }
    expect(iscoRadius(1, 0, 'prograde')).toBeCloseTo(6, 12);
    expect(iscoRadius(1, 0, 'retrograde')).toBeCloseTo(6, 12);
    expect(iscoRadius(1, 1, 'prograde')).toBeCloseTo(1, 9);
    expect(iscoRadius(1, 1, 'retrograde')).toBeCloseTo(9, 9);
  });

  it('keeps the ISCO outside the photon orbit, which is outside the horizon', () => {
    for (const a of SPINS) {
      for (const sense of ['prograde', 'retrograde'] as const) {
        expect(iscoRadius(1, a, sense)).toBeGreaterThanOrEqual(equatorialPhotonOrbitRadius(1, a, sense) - 1e-12);
        expect(equatorialPhotonOrbitRadius(1, a, sense)).toBeGreaterThanOrEqual(outerHorizonRadius(1, a) - 1e-12);
      }
    }
  });

  it('approaches the Lense-Thirring form omega = 2 M a / r^3 far from the hole', () => {
    const a = 0.7;
    let previous = Number.POSITIVE_INFINITY;
    for (const r of [50, 500, 5000]) {
      const miss = relativeMiss(frameDraggingOmega(1, a, r, Math.PI / 2), (2 * a) / r ** 3);
      expect(miss).toBeLessThan(previous);
      previous = miss;
    }
    expect(previous).toBeLessThan(1e-6);
  });
});

describe('the Carter constant, in this project’s convention', () => {
  it('vanishes for any equatorial trajectory', () => {
    // Q = p_theta^2 + cos^2(theta)[a^2(mu^2 - E^2) + L_z^2/sin^2(theta)], which is
    // identically 0 in the equatorial plane with no polar motion. That is the whole
    // reason for choosing this convention over K = Q + (L_z - a E)^2.
    const model = kerr(1, 0.8);
    for (const r of [4, 10, 40]) {
      const x: Vec4 = [0, r, Math.PI / 2, 0];
      const g = new Float64Array(16);
      model.metricInto(x, g);
      // A photon with some radial and azimuthal motion, normalized to be null.
      const kPhi = 0.02;
      const kR = 0.5;
      // Solve g_tt (k^t)^2 + 2 g_t_phi k^t k^phi + g_phi_phi (k^phi)^2 + g_rr (k^r)^2 = 0.
      const A = g[0];
      const B = 2 * g[3] * kPhi;
      const C = g[15] * kPhi * kPhi + g[5] * kR * kR;
      const kT = (-B - Math.sqrt(B * B - 4 * A * C)) / (2 * A);
      const tangent: Vec4 = [kT, kR, 0, kPhi];
      expect(Math.abs(carterConstant(1, 0.8, x, tangent)), `r = ${r}`).toBeLessThan(1e-20);
    }
  });

  it('is positive when there is polar motion, and grows with it', () => {
    const model = kerr(1, 0.8);
    const x: Vec4 = [0, 12, Math.PI / 2, 0];
    const g = new Float64Array(16);
    model.metricInto(x, g);
    let previous = 0;
    for (const kTheta of [0.005, 0.01, 0.02]) {
      const A = g[0];
      const C = g[10] * kTheta * kTheta;
      const kT = Math.sqrt(-C / A);
      const Q = carterConstant(1, 0.8, x, [kT, 0, kTheta, 0]);
      expect(Q).toBeGreaterThan(previous);
      previous = Q;
    }
  });
});
