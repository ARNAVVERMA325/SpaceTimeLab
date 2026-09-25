import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { nullState } from '../../src/physics/core/phase-space.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { kerr, outerHorizonRadius } from '../../src/physics/spacetimes/kerr.js';
import {
  constantsFromImagePosition,
  shadowBoundaryPoint,
  shadowCurve,
  shadowExtent,
  sphericalPhotonOrbitConstants,
} from '../../src/physics/spacetimes/kerr-shadow.js';
import { criticalImpactParameter } from '../../src/physics/spacetimes/schwarzschild.js';
import { normalizationResidual } from '../../src/physics/validation/normalization.js';

/**
 * The Kerr shadow (ROADMAP.md 4A: "asymmetric shadow vs. reference values").
 *
 * Two independent things are checked here. First, the analytic critical curve against
 * results that have been in the literature since Bardeen (1973): it must become
 * Schwarzschild's circle of radius 3 sqrt(3) M as a -> 0, and for an extremal hole seen
 * edge-on it must run from alpha = -2M to alpha = +7M. Second — the part that actually
 * exercises the engine — traced null geodesics launched from a distant observer at image
 * positions just inside and just outside that curve must be captured and escape
 * respectively. The curve is never used to decide a ray's fate; the integration is.
 */

const M = 1;
const OBSERVER_RADIUS = 2000;

function tight(): Dopri5Integrator {
  return new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-11, relative: 1e-11 } });
}

/**
 * Initial data for the ray that reaches image position (alpha, beta), launched inward
 * from a distant observer and traced backward.
 */
function rayFromImagePosition(
  a: number,
  observerPolarAngle: number,
  alpha: number,
  beta: number,
): ReturnType<typeof nullState> {
  const model = kerr(M, a);
  const x: Vec4 = [0, OBSERVER_RADIUS, observerPolarAngle, 0];
  const { xi } = constantsFromImagePosition(a, observerPolarAngle, alpha, beta);

  const gInv = new Float64Array(16);
  model.inverseMetricInto(x, gInv);

  // Theta(theta_o) = eta + a^2 cos^2 - xi^2 cot^2 is beta^2 by construction; either sign
  // of p_theta gives a point of the same |beta|, and the shadow is symmetric in beta.
  const p_theta = -beta;
  const p_t = -1;
  const p_phi = xi;
  const rest =
    gInv[0] * p_t * p_t + 2 * gInv[3] * p_t * p_phi + gInv[15] * p_phi * p_phi + gInv[10] * p_theta * p_theta;
  const radialSquared = -rest / gInv[5];
  if (!(radialSquared >= 0)) {
    throw new RangeError(`rayFromImagePosition: no null ray reaches (${alpha}, ${beta}) from r = ${OBSERVER_RADIUS}.`);
  }
  // Ingoing: k^r = g^rr p_r must be negative, and g^rr > 0 outside the horizon.
  const p_r = -Math.sqrt(radialSquared);

  const k: [number, number, number, number] = [0, 0, 0, 0];
  const p = [p_t, p_r, p_theta, p_phi];
  for (let mu = 0; mu < 4; mu += 1) {
    let sum = 0;
    for (let nu = 0; nu < 4; nu += 1) sum += gInv[mu * 4 + nu] * p[nu];
    k[mu] = sum;
  }
  return nullState(x, k);
}

/** Whether the ray that reaches (alpha, beta) falls into the hole. */
function isCaptured(a: number, observerPolarAngle: number, alpha: number, beta: number): boolean {
  const model = kerr(M, a);
  const horizon = outerHorizonRadius(M, a);
  const initial = rayFromImagePosition(a, observerPolarAngle, alpha, beta);
  expect(Math.abs(normalizationResidual(model, initial))).toBeLessThan(1e-9);

  let captured = false;
  integrateGeodesic({
    model,
    integrator: tight(),
    initial,
    limits: { initialStep: 1e-2, parameterMax: 200_000, maxSteps: 400_000, maxStep: 50 },
    events: [{ id: 'escape', value: (x) => x[1] - OBSERVER_RADIUS * 1.5, direction: 1, terminal: true }],
    terminator: (x) => {
      // Stopped a hair outside the horizon: Boyer-Lindquist coordinates do not reach it,
      // and a ray this deep inside the photon region cannot come back out.
      captured = x[1] < horizon * 1.0001;
      return captured;
    },
  });
  return captured;
}

describe('the analytic critical curve', () => {
  it('becomes Schwarzschild’s circle of radius 3 sqrt(3) M as a -> 0', () => {
    const b_c = criticalImpactParameter(M);
    let previous = Number.POSITIVE_INFINITY;
    for (const a of [1e-2, 1e-3, 1e-4]) {
      let worst = 0;
      for (const point of shadowCurve(M, a, Math.PI / 3, 64)) {
        const radius = Math.hypot(point.alpha, point.beta);
        worst = Math.max(worst, Math.abs(radius / b_c - 1));
      }
      // The departure from a circle is first order in a, not second: the shadow shifts
      // before it deforms.
      expect(worst, `a = ${a}`).toBeLessThan(0.4 * a);
      expect(worst).toBeLessThan(previous);
      previous = worst;
    }
    expect(previous).toBeLessThan(1e-4);
    // The parametrization is singular at a = 0 and says so rather than dividing by zero.
    expect(() => shadowCurve(M, 0, Math.PI / 3)).toThrow(/singular at a = 0/);
  });

  it('runs from alpha = -2M to alpha = +7M for an extremal hole seen edge-on', () => {
    // The textbook figure (Bardeen 1973): the prograde side is flattened at -2M while the
    // retrograde side bulges to +7M. Both ends are exact values, not fits.
    const extent = shadowExtent(M, 1, Math.PI / 2, 40_000);
    expect(extent.alphaMin).toBeCloseTo(-2, 5);
    expect(extent.alphaMax).toBeCloseTo(7, 9);
    const atProgradeEnd = sphericalPhotonOrbitConstants(M, 1, 1 + 1e-9);
    expect(atProgradeEnd.xi).toBeCloseTo(2, 6);
  });

  it('keeps the same height at every spin while the width shrinks and the centre shifts', () => {
    // Seen edge-on, the shadow's extent along the spin axis stays 2 * 3 sqrt(3) M however
    // fast the hole spins; only the width and the displacement change. The asymmetry is
    // entirely in alpha.
    const height = 2 * criticalImpactParameter(M);
    let previousWidth = Number.POSITIVE_INFINITY;
    let previousCentre = Number.NEGATIVE_INFINITY;
    for (const a of [0.1, 0.3, 0.6, 0.9, 0.998]) {
      const extent = shadowExtent(M, a, Math.PI / 2, 40_000);
      expect(Math.abs(extent.height / height - 1), `height at a = ${a}`).toBeLessThan(1e-6);
      expect(extent.width).toBeLessThan(previousWidth);
      expect(extent.alphaCentre).toBeGreaterThan(previousCentre);
      previousWidth = extent.width;
      previousCentre = extent.alphaCentre;
    }
    expect(previousCentre).toBeGreaterThan(2);
  });

  it('loses its displacement as sin(theta_o), becoming a centred circle down the axis', () => {
    // Only the projected rotation displaces the shadow, so the displacement goes as
    // sin(theta_o) and the curve closes back into a circle looking along the axis. The
    // ratio below is what makes it a measurement rather than a limit: it stays constant
    // across three decades of inclination.
    const a = 0.9;
    const ratios: number[] = [];
    for (const theta of [1e-2, 1e-3, 1e-4]) {
      const extent = shadowExtent(M, a, theta, 20_000);
      ratios.push(extent.alphaCentre / (a * Math.sin(theta)));
      expect(Math.abs(extent.width / extent.height - 1), `theta_o = ${theta}`).toBeLessThan(10 * theta);
    }
    for (const ratio of ratios) expect(Math.abs(ratio / ratios[ratios.length - 1] - 1)).toBeLessThan(1e-3);
    expect(ratios[ratios.length - 1]).toBeGreaterThan(1.5);
  });
});

describe('traced rays against the analytic curve', () => {
  const cases = [
    { a: 0.9, inclination: Math.PI / 2, label: 'a = 0.9, edge-on' },
    { a: 0.9, inclination: Math.PI / 3, label: 'a = 0.9, 60 degrees' },
    { a: 0.5, inclination: Math.PI / 2, label: 'a = 0.5, edge-on' },
  ];

  for (const { a, inclination, label } of cases) {
    it(`captures inside and releases outside the curve: ${label}`, () => {
      const extent = shadowExtent(M, a, inclination, 8192);
      const centre = { alpha: extent.alphaCentre, beta: 0 };
      const curve = shadowCurve(M, a, inclination, 1024);
      // Sample the curve evenly rather than taking neighbouring points.
      const sampled = [0.08, 0.28, 0.5, 0.72, 0.92].map((f) => curve[Math.floor(f * (curve.length - 1))]);

      for (const point of sampled) {
        for (const [scale, shouldBeCaptured] of [[0.995, true], [1.005, false]] as const) {
          const alpha = centre.alpha + (point.alpha - centre.alpha) * scale;
          const beta = centre.beta + (point.beta - centre.beta) * scale;
          expect(
            isCaptured(a, inclination, alpha, beta),
            `${label}: (${alpha.toFixed(4)}, ${beta.toFixed(4)}) at ${scale} of the boundary`,
          ).toBe(shouldBeCaptured);
        }
      }
    });
  }

  it('locates the boundary by bisection to within 1e-4 of the analytic curve', () => {
    const a = 0.9;
    const inclination = Math.PI / 2;
    const extent = shadowExtent(M, a, inclination, 8192);
    const centreAlpha = extent.alphaCentre;

    for (const orbitFraction of [0.25, 0.6]) {
      const curve = shadowCurve(M, a, inclination, 512);
      const target = curve[Math.floor(orbitFraction * (curve.length - 1))];
      let lo = 0.9; // certainly inside
      let hi = 1.1; // certainly outside
      for (let i = 0; i < 22; i += 1) {
        const mid = 0.5 * (lo + hi);
        const alpha = centreAlpha + (target.alpha - centreAlpha) * mid;
        const beta = target.beta * mid;
        if (isCaptured(a, inclination, alpha, beta)) lo = mid;
        else hi = mid;
      }
      const boundary = 0.5 * (lo + hi);
      expect(Math.abs(boundary - 1), `boundary scale at orbit fraction ${orbitFraction}`).toBeLessThan(1e-4);
    }
  });

  it('shows the shadow displaced toward the retrograde side, as the curve says', () => {
    // A direct consequence of frame dragging rather than a restatement of the formula:
    // the point at alpha = -4M, beta = 0 is captured for a = 0.9 edge-on, while its
    // mirror image at alpha = +4M is not, because the shadow has shifted to +alpha.
    const a = 0.9;
    const inclination = Math.PI / 2;
    expect(isCaptured(a, inclination, 4, 0)).toBe(true);
    expect(isCaptured(a, inclination, -4, 0)).toBe(false);
    // And with the spin reversed, so is the asymmetry.
    expect(isCaptured(-a, inclination, -4, 0)).toBe(true);
    expect(isCaptured(-a, inclination, 4, 0)).toBe(false);
  });

  it('places the boundary point of a spherical photon orbit on a ray that neither falls nor escapes quickly', () => {
    // The orbits the curve is built from are real: a ray aimed exactly at the boundary
    // spends a long time near the photon region before its fate is decided, which is why
    // the edge of a shadow is sharp but infinitely fine-structured.
    const a = 0.9;
    const inclination = Math.PI / 2;
    const point = shadowBoundaryPoint(M, a, inclination, 2.5);
    expect(point).toBeDefined();
    if (!point) return;

    const model = kerr(M, a);
    const result = integrateGeodesic({
      model,
      integrator: tight(),
      initial: rayFromImagePosition(a, inclination, point.alpha, point.beta),
      limits: { initialStep: 1e-2, parameterMax: 200_000, maxSteps: 400_000, maxStep: 50 },
      events: [{ id: 'escape', value: (x) => x[1] - OBSERVER_RADIUS * 1.5, direction: 1, terminal: true }],
      terminator: (x) => x[1] < outerHorizonRadius(M, a) * 1.0001,
    });
    // However it ends, it has lingered: many times the steps a ray that misses the hole
    // by a wide margin needs.
    const wideMiss = integrateGeodesic({
      model,
      integrator: tight(),
      initial: rayFromImagePosition(a, inclination, 40, 0),
      limits: { initialStep: 1e-2, parameterMax: 200_000, maxSteps: 400_000, maxStep: 50 },
      events: [{ id: 'escape', value: (x) => x[1] - OBSERVER_RADIUS * 1.5, direction: 1, terminal: true }],
    });
    expect(result.steps).toBeGreaterThan(3 * wideMiss.steps);
  });
});
