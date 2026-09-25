import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { generateNullRay } from '../../src/physics/observer/observer.js';
import { zamoObserver } from '../../src/physics/observer/zamo.js';
import { kerr, outerHorizonRadius } from '../../src/physics/spacetimes/kerr.js';
import { constantsOfMotion, isCaptured, radialPotential } from '../../src/physics/spacetimes/kerr-rays.js';

/**
 * The Kerr capture test.
 *
 * It claims something exact: an inward-moving photon reaches the horizon precisely when
 * the radial potential R(r) has no root between r_+ and where it is. The claim is checked
 * two ways — against a brute-force scan of R, and against what the integrator actually
 * does with the ray.
 */

const M = 1;
const a = 0.9;
const model = kerr(M, a);
const horizon = outerHorizonRadius(M, a);

function rays(): { ray: ReturnType<typeof generateNullRay>; label: string }[] {
  const out: { ray: ReturnType<typeof generateNullRay>; label: string }[] = [];
  for (const r of [2.2, 4, 8, 25]) {
    for (const theta of [Math.PI / 2, 1.0, 2.2]) {
      const observer = zamoObserver(model, [0, r, theta, 0]);
      for (const raw of [
        [-1, 0, 0],
        [-0.96, 0.28, 0],
        [-0.9, -0.3, 0.32],
        [-0.8, 0.5, 0.33],
        [-0.6, 0.7, 0.39],
        [-0.3, 0.9, 0.32],
        [0.2, 0.9, 0.39],
      ] as const) {
        const norm = Math.hypot(...raw);
        out.push({
          ray: generateNullRay(observer, [raw[0] / norm, raw[1] / norm, raw[2] / norm]),
          label: `r = ${r}, theta = ${theta}, direction ${raw.join(', ')}`,
        });
      }
    }
  }
  return out;
}

describe('the exact capture test', () => {
  it('agrees with a brute-force scan of the radial potential', () => {
    for (const { ray, label } of rays()) {
      const r = ray.position_x[1];
      const inward = ray.tangent[1] < 0;
      const constants = constantsOfMotion(model, ray.position_x, ray.tangent);

      let rootBetween = false;
      const samples = 20_000;
      for (let i = 0; i <= samples; i += 1) {
        const sample = horizon + ((r - horizon) * i) / samples;
        if (radialPotential(M, a, constants, sample) <= 0) {
          rootBetween = true;
          break;
        }
      }
      const expected = inward && !rootBetween;
      expect(isCaptured(model, ray.position_x, ray.tangent), label).toBe(expected);
    }
  });

  it('never claims an outward-moving ray is captured', () => {
    for (const { ray } of rays()) {
      if (ray.tangent[1] <= 0) continue;
      expect(isCaptured(model, ray.position_x, ray.tangent)).toBe(false);
    }
  });

  it('matches what the integrator does with the same ray', () => {
    // The test that matters: stopping a ray early on this criterion must reach the same
    // verdict as following it all the way down to the horizon.
    const integrator = (): Dopri5Integrator =>
      new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-11, relative: 1e-11 } });
    const limits = { initialStep: 1e-3, parameterMax: 5000, maxSteps: 200_000, maxStep: 5 };
    const escape = { id: 'escape', value: (x: Vec4) => x[1] - 200, direction: 1 as const, terminal: true };

    for (const { ray, label } of rays()) {
      let earlyVerdict = false;
      integrateGeodesic({
        model,
        integrator: integrator(),
        initial: ray,
        limits,
        events: [escape],
        terminator: (x, k) => (earlyVerdict = isCaptured(model, x, k)),
      });

      let reachedHorizon = false;
      integrateGeodesic({
        model,
        integrator: integrator(),
        initial: ray,
        limits,
        events: [escape],
        terminator: (x) => (reachedHorizon = x[1] < horizon * 1.0005),
      });

      expect(earlyVerdict, label).toBe(reachedHorizon);
    }
  });

  it('saves the integrator from chasing a coordinate singularity', () => {
    // A captured ray stopped on the exact criterion takes far fewer steps than one
    // followed to within 5e-4 of r_+, where dphi/dlambda is diverging.
    const observer = zamoObserver(model, [0, 25, Math.PI / 2, 0]);
    const ray = generateNullRay(observer, [-1, 0, 0]);
    const limits = { initialStep: 1e-3, parameterMax: 5000, maxSteps: 400_000, maxStep: 5 };
    const early = integrateGeodesic({
      model,
      integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-11, relative: 1e-11 } }),
      initial: ray,
      limits,
      terminator: (x, k) => isCaptured(model, x, k),
    });
    const late = integrateGeodesic({
      model,
      integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-11, relative: 1e-11 } }),
      initial: ray,
      limits,
      terminator: (x) => x[1] < horizon * 1.0005,
    });
    expect(early.final.position_x[1]).toBeGreaterThan(late.final.position_x[1]);
    expect(early.steps).toBeLessThan(late.steps);
  });
});
