import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { nullState } from '../../src/physics/core/phase-space.js';
import { HAMILTONIAN } from '../../src/physics/geodesic/formulation.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { generateNullRay } from '../../src/physics/observer/observer.js';
import { orthonormalityResidual } from '../../src/physics/observer/tetrad.js';
import { zamoLapse, zamoObserver, zamoTetrad } from '../../src/physics/observer/zamo.js';
import {
  carterConstant,
  equatorialPhotonOrbitRadius,
  ergosphereRadius,
  frameDraggingOmega,
  kerr,
  outerHorizonRadius,
} from '../../src/physics/spacetimes/kerr.js';
import { schwarzschild } from '../../src/physics/spacetimes/schwarzschild.js';
import { normalizationResidual } from '../../src/physics/validation/normalization.js';

/**
 * Geodesics and observers in Kerr (ROADMAP.md 4A).
 *
 * Two things are worth separating here. E = -p_t and L_z = p_phi come from Killing
 * vectors, and the Hamiltonian formulation conserves them by construction, so their
 * conservation tests nothing about the integration. The Carter constant comes from a
 * Killing *tensor*: nothing in the update rule preserves it, so its drift is a genuine
 * measurement of the integration's quality, and it is the one conservation test in this
 * file that can actually fail.
 */

const M = 1;

function tight(): Dopri5Integrator {
  return new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-12, relative: 1e-12 } });
}

/**
 * The radial momentum that makes (p_t, p_r, p_theta, p_phi) null at this event.
 *
 * Solving g^{mu nu} p_mu p_nu = 0 for p_r, which is the only free component once E, L_z
 * and p_theta are chosen. The sign selects an ingoing or outgoing ray.
 */
function nullRadialMomentum(
  model: ReturnType<typeof kerr>,
  x: Vec4,
  p_t: number,
  p_theta: number,
  p_phi: number,
  sign: 1 | -1,
): number {
  const gInv = new Float64Array(16);
  model.inverseMetricInto(x, gInv);
  const rest =
    gInv[0] * p_t * p_t + 2 * gInv[3] * p_t * p_phi + gInv[15] * p_phi * p_phi + gInv[10] * p_theta * p_theta;
  const squared = -rest / gInv[5];
  if (!(squared >= 0)) {
    throw new RangeError(`nullRadialMomentum: no null ray with these constants at r = ${x[1]}.`);
  }
  return sign * Math.sqrt(squared);
}

/** Raise p_mu to k^mu with the model's inverse metric. */
function raise(model: ReturnType<typeof kerr>, x: Vec4, p: Vec4): Vec4 {
  const gInv = new Float64Array(16);
  model.inverseMetricInto(x, gInv);
  const k: [number, number, number, number] = [0, 0, 0, 0];
  for (let mu = 0; mu < 4; mu += 1) {
    let sum = 0;
    for (let nu = 0; nu < 4; nu += 1) sum += gInv[mu * 4 + nu] * p[nu];
    k[mu] = sum;
  }
  return k;
}

describe('the zero-angular-momentum observer', () => {
  it('is orthonormal, including inside the ergosphere where no static observer exists', () => {
    const a = 0.9;
    const model = kerr(M, a);
    const horizon = outerHorizonRadius(M, a);
    for (const [r, theta] of [[1.6, Math.PI / 2], [3, 0.7], [12, 1.4], [400, 2.2]] as const) {
      expect(r).toBeGreaterThan(horizon);
      const metric = model.metricAt([0, r, theta, 0]);
      const residual = orthonormalityResidual(metric, zamoTetrad(model, [0, r, theta, 0]));
      expect(residual, `r = ${r}, theta = ${theta}`).toBeLessThan(1e-14);
    }
    // r = 1.6 is inside the equatorial static limit r_E = 2M, so g_tt > 0 there: the
    // timelike Killing vector has turned spacelike and hovering is impossible.
    expect(ergosphereRadius(M, a, Math.PI / 2)).toBeCloseTo(2, 12);
    expect(model.metricAt([0, 1.6, Math.PI / 2, 0]).g_mu_nu[0]).toBeGreaterThan(0);
  });

  it('has exactly zero angular momentum, which is what names it', () => {
    const model = kerr(M, 0.7);
    for (const [r, theta] of [[2.5, Math.PI / 2], [6, 1.0], [80, 2.0]] as const) {
      const x: Vec4 = [0, r, theta, 0];
      const g = new Float64Array(16);
      model.metricInto(x, g);
      const u = zamoTetrad(model, x).four_velocity_u();
      let p_phi = 0;
      for (let mu = 0; mu < 4; mu += 1) p_phi += g[3 * 4 + mu] * u[mu];
      expect(Math.abs(p_phi), `r = ${r}`).toBeLessThan(1e-15);
      // And it is swept around all the same: omega is not zero.
      expect(u[3] / u[0]).toBeCloseTo(frameDraggingOmega(M, 0.7, r, theta), 14);
    }
  });

  it('has a lapse that falls to zero at the horizon and rises to one far away', () => {
    const a = 0.9;
    const model = kerr(M, a);
    const horizon = outerHorizonRadius(M, a);
    // alpha -> 1 - M/r far away, so (1 - alpha) r -> M rather than alpha reaching 1.
    for (const r of [1e4, 1e6, 1e8]) {
      expect(Math.abs((1 - zamoLapse(model, [0, r, 1.2, 0])) * r - M), `r = ${r}`).toBeLessThan(2e-3);
    }
    let previous = 1;
    for (const offset of [1, 1e-2, 1e-4, 1e-6]) {
      const alpha = zamoLapse(model, [0, horizon + offset, 1.2, 0]);
      expect(alpha).toBeLessThan(previous);
      previous = alpha;
    }
    expect(previous).toBeLessThan(1e-3);
  });

  it('generates genuinely null rays', () => {
    const model = kerr(M, 0.9);
    const observer = zamoObserver(model, [0, 10, Math.PI / 2.2, 0]);
    for (const raw of [
      [1, 0, 0],
      [-1, 0, 0],
      [0.3, -0.6, 0.74],
      [-0.5, 0.5, 0.707],
    ] as const) {
      const norm = Math.hypot(...raw);
      const ray = generateNullRay(observer, [raw[0] / norm, raw[1] / norm, raw[2] / norm]);
      expect(Math.abs(normalizationResidual(model, ray)), `direction ${raw.join(', ')}`).toBeLessThan(1e-14);
    }
  });
});

describe('conserved quantities along a Kerr geodesic', () => {
  it('holds E and L_z exactly, and the Carter constant to 1e-9 over a long trace', () => {
    const a = 0.9;
    const model = kerr(M, a);
    const session = HAMILTONIAN.bind(model);
    const observer = zamoObserver(model, [0, 15, 1.2, 0]);

    for (const direction of [
      [-1, 0.25, 0.4],
      [-0.9, -0.3, 0.32],
      [-0.7, 0.5, -0.51],
    ] as const) {
      const norm = Math.hypot(...direction);
      const ray = generateNullRay(observer, [direction[0] / norm, direction[1] / norm, direction[2] / norm]);
      const initialQ = carterConstant(M, a, ray.position_x, ray.tangent);

      const result = integrateGeodesic({
        model,
        integrator: tight(),
        initial: ray,
        session,
        limits: { initialStep: 1e-3, parameterMax: 400, maxSteps: 200_000, maxStep: 5 },
        events: [{ id: 'escape', value: (x) => x[1] - 300, direction: 1, terminal: true }],
        terminator: (x) => x[1] < outerHorizonRadius(M, a) * 1.001,
        recordPath: true,
      });

      const packedStart = session.pack(ray, new Float64Array(STATE_DIM));
      // Momenta conjugate to t and phi: conserved by construction in this formulation.
      expect(Math.abs(result.packed[4] - packedStart[4]), 'E').toBeLessThan(1e-15);
      expect(Math.abs(result.packed[7] - packedStart[7]), 'L_z').toBeLessThan(1e-15);

      const path = result.path;
      expect(path).toBeDefined();
      let worstQ = 0;
      for (const state of path ?? []) {
        const Q = carterConstant(M, a, state.position_x, state.tangent);
        worstQ = Math.max(worstQ, Math.abs(Q - initialQ) / Math.max(Math.abs(initialQ), 1e-3));
      }
      expect(worstQ, `Carter Q drift for direction ${direction.join(', ')}`).toBeLessThan(1e-9);
      expect(initialQ).toBeGreaterThan(0);
    }
  });

  it('keeps an equatorial ray equatorial, because its Carter constant is zero', () => {
    const a = 0.8;
    const model = kerr(M, a);
    const x: Vec4 = [0, 20, Math.PI / 2, 0];
    // E = 1, L_z = 5, p_theta = 0, and p_r fixed by the null condition: Q = 0, so there is
    // no polar motion to start and none can develop. Nothing in the integrator enforces it.
    const p_r = nullRadialMomentum(model, x, -1, 0, 5, -1);
    const k = raise(model, x, [-1, p_r, 0, 5]);
    const initial = nullState(x, k);
    expect(Math.abs(normalizationResidual(model, initial))).toBeLessThan(1e-13);
    expect(Math.abs(carterConstant(M, a, x, k))).toBeLessThan(1e-20);

    const result = integrateGeodesic({
      model,
      integrator: tight(),
      initial,
      limits: { initialStep: 1e-3, parameterMax: 500, maxSteps: 100_000, maxStep: 5 },
      events: [{ id: 'escape', value: (y) => y[1] - 200, direction: 1, terminal: true }],
      recordPath: true,
    });
    for (const state of result.path ?? []) {
      expect(Math.abs(state.position_x[2] - Math.PI / 2)).toBeLessThan(1e-13);
    }
  });

  it('drags a photon with no angular momentum around the hole', () => {
    // L_z = p_phi = 0, yet dphi/dlambda = g^{phi t} p_t is not zero: the photon is swept
    // in the direction of the spin. This is frame dragging with nothing else mixed in,
    // since the photon carries no angular momentum of its own to confuse it with.
    const a = 0.9;
    const model = kerr(M, a);
    const x: Vec4 = [0, 8, Math.PI / 2, 0];
    const k = raise(model, x, [-1, 0, 0, 0]);
    expect(k[3]).toBeGreaterThan(0); // dragged toward +phi, the direction of the spin

    const result = integrateGeodesic({
      model,
      integrator: tight(),
      initial: nullState(x, k),
      limits: { initialStep: 1e-4, parameterMax: 200, maxSteps: 200_000 },
      terminator: (y) => y[1] < outerHorizonRadius(M, a) * 1.0005,
    });
    // It falls in — a radially infalling photon must — having been carried around in phi.
    expect(result.final.position_x[1]).toBeLessThan(2);
    expect(result.final.position_x[3]).toBeGreaterThan(0.5);

    // With no spin the same initial data falls straight in with no change in phi at all.
    const still = kerr(M, 0);
    const kStill = raise(still, x, [-1, 0, 0, 0]);
    expect(Math.abs(kStill[3])).toBeLessThan(1e-18);
  });
});

describe('circular photon orbits', () => {
  it('locks to the prograde and retrograde equatorial photon orbits', () => {
    const a = 0.9;
    const model = kerr(M, a);
    for (const sense of ['prograde', 'retrograde'] as const) {
      const r = equatorialPhotonOrbitRadius(M, a, sense);
      // The spherical-photon-orbit constants at this radius, from R(r) = R'(r) = 0:
      //   xi = [M(r^2 - a^2) - r Delta] / [a (r - M)],  eta = 0 in the equatorial plane.
      const delta = model.deltaFunction(r);
      const xi = (M * (r * r - a * a) - r * delta) / (a * (r - M));
      const x: Vec4 = [0, r, Math.PI / 2, 0];
      const k = raise(model, x, [-1, 0, 0, xi]);
      const initial = nullState(x, k);
      expect(Math.abs(normalizationResidual(model, initial)), `${sense} null`).toBeLessThan(1e-13);

      // Two turns, no more. A circular photon orbit is unstable — that is what makes it
      // the edge of the shadow — so the rounding in the initial data grows exponentially
      // along it, at the orbit's Lyapunov rate. Holding it for two turns measures the
      // integration; asking for twenty would measure the instability instead.
      const result = integrateGeodesic({
        model,
        integrator: tight(),
        initial,
        limits: { initialStep: 1e-4, parameterMax: 4000, maxSteps: 400_000, maxStep: 1 },
        events: [{ id: 'two-turns', value: (y) => Math.abs(y[3]) - 4 * Math.PI, direction: 1, terminal: true }],
        recordPath: true,
      });

      let worst = 0;
      for (const state of result.path ?? []) worst = Math.max(worst, Math.abs(state.position_x[1] - r));
      const sweptTurns = Math.abs(result.final.position_x[3]) / (2 * Math.PI);
      expect(sweptTurns, `${sense} turns completed`).toBeGreaterThan(1.9);
      expect(worst, `${sense} radial wander over ${sweptTurns.toFixed(1)} turns`).toBeLessThan(1e-8);
    }
  });

  it('puts the prograde orbit inside the retrograde one, both closer in than Schwarzschild', () => {
    const a = 0.9;
    const prograde = equatorialPhotonOrbitRadius(M, a, 'prograde');
    const retrograde = equatorialPhotonOrbitRadius(M, a, 'retrograde');
    expect(prograde).toBeLessThan(3);
    expect(retrograde).toBeGreaterThan(3);
    expect(prograde).toBeLessThan(retrograde);
  });
});

describe('the Schwarzschild limit of a traced geodesic', () => {
  it('reproduces a Schwarzschild ray as a -> 0', () => {
    const x: Vec4 = [0, 25, Math.PI / 2, 0];
    const trace = (model: ReturnType<typeof kerr> | ReturnType<typeof schwarzschild>): number => {
      const gInv = new Float64Array(16);
      model.inverseMetricInto(x, gInv);
      const p: Vec4 = [-1, 0, 0, 6];
      const k: [number, number, number, number] = [0, 0, 0, 0];
      for (let mu = 0; mu < 4; mu += 1) {
        let sum = 0;
        for (let nu = 0; nu < 4; nu += 1) sum += gInv[mu * 4 + nu] * p[nu];
        k[mu] = sum;
      }
      const result = integrateGeodesic({
        model,
        integrator: tight(),
        initial: nullState(x, k),
        limits: { initialStep: 1e-3, parameterMax: 2000, maxSteps: 200_000, maxStep: 5 },
        events: [{ id: 'escape', value: (y) => y[1] - 400, direction: 1, terminal: true }],
      });
      return result.final.position_x[3];
    };

    const reference = trace(schwarzschild(M));
    let previous = Number.POSITIVE_INFINITY;
    for (const a of [1e-2, 1e-4, 1e-6]) {
      const miss = Math.abs(trace(kerr(M, a)) - reference);
      expect(miss, `a = ${a}`).toBeLessThan(previous);
      previous = miss;
    }
    expect(previous).toBeLessThan(1e-6);
  });
});
