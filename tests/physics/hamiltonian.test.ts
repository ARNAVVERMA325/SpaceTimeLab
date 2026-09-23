import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { nullState, timelikeState } from '../../src/physics/core/phase-space.js';
import { HAMILTONIAN, LAGRANGIAN } from '../../src/physics/geodesic/formulation.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { inverseResidual } from '../../src/physics/geometry/invert-metric.js';
import { minkowski } from '../../src/physics/spacetimes/minkowski.js';
import { schwarzschild } from '../../src/physics/spacetimes/schwarzschild.js';
import {
  deflectionAngleExact,
  impactParameterForTurningPoint,
} from '../../src/physics/spacetimes/schwarzschild-analytic.js';
import { equatorialNullRay, measureDeflection } from '../../src/physics/spacetimes/schwarzschild-rays.js';
import type { SpacetimeModel } from '../../src/physics/spacetimes/spacetime-model.js';
import { normalizationResidual } from '../../src/physics/validation/normalization.js';

/**
 * The Hamiltonian formulation (CLAUDE.md §2: "or preferably, where appropriate, in
 * Hamiltonian first-order form").
 *
 * Symbolic check, recorded here because the test suite cannot run a CAS: feeding the
 * Schwarzschild d_alpha g^{mu nu} through Hamilton's equations and subtracting
 * d/dlambda (g_{mu nu} k^nu) along the Christoffel flow gives identically zero for every
 * component (sympy 1.14). The two formulations integrate the same curves; these tests
 * establish that numerically, and then measure what differs between them.
 */

const schw = schwarzschild(1);
const MODELS: readonly SpacetimeModel[] = [minkowski, schw];

const EVENTS: Readonly<Record<string, readonly Vec4[]>> = {
  minkowski: [
    [0, 0, 0, 0],
    [3, -2, 7, 1],
  ],
  schwarzschild: [
    [0, 3, Math.PI / 2, 0],
    [0, 7.5, Math.PI / 3, 1],
    [5, 40, 0.4, 2],
    [0, 2.2, 2.5, 0],
  ],
};

function tight(): Dopri5Integrator {
  return new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-12, relative: 1e-12 } });
}

describe('model inputs to the Hamiltonian', () => {
  for (const model of MODELS) {
    it(`${model.id}: g^{mu nu} from inverseMetricInto inverts g_{mu nu} from metricInto`, () => {
      const g = new Float64Array(16);
      const gInv = new Float64Array(16);
      for (const x of EVENTS[model.id]) {
        model.metricInto(x, g);
        model.inverseMetricInto(x, gInv);
        expect(inverseResidual(g, gInv), `at ${x}`).toBeLessThan(1e-13);
        // And agree with the allocating accessor.
        const reference = model.metricAt(x);
        for (let i = 0; i < 16; i += 1) {
          expect(g[i]).toBe(reference.g_mu_nu[i]);
          expect(gInv[i]).toBe(reference.g_inv_mu_nu[i]);
        }
      }
    });

    it(`${model.id}: analytic d_alpha g^{mu nu} match central differences`, () => {
      const d = new Float64Array(64);
      const plus = new Float64Array(16);
      const minus = new Float64Array(16);
      for (const x of EVENTS[model.id]) {
        model.inverseMetricDerivativesInto(x, d);
        for (let alpha = 0; alpha < 4; alpha += 1) {
          const h = 1e-5 * Math.max(1, Math.abs(x[alpha]));
          const xp = [...x] as [number, number, number, number];
          const xm = [...x] as [number, number, number, number];
          xp[alpha] += h;
          xm[alpha] -= h;
          model.inverseMetricInto(xp, plus);
          model.inverseMetricInto(xm, minus);
          for (let i = 0; i < 16; i += 1) {
            const numeric = (plus[i] - minus[i]) / (xp[alpha] - xm[alpha]);
            const analytic = d[alpha * 16 + i];
            const scale = Math.max(1, Math.abs(analytic));
            expect(
              Math.abs(numeric - analytic) / scale,
              `d_${alpha} g^[${i}] at ${x}: analytic ${analytic}, numeric ${numeric}`,
            ).toBeLessThan(1e-7);
          }
        }
      }
    });
  }
});

describe('Hamiltonian and Lagrangian formulations integrate the same curves', () => {
  it('round-trips a state through pack and unpack', () => {
    const session = HAMILTONIAN.bind(schw);
    const state = nullState([0, 12, 1.1, 0.3], [1.2, -0.4, 0.01, 0.02]);
    const packed = session.pack(state, new Float64Array(STATE_DIM));
    const back = session.unpack(packed, 'null', 0);
    for (let mu = 0; mu < 4; mu += 1) {
      expect(back.position_x[mu]).toBe(state.position_x[mu]);
      expect(Math.abs(back.tangent[mu] - state.tangent[mu])).toBeLessThan(1e-15 * Math.max(1, Math.abs(state.tangent[mu])));
    }
  });

  it('agrees on a strongly bent Schwarzschild ray to integrator tolerance', () => {
    const initial = equatorialNullRay(schw, 50, 7, 'ingoing');
    const run = (formulation: typeof HAMILTONIAN) =>
      integrateGeodesic({
        model: schw,
        integrator: tight(),
        initial,
        session: formulation.bind(schw),
        limits: { initialStep: 1e-3, parameterMax: 90, maxSteps: 200_000 },
      });
    const h = run(HAMILTONIAN);
    const l = run(LAGRANGIAN);
    expect(h.reason).toBe('parameter-limit');
    expect(l.reason).toBe('parameter-limit');
    for (let mu = 0; mu < 4; mu += 1) {
      const scale = Math.max(1, Math.abs(h.final.position_x[mu]));
      expect(Math.abs(h.final.position_x[mu] - l.final.position_x[mu]) / scale).toBeLessThan(1e-9);
      expect(Math.abs(h.final.tangent[mu] - l.final.tangent[mu])).toBeLessThan(1e-9);
    }
  });

  it('agrees on a timelike worldline', () => {
    const r = 12;
    const f = 1 - 2 / r;
    // A bound, non-circular orbit: E and L chosen below the circular values at r = 12.
    const E = 0.97;
    const L = 4.1;
    const kr2 = E * E - f * (1 + (L * L) / (r * r));
    const initial = timelikeState([0, r, Math.PI / 2, 0], [E / f, Math.sqrt(kr2), 0, L / (r * r)]);
    expect(Math.abs(normalizationResidual(schw, initial))).toBeLessThan(1e-14);

    const run = (formulation: typeof HAMILTONIAN) =>
      integrateGeodesic({
        model: schw,
        integrator: tight(),
        initial,
        session: formulation.bind(schw),
        limits: { initialStep: 1e-3, parameterMax: 400, maxSteps: 400_000 },
      });
    const h = run(HAMILTONIAN);
    const l = run(LAGRANGIAN);
    expect(Math.abs(h.final.position_x[1] - l.final.position_x[1])).toBeLessThan(1e-7);
    expect(Math.abs(h.final.position_x[3] - l.final.position_x[3])).toBeLessThan(1e-7);
  });
});

describe('what the Hamiltonian formulation buys', () => {
  it('conserves p_t and p_phi bit for bit, with any integrator', () => {
    // The metric does not depend on t or phi, so d_t g and d_phi g are identically zero,
    // dp_t/dlambda and dp_phi/dlambda are identically zero, and every Runge-Kutta stage
    // adds exactly nothing. This is exact conservation, not small drift.
    const initial = equatorialNullRay(schw, 200, 6, 'ingoing');
    const session = HAMILTONIAN.bind(schw);
    const before = session.pack(initial, new Float64Array(STATE_DIM));
    const result = integrateGeodesic({
      model: schw,
      integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-8, relative: 1e-8 } }),
      initial,
      session,
      limits: { initialStep: 1e-3, parameterMax: 600, maxSteps: 200_000 },
    });
    expect(result.steps).toBeGreaterThan(50);
    expect(result.packed[4]).toBe(before[4]); // p_t, so energy_E = -p_t
    expect(result.packed[7]).toBe(before[7]); // p_phi = angular_momentum_Lz
  });

  it('is more accurate than the Christoffel form at equal tolerance, for fewer steps', () => {
    // Measured at tolerance 1e-10 on the r_0 = 10M deflection: Hamiltonian 1.4e-9 relative
    // error in 194 steps, Lagrangian 1.7e-8 in 274. Asserted with margin.
    const b = impactParameterForTurningPoint(1, 10);
    const exact = deflectionAngleExact(1, 10);
    const run = (formulation: typeof HAMILTONIAN) =>
      measureDeflection({
        model: schw,
        integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-10, relative: 1e-10 } }),
        impactParameter: b,
        startRadius: 1000,
        formulation,
      });
    const h = run(HAMILTONIAN);
    const l = run(LAGRANGIAN);
    const errorH = Math.abs(h.deflectionAngle! - exact) / exact;
    const errorL = Math.abs(l.deflectionAngle! - exact) / exact;
    expect(errorH).toBeLessThan(errorL / 4);
    expect(h.steps).toBeLessThan(l.steps);
    expect(h.angularMomentumDrift).toBe(0);
    expect(l.angularMomentumDrift).toBeGreaterThan(1e-10);
  });

  it('leaves the mass shell as the one invariant that drifts, and reports it', () => {
    // The trade, stated rather than hidden: with p_t and p_phi pinned, the null
    // constraint H = 0 is what absorbs the discretization error. Measured at 4.4e-10 for
    // the Hamiltonian against 1.3e-10 for the Lagrangian form on the same trace.
    const b = impactParameterForTurningPoint(1, 10);
    const m = measureDeflection({
      model: schw,
      integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-10, relative: 1e-10 } }),
      impactParameter: b,
      startRadius: 1000,
      formulation: HAMILTONIAN,
    });
    expect(m.nullResidual).toBeGreaterThan(0);
    expect(m.nullResidual).toBeLessThan(1e-9);
  });
});
