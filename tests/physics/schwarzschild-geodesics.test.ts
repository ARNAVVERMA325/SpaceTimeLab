import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { nullState } from '../../src/physics/core/phase-space.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import type { ErrorTolerance } from '../../src/physics/geodesic/integrators/integrator.js';
import { RKF45Integrator } from '../../src/physics/geodesic/integrators/rkf45.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import {
  criticalImpactParameter,
  photonSphereRadius,
  schwarzschild,
} from '../../src/physics/spacetimes/schwarzschild.js';
import {
  deflectionAngleExact,
  deflectionAngleWeakField,
  impactParameterForTurningPoint,
} from '../../src/physics/spacetimes/schwarzschild-analytic.js';
import {
  equatorialNullRay,
  isCaptured,
  measureDeflection,
} from '../../src/physics/spacetimes/schwarzschild-rays.js';
import { conservedQuantityDrift } from '../../src/physics/validation/conserved.js';
import { normalizationResidual } from '../../src/physics/validation/normalization.js';
import {
  CONSERVED_QUANTITY_DRIFT_SCHWARZSCHILD,
  CRITICAL_IMPACT_PARAMETER,
  DEFLECTION_VS_EXACT,
  NULL_NORMALIZATION_TRACED,
  PHOTON_SPHERE_LOCK,
  checkTolerance,
} from '../../src/physics/validation/tolerances.js';

/**
 * ROADMAP.md 2A.4 — the Milestone 2A physics validation gate.
 *
 * Three checks are named by the roadmap: the photon sphere locking at r = 3M, the
 * conserved quantities E and L_z holding to better than 1e-6 relative, and the
 * weak-field deflection matching alpha ~ 4M/b. All three are here, together with a
 * comparison against the exact orbit-equation quadrature that is stronger than any of
 * them individually.
 */

const M = 1;
const model = schwarzschild(M);
const TIGHT: ErrorTolerance = { absolute: 1e-12, relative: 1e-12 };

function integrator(tolerance: ErrorTolerance = TIGHT): RKF45Integrator {
  return new RKF45Integrator(STATE_DIM, { tolerance });
}

describe('photon sphere (ROADMAP.md 2A.4)', () => {
  const rPhoton = photonSphereRadius(M);

  it('puts the null effective potential maximum at r = 3M', () => {
    // V(r) = f / r^2 for a null geodesic. V'(r) = (6M - 2r)/r^4 vanishes at r = 3M,
    // and the second derivative is negative there, so the circular orbit is unstable.
    const V = (r: number): number => (1 - (2 * M) / r) / (r * r);
    expect(V(rPhoton)).toBeCloseTo(1 / 27, 15);
    expect(V(rPhoton)).toBeGreaterThan(V(rPhoton * 0.99));
    expect(V(rPhoton)).toBeGreaterThan(V(rPhoton * 1.01));

    // b_c = 1/sqrt(V(3M)) = 3 sqrt(3) M.
    expect(1 / Math.sqrt(V(rPhoton))).toBeCloseTo(criticalImpactParameter(M), 13);
  });

  it('locks a circular null orbit at r = 3M over several complete orbits', () => {
    const f = model.lapseFunction(rPhoton);
    const angular_momentum_Lz = criticalImpactParameter(M);
    // k^r = k^theta = 0; E = 1 fixes k^t, and L = b_c fixes k^phi.
    const k: Vec4 = [1 / f, 0, 0, angular_momentum_Lz / (rPhoton * rPhoton)];
    const initial = nullState([0, rPhoton, Math.PI / 2, 0], k);

    // The initial data really is null, before anything is integrated.
    expect(Math.abs(normalizationResidual(model, initial))).toBeLessThan(1e-15);

    const result = integrateGeodesic({
      model,
      integrator: integrator({ absolute: 1e-14, relative: 1e-14 }),
      initial,
      limits: { initialStep: 1e-4, parameterMax: 50, maxSteps: 2_000_000, maxStep: 0.05 },
    });

    expect(result.reason).toBe('parameter-limit');
    const check = checkTolerance(PHOTON_SPHERE_LOCK, result.final.position_x[1] - rPhoton);
    expect(check.withinTolerance, check.message).toBe(true);

    // It must genuinely have gone round, not merely failed to move.
    const turns = result.final.position_x[3] / (2 * Math.PI);
    expect(turns).toBeGreaterThan(4);
  });

  it('shows the orbit is unstable, as the potential maximum requires', () => {
    // Guards the lock test above: an orbit that stayed at 3M forever would mean the
    // integrator was not following the physics, since this equilibrium is a maximum.
    const f = model.lapseFunction(rPhoton);
    const k: Vec4 = [1 / f, 0, 0, criticalImpactParameter(M) / (rPhoton * rPhoton)];

    const unperturbed = integrateGeodesic({
      model,
      integrator: integrator({ absolute: 1e-14, relative: 1e-14 }),
      initial: nullState([0, rPhoton, Math.PI / 2, 0], k),
      limits: { initialStep: 1e-4, parameterMax: 150, maxSteps: 2_000_000, maxStep: 0.05 },
    });
    // By lambda = 150M, rounding in the initial data alone has grown to order unity.
    expect(Math.abs(unperturbed.final.position_x[1] - rPhoton)).toBeGreaterThan(1);

    // An explicit outward perturbation escapes outward.
    const perturbed = integrateGeodesic({
      model,
      integrator: integrator({ absolute: 1e-14, relative: 1e-14 }),
      initial: nullState([0, rPhoton * (1 + 1e-8), Math.PI / 2, 0], k),
      limits: { initialStep: 1e-4, parameterMax: 300, maxSteps: 2_000_000, maxStep: 0.05 },
    });
    expect(perturbed.final.position_x[1]).toBeGreaterThan(100);
  });

  it('puts the capture threshold at b_c = 3 sqrt(3) M', () => {
    // The strongest check in the milestone: it exercises the metric, the connection, the
    // integrator and the capture condition at once, and compares against a closed-form
    // constant that none of them knows about.
    const captured = (b: number): boolean =>
      measureDeflection({
        model,
        integrator: integrator(),
        impactParameter: b,
        startRadius: 1000,
      }).captured;

    let lo = 5.0;
    let hi = 5.5;
    expect(captured(lo)).toBe(true);
    expect(captured(hi)).toBe(false);

    for (let i = 0; i < 45 && hi - lo > 1e-13; i += 1) {
      const mid = 0.5 * (lo + hi);
      if (captured(mid)) lo = mid;
      else hi = mid;
    }

    const bisected = 0.5 * (lo + hi);
    const exact = criticalImpactParameter(M);
    const check = checkTolerance(CRITICAL_IMPACT_PARAMETER, bisected - exact, exact);
    expect(check.withinTolerance, `${check.message} (bisected ${bisected}, exact ${exact})`).toBe(
      true,
    );
  });

  it('applies the capture condition only where it is exact', () => {
    // Inward-moving below the photon sphere: captured, with no escape possible.
    expect(isCaptured(model, [0, 2.5, Math.PI / 2, 0], [1, -1, 0, 0])).toBe(true);
    // Outward-moving below the photon sphere: undecided, so not declared captured.
    expect(isCaptured(model, [0, 2.5, Math.PI / 2, 0], [1, 1, 0, 0])).toBe(false);
    // Above the photon sphere: never captured by this test, whichever way it moves.
    expect(isCaptured(model, [0, 5, Math.PI / 2, 0], [1, -1, 0, 0])).toBe(false);
  });
});

describe('conserved quantities (ROADMAP.md 2A.4, CLAUDE.md §16)', () => {
  it('declares only the symmetries Schwarzschild actually has', () => {
    const ids = model.killingVectors.map((k) => k.id);
    expect(ids).toEqual(['d_dt', 'd_dphi']);
    const names = model.killingVectors.map((k) => k.conservedQuantityName);
    expect(names).toEqual(['energy_E', 'angular_momentum_Lz']);
  });

  it('holds E and L_z along traced geodesics, well inside the roadmap budget', () => {
    for (const r0 of [3.2, 4, 6, 10, 100, 1000]) {
      const b = impactParameterForTurningPoint(M, r0);
      const measurement = measureDeflection({
        model,
        integrator: integrator(),
        impactParameter: b,
        startRadius: Math.max(1000, 20 * r0),
      });
      expect(measurement.captured, `b = ${b} should escape`).toBe(false);

      const energy = checkTolerance(
        CONSERVED_QUANTITY_DRIFT_SCHWARZSCHILD,
        measurement.energyDrift,
        1,
      );
      expect(energy.withinTolerance, `energy_E at r_0 = ${r0}M: ${energy.message}`).toBe(true);

      const angular = checkTolerance(
        CONSERVED_QUANTITY_DRIFT_SCHWARZSCHILD,
        measurement.angularMomentumDrift,
        1,
      );
      expect(
        angular.withinTolerance,
        `angular_momentum_Lz at r_0 = ${r0}M: ${angular.message}`,
      ).toBe(true);

      // ROADMAP.md 2A.4's own figure, asserted explicitly so the stated gate is visibly met.
      expect(measurement.energyDrift).toBeLessThan(1e-6);
      expect(measurement.angularMomentumDrift).toBeLessThan(1e-6);

      const nullCheck = checkTolerance(NULL_NORMALIZATION_TRACED, measurement.nullResidual);
      expect(nullCheck.withinTolerance, `r_0 = ${r0}M: ${nullCheck.message}`).toBe(true);
    }
  });

  it('agrees with the generic Killing-vector machinery', () => {
    // The deflection helper reads E and L_z straight off the tangent. The validation
    // layer computes them from the Killing vectors and the metric. They must agree.
    const initial = equatorialNullRay(model, 100, 20, 'ingoing');
    const result = integrateGeodesic({
      model,
      integrator: integrator(),
      initial,
      limits: { initialStep: 1e-3, parameterMax: 50, maxSteps: 200_000 },
    });

    const drifts = conservedQuantityDrift(model, initial, result.final);
    const energy = drifts.find((d) => d.name === 'energy_E');
    const angular = drifts.find((d) => d.name === 'angular_momentum_Lz');

    expect(energy?.initial).toBeCloseTo(1, 12);
    expect(angular?.initial).toBeCloseTo(20, 10);
    expect(energy!.relativeDrift).toBeLessThan(1e-7);
    expect(angular!.relativeDrift).toBeLessThan(1e-7);
  });
});

describe('deflection angle (ROADMAP.md 2A.4)', () => {
  it('matches the exact orbit-equation quadrature from weak to strong field', () => {
    for (const r0 of [3.2, 4, 6, 10, 20, 100, 1000, 1e4]) {
      const b = impactParameterForTurningPoint(M, r0);
      const expected = deflectionAngleExact(M, r0);
      const measurement = measureDeflection({
        model,
        integrator: integrator(),
        impactParameter: b,
        startRadius: Math.max(1000, 20 * r0),
      });

      expect(measurement.deflectionAngle, `r_0 = ${r0}M was captured`).toBeDefined();
      const check = checkTolerance(
        DEFLECTION_VS_EXACT,
        measurement.deflectionAngle! - expected,
        expected,
      );
      expect(
        check.withinTolerance,
        `r_0 = ${r0}M: traced ${measurement.deflectionAngle}, exact ${expected}. ${check.message}`,
      ).toBe(true);
    }
  });

  it('converges on the weak-field limit 4M/b as the impact parameter grows', () => {
    // The roadmap's named benchmark. Stated as convergence rather than equality: 4M/b is
    // the leading term of an expansion, and the correction falls off as M/b.
    const ratios: number[] = [];
    for (const r0 of [100, 1000, 1e4, 1e5]) {
      const b = impactParameterForTurningPoint(M, r0);
      const measurement = measureDeflection({
        model,
        integrator: integrator(),
        impactParameter: b,
        startRadius: Math.max(1000, 20 * r0),
      });
      ratios.push(measurement.deflectionAngle! / deflectionAngleWeakField(M, b));
    }

    // Each step of ten in b cuts the discrepancy by about ten: the correction is O(M/b).
    for (let i = 0; i < ratios.length; i += 1) {
      expect(ratios[i], `ratio ${ratios[i]} should exceed 1`).toBeGreaterThan(1);
      if (i > 0) {
        const previousExcess = ratios[i - 1] - 1;
        const excess = ratios[i] - 1;
        expect(excess).toBeLessThan(previousExcess);
        expect(previousExcess / excess).toBeGreaterThan(5);
        expect(previousExcess / excess).toBeLessThan(20);
      }
    }
    // At b ~ 10^5 M the weak-field formula is good to a part in 10^4.
    expect(Math.abs(ratios[ratios.length - 1] - 1)).toBeLessThan(1e-4);
  });

  it('shows the weak-field formula failing in the strong field', () => {
    // Without this, the benchmark above could be satisfied by a renderer that bent
    // light by 4M/b everywhere and got the strong field badly wrong.
    const b = impactParameterForTurningPoint(M, 4);
    const measurement = measureDeflection({
      model,
      integrator: integrator(),
      impactParameter: b,
      startRadius: 1000,
    });
    expect(measurement.deflectionAngle!).toBeGreaterThan(2);
    expect(deflectionAngleWeakField(M, b)).toBeLessThan(0.75);
  });

  it('bends light by more than a half turn just outside the photon sphere', () => {
    const b = impactParameterForTurningPoint(M, 3.1);
    const measurement = measureDeflection({
      model,
      integrator: integrator(),
      impactParameter: b,
      startRadius: 1000,
    });
    expect(measurement.captured).toBe(false);
    expect(measurement.deflectionAngle!).toBeGreaterThan(Math.PI);
    expect(measurement.minimumRadius).toBeGreaterThan(photonSphereRadius(M));
  });

  it('captures rays inside the critical impact parameter', () => {
    for (const b of [0.5, 2, 4, 5.19]) {
      const measurement = measureDeflection({
        model,
        integrator: integrator(),
        impactParameter: b,
        startRadius: 1000,
      });
      expect(measurement.captured, `b = ${b}M should be captured`).toBe(true);
      expect(measurement.deflectionAngle).toBeUndefined();
    }
  });

  it('refuses a ray set up inside its own turning point', () => {
    // b = 100 has a turning point near r = 100, so no such geodesic passes through r = 10.
    expect(() => equatorialNullRay(model, 10, 100, 'ingoing')).toThrow(RangeError);
  });
});
