import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { nullState } from '../../src/physics/core/phase-space.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import type { ErrorTolerance } from '../../src/physics/geodesic/integrators/integrator.js';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
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

  it('grows perturbations at the analytic Lyapunov exponent 1/(3 sqrt(3) M)', () => {
    // The instability, measured rather than merely observed. Linearizing
    // d^2r/dlambda^2 = -1/2 L^2 V'(r) about r = 3M gives d^2 delta/dlambda^2 = kappa^2 delta
    // with kappa = E / (sqrt(3) M) per unit affine parameter; with dt/dlambda = 3E at the
    // photon sphere that is lambda_t = 1/(3 sqrt(3) M) in coordinate time — the result of
    // Cardoso et al., Phys. Rev. D 79, 064016 (2009), and equal to the orbital frequency,
    // so a perturbation grows by exactly e^pi per half orbit.
    //
    // An earlier version of this test only asserted that an *unperturbed* orbit had
    // wandered off by lambda = 150M. That was testing an artifact: in the Christoffel form,
    // rounding seeded the instability, while in the Hamiltonian form dp_r/dlambda at r = 3M
    // evaluates to exactly zero and the circular orbit is a bit-exact fixed point — which
    // is the correct solution. The physics is in the growth rate, so that is what is
    // measured.
    const b = criticalImpactParameter(M);
    const epsilon = 3e-9;
    const r = rPhoton + epsilon;
    const f = model.lapseFunction(r);
    // The pure growing mode. With E = 1 and b = b_c the radial null condition factors as
    // (k^r)^2 = (9 eps^2 + eps^3) / r^3, which avoids the catastrophic cancellation in
    // E^2 - f L^2 / r^2: that difference of two O(1) numbers is about 3e-18 here, far
    // below the rounding of either term.
    const kr = Math.sqrt((9 * epsilon * epsilon + epsilon ** 3) / r ** 3);
    const initial = nullState([0, r, Math.PI / 2, 0], [1 / f, kr, 0, b / (r * r)]);
    expect(Math.abs(normalizationResidual(model, initial))).toBeLessThan(1e-15);

    const result = integrateGeodesic({
      model,
      integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-13, relative: 1e-13 } }),
      initial,
      recordPath: true,
      limits: { initialStep: 1e-3, parameterMax: 40, maxSteps: 1_000_000, maxStep: 0.05 },
    });

    // Least-squares slope of ln(r - 3M) against lambda, in the linear regime.
    const points = result.path!
      .map((state) => [state.parameter, Math.log(state.position_x[1] - rPhoton)] as const)
      .filter(([, logDelta]) => logDelta > Math.log(3e-7) && logDelta < Math.log(3e-4));
    expect(points.length).toBeGreaterThan(100);
    const meanX = points.reduce((sum, [x]) => sum + x, 0) / points.length;
    const meanY = points.reduce((sum, [, y]) => sum + y, 0) / points.length;
    let sxx = 0;
    let sxy = 0;
    for (const [x, y] of points) {
      sxx += (x - meanX) ** 2;
      sxy += (x - meanX) * (y - meanY);
    }
    const measuredKappa = sxy / sxx;
    const analyticKappa = 1 / (Math.sqrt(3) * M);

    // Measured 1.2e-5, independent of the integration tolerance and of epsilon: the
    // residual is the O(delta) nonlinear term inside the fit window, not numerical error.
    expect(Math.abs(measuredKappa / analyticKappa - 1)).toBeLessThan(1e-4);

    // And the growth per half orbit is e^pi: kappa * (half-period in lambda).
    const halfPeriodAffine = Math.PI / (criticalImpactParameter(M) / (rPhoton * rPhoton));
    expect(measuredKappa * halfPeriodAffine).toBeCloseTo(Math.PI, 3);
  });

  it('holds the unperturbed circular orbit exactly in the Hamiltonian formulation', () => {
    const f = model.lapseFunction(rPhoton);
    const k: Vec4 = [1 / f, 0, 0, criticalImpactParameter(M) / (rPhoton * rPhoton)];
    const result = integrateGeodesic({
      model,
      integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-12, relative: 1e-12 } }),
      initial: nullState([0, rPhoton, Math.PI / 2, 0], k),
      limits: { initialStep: 1e-3, parameterMax: 300, maxSteps: 1_000_000, maxStep: 0.1 },
    });
    // Over some 27 orbits the radius stays at 3M: the equilibrium is exact, and nothing
    // in the Hamiltonian discretization perturbs it.
    const check = checkTolerance(PHOTON_SPHERE_LOCK, result.final.position_x[1] - rPhoton);
    expect(check.withinTolerance, check.message).toBe(true);
    expect(result.final.position_x[3] / (2 * Math.PI)).toBeGreaterThan(25);
  });

  it('lets an explicit outward perturbation escape to large radius', () => {
    const f = model.lapseFunction(rPhoton);
    const k: Vec4 = [1 / f, 0, 0, criticalImpactParameter(M) / (rPhoton * rPhoton)];
    const perturbed = integrateGeodesic({
      model,
      integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-12, relative: 1e-12 } }),
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
