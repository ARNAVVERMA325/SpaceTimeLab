import { describe, expect, it } from 'vitest';
import { timelikeState } from '../../src/physics/core/phase-space.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { GaussLegendre4Integrator } from '../../src/physics/geodesic/integrators/gauss-legendre.js';
import { RK4Integrator } from '../../src/physics/geodesic/integrators/rk4.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { schwarzschild } from '../../src/physics/spacetimes/schwarzschild.js';
import {
  circularOrbit,
  eccentricOrbitConstants,
  ellipticK,
  equatorialTimelikeState,
  iscoRadius,
  periapsisPrecessionExact,
  periapsisPrecessionWeakField,
  radialEpicyclicFrequency,
  radialInfallCoordinateTime,
  radialInfallProperTime,
} from '../../src/physics/spacetimes/schwarzschild-orbits.js';
import { normalizationResidual } from '../../src/physics/validation/normalization.js';

/**
 * Timelike Schwarzschild geodesics (CLAUDE.md §16: "Schwarzschild ISCO").
 *
 * Reference values for the precession were computed independently in mpmath three ways —
 * elliptic closed form, quadrature in the relativistic anomaly, direct quadrature in r —
 * which agreed to 22 digits.
 */

const M = 1;
const model = schwarzschild(M);

function tight(tolerance = 1e-13): Dopri5Integrator {
  return new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: tolerance, relative: tolerance } });
}

/** Periapsis advance per radial period, from mpmath (22 digits). */
const PRECESSION: readonly (readonly [number, number, number])[] = [
  [10, 0.1, 3.656066174600919439965],
  [20, 0.5, 1.233861806265436012065],
  [50, 0.3, 0.4149462859202380195607],
  [100, 0.2, 0.1974475160236572840877],
  [1000, 0.5, 0.01893600123348775831921],
  [12, 0.6, 2.670722986569505392589],
];

describe('closed forms', () => {
  it('computes K(m) by the arithmetic-geometric mean', () => {
    expect(ellipticK(0)).toBeCloseTo(Math.PI / 2, 15);
    // K(1/2) = Gamma(1/4)^2 / (4 sqrt(pi)).
    expect(ellipticK(0.5)).toBeCloseTo(1.8540746773013719, 14);
    expect(() => ellipticK(1)).toThrow(RangeError);
  });

  it('reproduces the mpmath periapsis advance', () => {
    for (const [p, e, expected] of PRECESSION) {
      expect(Math.abs(periapsisPrecessionExact(M, p, e) - expected) / expected, `p=${p} e=${e}`).toBeLessThan(1e-13);
    }
  });

  it('puts the ISCO at 6M with E = sqrt(8/9) and L = 2 sqrt(3) M', () => {
    expect(iscoRadius(M)).toBe(6);
    const isco = circularOrbit(M, 6);
    expect(isco.energy_E).toBeCloseTo(Math.sqrt(8 / 9), 15);
    expect(isco.angular_momentum_Lz).toBeCloseTo(2 * Math.sqrt(3), 14);
    // E and L are minimized along the circular family at the ISCO.
    for (const r of [5.9, 6.1]) {
      expect(circularOrbit(M, r).energy_E).toBeGreaterThan(isco.energy_E);
      expect(circularOrbit(M, r).angular_momentum_Lz).toBeGreaterThan(isco.angular_momentum_Lz);
    }
  });

  it('refuses orbits that do not exist', () => {
    expect(() => circularOrbit(M, 3)).toThrow(/no timelike circular orbit/);
    expect(() => radialEpicyclicFrequency(M, 6)).toThrow(RangeError);
    expect(() => eccentricOrbitConstants(M, 7, 0.6)).toThrow(/separatrix/);
  });
});

describe('circular orbits', () => {
  for (const r of [7, 10, 20, 100]) {
    it(`obeys Kepler's third law exactly in coordinate time at r = ${r}M`, () => {
      const orbit = circularOrbit(M, r);
      const initial = equatorialTimelikeState(model, r, orbit.energy_E, orbit.angular_momentum_Lz, 'turning-point');
      expect(Math.abs(normalizationResidual(model, initial))).toBeLessThan(1e-14);

      const period = (2 * Math.PI) / orbit.omega / orbit.ut; // in proper time
      const result = integrateGeodesic({
        model,
        integrator: tight(),
        initial,
        limits: { initialStep: 1e-2, parameterMax: 5 * period, maxSteps: 1_000_000 },
      });
      expect(Math.abs(result.final.position_x[1] - r) / r).toBeLessThan(1e-10);
      // dphi/dt = sqrt(M / r^3), with no relativistic correction in this time coordinate.
      const omega = result.final.position_x[3] / result.final.position_x[0];
      expect(Math.abs(omega - Math.sqrt(M / r ** 3)) / omega).toBeLessThan(1e-11);
      // And dt/dtau = 1 / sqrt(1 - 3M/r).
      expect(result.final.position_x[0] / result.final.parameter).toBeCloseTo(orbit.ut, 11);
    });
  }
});

describe('the ISCO, measured from both sides', () => {
  it('oscillates radially at omega_r = Omega sqrt(1 - 6M/r) outside it', () => {
    // A nearly circular orbit (e = 1e-4) about r = p. Its radial period in coordinate
    // time, periapsis to periapsis, must match the epicyclic frequency, with corrections
    // of order e^2 ~ 1e-8.
    for (const p of [8, 10, 20]) {
      const e = 1e-4;
      const c = eccentricOrbitConstants(M, p, e);
      const result = integrateGeodesic({
        model,
        integrator: tight(),
        initial: equatorialTimelikeState(model, c.periapsis, c.energy_E, c.angular_momentum_Lz, 'turning-point'),
        limits: { initialStep: 1e-2, parameterMax: 1e7, maxSteps: 1_000_000 },
        events: [{ id: 'periapsis', value: (_x, u) => u[1], direction: 1, terminal: true }],
      });
      expect(result.reason).toBe('event');
      const measured = (2 * Math.PI) / result.final.position_x[0];
      const expected = radialEpicyclicFrequency(M, p);
      expect(Math.abs(measured / expected - 1), `p = ${p}`).toBeLessThan(1e-6);
    }
  });

  it('is unstable inside it, with growth rate sqrt(M (6M - r)) / r^2', () => {
    // omega_r^2 = M (r - 6M) / r^4 turns negative inside the ISCO, so a perturbed circular
    // orbit departs exponentially at |omega_r| per unit coordinate time. Measured by the
    // slope of ln |r - r_c| against t.
    const rc = 5;
    const orbit = circularOrbit(M, rc);
    const f = model.lapseFunction(rc);
    const delta = 1e-9;
    // Energy nudged above the potential maximum by a tiny amount; start moving outward.
    const initial = timelikeState(
      [0, rc, Math.PI / 2, 0],
      [orbit.energy_E / f, delta, 0, orbit.angular_momentum_Lz / (rc * rc)],
    );
    const result = integrateGeodesic({
      model,
      integrator: tight(),
      initial,
      recordPath: true,
      limits: { initialStep: 1e-2, parameterMax: 2000, maxSteps: 1_000_000, maxStep: 0.5 },
      events: [{ id: 'escaped', value: (x) => x[1] - rc * 1.05, direction: 1, terminal: true }],
    });
    const points = result.path!
      .map((s) => [s.position_x[0], Math.log(Math.abs(s.position_x[1] - rc))] as const)
      .filter(([, l]) => l > Math.log(1e-6) && l < Math.log(1e-3));
    expect(points.length).toBeGreaterThan(20);
    const mx = points.reduce((a, [x]) => a + x, 0) / points.length;
    const my = points.reduce((a, [, y]) => a + y, 0) / points.length;
    let sxx = 0;
    let sxy = 0;
    for (const [x, y] of points) {
      sxx += (x - mx) ** 2;
      sxy += (x - mx) * (y - my);
    }
    const slope = sxy / sxx;
    const expected = Math.sqrt(M * (6 * M - rc)) / (rc * rc);
    expect(Math.abs(slope / expected - 1)).toBeLessThan(1e-3);
  });
});

describe('periapsis precession', () => {
  for (const [p, e, expected] of PRECESSION) {
    it(`advances by the exact amount per radial period at p = ${p}M, e = ${e}`, () => {
      const c = eccentricOrbitConstants(M, p, e);
      const result = integrateGeodesic({
        model,
        integrator: tight(),
        initial: equatorialTimelikeState(model, c.periapsis, c.energy_E, c.angular_momentum_Lz, 'turning-point'),
        limits: { initialStep: 1e-3, parameterMax: 1e9, maxSteps: 5_000_000 },
        events: [
          { id: 'apoapsis', value: (_x, u) => u[1], direction: -1, terminal: false },
          { id: 'periapsis', value: (_x, u) => u[1], direction: 1, terminal: true },
        ],
      });
      expect(result.reason).toBe('event');
      const advance = result.final.position_x[3] - 2 * Math.PI;
      // Measured 2e-14 to 4e-13 relative.
      expect(Math.abs(advance - expected) / expected).toBeLessThan(5e-12);
      // The turning points land where the orbit constants put them.
      const apoapsis = result.events.find((h) => h.id === 'apoapsis')!.state.position_x[1];
      expect(Math.abs(apoapsis - c.apoapsis) / c.apoapsis).toBeLessThan(1e-11);
      expect(Math.abs(result.final.position_x[1] - c.periapsis) / c.periapsis).toBeLessThan(1e-11);
    });
  }

  it("converges on Einstein's 6 pi M / p in the weak field", () => {
    const ratios = [20, 100, 1000, 10000].map(
      (p) => periapsisPrecessionExact(M, p, 0.3) / periapsisPrecessionWeakField(M, p),
    );
    for (let i = 1; i < ratios.length; i += 1) expect(ratios[i]).toBeLessThan(ratios[i - 1]);
    expect(ratios[3] - 1).toBeLessThan(1e-3);
    expect(ratios[3]).toBeGreaterThan(1);
  });
});

describe('radial free fall (CLAUDE.md §18: which clock?)', () => {
  it('takes the analytic proper time to fall from 50M to 2.01M', () => {
    const r0 = 50;
    const rEnd = 2.01;
    const result = integrateGeodesic({
      model,
      integrator: tight(),
      // From rest at infinity: E = 1, L = 0, so u^r = -sqrt(2M/r).
      initial: equatorialTimelikeState(model, r0, 1, 0, 'inward'),
      limits: { initialStep: 1e-2, parameterMax: 1e4, maxSteps: 1_000_000 },
      events: [{ id: 'near-horizon', value: (x) => x[1] - rEnd, direction: -1, terminal: true }],
    });
    expect(result.reason).toBe('event');
    const tau = result.final.parameter;
    expect(Math.abs(tau - radialInfallProperTime(M, r0, rEnd)) / tau).toBeLessThan(1e-11);
    const t = result.final.position_x[0];
    const expectedT = radialInfallCoordinateTime(M, rEnd) - radialInfallCoordinateTime(M, r0);
    expect(Math.abs(t - expectedT) / expectedT).toBeLessThan(1e-10);
  });

  it('keeps proper time finite while coordinate time diverges at the horizon', () => {
    // The infalling clock reads almost the same at r = 2.01M and r = 2(1 + 1e-9)M; the
    // Schwarzschild time coordinate keeps growing like -2M ln(r - 2M). Two clocks, two
    // statements, neither of them "time stops".
    const r0 = 50;
    const tauFar = radialInfallProperTime(M, r0, 2.01);
    const tauNear = radialInfallProperTime(M, r0, 2 * (1 + 1e-9));
    const tFar = radialInfallCoordinateTime(M, 2.01) - radialInfallCoordinateTime(M, r0);
    const tNear = radialInfallCoordinateTime(M, 2 * (1 + 1e-9)) - radialInfallCoordinateTime(M, r0);
    expect(tauNear - tauFar).toBeLessThan(0.02);
    expect(tNear - tFar).toBeGreaterThan(30);
  });
});

describe('symplectic integration of long-lived orbits (CLAUDE.md §7.2, §7.3)', () => {
  it('keeps the mass-shell error bounded where explicit Runge-Kutta drifts', () => {
    // 300 radial periods of a p = 20M, e = 0.5 orbit at a fixed step of 2M. Measured:
    // RK4 grows linearly from 2.8e-10 to 2.1e-8; Gauss-Legendre stays at 1.0e-10.
    const c = eccentricOrbitConstants(M, 20, 0.5);
    const initial = equatorialTimelikeState(model, c.periapsis, c.energy_E, c.angular_momentum_Lz, 'turning-point');
    const span = 300 * 2 * Math.PI * Math.pow(20 / 0.75, 1.5);

    const errorProfile = (integrator: RK4Integrator | GaussLegendre4Integrator): number[] => {
      const result = integrateGeodesic({
        model,
        integrator,
        initial,
        recordPath: true,
        limits: { initialStep: 2, maxStep: 2, parameterMax: span, maxSteps: 10_000_000 },
      });
      const path = result.path!;
      return [0.25, 0.5, 1].map((q) =>
        Math.abs(normalizationResidual(model, path[Math.floor(q * (path.length - 1))])),
      );
    };

    const rk4 = errorProfile(new RK4Integrator(STATE_DIM));
    const gl4 = errorProfile(new GaussLegendre4Integrator(STATE_DIM));

    // RK4: secular, roughly linear growth — four times the run, about four times the error.
    expect(rk4[2] / rk4[0]).toBeGreaterThan(3);
    // Gauss-Legendre: bounded — no growth from a quarter of the run to the end.
    expect(gl4[2] / gl4[0]).toBeLessThan(1.5);
    // And by the end it is two orders of magnitude better.
    expect(gl4[2]).toBeLessThan(rk4[2] / 100);
  });

  it('reports itself as symplectic, implicit and fixed-step', () => {
    const gl = new GaussLegendre4Integrator(STATE_DIM);
    expect(gl.symplectic).toBe(true);
    expect(gl.adaptive).toBe(false);
    expect(gl.order).toBe(4);
  });

  it('rejects a step whose implicit stages do not converge, rather than accepting it', () => {
    const gl = new GaussLegendre4Integrator(2, { maxIterations: 2 });
    const y = new Float64Array([1, 0]);
    // Stiff enough that two fixed-point iterations cannot converge at h = 1.
    const result = gl.step((_t, s, d) => {
      d[0] = -50 * s[0];
      d[1] = s[0];
    }, 0, y, 1);
    expect(result.accepted).toBe(false);
    expect(Array.from(y)).toEqual([1, 0]);
  });
});
