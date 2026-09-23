import { describe, expect, it } from 'vitest';
import type { Vec4 } from '../../src/physics/core/indices.js';
import { HAMILTONIAN } from '../../src/physics/geodesic/formulation.js';
import { integrateGeodesic } from '../../src/physics/geodesic/integrate.js';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import {
  frequencyRatio,
  keplerianEmitter,
  STATIC_EMITTER,
} from '../../src/physics/observer/frequency-shift.js';
import {
  freeFallingObserver,
  generateNullRay,
  staticObserver,
} from '../../src/physics/observer/observer.js';
import { boostTetrad, orthonormalityResidual, Tetrad } from '../../src/physics/observer/tetrad.js';
import { gaussLegendreUnitInterval } from '../../src/physics/spacetimes/schwarzschild-analytic.js';
import { criticalImpactParameter, schwarzschild } from '../../src/physics/spacetimes/schwarzschild.js';
import { circularOrbit, equatorialTimelikeState } from '../../src/physics/spacetimes/schwarzschild-orbits.js';
import { isCaptured } from '../../src/physics/spacetimes/schwarzschild-rays.js';
import {
  diskOrbitEnergy,
  NOVIKOV_THORNE_ASYMPTOTIC_C,
  NOVIKOV_THORNE_EFFICIENCY,
  NOVIKOV_THORNE_PEAK_RADIUS,
  novikovThorneFluxPerAccretionRate,
  thinDiskScale,
} from '../../src/physics/spacetimes/novikov-thorne.js';

/**
 * Milestone 3 physics: the Novikov-Thorne disk, the frequency shift, and observers in
 * motion. Reference values from mpmath (disk flux) and an independent Python calculation
 * with CODATA 2018 constants (physical temperatures).
 */

const M = 1;
const model = schwarzschild(M);
const tight = (): Dopri5Integrator =>
  new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-12, relative: 1e-12 } });

/** F M^2 / Mdot at r / M, mpmath, 18 digits. */
const FLUX: readonly (readonly [number, number])[] = [
  [6.5, 2.36305761017507855e-6],
  [8, 1.1447491074846865e-5],
  [10, 1.3568754638046947e-5],
  [20, 4.72514406575524364e-6],
  [50, 5.06867951247639854e-7],
  [100, 7.78585359203804547e-8],
];

describe('Novikov-Thorne flux (ROADMAP.md 3.4)', () => {
  it('matches the mpmath closed form', () => {
    for (const [r, expected] of FLUX) {
      expect(Math.abs(novikovThorneFluxPerAccretionRate(M, r) / expected - 1), `r = ${r}M`).toBeLessThan(1e-12);
    }
  });

  it('vanishes at the ISCO and is never negative', () => {
    expect(novikovThorneFluxPerAccretionRate(M, 6)).toBe(0);
    expect(novikovThorneFluxPerAccretionRate(M, 5)).toBe(0);
    for (let r = 6.0001; r < 200; r *= 1.01) {
      expect(novikovThorneFluxPerAccretionRate(M, r)).toBeGreaterThanOrEqual(0);
    }
  });

  it('peaks at r = 9.550928 M', () => {
    let best = 0;
    let bestR = 0;
    for (let r = 8; r < 12; r += 1e-6) {
      const F = novikovThorneFluxPerAccretionRate(M, r);
      if (F > best) {
        best = F;
        bestR = r;
      }
    }
    expect(Math.abs(bestR - NOVIKOV_THORNE_PEAK_RADIUS)).toBeLessThan(2e-6);
  });

  it('radiates exactly the ISCO binding energy: L_inf = Mdot (1 - sqrt(8/9))', () => {
    // The model's own energy balance, and the strongest check on the flux formula: the
    // luminosity reaching infinity from both faces, Integral 4 pi r E(r) F(r) dr, must equal
    // the rest-mass energy the gas loses spiralling in from infinity to the ISCO.
    const { nodes, weights } = gaussLegendreUnitInterval(64);
    const integrand = (r: number): number =>
      4 * Math.PI * r * diskOrbitEnergy(M, r) * novikovThorneFluxPerAccretionRate(M, r);
    let luminosity = 0;
    const segments: readonly (readonly [number, number])[] = [
      [6, 7], [7, 10], [10, 30], [30, 100], [100, 1000],
    ];
    for (const [a, b] of segments) {
      for (let i = 0; i < nodes.length; i += 1) {
        luminosity += weights[i] * (b - a) * integrand(a + (b - a) * nodes[i]);
      }
    }
    // Tail to infinity with r = 1000 / w^2, dr = -2000 dw / w^3. Squaring matters: far out
    // the flux carries a sqrt(M/r) term, so in u = 1000/r the integrand behaves like sqrt(u)
    // at u = 0 and Gauss-Legendre converges only algebraically. In w = sqrt(u) it is smooth.
    for (let i = 0; i < nodes.length; i += 1) {
      const w = nodes[i];
      luminosity += weights[i] * integrand(1000 / (w * w)) * (2000 / (w * w * w));
    }
    expect(NOVIKOV_THORNE_EFFICIENCY).toBeCloseTo(0.057190958417936634, 15);
    expect(Math.abs(luminosity / NOVIKOV_THORNE_EFFICIENCY - 1)).toBeLessThan(1e-9);
  });

  it('approaches (3 M Mdot / 8 pi r^3)(1 - C sqrt(M/r)) far out, with C = sqrt 6 + sqrt 3 ln(1 + sqrt 2)', () => {
    // Not the naive Newtonian 1 - sqrt(6M/r): the logarithms in the relativistic flux change
    // the inner-boundary term, giving C = 3.976 rather than sqrt 6 = 2.449.
    expect(NOVIKOV_THORNE_ASYMPTOTIC_C).toBeCloseTo(3.9761, 4);
    let previous = Infinity;
    for (const r of [1e4, 1e6, 1e8]) {
      const scaled = (novikovThorneFluxPerAccretionRate(M, r) * 8 * Math.PI * r ** 3) / (3 * M);
      const correction = (1 - scaled) / (NOVIKOV_THORNE_ASYMPTOTIC_C * Math.sqrt(M / r));
      const miss = Math.abs(correction - 1);
      expect(miss, `r = ${r}`).toBeLessThan(previous);
      previous = miss;
    }
    expect(previous).toBeLessThan(1e-3);
  });
});

describe('disk temperature in physical units', () => {
  it('reproduces an independent CODATA 2018 calculation', () => {
    // From separate Python, same constants: T_max for (M, f_Edd).
    const cases: readonly (readonly [number, number, number])[] = [
      [1e9, 0.1, 39488.99552707843],
      [1e8, 0.3, 92417.96482446574],
      [4e10, 0.01, 8830.007836173249],
    ];
    for (const [massSolar, eddingtonFraction, expected] of cases) {
      const T = thinDiskScale({ massSolar, eddingtonFraction }).peakTemperatureK;
      expect(Math.abs(T / expected - 1), `M = ${massSolar}`).toBeLessThan(1e-9);
    }
  });

  it('scales as (f_Edd / M)^(1/4)', () => {
    // sigma T^4 ~ Mdot / r_g^2 ~ f_Edd M / M^2.
    const base = thinDiskScale({ massSolar: 1e9, eddingtonFraction: 0.1 }).peakTemperatureK;
    const heavier = thinDiskScale({ massSolar: 16e9, eddingtonFraction: 0.1 }).peakTemperatureK;
    const faster = thinDiskScale({ massSolar: 1e9, eddingtonFraction: 1.6 }).peakTemperatureK;
    expect(heavier / base).toBeCloseTo(0.5, 12);
    expect(faster / base).toBeCloseTo(2, 12);
  });

  it('is hot for any real black hole: a physical thin disk is not orange', () => {
    // A stellar-mass hole's disk peaks in X-rays; even the most massive known holes at low
    // accretion rates stay near 10^4 K. The colour on screen is whatever this gives.
    expect(thinDiskScale({ massSolar: 10, eddingtonFraction: 0.1 }).peakTemperatureK).toBeGreaterThan(1e6);
    expect(thinDiskScale({ massSolar: 4e10, eddingtonFraction: 0.01 }).peakTemperatureK).toBeGreaterThan(8000);
  });
});

describe('frequency shift (ROADMAP.md 3.3)', () => {
  it('gives the gravitational redshift sqrt(f_emit / f_obs) between static observers', () => {
    // A static emitter at r = 8M seen by a static observer at r = 20M, along a radial ray.
    const observer = staticObserver(model, [0, 20, Math.PI / 2, 0]);
    const result = integrateGeodesic({
      model,
      integrator: tight(),
      initial: generateNullRay(observer, [-1, 0, 0]),
      limits: { initialStep: 1e-3, parameterMax: 100, maxSteps: 100_000 },
      events: [{ id: 'emitter', value: (x) => x[1] - 8, direction: -1, terminal: true }],
    });
    const f8 = model.lapseFunction(8);
    const photon = { energy: Math.sqrt(f8) * result.final.tangent[0], momentum: [0, 0, 0] as const };
    const g = frequencyRatio(photon, STATIC_EMITTER);
    expect(Math.abs(g - Math.sqrt(f8 / model.lapseFunction(20)))).toBeLessThan(1e-12);
    expect(g).toBeLessThan(1); // climbing out: redshift
  });

  it('agrees with the covariant formula for rays hitting a Keplerian disk', () => {
    // Two independent routes to g. Covariant: k . u_emit = u^t (p_t + Omega p_phi), with p_t
    // and p_phi the conserved covariant momenta of the traced ray. Static-frame: boost the
    // emitter's velocity against the photon's static-frame energy and direction — what the
    // renderer does. Integrated here in the full world chart, with no orbital-plane
    // reduction, so nothing is shared with the renderer's path.
    const camera: Vec4 = [0, 30, (60 * Math.PI) / 180, 0];
    const observer = staticObserver(model, camera);
    const session = HAMILTONIAN.bind(model);
    let checked = 0;
    for (const [a, b] of [[0.1, 0.15], [-0.12, 0.2], [0.05, -0.18], [0.2, 0.05], [-0.2, -0.1]] as const) {
      const n = [-1, a, b];
      const norm = Math.hypot(...n);
      const result = integrateGeodesic({
        model,
        integrator: tight(),
        initial: generateNullRay(observer, [n[0] / norm, n[1] / norm, n[2] / norm]),
        session,
        limits: { initialStep: 1e-3, parameterMax: 500, maxSteps: 200_000 },
        events: [
          {
            id: 'disk',
            value: (x) => Math.cos(x[2]),
            direction: 0,
            terminal: (s) => s.position_x[1] >= 6 && s.position_x[1] <= 20,
          },
        ],
        terminator: (x, k) => isCaptured(model, x, k),
      });
      if (result.reason !== 'event') continue;
      checked += 1;
      const x = result.final.position_x;
      const r = x[1];
      const pt = result.packed[4];
      const pphi = result.packed[7];
      const orbit = circularOrbit(M, r);
      const gCovariant = 1 / (orbit.ut * (pt + orbit.omega * pphi));

      const f = model.lapseFunction(r);
      const photon = {
        energy: Math.sqrt(f) * result.final.tangent[0],
        momentum: model.geometry.toCartesianDirection(x, result.final.tangent),
      };
      const gStatic = frequencyRatio(photon, keplerianEmitter(M, model.geometry.toCartesianPosition(x)));
      expect(Math.abs(gStatic / gCovariant - 1), `ray (${a}, ${b})`).toBeLessThan(1e-12);
    }
    expect(checked).toBeGreaterThanOrEqual(3);
  });

  it('gives g = sqrt(1 - 3M/r) / sqrt(f_obs) for rays with no angular momentum about the disk axis', () => {
    // A ray confined to a meridional plane has p_phi = 0, so only the gravitational and
    // transverse-Doppler parts of the shift survive: g = 1 / (u^t sqrt(f_obs)).
    const camera: Vec4 = [0, 30, (40 * Math.PI) / 180, 0];
    const observer = staticObserver(model, camera);
    // Tilts toward the equator chosen to meet the disk plane between about 6M and 25M. A
    // ray aimed almost straight at the hole would be captured before reaching the plane.
    for (const tilt of [0.3, 0.4, 0.5]) {
      // Legs (1) radial and (2) polar only: no azimuthal component, so p_phi = 0 exactly.
      const n = [-Math.cos(tilt), Math.sin(tilt), 0] as const;
      const result = integrateGeodesic({
        model,
        integrator: tight(),
        initial: generateNullRay(observer, n),
        limits: { initialStep: 1e-3, parameterMax: 500, maxSteps: 200_000 },
        events: [{ id: 'disk', value: (x) => Math.cos(x[2]), direction: 0, terminal: true }],
        terminator: (x, k) => isCaptured(model, x, k),
      });
      expect(result.reason, `tilt ${tilt}`).toBe('event');
      expect(result.packed[7]).toBe(0);
      const r = result.final.position_x[1];
      const f = model.lapseFunction(r);
      const photon = {
        energy: Math.sqrt(f) * result.final.tangent[0],
        momentum: model.geometry.toCartesianDirection(result.final.position_x, result.final.tangent),
      };
      const g = frequencyRatio(photon, keplerianEmitter(M, model.geometry.toCartesianPosition(result.final.position_x)));
      const expected = Math.sqrt(1 - (3 * M) / r) / Math.sqrt(model.lapseFunction(30));
      expect(Math.abs(g / expected - 1), `tilt ${tilt}, r = ${r.toFixed(3)}`).toBeLessThan(1e-10);
    }
  });
});

describe('observers in motion (ROADMAP.md 3.1)', () => {
  it('keeps a boosted frame orthonormal', () => {
    const metric = model.metricAt([0, 12, 1.1, 0.4]);
    const boosted = boostTetrad(Tetrad.diagonalStatic(metric), [0.3, -0.5, 0.2]);
    expect(orthonormalityResidual(metric, boosted)).toBeLessThan(1e-14);
  });

  it('gives the free-faller the radial-infall four-velocity, which is a geodesic', () => {
    for (const r of [4, 10, 20, 100]) {
      const faller = freeFallingObserver(model, [0, r, Math.PI / 2, 0], M);
      const u = faller.tetrad.four_velocity_u();
      const f = model.lapseFunction(r);
      expect(u[0]).toBeCloseTo(1 / f, 13);
      expect(u[1]).toBeCloseTo(-Math.sqrt((2 * M) / r), 14);
      // The same four-velocity as a timelike geodesic falling from rest at infinity.
      const geodesic = equatorialTimelikeState(model, r, 1, 0, 'inward');
      for (let mu = 0; mu < 4; mu += 1) expect(u[mu]).toBeCloseTo(geodesic.tangent[mu], 13);
      expect(orthonormalityResidual(model.metricAt([0, r, Math.PI / 2, 0]), faller.tetrad)).toBeLessThan(1e-14);
    }
  });

  it('sees the sky behind it redshifted by exactly 1 / (1 + sqrt(2M/r))', () => {
    // Looking straight up, away from the hole, at light from a static source at infinity.
    // Gravitational blueshift 1/sqrt(f) times the Doppler redshift of receding from the
    // source, sqrt((1 - v)/(1 + v)), with f = 1 - v^2, collapses to 1 / (1 + v).
    for (const r of [5, 20, 100]) {
      const faller = freeFallingObserver(model, [0, r, Math.PI / 2, 0], M);
      const ray = generateNullRay(faller, [1, 0, 0]);
      const packed = HAMILTONIAN.bind(model).pack(ray, new Float64Array(STATE_DIM));
      // At infinity a static source has k . u = p_t, conserved exactly along the ray.
      const g = 1 / packed[4];
      expect(Math.abs(g - 1 / (1 + Math.sqrt((2 * M) / r))), `r = ${r}`).toBeLessThan(1e-14);
    }
  });

  it('sees the shadow shrunk by relativistic aberration', () => {
    // A static observer at r sees the shadow edge at sin(psi_s) = b_c sqrt(f) / r. The
    // infaller, moving inward at v = sqrt(2M/r), sees it at cos(psi') = (cos psi_s + v) /
    // (1 + v cos psi_s): smaller, because aberration crowds what lies ahead toward the
    // direction of motion. Measured by bisecting on capture in the infaller's own frame.
    const r = 20;
    const faller = freeFallingObserver(model, [0, r, Math.PI / 2, 0], M);
    const v = Math.sqrt((2 * M) / r);
    const psiStatic = Math.asin((criticalImpactParameter(M) * Math.sqrt(model.lapseFunction(r))) / r);
    const predicted = Math.acos((Math.cos(psiStatic) + v) / (1 + v * Math.cos(psiStatic)));

    const captured = (psi: number): boolean => {
      let hit = false;
      integrateGeodesic({
        model,
        integrator: tight(),
        initial: generateNullRay(faller, [-Math.cos(psi), 0, -Math.sin(psi)]),
        limits: { initialStep: 1e-3, parameterMax: 5000, maxSteps: 400_000 },
        events: [{ id: 'out', value: (x) => x[1] - 500, direction: 1, terminal: true }],
        terminator: (x, k) => (hit = isCaptured(model, x, k)),
      });
      return hit;
    };
    let lo = predicted * 0.5;
    let hi = predicted * 1.5;
    expect(captured(lo)).toBe(true);
    expect(captured(hi)).toBe(false);
    for (let i = 0; i < 45; i += 1) {
      const mid = 0.5 * (lo + hi);
      if (captured(mid)) lo = mid;
      else hi = mid;
    }
    const measured = 0.5 * (lo + hi);
    expect(Math.abs(measured / predicted - 1)).toBeLessThan(1e-8);
    expect(measured).toBeLessThan(psiStatic);
  });
});
