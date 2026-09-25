import { CONVENTIONS } from '../conventions.js';
import { ChristoffelSymbols } from '../core/christoffel.js';
import type { Vec4 } from '../core/indices.js';
import { MetricTensor } from '../core/metric-tensor.js';
import { christoffelFromMetricDerivatives } from '../geometry/christoffel-from-derivatives.js';
import {
  IN_DOMAIN,
  type CartesianVec3,
  type ChartGeometry,
  type CoordinateChart,
  type DomainStatus,
  type KillingVector,
  type SpacetimeModel,
  type Symmetries,
} from './spacetime-model.js';

/**
 * Kerr spacetime in Boyer-Lindquist coordinates (ROADMAP.md 4A).
 *
 * Chart: (t, r, theta, phi), covering the stationary exterior r > r_+ only.
 *
 *   ds^2 = -(1 - 2Mr/Sigma) dt^2 - (4 M a r sin^2(theta) / Sigma) dt dphi
 *          + (Sigma/Delta) dr^2 + Sigma dtheta^2
 *          + (A sin^2(theta) / Sigma) dphi^2
 *
 * with
 *
 *   Sigma = r^2 + a^2 cos^2(theta)
 *   Delta = r^2 - 2Mr + a^2
 *   A     = (r^2 + a^2)^2 - a^2 Delta sin^2(theta)
 *
 * in geometric units (G = c = 1) and the (-,+,+,+) signature. The spin parameter a has
 * the dimensions of M; a > 0 spins about +z, so a prograde orbit runs toward +phi.
 *
 * An exact vacuum solution: the Ricci tensor vanishes identically, which was confirmed
 * symbolically from these components before they were used, along with
 * g_{mu nu} g^{nu sigma} = delta^sigma_mu, det g = -Sigma^2 sin^2(theta), and the
 * Kretschmann scalar below. Integrating geodesics here evaluates the consequences of a
 * specified geometry and is not solving the field equations dynamically (CLAUDE.md §1.2).
 *
 * Two structural differences from Schwarzschild drive most of what follows:
 *
 * - The metric has a g_{t phi} cross term, so there is no diagonal static tetrad and the
 *   natural local frame is the zero-angular-momentum observer's (see `observer/zamo.ts`).
 * - Kerr is axisymmetric but *not* spherically symmetric. Its geodesics do not lie in
 *   planes through the centre, so the renderer's orbital-plane reduction does not apply
 *   and the model does not declare spherical symmetry. The reduction refuses a model that
 *   does not, which is why this cannot be reached for by accident.
 */

export const KERR_BOYER_LINDQUIST_CHART: CoordinateChart = Object.freeze<CoordinateChart>({
  id: 'kerr-boyer-lindquist',
  displayName: 'Kerr, Boyer-Lindquist (t, r, theta, phi)',
  kind: 'spherical',
  coordinateNames: ['t', 'r', 'theta', 'phi'],
  horizonPenetrating: false,
  notes:
    'Stationary, axisymmetric exterior chart, valid for r > r_+ = M + sqrt(M^2 - a^2). ' +
    'Delta vanishes at r_+, which is a coordinate singularity and not a curvature one: ' +
    'the Kretschmann scalar is finite there. The chart also degenerates on the polar axis ' +
    'sin(theta) = 0, an artifact of spherical-type coordinates. It does not reduce to a ' +
    'rotating frame at infinity: t and phi are tied to a static observer at infinity. ' +
    'Continuing an integration through r_+ requires horizon-penetrating coordinates such ' +
    'as Kerr-Schild, which this chart is not.',
});

const KERR_SYMMETRIES: Symmetries = Object.freeze({
  stationary: true,
  axisymmetric: true,
  // Not spherically symmetric: a picks out an axis. Geodesics are confined to a cone of
  // constant Carter constant, not to a plane through the centre.
  sphericallySymmetric: false,
});

/**
 * How close to the polar axis the chart is treated as broken down.
 *
 * On the axis g_{phi phi} vanishes, so g^{phi phi} = (Delta - a^2 sin^2 theta) /
 * (Sigma Delta sin^2 theta) diverges. Detected from that component, never from the
 * determinant, which CLAUDE.md §6.1 forbids as a horizon or breakdown test.
 *
 * A geodesic with L_z != 0 never reaches the axis — the Carter constant bounds theta away
 * from 0 and pi — so this is a backstop for nearly axial rays rather than a routine
 * termination.
 */
export const KERR_POLAR_AXIS_SIN_THETA_FLOOR = 1e-7;

/** The outer event horizon r_+ = M + sqrt(M^2 - a^2). */
export function outerHorizonRadius(M: number, a: number): number {
  return M + Math.sqrt(M * M - a * a);
}

/** The inner (Cauchy) horizon r_- = M - sqrt(M^2 - a^2). Outside this chart's domain. */
export function innerHorizonRadius(M: number, a: number): number {
  return M - Math.sqrt(M * M - a * a);
}

/**
 * The static limit (outer ergosurface) r_E(theta) = M + sqrt(M^2 - a^2 cos^2 theta).
 *
 * Inside it g_{tt} > 0: no observer can remain at fixed (r, theta, phi), because that
 * worldline is spacelike. It touches the horizon on the axis and reaches 2M at the
 * equator. Being unable to stand still there is a statement about worldlines, not about
 * the horizon, which lies further in.
 */
export function ergosphereRadius(M: number, a: number, theta: number): number {
  const c = Math.cos(theta);
  return M + Math.sqrt(M * M - a * a * c * c);
}

/**
 * The Kretschmann scalar
 *
 *   K = 48 M^2 (r^2 - a^2 cos^2 theta)
 *       [ (r^2 + a^2 cos^2 theta)^2 - 16 r^2 a^2 cos^2 theta ] / Sigma^6
 *
 * A genuine curvature invariant (CLAUDE.md §5.3). Checked against the full contraction
 * R_{abcd} R^{abcd} computed symbolically from these metric components, agreeing to
 * 1e-15 relative, and reducing to 48 M^2 / r^6 when a = 0.
 *
 * Finite at both horizons and divergent only on the ring r = 0, theta = pi/2: the
 * horizons are features of the chart, the ring is a physical curvature singularity.
 */
export function kretschmann(M: number, a: number, r: number, theta: number): number {
  const c = Math.cos(theta);
  const ac2 = a * a * c * c;
  const r2 = r * r;
  const sigma = r2 + ac2;
  return (48 * M * M * (r2 - ac2) * ((r2 + ac2) * (r2 + ac2) - 16 * r2 * ac2)) / Math.pow(sigma, 6);
}

/**
 * The radius of the equatorial circular photon orbit,
 *
 *   r_ph = 2M { 1 + cos[ (2/3) arccos( -+ a/M ) ] }
 *
 * with the upper sign for a prograde orbit. Bardeen, Press & Teukolsky (1972). It runs
 * from 3M at a = 0 to M (prograde) and 4M (retrograde) at a = M, and each value was
 * checked to satisfy the defining cubic r^3 - 6Mr^2 + 9M^2 r - 4Ma^2 = 0 to 1e-40.
 */
export function equatorialPhotonOrbitRadius(M: number, a: number, sense: OrbitSense): number {
  const sign = sense === 'prograde' ? -1 : 1;
  return 2 * M * (1 + Math.cos((2 / 3) * Math.acos((sign * a) / M)));
}

export type OrbitSense = 'prograde' | 'retrograde';

/**
 * The equatorial ISCO radius (Bardeen, Press & Teukolsky 1972),
 *
 *   Z1 = 1 + (1 - a^2/M^2)^(1/3) [ (1 + a/M)^(1/3) + (1 - a/M)^(1/3) ]
 *   Z2 = sqrt(3 a^2/M^2 + Z1^2)
 *   r_isco = M { 3 + Z2 -+ sqrt[ (3 - Z1)(3 + Z1 + 2 Z2) ] }
 *
 * with the upper sign for a prograde orbit: 6M at a = 0, falling to M for a prograde
 * orbit around an extremal hole and rising to 9M for a retrograde one. Each value was
 * checked against the marginal-stability condition
 * r^2 - 6Mr +- 8a sqrt(Mr) - 3a^2 = 0 to 1e-40.
 */
export function iscoRadius(M: number, a: number, sense: OrbitSense): number {
  const chi = a / M;
  const z1 = 1 + Math.cbrt(1 - chi * chi) * (Math.cbrt(1 + chi) + Math.cbrt(1 - chi));
  const z2 = Math.sqrt(3 * chi * chi + z1 * z1);
  const sign = sense === 'prograde' ? -1 : 1;
  return M * (3 + z2 + sign * Math.sqrt((3 - z1) * (3 + z1 + 2 * z2)));
}

/**
 * The frame-dragging angular velocity omega = -g_{t phi} / g_{phi phi} = 2 M a r / A.
 *
 * The coordinate angular velocity of a zero-angular-momentum observer: the rate at which
 * a locally non-rotating frame is carried around in phi, as seen from infinity. It is a
 * property of the chart's time slicing as much as of the spacetime, and falls off as
 * 2 M a / r^3 far away — the Lense-Thirring form.
 */
export function frameDraggingOmega(M: number, a: number, r: number, theta: number): number {
  const s = Math.sin(theta);
  const r2 = r * r;
  const a2 = a * a;
  const delta = r2 - 2 * M * r + a2;
  const A = (r2 + a2) * (r2 + a2) - a2 * delta * s * s;
  return (2 * M * a * r) / A;
}

/**
 * The Carter constant, in one convention, fixed here and used nowhere else.
 *
 *   Q = p_theta^2 + cos^2(theta) [ a^2 (mu^2 - E^2) + L_z^2 / sin^2(theta) ]
 *
 * with mu^2 = -g_{mu nu} p^mu p^nu, so mu = 1 for a timelike worldline parametrized by
 * proper time and mu = 0 for a null one. In this convention Q vanishes exactly for an
 * equatorial orbit, which makes "is this orbit equatorial?" a question about a number
 * being zero.
 *
 * CLAUDE.md §16 requires exactly one documented convention. The common alternative is
 * the Carter-Walker constant K = Q + (L_z - a E)^2, which does not vanish in the
 * equatorial plane; it is not used anywhere in this project.
 *
 * Unlike E and L_z, Q is not generated by a Killing vector but by a Killing tensor, so
 * the Hamiltonian formulation does not conserve it by construction. That makes it a
 * genuine independent check on a Kerr integration rather than a restatement of the
 * update rule.
 */
export function carterConstant(M: number, a: number, x: Vec4, tangent: Vec4): number {
  const r = x[1];
  const theta = x[2];
  const s = Math.sin(theta);
  const c = Math.cos(theta);
  const sigma = r * r + a * a * c * c;

  const g = new Float64Array(16);
  writeMetric(M, a, x, g);

  const p_t = g[0] * tangent[0] + g[3] * tangent[3];
  const p_phi = g[12] * tangent[0] + g[15] * tangent[3];
  const p_theta = sigma * tangent[2];

  let normSquared = 0;
  for (let mu = 0; mu < 4; mu += 1) {
    for (let nu = 0; nu < 4; nu += 1) {
      const component = g[mu * 4 + nu];
      if (component !== 0) normSquared += component * tangent[mu] * tangent[nu];
    }
  }
  const muSquared = -normSquared;
  const energy_E = -p_t;

  return p_theta * p_theta + c * c * (a * a * (muSquared - energy_E * energy_E) + (p_phi * p_phi) / (s * s));
}

/**
 * The repeated combinations, computed once per evaluation.
 *
 * Sigma, Delta and A appear in every component and in every derivative, so they are
 * built here with their r and theta derivatives rather than being recomputed. Nothing in
 * this struct is a physical quantity on its own; they are the standard abbreviations of
 * the Kerr line element.
 */
interface KerrTerms {
  readonly s: number;
  readonly c: number;
  readonly s2: number;
  readonly sigma: number;
  readonly delta: number;
  readonly A: number;
  readonly dSigma_dr: number;
  readonly dSigma_dtheta: number;
  readonly dDelta_dr: number;
  readonly dA_dr: number;
  readonly dA_dtheta: number;
}

/**
 * Delta = (r - r_+)(r - r_-), which is the same polynomial as r^2 - 2Mr + a^2.
 *
 * The factored form is used because the expanded one loses almost all of its relative
 * accuracy near the horizon: at r = r_+ + 1e-9 the three terms are each of order 1 and
 * cancel to 1e-9, so binary64 leaves only about seven correct digits, and g_rr = Sigma /
 * Delta inherits that. The roots are exact in the same arithmetic (r_+ + r_- = 2M and
 * r_+ r_- = a^2 hold to rounding), so the factored form keeps full relative accuracy
 * right up to the horizon. This is an identity, not an approximation.
 */
function deltaAt(M: number, a: number, r: number): number {
  const root = Math.sqrt(M * M - a * a);
  return (r - (M + root)) * (r - (M - root));
}

function terms(M: number, a: number, r: number, theta: number): KerrTerms {
  const s = Math.sin(theta);
  const c = Math.cos(theta);
  const s2 = s * s;
  const r2 = r * r;
  const a2 = a * a;
  const sigma = r2 + a2 * c * c;
  const delta = deltaAt(M, a, r);
  const A = (r2 + a2) * (r2 + a2) - a2 * delta * s2;
  const dDelta_dr = 2 * r - 2 * M;
  return {
    s,
    c,
    s2,
    sigma,
    delta,
    A,
    dSigma_dr: 2 * r,
    dSigma_dtheta: -2 * a2 * s * c,
    dDelta_dr,
    dA_dr: 4 * r * (r2 + a2) - a2 * s2 * dDelta_dr,
    dA_dtheta: -2 * a2 * delta * s * c,
  };
}

function writeMetric(M: number, a: number, x: Vec4, out: Float64Array): void {
  const r = x[1];
  const { s2, sigma, delta, A } = terms(M, a, r, x[2]);
  out.fill(0);
  out[0] = -(1 - (2 * M * r) / sigma);
  out[3] = out[12] = (-2 * M * a * r * s2) / sigma;
  out[5] = sigma / delta;
  out[10] = sigma;
  out[15] = (A * s2) / sigma;
}

function writeInverseMetric(M: number, a: number, x: Vec4, out: Float64Array): void {
  const r = x[1];
  const { s2, sigma, delta, A } = terms(M, a, r, x[2]);
  const a2 = a * a;
  out.fill(0);
  out[0] = -A / (sigma * delta);
  out[3] = out[12] = (-2 * M * a * r) / (sigma * delta);
  out[5] = delta / sigma;
  out[10] = 1 / sigma;
  out[15] = (delta - a2 * s2) / (sigma * delta * s2);
}

/**
 * d_alpha g_{mu nu}, index alpha * 16 + mu * 4 + nu.
 *
 * Only the r and theta derivatives are non-zero: the metric is stationary and
 * axisymmetric. Each line is the quotient rule applied to one component of the line
 * element above, with the abbreviations from `terms`.
 */
function writeMetricDerivatives(M: number, a: number, x: Vec4, out: Float64Array): void {
  const r = x[1];
  const T = terms(M, a, r, x[2]);
  const { s, c, s2, sigma, delta, A, dSigma_dr, dSigma_dtheta, dDelta_dr, dA_dr, dA_dtheta } = T;
  const sigma2 = sigma * sigma;
  out.fill(0);

  const R = 16; // d_r block
  const TH = 32; // d_theta block

  // g_tt = -1 + 2Mr/Sigma
  out[R + 0] = (2 * M * (sigma - r * dSigma_dr)) / sigma2;
  out[TH + 0] = (-2 * M * r * dSigma_dtheta) / sigma2;

  // g_t_phi = -2 M a r sin^2(theta) / Sigma
  const dg_tphi_dr = (-2 * M * a * s2 * (sigma - r * dSigma_dr)) / sigma2;
  const dg_tphi_dtheta =
    (-2 * M * a * r * (2 * s * c * sigma - s2 * dSigma_dtheta)) / sigma2;
  out[R + 3] = out[R + 12] = dg_tphi_dr;
  out[TH + 3] = out[TH + 12] = dg_tphi_dtheta;

  // g_rr = Sigma / Delta
  out[R + 5] = (dSigma_dr * delta - sigma * dDelta_dr) / (delta * delta);
  out[TH + 5] = dSigma_dtheta / delta;

  // g_theta_theta = Sigma
  out[R + 10] = dSigma_dr;
  out[TH + 10] = dSigma_dtheta;

  // g_phi_phi = A sin^2(theta) / Sigma
  out[R + 15] = (s2 * (dA_dr * sigma - A * dSigma_dr)) / sigma2;
  out[TH + 15] =
    ((dA_dtheta * s2 + 2 * A * s * c) * sigma - A * s2 * dSigma_dtheta) / sigma2;
}

/**
 * d_alpha g^{mu nu}, index alpha * 16 + mu * 4 + nu.
 *
 * The only geometric input the Hamiltonian equations need besides g^{mu nu} itself.
 * Written out by the quotient rule from the closed-form inverse rather than by inverting
 * and differencing, and checked against both central differences of `inverseMetricInto`
 * and an independent SymPy evaluation.
 */
function writeInverseMetricDerivatives(M: number, a: number, x: Vec4, out: Float64Array): void {
  const r = x[1];
  const T = terms(M, a, r, x[2]);
  const { s, c, s2, sigma, delta, A, dSigma_dr, dSigma_dtheta, dDelta_dr, dA_dr, dA_dtheta } = T;
  const a2 = a * a;
  const sigmaDelta = sigma * delta;
  const sigmaDelta2 = sigmaDelta * sigmaDelta;
  const sigma2 = sigma * sigma;
  out.fill(0);

  const R = 16;
  const TH = 32;

  // g^tt = -A / (Sigma Delta)
  out[R + 0] = -(dA_dr * sigmaDelta - A * (dSigma_dr * delta + sigma * dDelta_dr)) / sigmaDelta2;
  out[TH + 0] = -(dA_dtheta * sigmaDelta - A * dSigma_dtheta * delta) / sigmaDelta2;

  // g^t_phi = -2 M a r / (Sigma Delta)
  const dgi_tphi_dr =
    (-2 * M * a * (sigmaDelta - r * (dSigma_dr * delta + sigma * dDelta_dr))) / sigmaDelta2;
  const dgi_tphi_dtheta = (2 * M * a * r * dSigma_dtheta * delta) / sigmaDelta2;
  out[R + 3] = out[R + 12] = dgi_tphi_dr;
  out[TH + 3] = out[TH + 12] = dgi_tphi_dtheta;

  // g^rr = Delta / Sigma
  out[R + 5] = (dDelta_dr * sigma - delta * dSigma_dr) / sigma2;
  out[TH + 5] = (-delta * dSigma_dtheta) / sigma2;

  // g^theta_theta = 1 / Sigma
  out[R + 10] = -dSigma_dr / sigma2;
  out[TH + 10] = -dSigma_dtheta / sigma2;

  // g^phi_phi = (Delta - a^2 sin^2 theta) / (Sigma Delta sin^2 theta)
  const numerator = delta - a2 * s2;
  const denominator = sigmaDelta * s2;
  const dDen_dr = (dSigma_dr * delta + sigma * dDelta_dr) * s2;
  const dDen_dtheta = dSigma_dtheta * delta * s2 + sigmaDelta * 2 * s * c;
  out[R + 15] = (dDelta_dr * denominator - numerator * dDen_dr) / (denominator * denominator);
  out[TH + 15] =
    (-2 * a2 * s * c * denominator - numerator * dDen_dtheta) / (denominator * denominator);
}

function killingVectors(): readonly KillingVector[] {
  return Object.freeze([
    {
      id: 'd_dt',
      displayName: 'Time translation d/dt (stationarity)',
      conservedQuantityName: 'energy_E',
      sign: -1 as const,
      at: (): Vec4 => [1, 0, 0, 0],
    },
    {
      id: 'd_dphi',
      displayName: 'Azimuthal rotation d/dphi (axisymmetry)',
      conservedQuantityName: 'angular_momentum_Lz',
      sign: 1 as const,
      at: (): Vec4 => [0, 0, 0, 1],
    },
  ]);
}

/**
 * The chart's relation to the auxiliary Cartesian visualization axes.
 *
 * Boyer-Lindquist r is not a spherical radius: surfaces of constant r are confocal
 * oblate spheroids, and the standard relation to Kerr-Schild Cartesian coordinates is
 *
 *   x = sqrt(r^2 + a^2) sin(theta) cos(phi),  y = sqrt(r^2 + a^2) sin(theta) sin(phi),
 *   z = r cos(theta)
 *
 * which is used here so that r = 0 maps to a ring of radius a in the equatorial plane,
 * as it should. The difference from a naive spherical mapping is O(a^2 / r^2) and
 * vanishes at the radii where the background is sampled.
 *
 * Directions are reported in the zero-angular-momentum observer's orthonormal spatial
 * triad. A static frame would do for Schwarzschild but does not exist inside the
 * ergosphere, and the ZAMO frame is asymptotically inertial, so far from the hole this
 * is the direction the ray is really travelling. Close in it is that observer's view and
 * nothing more (CLAUDE.md §1.3).
 */
function geometryFor(M: number, a: number): ChartGeometry {
  return {
    spatialRadius: (x: Vec4): number => x[1],

    toCartesianPosition: (x: Vec4): CartesianVec3 => {
      const [, r, theta, phi] = x;
      const rho = Math.sqrt(r * r + a * a) * Math.sin(theta);
      return [rho * Math.cos(phi), rho * Math.sin(phi), r * Math.cos(theta)];
    },

    toCartesianDirection: (x: Vec4, tangent: Vec4): CartesianVec3 => {
      const [, r, theta, phi] = x;
      const { s, c, sigma, delta, A } = terms(M, a, r, theta);
      const omega = (2 * M * a * r) / A;

      const vRadial = delta > 0 ? tangent[1] * Math.sqrt(sigma / delta) : tangent[1];
      const vPolar = Math.sqrt(sigma) * tangent[2];
      const vAzimuthal = Math.sqrt(A / sigma) * s * (tangent[3] - omega * tangent[0]);

      const cosPhi = Math.cos(phi);
      const sinPhi = Math.sin(phi);
      const rHat: CartesianVec3 = [s * cosPhi, s * sinPhi, c];
      const thetaHat: CartesianVec3 = [c * cosPhi, c * sinPhi, -s];
      const phiHat: CartesianVec3 = [-sinPhi, cosPhi, 0];

      return [
        vRadial * rHat[0] + vPolar * thetaHat[0] + vAzimuthal * phiHat[0],
        vRadial * rHat[1] + vPolar * thetaHat[1] + vAzimuthal * phiHat[1],
        vRadial * rHat[2] + vPolar * thetaHat[2] + vAzimuthal * phiHat[2],
      ];
    },
  };
}

class KerrSpacetime implements SpacetimeModel {
  readonly id: string;
  readonly displayName: string;
  readonly classification = 'exact-analytical' as const;
  readonly chart = KERR_BOYER_LINDQUIST_CHART;
  readonly conventions = CONVENTIONS;
  readonly parameters: Readonly<Record<string, number>>;
  readonly killingVectors = killingVectors();
  readonly symmetries = KERR_SYMMETRIES;
  readonly geometry: ChartGeometry;
  readonly description =
    'Exact vacuum solution of the Einstein field equations: stationary, axisymmetric and ' +
    'rotating. Ricci-flat everywhere it is defined. Axisymmetric but not spherically ' +
    'symmetric, so geodesics do not lie in planes through the centre. Integrating ' +
    'geodesics here evaluates the consequences of a specified geometry, not a dynamical ' +
    'solution of the field equations (CLAUDE.md §1.2).';

  readonly M: number;
  readonly a: number;

  private readonly scratchInverse = new Float64Array(16);
  private readonly scratchDerivatives = new Float64Array(64);

  constructor(M: number, a: number) {
    if (!(M > 0) || !Number.isFinite(M)) {
      throw new RangeError(`Kerr: the mass parameter must be finite and positive, received ${M}.`);
    }
    if (!Number.isFinite(a) || Math.abs(a) > M) {
      throw new RangeError(
        `Kerr: the spin parameter must satisfy |a| <= M, received a = ${a} with M = ${M}. ` +
          'Beyond that there is no horizon and the ring singularity is naked, which is a ' +
          'different spacetime and outside this project’s scope.',
      );
    }
    this.M = M;
    this.a = a;
    this.id = 'kerr';
    this.displayName = `Kerr spacetime (M = ${M}, a = ${a})`;
    this.parameters = Object.freeze({ M, a });
    this.geometry = geometryFor(M, a);
  }

  /** The outer horizon radius for this model. */
  get horizonRadius(): number {
    return outerHorizonRadius(this.M, this.a);
  }

  /** Delta = (r - r_+)(r - r_-) = r^2 - 2Mr + a^2, which vanishes at both horizons. */
  deltaFunction(r: number): number {
    return deltaAt(this.M, this.a, r);
  }

  metricAt(x: Vec4): MetricTensor {
    this.assertInDomain(x, 'metricAt');
    const g = new Float64Array(16);
    const gInv = new Float64Array(16);
    writeMetric(this.M, this.a, x, g);
    writeInverseMetric(this.M, this.a, x, gInv);
    return new MetricTensor(g, gInv);
  }

  christoffelAt(x: Vec4): ChristoffelSymbols {
    const out = new Float64Array(64);
    this.christoffelInto(x, out);
    return new ChristoffelSymbols(out);
  }

  christoffelInto(x: Vec4, out: Float64Array): void {
    writeInverseMetric(this.M, this.a, x, this.scratchInverse);
    writeMetricDerivatives(this.M, this.a, x, this.scratchDerivatives);
    christoffelFromMetricDerivatives(this.scratchInverse, this.scratchDerivatives, out);
  }

  metricInto(x: Vec4, out: Float64Array): void {
    writeMetric(this.M, this.a, x, out);
  }

  inverseMetricInto(x: Vec4, out: Float64Array): void {
    writeInverseMetric(this.M, this.a, x, out);
  }

  inverseMetricDerivativesInto(x: Vec4, out: Float64Array): void {
    writeInverseMetricDerivatives(this.M, this.a, x, out);
  }

  /** Analytic d_alpha g_{mu nu}, exposed for cross-validation. */
  metricDerivativesInto(x: Vec4, out: Float64Array): void {
    writeMetricDerivatives(this.M, this.a, x, out);
  }

  domainCheck(x: Vec4): DomainStatus {
    const r = x[1];
    const horizon = this.horizonRadius;
    if (!(r > horizon)) {
      return {
        inDomain: false,
        code: 'coordinate-breakdown',
        reason:
          `Delta = r^2 - 2Mr + a^2 = ${this.deltaFunction(r)} at r = ${r}, at or inside the ` +
          `outer horizon r_+ = ${horizon}. Boyer-Lindquist coordinates cover the exterior ` +
          'only; the breakdown is in the chart, not in the geometry, whose curvature ' +
          'invariants are finite there.',
      };
    }
    const sinTheta = Math.abs(Math.sin(x[2]));
    if (!(sinTheta > KERR_POLAR_AXIS_SIN_THETA_FLOOR)) {
      return {
        inDomain: false,
        code: 'coordinate-breakdown',
        reason:
          `sin(theta) = ${sinTheta} is at the polar axis, where g^{phi phi} diverges. This is ` +
          'an artifact of spherical-type coordinates rather than anything physical.',
      };
    }
    return IN_DOMAIN;
  }

  private assertInDomain(x: Vec4, caller: string): void {
    const status = this.domainCheck(x);
    if (!status.inDomain) throw new RangeError(`Kerr.${caller}: ${status.reason}`);
  }
}

/** Kerr spacetime with mass M and spin parameter a, in Boyer-Lindquist coordinates. */
export function kerr(M: number, a: number): SpacetimeModel & {
  readonly M: number;
  readonly a: number;
  readonly horizonRadius: number;
  deltaFunction(r: number): number;
  metricDerivativesInto(x: Vec4, out: Float64Array): void;
} {
  return new KerrSpacetime(M, a);
}
