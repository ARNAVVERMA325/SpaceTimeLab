import { CONVENTIONS } from '../conventions.js';
import { ChristoffelSymbols } from '../core/christoffel.js';
import { MetricTensor } from '../core/metric-tensor.js';
import type { Vec4 } from '../core/indices.js';
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
 * Schwarzschild spacetime in standard Schwarzschild coordinates (ROADMAP.md 2A.1).
 *
 * Chart: (t, r, theta, phi), covering the static exterior region r > 2M only.
 *
 *   ds^2 = -f dt^2 + f^-1 dr^2 + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2,   f = 1 - 2M/r
 *
 * in geometric units (G = c = 1) and the (-,+,+,+) signature.
 *
 * This is an exact vacuum solution: the Ricci tensor vanishes identically. Integrating
 * geodesics through it evaluates the consequences of a specified geometry and is not
 * solving the Einstein field equations dynamically (CLAUDE.md §1.2).
 *
 * The Christoffel symbols below were derived symbolically from the definition in
 * CLAUDE.md §2 and cross-checked by confirming that the resulting Riemann tensor gives
 * R_mu_nu = 0 and the Kretschmann scalar K = 48 M^2 / r^6. The non-zero symbols are:
 *
 *   Gamma^t_{t r}         =  M / (r^2 f)
 *   Gamma^r_{t t}         =  M f / r^2
 *   Gamma^r_{r r}         = -M / (r^2 f)
 *   Gamma^r_{theta theta} = -r f
 *   Gamma^r_{phi phi}     = -r f sin^2(theta)
 *   Gamma^theta_{r theta} =  1 / r
 *   Gamma^theta_{phi phi} = -sin(theta) cos(theta)
 *   Gamma^phi_{r phi}     =  1 / r
 *   Gamma^phi_{theta phi} =  cot(theta)
 *
 * all others zero, with the lower index pair symmetric.
 */

export const SCHWARZSCHILD_CHART: CoordinateChart = Object.freeze<CoordinateChart>({
  id: 'schwarzschild-spherical',
  displayName: 'Schwarzschild (t, r, theta, phi)',
  kind: 'spherical',
  coordinateNames: ['t', 'r', 'theta', 'phi'],
  horizonPenetrating: false,
  notes:
    'Static exterior chart, valid for r > 2M. The chart breaks down at r = 2M, which is ' +
    'a coordinate singularity and not a curvature singularity: the Kretschmann scalar ' +
    'K = 48 M^2 / r^6 is finite there. It also degenerates on the polar axis ' +
    'sin(theta) = 0, which is an artifact of spherical coordinates rather than anything ' +
    'physical. Continuing an integration through r = 2M requires horizon-penetrating ' +
    'coordinates, which this chart is not.',
});

const SCHWARZSCHILD_SYMMETRIES: Symmetries = Object.freeze({
  stationary: true,
  axisymmetric: true,
  sphericallySymmetric: true,
});

/**
 * How close to the polar axis the chart is treated as broken down.
 *
 * On the axis g_phi_phi = r^2 sin^2(theta) vanishes, so g^phi_phi diverges and
 * Gamma^phi_{theta phi} = cot(theta) is unbounded. Per CLAUDE.md §6.1 the breakdown is
 * detected from those components rather than from the metric determinant, which is
 * -r^4 sin^2(theta) and vanishes on the axis for the same reason but also happens to
 * vanish there for a chart that is otherwise perfectly well behaved.
 *
 * At sin(theta) = 1e-6 the connection is already of order 1e6 and an error-controlled
 * integrator has effectively stalled. Rays are kept off the axis entirely by the
 * orbital-plane reduction the renderer uses, so this bound is a backstop rather than a
 * routine termination.
 */
export const POLAR_AXIS_SIN_THETA_FLOOR = 1e-6;

/**
 * The Kretschmann scalar for Schwarzschild, K = R_abcd R^abcd = 48 M^2 / r^6.
 *
 * A genuine curvature invariant (CLAUDE.md §5.3), not a coordinate-dependent quantity.
 * It is finite at r = 2M and diverges only as r -> 0, which is the statement that the
 * horizon is a feature of the chart while r = 0 is a physical curvature singularity.
 */
export function kretschmann(M: number, r: number): number {
  return (48 * M * M) / Math.pow(r, 6);
}

/** The event-horizon radius r = 2M in geometric units. */
export function horizonRadius(M: number): number {
  return 2 * M;
}

/**
 * The photon-sphere radius r = 3M.
 *
 * The circular null orbit is where the null effective potential V(r) = f / r^2 is
 * stationary: V'(r) = (6M - 2r) / r^4 vanishes at r = 3M, and V''(3M) < 0, so the orbit
 * is a maximum of the potential and therefore unstable.
 */
export function photonSphereRadius(M: number): number {
  return 3 * M;
}

/**
 * The critical impact parameter b_c = 3 sqrt(3) M.
 *
 * At r = 3M a circular null orbit has (E/L)^2 = f / r^2 = 1 / (27 M^2), so
 * b = L / E = 3 sqrt(3) M ~ 5.196 M. Rays with b < b_c are captured; rays with
 * b > b_c escape after turning at some r_0 > 3M.
 */
export function criticalImpactParameter(M: number): number {
  return 3 * Math.sqrt(3) * M;
}

function killingVectors(): readonly KillingVector[] {
  return Object.freeze([
    {
      id: 'd_dt',
      displayName: 'Time translation d/dt (stationarity)',
      conservedQuantityName: 'energy_E',
      // E = -p_t, so the conserved quantity is -(g_mu_nu xi^mu t^nu) for xi = d/dt.
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
 * The spherical chart's relation to the auxiliary Cartesian visualization axes.
 *
 * The direction uses the tangent's orthonormal spatial components as measured by a
 * static observer, obtained for this diagonal metric as v^(i) = sqrt(g_ii) k^i:
 *
 *   v_radial    = k^r / sqrt(f)
 *   v_polar     = r k^theta
 *   v_azimuthal = r sin(theta) k^phi
 *
 * combined with the usual spherical basis vectors. At large r the static frame is
 * asymptotically inertial, so this is the direction the ray is really travelling.
 */
function geometryFor(M: number): ChartGeometry {
  return {
    spatialRadius: (x: Vec4): number => x[1],

    toCartesianPosition: (x: Vec4): CartesianVec3 => {
      const [, r, theta, phi] = x;
      const sinTheta = Math.sin(theta);
      return [r * sinTheta * Math.cos(phi), r * sinTheta * Math.sin(phi), r * Math.cos(theta)];
    },

    toCartesianDirection: (x: Vec4, tangent: Vec4): CartesianVec3 => {
      const [, r, theta, phi] = x;
      const f = 1 - (2 * M) / r;
      const sinTheta = Math.sin(theta);
      const cosTheta = Math.cos(theta);
      const cosPhi = Math.cos(phi);
      const sinPhi = Math.sin(phi);

      // Orthonormal spatial components in the static frame.
      const vRadial = f > 0 ? tangent[1] / Math.sqrt(f) : tangent[1];
      const vPolar = r * tangent[2];
      const vAzimuthal = r * sinTheta * tangent[3];

      // Spherical basis vectors on the Cartesian axes.
      const rHat: CartesianVec3 = [sinTheta * cosPhi, sinTheta * sinPhi, cosTheta];
      const thetaHat: CartesianVec3 = [cosTheta * cosPhi, cosTheta * sinPhi, -sinTheta];
      const phiHat: CartesianVec3 = [-sinPhi, cosPhi, 0];

      return [
        vRadial * rHat[0] + vPolar * thetaHat[0] + vAzimuthal * phiHat[0],
        vRadial * rHat[1] + vPolar * thetaHat[1] + vAzimuthal * phiHat[1],
        vRadial * rHat[2] + vPolar * thetaHat[2] + vAzimuthal * phiHat[2],
      ];
    },
  };
}

class SchwarzschildSpacetime implements SpacetimeModel {
  readonly id: string;
  readonly displayName: string;
  readonly classification = 'exact-analytical' as const;
  readonly chart = SCHWARZSCHILD_CHART;
  readonly conventions = CONVENTIONS;
  readonly parameters: Readonly<Record<string, number>>;
  readonly killingVectors = killingVectors();
  readonly symmetries = SCHWARZSCHILD_SYMMETRIES;
  readonly geometry: ChartGeometry;
  readonly description =
    'Exact vacuum solution of the Einstein field equations, spherically symmetric and ' +
    'static. Ricci-flat everywhere it is defined. Integrating geodesics here evaluates ' +
    'the consequences of a specified geometry, not a dynamical solution of the field ' +
    'equations (CLAUDE.md §1.2).';

  /** The mass parameter M, in geometric units. */
  readonly M: number;

  constructor(M: number) {
    if (!(M > 0) || !Number.isFinite(M)) {
      throw new RangeError(
        `Schwarzschild: the mass parameter must be finite and positive, received ${M}. ` +
          'For the flat-space limit use a small positive M and check convergence, or use ' +
          'the Minkowski model directly.',
      );
    }
    this.M = M;
    this.id = 'schwarzschild';
    this.displayName = `Schwarzschild spacetime (M = ${M})`;
    this.parameters = Object.freeze({ M });
    this.geometry = geometryFor(M);
  }

  /** The metric function f = 1 - 2M/r. Negative inside the horizon. */
  lapseFunction(r: number): number {
    return 1 - (2 * this.M) / r;
  }

  metricAt(x: Vec4): MetricTensor {
    const r = x[1];
    const theta = x[2];
    const f = this.lapseFunction(r);

    if (!(f > 0)) {
      throw new RangeError(
        `Schwarzschild.metricAt: f = 1 - 2M/r = ${f} at r = ${r}, so the event is at or ` +
          `inside the horizon r = ${horizonRadius(this.M)}. This chart covers the exterior ` +
          'only; continuing through the horizon needs horizon-penetrating coordinates.',
      );
    }

    const sinTheta = Math.sin(theta);
    const g = new Float64Array(16);
    const gInv = new Float64Array(16);

    g[0] = -f;
    g[5] = 1 / f;
    g[10] = r * r;
    g[15] = r * r * sinTheta * sinTheta;

    // Diagonal metric: the inverse is the componentwise reciprocal.
    gInv[0] = -1 / f;
    gInv[5] = f;
    gInv[10] = 1 / (r * r);
    gInv[15] = 1 / (r * r * sinTheta * sinTheta);

    return new MetricTensor(g, gInv);
  }

  christoffelAt(x: Vec4): ChristoffelSymbols {
    const out = new Float64Array(64);
    this.christoffelInto(x, out);
    return new ChristoffelSymbols(out);
  }

  christoffelInto(x: Vec4, out: Float64Array): void {
    if (out.length !== 64) {
      throw new RangeError('christoffelInto expects a 64-entry buffer.');
    }
    out.fill(0);

    const r = x[1];
    const theta = x[2];
    const M = this.M;
    const f = 1 - (2 * M) / r;
    const sinTheta = Math.sin(theta);
    const cosTheta = Math.cos(theta);

    // Index layout: mu * 16 + alpha * 4 + beta, with (t, r, theta, phi) = (0, 1, 2, 3).
    const m_over_r2f = M / (r * r * f);

    // Gamma^t_{t r} = Gamma^t_{r t} = M / (r^2 f)
    out[0 * 16 + 0 * 4 + 1] = m_over_r2f;
    out[0 * 16 + 1 * 4 + 0] = m_over_r2f;

    // Gamma^r_{t t} = M f / r^2
    out[1 * 16 + 0 * 4 + 0] = (M * f) / (r * r);
    // Gamma^r_{r r} = -M / (r^2 f)
    out[1 * 16 + 1 * 4 + 1] = -m_over_r2f;
    // Gamma^r_{theta theta} = -r f
    out[1 * 16 + 2 * 4 + 2] = -r * f;
    // Gamma^r_{phi phi} = -r f sin^2(theta)
    out[1 * 16 + 3 * 4 + 3] = -r * f * sinTheta * sinTheta;

    // Gamma^theta_{r theta} = Gamma^theta_{theta r} = 1 / r
    out[2 * 16 + 1 * 4 + 2] = 1 / r;
    out[2 * 16 + 2 * 4 + 1] = 1 / r;
    // Gamma^theta_{phi phi} = -sin(theta) cos(theta)
    out[2 * 16 + 3 * 4 + 3] = -sinTheta * cosTheta;

    // Gamma^phi_{r phi} = Gamma^phi_{phi r} = 1 / r
    out[3 * 16 + 1 * 4 + 3] = 1 / r;
    out[3 * 16 + 3 * 4 + 1] = 1 / r;
    // Gamma^phi_{theta phi} = Gamma^phi_{phi theta} = cot(theta)
    const cotTheta = cosTheta / sinTheta;
    out[3 * 16 + 2 * 4 + 3] = cotTheta;
    out[3 * 16 + 3 * 4 + 2] = cotTheta;
  }

  /**
   * Whether an event lies in the static exterior chart.
   *
   * Detected from the metric components actually in use: f = 1 - 2M/r appears directly
   * in g_t_t and as the reciprocal in g_r_r, and sin(theta) appears in g_phi_phi. The
   * metric determinant is deliberately not consulted — CLAUDE.md §6.1 is explicit that
   * the horizon must not be found that way, and for this metric the determinant is
   * -r^4 sin^2(theta), which is perfectly finite and non-zero at r = 2M.
   */
  domainCheck(x: Vec4): DomainStatus {
    const r = x[1];
    const theta = x[2];

    if (!Number.isFinite(r) || !Number.isFinite(theta)) {
      return {
        inDomain: false,
        code: 'outside-chart',
        reason: `Non-finite coordinate: r = ${r}, theta = ${theta}.`,
      };
    }

    const f = 1 - (2 * this.M) / r;
    if (!(f > 0)) {
      return {
        inDomain: false,
        code: 'coordinate-breakdown',
        reason:
          `f = 1 - 2M/r = ${f.toExponential(3)} at r = ${r}, so the event is at or inside ` +
          `the horizon r = ${horizonRadius(this.M)}. g_t_t vanishes and g_r_r diverges ` +
          'there: the Schwarzschild chart breaks down. This is a coordinate singularity, ' +
          `not a curvature singularity — the Kretschmann scalar is ` +
          `${kretschmann(this.M, r).toExponential(3)} and finite.`,
      };
    }

    const sinTheta = Math.abs(Math.sin(theta));
    if (sinTheta < POLAR_AXIS_SIN_THETA_FLOOR) {
      return {
        inDomain: false,
        code: 'coordinate-breakdown',
        reason:
          `sin(theta) = ${sinTheta.toExponential(3)} is below the polar-axis floor ` +
          `${POLAR_AXIS_SIN_THETA_FLOOR}. g_phi_phi vanishes and Gamma^phi_{theta phi} = ` +
          'cot(theta) diverges on the axis. This is an artifact of spherical coordinates, ' +
          'not a physical boundary; the spacetime itself is perfectly regular there.',
      };
    }

    return IN_DOMAIN;
  }
}

export type SchwarzschildModel = SchwarzschildSpacetime;

/** Build a Schwarzschild model with mass parameter `M` in geometric units. */
export function schwarzschild(M = 1): SchwarzschildModel {
  return new SchwarzschildSpacetime(M);
}
