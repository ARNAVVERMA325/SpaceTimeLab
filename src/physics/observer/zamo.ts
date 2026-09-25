import type { Vec4 } from '../core/indices.js';
import type { SpacetimeModel } from '../spacetimes/spacetime-model.js';
import type { Observer } from './observer.js';
import { Tetrad } from './tetrad.js';

/**
 * The zero-angular-momentum observer (ZAMO), also called the locally non-rotating frame.
 *
 * Bardeen, Press & Teukolsky, ApJ 178, 347 (1972). In a stationary axisymmetric metric
 * with a g_{t phi} cross term there is no diagonal static tetrad, and inside the
 * ergosphere there is no static observer at all: the timelike Killing vector has turned
 * spacelike, so no worldline of fixed (r, theta, phi) is physical. The ZAMO exists
 * everywhere outside the horizon. It is the observer whose angular momentum vanishes,
 *
 *   p_phi = g_{phi mu} u^mu = 0,
 *
 * which fixes its coordinate angular velocity to
 *
 *   omega = -g_{t phi} / g_{phi phi} = 2 M a r / A.
 *
 * This observer is dragged around the hole even though it is "not rotating" by its own
 * local measurements. That is not a paradox and not a choice of convenience: it is what
 * frame dragging means. A ZAMO is *not* a geodesic observer — it accelerates — and its
 * frame is tied to the chart's time slicing, so quantities it measures are that
 * observer's measurements and not invariants (CLAUDE.md §1.3, §4).
 *
 * The legs are
 *
 *   e_(0) = (1/alpha) (d_t + omega d_phi),   alpha = sqrt(Sigma Delta / A)
 *   e_(1) = sqrt(Delta / Sigma) d_r
 *   e_(2) = (1/sqrt(Sigma)) d_theta
 *   e_(3) = (1/varpi) d_phi,                 varpi = sqrt(A / Sigma) sin(theta)
 *
 * built here from the metric components themselves rather than from those closed forms,
 * so the construction works for any stationary axisymmetric metric whose only cross term
 * is g_{t phi} — Kerr in Boyer-Lindquist today, another such chart later.
 */
export function zamoTetrad(model: SpacetimeModel, position_x: Vec4): Tetrad {
  const g = new Float64Array(16);
  model.metricInto(position_x, g);

  const g_tt = g[0];
  const g_tphi = g[3];
  const g_rr = g[5];
  const g_thth = g[10];
  const g_phiphi = g[15];

  for (let mu = 0; mu < 4; mu += 1) {
    for (let nu = 0; nu < 4; nu += 1) {
      if (mu === nu) continue;
      const isTimeAzimuth = (mu === 0 && nu === 3) || (mu === 3 && nu === 0);
      if (!isTimeAzimuth && g[mu * 4 + nu] !== 0) {
        throw new RangeError(
          `zamoTetrad: the metric has an off-diagonal component g_${mu}_${nu} = ` +
            `${g[mu * 4 + nu]} other than g_t_phi. This construction handles a stationary ` +
            'axisymmetric chart whose only cross term is between t and phi.',
        );
      }
    }
  }
  if (!(g_phiphi > 0) || !(g_rr > 0) || !(g_thth > 0)) {
    throw new RangeError(
      `zamoTetrad: the spatial metric has degenerated at r = ${position_x[1]}, ` +
        `theta = ${position_x[2]} (g_rr = ${g_rr}, g_theta_theta = ${g_thth}, ` +
        `g_phi_phi = ${g_phiphi}). On the polar axis or at the horizon this chart does not ` +
        'support an orthonormal frame.',
    );
  }

  const omega = -g_tphi / g_phiphi;
  // g(u, u) = g_tt - g_t_phi^2 / g_phi_phi, which must be negative for a timelike observer.
  const normSquared = g_tt - (g_tphi * g_tphi) / g_phiphi;
  if (!(normSquared < 0)) {
    throw new RangeError(
      `zamoTetrad: g_tt - g_t_phi^2 / g_phi_phi = ${normSquared} is not negative at ` +
        `r = ${position_x[1]}. The observer's worldline would not be timelike; this happens ` +
        'at and inside the horizon, which this chart does not cover.',
    );
  }
  const lapse = Math.sqrt(-normSquared);

  const components = new Float64Array(16);
  components[0] = 1 / lapse;
  components[3] = omega / lapse;
  components[4 + 1] = 1 / Math.sqrt(g_rr);
  components[8 + 2] = 1 / Math.sqrt(g_thth);
  components[12 + 3] = 1 / Math.sqrt(g_phiphi);
  return new Tetrad(components);
}

/**
 * The lapse alpha = sqrt(-(g_tt - g_t_phi^2 / g_phi_phi)) of the ZAMO slicing.
 *
 * The rate of the observer's proper time against coordinate time t, which is normalized
 * to a static observer at infinity. CLAUDE.md §18 applies: this is a statement about one
 * family of observers in one slicing, not a universal rate at which "time passes".
 */
export function zamoLapse(model: SpacetimeModel, position_x: Vec4): number {
  const g = new Float64Array(16);
  model.metricInto(position_x, g);
  return Math.sqrt(-(g[0] - (g[3] * g[3]) / g[15]));
}

/** A zero-angular-momentum observer at an event, as a physical observer. */
export function zamoObserver(model: SpacetimeModel, position_x: Vec4, id = 'zamo'): Observer {
  const tetrad = zamoTetrad(model, position_x);
  const u = tetrad.four_velocity_u();
  const omega = u[3] / u[0];
  return {
    id,
    displayName: `Zero-angular-momentum observer at r = ${position_x[1]}`,
    kind: 'physical-observer',
    position_x,
    tetrad,
    description:
      'The locally non-rotating frame of Bardeen, Press & Teukolsky (1972): the observer ' +
      'whose angular momentum p_phi vanishes. It is carried around the hole at coordinate ' +
      `angular velocity omega = ${omega.toExponential(3)} by frame dragging, and it ` +
      'accelerates, so it is not a freely falling frame. Unlike a static observer it exists ' +
      'inside the ergosphere, where no one can hold station at fixed phi.',
  };
}
