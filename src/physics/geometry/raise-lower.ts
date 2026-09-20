import type { MetricTensor } from '../core/metric-tensor.js';
import { contravariant, covariant, requireVariance, type FourVector } from '../core/four-vector.js';
import type { Vec4 } from '../core/indices.js';

/**
 * Differential-geometry layer (CLAUDE.md §20): moving indices with the metric.
 *
 * Both routines check the input's index position first, because raising an already-upper
 * index is the kind of silent error CLAUDE.md §24 warns about.
 */

/** V_mu = g_{mu nu} V^nu. */
export function lowerIndex(metric: MetricTensor, v: FourVector): FourVector {
  const up = requireVariance(v, 'contravariant', 'lowerIndex');
  const out: [number, number, number, number] = [0, 0, 0, 0];
  for (let mu = 0; mu < 4; mu += 1) {
    let sum = 0;
    for (let nu = 0; nu < 4; nu += 1) {
      sum += metric.g_mu_nu[mu * 4 + nu] * up[nu];
    }
    out[mu] = sum;
  }
  return covariant(out);
}

/** V^mu = g^{mu nu} V_nu. */
export function raiseIndex(metric: MetricTensor, v: FourVector): FourVector {
  const down = requireVariance(v, 'covariant', 'raiseIndex');
  const out: [number, number, number, number] = [0, 0, 0, 0];
  for (let mu = 0; mu < 4; mu += 1) {
    let sum = 0;
    for (let nu = 0; nu < 4; nu += 1) {
      sum += metric.g_inv_mu_nu[mu * 4 + nu] * down[nu];
    }
    out[mu] = sum;
  }
  return contravariant(out);
}

/**
 * The covariant four-momentum p_mu = g_{mu nu} t^nu for a geodesic tangent t^nu.
 *
 * This is the object the conserved quantities of CLAUDE.md §16 are built from:
 * E = -p_t for a stationary spacetime, L_z = p_phi for an axisymmetric one.
 */
export function four_momentum_p(metric: MetricTensor, tangent: Vec4): Vec4 {
  return lowerIndex(metric, contravariant(tangent)).components;
}
