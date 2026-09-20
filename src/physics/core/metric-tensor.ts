import { idx2, type Index4, type Vec4 } from './indices.js';

/**
 * The metric tensor at one event: covariant components and their inverse.
 *
 * Naming follows CLAUDE.md §19 exactly. `g_inv_mu_nu` holds g^{mu nu}; the name is the
 * one the project rules mandate, so it is used verbatim rather than "corrected" here.
 *
 * Components are stored flat (row-major, 16 entries) so the hot integration path can
 * work on typed arrays without allocating per step.
 */
export class MetricTensor {
  readonly g_mu_nu: Float64Array;
  readonly g_inv_mu_nu: Float64Array;

  constructor(g_mu_nu: Float64Array, g_inv_mu_nu: Float64Array) {
    if (g_mu_nu.length !== 16 || g_inv_mu_nu.length !== 16) {
      throw new RangeError('MetricTensor expects 16 components for each of g_mu_nu and g_inv_mu_nu.');
    }
    this.g_mu_nu = g_mu_nu;
    this.g_inv_mu_nu = g_inv_mu_nu;
  }

  /** g_{mu nu} */
  g(mu: Index4, nu: Index4): number {
    return this.g_mu_nu[idx2(mu, nu)];
  }

  /** g^{mu nu} */
  gInv(mu: Index4, nu: Index4): number {
    return this.g_inv_mu_nu[idx2(mu, nu)];
  }

  /** Build a metric from a diagonal, inverting componentwise. */
  static diagonal(d0: number, d1: number, d2: number, d3: number): MetricTensor {
    const g = new Float64Array(16);
    const gInv = new Float64Array(16);
    const d = [d0, d1, d2, d3];
    for (let i = 0; i < 4; i += 1) {
      if (d[i] === 0) {
        throw new RangeError(`MetricTensor.diagonal: component ${i} is zero; the metric is degenerate.`);
      }
      g[i * 4 + i] = d[i];
      gInv[i * 4 + i] = 1 / d[i];
    }
    return new MetricTensor(g, gInv);
  }
}

/**
 * The scalar g_mu_nu a^mu b^nu, for two contravariant four-vectors.
 *
 * This is the contraction behind every normalization check in CLAUDE.md §16.
 */
export function contractLower(metric: MetricTensor, a: Vec4, b: Vec4): number {
  let sum = 0;
  for (let mu = 0; mu < 4; mu += 1) {
    for (let nu = 0; nu < 4; nu += 1) {
      const g = metric.g_mu_nu[mu * 4 + nu];
      if (g !== 0) sum += g * a[mu] * b[nu];
    }
  }
  return sum;
}

/** The scalar g^{mu nu} a_mu b_nu, for two covariant four-vectors. */
export function contractUpper(metric: MetricTensor, a: Vec4, b: Vec4): number {
  let sum = 0;
  for (let mu = 0; mu < 4; mu += 1) {
    for (let nu = 0; nu < 4; nu += 1) {
      const gInv = metric.g_inv_mu_nu[mu * 4 + nu];
      if (gInv !== 0) sum += gInv * a[mu] * b[nu];
    }
  }
  return sum;
}
