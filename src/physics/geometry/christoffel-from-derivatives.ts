/**
 * Christoffel symbols from a metric and its analytic derivatives.
 *
 * CLAUDE.md §2 defines
 *
 *   Gamma^mu_{alpha beta}
 *     = 1/2 g^{mu sigma} ( d_alpha g_{sigma beta} + d_beta g_{sigma alpha}
 *                          - d_sigma g_{alpha beta} )
 *
 * and this evaluates exactly that. Schwarzschild's nine non-zero symbols are short
 * enough to write out by hand; Kerr's are not, and hand-transcribing two dozen rational
 * functions of (r, theta) is a reliable way to introduce a silent error. Supplying the
 * metric derivatives in closed form and contracting them here keeps the hand-written
 * input small enough to check term by term, and the result is still analytic — no
 * differencing, no step size, no truncation error.
 *
 * Index layout matches the rest of the engine: `dg` holds d_alpha g_{mu nu} at
 * alpha * 16 + mu * 4 + nu, `gInv` holds g^{mu nu} at mu * 4 + nu, and `out` receives
 * Gamma^mu_{alpha beta} at mu * 16 + alpha * 4 + beta.
 */
export function christoffelFromMetricDerivatives(
  gInv: Float64Array,
  dg: Float64Array,
  out: Float64Array,
): void {
  if (gInv.length !== 16) throw new RangeError('christoffelFromMetricDerivatives: gInv needs 16 entries.');
  if (dg.length !== 64) throw new RangeError('christoffelFromMetricDerivatives: dg needs 64 entries.');
  if (out.length !== 64) throw new RangeError('christoffelFromMetricDerivatives: out needs 64 entries.');

  for (let mu = 0; mu < 4; mu += 1) {
    for (let alpha = 0; alpha < 4; alpha += 1) {
      for (let beta = alpha; beta < 4; beta += 1) {
        let sum = 0;
        for (let sigma = 0; sigma < 4; sigma += 1) {
          const inverse = gInv[mu * 4 + sigma];
          if (inverse === 0) continue;
          sum +=
            inverse *
            (dg[alpha * 16 + sigma * 4 + beta] +
              dg[beta * 16 + sigma * 4 + alpha] -
              dg[sigma * 16 + alpha * 4 + beta]);
        }
        const value = 0.5 * sum;
        // Symmetric in the lower pair for a Levi-Civita connection, so it is computed
        // once and written to both slots rather than computed twice.
        out[mu * 16 + alpha * 4 + beta] = value;
        out[mu * 16 + beta * 4 + alpha] = value;
      }
    }
  }
}
