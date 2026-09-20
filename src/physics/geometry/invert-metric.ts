/**
 * Differential-geometry layer (CLAUDE.md §20): inverting g_mu_nu to obtain g^{mu nu}.
 *
 * Gauss-Jordan elimination with partial pivoting on the 4x4 system. A general routine
 * is used rather than a closed-form determinant expression because later models
 * (Kerr in Boyer-Lindquist, Kerr-Schild) have off-diagonal metrics, and partial
 * pivoting keeps the elimination stable when a diagonal entry is small.
 *
 * Note that a singular result here means the metric *components in this chart* have
 * degenerated, which is a statement about the chart. CLAUDE.md §6.1 forbids using a
 * determinant test as a horizon detector; the horizon check belongs to each model's
 * `domainCheck`.
 */
export function invertMetric4(g: Float64Array): Float64Array {
  if (g.length !== 16) {
    throw new RangeError('invertMetric4 expects 16 components.');
  }

  // Augmented [g | I], 4 rows of 8.
  const a = new Float64Array(32);
  for (let row = 0; row < 4; row += 1) {
    for (let col = 0; col < 4; col += 1) {
      a[row * 8 + col] = g[row * 4 + col];
    }
    a[row * 8 + 4 + row] = 1;
  }

  for (let col = 0; col < 4; col += 1) {
    let pivotRow = col;
    let pivotMag = Math.abs(a[col * 8 + col]);
    for (let row = col + 1; row < 4; row += 1) {
      const mag = Math.abs(a[row * 8 + col]);
      if (mag > pivotMag) {
        pivotMag = mag;
        pivotRow = row;
      }
    }

    if (pivotMag === 0 || !Number.isFinite(pivotMag)) {
      throw new RangeError(
        `invertMetric4: the metric components are singular or non-finite at column ${col}; ` +
          'g_mu_nu cannot be inverted in this chart at this event.',
      );
    }

    if (pivotRow !== col) {
      for (let k = 0; k < 8; k += 1) {
        const tmp = a[col * 8 + k];
        a[col * 8 + k] = a[pivotRow * 8 + k];
        a[pivotRow * 8 + k] = tmp;
      }
    }

    const pivot = a[col * 8 + col];
    for (let k = 0; k < 8; k += 1) a[col * 8 + k] /= pivot;

    for (let row = 0; row < 4; row += 1) {
      if (row === col) continue;
      const factor = a[row * 8 + col];
      if (factor === 0) continue;
      for (let k = 0; k < 8; k += 1) {
        a[row * 8 + k] -= factor * a[col * 8 + k];
      }
    }
  }

  const inverse = new Float64Array(16);
  for (let row = 0; row < 4; row += 1) {
    for (let col = 0; col < 4; col += 1) {
      inverse[row * 4 + col] = a[row * 8 + 4 + col];
    }
  }
  return inverse;
}

/**
 * Largest |g^{mu sigma} g_{sigma nu} - delta^mu_nu| over all components.
 *
 * The direct check that an inverse really is one, used by the validation layer.
 */
export function inverseResidual(g: Float64Array, gInv: Float64Array): number {
  let worst = 0;
  for (let mu = 0; mu < 4; mu += 1) {
    for (let nu = 0; nu < 4; nu += 1) {
      let sum = 0;
      for (let sigma = 0; sigma < 4; sigma += 1) {
        sum += gInv[mu * 4 + sigma] * g[sigma * 4 + nu];
      }
      const delta = mu === nu ? 1 : 0;
      const residual = Math.abs(sum - delta);
      if (residual > worst) worst = residual;
    }
  }
  return worst;
}
