import type { Vec4 } from '../core/indices.js';
import type { SpacetimeModel } from '../spacetimes/spacetime-model.js';
import { POSITION_OFFSET, STATE_DIM, TANGENT_OFFSET } from './state-vector.js';

/**
 * A first-order ODE right-hand side, written in place.
 *
 * `dydp` is filled with dy/dparam at parameter `p`. In place because the integrator
 * calls this several times per step and per ray; allocating here would dominate the
 * cost of a full image.
 */
export type DerivativeFn = (p: number, y: Float64Array, dydp: Float64Array) => void;

/**
 * The geodesic equation in first-order form (CLAUDE.md §2, ROADMAP.md 1.3).
 *
 * The second-order geodesic equation
 *
 *   d^2 x^mu / dparam^2 + Gamma^mu_{alpha beta} (dx^alpha/dparam)(dx^beta/dparam) = 0
 *
 * is split into the pair
 *
 *   dx^mu / dparam = t^mu
 *   dt^mu / dparam = -Gamma^mu_{alpha beta} t^alpha t^beta
 *
 * with t^mu the tangent: k^mu against affine parameter lambda for null worldlines,
 * u^mu against proper time tau for timelike ones. The equation itself is identical in
 * both cases; only the normalization of the tangent and the name of the parameter
 * differ, which is why the worldline kind lives on `PhaseSpaceState` rather than here.
 *
 * The returned closure owns a Christoffel scratch buffer, so it is stateful and must
 * not be shared across concurrent integrations.
 */
export function geodesicDerivative(model: SpacetimeModel): DerivativeFn {
  const christoffel = new Float64Array(64);
  const position: [number, number, number, number] = [0, 0, 0, 0];

  return function geodesicRhs(_p: number, y: Float64Array, dydp: Float64Array): void {
    if (y.length !== STATE_DIM || dydp.length !== STATE_DIM) {
      throw new RangeError(`geodesicDerivative expects ${STATE_DIM}-entry state buffers.`);
    }

    position[0] = y[POSITION_OFFSET];
    position[1] = y[POSITION_OFFSET + 1];
    position[2] = y[POSITION_OFFSET + 2];
    position[3] = y[POSITION_OFFSET + 3];

    model.christoffelInto(position as Vec4, christoffel);

    // dx^mu / dparam = t^mu
    dydp[POSITION_OFFSET] = y[TANGENT_OFFSET];
    dydp[POSITION_OFFSET + 1] = y[TANGENT_OFFSET + 1];
    dydp[POSITION_OFFSET + 2] = y[TANGENT_OFFSET + 2];
    dydp[POSITION_OFFSET + 3] = y[TANGENT_OFFSET + 3];

    // dt^mu / dparam = -Gamma^mu_{alpha beta} t^alpha t^beta
    for (let mu = 0; mu < 4; mu += 1) {
      let acc = 0;
      const base = mu * 16;
      for (let alpha = 0; alpha < 4; alpha += 1) {
        const tAlpha = y[TANGENT_OFFSET + alpha];
        if (tAlpha === 0) continue;
        const rowBase = base + alpha * 4;
        for (let beta = 0; beta < 4; beta += 1) {
          const gamma = christoffel[rowBase + beta];
          if (gamma !== 0) acc += gamma * tAlpha * y[TANGENT_OFFSET + beta];
        }
      }
      dydp[TANGENT_OFFSET + mu] = -acc;
    }
  };
}
