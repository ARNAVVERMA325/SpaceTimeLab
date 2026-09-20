import type { Vec4 } from '../core/indices.js';
import type { PhaseSpaceState, WorldlineKind } from '../core/phase-space.js';

/**
 * Packed representation of a phase-space state for the integrator core.
 *
 * Layout: [ x^0, x^1, x^2, x^3, t^0, t^1, t^2, t^3 ] where t^mu is the tangent
 * (k^mu for null worldlines, u^mu for timelike ones).
 *
 * The integrators work on this flat buffer so that a step allocates nothing.
 * `PhaseSpaceState` remains the readable object the rest of the engine passes around.
 */
export const STATE_DIM = 8;
export const POSITION_OFFSET = 0;
export const TANGENT_OFFSET = 4;

export function packState(state: PhaseSpaceState, out = new Float64Array(STATE_DIM)): Float64Array {
  if (out.length !== STATE_DIM) {
    throw new RangeError(`packState expects a ${STATE_DIM}-entry buffer.`);
  }
  for (let mu = 0; mu < 4; mu += 1) {
    out[POSITION_OFFSET + mu] = state.position_x[mu];
    out[TANGENT_OFFSET + mu] = state.tangent[mu];
  }
  return out;
}

export function unpackState(y: Float64Array, kind: WorldlineKind, parameter: number): PhaseSpaceState {
  if (y.length !== STATE_DIM) {
    throw new RangeError(`unpackState expects a ${STATE_DIM}-entry buffer.`);
  }
  const position_x: Vec4 = [y[0], y[1], y[2], y[3]];
  const tangent: Vec4 = [y[4], y[5], y[6], y[7]];
  return { kind, parameter, position_x, tangent };
}

export function positionOf(y: Float64Array): Vec4 {
  return [y[0], y[1], y[2], y[3]];
}

export function tangentOf(y: Float64Array): Vec4 {
  return [y[4], y[5], y[6], y[7]];
}
