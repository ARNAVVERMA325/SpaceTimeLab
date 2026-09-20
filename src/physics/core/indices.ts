/** Spacetime index running 0..3, in the declared coordinate ordering. */
export type Index4 = 0 | 1 | 2 | 3;

/** All four spacetime indices, for iteration. */
export const INDICES: readonly Index4[] = Object.freeze([0, 1, 2, 3]) as readonly Index4[];

/** Four components in the declared coordinate ordering of the active chart. */
export type Vec4 = readonly [number, number, number, number];

/** Row-major offset into a rank-2 4x4 array stored flat. */
export function idx2(mu: Index4, nu: Index4): number {
  return mu * 4 + nu;
}

/** Row-major offset into a rank-3 4x4x4 array stored flat. */
export function idx3(mu: Index4, alpha: Index4, beta: Index4): number {
  return mu * 16 + alpha * 4 + beta;
}

export function vec4(a: number, b: number, c: number, d: number): Vec4 {
  return [a, b, c, d];
}

/** Copy a Vec4 out of a flat buffer starting at `offset`. */
export function readVec4(buffer: Float64Array, offset: number): Vec4 {
  return [buffer[offset], buffer[offset + 1], buffer[offset + 2], buffer[offset + 3]];
}
