import type { Vec4 } from './indices.js';

/**
 * What kind of worldline a phase-space state describes.
 *
 * CLAUDE.md §3 requires the engine to distinguish massive test-particle worldlines,
 * photon / null-ray trajectories, and observer worldlines. The first two differ in
 * their normalization and in their integration parameter; observer worldlines are
 * timelike and are modelled in the observer layer.
 */
export type WorldlineKind = 'null' | 'timelike';

/**
 * A point in geodesic phase space: an event plus the tangent to the worldline there.
 *
 * This is the state the first-order geodesic equations advance (ROADMAP.md 1.3):
 *
 *   dx^mu / dparam = t^mu
 *   dt^mu / dparam = -Gamma^mu_{alpha beta} t^alpha t^beta
 *
 * `tangent` is deliberately named neutrally. Reach for it through
 * `null_wavevector_k` or `four_velocity_u`, which check the worldline kind first:
 * CLAUDE.md §3 and §24 both forbid calling a null vector a four-velocity, and a
 * neutral field name with checked accessors makes that mistake fail loudly rather
 * than silently produce a plausible number.
 */
export interface PhaseSpaceState {
  readonly kind: WorldlineKind;
  /** Affine parameter lambda for null worldlines; proper time tau for timelike ones. */
  readonly parameter: number;
  /** The event x^mu, in the active chart's coordinate ordering. */
  readonly position_x: Vec4;
  /** The tangent: k^mu if `kind` is null, u^mu if timelike. */
  readonly tangent: Vec4;
}

/** The name of the integration parameter for a worldline kind (CLAUDE.md §2). */
export function parameterName(kind: WorldlineKind): 'lambda' | 'tau' {
  return kind === 'null' ? 'lambda' : 'tau';
}

/**
 * The null wavevector k^mu = dx^mu/dlambda.
 *
 * Throws for a timelike state. CLAUDE.md §3: "Do not call k^mu a four-velocity."
 */
export function null_wavevector_k(state: PhaseSpaceState): Vec4 {
  if (state.kind !== 'null') {
    throw new TypeError(
      'null_wavevector_k: this state describes a timelike worldline. Its tangent is a ' +
        'four-velocity u^mu with g_mu_nu u^mu u^nu = -1, not a null wavevector. ' +
        'Use four_velocity_u instead.',
    );
  }
  return state.tangent;
}

/**
 * The four-velocity u^mu = dx^mu/dtau.
 *
 * Throws for a null state: a null tangent is not a four-velocity and cannot be
 * normalized to -1.
 */
export function four_velocity_u(state: PhaseSpaceState): Vec4 {
  if (state.kind !== 'timelike') {
    throw new TypeError(
      'four_velocity_u: this state describes a null worldline. Its tangent is a null ' +
        'wavevector k^mu with g_mu_nu k^mu k^nu = 0, and a null vector is not a ' +
        'four-velocity. Use null_wavevector_k instead.',
    );
  }
  return state.tangent;
}

export function nullState(position_x: Vec4, null_wavevector_k: Vec4, lambda = 0): PhaseSpaceState {
  return { kind: 'null', parameter: lambda, position_x, tangent: null_wavevector_k };
}

export function timelikeState(position_x: Vec4, four_velocity_u: Vec4, tau = 0): PhaseSpaceState {
  return { kind: 'timelike', parameter: tau, position_x, tangent: four_velocity_u };
}
