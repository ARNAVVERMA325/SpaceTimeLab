import type { Vec4 } from '../core/indices.js';
import type { PhaseSpaceState } from '../core/phase-space.js';
import { parameterName } from '../core/phase-space.js';
import type { SpacetimeModel } from '../spacetimes/spacetime-model.js';
import { firstNonFiniteIndex } from '../validation/numeric-health.js';
import { geodesicDerivative, type DerivativeFn } from './geodesic-system.js';
import type { Integrator } from './integrators/integrator.js';
import { packState, STATE_DIM, unpackState } from './state-vector.js';

/**
 * Geodesic / dynamical layer driver (CLAUDE.md §20).
 *
 * The driver owns termination policy, step bookkeeping and numerical-health checks.
 * The spacetime model owns the geometry and the integrator owns the stepping; neither
 * decides when to stop.
 */

export type TerminationReason =
  /** Reached the requested parameter bound. */
  | 'parameter-limit'
  /** Hit the step budget before any other condition. */
  | 'step-limit'
  /** The model reported the event is outside its chart's valid domain. */
  | 'domain-exit'
  /** A caller-supplied terminator fired, e.g. a ray leaving the computational domain. */
  | 'terminator'
  /** Adaptive control drove the step below the floor without accepting. */
  | 'step-underflow'
  /** A NaN or infinity appeared in the state. */
  | 'numerical-failure';

/**
 * Diagnostics for a failed or truncated integration (CLAUDE.md §17).
 *
 * CLAUDE.md §17 requires that a numerical failure expose the integrator, model,
 * coordinate system, current parameter, step size and failed quantity, rather than
 * being silently replaced by a plausible-looking result.
 */
export interface IntegrationDiagnostics {
  readonly integrator: string;
  readonly model: string;
  readonly coordinateSystem: string;
  readonly parameterName: 'lambda' | 'tau';
  readonly parameterValue: number;
  readonly stepSize: number;
  readonly failedQuantity: string;
  readonly detail: string;
}

export interface IntegrationLimits {
  /** Upper bound on |parameter - parameter_0|. */
  readonly parameterMax: number;
  readonly maxSteps: number;
  /** Initial step. For a fixed-step integrator this is *the* step. */
  readonly initialStep: number;
  /** Adaptive control gives up below this magnitude. */
  readonly minStep: number;
  /** Optional cap on the adaptive step. */
  readonly maxStep?: number;
}

export const DEFAULT_LIMITS: IntegrationLimits = Object.freeze({
  parameterMax: 1e3,
  maxSteps: 100_000,
  initialStep: 1e-2,
  minStep: 1e-12,
});

/**
 * A caller-supplied stopping condition evaluated after each accepted step.
 *
 * Receives the tangent as well as the position, because a physically exact stopping
 * condition often needs both: a photon inside the Schwarzschild photon sphere is
 * captured if and only if it is moving inward, which no position alone can express.
 */
export type Terminator = (position_x: Vec4, tangent: Vec4, parameter: number) => boolean;

export interface IntegrateOptions {
  readonly model: SpacetimeModel;
  readonly integrator: Integrator;
  readonly initial: PhaseSpaceState;
  readonly limits?: Partial<IntegrationLimits>;
  readonly terminator?: Terminator;
  /** Reuse a derivative closure across rays; one is built per call otherwise. */
  readonly derivative?: DerivativeFn;
  /** Record every accepted state. Off by default: a full image would exhaust memory. */
  readonly recordPath?: boolean;
}

export interface IntegrationResult {
  readonly final: PhaseSpaceState;
  readonly reason: TerminationReason;
  readonly steps: number;
  readonly rejectedSteps: number;
  readonly finalStep: number;
  /** Largest scaled local error norm seen on an accepted step; NaN for fixed-step methods. */
  readonly maxErrorNorm: number;
  readonly diagnostics?: IntegrationDiagnostics;
  /** Present only when `recordPath` was set. */
  readonly path?: readonly PhaseSpaceState[];
}

const STATE_COMPONENT_NAMES = [
  'x^0', 'x^1', 'x^2', 'x^3',
  'tangent^0', 'tangent^1', 'tangent^2', 'tangent^3',
] as const;

/**
 * Integrate a geodesic until a termination condition fires.
 *
 * NaN and infinity are detected explicitly after every accepted step (CLAUDE.md §17),
 * and the last finite state is returned together with diagnostics. The caller is told
 * what happened; nothing is quietly substituted.
 */
export function integrateGeodesic(options: IntegrateOptions): IntegrationResult {
  const { model, integrator, initial, terminator } = options;
  const limits: IntegrationLimits = { ...DEFAULT_LIMITS, ...options.limits };
  const derivative = options.derivative ?? geodesicDerivative(model);
  const paramName = parameterName(initial.kind);

  const y = packState(initial);
  const previous = new Float64Array(STATE_DIM);

  // Scratch tuples reused for every step's domain and terminator checks.
  //
  // A full image traces hundreds of thousands of rays of a few hundred steps each, and
  // allocating a fresh position and tangent per step dominated the render time. Both
  // callbacks receive readonly tuples and only read them, so one buffer each is enough.
  // They are transient: a callback that stores one will see it change under it.
  const positionScratch: [number, number, number, number] = [0, 0, 0, 0];
  const tangentScratch: [number, number, number, number] = [0, 0, 0, 0];

  const readPosition = (): Vec4 => {
    positionScratch[0] = y[0];
    positionScratch[1] = y[1];
    positionScratch[2] = y[2];
    positionScratch[3] = y[3];
    return positionScratch;
  };

  const readTangent = (): Vec4 => {
    tangentScratch[0] = y[4];
    tangentScratch[1] = y[5];
    tangentScratch[2] = y[6];
    tangentScratch[3] = y[7];
    return tangentScratch;
  };

  let parameter = initial.parameter;
  const parameterStart = parameter;
  let step = limits.initialStep;
  let steps = 0;
  let rejectedSteps = 0;
  let maxErrorNorm = Number.NaN;

  const path: PhaseSpaceState[] | undefined = options.recordPath ? [initial] : undefined;

  const finish = (
    reason: TerminationReason,
    diagnostics?: IntegrationDiagnostics,
  ): IntegrationResult => ({
    final: unpackState(y, initial.kind, parameter),
    reason,
    steps,
    rejectedSteps,
    finalStep: step,
    maxErrorNorm,
    ...(diagnostics ? { diagnostics } : {}),
    ...(path ? { path } : {}),
  });

  const diagnose = (failedQuantity: string, detail: string): IntegrationDiagnostics => ({
    integrator: integrator.id,
    model: model.id,
    coordinateSystem: model.chart.id,
    parameterName: paramName,
    parameterValue: parameter,
    stepSize: step,
    failedQuantity,
    detail,
  });

  // The starting event must itself be in the chart's domain.
  const initialDomain = model.domainCheck(initial.position_x);
  if (!initialDomain.inDomain) {
    return finish('domain-exit', diagnose('position_x', initialDomain.reason));
  }

  while (steps < limits.maxSteps) {
    const remaining = limits.parameterMax - Math.abs(parameter - parameterStart);
    if (remaining <= 0) return finish('parameter-limit');

    let attempt = Math.min(step, remaining);
    if (limits.maxStep !== undefined) attempt = Math.min(attempt, limits.maxStep);

    previous.set(y);
    const result = integrator.step(derivative, parameter, y, attempt);

    if (!result.accepted) {
      rejectedSteps += 1;
      // The integrator leaves `y` untouched on rejection, but restore defensively so a
      // partially written buffer from a future integrator cannot leak into the result.
      y.set(previous);
      step = result.nextStep;
      if (!(Math.abs(step) >= limits.minStep) || !Number.isFinite(step)) {
        return finish(
          'step-underflow',
          diagnose(
            'step size',
            `Adaptive control drove the step to ${step}, below the floor ${limits.minStep}, ` +
              `with a scaled local error norm of ${result.errorNorm}. The requested tolerance ` +
              'cannot be met here.',
          ),
        );
      }
      continue;
    }

    const bad = firstNonFiniteIndex(y);
    if (bad >= 0) {
      const component = STATE_COMPONENT_NAMES[bad];
      const badValue = y[bad];
      y.set(previous);
      return finish(
        'numerical-failure',
        diagnose(
          component,
          `Component ${component} became ${badValue} after an accepted step of ${attempt}. ` +
            'The last finite state is returned; no substitute value was invented.',
        ),
      );
    }

    steps += 1;
    parameter += attempt;
    step = result.nextStep;
    if (Number.isFinite(result.errorNorm)) {
      maxErrorNorm = Number.isNaN(maxErrorNorm)
        ? result.errorNorm
        : Math.max(maxErrorNorm, result.errorNorm);
    }

    if (path) path.push(unpackState(y, initial.kind, parameter));

    const domain = model.domainCheck(readPosition());
    if (!domain.inDomain) {
      return finish('domain-exit', diagnose('position_x', domain.reason));
    }

    if (terminator && terminator(readPosition(), readTangent(), parameter)) {
      return finish('terminator');
    }
  }

  return finish('step-limit');
}
