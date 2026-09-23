import type { Vec4 } from '../core/indices.js';
import type { PhaseSpaceState } from '../core/phase-space.js';
import { parameterName } from '../core/phase-space.js';
import type { SpacetimeModel } from '../spacetimes/spacetime-model.js';
import { firstNonFiniteIndex } from '../validation/numeric-health.js';
import { HAMILTONIAN, type FormulationSession } from './formulation.js';
import type { Integrator } from './integrators/integrator.js';
import { STATE_DIM } from './state-vector.js';

/**
 * Geodesic / dynamical layer driver (CLAUDE.md §20).
 *
 * The driver owns termination policy, event location, step bookkeeping and
 * numerical-health checks. The model owns the geometry, the formulation owns how the
 * geodesic equation is written, and the integrator owns the stepping; none of them
 * decides when to stop.
 */

export type TerminationReason =
  /** Reached the requested parameter bound. */
  | 'parameter-limit'
  /** Hit the step budget before any other condition. */
  | 'step-limit'
  /** The model reported the event is outside its chart's valid domain. */
  | 'domain-exit'
  /** A caller-supplied boolean terminator fired after a step. */
  | 'terminator'
  /** A terminal event was located, and the state was placed exactly on it. */
  | 'event'
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
  readonly formulation: string;
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
 * A boolean stopping condition evaluated after each accepted step.
 *
 * Suitable for conditions that are *regions* rather than surfaces — "inside the photon
 * sphere and moving inward" — where there is no crossing to locate. For a surface, use
 * an event instead, which places the state on it exactly.
 */
export type Terminator = (position_x: Vec4, tangent: Vec4, parameter: number) => boolean;

/**
 * A surface in phase space, defined as the zero set of a scalar function.
 *
 * After every accepted step the driver checks each event function for a sign change,
 * and locates any crossing to near machine precision in the step parameter using the
 * integrator's continuous extension and Brent's method. This is the same approach as
 * SciPy's solve_ivp events. The reported state then lies on the surface to the accuracy
 * of the continuous extension, instead of wherever the step happened to land.
 */
export interface GeodesicEvent {
  readonly id: string;
  value(position_x: Vec4, tangent: Vec4): number;
  /** +1: only crossings from negative to positive; -1: the reverse; 0: either. */
  readonly direction: -1 | 0 | 1;
  /**
   * Whether reaching the event stops the integration. A function allows the decision
   * to depend on where the crossing happened — a ray crossing the equatorial plane
   * stops if it hits the disk annulus and passes straight through elsewhere.
   */
  readonly terminal: boolean | ((state: PhaseSpaceState) => boolean);
}

export interface EventHit {
  readonly id: string;
  readonly parameter: number;
  readonly state: PhaseSpaceState;
  readonly terminal: boolean;
}

export interface IntegrateOptions {
  readonly model: SpacetimeModel;
  readonly integrator: Integrator;
  readonly initial: PhaseSpaceState;
  readonly limits?: Partial<IntegrationLimits>;
  readonly terminator?: Terminator;
  readonly events?: readonly GeodesicEvent[];
  /**
   * A formulation already bound to `model`, reused across many integrations.
   * Defaults to a fresh Hamiltonian session.
   */
  readonly session?: FormulationSession;
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
  /** Every event crossing located, in order, including the terminal one if any. */
  readonly events: readonly EventHit[];
  readonly formulation: string;
  /**
   * The raw final integrator state, [x^0..x^3, momentum_0..momentum_3]. For the
   * Hamiltonian formulation the momentum half is the covariant p_mu that was actually
   * integrated, which is where exact conservation of p_t and p_phi can be observed
   * directly — the contravariant tangent in `final` has been through a metric
   * contraction and carries its rounding.
   */
  readonly packed: readonly number[];
  readonly diagnostics?: IntegrationDiagnostics;
  /** Present only when `recordPath` was set. */
  readonly path?: readonly PhaseSpaceState[];
}

const STATE_COMPONENT_NAMES = [
  'x^0', 'x^1', 'x^2', 'x^3',
  'momentum^0', 'momentum^1', 'momentum^2', 'momentum^3',
] as const;

/**
 * Brent's method for a root of g on [a, b], given g(a) and g(b) of opposite sign.
 *
 * Brent, "Algorithms for Minimization without Derivatives" (1973), ch. 4: inverse
 * quadratic interpolation with bisection safeguarding, so it converges superlinearly on
 * smooth functions and never worse than bisection.
 */
export function brentRoot(
  g: (x: number) => number,
  a: number,
  b: number,
  ga: number,
  gb: number,
  tolerance = 4 * Number.EPSILON,
  maxIterations = 100,
): number {
  if (ga === 0) return a;
  if (gb === 0) return b;
  if (ga * gb > 0) {
    throw new RangeError('brentRoot: the bracket does not straddle a sign change.');
  }

  let c = a;
  let gc = ga;
  let d = b - a;
  let e = d;

  for (let iteration = 0; iteration < maxIterations; iteration += 1) {
    if (gb * gc > 0) {
      c = a;
      gc = ga;
      d = b - a;
      e = d;
    }
    if (Math.abs(gc) < Math.abs(gb)) {
      a = b;
      b = c;
      c = a;
      ga = gb;
      gb = gc;
      gc = ga;
    }

    const tol = 2 * tolerance * Math.abs(b) + 0.5 * tolerance;
    const m = 0.5 * (c - b);
    if (Math.abs(m) <= tol || gb === 0) return b;

    if (Math.abs(e) >= tol && Math.abs(ga) > Math.abs(gb)) {
      const s = gb / ga;
      let pNum: number;
      let qDen: number;
      if (a === c) {
        pNum = 2 * m * s;
        qDen = 1 - s;
      } else {
        const q = ga / gc;
        const r = gb / gc;
        pNum = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
        qDen = (q - 1) * (r - 1) * (s - 1);
      }
      if (pNum > 0) qDen = -qDen;
      else pNum = -pNum;
      if (2 * pNum < Math.min(3 * m * qDen - Math.abs(tol * qDen), Math.abs(e * qDen))) {
        e = d;
        d = pNum / qDen;
      } else {
        d = m;
        e = m;
      }
    } else {
      d = m;
      e = m;
    }

    a = b;
    ga = gb;
    b += Math.abs(d) > tol ? d : m > 0 ? tol : -tol;
    gb = g(b);
  }
  return b;
}

function crosses(before: number, after: number, direction: -1 | 0 | 1): boolean {
  if (!Number.isFinite(before) || !Number.isFinite(after)) return false;
  if (before === 0) return false;
  const signChange = before < 0 ? after >= 0 : after <= 0;
  if (!signChange) return false;
  if (direction === 0) return true;
  return direction > 0 ? before < 0 : before > 0;
}

/**
 * Integrate a geodesic until a termination condition fires.
 *
 * NaN and infinity are detected explicitly after every accepted step (CLAUDE.md §17),
 * and the last finite state is returned together with diagnostics. The caller is told
 * what happened; nothing is quietly substituted.
 */
export function integrateGeodesic(options: IntegrateOptions): IntegrationResult {
  const { model, integrator, initial, terminator } = options;
  const events = options.events ?? [];
  const limits: IntegrationLimits = { ...DEFAULT_LIMITS, ...options.limits };
  const session = options.session ?? HAMILTONIAN.bind(model);
  if (session.model !== model) {
    throw new RangeError('integrateGeodesic: the formulation session is bound to a different model.');
  }
  const derivative = session.derivative;
  const paramName = parameterName(initial.kind);

  const y = session.pack(initial, new Float64Array(STATE_DIM));
  const previous = new Float64Array(STATE_DIM);
  const probe = new Float64Array(STATE_DIM);

  let parameter = initial.parameter;
  const parameterStart = parameter;
  let step = limits.initialStep;
  let steps = 0;
  let rejectedSteps = 0;
  let maxErrorNorm = Number.NaN;
  const hits: EventHit[] = [];

  const path: PhaseSpaceState[] | undefined = options.recordPath ? [initial] : undefined;

  // Scratch reused for every callback. Callbacks receive readonly tuples and only read
  // them; one that stores a reference will see it change under it.
  const positionScratch: [number, number, number, number] = [0, 0, 0, 0];
  const tangentScratch: [number, number, number, number] = [0, 0, 0, 0];

  const readPosition = (buffer: Float64Array): Vec4 => {
    positionScratch[0] = buffer[0];
    positionScratch[1] = buffer[1];
    positionScratch[2] = buffer[2];
    positionScratch[3] = buffer[3];
    return positionScratch;
  };
  const readTangent = (buffer: Float64Array): Vec4 => {
    session.tangentInto(buffer, tangentScratch);
    return tangentScratch;
  };
  const eventValue = (event: GeodesicEvent, buffer: Float64Array): number =>
    event.value(readPosition(buffer), readTangent(buffer));

  const finish = (reason: TerminationReason, diagnostics?: IntegrationDiagnostics): IntegrationResult => ({
    final: session.unpack(y, initial.kind, parameter),
    reason,
    steps,
    rejectedSteps,
    finalStep: step,
    maxErrorNorm,
    events: hits,
    formulation: session.formulation.id,
    packed: Array.from(y),
    ...(diagnostics ? { diagnostics } : {}),
    ...(path ? { path } : {}),
  });

  const diagnose = (failedQuantity: string, detail: string): IntegrationDiagnostics => ({
    integrator: integrator.id,
    formulation: session.formulation.id,
    model: model.id,
    coordinateSystem: model.chart.id,
    parameterName: paramName,
    parameterValue: parameter,
    stepSize: step,
    failedQuantity,
    detail,
  });

  const initialDomain = model.domainCheck(initial.position_x);
  if (!initialDomain.inDomain) {
    return finish('domain-exit', diagnose('position_x', initialDomain.reason));
  }

  let eventValuesBefore = events.map((event) => eventValue(event, y));

  // Continuous extension over the last accepted step: the integrator's own if it has
  // one, otherwise cubic Hermite from the derivatives at both ends (third order, which
  // is enough to bracket and refine a crossing for the fixed-step methods).
  const hermiteStart = new Float64Array(STATE_DIM);
  const hermiteEnd = new Float64Array(STATE_DIM);
  const makeInterpolant = (
    yStart: Float64Array,
    yEnd: Float64Array,
    pStart: number,
    h: number,
  ): ((theta: number, out: Float64Array) => void) => {
    if (integrator.interpolate) {
      return (theta, out) => integrator.interpolate!(theta, out);
    }
    derivative(pStart, yStart, hermiteStart);
    derivative(pStart + h, yEnd, hermiteEnd);
    return (theta, out) => {
      const t2 = theta * theta;
      const t3 = t2 * theta;
      const h00 = 2 * t3 - 3 * t2 + 1;
      const h10 = t3 - 2 * t2 + theta;
      const h01 = -2 * t3 + 3 * t2;
      const h11 = t3 - t2;
      for (let i = 0; i < STATE_DIM; i += 1) {
        out[i] =
          h00 * yStart[i] + h10 * h * hermiteStart[i] + h01 * yEnd[i] + h11 * h * hermiteEnd[i];
      }
    };
  };

  while (steps < limits.maxSteps) {
    const remaining = limits.parameterMax - Math.abs(parameter - parameterStart);
    if (remaining <= 0) return finish('parameter-limit');

    let attempt = Math.min(step, remaining);
    if (limits.maxStep !== undefined) attempt = Math.min(attempt, limits.maxStep);

    previous.set(y);
    const result = integrator.step(derivative, parameter, y, attempt);

    if (!result.accepted) {
      rejectedSteps += 1;
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

    const stepStart = parameter;
    steps += 1;
    parameter += attempt;
    step = result.nextStep;
    if (Number.isFinite(result.errorNorm)) {
      maxErrorNorm = Number.isNaN(maxErrorNorm)
        ? result.errorNorm
        : Math.max(maxErrorNorm, result.errorNorm);
    }

    // Event location over the step just taken.
    if (events.length > 0) {
      const eventValuesAfter = events.map((event) => eventValue(event, y));
      const crossings: { index: number; theta: number }[] = [];
      let interpolate: ((theta: number, out: Float64Array) => void) | undefined;

      for (let i = 0; i < events.length; i += 1) {
        if (!crosses(eventValuesBefore[i], eventValuesAfter[i], events[i].direction)) continue;
        interpolate ??= makeInterpolant(previous, y, stepStart, attempt);
        const g = (theta: number): number => {
          interpolate!(theta, probe);
          return eventValue(events[i], probe);
        };
        const theta = brentRoot(g, 0, 1, eventValuesBefore[i], eventValuesAfter[i]);
        crossings.push({ index: i, theta });
      }

      crossings.sort((a, b) => a.theta - b.theta);
      for (const crossing of crossings) {
        interpolate!(crossing.theta, probe);
        const event = events[crossing.index];
        const hitParameter = stepStart + crossing.theta * attempt;
        const interpolated = session.unpack(probe, initial.kind, hitParameter);
        const terminal =
          typeof event.terminal === 'function' ? event.terminal(interpolated) : event.terminal;

        if (!terminal) {
          // A non-terminal crossing is recorded from the continuous extension; the
          // integration itself carries on from the step's true endpoint, untouched.
          hits.push({ id: event.id, parameter: hitParameter, state: interpolated, terminal });
          continue;
        }

        // A terminal crossing becomes the final state, so it is worth landing on with a
        // genuine integrator step rather than trusting the interpolant. The continuous
        // extension is used only to find *where* the event is; the state reported there
        // then carries the integrator's full order. Without this, an integrator lacking
        // its own dense output falls back to third-order Hermite interpolation, which is
        // less accurate than the method itself and visibly degraded the Milestone 2A
        // deflection benchmark before this refinement was added.
        //
        // The root found on the interpolant and the endpoint of a real step to it differ
        // by the interpolant's own error, so the landed state can sit a tolerance-sized
        // distance off the surface. A short Newton iteration on theta, driven by the event
        // value at genuinely integrated states and the slope from the interpolant, removes
        // that residual. It runs only at terminal events — once per ray — so its cost is
        // negligible.
        const landed = new Float64Array(STATE_DIM);
        const g = (buffer: Float64Array): number => eventValue(event, buffer);
        const eventScale =
          Math.abs(eventValuesAfter[crossing.index] - eventValuesBefore[crossing.index]) || 1;

        // Everything needed from the interpolant is taken now, before any re-step: an
        // accepted re-step replaces the integrator's continuous extension with its own,
        // so the original step's interpolant is gone afterwards.
        const theta0 = crossing.theta;
        const fallback = new Float64Array(STATE_DIM);
        interpolate!(theta0, fallback);
        const delta = 1e-6;
        const thetaHi = Math.min(1, theta0 + delta);
        const thetaLo = Math.max(0, theta0 - delta);
        interpolate!(thetaHi, probe);
        const gHi = g(probe);
        interpolate!(thetaLo, probe);
        const gLo = g(probe);
        const slope = (gHi - gLo) / (thetaHi - thetaLo);

        const stepTo = (theta: number): boolean => {
          landed.set(previous);
          if (theta <= 0) return true;
          const r = integrator.step(derivative, stepStart, landed, theta * attempt);
          return r.accepted && firstNonFiniteIndex(landed) < 0;
        };

        let theta = theta0;
        let landedOk = stepTo(theta);
        if (Number.isFinite(slope) && slope !== 0) {
          for (let iteration = 0; landedOk && iteration < 4; iteration += 1) {
            const residual = g(landed);
            if (Math.abs(residual) <= 1e-14 * eventScale) break;
            const next = Math.min(1, Math.max(0, theta - residual / slope));
            if (next === theta) break;
            theta = next;
            landedOk = stepTo(theta);
          }
        }

        if (landedOk) {
          y.set(landed);
        } else {
          // The shorter step should always be accepted, since its local error is smaller.
          // If it somehow is not, the interpolated state is still a valid answer.
          y.set(fallback);
          theta = theta0;
        }
        parameter = stepStart + theta * attempt;
        const state = session.unpack(y, initial.kind, parameter);
        hits.push({ id: event.id, parameter, state, terminal });
        if (path) path.push(state);
        return finish('event');
      }
      eventValuesBefore = eventValuesAfter;
    }

    if (path) path.push(session.unpack(y, initial.kind, parameter));

    const domain = model.domainCheck(readPosition(y));
    if (!domain.inDomain) {
      return finish('domain-exit', diagnose('position_x', domain.reason));
    }

    if (terminator && terminator(readPosition(y), readTangent(y), parameter)) {
      return finish('terminator');
    }
  }

  return finish('step-limit');
}
