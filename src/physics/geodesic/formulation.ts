import type { Vec4 } from '../core/indices.js';
import type { PhaseSpaceState, WorldlineKind } from '../core/phase-space.js';
import type { SpacetimeModel } from '../spacetimes/spacetime-model.js';
import { geodesicDerivative, type DerivativeFn } from './geodesic-system.js';
import { STATE_DIM } from './state-vector.js';

/**
 * How the geodesic equation is written for the integrator (CLAUDE.md §2).
 *
 * Two mathematically equivalent formulations, kept side by side so that each can be
 * checked against the other (CLAUDE.md §7.3: "when two methods are available, benchmark
 * them"). They integrate the same curves; they differ in which quantities the
 * discretization preserves.
 *
 * Lagrangian (second-order, via Christoffel symbols), state (x^mu, k^mu):
 *
 *   dx^mu / dlambda = k^mu
 *   dk^mu / dlambda = -Gamma^mu_{alpha beta} k^alpha k^beta
 *
 * Hamiltonian, state (x^mu, p_mu) with H = 1/2 g^{mu nu} p_mu p_nu:
 *
 *   dx^mu / dlambda =  g^{mu nu} p_nu
 *   dp_mu / dlambda = -1/2 (d_mu g^{alpha beta}) p_alpha p_beta
 *
 * The two flows are identical — `tests/physics/hamiltonian.test.ts` records the symbolic
 * check — and the affine parameter is the same, since dx^mu/dlambda = g^{mu nu} p_nu is
 * exactly k^mu.
 *
 * The Hamiltonian form is preferred wherever it applies, for a structural reason. If the
 * metric is independent of a coordinate x^c, then d_c g^{alpha beta} is identically zero,
 * so dp_c/dlambda is identically zero, and every Runge-Kutta stage adds exactly nothing
 * to p_c. The conserved momentum is preserved to the last bit by *any* integrator. For
 * Schwarzschild that makes energy_E = -p_t and angular_momentum_Lz = p_phi exact, where
 * the Lagrangian form lets them drift at the integrator's tolerance. The only invariant
 * left to drift is the mass shell H itself, which becomes the honest diagnostic.
 */
export type FormulationId = 'lagrangian' | 'hamiltonian';

/**
 * A formulation bound to one model, owning the scratch buffers its hot path needs.
 *
 * Stateful and not shareable across concurrent integrations; build one per worker.
 */
export interface FormulationSession {
  readonly formulation: GeodesicFormulation;
  readonly model: SpacetimeModel;
  readonly derivative: DerivativeFn;
  /** Write the integrator state for a phase-space state. */
  pack(state: PhaseSpaceState, out: Float64Array): Float64Array;
  /** Read a phase-space state back out. The tangent is always contravariant. */
  unpack(y: Float64Array, kind: WorldlineKind, parameter: number): PhaseSpaceState;
  /** Write the contravariant tangent k^mu of a packed state into `out`. */
  tangentInto(y: Float64Array, out: [number, number, number, number]): void;
}

export interface GeodesicFormulation {
  readonly id: FormulationId;
  readonly displayName: string;
  readonly description: string;
  bind(model: SpacetimeModel): FormulationSession;
}

function requireStateBuffer(buffer: Float64Array, context: string): void {
  if (buffer.length !== STATE_DIM) {
    throw new RangeError(`${context}: expected a ${STATE_DIM}-entry state buffer.`);
  }
}

export const LAGRANGIAN: GeodesicFormulation = Object.freeze({
  id: 'lagrangian' as const,
  displayName: 'Geodesic equation (Christoffel form)',
  description:
    'Second-order geodesic equation split into first-order form, state (x^mu, k^mu). ' +
    'Conserved quantities are not built in and drift at the integrator tolerance.',

  bind(model: SpacetimeModel): FormulationSession {
    const derivative = geodesicDerivative(model);
    return {
      formulation: LAGRANGIAN,
      model,
      derivative,
      pack(state, out) {
        requireStateBuffer(out, 'LAGRANGIAN.pack');
        for (let mu = 0; mu < 4; mu += 1) {
          out[mu] = state.position_x[mu];
          out[4 + mu] = state.tangent[mu];
        }
        return out;
      },
      unpack(y, kind, parameter) {
        requireStateBuffer(y, 'LAGRANGIAN.unpack');
        return {
          kind,
          parameter,
          position_x: [y[0], y[1], y[2], y[3]],
          tangent: [y[4], y[5], y[6], y[7]],
        };
      },
      tangentInto(y, out) {
        out[0] = y[4];
        out[1] = y[5];
        out[2] = y[6];
        out[3] = y[7];
      },
    };
  },
});

/** Hamilton's equations for H = 1/2 g^{mu nu} p_mu p_nu, written in place. */
function hamiltonianDerivative(model: SpacetimeModel): DerivativeFn {
  const gInv = new Float64Array(16);
  const dgInv = new Float64Array(64);
  const position: [number, number, number, number] = [0, 0, 0, 0];

  return function hamiltonRhs(_p: number, y: Float64Array, dydp: Float64Array): void {
    position[0] = y[0];
    position[1] = y[1];
    position[2] = y[2];
    position[3] = y[3];

    model.inverseMetricInto(position, gInv);
    model.inverseMetricDerivativesInto(position, dgInv);

    // dx^mu / dlambda = g^{mu nu} p_nu
    for (let mu = 0; mu < 4; mu += 1) {
      let sum = 0;
      const row = mu * 4;
      for (let nu = 0; nu < 4; nu += 1) {
        const g = gInv[row + nu];
        if (g !== 0) sum += g * y[4 + nu];
      }
      dydp[mu] = sum;
    }

    // dp_mu / dlambda = -1/2 (d_mu g^{alpha beta}) p_alpha p_beta
    for (let mu = 0; mu < 4; mu += 1) {
      let sum = 0;
      const block = mu * 16;
      for (let alpha = 0; alpha < 4; alpha += 1) {
        const pAlpha = y[4 + alpha];
        if (pAlpha === 0) continue;
        const row = block + alpha * 4;
        for (let beta = 0; beta < 4; beta += 1) {
          const d = dgInv[row + beta];
          if (d !== 0) sum += d * pAlpha * y[4 + beta];
        }
      }
      dydp[4 + mu] = -0.5 * sum;
    }
  };
}

export const HAMILTONIAN: GeodesicFormulation = Object.freeze({
  id: 'hamiltonian' as const,
  displayName: "Hamilton's equations, H = 1/2 g^{mu nu} p_mu p_nu",
  description:
    'First-order Hamiltonian form, state (x^mu, p_mu). A momentum conjugate to a ' +
    'coordinate the metric does not depend on has an identically zero equation of ' +
    'motion, so it is conserved exactly by any integrator; only the mass shell ' +
    'H = 1/2 g^{mu nu} p_mu p_nu is left to drift.',

  bind(model: SpacetimeModel): FormulationSession {
    const derivative = hamiltonianDerivative(model);
    const g = new Float64Array(16);
    const gInv = new Float64Array(16);
    const position: [number, number, number, number] = [0, 0, 0, 0];

    const raise = (y: Float64Array, out: [number, number, number, number]): void => {
      position[0] = y[0];
      position[1] = y[1];
      position[2] = y[2];
      position[3] = y[3];
      model.inverseMetricInto(position, gInv);
      for (let mu = 0; mu < 4; mu += 1) {
        let sum = 0;
        for (let nu = 0; nu < 4; nu += 1) {
          const v = gInv[mu * 4 + nu];
          if (v !== 0) sum += v * y[4 + nu];
        }
        out[mu] = sum;
      }
    };

    return {
      formulation: HAMILTONIAN,
      model,
      derivative,
      pack(state, out) {
        requireStateBuffer(out, 'HAMILTONIAN.pack');
        model.metricInto(state.position_x, g);
        for (let mu = 0; mu < 4; mu += 1) {
          out[mu] = state.position_x[mu];
          let sum = 0;
          for (let nu = 0; nu < 4; nu += 1) {
            const v = g[mu * 4 + nu];
            if (v !== 0) sum += v * state.tangent[nu];
          }
          // p_mu = g_{mu nu} k^nu
          out[4 + mu] = sum;
        }
        return out;
      },
      unpack(y, kind, parameter) {
        requireStateBuffer(y, 'HAMILTONIAN.unpack');
        const tangent: [number, number, number, number] = [0, 0, 0, 0];
        raise(y, tangent);
        return {
          kind,
          parameter,
          position_x: [y[0], y[1], y[2], y[3]] as Vec4,
          tangent,
        };
      },
      tangentInto(y, out) {
        raise(y, out);
      },
    };
  },
});

export const FORMULATIONS: Readonly<Record<FormulationId, GeodesicFormulation>> = Object.freeze({
  lagrangian: LAGRANGIAN,
  hamiltonian: HAMILTONIAN,
});
