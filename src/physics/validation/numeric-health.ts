import type { Vec4 } from '../core/indices.js';

/**
 * Numerical-validation layer (CLAUDE.md §20, §17): explicit NaN / infinity detection.
 *
 * CLAUDE.md §17 requires NaN and infinity detection to be explicit and forbids
 * silently replacing a failed physical calculation with a visually plausible result.
 * These helpers are the shared predicates the driver and the renderer both use.
 */

export function isFiniteVec4(v: Vec4): boolean {
  return (
    Number.isFinite(v[0]) && Number.isFinite(v[1]) && Number.isFinite(v[2]) && Number.isFinite(v[3])
  );
}

export function isFiniteBuffer(buffer: Float64Array): boolean {
  for (let i = 0; i < buffer.length; i += 1) {
    if (!Number.isFinite(buffer[i])) return false;
  }
  return true;
}

/** Index of the first non-finite entry, or -1 if all entries are finite. */
export function firstNonFiniteIndex(buffer: Float64Array): number {
  for (let i = 0; i < buffer.length; i += 1) {
    if (!Number.isFinite(buffer[i])) return i;
  }
  return -1;
}

export type HealthLevel = 'ok' | 'degraded' | 'failed';

export interface NumericalHealth {
  readonly level: HealthLevel;
  readonly messages: readonly string[];
}

/**
 * Combine tolerance outcomes into a single indicator for the UI.
 *
 * `degraded` means a checked quantity exceeded its validated tolerance; `failed` means
 * a non-finite value appeared. Both are reported, never smoothed over.
 */
export function summarizeHealth(
  checks: readonly { readonly withinTolerance: boolean; readonly message: string }[],
  sawNonFinite: boolean,
): NumericalHealth {
  const failures = checks.filter((c) => !c.withinTolerance).map((c) => c.message);
  if (sawNonFinite) {
    return {
      level: 'failed',
      messages: ['A NaN or infinity appeared during integration.', ...failures],
    };
  }
  if (failures.length > 0) return { level: 'degraded', messages: failures };
  return { level: 'ok', messages: [] };
}
