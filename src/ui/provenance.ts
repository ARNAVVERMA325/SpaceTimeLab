import { HAMILTONIAN, type GeodesicFormulation } from '../physics/geodesic/formulation.js';
import type { Integrator } from '../physics/geodesic/integrators/integrator.js';
import type { Observer } from '../physics/observer/observer.js';
import type { SpacetimeModel } from '../physics/spacetimes/spacetime-model.js';
import type { RenderDiagnostics } from '../visualization/raytracer.js';
import { ALL_TOLERANCES } from '../physics/validation/tolerances.js';

/**
 * UI / educational layer (CLAUDE.md §20, §22).
 *
 * CLAUDE.md §22 requires every major visualization to let the user inspect the model,
 * metric, coordinates, units, observer, numerical method, approximation level,
 * external-data provenance and validation status. This module assembles exactly that
 * list from the objects actually used in the render, so the panel cannot drift out of
 * step with what was computed.
 *
 * The wording follows §22's guidance: no "gravity pulls the light", no "this is what a
 * black hole really looks like".
 */

export interface ProvenanceEntry {
  readonly label: string;
  readonly value: string;
  readonly note?: string;
  /** Semantic state, so a verdict reads as one at a glance rather than as more prose. */
  readonly tone?: 'ok' | 'warn' | 'fail';
}

export interface ProvenanceReport {
  readonly headline: string;
  readonly entries: readonly ProvenanceEntry[];
  readonly validation: readonly ProvenanceEntry[];
  readonly dataHierarchy: string;
}

const CLASSIFICATION_LABEL: Record<SpacetimeModel['classification'], string> = {
  'exact-analytical': 'A. Exact analytical — a mathematically defined exact solution',
  'approximate-perturbative':
    'B. Approximate / perturbative — not an exact solution of the full nonlinear equations',
  'theoretical-metric-model':
    'A specified theoretical metric model / ansatz. No observational confirmation implied.',
  'imported-numerical-relativity':
    'C. Imported numerical-relativity result — precomputed, not a live browser calculation',
};

export function buildProvenanceReport(options: {
  readonly model: SpacetimeModel;
  readonly observer: Observer;
  readonly integrator: Integrator;
  readonly diagnostics: RenderDiagnostics;
  readonly backgroundRadius: number;
  /** The formulation integrated; defaults to the renderer's default, Hamiltonian. */
  readonly formulation?: GeodesicFormulation;
  /** Scene-specific entries, appended after the shared ones. */
  readonly sceneEntries?: readonly ProvenanceEntry[];
  /** Whether rays can end on an emitting disk as well as on the background. */
  readonly emitter?: 'none' | 'thin-disk';
  /** How the render was executed; scheduling only, recorded so timings can be compared. */
  readonly execution?: { readonly threads: number; readonly elapsedMs: number };
}): ProvenanceReport {
  const { model, observer, integrator, diagnostics, backgroundRadius } = options;
  const emitter = options.emitter ?? 'none';
  const samples = diagnostics.samplesPerPixel;
  const coords = model.chart.coordinateNames.join(', ');

  const entries: ProvenanceEntry[] = [
    {
      label: 'Spacetime model',
      value: model.displayName,
      note: model.description,
    },
    {
      label: 'Model classification',
      value: CLASSIFICATION_LABEL[model.classification],
    },
    {
      label: 'Metric',
      value: `g_mu_nu in (${coords})`,
      note:
        'Geodesics are integrated through this specified geometry. That is evaluating the ' +
        'consequences of a known spacetime, not solving the Einstein field equations ' +
        'dynamically.',
    },
    {
      label: 'Coordinates',
      value: model.chart.displayName,
      note: model.chart.notes,
    },
    {
      label: 'Signature',
      value: model.conventions.signature,
    },
    {
      label: 'Units',
      value: model.conventions.unitSystem,
      note: model.conventions.unitSystemNote,
    },
    {
      label: 'Observer',
      value: `${observer.displayName} (${
        observer.kind === 'physical-observer'
          ? 'physical observer'
          : 'visualization coordinate system, not a physical observer'
      })`,
      note: observer.description,
    },
    {
      label: 'Geodesic formulation',
      value: (options.formulation ?? HAMILTONIAN).displayName,
      note: (options.formulation ?? HAMILTONIAN).description,
    },
    {
      label: 'Numerical method',
      value: integrator.displayName,
      note: integrator.description,
    },
    {
      label: 'Floating-point precision',
      value: model.conventions.floatingPoint,
    },
    {
      label: 'Rendering method',
      value: `Backward null-geodesic ray tracing, ${samples} ray${samples === 1 ? '' : 's'} per pixel`,
      note:
        'For each sample the observer-frame viewing direction is converted into a null ' +
        'wavevector k^mu and integrated backward in time until the ray ' +
        (emitter === 'thin-disk' ? 'crosses the disk, ' : '') +
        `is captured, or reaches coordinate radius ${backgroundRadius}, where the ` +
        'background is sampled. The photon follows a null geodesic of the chosen spacetime.',
    },
  ];

  const nullResidualText = Number.isFinite(diagnostics.maxNullResidual)
    ? diagnostics.maxNullResidual.toExponential(3)
    : String(diagnostics.maxNullResidual);

  const validation: ProvenanceEntry[] = [
    {
      label: 'Validation status',
      value:
        diagnostics.health.level === 'ok'
          ? 'All checked invariants within validated tolerance'
          : diagnostics.health.messages.join(' '),
      tone: diagnostics.health.level === 'ok' ? 'ok' : diagnostics.health.level === 'degraded' ? 'warn' : 'fail',
    },
    {
      label: 'Null normalization',
      value: `max |g_mu_nu k^mu k^nu| = ${nullResidualText}`,
      tone: diagnostics.maxNullResidual <= diagnostics.residualTolerance.value ? 'ok' : 'fail',
      note:
        'Target is exactly 0 for a null geodesic. Judged against the ' +
        `"${diagnostics.residualTolerance.id}" tolerance of ` +
        `${diagnostics.residualTolerance.value.toExponential(0)}. ` +
        diagnostics.residualTolerance.justification,
    },
    {
      label: 'Rays traced',
      value:
        `${diagnostics.raysTraced} total at ${diagnostics.samplesPerPixel} per pixel — ` +
        `${diagnostics.raysReachingBackground} reached ` +
        `the background, ${diagnostics.raysCaptured} were captured, ` +
        `${diagnostics.raysFailed} failed numerically`,
    },
    {
      label: 'Integration steps',
      value: diagnostics.totalSteps.toLocaleString('en-US'),
    },
    ...diskValidation(diagnostics, emitter),
    ...(diagnostics.pixelsWithFailures > 0
      ? [
          {
            label: 'Failed pixels',
            value: `${diagnostics.pixelsWithFailures} pixel(s) drawn in magenta`,
            tone: 'fail' as const,
            note:
              'A pixel containing any numerically failed sample is painted in the failure ' +
              'colour instead of being averaged into a plausible tint (CLAUDE.md §17).',
          },
        ]
      : []),
    ...(options.execution
      ? [
          {
            label: 'Execution',
            value:
              `${(options.execution.elapsedMs / 1000).toFixed(1)} s on ` +
              `${options.execution.threads} thread${options.execution.threads === 1 ? '' : 's'}`,
            note:
              'Rows are dealt to parallel workers, each tracing its rays independently. The ' +
              'assembled image is bit-identical to a single-threaded render, which the test ' +
              'suite checks; the thread count changes the time taken, never the result.',
          },
        ]
      : []),
    ...ALL_TOLERANCES.map((t) => ({
      label: `Tolerance: ${t.id}`,
      value: `${t.kind} ${t.value.toExponential(0)} on ${t.quantity}`,
      note: t.justification,
    })),
  ];

  if (options.sceneEntries) entries.push(...options.sceneEntries);

  return {
    headline:
      'Computed appearance for the selected spacetime, observer, background model and ' +
      'rendering assumptions.',
    entries,
    validation,
    dataHierarchy:
      'Everything shown here is category B: a numerical approximation obtained by ' +
      'integrating the geodesic equations of a category A exact analytical metric' +
      (emitter === 'thin-disk'
        ? ', lit by an analytical emission model (Novikov-Thorne) whose physical scale uses ' +
          'CODATA 2018 constants. The CIE 1931 colour-matching functions are a standard ' +
          'colorimetric table, used only to turn a computed spectrum into display colour.'
        : '.') +
      ' No observational data and no imported numerical-relativity data are used.',
  };
}

function formatRange([low, high]: readonly [number, number], digits: number): string {
  return `${low.toFixed(digits)} to ${high.toFixed(digits)}`;
}

/**
 * What a disk render measured about its own light: the range of frequency shift and the
 * range of observed temperature over every sample that landed on the disk.
 */
function diskValidation(diagnostics: RenderDiagnostics, emitter: 'none' | 'thin-disk'): ProvenanceEntry[] {
  if (emitter !== 'thin-disk') return [];
  if (diagnostics.raysHittingDisk === 0) {
    return [{ label: 'Disk samples', value: 'No sample landed on the disk in this view' }];
  }
  const [tLow, tHigh] = diagnostics.observedTemperatureRange;
  return [
    {
      label: 'Disk samples',
      value: `${diagnostics.raysHittingDisk.toLocaleString('en-US')} rays ended on the disk`,
    },
    {
      label: 'Frequency shift over the disk',
      value: `g = nu_obs / nu_emit from ${formatRange(diagnostics.frequencyRatioRange, 3)}`,
      note:
        'Measured on this render, from the traced wavevector and the observer and emitter ' +
        'four-velocities. g < 1 is a net redshift, g > 1 a net blueshift. The same code is ' +
        'tested against the static-to-static redshift sqrt(f_emit / f_obs) (to 1e-12), the ' +
        'closed form sqrt(1 - 3M/r) / sqrt(f_obs) for disk rays with no angular momentum about ' +
        'the axis (1e-10), and an independent covariant evaluation on traced rays (1e-12).',
    },
    {
      label: 'Observed temperature',
      value:
        `${Math.round(tLow).toLocaleString('en-US')} K to ` +
        `${Math.round(tHigh).toLocaleString('en-US')} K`,
      note:
        'The colour temperature of the observed blackbody, g T(r), over every disk sample. ' +
        'A blackbody at T seen with frequency ratio g is exactly a blackbody at g T.',
    },
  ];
}
