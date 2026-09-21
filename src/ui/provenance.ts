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
  /** Scene-specific entries, appended after the shared ones. */
  readonly sceneEntries?: readonly ProvenanceEntry[];
}): ProvenanceReport {
  const { model, observer, integrator, diagnostics, backgroundRadius } = options;
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
      value: 'Backward null-geodesic ray tracing, one ray per pixel',
      note:
        'For each pixel the observer-frame viewing direction is converted into a null ' +
        `wavevector k^mu and integrated until the ray reaches coordinate radius ` +
        `${backgroundRadius}, where the background grid is sampled. The photon follows a ` +
        'null geodesic of the chosen spacetime.',
    },
    {
      label: 'Omitted physics',
      value: 'Emission, absorption, plasma, radiative transfer, frequency shift',
      note:
        'There is no emitting matter and no relative motion between emitter and observer ' +
        'in this scene, so there is nothing to shift. Redshift, Doppler boosting, ' +
        'relativistic beaming and an accretion disk arrive with Milestone 3 rather than ' +
        'being approximated now. Nothing here is a prediction of an observed image.',
    },
    {
      label: 'Background',
      value: 'A latitude/longitude grid on a sphere at fixed coordinate radius',
      note:
        'A visualization mapping, not a physical emitting surface. It carries no emission ' +
        'model and no spectrum.',
    },
    {
      label: 'Image projection',
      value: 'Rectilinear pinhole projection of the celestial sphere',
      note:
        'Curvature visible in the grid lines is projection geometry, not light deflection. ' +
        'A rectilinear camera maps great circles to straight lines, so meridians (constant ' +
        'azimuth) appear straight, while parallels (constant polar angle) are small circles ' +
        'and appear curved -- the equator excepted, since it is a great circle. In flat ' +
        'spacetime this image is identical, pixel for pixel, to sampling the grid along ' +
        'each pixel\'s initial viewing direction with no integration at all, which is what ' +
        'the M1 validation suite asserts. Any deviation from that reference would be a ' +
        'defect, not lensing.',
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
    },
    {
      label: 'Null normalization',
      value: `max |g_mu_nu k^mu k^nu| = ${nullResidualText}`,
      note:
        'Target is exactly 0 for a null geodesic. Judged against the ' +
        `"${diagnostics.residualTolerance.id}" tolerance of ` +
        `${diagnostics.residualTolerance.value.toExponential(0)}. ` +
        diagnostics.residualTolerance.justification,
    },
    {
      label: 'Rays traced',
      value:
        `${diagnostics.raysTraced} total — ${diagnostics.raysReachingBackground} reached ` +
        `the background, ${diagnostics.raysCaptured} were captured, ` +
        `${diagnostics.raysFailed} failed numerically`,
    },
    {
      label: 'Integration steps',
      value: diagnostics.totalSteps.toLocaleString('en-US'),
    },
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
      'integrating the geodesic equations of a category A exact analytical metric. No ' +
      'observational data and no imported numerical-relativity data are used.',
  };
}
