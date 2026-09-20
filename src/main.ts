import { RKF45Integrator } from './physics/geodesic/integrators/rkf45.js';
import { STATE_DIM } from './physics/geodesic/state-vector.js';
import { defaultScreen, staticMinkowskiObserver } from './physics/observer/observer.js';
import { minkowski } from './physics/spacetimes/minkowski.js';
import { buildProvenanceReport, type ProvenanceEntry } from './ui/provenance.js';
import { renderToCanvas } from './visualization/canvas-renderer.js';
import { DEFAULT_CELESTIAL_GRID } from './visualization/celestial-grid.js';
import type { TraceConfig } from './visualization/raytracer.js';

/**
 * Milestone 1 application entry point.
 *
 * ROADMAP.md 1.5 calls this visualizer deliberately undramatic: flat spacetime produces
 * an undistorted grid, and the point is proving the pipeline rather than the picture.
 * Every stage the later milestones need is already here and in the right order —
 * observer frame, ray generation, null-geodesic integration, background intersection —
 * so Milestone 2A changes the metric, not the architecture.
 */

const OBSERVER_POSITION = [0, 0, 0, 0] as const;
const HORIZONTAL_FOV = Math.PI / 2;

function element<T extends HTMLElement>(id: string): T {
  const found = document.getElementById(id);
  if (!found) throw new Error(`main: expected an element with id "${id}".`);
  return found as T;
}

function renderEntries(container: HTMLElement, entries: readonly ProvenanceEntry[]): void {
  container.replaceChildren(
    ...entries.map((entry) => {
      const row = document.createElement('div');
      row.className = 'entry';

      const label = document.createElement('div');
      label.className = 'entry-label';
      label.textContent = entry.label;

      const value = document.createElement('div');
      value.className = 'entry-value';
      value.textContent = entry.value;

      row.append(label, value);

      if (entry.note) {
        const note = document.createElement('div');
        note.className = 'entry-note';
        note.textContent = entry.note;
        row.append(note);
      }
      return row;
    }),
  );
}

async function main(): Promise<void> {
  const canvas = element<HTMLCanvasElement>('view');
  const status = element<HTMLElement>('status');
  const headline = element<HTMLElement>('headline');
  const modelEntries = element<HTMLElement>('model-entries');
  const validationEntries = element<HTMLElement>('validation-entries');
  const hierarchy = element<HTMLElement>('data-hierarchy');

  const observer = staticMinkowskiObserver([...OBSERVER_POSITION]);
  const grid = DEFAULT_CELESTIAL_GRID;

  // The adaptive method is the right choice here and later: in flat space the embedded
  // pair agrees exactly and the controller strides across the empty region, while in a
  // curved spacetime the same controller will shorten its step where curvature demands
  // it. The step cap keeps a ray from overshooting the background sphere by a wide
  // margin, which will matter once the geometry actually bends the ray.
  const integrator = new RKF45Integrator(STATE_DIM);

  const config: TraceConfig = {
    model: minkowski,
    integrator,
    observer,
    grid,
    limits: {
      initialStep: 0.5,
      maxStep: grid.radius / 10,
      parameterMax: grid.radius * 10,
      maxSteps: 10_000,
    },
  };

  const screen = defaultScreen(
    Math.min(720, Math.max(240, Math.floor(canvas.clientWidth || 640))),
    Math.min(540, Math.max(180, Math.floor((canvas.clientWidth || 640) * 0.75))),
    HORIZONTAL_FOV,
  );

  status.textContent = `Tracing ${screen.widthPx * screen.heightPx} null geodesics...`;

  const started = performance.now();
  const result = await renderToCanvas(canvas, {
    config,
    screen,
    rowsPerBand: 12,
    onProgress: (rows, total) => {
      status.textContent = `Tracing null geodesics: ${rows} / ${total} rows`;
    },
  });
  const elapsedMs = performance.now() - started;

  const report = buildProvenanceReport({
    model: minkowski,
    observer,
    integrator,
    diagnostics: result.diagnostics,
    backgroundRadius: grid.radius,
  });

  headline.textContent = report.headline;
  renderEntries(modelEntries, report.entries);
  renderEntries(validationEntries, report.validation);
  hierarchy.textContent = report.dataHierarchy;

  const health = result.diagnostics.health;
  status.textContent =
    `${result.diagnostics.raysTraced} rays, ` +
    `${result.diagnostics.totalSteps.toLocaleString('en-US')} integration steps, ` +
    `${elapsedMs.toFixed(0)} ms. ` +
    (health.level === 'ok'
      ? 'All checked invariants within validated tolerance.'
      : health.messages.join(' '));
  status.dataset.health = health.level;
}

main().catch((error: unknown) => {
  // A failed physical calculation is reported, never replaced by a plausible-looking
  // image (CLAUDE.md §17, §24).
  const status = document.getElementById('status');
  const message = error instanceof Error ? error.message : String(error);
  if (status) {
    status.textContent = `Render failed: ${message}`;
    status.dataset.health = 'failed';
  }
  console.error(error);
});
