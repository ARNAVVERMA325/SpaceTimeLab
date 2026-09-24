import { buildProvenanceReport, type ProvenanceEntry } from './ui/provenance.js';
import { RenderPool, renderSerially } from './visualization/parallel-renderer.js';
import type { ImageAssembly } from './visualization/render-rows.js';
import { buildScene, type ObserverChoice, type SceneDescription } from './visualization/scene.js';

/**
 * Spacetime Lab application entry point.
 *
 * Three scenes share one pipeline: observer tetrad, ray generation, null-geodesic
 * integration, and an intersection test against whatever the scene contains. Changing
 * scene changes the metric, the observer and the emission model, and nothing else —
 * which is the point. Flat spacetime bends nothing and the sky grid arrives undistorted;
 * Schwarzschild bends light into a shadow and an Einstein ring; the thin disk adds an
 * emitting surface whose observed colour follows from the frequency shift of the traced
 * rays rather than from any artistic choice.
 *
 * Rendering happens in Web Workers, which is scheduling only: the image assembled from
 * row blocks is bit-identical to a single-threaded render, and the test suite checks that
 * on every scene kind. This is still the CPU reference path. CLAUDE.md §21 is explicit
 * that 60 FPS is a rendering goal rather than a scientific-validity requirement; the GPU
 * port is its own milestone with its own cross-validation gate.
 */

type SceneKind = SceneDescription['kind'];

const MAX_THREADS = 8;

function element<T extends HTMLElement>(id: string): T {
  const found = document.getElementById(id);
  if (!found) throw new Error(`main: expected an element with id "${id}".`);
  return found as T;
}

function threadCount(): number {
  if (typeof Worker === 'undefined') return 1;
  const cores = typeof navigator === 'object' ? navigator.hardwareConcurrency : undefined;
  return Math.max(1, Math.min(MAX_THREADS, cores ?? 4));
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
  const renderButton = element<HTMLButtonElement>('render');

  const sceneSelect = element<HTMLSelectElement>('scene');
  const widthSelect = element<HTMLSelectElement>('resolution');
  const samplesSelect = element<HTMLSelectElement>('samples');
  const observerSelect = element<HTMLSelectElement>('observer');
  const inclinationInput = element<HTMLInputElement>('inclination');
  const inclinationValue = element<HTMLElement>('inclination-value');
  const massSelect = element<HTMLSelectElement>('mass');
  const eddingtonSelect = element<HTMLSelectElement>('eddington');
  const exposureInput = element<HTMLInputElement>('exposure');
  const exposureValue = element<HTMLElement>('exposure-value');
  const skyGridInput = element<HTMLInputElement>('sky-grid');

  const threads = threadCount();
  const pool = typeof Worker === 'undefined' ? undefined : new RenderPool(threads);

  function describeScene(): SceneDescription {
    const kind = sceneSelect.value as SceneKind;
    const widthPx = Number.parseInt(widthSelect.value, 10);
    const heightPx = Math.round((widthPx * 3) / 4);
    const common = {
      widthPx,
      heightPx,
      samplesPerAxis: Number.parseInt(samplesSelect.value, 10),
      seed: 1,
    };
    const observer = observerSelect.value as ObserverChoice;

    if (kind === 'minkowski') return { ...common, kind };
    if (kind === 'schwarzschild-sky') {
      return { ...common, kind, cameraRadius: 20, observer };
    }
    return {
      ...common,
      kind: 'schwarzschild-disk',
      cameraRadius: 30,
      inclinationDeg: Number.parseFloat(inclinationInput.value),
      observer,
      massSolar: Number.parseFloat(massSelect.value),
      eddingtonFraction: Number.parseFloat(eddingtonSelect.value),
      exposureStops: Number.parseFloat(exposureInput.value),
      skyGrid: skyGridInput.checked,
    };
  }

  /** Controls that only mean something for some scenes are hidden for the others. */
  function syncControlVisibility(): void {
    const kind = sceneSelect.value as SceneKind;
    for (const field of document.querySelectorAll<HTMLElement>('[data-scenes]')) {
      const scenes = (field.dataset.scenes ?? '').split(' ');
      field.hidden = !scenes.includes(kind);
    }
    inclinationValue.textContent = `${Number.parseFloat(inclinationInput.value).toFixed(0)}°`;
    const stops = Number.parseFloat(exposureInput.value);
    exposureValue.textContent = `${stops > 0 ? '+' : ''}${stops.toFixed(1)}`;
  }

  let controller: AbortController | undefined;

  async function run(): Promise<void> {
    controller?.abort();
    pool?.cancel();
    const local = new AbortController();
    controller = local;

    const description = describeScene();
    const scene = buildScene(description);
    const context = canvas.getContext('2d');
    if (!context) throw new Error('main: could not acquire a 2D rendering context.');

    canvas.width = description.widthPx;
    canvas.height = description.heightPx;
    context.clearRect(0, 0, canvas.width, canvas.height);

    const rays = description.widthPx * description.heightPx * description.samplesPerAxis ** 2;
    renderButton.disabled = true;
    status.dataset.health = '';
    status.textContent = `Tracing ${rays.toLocaleString('en-US')} null geodesics on ${threads} thread${threads === 1 ? '' : 's'}...`;
    headline.textContent = scene.caption;

    const drawRows = (assembly: ImageAssembly, rows: readonly number[]): void => {
      const width = description.widthPx;
      for (const row of rows) {
        const line = assembly.pixels.subarray(row * width * 4, (row + 1) * width * 4);
        context.putImageData(new ImageData(line, width, 1), 0, row);
      }
      status.textContent =
        `Tracing null geodesics: ${assembly.rowsCompleted} / ${description.heightPx} rows ` +
        `(${threads} thread${threads === 1 ? '' : 's'})`;
    };

    const started = performance.now();
    try {
      const options = { onRows: drawRows, signal: local.signal };
      const result = pool
        ? await pool.render(description, options)
        : await renderSerially(description, options);
      const elapsedMs = performance.now() - started;

      const report = buildProvenanceReport({
        model: scene.model,
        observer: scene.observer,
        integrator: scene.config.integrator,
        formulation: scene.config.formulation,
        diagnostics: result.diagnostics,
        backgroundRadius: scene.config.grid.radius,
        sceneEntries: scene.entries,
        emitter: description.kind === 'schwarzschild-disk' ? 'thin-disk' : 'none',
        execution: { threads: pool ? threads : 1, elapsedMs },
      });

      renderEntries(modelEntries, report.entries);
      renderEntries(validationEntries, report.validation);
      hierarchy.textContent = report.dataHierarchy;

      const { health } = result.diagnostics;
      status.textContent =
        `${result.diagnostics.raysTraced.toLocaleString('en-US')} rays, ` +
        `${result.diagnostics.totalSteps.toLocaleString('en-US')} integration steps, ` +
        `${(elapsedMs / 1000).toFixed(1)} s on ${pool ? threads : 1} thread${pool && threads > 1 ? 's' : ''}. ` +
        (health.level === 'ok'
          ? 'All checked invariants within validated tolerance.'
          : health.messages.join(' '));
      status.dataset.health = health.level;
    } catch (error: unknown) {
      if (local.signal.aborted) return;
      throw error;
    } finally {
      if (controller === local) renderButton.disabled = false;
    }
  }

  sceneSelect.addEventListener('change', () => {
    syncControlVisibility();
    void run();
  });
  for (const control of [widthSelect, samplesSelect, observerSelect, massSelect, eddingtonSelect]) {
    control.addEventListener('change', () => void run());
  }
  skyGridInput.addEventListener('change', () => void run());
  for (const slider of [inclinationInput, exposureInput]) {
    slider.addEventListener('input', syncControlVisibility);
    slider.addEventListener('change', () => void run());
  }
  renderButton.addEventListener('click', () => void run());

  syncControlVisibility();
  await run();
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
