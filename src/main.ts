import { buildProvenanceReport, type ProvenanceEntry } from './ui/provenance.js';
import { inspectPixel, type InspectorRow, type Tone } from './ui/ray-inspector.js';
import { RenderPool, renderSerially } from './visualization/parallel-renderer.js';
import type { ImageAssembly } from './visualization/render-rows.js';
import { buildScene, type BuiltScene, type ObserverChoice, type SceneDescription } from './visualization/scene.js';

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
 * Clicking the image re-traces that one ray on its own. That is the point of the whole
 * interface: CLAUDE.md §22 asks that a viewer be able to inspect the model, the observer,
 * the numerical method and the validation status behind what they are looking at, and a
 * picture with a panel beside it only half answers that. A single line of sight, followed
 * from the camera to wherever it ends, answers it for the pixel actually in question.
 *
 * Rendering happens in Web Workers, which is scheduling only: the image assembled from
 * row blocks is bit-identical to a single-threaded render, and the test suite checks that
 * on every scene kind. This is still the CPU reference path. CLAUDE.md §21 is explicit
 * that 60 FPS is a rendering goal rather than a scientific-validity requirement; the GPU
 * port is its own milestone with its own cross-validation gate.
 */

type SceneKind = SceneDescription['kind'];

interface DisplayEntry {
  readonly label: string;
  readonly value: string;
  readonly note?: string;
  readonly tone?: Tone;
  readonly swatch?: string;
}

const MAX_THREADS = 8;
const TABS = ['model', 'validation', 'ray', 'data'] as const;

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

function plural(n: number, word: string): string {
  return `${n} ${word}${n === 1 ? '' : 's'}`;
}

/**
 * One label / value / note row.
 *
 * The note carries the disclosure CLAUDE.md §22 requires — assumptions, frames, what was
 * omitted — and there is a lot of it. Keeping it one disclosure click away leaves the
 * panel scannable without hiding anything: nothing is dropped, and "Expand all" opens the
 * lot.
 */
function entryElement(entry: DisplayEntry, numeric: boolean): HTMLElement {
  const row = document.createElement('div');
  row.className = entry.tone && entry.tone !== 'neutral' ? `entry ${entry.tone}` : 'entry';

  const label = document.createElement('div');
  label.className = 'entry-label';
  label.textContent = entry.label;

  const value = document.createElement('div');
  value.className = numeric ? 'entry-value numeric' : 'entry-value';
  if (entry.swatch) {
    const swatch = document.createElement('span');
    swatch.className = 'swatch';
    swatch.style.background = entry.swatch;
    value.append(swatch);
  }
  value.append(document.createTextNode(entry.value));

  row.append(label, value);

  if (entry.note) {
    const details = document.createElement('details');
    const summary = document.createElement('summary');
    summary.textContent = 'Why';
    const note = document.createElement('div');
    note.className = 'note';
    note.textContent = entry.note;
    details.append(summary, note);
    row.append(details);
  }
  return row;
}

function renderEntries(container: HTMLElement, entries: readonly DisplayEntry[], numeric = false): void {
  container.replaceChildren(...entries.map((entry) => entryElement(entry, numeric)));
}

async function main(): Promise<void> {
  const canvas = element<HTMLCanvasElement>('view');
  const frame = element<HTMLElement>('frame');
  const rail = element<HTMLElement>('rail');
  const railFill = element<HTMLElement>('rail-fill');
  const status = element<HTMLElement>('status');
  const verdict = element<HTMLElement>('verdict');
  const headline = element<HTMLElement>('headline');
  const modelEntries = element<HTMLElement>('model-entries');
  const validationEntries = element<HTMLElement>('validation-entries');
  const rayEntries = element<HTMLElement>('ray-entries');
  const rayHeadline = element<HTMLElement>('ray-headline');
  const hierarchy = element<HTMLElement>('data-hierarchy');
  const renderButton = element<HTMLButtonElement>('render');
  const stopButton = element<HTMLButtonElement>('stop');

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
  const spinInput = element<HTMLInputElement>('spin');
  const spinValue = element<HTMLElement>('spin-value');

  const threads = threadCount();
  const pool = typeof Worker === 'undefined' ? undefined : new RenderPool(threads);

  let scene: BuiltScene | undefined;
  let selectedPixel: readonly [number, number] | undefined;
  let crosshair: HTMLElement | undefined;

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
    if (kind === 'schwarzschild-sky') return { ...common, kind, cameraRadius: 20, observer };
    if (kind === 'kerr-sky') {
      return {
        ...common,
        kind,
        cameraRadius: 60,
        inclinationDeg: Number.parseFloat(inclinationInput.value),
        spin: Number.parseFloat(spinInput.value),
      };
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

  /** Controls that mean nothing for the chosen scene are hidden rather than left to confuse. */
  function syncControls(): void {
    const kind = sceneSelect.value as SceneKind;
    for (const field of document.querySelectorAll<HTMLElement>('[data-scenes]')) {
      field.hidden = !(field.dataset.scenes ?? '').split(' ').includes(kind);
    }
    inclinationValue.textContent = `${Number.parseFloat(inclinationInput.value).toFixed(0)}°`;
    const spin = Number.parseFloat(spinInput.value);
    spinValue.textContent = `${spin > 0 ? '+' : ''}${spin.toFixed(3)}`;
    const stops = Number.parseFloat(exposureInput.value);
    exposureValue.textContent = `${stops > 0 ? '+' : ''}${stops.toFixed(1)} EV`;
  }

  function selectTab(name: (typeof TABS)[number]): void {
    for (const other of TABS) {
      const tab = element<HTMLElement>(`tab-${other}`);
      const pane = element<HTMLElement>(`pane-${other}`);
      const chosen = other === name;
      tab.setAttribute('aria-selected', String(chosen));
      pane.hidden = !chosen;
    }
  }

  function showInspection(i: number, j: number): void {
    if (!scene) return;
    const inspection = inspectPixel(scene, i, j);
    rayHeadline.textContent = inspection.headline;
    renderEntries(rayEntries, inspection.rows as readonly InspectorRow[], true);
  }

  function placeCrosshair(i: number, j: number): void {
    if (!scene) return;
    crosshair ??= frame.appendChild(Object.assign(document.createElement('div'), { className: 'crosshair' }));
    crosshair.style.left = `${((i + 0.5) / scene.description.widthPx) * 100}%`;
    crosshair.style.top = `${((j + 0.5) / scene.description.heightPx) * 100}%`;
    crosshair.hidden = false;
  }

  function inspect(i: number, j: number): void {
    selectedPixel = [i, j];
    placeCrosshair(i, j);
    showInspection(i, j);
  }

  let controller: AbortController | undefined;

  async function run(): Promise<void> {
    controller?.abort();
    pool?.cancel();
    const local = new AbortController();
    controller = local;

    const description = describeScene();
    scene = buildScene(description);
    const context = canvas.getContext('2d');
    if (!context) throw new Error('main: could not acquire a 2D rendering context.');

    canvas.width = description.widthPx;
    canvas.height = description.heightPx;
    context.clearRect(0, 0, canvas.width, canvas.height);
    if (crosshair) crosshair.hidden = true;

    const rays = description.widthPx * description.heightPx * description.samplesPerAxis ** 2;
    renderButton.disabled = true;
    stopButton.disabled = false;
    rail.dataset.state = 'running';
    railFill.style.width = '0%';
    verdict.textContent = '';
    verdict.className = 'verdict';
    status.textContent = `Tracing ${rays.toLocaleString('en-US')} null geodesics on ${plural(threads, 'thread')}…`;
    headline.textContent = scene.caption;

    const drawRows = (assembly: ImageAssembly, rows: readonly number[]): void => {
      const width = description.widthPx;
      for (const row of rows) {
        const line = assembly.pixels.subarray(row * width * 4, (row + 1) * width * 4);
        context.putImageData(new ImageData(line, width, 1), 0, row);
      }
      const fraction = assembly.rowsCompleted / description.heightPx;
      railFill.style.width = `${(fraction * 100).toFixed(1)}%`;
      status.textContent = `Tracing: ${assembly.rowsCompleted} / ${description.heightPx} rows`;
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

      renderEntries(modelEntries, report.entries as readonly ProvenanceEntry[]);
      renderEntries(validationEntries, report.validation as readonly ProvenanceEntry[], true);
      hierarchy.textContent = report.dataHierarchy;

      const { health } = result.diagnostics;
      railFill.style.width = '100%';
      rail.dataset.state = health.level === 'ok' ? 'done' : 'failed';
      status.textContent =
        `${result.diagnostics.raysTraced.toLocaleString('en-US')} rays · ` +
        `${result.diagnostics.totalSteps.toLocaleString('en-US')} steps · ` +
        `${(elapsedMs / 1000).toFixed(1)} s · ${plural(pool ? threads : 1, 'thread')}`;
      verdict.textContent =
        health.level === 'ok' ? 'All checked invariants within tolerance' : health.messages.join(' ');
      verdict.className = `verdict ${health.level}`;

      if (selectedPixel) {
        const [i, j] = selectedPixel;
        if (i < description.widthPx && j < description.heightPx) inspect(i, j);
        else selectedPixel = undefined;
      }
    } catch (error: unknown) {
      if (local.signal.aborted) return;
      throw error;
    } finally {
      if (controller === local) {
        renderButton.disabled = false;
        stopButton.disabled = true;
      }
    }
  }

  canvas.addEventListener('click', (event) => {
    if (!scene) return;
    const box = canvas.getBoundingClientRect();
    const i = Math.min(
      scene.description.widthPx - 1,
      Math.max(0, Math.floor(((event.clientX - box.left) / box.width) * scene.description.widthPx)),
    );
    const j = Math.min(
      scene.description.heightPx - 1,
      Math.max(0, Math.floor(((event.clientY - box.top) / box.height) * scene.description.heightPx)),
    );
    inspect(i, j);
    selectTab('ray');
  });

  for (const name of TABS) {
    element<HTMLElement>(`tab-${name}`).addEventListener('click', () => selectTab(name));
  }

  for (const [buttonId, containerId] of [
    ['expand-model', 'model-entries'],
    ['expand-validation', 'validation-entries'],
  ] as const) {
    const button = element<HTMLButtonElement>(buttonId);
    button.addEventListener('click', () => {
      const container = element<HTMLElement>(containerId);
      const details = [...container.querySelectorAll('details')];
      const open = !details.every((d) => d.open);
      for (const d of details) d.open = open;
      button.textContent = open ? 'Collapse all' : 'Expand all';
    });
  }

  sceneSelect.addEventListener('change', () => {
    selectedPixel = undefined;
    syncControls();
    void run();
  });
  for (const control of [widthSelect, samplesSelect, observerSelect, massSelect, eddingtonSelect]) {
    control.addEventListener('change', () => void run());
  }
  skyGridInput.addEventListener('change', () => void run());
  for (const slider of [inclinationInput, exposureInput, spinInput]) {
    slider.addEventListener('input', syncControls);
    slider.addEventListener('change', () => void run());
  }
  renderButton.addEventListener('click', () => void run());
  stopButton.addEventListener('click', () => {
    controller?.abort();
    pool?.cancel();
    rail.dataset.state = 'failed';
    status.textContent = 'Stopped. The rows already traced are the ones shown.';
    renderButton.disabled = false;
    stopButton.disabled = true;
  });

  syncControls();
  await run();
}

main().catch((error: unknown) => {
  // A failed physical calculation is reported, never replaced by a plausible-looking
  // image (CLAUDE.md §17, §24).
  const status = document.getElementById('status');
  const rail = document.getElementById('rail');
  const message = error instanceof Error ? error.message : String(error);
  if (status) status.textContent = `Render failed: ${message}`;
  if (rail) rail.dataset.state = 'failed';
  console.error(error);
});
