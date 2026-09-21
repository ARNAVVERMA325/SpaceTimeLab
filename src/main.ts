import type { Vec4 } from './physics/core/indices.js';
import type { Integrator } from './physics/geodesic/integrators/integrator.js';
import { RKF45Integrator } from './physics/geodesic/integrators/rkf45.js';
import { STATE_DIM } from './physics/geodesic/state-vector.js';
import {
  defaultScreen,
  inwardFacingScreen,
  staticMinkowskiObserver,
  staticObserver,
  type Observer,
  type PinholeScreen,
} from './physics/observer/observer.js';
import { minkowski } from './physics/spacetimes/minkowski.js';
import {
  criticalImpactParameter,
  photonSphereRadius,
  schwarzschild,
} from './physics/spacetimes/schwarzschild.js';
import { isCaptured } from './physics/spacetimes/schwarzschild-rays.js';
import type { SpacetimeModel } from './physics/spacetimes/spacetime-model.js';
import {
  NULL_NORMALIZATION_PREVIEW,
  NULL_NORMALIZATION_TRACED,
} from './physics/validation/tolerances.js';
import { buildProvenanceReport, type ProvenanceEntry } from './ui/provenance.js';
import { renderToCanvas } from './visualization/canvas-renderer.js';
import { DEFAULT_CELESTIAL_GRID } from './visualization/celestial-grid.js';
import type { TraceConfig } from './visualization/raytracer.js';

/**
 * Spacetime Lab application entry point.
 *
 * Two scenes share one pipeline: observer tetrad, ray generation, null-geodesic
 * integration, background intersection. Switching between them changes the metric and
 * nothing else, which is the point — flat spacetime bends nothing and the sky grid
 * arrives undistorted, while Schwarzschild bends light into a shadow and an Einstein
 * ring using the same code.
 *
 * The interactive render is the CPU reference path. It is not fast, and CLAUDE.md §21
 * is explicit that 60 FPS is a rendering goal rather than a scientific-validity
 * requirement; the GPU port is its own milestone with its own cross-validation gate.
 */

interface Scene {
  readonly id: string;
  readonly label: string;
  readonly model: SpacetimeModel;
  readonly observer: Observer;
  readonly screen: (widthPx: number, heightPx: number) => PinholeScreen;
  readonly config: Omit<TraceConfig, 'observer'>;
  readonly entries: readonly ProvenanceEntry[];
  readonly caption: string;
}

/**
 * Preview integration tolerance.
 *
 * Looser than the test suite's reference setting, and paired with
 * NULL_NORMALIZATION_PREVIEW so the reported health describes the preview rather than
 * failing against a gate it was never asked to meet. This lowers a numerical budget,
 * not the physical model (CLAUDE.md §1.1, §21).
 */
const PREVIEW_TOLERANCE = { absolute: 1e-10, relative: 1e-10 };

const BLACK_HOLE_MASS = 1;
const CAMERA_RADIUS = 20;

function buildScenes(): readonly Scene[] {
  const integrator: Integrator = new RKF45Integrator(STATE_DIM, {
    tolerance: PREVIEW_TOLERANCE,
  });

  const flatObserver = staticMinkowskiObserver([0, 0, 0, 0]);
  const flatGrid = { ...DEFAULT_CELESTIAL_GRID, radius: 100 };

  const holeModel = schwarzschild(BLACK_HOLE_MASS);
  const holeObserver = staticObserver(holeModel, [0, CAMERA_RADIUS, Math.PI / 2, 0]);
  const holeGrid = { ...DEFAULT_CELESTIAL_GRID, radius: 200 };

  // The shadow's angular radius for a static observer at r: a photon arriving at angle
  // psi from the inward radial direction has b = r sin(psi) / sqrt(f), so the capture
  // boundary b = b_c sits at sin(psi) = b_c sqrt(f) / r.
  const f = 1 - (2 * BLACK_HOLE_MASS) / CAMERA_RADIUS;
  const shadowAngle = Math.asin(
    (criticalImpactParameter(BLACK_HOLE_MASS) * Math.sqrt(f)) / CAMERA_RADIUS,
  );

  return [
    {
      id: 'schwarzschild',
      label: 'Schwarzschild black hole',
      model: holeModel,
      observer: holeObserver,
      screen: (w, h) => inwardFacingScreen(w, h, 3.2 * shadowAngle),
      config: {
        model: holeModel,
        integrator,
        grid: holeGrid,
        limits: { initialStep: 1e-3, parameterMax: 20_000, maxSteps: 200_000, maxStep: 5 },
        captureTest: (position_x: Vec4, tangent: Vec4) => isCaptured(holeModel, position_x, tangent),
        orbitalPlaneReduction: true,
        residualTolerance: NULL_NORMALIZATION_PREVIEW,
      },
      caption:
        'Computed appearance of the background grid for the selected spacetime, observer ' +
        'and rendering assumptions. The dark disc is the black-hole shadow: lines of ' +
        'sight along which no background light reaches the observer.',
      entries: [
        {
          label: 'Mass parameter',
          value: `M = ${BLACK_HOLE_MASS} (geometric units)`,
          note: `Horizon at r = ${2 * BLACK_HOLE_MASS}, photon sphere at r = ${photonSphereRadius(BLACK_HOLE_MASS)}.`,
        },
        {
          label: 'Camera position',
          value: `Static observer at r = ${CAMERA_RADIUS}M, theta = pi/2`,
          note:
            'Hovering at fixed radius, which requires proper acceleration. A freely-falling ' +
            'observer at the same event would see a different image; that comparison is ' +
            'Milestone 3.',
        },
        {
          label: 'Shadow',
          value: `Angular radius ${(shadowAngle * (180 / Math.PI)).toFixed(2)} degrees`,
          note:
            'Predicted in closed form by sin(psi) = b_c sqrt(f) / r with ' +
            `b_c = 3 sqrt(3) M ~ ${criticalImpactParameter(BLACK_HOLE_MASS).toFixed(4)}M, and ` +
            'matched by the traced image to better than a part in a million. The shadow is ' +
            'larger than the horizon and is not a picture of it: it is the set of directions ' +
            'whose backward-traced rays end on the hole.',
        },
        {
          label: 'Ray termination',
          value: 'Captured when r < 3M with k^r < 0',
          note:
            'Exact rather than a tuned cutoff. The null effective potential f/r^2 increases ' +
            'inward of r = 3M, so a photon moving inward there can never turn around. Rays ' +
            'are stopped while the chart is still well behaved rather than integrated toward ' +
            'r = 2M, where these coordinates break down.',
        },
        {
          label: 'Sampling and aliasing',
          value: 'One ray per pixel, no anti-aliasing',
          note:
            'The stippling in the fine bands hugging the shadow is aliasing, and it is ' +
            'showing you something real. Approaching the capture boundary, the lensing map ' +
            'compresses an unbounded sequence of images of the whole sky into a vanishing ' +
            'angular width, so no finite ray count can resolve it. That is a sampling limit ' +
            'of this render, not an error in the trajectories. No smoothing is applied: ' +
            'anti-aliasing is a rendering operation and must not be allowed to stand in for ' +
            'resolving the structure (CLAUDE.md §9).',
        },
        {
          label: 'Integration domain',
          value: 'Exterior only, r > 2M',
          note:
            'Schwarzschild coordinates do not cover the horizon. That is a property of the ' +
            'chart, not of the spacetime: the Kretschmann scalar K = 48 M^2 / r^6 is finite ' +
            'at r = 2M. Continuing through the horizon needs horizon-penetrating coordinates ' +
            'and is a later milestone.',
        },
      ],
    },
    {
      id: 'minkowski',
      label: 'Minkowski (flat) baseline',
      model: minkowski,
      observer: flatObserver,
      screen: (w, h) => defaultScreen(w, h, Math.PI / 2),
      config: {
        model: minkowski,
        integrator,
        grid: flatGrid,
        limits: { initialStep: 0.5, maxStep: 10, parameterMax: 1000, maxSteps: 10_000 },
        residualTolerance: NULL_NORMALIZATION_TRACED,
      },
      caption:
        'Flat spacetime: the pipeline bends nothing, so the sky grid arrives exactly as a ' +
        'pinhole camera projects it. Curved grid lines here are rectilinear projection of a ' +
        'sphere, not light deflection.',
      entries: [
        {
          label: 'Image projection',
          value: 'Rectilinear pinhole projection of the celestial sphere',
          note:
            'A rectilinear camera maps great circles to straight lines, so meridians appear ' +
            'straight while parallels do not. In flat spacetime this image is identical, pixel ' +
            'for pixel, to sampling the grid along each initial viewing direction with no ' +
            'integration at all, which is what the validation suite asserts.',
        },
      ],
    },
  ];
}

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
  const sceneSelect = element<HTMLSelectElement>('scene');
  const resolutionSelect = element<HTMLSelectElement>('resolution');
  const renderButton = element<HTMLButtonElement>('render');

  const scenes = buildScenes();
  sceneSelect.replaceChildren(
    ...scenes.map((scene) => new Option(scene.label, scene.id)),
  );

  let controller: AbortController | undefined;

  async function run(): Promise<void> {
    controller?.abort();
    const local = new AbortController();
    controller = local;

    const scene = scenes.find((s) => s.id === sceneSelect.value) ?? scenes[0];
    const width = Number.parseInt(resolutionSelect.value, 10);
    const screen = scene.screen(width, Math.round((width * 3) / 4));
    const config: TraceConfig = { ...scene.config, observer: scene.observer };

    renderButton.disabled = true;
    status.dataset.health = '';
    status.textContent = `Tracing ${screen.widthPx * screen.heightPx} null geodesics...`;
    headline.textContent = scene.caption;

    const started = performance.now();
    try {
      const result = await renderToCanvas(canvas, {
        config,
        screen,
        rowsPerBand: 4,
        signal: local.signal,
        onProgress: (rows, total) => {
          status.textContent = `Tracing null geodesics: ${rows} / ${total} rows`;
        },
      });
      const elapsedMs = performance.now() - started;

      const report = buildProvenanceReport({
        model: scene.model,
        observer: scene.observer,
        integrator: scene.config.integrator,
        diagnostics: result.diagnostics,
        backgroundRadius: scene.config.grid.radius,
        sceneEntries: scene.entries,
      });

      renderEntries(modelEntries, report.entries);
      renderEntries(validationEntries, report.validation);
      hierarchy.textContent = report.dataHierarchy;

      const health = result.diagnostics.health;
      status.textContent =
        `${result.diagnostics.raysTraced} rays, ` +
        `${result.diagnostics.totalSteps.toLocaleString('en-US')} integration steps, ` +
        `${(elapsedMs / 1000).toFixed(1)} s. ` +
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

  sceneSelect.addEventListener('change', () => void run());
  resolutionSelect.addEventListener('change', () => void run());
  renderButton.addEventListener('click', () => void run());

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
