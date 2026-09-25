import { NULL_NORMALIZATION_TRACED } from '../physics/validation/tolerances.js';
import {
  createAccumulator,
  diagnosticsFromStats,
  emptyStats,
  mergeStats,
  renderBand,
  type RenderResult,
  type RenderStats,
} from './raytracer.js';
import { buildScene, type BuiltScene, type SceneDescription } from './scene.js';

/**
 * Rendering an arbitrary set of rows, and reassembling an image from such sets.
 *
 * This is the whole of the parallel renderer's correctness argument. Each ray depends only
 * on its pixel, its sample index and the seed; rows share nothing. So rendering any
 * partition of the rows and reassembling the pieces must reproduce a single-pass render
 * exactly, bit for bit — and `tests/visualization/parallel-render.test.ts` checks that it
 * does, on every kind of scene. The workers then only move these blocks between threads.
 */

export interface RowBlock {
  readonly rows: readonly number[];
  /** RGBA8 for each listed row, in the order listed. */
  readonly pixels: Uint8ClampedArray<ArrayBuffer>;
  /** Linear RGB for each listed row, in the order listed. */
  readonly linear: Float64Array<ArrayBuffer>;
  readonly stats: RenderStats;
}

export function renderRows(scene: BuiltScene, rows: readonly number[]): RowBlock {
  const { config, screen } = scene;
  const accumulator = createAccumulator(config, screen);
  for (const row of rows) renderBand(config, screen, accumulator, row, row + 1);

  const width = screen.widthPx;
  const pixels = new Uint8ClampedArray(rows.length * width * 4);
  const linear = new Float64Array(rows.length * width * 3);
  rows.forEach((row, index) => {
    pixels.set(accumulator.pixels.subarray(row * width * 4, (row + 1) * width * 4), index * width * 4);
    linear.set(accumulator.linear.subarray(row * width * 3, (row + 1) * width * 3), index * width * 3);
  });

  const stats = emptyStats();
  mergeStats(stats, accumulator);
  return { rows, pixels, linear, stats };
}

/** An image being assembled from row blocks as they arrive. */
export class ImageAssembly {
  readonly pixels: Uint8ClampedArray<ArrayBuffer>;
  readonly linear: Float64Array<ArrayBuffer>;
  readonly stats = emptyStats();
  private rowsDone = 0;

  constructor(readonly description: SceneDescription) {
    this.pixels = new Uint8ClampedArray(description.widthPx * description.heightPx * 4);
    this.linear = new Float64Array(description.widthPx * description.heightPx * 3);
  }

  get rowsCompleted(): number {
    return this.rowsDone;
  }

  add(block: RowBlock): void {
    const width = this.description.widthPx;
    block.rows.forEach((row, index) => {
      this.pixels.set(block.pixels.subarray(index * width * 4, (index + 1) * width * 4), row * width * 4);
      this.linear.set(block.linear.subarray(index * width * 3, (index + 1) * width * 3), row * width * 3);
    });
    mergeStats(this.stats, block.stats);
    this.rowsDone += block.rows.length;
  }

  result(): RenderResult {
    const { widthPx, heightPx, samplesPerAxis } = this.description;
    return {
      widthPx,
      heightPx,
      pixels: this.pixels,
      linear: this.linear,
      diagnostics: diagnosticsFromStats(
        this.stats,
        widthPx * heightPx,
        samplesPerAxis,
        NULL_NORMALIZATION_TRACED,
        // Only the Schwarzschild scenes carry the exact asymptotic correction; Kerr has no
        // closed-form tail integral, and flat spacetime needs none.
        this.description.kind === 'schwarzschild-sky' || this.description.kind === 'schwarzschild-disk',
      ),
    };
  }
}

/**
 * Deal rows out round-robin in small chunks.
 *
 * Interleaving rather than contiguous bands, because the cost of a row is far from
 * uniform: rays near the shadow edge loop around the photon sphere and take many times
 * more steps than rays that miss the hole. Contiguous bands would hand the whole expensive
 * middle of the image to one or two workers.
 */
export function interleaveRows(heightPx: number, workers: number, chunk = 2): number[][] {
  const assignments: number[][] = Array.from({ length: workers }, () => []);
  for (let start = 0, c = 0; start < heightPx; start += chunk, c += 1) {
    for (let row = start; row < Math.min(heightPx, start + chunk); row += 1) {
      assignments[c % workers].push(row);
    }
  }
  return assignments;
}

export { buildScene };
