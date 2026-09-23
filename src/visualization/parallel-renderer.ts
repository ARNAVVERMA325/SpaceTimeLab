import type { RenderJob, WorkerMessage } from '../workers/render-worker.js';
import type { RenderResult } from './raytracer.js';
import { buildScene, ImageAssembly, interleaveRows, renderRows } from './render-rows.js';
import type { SceneDescription } from './scene.js';

/**
 * Rendering across all available cores.
 *
 * Rows are dealt round-robin to a pool of Web Workers; each posts finished rows back as it
 * goes, and they are drawn immediately. Parallelism changes no trajectory: the assembled
 * image is bit-identical to a single-threaded one, which the test suite checks directly on
 * the same row-partitioning code the workers run.
 *
 * Cancelling a render terminates the workers rather than asking them to stop, since a
 * worker in the middle of a long ray cannot be interrupted any other way; the pool is
 * rebuilt on the next render.
 */

export interface ParallelRenderOptions {
  readonly onRows?: (assembly: ImageAssembly, rows: readonly number[]) => void;
  readonly signal?: AbortSignal;
}

export class RenderPool {
  private workers: Worker[] = [];
  private nextJob = 1;

  constructor(readonly size: number) {}

  private ensureWorkers(): Worker[] {
    if (this.workers.length === 0) {
      this.workers = Array.from(
        { length: this.size },
        () => new Worker(new URL('../workers/render-worker.ts', import.meta.url), { type: 'module' }),
      );
    }
    return this.workers;
  }

  /** Stop everything in flight. The next render starts a fresh pool. */
  cancel(): void {
    for (const worker of this.workers) worker.terminate();
    this.workers = [];
  }

  render(description: SceneDescription, options: ParallelRenderOptions = {}): Promise<RenderResult> {
    const workers = this.ensureWorkers();
    const jobId = this.nextJob++;
    const assembly = new ImageAssembly(description);
    const assignments = interleaveRows(description.heightPx, workers.length);

    return new Promise<RenderResult>((resolve, reject) => {
      let remaining = workers.length;
      const onAbort = (): void => {
        this.cancel();
        reject(new DOMException('Render aborted.', 'AbortError'));
      };
      options.signal?.addEventListener('abort', onAbort, { once: true });

      workers.forEach((worker, index) => {
        worker.onmessage = (event: MessageEvent<WorkerMessage>) => {
          const message = event.data;
          if (message.jobId !== jobId) return;
          if (message.type === 'block') {
            assembly.add(message.block);
            options.onRows?.(assembly, message.block.rows);
          } else if (message.type === 'error') {
            options.signal?.removeEventListener('abort', onAbort);
            this.cancel();
            reject(new Error(`Render worker failed: ${message.message}`));
          } else if (--remaining === 0) {
            options.signal?.removeEventListener('abort', onAbort);
            resolve(assembly.result());
          }
        };
        worker.onerror = (event) => {
          options.signal?.removeEventListener('abort', onAbort);
          this.cancel();
          reject(new Error(`Render worker crashed: ${event.message}`));
        };
        const job: RenderJob = { jobId, description, rows: assignments[index], chunk: 2 };
        worker.postMessage(job);
      });
    });
  }
}

/**
 * Single-threaded fallback with the same assembly path, for environments without Workers.
 * Yields between chunks so a page stays responsive.
 */
export async function renderSerially(
  description: SceneDescription,
  options: ParallelRenderOptions = {},
): Promise<RenderResult> {
  const scene = buildScene(description);
  const assembly = new ImageAssembly(description);
  for (let row = 0; row < description.heightPx; row += 4) {
    if (options.signal?.aborted) throw new DOMException('Render aborted.', 'AbortError');
    const rows = Array.from({ length: Math.min(4, description.heightPx - row) }, (_, i) => row + i);
    assembly.add(renderRows(scene, rows));
    options.onRows?.(assembly, rows);
    await new Promise((resolve) => setTimeout(resolve, 0));
  }
  return assembly.result();
}
