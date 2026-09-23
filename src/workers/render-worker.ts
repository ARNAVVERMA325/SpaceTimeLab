import { buildScene, renderRows, type RowBlock } from '../visualization/render-rows.js';
import type { BuiltScene, SceneDescription } from '../visualization/scene.js';

/**
 * A render worker: builds the scene from its plain-data description and traces the rows
 * it is dealt, posting each finished chunk back as soon as it is done.
 *
 * It holds no state that could make it disagree with the main thread: the scene is
 * rebuilt from the same description by the same function, and every ray depends only on
 * its pixel, sample and seed.
 */

export interface RenderJob {
  readonly jobId: number;
  readonly description: SceneDescription;
  readonly rows: readonly number[];
  /** Rows per posted block: small enough to show progress, large enough to amortize messaging. */
  readonly chunk: number;
}

export type WorkerMessage =
  | { readonly jobId: number; readonly type: 'block'; readonly block: RowBlock }
  | { readonly jobId: number; readonly type: 'done' }
  | { readonly jobId: number; readonly type: 'error'; readonly message: string };

interface WorkerScope {
  onmessage: ((event: MessageEvent<RenderJob>) => void) | null;
  postMessage(message: WorkerMessage, transfer?: Transferable[]): void;
}

const scope = self as unknown as WorkerScope;
let cached: { readonly key: string; readonly scene: BuiltScene } | undefined;

scope.onmessage = (event) => {
  const { jobId, description, rows, chunk } = event.data;
  try {
    const key = JSON.stringify(description);
    if (!cached || cached.key !== key) cached = { key, scene: buildScene(description) };
    for (let i = 0; i < rows.length; i += chunk) {
      const block = renderRows(cached.scene, rows.slice(i, i + chunk));
      scope.postMessage({ jobId, type: 'block', block }, [block.pixels.buffer, block.linear.buffer]);
    }
    scope.postMessage({ jobId, type: 'done' });
  } catch (error) {
    scope.postMessage({ jobId, type: 'error', message: error instanceof Error ? error.message : String(error) });
  }
};
