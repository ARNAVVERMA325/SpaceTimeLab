import type { PinholeScreen } from '../physics/observer/observer.js';
import {
  createAccumulator,
  finalizeRender,
  renderBand,
  type RenderResult,
  type TraceConfig,
} from './raytracer.js';

/**
 * Visualization layer (CLAUDE.md §20): drawing a traced image onto a Canvas2D surface.
 *
 * Rendering proceeds in row bands with a yield between them, so the page stays
 * responsive while a full image is traced. Banding is a scheduling detail only: each
 * ray is traced independently and the physical result is identical to a single pass.
 *
 * CLAUDE.md §9 is explicit that graphics techniques must not alter the underlying
 * physical trajectory. Nothing here touches the geodesics; the canvas receives the
 * traced colours unmodified, with no smoothing, interpolation or denoising.
 */

export interface CanvasRenderOptions {
  readonly config: TraceConfig;
  readonly screen: PinholeScreen;
  readonly rowsPerBand?: number;
  readonly onProgress?: (rowsCompleted: number, totalRows: number) => void;
  readonly signal?: AbortSignal;
}

function nextFrame(): Promise<void> {
  return new Promise((resolve) => {
    if (typeof requestAnimationFrame === 'function') {
      requestAnimationFrame(() => resolve());
    } else {
      setTimeout(resolve, 0);
    }
  });
}

/**
 * Trace an image band by band, blitting each completed band to the canvas.
 *
 * Resolves with the finished render and its diagnostics, or rejects if the caller
 * aborts.
 */
export async function renderToCanvas(
  canvas: HTMLCanvasElement,
  options: CanvasRenderOptions,
): Promise<RenderResult> {
  const { config, screen, onProgress, signal } = options;
  const rowsPerBand = Math.max(1, options.rowsPerBand ?? 16);

  canvas.width = screen.widthPx;
  canvas.height = screen.heightPx;

  const context = canvas.getContext('2d');
  if (!context) {
    throw new Error('renderToCanvas: could not acquire a 2D rendering context.');
  }

  const accumulator = createAccumulator(config, screen);

  for (let rowStart = 0; rowStart < screen.heightPx; rowStart += rowsPerBand) {
    if (signal?.aborted) {
      throw new DOMException('Render aborted.', 'AbortError');
    }

    const rowEnd = Math.min(rowStart + rowsPerBand, screen.heightPx);
    renderBand(config, screen, accumulator, rowStart, rowEnd);

    const bandHeight = rowEnd - rowStart;
    const bandPixels = accumulator.pixels.subarray(
      rowStart * screen.widthPx * 4,
      rowEnd * screen.widthPx * 4,
    );
    context.putImageData(new ImageData(bandPixels, screen.widthPx, bandHeight), 0, rowStart);

    onProgress?.(rowEnd, screen.heightPx);
    await nextFrame();
  }

  return finalizeRender(accumulator, config);
}
