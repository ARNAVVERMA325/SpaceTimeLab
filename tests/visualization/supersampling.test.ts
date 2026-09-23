import { describe, expect, it } from 'vitest';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import {
  defaultScreen,
  localRayDirectionAt,
  staticMinkowskiObserver,
} from '../../src/physics/observer/observer.js';
import { minkowski } from '../../src/physics/spacetimes/minkowski.js';
import { DEFAULT_CELESTIAL_GRID, sampleCelestialGrid } from '../../src/visualization/celestial-grid.js';
import { decodeSrgb, hashToUnit, linearToSrgb8, srgbToLinear } from '../../src/visualization/color.js';
import { FAILED_RAY_COLOR, renderImage, type TraceConfig } from '../../src/visualization/raytracer.js';

/**
 * Supersampling (CLAUDE.md §9 step 7), and the linear-light accumulation it depends on.
 */

function flatConfig(samplesPerAxis: number, seed = 0): TraceConfig {
  return {
    model: minkowski,
    integrator: new Dopri5Integrator(STATE_DIM),
    observer: staticMinkowskiObserver([0, 0, 0, 0]),
    grid: { ...DEFAULT_CELESTIAL_GRID, radius: 100 },
    limits: { initialStep: 0.5, maxStep: 10, parameterMax: 1000, maxSteps: 10_000 },
    sampling: { samplesPerAxis, seed },
  };
}

describe('sRGB transfer functions', () => {
  it('round-trips every 8-bit value exactly', () => {
    for (let c = 0; c < 256; c += 1) expect(linearToSrgb8(srgbToLinear(c))).toBe(c);
  });

  it('averages in linear light, not in encoded values', () => {
    // Half black and half white is 50% luminance, which encodes to sRGB 188 — not the
    // 128 that averaging the encoded values would give.
    const mean = (srgbToLinear(0) + srgbToLinear(255)) / 2;
    expect(linearToSrgb8(mean)).toBe(188);
  });
});

describe('deterministic jitter', () => {
  it('is reproducible and lies in [0, 1)', () => {
    for (let i = 0; i < 1000; i += 1) {
      const u = hashToUnit(i, 7, 3, 42);
      expect(u).toBe(hashToUnit(i, 7, 3, 42));
      expect(u).toBeGreaterThanOrEqual(0);
      expect(u).toBeLessThan(1);
    }
  });

  it('is roughly uniform', () => {
    const bins = new Array(10).fill(0);
    const n = 20_000;
    for (let i = 0; i < n; i += 1) bins[Math.floor(hashToUnit(i, i >> 3, i & 7, 1) * 10)] += 1;
    // Each decile within 5% of n/10: loose, since this only needs to be a decent jitter.
    for (const count of bins) expect(Math.abs(count - n / 10) / (n / 10)).toBeLessThan(0.05);
  });
});

describe('stratified supersampling', () => {
  it('produces bit-identical images for the same seed', () => {
    const screen = defaultScreen(16, 12, Math.PI / 2);
    const a = renderImage(flatConfig(3, 9), screen);
    const b = renderImage(flatConfig(3, 9), screen);
    expect(Array.from(a.pixels)).toEqual(Array.from(b.pixels));
    const c = renderImage(flatConfig(3, 10), screen);
    expect(Array.from(a.pixels)).not.toEqual(Array.from(c.pixels));
  });

  it('counts every sample as a traced ray', () => {
    const screen = defaultScreen(10, 8, Math.PI / 2);
    const result = renderImage(flatConfig(4), screen);
    expect(result.diagnostics.samplesPerPixel).toBe(16);
    expect(result.diagnostics.raysTraced).toBe(10 * 8 * 16);
    expect(result.diagnostics.raysReachingBackground).toBe(10 * 8 * 16);
  });

  it('converges to the exact pixel average at the stratified-sampling rate', () => {
    // The reference is independent of the renderer: in flat space the traced image is
    // already proven to equal the grid sampled along each initial viewing direction, so
    // the exact pixel value is a quadrature of the grid over the pixel footprint, done
    // here on a dense 48 x 48 lattice per pixel without tracing a single ray.
    //
    // For images with edges, stratified jitter converges as O(N^-3/4) in samples per pixel
    // N, so each doubling of samples per axis should cut the error by about 2^1.5 ~ 2.8
    // (Mitchell, "Consequences of stratified sampling in graphics", SIGGRAPH 1996).
    const screen = defaultScreen(32, 24, Math.PI / 2);
    const lattice = 48;
    const reference = new Float64Array(screen.widthPx * screen.heightPx * 3);
    for (let j = 0; j < screen.heightPx; j += 1) {
      for (let i = 0; i < screen.widthPx; i += 1) {
        let r = 0;
        let g = 0;
        let b = 0;
        for (let sy = 0; sy < lattice; sy += 1) {
          for (let sx = 0; sx < lattice; sx += 1) {
            const d = localRayDirectionAt(screen, i + (sx + 0.5) / lattice, j + (sy + 0.5) / lattice);
            const color = decodeSrgb(sampleCelestialGrid(DEFAULT_CELESTIAL_GRID, d[0], d[1], d[2]));
            r += color.r;
            g += color.g;
            b += color.b;
          }
        }
        const k = (j * screen.widthPx + i) * 3;
        const count = lattice * lattice;
        reference[k] = r / count;
        reference[k + 1] = g / count;
        reference[k + 2] = b / count;
      }
    }

    const errors = [1, 2, 4, 8].map((n) => {
      const image = renderImage(flatConfig(n), screen).linear;
      let sum = 0;
      for (let i = 0; i < reference.length; i += 1) sum += Math.abs(image[i] - reference[i]);
      return sum / reference.length;
    });
    for (let i = 1; i < errors.length; i += 1) {
      const ratio = errors[i - 1] / errors[i];
      expect(ratio, `errors ${errors.map((e) => e.toExponential(2)).join(', ')}`).toBeGreaterThan(2);
      expect(ratio).toBeLessThan(4);
    }
    expect(errors[3]).toBeLessThan(errors[0] / 12);
  });

  it('paints a pixel with any failed sample in the failure colour, never a blend', () => {
    // A budget too small to reach the background makes every sample fail.
    const starved: TraceConfig = {
      ...flatConfig(2),
      limits: { initialStep: 0.5, maxStep: 0.5, parameterMax: 1, maxSteps: 10 },
    };
    const result = renderImage(starved, defaultScreen(4, 4, Math.PI / 2));
    expect(result.diagnostics.pixelsWithFailures).toBe(16);
    expect(result.diagnostics.health.level).toBe('failed');
    for (let p = 0; p < result.pixels.length; p += 4) {
      expect([result.pixels[p], result.pixels[p + 1], result.pixels[p + 2]]).toEqual([
        FAILED_RAY_COLOR.r,
        FAILED_RAY_COLOR.g,
        FAILED_RAY_COLOR.b,
      ]);
    }
  });
});
