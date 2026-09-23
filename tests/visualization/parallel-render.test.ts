import { describe, expect, it } from 'vitest';
import { renderImage } from '../../src/visualization/raytracer.js';
import {
  buildScene,
  ImageAssembly,
  interleaveRows,
  renderRows,
} from '../../src/visualization/render-rows.js';
import type { SceneDescription } from '../../src/visualization/scene.js';

/**
 * The parallel renderer's correctness argument, tested directly: rendering any partition
 * of the rows and reassembling must reproduce a single-pass render bit for bit. The Web
 * Workers only carry these row blocks between threads, so this is what makes a parallel
 * image trustworthy — not that it looks the same, but that it is the same.
 */

const common = { widthPx: 40, heightPx: 24, samplesPerAxis: 2, seed: 7 } as const;
const SCENES: readonly SceneDescription[] = [
  { ...common, kind: 'minkowski' },
  { ...common, kind: 'schwarzschild-sky', cameraRadius: 20, observer: 'static' },
  { ...common, kind: 'schwarzschild-sky', cameraRadius: 20, observer: 'free-fall' },
  {
    ...common,
    kind: 'schwarzschild-disk',
    cameraRadius: 60,
    inclinationDeg: 75,
    observer: 'static',
    massSolar: 4e10,
    eddingtonFraction: 0.01,
    exposureStops: 0,
    skyGrid: false,
  },
];

describe('parallel rendering reproduces a serial render exactly', () => {
  for (const description of SCENES) {
    it(`${description.kind}${'observer' in description ? ` (${description.observer})` : ''}`, () => {
      const serial = renderImage(buildScene(description).config, buildScene(description).screen);

      // Four simulated workers, each with its own freshly built scene, as in the browser.
      const assembly = new ImageAssembly(description);
      for (const rows of interleaveRows(description.heightPx, 4)) {
        assembly.add(renderRows(buildScene(description), rows));
      }
      const parallel = assembly.result();

      expect(assembly.rowsCompleted).toBe(description.heightPx);
      expect(Array.from(parallel.pixels)).toEqual(Array.from(serial.pixels));
      expect(Array.from(parallel.linear)).toEqual(Array.from(serial.linear));
      const a = parallel.diagnostics;
      const b = serial.diagnostics;
      expect(a.raysTraced).toBe(b.raysTraced);
      expect(a.raysCaptured).toBe(b.raysCaptured);
      expect(a.raysHittingDisk).toBe(b.raysHittingDisk);
      expect(a.raysReachingBackground).toBe(b.raysReachingBackground);
      expect(a.totalSteps).toBe(b.totalSteps);
      expect(a.maxNullResidual).toBe(b.maxNullResidual);
      expect(a.health.level).toBe(b.health.level);
    });
  }
});

describe('row interleaving', () => {
  it('assigns every row exactly once', () => {
    for (const [height, workers] of [[24, 4], [37, 3], [5, 8], [1, 2]] as const) {
      const all = interleaveRows(height, workers).flat().sort((x, y) => x - y);
      expect(all).toEqual(Array.from({ length: height }, (_, i) => i));
    }
  });

  it('spreads the expensive middle rows across workers rather than giving them to one', () => {
    const assignments = interleaveRows(120, 4);
    const middle = new Set(Array.from({ length: 20 }, (_, i) => 50 + i));
    const counts = assignments.map((rows) => rows.filter((r) => middle.has(r)).length);
    expect(Math.max(...counts) - Math.min(...counts)).toBeLessThanOrEqual(2);
  });
});
