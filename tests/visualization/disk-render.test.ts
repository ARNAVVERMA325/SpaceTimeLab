import { describe, expect, it } from 'vitest';
import { Dopri5Integrator } from '../../src/physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../../src/physics/geodesic/state-vector.js';
import { inwardFacingScreen, staticObserver } from '../../src/physics/observer/observer.js';
import { thinDiskScale } from '../../src/physics/spacetimes/novikov-thorne.js';
import { schwarzschild } from '../../src/physics/spacetimes/schwarzschild.js';
import { asymptoticDirection, isCaptured } from '../../src/physics/spacetimes/schwarzschild-rays.js';
import { DEFAULT_CELESTIAL_GRID } from '../../src/visualization/celestial-grid.js';
import { diskReferenceLuminance, displayGain, toneMapReinhard } from '../../src/visualization/disk-emission.js';
import { renderImage, type RenderResult } from '../../src/visualization/raytracer.js';

/**
 * ROADMAP.md 3.4 — the thin accretion disk as rendered.
 *
 * The image features asserted here are not drawn: each follows from the geodesics, the
 * emitter's motion and the blackbody. The far side of the disk appears above the shadow
 * because rays passing over the hole are bent down onto it; the approaching side is
 * brighter because g > 1 there.
 */

const M = 1;
const model = schwarzschild(M);

function renderDisk(inclinationDeg: number, massSolar: number, eddingtonFraction: number, width = 96): RenderResult {
  const disk = { mass: M, innerRadius: 6, outerRadius: 20, scale: thinDiskScale({ massSolar, eddingtonFraction }) };
  const observer = staticObserver(model, [0, 60, (inclinationDeg * Math.PI) / 180, 0]);
  return renderImage(
    {
      model,
      observer,
      integrator: new Dopri5Integrator(STATE_DIM, { tolerance: { absolute: 1e-10, relative: 1e-10 } }),
      grid: { ...DEFAULT_CELESTIAL_GRID, radius: 200 },
      limits: { initialStep: 1e-3, parameterMax: 20_000, maxSteps: 200_000, maxStep: 5 },
      captureTest: (x, k) => isCaptured(model, x, k),
      orbitalPlaneReduction: true,
      asymptoticDirection: (p, d) => asymptoticDirection(model, p, d),
      disk,
      display: { referenceLuminance: diskReferenceLuminance(disk), exposureStops: 0, toneMap: 'none' },
      background: 'black',
    },
    inwardFacingScreen(width, Math.round(width * 0.5625), 0.9),
  );
}

/** Mean display-linear luminance over a pixel region, counting only lit pixels. */
function meanLuminance(result: RenderResult, x0: number, x1: number, y0: number, y1: number): number {
  let sum = 0;
  let count = 0;
  for (let y = y0; y < y1; y += 1) {
    for (let x = x0; x < x1; x += 1) {
      const k = (y * result.widthPx + x) * 3;
      const Y = 0.2126 * result.linear[k] + 0.7152 * result.linear[k + 1] + 0.0722 * result.linear[k + 2];
      if (Y > 1e-6) {
        sum += Y;
        count += 1;
      }
    }
  }
  return count > 0 ? sum / count : 0;
}

describe('rendered thin disk (ROADMAP.md 3.4)', () => {
  const edgeOn = renderDisk(80, 4e10, 0.01);

  it('renders within the validated tolerance, with no failed rays', () => {
    expect(edgeOn.diagnostics.raysFailed).toBe(0);
    expect(edgeOn.diagnostics.raysHittingDisk).toBeGreaterThan(0);
    expect(edgeOn.diagnostics.health.level).toBe('ok');
  });

  it('shows both redshift and blueshift across a steeply inclined disk', () => {
    const [gMin, gMax] = edgeOn.diagnostics.frequencyRatioRange;
    expect(gMin).toBeLessThan(0.6);
    expect(gMax).toBeGreaterThan(1.3);
  });

  it('makes the approaching side brighter: Doppler beaming', () => {
    // Prograde rotation about +z, camera at phi = 0 looking in: the gas on the image's
    // right moves toward the camera. At 8830 K the visible band is near the Wien peak and
    // the asymmetry is strong — measured about 3x in display-linear luminance.
    const w = edgeOn.widthPx;
    const h = edgeOn.heightPx;
    const band = [Math.floor(h * 0.55), Math.floor(h * 0.68)] as const;
    const left = meanLuminance(edgeOn, 0, Math.floor(w * 0.4), band[0], band[1]);
    const right = meanLuminance(edgeOn, Math.ceil(w * 0.6), w, band[0], band[1]);
    expect(right / left).toBeGreaterThan(2);
  });

  it('shows the far side of the disk lensed over the top of the shadow', () => {
    // The signature of a lensed disk (Luminet 1979). Directly above the shadow in the
    // image there is disk light that, without lensing, would be hidden behind the hole.
    const w = edgeOn.widthPx;
    const h = edgeOn.heightPx;
    const cx = Math.floor(w / 2);
    // Find the shadow: the dark run through the centre column.
    let shadowTop = -1;
    for (let y = Math.floor(h / 2); y > 0; y -= 1) {
      const k = (y * w + cx) * 3;
      if (edgeOn.linear[k] + edgeOn.linear[k + 1] + edgeOn.linear[k + 2] > 1e-6) {
        shadowTop = y;
        break;
      }
    }
    expect(shadowTop, 'no lit pixel above the shadow centre').toBeGreaterThan(0);
    expect(meanLuminance(edgeOn, cx - 1, cx + 2, Math.max(0, shadowTop - 3), shadowTop + 1)).toBeGreaterThan(0.01);
  });

  it('is nearly symmetric and entirely redshifted seen face-on', () => {
    // With almost no line-of-sight orbital velocity the shift is gravitational and
    // transverse-Doppler only, g = sqrt(1 - 3M/r) / sqrt(f_obs) < 1 everywhere on the disk.
    const faceOn = renderDisk(2, 4e10, 0.01, 64);
    const [, gMax] = faceOn.diagnostics.frequencyRatioRange;
    expect(gMax).toBeLessThan(1);
    const w = faceOn.widthPx;
    const h = faceOn.heightPx;
    const left = meanLuminance(faceOn, 0, Math.floor(w / 2), 0, h);
    const right = meanLuminance(faceOn, Math.ceil(w / 2), w, 0, h);
    expect(Math.abs(right / left - 1)).toBeLessThan(0.1);
  });
});

describe('display mapping (a visualization choice, disclosed in the UI)', () => {
  it('doubles the display value per stop of exposure', () => {
    const base = displayGain({ referenceLuminance: 5e8, exposureStops: 0, toneMap: 'none' });
    expect(displayGain({ referenceLuminance: 5e8, exposureStops: 1, toneMap: 'none' }) / base).toBeCloseTo(2, 14);
    expect(base * 5e8).toBeCloseTo(1, 14);
  });

  it("compresses luminance with Reinhard's curve while keeping chromaticity", () => {
    const input = { r: 3, g: 1.5, b: 0.6 };
    const out = toneMapReinhard(input);
    expect(out.r / out.g).toBeCloseTo(input.r / input.g, 14);
    expect(out.g / out.b).toBeCloseTo(input.g / input.b, 14);
    const L = 0.2126 * input.r + 0.7152 * input.g + 0.0722 * input.b;
    const Lout = 0.2126 * out.r + 0.7152 * out.g + 0.0722 * out.b;
    expect(Lout).toBeCloseTo(L / (1 + L), 14);
  });

  it('clamps negative, out-of-gamut components to zero', () => {
    const out = toneMapReinhard({ r: -0.2, g: 0.5, b: 0.5 });
    expect(out.r).toBe(0);
  });
});
