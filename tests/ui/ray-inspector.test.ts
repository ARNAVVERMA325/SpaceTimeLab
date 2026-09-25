import { describe, expect, it } from 'vitest';
import { generateNullRay, localRayDirectionAt } from '../../src/physics/observer/observer.js';
import { criticalImpactParameter, schwarzschild } from '../../src/physics/spacetimes/schwarzschild.js';
import { impactParameter, inspectPixel, temperatureSwatch } from '../../src/ui/ray-inspector.js';
import { buildScene, type SceneDescription } from '../../src/visualization/scene.js';

/**
 * The pixel inspector (CLAUDE.md §22: the viewer must be able to inspect what produced
 * what they are looking at).
 *
 * The inspector's job is to report the ray behind a pixel, so these tests check that what
 * it reports is the same physics the renderer used and the same physics the closed forms
 * predict — not that it produces some plausible-looking text.
 */

const M = 1;
const model = schwarzschild(M);
const b_c = criticalImpactParameter(M);

function sky(widthPx: number, heightPx: number): SceneDescription {
  return { kind: 'schwarzschild-sky', widthPx, heightPx, samplesPerAxis: 1, seed: 1, cameraRadius: 20, observer: 'static' };
}

function rowValue(rows: readonly { label: string; value: string }[], label: string): string {
  const row = rows.find((r) => r.label === label);
  if (!row) throw new Error(`no row labelled "${label}" in: ${rows.map((r) => r.label).join(', ')}`);
  return row.value;
}

describe('impact parameter from a ray', () => {
  it('agrees with r sin(psi) / sqrt(f) for a static observer', () => {
    // The closed form the shadow prediction uses, from the angle psi between the ray and
    // the inward radial direction in the observer's frame. The inspector instead builds
    // b = L / E from the conserved quantities, so the two routes share no arithmetic.
    const scene = buildScene(sky(41, 31));
    const r = 20;
    const f = model.lapseFunction(r);
    for (const [i, j] of [[0, 0], [20, 15], [40, 30], [10, 22], [33, 4]] as const) {
      const direction = localRayDirectionAt(scene.screen, i + 0.5, j + 0.5);
      const ray = generateNullRay(scene.observer, direction);
      const psi = Math.acos(Math.min(1, Math.max(-1, -direction[0])));
      const closedForm = (r * Math.sin(psi)) / Math.sqrt(f);
      const measured = impactParameter(model, ray.position_x, ray.tangent);
      // Absolute at the image centre, where b is 0 and a ratio would be meaningless.
      expect(Math.abs(measured - closedForm), `pixel (${i}, ${j})`).toBeLessThan(1e-12 * Math.max(1, closedForm));
    }
  });

  it('decides capture exactly at b = 3 sqrt(3) M, on every pixel it is asked about', () => {
    // b < b_c must fall in and b > b_c must escape: a closed-form threshold against a
    // traced outcome, checked across the image rather than at one chosen ray.
    const scene = buildScene(sky(37, 27));
    let captured = 0;
    let escaped = 0;
    for (let j = 0; j < 27; j += 3) {
      for (let i = 0; i < 37; i += 3) {
        const direction = localRayDirectionAt(scene.screen, i + 0.5, j + 0.5);
        const ray = generateNullRay(scene.observer, direction);
        const b = impactParameter(model, ray.position_x, ray.tangent);
        // Rays within a part in a million of the threshold are genuinely ill-conditioned;
        // the renderer's own tolerance cannot resolve which side they fall on.
        if (Math.abs(b / b_c - 1) < 1e-6) continue;
        const outcome = inspectPixel(scene, i, j).outcome;
        if (b < b_c) {
          expect(outcome, `pixel (${i}, ${j}), b = ${b}`).toBe('captured');
          captured += 1;
        } else {
          expect(outcome, `pixel (${i}, ${j}), b = ${b}`).toBe('background');
          escaped += 1;
        }
      }
    }
    expect(captured).toBeGreaterThan(10);
    expect(escaped).toBeGreaterThan(10);
  });
});

describe('what the inspector reports', () => {
  it('reports no turn at all in flat spacetime', () => {
    const scene = buildScene({ kind: 'minkowski', widthPx: 21, heightPx: 15, samplesPerAxis: 1, seed: 1 });
    for (const [i, j] of [[0, 0], [10, 7], [20, 14]] as const) {
      const inspection = inspectPixel(scene, i, j);
      expect(inspection.outcome).toBe('background');
      const turnDeg = Number.parseFloat(rowValue(inspection.rows, 'Total sweep'));
      expect(Math.abs(turnDeg), `pixel (${i}, ${j})`).toBeLessThan(1e-6);
    }
    // And no impact-parameter row: there is no centre to have an impact parameter about.
    expect(inspectPixel(scene, 10, 7).rows.some((r) => r.label === 'Impact parameter')).toBe(false);
  });

  it('turns a ray that grazes the photon sphere much further than a distant one', () => {
    // Strong deflection is the whole content of the lensed image: alpha diverges
    // logarithmically as b approaches b_c, so the pixel just outside the shadow must turn
    // far more than one at the edge of the frame. Columns are chosen by their own b rather
    // than by a guessed pixel index.
    const scene = buildScene(sky(101, 75));
    const middleRow = 37;
    const turns = new Map<number, { b: number; turn: number }>();
    for (let i = 51; i < 101; i += 1) {
      const direction = localRayDirectionAt(scene.screen, i + 0.5, middleRow + 0.5);
      const ray = generateNullRay(scene.observer, direction);
      const b = impactParameter(model, ray.position_x, ray.tangent);
      if (b < b_c * 1.0005) continue;
      const inspection = inspectPixel(scene, i, middleRow);
      if (inspection.outcome !== 'background') continue;
      turns.set(i, { b, turn: Number.parseFloat(rowValue(inspection.rows, 'Total sweep')) });
    }
    const byB = [...turns.values()].sort((a, b) => a.b - b.b);
    expect(byB.length).toBeGreaterThan(5);
    const grazing = byB[0];
    const distant = byB[byB.length - 1];
    expect(grazing.turn).toBeGreaterThan(distant.turn);
    expect(distant.turn).toBeGreaterThan(0);
    // Monotone in b: every step outward from the hole turns the ray less than the last.
    for (let n = 1; n < byB.length; n += 1) {
      expect(byB[n].turn, `b = ${byB[n].b}`).toBeLessThan(byB[n - 1].turn);
    }
  });

  it('splits a disk ray’s shift into two factors whose product is g', () => {
    const scene = buildScene({
      kind: 'schwarzschild-disk',
      widthPx: 61,
      heightPx: 45,
      samplesPerAxis: 1,
      seed: 1,
      cameraRadius: 30,
      inclinationDeg: 80,
      observer: 'static',
      massSolar: 4e10,
      eddingtonFraction: 0.01,
      exposureStops: 0,
      skyGrid: false,
    });

    let checked = 0;
    for (let i = 0; i < 61; i += 5) {
      const inspection = inspectPixel(scene, i, 30);
      if (inspection.outcome !== 'disk') continue;
      checked += 1;
      const g = Number.parseFloat(rowValue(inspection.rows, 'Frequency shift').replace('g = ', ''));
      const [staticPart, dopplerPart] = rowValue(inspection.rows, 'Split into two factors')
        .split('x')
        .map((part) => Number.parseFloat(part.trim()));
      expect(Math.abs(staticPart * dopplerPart - g), `pixel ${i}`).toBeLessThan(1e-4);

      const temperatures = rowValue(inspection.rows, 'Temperature');
      expect(temperatures).toMatch(/K emitted, seen as .* K/);
    }
    expect(checked).toBeGreaterThanOrEqual(3);
  });
});

describe('temperature swatches', () => {
  it('runs red at low temperature and blue at high, and is always a valid colour', () => {
    const parse = (css: string): [number, number, number] => {
      const match = /^rgb\((\d+), (\d+), (\d+)\)$/.exec(css);
      if (!match) throw new Error(`not a colour: ${css}`);
      return [Number(match[1]), Number(match[2]), Number(match[3])];
    };
    const [rCool, , bCool] = parse(temperatureSwatch(2000));
    const [rHot, , bHot] = parse(temperatureSwatch(30000));
    expect(rCool).toBeGreaterThan(bCool);
    expect(bHot).toBeGreaterThan(rHot);
    // Normalized to full brightness, so the brightest channel is always at the top.
    expect(Math.max(rCool, bCool)).toBe(255);
    expect(temperatureSwatch(0)).toBe('#000000');
  });
});
