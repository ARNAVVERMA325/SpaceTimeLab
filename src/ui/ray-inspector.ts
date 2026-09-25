import type { Vec4 } from '../physics/core/indices.js';
import { generateNullRay, localRayDirectionAt } from '../physics/observer/observer.js';
import { frequencyRatio, STATIC_EMITTER } from '../physics/observer/frequency-shift.js';
import { blackbodyXYZ, xyzToLinearSrgb } from '../physics/radiation/blackbody.js';
import { asymptoticSweepTail } from '../physics/spacetimes/schwarzschild-analytic.js';
import { criticalImpactParameter } from '../physics/spacetimes/schwarzschild.js';
import type { PhaseSpaceState } from '../physics/core/phase-space.js';
import type { CartesianVec3, SpacetimeModel } from '../physics/spacetimes/spacetime-model.js';
import { checkTolerance, NULL_NORMALIZATION_TRACED } from '../physics/validation/tolerances.js';
import { linearToSrgb8 } from '../visualization/color.js';
import { liftDirection, liftPosition, reduceToOrbitalPlane } from '../visualization/orbital-plane.js';
import { traceRay, type RayOutcome } from '../visualization/raytracer.js';
import type { BuiltScene } from '../visualization/scene.js';

/**
 * One pixel, traced on its own and reported in full.
 *
 * The image says what the sky looks like; this says why a particular line of sight looks
 * that way. It re-traces the ray through the pixel centre with the scene's own
 * configuration — same metric, same integrator, same tolerance — so what it reports is the
 * ray behind that pixel and not a second, more convenient calculation.
 *
 * Every quantity here is either read off the traced ray or computed from it in closed
 * form. Where a number depends on a frame, the row says which frame (CLAUDE.md §1.3).
 */

export type Tone = 'neutral' | 'ok' | 'warn' | 'fail';

export interface InspectorRow {
  readonly label: string;
  readonly value: string;
  readonly note?: string;
  readonly tone?: Tone;
  /** A colour swatch to show beside the value, as a CSS colour. */
  readonly swatch?: string;
}

export interface RayInspection {
  readonly pixel: readonly [number, number];
  readonly outcome: RayOutcome;
  readonly headline: string;
  readonly rows: readonly InspectorRow[];
}

const DEG = 180 / Math.PI;

function angleBetween(a: CartesianVec3, b: CartesianVec3): number {
  const dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross: CartesianVec3 = [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
  // atan2 of |a x b| against a . b, which stays accurate at small and large angles alike.
  return Math.atan2(Math.hypot(cross[0], cross[1], cross[2]), dot);
}

/**
 * The colour of a blackbody at T, normalized to full brightness.
 *
 * Chromaticity only: the swatch shows the hue of the observed spectrum, not how bright it
 * is. Brightness in the image comes from the exposure and tone curve, which this deliberately
 * drops so that a dim, strongly redshifted ring still shows its colour.
 */
export function temperatureSwatch(temperatureK: number): string {
  if (!(temperatureK > 0)) return '#000000';
  const rgb = xyzToLinearSrgb(blackbodyXYZ(temperatureK));
  const peak = Math.max(rgb.r, rgb.g, rgb.b, Number.MIN_VALUE);
  const r = linearToSrgb8(rgb.r / peak);
  const g = linearToSrgb8(rgb.g / peak);
  const b = linearToSrgb8(rgb.b / peak);
  return `rgb(${r}, ${g}, ${b})`;
}

/**
 * Impact parameter b = L / E of a ray, from its initial phase-space state.
 *
 * A property of the geodesic rather than of whoever launched it: E = -p_t and the total
 * angular momentum L = r^2 sqrt((k^theta)^2 + sin^2(theta) (k^phi)^2) are both conserved in
 * a static, spherically symmetric spacetime, so two observers in relative motion at the
 * same event disagree about the viewing angle but agree about b.
 *
 * Magnitudes, because the renderer traces backward in time: its rays are past-directed and
 * carry E < 0, while the ray's geometry is the same either way.
 */
export function impactParameter(model: SpacetimeModel, position_x: Vec4, tangent: Vec4): number {
  const r = position_x[1];
  const sinTheta = Math.sin(position_x[2]);
  const metric = model.metricAt(position_x);
  // E = -p_t, general for a stationary metric; L is the total angular momentum, which is
  // conserved because the spacetime is also spherically symmetric.
  let energy_E = 0;
  for (let mu = 0; mu < 4; mu += 1) energy_E -= metric.g(0, mu as 0 | 1 | 2 | 3) * tangent[mu];
  const angularMomentum = r * r * Math.hypot(tangent[2], sinTheta * tangent[3]);
  return angularMomentum / Math.abs(energy_E);
}

function residualRow(nullResidual: number): InspectorRow {
  const check = checkTolerance(NULL_NORMALIZATION_TRACED, nullResidual);
  return {
    label: 'Null normalization',
    value: `g_mu_nu k^mu k^nu = ${nullResidual.toExponential(3)}`,
    tone: check.withinTolerance ? 'ok' : 'fail',
    note:
      `Exactly 0 for a null geodesic; judged against ${NULL_NORMALIZATION_TRACED.value.toExponential(0)} ` +
      `(${NULL_NORMALIZATION_TRACED.id}). This ray: ${check.withinTolerance ? 'within tolerance' : 'OUTSIDE tolerance'}.`,
  };
}

/**
 * The angle the ray sweeps about the centre, camera to infinity, in degrees.
 *
 * Not the angle between the launch and exit directions: that saturates at 180 degrees, so
 * a ray looping the photon sphere twice would report less turn than one looping once. The
 * sweep is read from the ray's own orbital plane, where phi advances monotonically, and
 * the sweep past the integration radius is added from the same tail integral the
 * background lookup uses.
 *
 * Undefined when the scene does not reduce to an orbital plane, as in flat spacetime,
 * where the caller falls back to the angle between the two directions.
 */
function totalSweepDeg(
  scene: BuiltScene,
  initial: PhaseSpaceState,
  ray: { readonly frame?: unknown; readonly final: PhaseSpaceState },
): number | undefined {
  if (!scene.config.orbitalPlaneReduction || !ray.frame) return undefined;
  const reduced = reduceToOrbitalPlane(scene.model, initial);
  const swept = Math.abs(ray.final.position_x[3] - reduced.initial.position_x[3]);
  const b = impactParameter(scene.model, initial.position_x, initial.tangent);
  return (swept + asymptoticSweepTail(1, b, ray.final.position_x[1])) * DEG;
}

export function inspectPixel(scene: BuiltScene, i: number, j: number): RayInspection {
  const { model, observer, config, screen } = scene;
  const localDirection = localRayDirectionAt(screen, i + 0.5, j + 0.5);
  const initial = generateNullRay(observer, localDirection);
  const ray = traceRay(config, initial);

  const launchDirection = model.geometry.toCartesianDirection(initial.position_x, initial.tangent);
  const offAxisDeg = angleBetween(localDirection, screen.forward) * DEG;

  const rows: InspectorRow[] = [
    {
      label: 'Line of sight',
      value: `pixel (${i}, ${j}), ${offAxisDeg.toFixed(2)} degrees off the image centre`,
      note:
        'The direction is set in the observer’s own orthonormal frame, so it is an angle ' +
        'this observer actually measures rather than a coordinate difference.',
    },
  ];

  const curved = scene.description.kind !== 'minkowski';
  if (curved) {
    const b = impactParameter(model, initial.position_x, initial.tangent);
    const bc = criticalImpactParameter(1);
    const margin = (b / bc - 1) * 100;
    rows.push({
      label: 'Impact parameter',
      value: `b = ${b.toFixed(4)} M  (b_c = 3 sqrt(3) M = ${bc.toFixed(4)} M)`,
      tone: b < bc ? 'warn' : 'neutral',
      note:
        `${Math.abs(margin).toFixed(2)}% ${margin < 0 ? 'below' : 'above'} the capture threshold. ` +
        'b = L / E is a property of the geodesic, not of the camera: a hovering and a falling ' +
        'observer at this event point their cameras differently to launch this same ray, and ' +
        'both get this b. Below b_c the ray must fall in; above it, it must escape.',
    });
  }

  let headline: string;
  switch (ray.outcome) {
    case 'captured':
      headline =
        'This line of sight ends on the black hole. No background light reaches the camera ' +
        'along it, which is what makes the shadow dark.';
      break;
    case 'disk': {
      const sample = ray.disk;
      headline = sample
        ? `This line of sight ends on the disk at r = ${(sample.radius).toFixed(2)} M, seen at g = ${sample.frequencyRatio.toFixed(4)}.`
        : 'This line of sight ends on the disk.';
      if (sample) {
        const position = ray.frame
          ? liftPosition(model, ray.frame, ray.final)
          : model.geometry.toCartesianPosition(ray.final.position_x);
        const direction = ray.frame
          ? liftDirection(model, ray.frame, ray.final)
          : model.geometry.toCartesianDirection(ray.final.position_x, ray.final.tangent);
        const f = 1 - 2 / sample.radius;
        const staticRatio = frequencyRatio(
          { energy: Math.sqrt(f) * ray.final.tangent[0], momentum: direction },
          STATIC_EMITTER,
        );
        const orbitalSpeed = Math.sqrt(1 / sample.radius) / Math.sqrt(f);
        rows.push(
          {
            label: 'Emission point',
            value: `r = ${sample.radius.toFixed(3)} M, |z| = ${Math.abs(position[2]).toExponential(1)} M off the plane`,
            note:
              'The disk is geometrically thin, so emission comes from the equatorial plane itself. ' +
              'That |z| is how precisely the crossing was located, not a thickness: the event is ' +
              'found by bisection on the dense output and then refined, rather than at whichever ' +
              'step happened to straddle the plane.',
          },
          {
            label: 'Gas velocity there',
            value: `${orbitalSpeed.toFixed(4)} c, prograde, measured by a static observer at that point`,
            note:
              'A circular geodesic: Omega = sqrt(M / r^3) in these coordinates, which a local ' +
              'static observer sees as v = r Omega / sqrt(f). At the ISCO it is exactly 1/2.',
          },
          {
            label: 'Frequency shift',
            value: `g = ${sample.frequencyRatio.toFixed(5)}  (${sample.frequencyRatio >= 1 ? 'blueshift' : 'redshift'})`,
            tone: 'neutral',
            note:
              'One contraction, (k . u_obs) / (k . u_emit), carrying gravitational redshift, ' +
              'transverse Doppler and line-of-sight Doppler together.',
          },
          {
            label: 'Split into two factors',
            value: `${staticRatio.toFixed(5)} without the orbital motion, x ${(sample.frequencyRatio / staticRatio).toFixed(5)} from it`,
            note:
              'A decomposition relative to the static frame at the emission event, not an ' +
              'invariant one: the first factor is what a static emitter there would show, the ' +
              'second is the Doppler factor of the orbiting gas against it.',
          },
          {
            label: 'Temperature',
            value: `${Math.round(sample.emittedTemperatureK).toLocaleString('en-US')} K emitted, seen as ${Math.round(sample.observedTemperatureK).toLocaleString('en-US')} K`,
            swatch: temperatureSwatch(sample.observedTemperatureK),
            note:
              'I_nu / nu^3 is invariant, so a blackbody at T seen with ratio g is exactly a ' +
              'blackbody at g T. The swatch is the hue of that spectrum at full brightness; how ' +
              'bright the pixel is depends on the exposure and tone curve.',
          },
        );
      }
      break;
    }
    case 'background': {
      const turned = totalSweepDeg(scene, initial, ray) ?? angleBetween(launchDirection, ray.exitDirection) * DEG;
      const loops = turned / 360;
      headline =
        `This line of sight reaches the background sphere. The ray swept ${turned.toFixed(2)} degrees ` +
        (loops >= 1
          ? `around the hole \u2014 ${loops.toFixed(2)} full loops \u2014 before escaping, so this pixel shows a part of the sky the camera is nowhere near pointing at.`
          : 'between the camera and infinity, so it shows a part of the sky that far from where the camera is pointing.');
      rows.push({
        label: 'Total sweep',
        value: `${turned.toFixed(4)} degrees${loops >= 1 ? `  (${loops.toFixed(2)} loops)` : ''}`,
        note:
          'The angle the ray sweeps about the hole in its own orbital plane, from the camera out ' +
          'to infinity. Exactly 0 in flat spacetime, and unbounded as b approaches b_c: a ray can ' +
          'circle the photon sphere any number of times, which is what stacks the Einstein rings ' +
          'against the shadow edge. ' +
          (ray.asymptoticallyCorrected
            ? 'The sweep beyond the integration radius is added from the orbit-equation tail integral, so it does not depend on where the integration stopped.'
            : 'Read off where the integration stopped, with no tail correction.'),
      });
      break;
    }
    default:
      headline =
        `This ray did not finish: ${ray.outcome} (${ray.terminationReason}). It is drawn in magenta ` +
        'rather than being blended into the picture.';
  }

  rows.push(
    {
      label: 'Integration',
      value: `${ray.steps.toLocaleString('en-US')} steps, ended by ${ray.terminationReason}`,
      note: `${config.integrator.displayName}, ${(config.formulation ?? { displayName: 'Hamiltonian' }).displayName}.`,
    },
    residualRow(Math.abs(ray.nullResidual)),
  );

  return { pixel: [i, j], outcome: ray.outcome, headline, rows };
}
