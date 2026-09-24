import type { Vec4 } from '../physics/core/indices.js';
import { HAMILTONIAN } from '../physics/geodesic/formulation.js';
import { Dopri5Integrator } from '../physics/geodesic/integrators/dopri5.js';
import { STATE_DIM } from '../physics/geodesic/state-vector.js';
import {
  defaultScreen,
  freeFallingObserver,
  inwardFacingScreen,
  staticMinkowskiObserver,
  staticObserver,
  type Observer,
  type PinholeScreen,
} from '../physics/observer/observer.js';
import { minkowski } from '../physics/spacetimes/minkowski.js';
import { NOVIKOV_THORNE_EFFICIENCY, thinDiskScale } from '../physics/spacetimes/novikov-thorne.js';
import { criticalImpactParameter, photonSphereRadius, schwarzschild } from '../physics/spacetimes/schwarzschild.js';
import { asymptoticDirection, isCaptured } from '../physics/spacetimes/schwarzschild-rays.js';
import type { SpacetimeModel } from '../physics/spacetimes/spacetime-model.js';
import type { ProvenanceEntry } from '../ui/provenance.js';
import { DEFAULT_CELESTIAL_GRID } from './celestial-grid.js';
import { diskReferenceLuminance, type ThinDiskEmitter } from './disk-emission.js';
import type { TraceConfig } from './raytracer.js';

/**
 * Scenes as plain data.
 *
 * A render worker cannot be handed functions, a model or an integrator, only
 * structured-cloneable data. So a scene is described by the parameters a person chooses,
 * and `buildScene` turns that description into a live configuration — identically on the
 * main thread and in every worker, which is what makes a parallel render reproduce a
 * serial one bit for bit.
 */

export type ObserverChoice = 'static' | 'free-fall';

interface SceneCommon {
  readonly widthPx: number;
  readonly heightPx: number;
  readonly samplesPerAxis: number;
  readonly seed: number;
}

export type SceneDescription =
  | (SceneCommon & { readonly kind: 'minkowski' })
  | (SceneCommon & {
      readonly kind: 'schwarzschild-sky';
      readonly cameraRadius: number;
      readonly observer: ObserverChoice;
    })
  | (SceneCommon & {
      readonly kind: 'schwarzschild-disk';
      readonly cameraRadius: number;
      readonly inclinationDeg: number;
      readonly observer: ObserverChoice;
      readonly massSolar: number;
      readonly eddingtonFraction: number;
      readonly exposureStops: number;
      readonly skyGrid: boolean;
    });

export interface BuiltScene {
  readonly description: SceneDescription;
  readonly model: SpacetimeModel;
  readonly observer: Observer;
  readonly config: TraceConfig;
  readonly screen: PinholeScreen;
  readonly caption: string;
  readonly entries: readonly ProvenanceEntry[];
}

/**
 * Integration tolerance for every rendered scene. With the Hamiltonian formulation and
 * Dormand-Prince it holds the worst null residual over a whole image at 5.1e-10 —
 * inside the same 1e-9 gate the test suite uses.
 */
export const RENDER_TOLERANCE = 1e-10;

const MASS = 1;

function integrator(): Dopri5Integrator {
  return new Dopri5Integrator(STATE_DIM, {
    tolerance: { absolute: RENDER_TOLERANCE, relative: RENDER_TOLERANCE },
  });
}

function schwarzschildObserver(choice: ObserverChoice, position: Vec4): Observer {
  const model = schwarzschild(MASS);
  return choice === 'free-fall' ? freeFallingObserver(model, position, MASS) : staticObserver(model, position);
}

function observerEntry(choice: ObserverChoice, r: number): ProvenanceEntry {
  if (choice === 'free-fall') {
    const v = Math.sqrt((2 * MASS) / r);
    return {
      label: 'Camera motion',
      value: `Falling radially from rest at infinity, at ${v.toFixed(3)}c relative to a static observer`,
      note:
        'A geodesic observer. Aberration crowds the view ahead toward the direction of fall, ' +
        'so the shadow appears smaller than a hovering observer at the same event sees it, ' +
        'and the Doppler shift brightens and blueshifts what lies ahead.',
    };
  }
  return {
    label: 'Camera motion',
    value: `Static: hovering at fixed r = ${r}M`,
    note:
      'Holding station requires proper acceleration; this is not a freely-falling frame. ' +
      'A freely-falling observer passing the same event would see a different image.',
  };
}

export function buildScene(description: SceneDescription): BuiltScene {
  const { widthPx, heightPx, samplesPerAxis, seed } = description;
  const sampling = { samplesPerAxis, seed };

  if (description.kind === 'minkowski') {
    const observer = staticMinkowskiObserver([0, 0, 0, 0]);
    return {
      description,
      model: minkowski,
      observer,
      screen: defaultScreen(widthPx, heightPx, Math.PI / 2),
      config: {
        model: minkowski,
        integrator: integrator(),
        formulation: HAMILTONIAN,
        observer,
        grid: { ...DEFAULT_CELESTIAL_GRID, radius: 100 },
        limits: { initialStep: 0.5, maxStep: 10, parameterMax: 1000, maxSteps: 10_000 },
        sampling,
      },
      caption:
        'Flat spacetime: nothing is bent, so the sky grid arrives exactly as a pinhole camera ' +
        'projects it. Curved grid lines here are rectilinear projection of a sphere, not light deflection.',
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
    };
  }

  const model = schwarzschild(MASS);
  const common = {
    model,
    integrator: integrator(),
    formulation: HAMILTONIAN,
    captureTest: (x: Vec4, k: Vec4) => isCaptured(model, x, k),
    orbitalPlaneReduction: true,
    asymptoticDirection: (p: readonly [number, number, number], d: readonly [number, number, number]) =>
      asymptoticDirection(model, p, d),
    sampling,
  };
  const sharedEntries: ProvenanceEntry[] = [
    {
      label: 'Mass parameter',
      value: `M = ${MASS} (geometric units)`,
      note: `Horizon at r = ${2 * MASS}M, photon sphere at r = ${photonSphereRadius(MASS)}M.`,
    },
    {
      label: 'Ray termination',
      value: 'Captured when r < 3M with k^r < 0',
      note:
        'Exact rather than a tuned cutoff: the null effective potential f/r^2 increases ' +
        'inward of r = 3M, so a photon moving inward there can never turn around. Rays stop ' +
        'while the chart is still well behaved rather than being integrated toward r = 2M, ' +
        'where these coordinates break down.',
    },
    {
      label: 'Background direction',
      value: 'Exact direction at infinity',
      note:
        'A ray is still being bent where the integration stops, so reading the sky off its ' +
        'local direction would bias the image by the deflection still to come. The remaining ' +
        'sweep to infinity is added exactly from the orbit-equation tail integral.',
    },
    {
      label: 'Sampling',
      value: `${samplesPerAxis * samplesPerAxis} stratified sample(s) per pixel, averaged in linear light`,
      note:
        'Every sample is an independent, fully integrated geodesic; accumulation only averages ' +
        'their results, so no trajectory is altered (CLAUDE.md §9). The finest photon-ring ' +
        'bands still alias at any sample count: each is e^pi ~ 23 times thinner than the last.',
    },
  ];

  if (description.kind === 'schwarzschild-sky') {
    const position: Vec4 = [0, description.cameraRadius, Math.PI / 2, 0];
    const observer = schwarzschildObserver(description.observer, position);
    const f = 1 - (2 * MASS) / description.cameraRadius;
    const shadowStatic = Math.asin((criticalImpactParameter(MASS) * Math.sqrt(f)) / description.cameraRadius);
    return {
      description,
      model,
      observer,
      screen: inwardFacingScreen(widthPx, heightPx, 3.2 * shadowStatic),
      config: {
        ...common,
        observer,
        grid: { ...DEFAULT_CELESTIAL_GRID, radius: 200 },
        limits: { initialStep: 1e-3, parameterMax: 20_000, maxSteps: 200_000, maxStep: 5 },
      },
      caption:
        'Computed appearance of the background grid for the selected spacetime, observer and ' +
        'rendering assumptions. The dark disc is the black-hole shadow: lines of sight along ' +
        'which no background light reaches the observer.',
      entries: [
        ...sharedEntries,
        observerEntry(description.observer, description.cameraRadius),
        {
          label: 'Shadow (static observer)',
          value: `Angular radius ${((shadowStatic * 180) / Math.PI).toFixed(2)} degrees`,
          note:
            'sin(psi) = b_c sqrt(f) / r with b_c = 3 sqrt(3) M, matched by traced rays to better ' +
            'than a part in a million. The shadow is larger than the horizon and is not a ' +
            'picture of it.',
        },
        {
          label: 'Background',
          value: 'A latitude/longitude grid on the celestial sphere',
          note: 'A visualization texture, not an emitting surface. It has no spectrum and is drawn unshifted.',
        },
      ],
    };
  }

  // Thin disk.
  const inclination = (Math.max(2, Math.min(89, description.inclinationDeg)) * Math.PI) / 180;
  const position: Vec4 = [0, description.cameraRadius, inclination, 0];
  const observer = schwarzschildObserver(description.observer, position);
  const scale = thinDiskScale({
    massSolar: description.massSolar,
    eddingtonFraction: description.eddingtonFraction,
  });
  const disk: ThinDiskEmitter = { mass: MASS, innerRadius: 6 * MASS, outerRadius: 20 * MASS, scale };
  const referenceLuminance = diskReferenceLuminance(disk);
  const fov = 2 * Math.atan((26 * MASS) / description.cameraRadius);

  return {
    description,
    model,
    observer,
    screen: inwardFacingScreen(widthPx, heightPx, fov * (widthPx / heightPx > 1.4 ? 1 : 1.2)),
    config: {
      ...common,
      observer,
      grid: { ...DEFAULT_CELESTIAL_GRID, radius: 200 },
      limits: { initialStep: 1e-3, parameterMax: 20_000, maxSteps: 200_000, maxStep: 5 },
      disk,
      display: { referenceLuminance, exposureStops: description.exposureStops, toneMap: 'reinhard' },
      background: description.skyGrid ? 'grid' : 'black',
    },
    caption:
      'Computed appearance of a thin accretion disk for the selected spacetime, observer, ' +
      'emission model and display mapping. Not a prediction of any observed image: real ' +
      'accretion flows add plasma, scattering and radiative transfer this model omits.',
    entries: [
      ...sharedEntries,
      observerEntry(description.observer, description.cameraRadius),
      {
        label: 'Viewing inclination',
        value: `${((inclination * 180) / Math.PI).toFixed(1)} degrees from the disk axis, at r = ${description.cameraRadius}M`,
      },
      {
        label: 'Emission model',
        value: 'Novikov-Thorne thin disk, 6M to 20M, local blackbody',
        note:
          'Geometrically thin, optically thick disk in the equatorial plane; gas on prograde ' +
          'Keplerian circular geodesics; zero torque at the ISCO; each face radiating as a ' +
          'blackbody at sigma T^4 = F with the Page-Thorne flux. Verified by its own energy ' +
          `balance: the radiated luminosity equals Mdot (1 - E_isco) = ${(NOVIKOV_THORNE_EFFICIENCY * 100).toFixed(3)}% of Mdot c^2.`,
      },
      {
        label: 'Physical scale',
        value:
          `M = ${description.massSolar.toExponential(1)} solar masses, accreting at ` +
          `${description.eddingtonFraction} of Eddington (${scale.accretionRateSolarPerYear.toPrecision(3)} solar masses per year)`,
        note:
          `Peak rest-frame temperature ${Math.round(scale.peakTemperatureK).toLocaleString('en-US')} K at ` +
          'r = 9.55M, from CODATA 2018 constants. Nothing here is tuned for appearance: a thin ' +
          'disk around any real black hole is hot, and the colour shown is what that temperature gives.',
      },
      {
        label: 'Frequency shift',
        value: 'g = nu_obs / nu_emit, from the traced wavevector',
        note:
          'Gravitational redshift, transverse Doppler and line-of-sight Doppler all come from one ' +
          'contraction, (k . u_obs) / (k . u_emit). The observed spectrum of a blackbody at T is a ' +
          'blackbody at g T (I_nu / nu^3 is invariant), so brightening and colour change are one effect.',
      },
      {
        label: 'Colour',
        value: 'CIE 1931 colorimetry of the observed blackbody, in sRGB',
        note:
          'Planck spectrum at g T integrated against the CIE 1931 2-degree colour-matching ' +
          'functions (CIE 015:2018), converted to sRGB. Out-of-gamut components are clamped to zero.',
      },
      {
        label: 'Display mapping',
        value: `Exposure ${description.exposureStops >= 0 ? '+' : ''}${description.exposureStops} stops, Reinhard tone curve`,
        note:
          'At 0 stops the hottest ring, seen at rest, has display luminance 1, which the tone ' +
          "curve L' = L / (1 + L) puts at 0.5. The curve compresses a range of radiance far too " +
          'wide for a screen while keeping chromaticity, so a blueshifted ring stays blue as it ' +
          'saturates. This changes how the picture looks, never what was computed: the linear ' +
          'radiance behind every pixel is kept unmodified.',
      },
      {
        label: 'Why the contrast depends on temperature',
        value: 'The visible band sits in a different part of the Planck curve for each disk',
        note:
          'The observed spectrum is a blackbody at g T, so the bolometric boost is exactly g^4 ' +
          'whatever the temperature. A colour image shows in-band luminance instead, and that ' +
          'depends on where the visible band falls on the Planck curve. Below about 10,000 K the ' +
          'band sits near the peak and the boost is steeper than g^4; far above it the band lies ' +
          'in the Rayleigh-Jeans tail, where B_nu is proportional to T and the boost falls to g. ' +
          'Measured on this renderer at 80 degrees, the approaching side outshines the receding ' +
          'one by 4.3x at 8,830 K but only 1.8x at 39,500 K, with identical kinematics. A hot ' +
          'disk looks almost uniformly blue-white however fast its gas moves.',
      },
      {
        label: 'Omitted physics',
        value: 'Disk atmosphere, electron scattering, limb darkening, returning radiation, plasma, jets',
        note:
          'No spectral hardening by electron scattering, no limb darkening, no self-irradiation or ' +
          'returning radiation, no emission from inside the ISCO, no optically thin or thick ' +
          'coronal flow. The disk does not affect the spacetime. Thin disks are a reasonable model ' +
          'only at moderate accretion rates; low-rate flows are thought to be hot and geometrically thick.',
      },
      {
        label: 'Sky',
        value: description.skyGrid ? 'Coordinate grid on the celestial sphere (visualization only)' : 'Empty',
        note: description.skyGrid
          ? 'A visualization texture with no spectrum, drawn unshifted. It is not physical starlight.'
          : 'No background sources are modelled.',
      },
    ],
  };
}
