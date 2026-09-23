import type { SpacetimeModel } from '../../src/physics/spacetimes/spacetime-model.js';

/**
 * A test-only model that behaves like `base` except where `overrides` says otherwise.
 *
 * Built with a Proxy so that it keeps working as the SpacetimeModel interface grows:
 * hand-written model literals had to be updated every time a method was added, and a
 * forgotten method would surface as `undefined is not a function` far from the cause.
 * Methods are bound to `base`, so delegated calls see the real model's state.
 */
export function overrideModel(
  base: SpacetimeModel,
  overrides: Partial<SpacetimeModel> & { readonly id: string },
): SpacetimeModel {
  return new Proxy(base, {
    get(target, property) {
      if (Object.prototype.hasOwnProperty.call(overrides, property)) {
        return (overrides as Record<PropertyKey, unknown>)[property];
      }
      const value = Reflect.get(target, property, target);
      return typeof value === 'function' ? value.bind(target) : value;
    },
  });
}
