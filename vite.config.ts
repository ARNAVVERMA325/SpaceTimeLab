import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    globals: true,
    environment: 'node',
    include: ['tests/**/*.test.ts'],
    // The flat-space propagation gate runs 10^6 integration steps (ROADMAP.md 1.4);
    // it needs more than Vitest's default per-test budget on a cold CI runner.
    testTimeout: 120_000,
  },
});
