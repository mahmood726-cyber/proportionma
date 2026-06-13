// node --test truth-recovery/test-truth-recovery.mjs
// Measured invariants for the proportionma truth-recovery yardstick. Every
// assertion is produced from seeded simulation -- no hand-entered numbers.
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';
import { generate, makeRng, SCENARIOS } from './dgp-prop.mjs';
import { runCell, runGrid, summarize } from './harness.mjs';

describe('proportion DGP', () => {
  it('is reproducible for a fixed seed', () => {
    const a = generate(0.2, 0.25, 10, 'none', makeRng(7));
    const b = generate(0.2, 0.25, 10, 'none', makeRng(7));
    assert.deepEqual(a.studies, b.studies);
  });
  it('returns exactly k studies with 0<=x<=n for every scenario', () => {
    const rng = makeRng(1);
    for (const scen of SCENARIOS) {
      const { studies } = generate(0.2, 0.25, 8, scen, rng);
      assert.equal(studies.length, 8);
      assert.ok(studies.every(s => s.events >= 0 && s.events <= s.total && s.total >= 5));
    }
  });
});

describe('Truth-recovery (measured)', () => {
  it('HKSJ recovers true prevalence at least as well as the z-interval (mean clean coverage)', () => {
    const grid = runGrid({ reps: 200 });
    const s = summarize(grid, c => c.scen === 'none');
    // The app ships HKSJ as an OPT-IN toggle (default off). This is the measured
    // justification for enabling it: it never recovers truth worse than z, and
    // closes the under-coverage gap under heterogeneity.
    assert.ok(s['logit+DL+HKSJ'].meanCoverage >= s['logit+DL+z'].meanCoverage,
      `HKSJ ${s['logit+DL+HKSJ'].meanCoverage} < z ${s['logit+DL+z'].meanCoverage}`);
    assert.ok(s['FT+DL+HKSJ'].meanCoverage >= s['FT+DL+z'].meanCoverage,
      `FT HKSJ ${s['FT+DL+HKSJ'].meanCoverage} < FT z ${s['FT+DL+z'].meanCoverage}`);
  });

  it('the z-interval UNDER-covers under strong heterogeneity; HKSJ restores it toward nominal', () => {
    const rng = makeRng(20260613);
    const cell = runCell(0.20, 1.0, 5, 'none', 400, rng);
    assert.ok(cell['logit+DL+z'].coverage < 0.92,
      `z coverage ${cell['logit+DL+z'].coverage} not under-covering`);
    assert.ok(cell['logit+DL+HKSJ'].coverage > cell['logit+DL+z'].coverage + 0.04,
      `HKSJ ${cell['logit+DL+HKSJ'].coverage} did not improve on z ${cell['logit+DL+z'].coverage}`);
  });

  it('FT is NOT worse than logit on coverage (honest: the real lever is the CI method, not the transform)', () => {
    const grid = runGrid({ reps: 200 });
    const s = summarize(grid, c => c.scen === 'none');
    // Both transforms land within 3pp of each other; logit carries slightly less bias.
    assert.ok(Math.abs(s['FT+DL+HKSJ'].meanCoverage - s['logit+DL+HKSJ'].meanCoverage) < 0.03,
      `FT ${s['FT+DL+HKSJ'].meanCoverage} vs logit ${s['logit+DL+HKSJ'].meanCoverage}`);
    assert.ok(s['logit+DL+HKSJ'].meanAbsBias <= s['FT+DL+HKSJ'].meanAbsBias + 1e-9,
      `logit bias ${s['logit+DL+HKSJ'].meanAbsBias} > FT bias ${s['FT+DL+HKSJ'].meanAbsBias}`);
  });
});
