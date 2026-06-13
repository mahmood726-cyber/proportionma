// ============================================================
// harness.mjs -- Truth-recovery yardstick for proportionma.
//
// Wires the app's OWN extracted engine (engine.mjs, copied verbatim from
// index.html) to the known-truth proportion DGP and measures how often each
// transform x CI-method interval covers the TRUE prevalence p0.
//
// Truth-first: every number printed comes from seeded simulation here.
// Run:  node truth-recovery/harness.mjs --reps 400
// ============================================================

import { transformStudies, poolDL, poolREML, backTransform, harmonicMean, qt, qnorm }
  from './engine.mjs';
import { generate, makeRng, SCENARIOS } from './dgp-prop.mjs';

const BASE_SEED = 20260613;

// Reproduces the EXACT CI path in index.html runAnalysis(): transform -> pool
// -> (HKSJ t-interval with max(1,.) floor | z-interval) -> back-transform.
// `piDfMode` lets us also test the t_{k-1} (Cochrane v6.5) PI variant vs the
// shipped t_{k-2}, but PI is reported separately; coverage here is for the CI.
function analyze(studies, { transform, poolMethod = 'DL', useHKSJ = false }) {
  const tr = transformStudies(studies, transform);
  const result = poolMethod === 'REML' ? poolREML(tr) : poolDL(tr);
  if (!result) return null;
  const alpha = 0.05;
  let mult, effSE;
  if (useHKSJ && result.k > 1) {
    const tVal = Math.abs(qt(alpha / 2, result.k - 1));
    const wRE = tr.map(s => 1 / (s.vi + result.tau2));
    const sumWRE = wRE.reduce((a, b) => a + b, 0);
    const qRE = wRE.reduce((s, wi, i) => s + wi * (tr[i].yi - result.mu) ** 2, 0);
    const adjFactor = Math.max(1, qRE / (result.k - 1));   // HKSJ floor
    effSE = Math.sqrt(adjFactor / sumWRE);
    mult = tVal;
  } else {
    mult = Math.abs(qnorm(alpha / 2));
    effSE = result.se;
  }
  const ciLo = result.mu - mult * effSE;
  const ciHi = result.mu + mult * effSE;
  const denoms = studies.map(s => (s.events === 0 || s.events === s.total) ? s.total + 1 : s.total);
  const nHarm = harmonicMean(denoms);
  return {
    pi: backTransform(result.mu, transform, nHarm),
    lo: backTransform(ciLo, transform, nHarm),
    hi: backTransform(ciHi, transform, nHarm),
    tau2: result.tau2,
  };
}

const METHODS = {
  'FT+DL+z':      (st) => analyze(st, { transform: 'FT', poolMethod: 'DL', useHKSJ: false }),
  'FT+DL+HKSJ':   (st) => analyze(st, { transform: 'FT', poolMethod: 'DL', useHKSJ: true }),
  'logit+DL+z':   (st) => analyze(st, { transform: 'logit', poolMethod: 'DL', useHKSJ: false }),
  'logit+DL+HKSJ':(st) => analyze(st, { transform: 'logit', poolMethod: 'DL', useHKSJ: true }),
  'logit+REML+HKSJ':(st)=> analyze(st, { transform: 'logit', poolMethod: 'REML', useHKSJ: true }),
};

const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;

export function runCell(p0, tau2, k, scenario, reps, rng) {
  const acc = {};
  for (const n of Object.keys(METHODS)) acc[n] = { cov: 0, biasSum: 0, sq: 0, wSum: 0, n: 0 };
  for (let r = 0; r < reps; r++) {
    const { studies } = generate(p0, tau2, k, scenario, rng);
    for (const [name, fn] of Object.entries(METHODS)) {
      let out; try { out = fn(studies); } catch { continue; }
      if (!out || !isFinite(out.pi) || !isFinite(out.lo) || !isFinite(out.hi)) continue;
      const a = acc[name];
      a.n++;
      a.biasSum += out.pi - p0;
      a.sq += (out.pi - p0) ** 2;
      a.wSum += out.hi - out.lo;
      if (out.lo <= p0 && p0 <= out.hi) a.cov++;
    }
  }
  const res = {};
  for (const [name, a] of Object.entries(acc)) {
    res[name] = {
      n: a.n,
      coverage: a.n ? +(a.cov / a.n).toFixed(4) : null,
      bias: a.n ? +(a.biasSum / a.n).toFixed(4) : null,
      rmse: a.n ? +Math.sqrt(a.sq / a.n).toFixed(4) : null,
      meanWidth: a.n ? +(a.wSum / a.n).toFixed(4) : null,
    };
  }
  return res;
}

export function runGrid({ reps = 400, ks = [5, 10, 20], p0 = 0.20,
                          tau2s = [0, 0.25, 1.0], scenarios = SCENARIOS } = {}) {
  const rng = makeRng(BASE_SEED);
  const grid = [];
  for (const scen of scenarios) {
    for (const tau2 of tau2s) {
      for (const k of ks) grid.push({ scen, tau2, k, results: runCell(p0, tau2, k, scen, reps, rng) });
    }
  }
  return grid;
}

export function summarize(grid, filter = () => true) {
  const names = Object.keys(grid[0].results);
  const out = {};
  for (const name of names) {
    const cov = [], ab = [], rm = [];
    for (const c of grid) {
      if (!filter(c)) continue;
      const r = c.results[name];
      if (r.coverage != null) cov.push(r.coverage);
      if (r.bias != null) ab.push(Math.abs(r.bias));
      if (r.rmse != null) rm.push(r.rmse);
    }
    out[name] = {
      meanCoverage: cov.length ? +mean(cov).toFixed(4) : null,
      meanAbsBias: ab.length ? +mean(ab).toFixed(4) : null,
      meanRmse: rm.length ? +mean(rm).toFixed(4) : null,
    };
  }
  return out;
}

const isMain = process.argv[1]?.endsWith('harness.mjs');
if (isMain) {
  const i = process.argv.indexOf('--reps');
  const reps = i >= 0 ? Number(process.argv[i + 1]) : 400;
  const t0 = Date.now();
  const grid = runGrid({ reps });
  console.log(`\n# Truth-recovery yardstick -- proportionma`);
  console.log(`reps=${reps}/cell  p0=0.20  seed=${BASE_SEED}\n`);
  console.log('## Mean coverage of TRUE prevalence, CLEAN data (scenario=none)\n');
  const clean = summarize(grid, c => c.scen === 'none');
  console.log('method             meanCov  meanAbsBias  meanRMSE');
  for (const [n, s] of Object.entries(clean))
    console.log(n.padEnd(18), String(s.meanCoverage).padStart(7), String(s.meanAbsBias).padStart(11), String(s.meanRmse).padStart(9));
  console.log('\n## Homogeneous (tau2=0), clean -- unambiguous truth=p0\n');
  console.log('method             meanCov');
  const homo = summarize(grid, c => c.scen === 'none' && c.tau2 === 0);
  for (const [n, s] of Object.entries(homo)) console.log(n.padEnd(18), String(s.meanCoverage).padStart(7));
  console.log('\n## Per-cell coverage (clean): logit+z vs logit+HKSJ\n');
  console.log('tau2     k    logit+z   logit+HKSJ   FT+z   FT+HKSJ');
  for (const c of grid.filter(c => c.scen === 'none')) {
    console.log(String(c.tau2).padEnd(8), String(c.k).padStart(3),
      String(c.results['logit+DL+z'].coverage).padStart(9),
      String(c.results['logit+DL+HKSJ'].coverage).padStart(11),
      String(c.results['FT+DL+z'].coverage).padStart(7),
      String(c.results['FT+DL+HKSJ'].coverage).padStart(8));
  }
  console.log(`\n(${(Date.now() - t0) / 1000}s)`);
}
