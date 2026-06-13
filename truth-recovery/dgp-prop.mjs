// ============================================================
// dgp-prop.mjs -- Known-truth DGP for PROPORTION / PREVALENCE meta-analysis.
//
// proportionma pools binomial proportions via FT / logit / raw transforms with
// DL/REML + optional HKSJ. R-metafor parity proves the FORMULAS are right; it
// does NOT tell you which transform x CI method actually RECOVERS the true
// prevalence under between-study heterogeneity. This DGP supplies that test.
//
// Generative model (logit-normal random effects -- the standard prevalence-MA
// model):  theta_i = logit(p0) + sqrt(tau2) * Z_i ,  p_i = expit(theta_i),
//          n_i ~ log-uniform[nLo,nHi],  x_i ~ Binomial(n_i, p_i).
//
// ESTIMAND (truth-first, stated explicitly): a logit random-effects pool targets
// the MEDIAN prevalence = expit(mean logit) = expit(logit(p0)) = p0. So the
// honest recovery target is p0 itself. (Under tau2=0 every transform shares this
// target unambiguously; under tau2>0 the FT back-transform targets a slightly
// different functional -- that divergence is part of what we MEASURE, not hide.)
//
// Optional small-study publication selection: studies whose 95% CI excludes a
// reference (here, prevalence far from a null) are more likely to be published.
// Prevalence MA rarely has classic p-value selection, so the headline regime is
// `none`; selection cells are reported as a stress test, not the main claim.
//
// Seeded -> reproducible. Standalone (no external deps).
// ============================================================

export const SCENARIOS = ['none', 'sel_weak', 'sel_strong'];

export function makeRng(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function randn(rng) {
  let u1 = rng(), u2 = rng();
  if (u1 < 1e-12) u1 = 1e-12;
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}
const expit = (x) => 1 / (1 + Math.exp(-x));
const logit = (p) => Math.log(p / (1 - p));

function drawN(rng, nLo, nHi) {
  const lo = Math.log(nLo), hi = Math.log(nHi);
  return Math.max(5, Math.round(Math.exp(lo + (hi - lo) * rng())));
}
function rbinom(rng, n, p) {
  let x = 0;
  for (let i = 0; i < n; i++) if (rng() < p) x++;
  return x;
}

// Selection weight: extreme proportions (away from 0.5) more likely published.
const SEL = { sel_weak: 0.6, sel_strong: 0.25 };

export function generate(p0, tau2, k, scenario, rng,
                         { nLo = 30, nHi = 800, maxFactor = 200 } = {}) {
  const muLogit = logit(p0);
  const sd = Math.sqrt(tau2);
  const studies = [];
  if (scenario === 'none') {
    for (let i = 0; i < k; i++) {
      const pi = expit(muLogit + sd * randn(rng));
      const n = drawN(rng, nLo, nHi);
      studies.push({ study: `s${i}`, events: rbinom(rng, n, pi), total: n });
    }
    return { studies, p0, info: { k, selFrac: 1 } };
  }
  const wMid = SEL[scenario];
  let nExamined = 0;
  const cap = maxFactor * k;
  while (studies.length < k && nExamined < cap) {
    const pi = expit(muLogit + sd * randn(rng));
    const n = drawN(rng, nLo, nHi);
    const x = rbinom(rng, n, pi);
    const phat = x / n;
    // studies near 0.5 (uninformative) down-weighted; extreme ones favoured
    const pub = (Math.abs(phat - 0.5) > 0.25) ? 1.0 : wMid;
    nExamined++;
    if (rng() < pub) studies.push({ study: `s${studies.length}`, events: x, total: n });
  }
  while (studies.length < k) {        // top-up guarantees k
    const pi = expit(muLogit + sd * randn(rng));
    const n = drawN(rng, nLo, nHi);
    studies.push({ study: `s${studies.length}`, events: rbinom(rng, n, pi), total: n });
  }
  return { studies, p0, info: { k, selFrac: k / Math.max(1, nExamined) } };
}
