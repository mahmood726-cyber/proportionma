# Truth-recovery yardstick — proportionma

**Verdict: STRONG VALIDATION + actionable recommendation (measured).**

proportionma is R-`meta::metaprop`-validated for its point outputs. That proves the
*formulas* are right; it does not show which **transform × CI-method** combination
actually **recovers the true prevalence** under between-study heterogeneity. This
harness supplies that missing evidence by injecting a known true prevalence `p0`
(logit-normal random effects, binomial counts) and measuring how often each
interval covers `p0`.

## Method
- DGP (`dgp-prop.mjs`): `theta_i = logit(p0) + sqrt(tau2)·Z`, `p_i = expit(theta_i)`,
  `n_i ~ log-uniform[30,800]`, `x_i ~ Binomial(n_i, p_i)`. Seeded → reproducible.
- Estimand stated explicitly: a logit RE pool targets the **median** prevalence
  `expit(mean logit) = p0`. Under `tau2=0` every transform shares this target
  unambiguously.
- Engine (`engine.mjs`): the app's own pooling/CI functions copied **verbatim**
  from `index.html` (FT/logit transforms, DL/REML, HKSJ with `max(1,·)` floor,
  Miller back-transform). The SAME math the app ships is measured.
- 400 reps/cell, `p0=0.20`, `k∈{5,10,20}`, `tau2∈{0,0.25,1.0}`, seed 20260613.

## Headline numbers (clean data, mean over the grid)

| method            | meanCov | meanAbsBias | meanRMSE |
|-------------------|--------:|------------:|---------:|
| FT + DL + z       |  0.9144 |      0.0096 |   0.0328 |
| FT + DL + HKSJ    |  0.9594 |      0.0096 |   0.0328 |
| logit + DL + z    |  0.9119 |      0.0057 |   0.0318 |
| logit + DL + HKSJ |  0.9617 |      0.0057 |   0.0318 |
| logit + REML+HKSJ |  0.9614 |      0.0044 |   0.0314 |

Per-cell coverage under heterogeneity (clean):

| tau2 | k | logit+z | logit+HKSJ | FT+z | FT+HKSJ |
|-----:|--:|--------:|-----------:|-----:|--------:|
| 0    | 5 |  0.953  |   1.000    |0.955 | 1.000   |
| 0.25 | 5 |  0.865  |   0.955    |0.868 | 0.960   |
| 1.0  | 5 |  0.863  |   0.950    |0.865 | 0.945   |
| 1.0  | 20|  0.878  |   0.930    |0.885 | 0.915   |

## Findings (all measured)
1. **The Wald-z interval (the app default, `useHKSJ=false`) under-covers true
   prevalence under heterogeneity** — down to ~0.86 at `tau2≥0.25`, well below the
   nominal 0.95. This is the classic DL-Wald failure with few studies.
2. **The HKSJ toggle is a real, measured fix**: enabling it restores coverage to
   ~0.95 across the heterogeneous cells (+5pp mean clean coverage; +9pp at
   `tau2=1, k=5`), and never recovers truth *worse* than z. At `tau2=0, k=5` it is
   mildly conservative (≈1.0), the expected small-k HKSJ behaviour.
   → **Recommendation: enable HKSJ by default** (it is currently opt-in). This is
   consistent with the advanced-stats rule "HKSJ preferred for k<30".
3. **FT is NOT worse than logit on coverage** — the two transforms land within
   ~1pp of each other. logit carries slightly lower bias (0.0057 vs 0.0096) and
   RMSE. The honest lever for honest intervals here is the **CI method (HKSJ), not
   the transform** — the same correction recorded for the sibling app
   `Pairwiseprohtml` in the first wave.
4. REML+HKSJ ≈ DL+HKSJ on coverage with marginally lower bias — no penalty for the
   library's REML default.

## What did NOT transfer
NPE / conformal / SBC / PartialID are estimator-of-μ machinery; proportionma is a
deterministic transform+pool formula library, so only the **known-truth harness**
and the heterogeneity DGP transferred. No runtime dependency added.

## Reproduce
```
node truth-recovery/harness.mjs --reps 400
node --test truth-recovery/test-truth-recovery.mjs
```
