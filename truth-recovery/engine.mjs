// engine.mjs -- pure numerical core EXTRACTED VERBATIM from index.html
// (proportionma). Only DOM-free functions are copied; exported for the
// truth-recovery harness so the SAME math the app ships is measured.

function qt(p, df) {
  // Returns quantile for probability p and degrees of freedom df
  // For meta-analysis we need qt(alpha/2, df) which is negative
  if (df <= 0) return NaN;
  if (df === 1) return Math.tan(Math.PI * (p - 0.5));
  if (df === 2) return (2 * p - 1) / Math.sqrt(2 * p * (1 - p));

  // Use Abramowitz & Stegun rational approximation via normal approx + correction
  const a = 1 / (df - 0.5);
  const b = 48 / (a * a);
  let z = qnorm(p);
  const y = z * z;

  // Cornish-Fisher expansion
  let x = z * (1 + (y - 3) / (4 * df) +
    ((5 * y * y - 56 * y + 75) / (96 * df * df)) +
    ((-3 * y * y * y + 105 * y * y - 735 * y + 945) / (384 * df * df * df)));

  // Newton-Raphson refinement
  for (let iter = 0; iter < 5; iter++) {
    const tp = tCDF(x, df) - p;
    const td = tPDF(x, df);
    if (Math.abs(td) < 1e-20) break;
    x = x - tp / td;
  }
  return x;
}

function qnorm(p) {
  // Rational approximation for the standard normal inverse CDF
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;

  const a = [
    -3.969683028665376e+01, 2.209460984245205e+02,
    -2.759285104469687e+02, 1.383577518672690e+02,
    -3.066479806614716e+01, 2.506628277459239e+00
  ];
  const b = [
    -5.447609879822406e+01, 1.615858368580409e+02,
    -1.556989798598866e+02, 6.680131188771972e+01,
    -1.328068155288572e+01
  ];
  const c = [
    -7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,
     4.374664141464968e+00, 2.938163982698783e+00
  ];
  const d = [
    7.784695709041462e-03, 3.224671290700398e-01,
    2.445134137142996e+00, 3.754408661907416e+00
  ];

  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  let q, r;

  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}

function tPDF(x, df) {
  const lnCoeff = lnGamma((df + 1) / 2) - lnGamma(df / 2) - 0.5 * Math.log(df * Math.PI);
  return Math.exp(lnCoeff - ((df + 1) / 2) * Math.log(1 + x * x / df));
}

function tCDF(x, df) {
  // Regularized incomplete beta function approach
  if (x === 0) return 0.5;
  const t2 = x * x;
  const ib = regIncBeta(df / 2, 0.5, df / (df + t2));
  return x > 0 ? 1 - 0.5 * ib : 0.5 * ib;
}

// --- Chi-squared CDF ---
function chiSqCDF(x, df) {
  if (x <= 0) return 0;
  return regGammaP(df / 2, x / 2);
}

// --- Log Gamma (Stirling) ---
function lnGamma(x) {
  if (x <= 0) return Infinity;
  // Lanczos approximation
  const g = 7;
  const coef = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
  ];
  if (x < 0.5) {
    return Math.log(Math.PI / Math.sin(Math.PI * x)) - lnGamma(1 - x);
  }
  x -= 1;
  let a = coef[0];
  const t = x + g + 0.5;
  for (let i = 1; i < coef.length; i++) {
    a += coef[i] / (x + i);
  }
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

function gamma(x) { return Math.exp(lnGamma(x)); }

// --- Regularized incomplete beta function ---
function regIncBeta(a, b, x) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;

  // Use symmetry relation for better convergence
  if (x > (a + 1) / (a + b + 2)) {
    return 1 - regIncBeta(b, a, 1 - x);
  }

  const lnBeta = lnGamma(a) + lnGamma(b) - lnGamma(a + b);
  const front = Math.exp(Math.log(x) * a + Math.log(1 - x) * b - lnBeta) / a;

  // Lentz's continued fraction
  let f = 1, c = 1, d = 1 - (a + b) * x / (a + 1);
  if (Math.abs(d) < 1e-30) d = 1e-30;
  d = 1 / d;
  f = d;

  for (let m = 1; m <= 200; m++) {
    // even step
    let numerator = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m));
    d = 1 + numerator * d;
    if (Math.abs(d) < 1e-30) d = 1e-30;
    c = 1 + numerator / c;
    if (Math.abs(c) < 1e-30) c = 1e-30;
    d = 1 / d;
    f *= c * d;

    // odd step
    numerator = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1));
    d = 1 + numerator * d;
    if (Math.abs(d) < 1e-30) d = 1e-30;
    c = 1 + numerator / c;
    if (Math.abs(c) < 1e-30) c = 1e-30;
    d = 1 / d;
    const delta = c * d;
    f *= delta;
    if (Math.abs(delta - 1) < 1e-12) break;
  }

  return front * f;
}

// --- Beta quantile (inverse CDF) via bisection + Newton ---
function qbeta(p, a, b) {
  if (p <= 0) return 0;
  if (p >= 1) return 1;
  if (a <= 0 || b <= 0) return NaN;

  // Initial guess using normal approximation
  let x = 0.5;
  // Better initial guess
  const mu = a / (a + b);
  x = Math.max(0.001, Math.min(0.999, mu));

  // Newton-Raphson
  for (let iter = 0; iter < 100; iter++) {
    const fx = regIncBeta(a, b, x) - p;
    if (Math.abs(fx) < 1e-12) break;

    // Beta PDF
    const lnBeta = lnGamma(a) + lnGamma(b) - lnGamma(a + b);
    const pdf = Math.exp((a - 1) * Math.log(Math.max(x, 1e-300)) + (b - 1) * Math.log(Math.max(1 - x, 1e-300)) - lnBeta);

    if (pdf < 1e-20) {
      // Bisection fallback
      if (fx < 0) x = x + (1 - x) / 2;
      else x = x / 2;
    } else {
      let step = fx / pdf;
      // Clamp step
      if (x - step <= 0) step = x * 0.9;
      if (x - step >= 1) step = -(1 - x) * 0.9;
      x = x - step;
      x = Math.max(1e-15, Math.min(1 - 1e-15, x));
    }
  }
  return x;
}

// --- Regularized lower incomplete gamma (for chi-sq CDF) ---
function regGammaP(a, x) {
  if (x <= 0) return 0;
  if (x < a + 1) {
    // Series expansion
    let sum = 1 / a, term = 1 / a;
    for (let n = 1; n < 200; n++) {
      term *= x / (a + n);
      sum += term;
      if (Math.abs(term) < Math.abs(sum) * 1e-12) break;
    }
    return sum * Math.exp(-x + a * Math.log(x) - lnGamma(a));
  } else {
    // Continued fraction
    return 1 - regGammaQ(a, x);
  }
}

function regGammaQ(a, x) {
  // Continued fraction for upper incomplete gamma
  let f = 1, c = 1, d = x + 1 - a;
  if (Math.abs(d) < 1e-30) d = 1e-30;
  d = 1 / d;
  f = d;

  for (let i = 1; i <= 200; i++) {
    const an = i * (a - i);
    const bn = x + 2 * i + 1 - a;
    d = bn + an * d;
    if (Math.abs(d) < 1e-30) d = 1e-30;
    c = bn + an / c;
    if (Math.abs(c) < 1e-30) c = 1e-30;
    d = 1 / d;
    const delta = c * d;
    f *= delta;
    if (Math.abs(delta - 1) < 1e-12) break;
  }

  return Math.exp(-x + a * Math.log(x) - lnGamma(a)) * f;
}


// ============================================================
// TRANSFORMATIONS
// ============================================================

function transformStudies(studies, method) {
  // Apply 0.5 continuity correction ONLY to studies with x=0 or x=n
  const corrected = studies.map(s => {
    let x = s.events, n = s.total;
    if (s.events === 0 || s.events === s.total) {
      x = x + 0.5;
      n = n + 1;
    }
    return { ...s, x_corr: x, n_corr: n };
  });

  return corrected.map(s => {
    const x = s.x_corr;
    const n = s.n_corr;
    const p = x / n;
    let yi, vi;

    if (method === 'FT') {
      // Freeman-Tukey double arcsine
      // Use corrected n (n_corr) which already includes +1 for zero-cell studies
      yi = Math.asin(Math.sqrt(x / (n + 1))) + Math.asin(Math.sqrt((x + 1) / (n + 1)));
      vi = 1 / (n + 0.5);
    } else if (method === 'logit') {
      const pClamped = Math.max(1e-10, Math.min(1 - 1e-10, p));
      yi = Math.log(pClamped / (1 - pClamped));
      vi = 1 / (n * pClamped * (1 - pClamped));
    } else {
      // raw
      yi = p;
      vi = p * (1 - p) / n;
      if (vi < 1e-20) vi = 1e-10; // avoid zero variance
    }

    return { ...s, yi, vi };
  });
}

// ============================================================
// POOLING
// ============================================================

function poolDL(transformed) {
  const k = transformed.length;
  if (k === 0) return null;

  const w = transformed.map(s => 1 / s.vi);
  const sumW = w.reduce((a, b) => a + b, 0);
  const yBar = w.reduce((s, wi, i) => s + wi * transformed[i].yi, 0) / sumW;

  if (k === 1) {
    return {
      mu: yBar, tau2: 0, se: Math.sqrt(transformed[0].vi),
      Q: 0, I2: 0, Qp: 1, k: 1, w, weights: [100]
    };
  }

  // Q statistic
  const Q = w.reduce((s, wi, i) => s + wi * Math.pow(transformed[i].yi - yBar, 2), 0);
  const df = k - 1;
  const C = sumW - w.reduce((s, wi) => s + wi * wi, 0) / sumW;
  let tau2 = Math.max(0, (Q - df) / C);

  // RE weights
  const wStar = transformed.map(s => 1 / (s.vi + tau2));
  const sumWStar = wStar.reduce((a, b) => a + b, 0);
  const mu = wStar.reduce((s, wi, i) => s + wi * transformed[i].yi, 0) / sumWStar;
  const se = Math.sqrt(1 / sumWStar);

  const I2 = Math.max(0, (Q - df) / Q) * 100;
  const Qp = 1 - chiSqCDF(Q, df);

  const totalW = wStar.reduce((a, b) => a + b, 0);
  const weights = wStar.map(wi => (wi / totalW) * 100);

  return { mu, tau2, se, Q, I2, Qp, k, w: wStar, weights };
}

function poolREML(transformed) {
  const k = transformed.length;
  if (k === 0) return null;
  if (k === 1) return poolDL(transformed);

  // Start from DL estimate
  const dlResult = poolDL(transformed);
  let tau2 = dlResult.tau2;

  // Fisher scoring iterations
  for (let iter = 0; iter < 100; iter++) {
    const wStar = transformed.map(s => 1 / (s.vi + tau2));
    const sumW = wStar.reduce((a, b) => a + b, 0);
    const mu = wStar.reduce((s, wi, i) => s + wi * transformed[i].yi, 0) / sumW;

    // First derivative of REML log-likelihood
    const dl = -0.5 * wStar.reduce((s, wi) => s + wi * wi / (wi), 0) / sumW +
      0.5 * wStar.reduce((s, wi, i) => {
        const r = transformed[i].yi - mu;
        return s + wi * wi * r * r;
      }, 0) / (sumW * sumW) * sumW;

    // Simplified: gradient and info
    let gradient = 0, info = 0;
    const wk = transformed.map(s => 1 / (s.vi + tau2));
    const sumWk = wk.reduce((a, b) => a + b, 0);
    const muK = wk.reduce((s, wi, i) => s + wi * transformed[i].yi, 0) / sumWk;

    for (let i = 0; i < k; i++) {
      const wi = wk[i];
      const ri = transformed[i].yi - muK;
      gradient += -wi * wi * (1 / wi - ri * ri);
      info += wi * wi;
    }
    gradient = -0.5 * (wk.reduce((s, wi) => s - wi, 0) + wk.reduce((s, wi, i) => {
      const r = transformed[i].yi - muK;
      return s + wi * wi * r * r;
    }, 0));

    // Actual REML score equation
    let score = 0;
    let fisher = 0;
    for (let i = 0; i < k; i++) {
      const wi2 = wk[i] * wk[i];
      score += wi2 * ((transformed[i].yi - muK) * (transformed[i].yi - muK) - (1 / wk[i])) +
        wi2 / sumWk;
    }
    score *= 0.5;

    for (let i = 0; i < k; i++) {
      fisher += wk[i] * wk[i];
    }
    fisher = 0.5 * (fisher - wk.reduce((s, wi) => s + wi * wi, 0) * wk.reduce((s, wi) => s + wi * wi, 0) / (sumWk * sumWk));

    // Simplified REML update: Paule-Mandel style iteration
    const Qk = wk.reduce((s, wi, i) => s + wi * (transformed[i].yi - muK) * (transformed[i].yi - muK), 0);
    const sumW2overW = wk.reduce((s, wi) => s + wi * wi, 0) / sumWk;
    const newTau2 = Math.max(0, tau2 + (Qk - (k - 1)) / (sumWk - sumW2overW));

    if (Math.abs(newTau2 - tau2) < 1e-8) {
      tau2 = newTau2;
      break;
    }
    tau2 = newTau2;
  }

  // Final estimates
  const wStar = transformed.map(s => 1 / (s.vi + tau2));
  const sumWStar = wStar.reduce((a, b) => a + b, 0);
  const mu = wStar.reduce((s, wi, i) => s + wi * transformed[i].yi, 0) / sumWStar;
  const se = Math.sqrt(1 / sumWStar);

  // Q using FE weights for heterogeneity test
  const wFE = transformed.map(s => 1 / s.vi);
  const sumWFE = wFE.reduce((a, b) => a + b, 0);
  const muFE = wFE.reduce((s, wi, i) => s + wi * transformed[i].yi, 0) / sumWFE;
  const Q = wFE.reduce((s, wi, i) => s + wi * Math.pow(transformed[i].yi - muFE, 2), 0);
  const df = k - 1;
  const I2 = Math.max(0, (Q - df) / Q) * 100;
  const Qp = 1 - chiSqCDF(Q, df);

  const totalW = wStar.reduce((a, b) => a + b, 0);
  const weights = wStar.map(wi => (wi / totalW) * 100);

  return { mu, tau2, se, Q, I2, Qp, k, w: wStar, weights };
}

// ============================================================
// BACK-TRANSFORMATION
// ============================================================

function backTransform(y, method, nHarm) {
  if (method === 'FT') {
    // Miller (1978) back-transformation using harmonic mean of denominators
    const sinVal = Math.sin(y / 2);
    let p = sinVal * sinVal;
    // Miller correction
    if (nHarm > 0) {
      p = p - (1 / nHarm) * (p * (1 - p) - y / (4 * nHarm));
    }
    return Math.max(0, Math.min(1, p));
  } else if (method === 'logit') {
    return 1 / (1 + Math.exp(-y));
  } else {
    return Math.max(0, Math.min(1, y));
  }
}

function harmonicMean(arr) {
  if (arr.length === 0) return 0;
  const sumInv = arr.reduce((s, n) => s + 1 / n, 0);
  return arr.length / sumInv;
}

// ============================================================
// CLOPPER-PEARSON EXACT CI
// ============================================================

function clopperPearsonCI(x, n, alpha) {
  if (n <= 0) return { lower: 0, upper: 1 };
  let lower, upper;

  if (x === 0) {
    lower = 0;
    upper = qbeta(1 - alpha / 2, x + 1, n - x);
  } else if (x === n) {
    lower = qbeta(alpha / 2, x, n - x + 1);
    upper = 1;
  } else {
    lower = qbeta(alpha / 2, x, n - x + 1);
    upper = qbeta(1 - alpha / 2, x + 1, n - x);
  }

  return {
    lower: Math.max(0, lower),
    upper: Math.min(1, upper)
  };
}

export { qt, qnorm, chiSqCDF, qbeta, transformStudies, poolDL, poolREML, backTransform, harmonicMean, clopperPearsonCI };
