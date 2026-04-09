# ProportionMA

Browser-based proportion/prevalence meta-analysis tool. No server, no dependencies -- runs entirely in a single HTML file.

## Features

- **Transformations**: Freeman-Tukey double arcsine, logit, raw proportion
- **Pooling methods**: DerSimonian-Laird (moment-based) and REML (iterative)
- **HKSJ adjustment**: Hartung-Knapp-Sidik-Jonkman with t-distribution and Q/(k-1) floor
- **Prediction interval**: Using t_{k-2} distribution (undefined for k < 3)
- **Per-study CIs**: Clopper-Pearson exact intervals
- **Zero handling**: 0.5 continuity correction applied only when needed
- **Subgroup analysis**: Stratified pooling with between-subgroup Q test
- **Heterogeneity**: Q statistic, I-squared, tau-squared
- **Forest plot**: SVG with proportional weight squares, summary diamond, optional PI line
- **Funnel plot**: SVG with pseudo-95% CI contour
- **Export**: CSV (study data + results), SVG (forest plot)

## Usage

1. Open `index.html` in any modern browser
2. Paste CSV data (Study, Events, Total, optional Subgroup) or use manual entry
3. Select transformation, method, and optional adjustments
4. Click Analyze

### CSV Format

```
Study,Events,Total,Subgroup
Smith 2020,15,100,Group A
Jones 2021,22,150,Group A
Lee 2019,8,50,Group B
```

## Statistical Details

### Back-transformation
- Freeman-Tukey: Miller (1978) correction using harmonic mean of denominators
- Logit: inverse logit (expit)
- Raw: identity (clamped to [0, 1])

### HKSJ
- Uses t-distribution with k-1 degrees of freedom (not normal approximation)
- Applies floor: max(1, Q/(k-1)) to prevent narrowing below DL

### Prediction Interval
- Uses t-distribution with k-2 degrees of freedom
- Undefined for k < 3 (displays message)

## Testing

```bash
cd C:\Models\ProportionMA
python -m pytest test_app.py -v
```

Requires: Python, pytest, selenium, Chrome/ChromeDriver.

29 tests covering: statistical correctness, edge cases, UI functionality, exports, input validation.

## Author

Mahmood Ahmad, Tahir Heart Institute
