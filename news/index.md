# Changelog

## fracreg 1.1.0

#### Output Standardization & Minor Fixes

- **Output Standardization:** Refined `print` and `summary` outputs for
  partial effects (`fracreg.pe`, `fracregpd.pe`, `fracreghet.pe`,
  `fracregridge.pe`, `fracregmlogit.pe`, `fracregmlogit.wtp`). Enforced
  a logical hierarchical model ordering (BIN0 -\> BIN1 -\> FRAC -\>
  OVERALL), adopted mathematical notation for component probabilities
  (e.g., `Pr(y > 0)`), and restructured `fracregmlogit.pe` outputs to
  universally present choices as columns and variables as rows. Added
  explanatory footnotes for omitted reduced-form parameters in
  `fracreghet.pe` outputs. Restored native
  [`print()`](https://rdrr.io/r/base/print.html) methods for
  point-estimates matrix displays while reserving full standard error
  tables for [`summary()`](https://rdrr.io/r/base/summary.html).
- **Plotting & WTP:** Updated
  [`wtp()`](https://sulmanolieko.github.io/fracreg/reference/wtp.md) and
  the respective plotting functions (`plot.fracregmlogit` and
  `plot.fracregmlogit.pe`) to natively support the newly transposed
  `fracregmlogit.pe` outputs, ensuring dimensional alignment across
  partial effects calculation pipelines. Added a specific check ensuring
  `indv.obs` arguments correctly validate the marginal effects method
  before allocating memory. Added stylistic headers (“Fractional
  multinomial logit regression”, “Willingness to Pay”) to the print
  output of
  [`wtp()`](https://sulmanolieko.github.io/fracreg/reference/wtp.md).
  Exported `broom::tidy()` and
  [`summary()`](https://rdrr.io/r/base/summary.html) for `.wtp` objects
  to the namespace.
- **Krinsky-Robb Simulation Bug Fix:** Fixed a critical bug in
  `fracregmlogit.pe` where Krinsky-Robb standard errors for
  `marg.type = "aveacr"` were incorrectly calculating point estimates
  across the single sample mean instead of iterating across all
  observations. Average Partial Effect standard errors are now correctly
  averaged across the full sample. Changed `se = TRUE` by default for
  `fracregmlogit.pe`.
- **gtsummary Intercept Output:** Introduced an explicit S3 method
  `tbl_regression.fracreg()` so that intercept terms are displayed by
  default when formatting `1P` and `2P` models with the `gtsummary`
  package, resolving a discrepancy where intercepts were omitted in
  non-mixture models.
- Removed redundant raw [`cat()`](https://rdrr.io/r/base/cat.html) print
  statements from the computation body of `fracregridge.pe`.
- Fixed an issue where `broom::tidy` outputs retained `frac_` prefixes
  in component columns for `fracregridge` models, ensuring clean
  integration with `tibble` and `gtsummary`.

## fracreg 1.0.1

CRAN release: 2026-08-05

#### New Features & Enhancements

- **Fractional Ridge Regression:** Added `fracregridge` and
  `fracregridge.pe` to implement fractional ridge regression and its
  corresponding partial effects.
- **Fractional Multinomial Logit & Krinsky-Robb:** Refactored
  `fracregmlogit.pe` to use
  [`MASS::mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) for
  Krinsky-Robb standard error simulations, resolving hardcoded
  constraints.
- **Panel Data Partial Effects:** Implemented `fracregpd.pe` for
  computing Average Partial Effects (APE) in panel data settings,
  including full support for Correlated Random Effects (CRE) models and
  Delta method standard errors.
- **Odds Ratios:** Added support for computing Odds Ratios and
  adjustable confidence intervals (`level`) in output summaries.
- **Offset Parameter:** Added support for an `offset` parameter in model
  specifications.
- **S3 Methods & Output Display:** Standardised console output by adding
  dedicated [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`print()`](https://rdrr.io/r/base/print.html) S3 methods across the
  package to align with R’s standard practices. Changed output tables to
  say “regression” rather than “model” for clarity.
- **Econometric Terminology & Standardization:** Renamed `Estimate`
  column headers to `dy/dx` universally across all partial effects
  tables (`fracreg.pe`, `fracreghet.pe`, `fracregmlogit.pe`,
  `fracregpd.pe`, `fracregridge.pe`), while enforcing `Coefficient` for
  non-partial effects outputs.
- **Standard Error Annotations:** Output tables now dynamically label
  column headers according to standard error type
  (e.g. `Robust Std.Err.`, `Cluster Std.Err.`) and cleanly append
  specific methodological notes to the tables (e.g., Delta Method,
  Krinsky-Robb, HC0, CRVE).

#### Documentation & Website

- **Roxygen2 Migration:** Fully automated package documentation by
  migrating all manual `.Rd` files into `roxygen2` comments embedded
  directly within the R source files, ensuring airtight synchronization
  between code and documentation.
- **Citations & References:** Added citations for the foundational
  `fmlogit` (Ji and Woodill) and `fracridge` (Rokem and Kay) packages
  across manual pages, README, and vignettes, formatting URLs properly
  for rendering.
- **Logo & Badges:** Improved the package hex logo and added status
  badges to the `README`.
- **Website (`pkgdown`):** Added a `sitemap.xml` and reordered the
  `README` and `vignettes` sections to properly showcase Fractional
  Multinomial Logit and Fractional Ridge Regression alongside Panel Data
  models. Added Willingness to Pay (WTP) plot visualizations in
  tutorials.
- **Help Files:** Updated `.Rd` files to properly reflect
  `table = FALSE` default arguments and fixed `\usage` widths. Updated
  the package title for clarity and consistency.

#### Bug Fixes & Under-the-Hood Improvements

- **Namespace Imports:** Added comprehensive `@importFrom` tags for base
  R `stats`, `graphics`, and `grid` functions, alongside `ggplot2`,
  eliminating all undefined global function NOTEs during checking.

- **S3 Method Registration:** Explicitly registered S3 methods (like
  `summary` and `plot`) using `@exportS3Method` to comply with strict
  `R CMD check` standards.

- **CRAN Compliance:** Wrapped computationally heavy examples (like
  `plot.fracregmlogit.pe`) in `\donttest{}` blocks to strictly respect
  CRAN’s 5-second execution limit.

- Resolved `length > 1` logical condition errors in `plot.fracregmlogit`
  standard error rendering by properly adopting
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

- Resolved WTP matrix subsetting issues by enforcing `drop = FALSE`.

- Cleaned up legacy `fmlogit` code dependencies and unified the
  namespace fully under `fracregmlogit`.

- Fixed `NA/NaN` evaluation errors in `nlminb` during `fracreghet`
  estimations.

- Resolved `devtools::check()` formatting issues, including S3 class
  inheritance checks
  ([`inherits()`](https://rdrr.io/r/base/class.html)),
  [`stats::approx`](https://rdrr.io/r/stats/approxfun.html) scoping, and
  fixing missing alt-text for accessibility. \# fracreg 1.0.0

- Initial CRAN release.
