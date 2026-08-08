# Changelog

## fracreg 1.1.0

#### Output Standardization & Aesthetics

- **Hierarchical Model Displays:** Refined `print` and `summary` outputs
  for partial effects (`fracreg.pe`, `fracregpd.pe`, `fracreghet.pe`,
  `fracregridge.pe`, `fracregmlogit.pe`, `fracregmlogit.wtp`). Enforced
  a logical hierarchical model ordering (BIN0 -\> BIN1 -\> FRAC -\>
  OVERALL) and adopted clear mathematical notation for component
  probabilities (e.g., `Pr(y > 0)`).
- **Matrix Transposition:** Restructured `fracregmlogit.pe` outputs to
  universally present discrete choices as columns and covariates as
  rows, standardizing the multidimensional matrix display.
- **Specification Tests Standardization:** Outputs for econometric
  specification tests (`fracreg.reset`, `fracreghet.reset`,
  `fracreg.ggoff`, `fracreg.ptest`) now display dynamic headers
  indicating precisely which parent model was tested (e.g. “RESET Test
  for Fractional logit regression”) rather than static generic text.
- **Willingness to Pay (WTP):** Updated
  [`wtp()`](https://sulmanolieko.github.io/fracreg/reference/wtp.md) and
  its plotting functions (`plot.fracregmlogit`, `plot.fracregmlogit.pe`)
  to natively support the newly transposed `fracregmlogit.pe` outputs.
  Added stylistic formatting headers to `wtp` console prints.
- **Native Print vs Summary:** Restored native
  [`print()`](https://rdrr.io/r/base/print.html) S3 methods across the
  package to output raw point-estimate matrices while reserving detailed
  statistical tables (with standard errors, p-values, and confidence
  intervals) for [`summary()`](https://rdrr.io/r/base/summary.html) S3
  methods.

#### Tidyverse & gtsummary Integration

- **Comprehensive Broom Support:** Exported `broom::tidy()` S3 methods
  for all package models, partial effects objects (`.pe`),
  willingness-to-pay objects (`.wtp`), and econometric tests (`.reset`,
  `.ggoff`, `.ptest`). Tidied tibbles correctly map test outputs to
  `term`, `statistic`, and `p.value` columns.
- **gtsummary Dispatch Resolution:** Re-engineered the
  [`gtsummary::tbl_regression()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_regression.html)
  integration. `tbl_regression.fracreg()` now correctly triggers S3
  dispatch via [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html)
  without inadvertently stripping class attributes, preventing
  `broom::tidy()` fallback failures and resolving `R CMD check` warnings
  related to unexported `gtsummary:::` namespaces.
- **Intercept Formatting:** The `gtsummary` integration now correctly
  lists `(Intercept)` terms by default for non-mixture (1P) and mixture
  (2P) fractional models.
- **Prefix Cleaning:** Fixed an issue where `broom::tidy` retained
  `frac_` internal prefixes for `fracregridge` components.

#### Bug Fixes & Code Quality Improvements

- **Krinsky-Robb Standard Errors:** Fixed a critical bug in
  `fracregmlogit.pe` where standard errors simulated via
  `marg.type = "aveacr"` calculated point estimates across a single
  sample mean instead of averaging across all individual observations.
  Average Partial Effect standard errors are now accurate across the
  full sample. Changed `se = TRUE` by default for `fracregmlogit.pe`.
- **R CMD Check Compliance:** Replaced all instances of `T` and `F` with
  `TRUE` and `FALSE` across R source files and `roxygen2` documentation
  blocks. Added `Depends: R (>= 3.5.0)` to DESCRIPTION to comply with
  serialized dataset formats.
- **Documentation Tags:** Systematically added missing `\value` tags to
  exported R methods (like `fitted`, `predict`, and `residuals`) and
  removed redundant introductory phrases like “This package…” from
  description fields to meet strict CRAN policies.
- **Example Scripts & Variances:** Fixed simulation bugs in
  `fracregpd.pe` example scripts where random error terms possessed
  static variances, restoring `R CMD check` validation.
- **GitHub Actions:** Resolved `403 Forbidden` errors during `pkgdown`
  automated documentation deployments by explicitly granting the
  `contents: write` permission to the GitHub Actions workflow.

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
