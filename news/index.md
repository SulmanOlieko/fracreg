# Changelog

## fracreg 1.0.1

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

#### Documentation & Website

- **Citations & References:** Added robust citations for the
  foundational `fmlogit` (Ji and Woodill) and `fracridge` (Rokem and
  Kay) packages across manual pages, README, and vignettes, formatting
  URLs properly for rendering.
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
  fixing missing alt-text for accessibility.

## fracreg 1.0.0

- Initial CRAN release.
