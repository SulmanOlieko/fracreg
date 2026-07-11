# fracreg 1.0.1

### New Features & Enhancements
* **Fractional Ridge Regression:** Added `fracregridge` and `fracregridge.pe` to implement fractional ridge regression and its corresponding partial effects.
* **Panel Data Partial Effects:** Implemented `fracregpd.pe` for computing Average Partial Effects (APE) in panel data settings, including full support for Correlated Random Effects (CRE) models and Delta method standard errors.
* **Odds Ratios:** Added support for computing Odds Ratios and adjustable confidence intervals (`level`) in output summaries.
* **Offset Parameter:** Added support for an `offset` parameter in model specifications.
* **S3 Methods & Output Display:** Standardised console output by adding dedicated `summary()` and `print()` S3 methods across the package to align with R's standard practices. Changed output tables to say "regression" rather than "model" for clarity.

### Documentation & Website
* **Logo & Badges:** Improved the package hex logo and added status badges to the `README`.
* **Website (`pkgdown`):** Added a `sitemap.xml` and reordered the `README` and `vignettes` sections to properly showcase Fractional Ridge Regression alongside Panel Data models.
* **Help Files:** Updated `.Rd` files to properly reflect `table = FALSE` default arguments and fixed `\usage` widths. Updated the package title for clarity and consistency.

### Bug Fixes & Under-the-Hood Improvements
* Fixed `NA/NaN` evaluation errors in `nlminb` during `fracreghet` estimations.
* Resolved `devtools::check()` formatting issues, including S3 class inheritance checks (`inherits()`), `stats::approx` scoping, and fixing missing alt-text for accessibility.

# fracreg 1.0.0

* Initial CRAN release.
