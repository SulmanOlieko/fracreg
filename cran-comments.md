## Resubmission

This is a resubmission. I have bumped the package to version 1.0.1 to include new features detailed in `NEWS.md`. In this version, I have also fully addressed the previous CRAN feedback:

* Added the required `\value` tag to `fitted.fracregmlogit.Rd` explicitly describing the output structure (`data.frame` class) and meaning for the exported `fitted`, `residuals`, and `predict` S3 methods on `fracregmlogit` objects.
* Removed the starting phrase "This package" and similar variations from the description field.
* Expanded the RESET acronym to "Regression Equation Specification Error Test" in the description field.
* Added method references to the description field in the required `authors (year) <doi:...>` format.
* Replaced all instances of `T` and `F` with `TRUE` and `FALSE` in function default arguments across the R source files and corresponding Rd documentation.

## Test environments

* local macOS install, R-devel
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Sulman Olieko Owili <oliekosulman@gmail.com>'
  New submission
