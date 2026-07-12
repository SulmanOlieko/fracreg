# Partial Effects for Fractional Panel Data Regression

Computes Average Partial Effects (APEs) for fractional panel data models
fit using `fracregpd`. In correlated random effects (CRE) models, the
time-averages of the covariates are used merely to control for
unobserved heterogeneity, and as such, they are automatically filtered
out when computing the structural partial effects. Standard errors are
computed using the Delta method and the bootstrapped parameter variance
matrix.

## Usage

``` r
fracregpd.pe(
  object,
  APE = TRUE,
  CPE = FALSE,
  at = NULL,
  which.x = NULL,
  variance = TRUE,
  table = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `fracregpd`.

- APE:

  logical. Compute Average Partial Effects?

- CPE:

  logical. Compute Conditional Partial Effects? (Not currently supported
  for panel data models).

- at:

  numeric vector. The values at which to evaluate the CPE.

- which.x:

  character vector. Variables for which to compute partial effects. By
  default, auxiliary CRE parameters (like `_mean` and `vhat`) are
  automatically excluded so that partial effects are only computed for
  the main structural parameters.

- variance:

  logical. Compute standard errors using the Delta method?

- table:

  logical. Print the resulting partial effects table?

- ...:

  further arguments passed to or from other methods.

## Value

An object of class `fracreg.pe` containing the standard coefficient
tables with Average Partial Effects.

## See also

[`fracregpd`](https://sulmanolieko.github.io/fracreg/reference/fracregpd.md),
[`fracreg.pe`](https://sulmanolieko.github.io/fracreg/reference/fracreg.pe.md)

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>
