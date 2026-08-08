# Partial Effects for Fractional Ridge Regression

Because Fractional Ridge Regression fits a linear model without a link
function, the partial effects are mathematically identical to the
estimated ridge coefficients. This function serves as a wrapper to
maintain API compatibility with the rest of the `fracreg` package,
printing a brief notification and returning the standard coefficient
tables.

## Usage

``` r
fracregridge.pe(
  object,
  APE = TRUE,
  CPE = FALSE,
  at = NULL,
  variance = TRUE,
  table = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `fracregridge`.

- APE:

  logical. Ignored for ridge regression.

- CPE:

  logical. Ignored for ridge regression.

- at:

  numeric vector. Ignored for ridge regression.

- variance:

  logical. Ignored for ridge regression.

- table:

  logical. Ignored for ridge regression.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class `fracreg.pe` containing the standard coefficient
tables.

## See also

[`fracregridge`](https://sulmanolieko.github.io/fracreg/reference/fracregridge.md),
[`fracreg.pe`](https://sulmanolieko.github.io/fracreg/reference/fracreg.pe.md)

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## Examples

``` r
# Generate random data
set.seed(123)
y <- rnorm(100)
X <- matrix(rnorm(1000), 100, 10)
colnames(X) <- paste0("X", 1:10)

# Fit Fractional Ridge Regression
mod <- fracregridge(y, X, fracs = c(0.3, 0.5))

# Compute Partial Effects (identical to coefficients)
pe <- fracregridge.pe(mod)
print(pe)
#> 
#> Fractional ridge regression
#> 
#> Average partial effects:
#> 
#> Note: Fractional ridge regression is a linear model without a link function.
#> Therefore, the partial effects are mathematically identical to the coefficients themselves.
#> 
#> Target Fraction: 0.3
#>                    dy/dx
#> (Intercept)  0.025969041
#> X1          -0.015190037
#> X2          -0.029634697
#> X3          -0.017972923
#> X4          -0.048194568
#> X5          -0.013807840
#> X6          -0.013152916
#> X7           0.048353290
#> X8          -0.006353284
#> X9           0.001567879
#> X10          0.020836124
#> 
#> Target Fraction: 0.5
#>                    dy/dx
#> (Intercept)  0.043431002
#> X1          -0.025277903
#> X2          -0.050625055
#> X3          -0.035289788
#> X4          -0.081684776
#> X5          -0.021406027
#> X6          -0.022143313
#> X7           0.078031858
#> X8          -0.014086532
#> X9           0.002529806
#> X10          0.031879315
#> 
```
