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
#> 
#> Note: Fractional Ridge Regression is a linear model without a link function.
#> Therefore, the partial effects are mathematically identical to the coefficients themselves.
#> 
print(pe)
#> 
#> Partial Effects for Fractional Ridge Regression
#> 
#> Call:
#> fracregridge.pe(object = mod)
#> 
```
