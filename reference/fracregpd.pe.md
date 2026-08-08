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

## Examples

``` r
# \donttest{
# Simulate Panel Data
N <- 50
T <- 5
id <- rep(1:N, each=T)
time <- rep(1:T, N)
x1 <- rnorm(N*T)
x2 <- rnorm(N*T)
z1 <- rnorm(N*T)
u <- rnorm(N*T)
y <- exp(x1 - x2 + u)/(1 + exp(x1 - x2 + u))
X <- cbind(x1 = x1, x2 = x2)
Z <- cbind(x1 = x1, z1 = z1)

# Estimate panel data model
mod_pd <- fracregpd(id=id, time=time, y=y, x=X, z=Z, 
                    x.exogenous=FALSE, type="GMMbgw", link="logit")

# Compute Average Partial Effects
pe_res <- fracregpd.pe(mod_pd)
summary(pe_res)
#> 
#> -------------------------------------------------------------------------------- 
#>                             Average partial effects 
#> -------------------------------------------------------------------------------- 
#>                      Panel data fractional logit regression 
#> -------------------------------------------------------------------------------- 
#> 
#> Note: Standard errors computed using the Delta method
#>         dy/dx Std. Error z value Pr(>|z|)
#> x1  9.893e-02  1.027e+07       0        1
#> x2 -2.073e-01  2.300e+07       0        1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-08-08 13:39:29 
#> -------------------------------------------------------------------------------- 
# }
```
