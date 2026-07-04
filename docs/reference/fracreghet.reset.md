# RESET Test for Fractional Regression Models under Neglected Heterogeneity

`fracreghet.reset` is used to test the specification of fractional
regression models estimated by GMMx or LINx.

## Usage

``` r
fracreghet.reset(object, lastpower.vec = 3, version = "Wald", table = T, ...)
```

## Arguments

- object:

  an object containing the results of an `fracreghet` command.

- lastpower.vec:

  a numeric vector containing the maximum powers of the linear
  predictors to be used in RESET tests.

- version:

  a vector containing the test versions to use. Available options:
  `Wald` (the default) and `LM` (only available for `GMMx`).

- table:

  a logical value indicating whether a summary table with the test
  results should be printed.

- ...:

  Arguments to pass to [nlminb](https://rdrr.io/r/stats/nlminb.html),
  which is used to estimate the model under the alternative hypothesis
  when `version` is equal to `"Wald"` and the null model was estimated
  by `GMMx`.

## Details

`fracreghet.reset` applies the RESET test statistic to fractional
regression models estimated via `fracreghet` using the options `GMMx` or
`LINx`. `fracreghet.reset` may be used to test simultaneously the
validity of the link specification and the transformation applied to the
response variable by each estimator. It is taken into account the option
that was chosen for computing standard errors in the model under
evaluation. See Ramalho and Ramalho (2017) for details.

## Value

`fracreghet.reset` returns a named vector with the test results.

## References

Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of
nonlinear regression models with boundary outcomes and endogeneity, with
applications to nonnegative and fractional responses", *Econometric
Reviews*, 36(4), 397-420.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreghet`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.md),
for fitting fractional regression models under unobserved
heterogeneity.  
[`fracreghet.pe`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.pe.md),
for computing partial effects.

## Examples

``` r
N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

Z <- cbind(rnorm(N),rnorm(N),rnorm(N))
dimnames(Z)[[2]] <- c("Z1","Z2","Z3")

y <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))

res <- fracreghet(y,X,type="GMMx",table=FALSE)

#LM and Wald versions of the RESET test, based on 1 or 2 fitted powers of xb
fracreghet.reset(res,2:3,c("Wald","LM"))
#> 
#> *** RESET test ***
#> Fractional logit regression model
#> 
#> H0:  Estimator: GMMx
#> 
#>  Version Statistic p-value 
#>    LM(2)     0.710   0.399 
#>  Wald(2)     0.453   0.501 
#>    LM(3)     2.433   0.296 
#>  Wald(3)     1.315   0.518 
#> 
```
