# Fractional Response Models under Unobserved Heterogeneity - Partial Effects

`fracreghet.pe` is used to compute average and/or conditional partial
effects in fractional response models under unobserved heterogeneity.

## Usage

``` r
fracreghet.pe(object, smearing = T, APE = T, CPE = F, at = NULL, 
              which.x = NULL, table = T, variance = T)
```

## Arguments

- object:

  an object containing the results of an `fracreghet` command.

- smearing:

  a logical value indicating whether the smearing correction is to be
  applied

- APE:

  a logical value indicating whether average partial effects are to be
  computed.

- CPE:

  a logical value indicating whether conditional partial effects are to
  be computed.

- at:

  a numeric vector containing the covariates' values at which the
  conditional partial effects are to be computed or the strings `"mean"`
  (the default) or `"median"`, in which cases the covariates are
  evaluated at their mean or median values (or mode, in case of dummy
  variables), respectively.

- which.x:

  a vector containing the names of the covariates to which the partial
  effects are to be computed.

- table:

  a logical value indicating whether a summary table with the results
  should be printed.

- variance:

  a logical value indicating whether the variance of the estimated
  partial effects should be calculated. Defaults to `TRUE` whenever
  `table = TRUE`.

## Details

`fracreghet.pe` calculates partial effects for fractional response
models estimated via `fracreghet`. `fracreghet.pe` may be used to
compute average or conditional partial effects. These partial effects
may be conditional only on observables, using the smearing estimator, or
also on unobservables, setting the error term to zero. For calculating
standard errors, it is taken into account the option that was previously
chosen for estimating the model. See Ramalho and Ramalho (2017) for
details on the computation of partial effects for fractional response
models under unobserved heterogeneity.

## Value

`fracreghet.pe` returns a list with the following element:

- PE.p:

  a named vector of partial effects.

If `variance = TRUE` or `table = TRUE`, the previous list also contains
the following element:

- PE.sd:

  a named vector of standard errors of the estimated partial effects.

When both average and conditional partial effects are requested, two
lists containing the previous elements are returned, indexed by the
prefixes `ape` and `cpe`.

## References

Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of
nonlinear regression models with boundary outcomes and endogeneity, with
applications to nonnegative and fractional responses", *Econometric
Reviews*, 36(4), 397-420.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreghet`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.md),
for fitting fractional response models under unobserved heterogeneity.  
[`fracreghet.reset`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.reset.md),
for the RESET test.  

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

#Smearing estimator of average partial effects for variable X1
fracreghet.pe(res,which.x="X1")
#> 
#> 
#> -------------------------------------------------------------------------------- 
#> Average partial effects (conditional only on observables, based on the smearing estimator) 
#> -------------------------------------------------------------------------------- 
#>                        Fractional logit regression model 
#>                                 Estimator: GMMx 
#> -------------------------------------------------------------------------------- 
#>    Estimate Std. Error z value Pr(>|z|)    
#> X1  0.16698    0.00957   17.45   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 17:25:27 
#> -------------------------------------------------------------------------------- 

#Naive estimator of conditional partial effects for all covariates,
#which are evaluated at X1=1 and X2=-1
fracreghet.pe(res,smearing=FALSE,APE=FALSE,CPE=TRUE,at=c(1,-1))
#> 
#> 
#> -------------------------------------------------------------------------------- 
#> Conditional partial effects (conditional on both observables and unobservables, with error term = 0) 
#> -------------------------------------------------------------------------------- 
#>                        Fractional logit regression model 
#>                                 Estimator: GMMx 
#> -------------------------------------------------------------------------------- 
#>    Estimate Std. Error z value Pr(>|z|)    
#> X1  0.22878    0.01322   17.31   <2e-16 ***
#> X2  0.20655    0.01447   14.28   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 17:25:27 
#> -------------------------------------------------------------------------------- 
#> 
#> Note: covariates evaluated at the following values:
#> 
#> X1 X2 
#>  1 -1 
```
