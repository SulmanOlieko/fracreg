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
also on unobservables, setting the error term to zero.

**Partial Effects under Unobserved Heterogeneity:** When unobserved
heterogeneity or endogeneity is present, calculating partial effects
requires dealing with the unobserved error \\v_i\\. Let the conditional
mean be \\E(y\|x, v) = G(x\beta + \gamma v)\\. - **Conditional on
Observables (Smearing):** The unobserved heterogeneity is integrated out
over its empirical distribution. The average partial effect for a
continuous variable \\x_k\\ is computed as: \$\$PE_k(x) = \frac{1}{N}
\sum\_{i=1}^N g(x\beta + \gamma \hat{v}\_i) \beta_k\$\$ - **Conditional
on Unobservables (Error = 0):** The partial effect is evaluated for an
individual with the mean level of unobserved heterogeneity (\\v = 0\\):
\$\$PE_k(x) = g(x\beta) \beta_k\$\$

For discrete variables, the partial effects are calculated as the
discrete differences evaluated using either the smearing approach or
setting the error term to zero.

For calculating standard errors, it is taken into account the option
that was previously chosen for estimating the model. See Ramalho and
Ramalho (2017) for details on the computation of partial effects for
fractional response models under unobserved heterogeneity.

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
### Empirical 401(k) Examples 
data("fracreg_k401k") 
y <- fracreg_k401k$prate 
X_het <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age)
 
# fracreghet estimators do not allow exact 1s or 0s
y_adj <- y
y_adj[y_adj == 1] <- 0.999

# Artificial instrumental variable for demonstration 
set.seed(42)
Z_emp <- cbind(X_het, z = fracreg_k401k$mrate * rnorm(length(y))) 
res_emp <- fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], type="QMLxv", link="logit", table=FALSE) 
#> Warning: NA/NaN function evaluation
fracreghet.pe(res_emp, which.x="mrate")
#> 
#> 
#> -------------------------------------------------------------------------------- 
#> Average partial effects (conditional only on observables, based on the smearing estimator) 
#> -------------------------------------------------------------------------------- 
#>                        Fractional logit regression model 
#>                                 Estimator: QMLxv 
#> -------------------------------------------------------------------------------- 
#>       Estimate Std. Error z value Pr(>|z|)
#> mrate   0.1081         NA      NA       NA
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-06 03:48:13 
#> -------------------------------------------------------------------------------- 
 
### Simulated Examples

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
#> X1   0.1364     0.0167   8.166 2.22e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-06 03:48:13 
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
#> X1  0.19889    0.02519   7.895 2.89e-15 ***
#> X2  0.25172    0.01495  16.842  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-06 03:48:13 
#> -------------------------------------------------------------------------------- 
#> 
#> Note: covariates evaluated at the following values:
#> 
#> X1 X2 
#>  1 -1 
```
