# Fitting Fractional Regression Models under Unobserved Heterogeneity

`fracreghet` is used to fit fractional regression models under
unobserved heterogeneity, i.e. regression models for proportions,
percentages or fractions that suffer from neglected heterogeneity and/or
endogeneity issues.

## Usage

``` r
fracreghet(y, x, z = x, var.endog, start, type = "GMMx", link = "logit", intercept = T,
  table = T, variance = T, var.type = "robust", var.cluster, adjust = 0, ...)
```

## Arguments

- y:

  a numeric vector containing the values of the response variable.

- x:

  a numeric matrix, with column names, containing the values of all
  covariates (exogenous and endogenous).

- z:

  a numeric matrix, with column names, containing the values of all
  exogenous variables (covariates and instrumental variables). Defaults
  to `x`.

- var.endog:

  a numeric vector containing the values of the endogenous covariate (or
  of some transformation of it), which will be used as dependent
  variable in the linear reduced form assumed for application of xv-type
  estimators.

- start:

  a numeric vector containing the initial values for the parameters to
  be optimized. Optional.

- type:

  a description of the estimator to compute: `GMMx` (the default),
  `GMMxv`, `GMMz`, `LINx`, `LINxv`, `LINz` or `QMLxv`.

- link:

  a description of the link function to use. Available options for all
  estimators: `logit` and `cloglog`. Additional available options for
  QML and LIN estimators: `probit`, `cauchit` and `loglog`.

- intercept:

  a logical value indicating whether the model should include a constant
  term or not.

- table:

  a logical value indicating whether a summary table with the regression
  results should be printed.

- variance:

  a logical value indicating whether the variance of the estimated
  parameters should be calculated. Defaults to `TRUE` whenever
  `table = TRUE`.

- var.type:

  a description of the type of variance of the estimated parameters to
  be calculated. Options are `robust`, the default, and `cluster`.

- var.cluster:

  a numeric vector containing the values of the variable that specifies
  to which cluster each observation belongs.

- adjust:

  the numeric value to be added to the response variable in case of
  boundary observations when the LIN estimators are applied or the
  string `drop`, which implies that the boundary observations are
  dropped.

- ...:

  Arguments to pass to [nlminb](https://rdrr.io/r/stats/nlminb.html).

## Details

`fracreghet` computes the GMM estimators proposed in Ramalho and Ramalho
(2016) for fractional regression models with unobserved heterogeneity:
GMMx, which allows for neglected heterogeneity but not for endogeneity;
GMMxv, which allows both issues and assumes a linear reduced form for
the endogeneous covariate (or for a transformation of it); and GMMz,
which also allows for both issues but does not require the assumption of
a reduced form for the endogenous covariate. In addition, `fracreghet`
also computes three linearized estimators (LINx, LINxv and LINz) that
have similar features to their GMM counterparts as well as a QML
estimator that allows for endogeneity but not for neglected
heterogeneity (QMLxv); see Ramalho and Ramalho (2016) for details on
each estimator. For overidentified models, `fracreghet` calculates
Hansen's J statistic. For `GMMx` and `LINx`, `fracreghet` stores the
information needed to implement the RESET test
([fracreghet.reset](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.reset.md)).
For all estimators, `fracreghet` stores the information needed to
calculate partial effects
([fracreghet.pe](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.pe.md)).

## Value

`fracreghet` returns a list with the following elements:

- class:

  "fracreghet".

- formula:

  the model formula.

- type:

  the name of the estimator computed.

- link:

  the name of the specified link.

- adjust:

  The value or the type of the adjustment applied to LIN estimators.

- p:

  a named vector of coefficients.

- Hy:

  the transformed values of the response variable when GMM or LIN
  estimators are computed or the values of the response variable in the
  QML case.

- xbhat:

  the fitted mean values of the linear predictor (for xv-type
  estimators, includes the term relative to the first-stage residual).

- converged:

  logical. Was the algorithm judged to have converged?

- x.names:

  a vector containing the names of the covariates.

In case of an overidentifying model, the following element is also
returned:

- J:

  the result of Hansen's J test of overidentifying moment conditions.

If `variance = TRUE` or `table = TRUE` and the algorithm converged
successfully, the previous list also contains the following elements:

- p.var:

  a named covariance matrix.

- var.type:

  covariance matrix type.

If `var.type = "cluster"`, the list also contains the following element:

- var.cluster:

  the variable that specifies to which cluster each observation belongs.

## References

Ramalho, E.A. and J.J.S. Ramalho (2016), "Moment-based estimation of
nonlinear regression models with boundary outcomes and endogeneity, with
applications to nonnegative and fractional responses", *Econometric
Reviews*, forthcoming (DOI: 10.1080/07474938.2014.976531).

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreghet.reset`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.reset.md),
for the RESET test.  
[`fracreghet.pe`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.pe.md),
for computing partial effects.  
[`fracreg`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md),
for fitting standard cross-sectional fractional regression models.  
[`fracregpd`](https://SulmanOlieko.github.io/fracreg/reference/fracregpd.md),
for fitting panel data fractional regression models.

## Examples

``` r
N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

Z <- cbind(rnorm(N),rnorm(N),rnorm(N))
dimnames(Z)[[2]] <- c("Z1","Z2","Z3")

y <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))

#Exogeneity, GMMx estimator
fracreghet(y,X,type="GMMx")
#> 
#> *** Fractional logit regression model ***
#> *** Estimator: GMMx
#> 
#>           Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT 0.383211   0.071031   5.395    0.000 ***
#> X1        0.915521   0.077771  11.772    0.000 ***
#> X2        1.080971   0.061189  17.666    0.000 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 250 
#> 
#> 

#Endogeneity, GMMz estimator
fracreghet(y,X,Z,type="GMMz")
#> 
#> *** Fractional logit regression model ***
#> *** Estimator: GMMz
#> 
#>           Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT 0.384353   0.092977   4.134    0.000 ***
#> X1        1.002693   0.423179   2.369    0.018 ** 
#> X2        1.006577   0.499748   2.014    0.044 ** 
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 250 
#> 
#> J test of overidentifying moment conditions: 0.1297948 (p-value: 0.7186449 )

#Endogeneity, GMMxv estimator
fracreghet(y,X,Z,X[,1],type="GMMxv")
#> 
#> *** Fractional logit regression model ***
#> *** Estimator: GMMxv
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT  0.361853   0.110813   3.265    0.001 ***
#> X1         1.133831   0.848603   1.336    0.182    
#> X2         1.083086   0.060107  18.019    0.000 ***
#> vhat      -0.221761   0.852539  -0.260    0.795    
#> 
#> Reduced form:
#>              Estimate Std. Error t value Pr(>|t|)    
#> Z_INTERCEPT  0.093652   0.058979   1.588    0.112    
#> Z_Z1        -0.007289   0.053477  -0.136    0.892    
#> Z_Z2         0.072459   0.051399   1.410    0.159    
#> Z_Z3        -0.042023   0.066847  -0.629    0.530    
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 250 
#> 
#> 
```
