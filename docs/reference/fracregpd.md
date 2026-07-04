# Fitting Panel Data Fractional Regression Models

`fracregpd` is used to fit panel data regression models when the
dependent variable has a bounded, fractional nature.

## Usage

``` r
fracregpd(id, time, y, x, z, var.endog, x.exogenous = T, lags, start, type,
  GMMww.cor = T, link = "logit", intercept = T, table = T, variance = T,
  var.type = "cluster", tdummies = F, bootstrap = F, B = 200, ...)
```

## Arguments

- id:

  a numeric vector identifying the cross-sectional units.

- time:

  a numeric vector identifying the time periods in which the
  cross-sectional units were observed.

- y:

  a numeric vector containing the values of the response variable.

- x:

  a numeric matrix, with column names, containing the values of all
  covariates (exogenous and endogenous).

- z:

  a numeric matrix, with column names, containing the values of all
  exogenous variables (covariates and external instrumental variables).
  Only required in case of endogenous explanatory variables.

- var.endog:

  a numeric vector containing the values of the endogenous covariate (or
  of some transformation of it), which will be used as dependent
  variable in the linear reduced form assumed for application of the
  `QMLcre` estimator. Only required for this estimator.

- x.exogenous:

  a logical value indicating whether all explanatory variables are
  assumed to be exogenous or not.

- lags:

  a logical value indicating whether the first lags of `x` or `z` should
  be used as instruments for `x`. Defaults to `TRUE` for the GMMww and
  GMMc estimators and to `FALSE` for the remaining estimators. The
  `GMMcre` and `QMLcre` estimators do not admit lagged instruments.

- start:

  a numeric vector containing the initial values for the parameters to
  be optimized. Optional.

- type:

  a description of the estimator to compute: `GMMww`, `GMMc`, `GMMbgw`,
  `GMMpfe`, `GMMcre`, `GMMpre` or `QMLcre`.

- GMMww.cor:

  a logical value indicating whether each explanatory variable should be
  transformed in deviations from its overall mean before computing the
  `GMMww` estimator.

- link:

  a description of the link function to use. Available options for all
  GMM estimators: `logit` and `cloglog`. Only option for the `QMLcre`
  estimator: `probit`.

- intercept:

  a logical value indicating whether the model should include a constant
  term or not. Only relevant for the `GMMpre` estimator.

- table:

  a logical value indicating whether a summary table with the regression
  results should be printed.

- variance:

  a logical value indicating whether the variance of the estimated
  parameters should be calculated. Defaults to `TRUE` whenever
  `table = TRUE`.

- var.type:

  a description of the type of variance of the estimated parameters to
  be calculated. Options are `cluster`, the default, and `robust`. In
  overidentified models, it also affects the parameter estimates via the
  GMM weighting matrix.

- tdummies:

  a logical value indicating whether time dummies should be included
  among the model explanatory variables.

- bootstrap:

  a logical value indicating whether bootstrap should be used in the
  estimation of the parameter standard errors.

- B:

  the number of bootstrap replications.

- ...:

  Arguments to pass to [nlminb](https://rdrr.io/r/stats/nlminb.html).

## Details

`fracregpd` computes the GMM estimators proposed in Ramalho and Ramalho
(2016) for panel data fractional regression models with both
time-variant and time-invariant unobserved heterogeneity and endogeneous
covariates: GMMww, GMMc, GMMbgw, GMMpfe, GMMcre and GMMpre. In addition,
`fracregpd` also computes QMLcre, which was proposed by Papke and
Wooldridge (2008) and Wooldridge (2010). For overidentified models,
`fracregpd` calculates Hansen's J statistic.

## Value

`fracregpd` returns a list with the following elements:

- type:

  the name of the estimator computed.

- link:

  the name of the specified link.

- p:

  a named vector of coefficients.

- Hy:

  the transformed values of the response variable when GMM estimators
  are computed or the values of the response variable in the QML case.

- converged:

  logical. Was the algorithm judged to have converged?

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

## References

Papke, L. and Wooldridge, J.M. (2008), "Panel data methods for
fractional response variables with an application to test pass rates",
*Journal of Econometrics*, 145(1-2), 121-233.

Ramalho, E.A. and J.J.S. Ramalho (2016), "Exponential regression of
fractional-response fixed-efects models with an application to firm
capital structure", mimeo.

Wooldridge, J.M. (2010), "Correlated random effects models with
unbalanced panels", mimeo.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md),
for fitting standard cross-sectional fractional regression models.  
[`fracreghet`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.md),
for fitting cross-sectional fractional regression models with unobserved
heterogeneity.

## Examples

``` r
id <- rep(1:50,each=5)
time <- rep(1:5,50)
NT <- 250

XBu <- rnorm(NT)
y <- exp(XBu)/(1+exp(XBu))

X <- cbind(rnorm(NT),rnorm(NT))
dimnames(X)[[2]] <- c("X1","X2")

Z <- cbind(rnorm(NT),rnorm(NT),rnorm(NT))
dimnames(Z)[[2]] <- c("Z1","Z2","Z3")

# Exogeneity, no lags, no time dummies, clustered standard errors, GMMbgw estimator
fracregpd(id,time,y,X,type="GMMbgw")
#> 
#> *** Fractional logit regression model ***
#> *** Estimator: GMMbgw
#> *** Exogeneity: TRUE
#> *** Use first lag of instruments: FALSE
#> 
#>     Estimate Std. Error t value Pr(>|t|)    
#> X1 -0.033483   0.077734  -0.431    0.667    
#> X2 -0.201427   0.072537  -2.777    0.005 ***
#> 
#> Note: cluster standard errors
#> 
#> Number of observations (initial): 250 
#> Number of observations (for estimation): 250 
#> Number of cross-sectional units (initial): 50 
#> Number of cross-sectional units (for estimation): 50 
#> Average number of time periods per cross-sectional unit (initial): 5 
#> Average number of time periods per cross-sectional unit (for estimation): 5 
#> 
#> 

# Lagged covariates and instruments, robust standard errors, GMMww estimator
fracregpd(id,time,y,X,lags=TRUE,type="GMMww",var.type="robust")
#> 
#> *** Fractional logit regression model ***
#> *** Estimator: GMMww
#> *** Exogeneity: TRUE
#> *** Use first lag of instruments: TRUE
#> 
#>     Estimate Std. Error t value Pr(>|t|)    
#> X1 -0.125826   0.088486  -1.422    0.155    
#> X2 -0.202961   0.071003  -2.858    0.004 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations (initial): 250 
#> Number of observations (for estimation): 200 
#> Number of cross-sectional units (initial): 50 
#> Number of cross-sectional units (for estimation): 50 
#> Average number of time periods per cross-sectional unit (initial): 5 
#> Average number of time periods per cross-sectional unit (for estimation): 4 
#> 
#> 

# Endogeneity, time dummies, GMMpfe estimator
fracregpd(id,time,y,X,Z,x.exogenous=FALSE,type="GMMpfe",tdummies=TRUE)
#> 
#> *** Fractional logit regression model ***
#> *** Estimator: GMMpfe
#> *** Exogeneity: FALSE
#> *** Use first lag of instruments: FALSE
#> 
#>         Estimate Std. Error t value Pr(>|t|) 
#> X1     -0.124274   0.435104  -0.286    0.775 
#> X2     -0.563168   0.744593  -0.756    0.449 
#> time.2  0.254132   0.401075   0.634    0.526 
#> time.3 -0.093777   0.299719  -0.313    0.754 
#> time.4 -0.266803   0.304761  -0.875    0.381 
#> time.5 -0.294691   0.329488  -0.894    0.371 
#> 
#> Note: cluster standard errors
#> 
#> Number of observations (initial): 250 
#> Number of observations (for estimation): 250 
#> Number of cross-sectional units (initial): 50 
#> Number of cross-sectional units (for estimation): 50 
#> Average number of time periods per cross-sectional unit (initial): 5 
#> Average number of time periods per cross-sectional unit (for estimation): 5 
#> 
#> J test of overidentifying moment conditions: 0.3391103 (p-value: 0.5603432 )

# Standard errors based on 100 bootstrap samples, QMLcre estimator (not run)
#fracregpd(id,time,y,X,Z,X[,1],x.exogenous=FALSE,type="QMLcre",link="probit",bootstrap=TRUE,B=100)
```
