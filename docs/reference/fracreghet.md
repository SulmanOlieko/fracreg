# Fitting Fractional Response Models under Unobserved Heterogeneity

`fracreghet` is used to fit fractional response models under unobserved
heterogeneity, i.e. regression models for proportions, percentages or
fractions that suffer from neglected heterogeneity and/or endogeneity
issues.

## Usage

``` r
fracreghet(y, x, z = x, var.endog, start, type = "GMMx", link = "logit", 
           intercept = T, table = T, variance = T, 
           var.type = "robust", var.cluster, adjust = 0, ...)
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
(2017) for fractional response models with unobserved heterogeneity:
GMMx, which allows for neglected heterogeneity but not for endogeneity;
GMMxv, which allows both issues and assumes a linear reduced form for
the endogeneous covariate (or for a transformation of it); and GMMz,
which also allows for both issues but does not require the assumption of
a reduced form for the endogenous covariate. In addition, `fracreghet`
also computes three linearized estimators (LINx, LINxv and LINz) that
have similar features to their GMM counterparts. It also provides a QML
estimator (QMLxv) that addresses endogeneity using a Control Function
(CF) approach, which includes the first-stage reduced-form residuals as
an additional regressor in the main fractional equation, providing a
Hausman-type test for endogeneity. See Ramalho and Ramalho (2017) for
details on each estimator. For overidentified models, `fracreghet`
calculates Hansen's J statistic. For `GMMx` and `LINx`, `fracreghet`
stores the information needed to implement the RESET test
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

Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of
nonlinear regression models with boundary outcomes and endogeneity, with
applications to nonnegative and fractional responses", *Econometric
Reviews*, 36(4), 397-420.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreghet.reset`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.reset.md),
for the RESET test.  
[`fracreghet.pe`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.pe.md),
for computing partial effects.  
[`fracreg`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md),
for fitting standard cross-sectional fractional response models.  
[`fracregpd`](https://SulmanOlieko.github.io/fracreg/reference/fracregpd.md),
for fitting panel data fractional response models.

## Examples

``` r
set.seed(123)
N <- 1000
x1 <- rnorm(N)

# Simulating an endogenous variable (var.endog) and an instrument (z1)
z1 <- rnorm(N)
u <- 0.5 * z1 + rnorm(N)
var.endog <- 0.8 * z1 + u
y_endog <- exp(0.5 * x1 + 1.2 * var.endog + u) / (1 + exp(0.5 * x1 + 1.2 * var.endog + u))

# Avoid exact 0 or 1 boundaries for some estimators
y_endog[y_endog <= 0] <- 0.01
y_endog[y_endog >= 1] <- 0.99

X <- cbind(x1 = x1, var.endog = var.endog)
Z <- cbind(x1 = x1, z1 = z1)

# Exogeneity (assuming var.endog is exogenous for comparison), GMMx estimator
fracreghet(y = y_endog, x = X, type = "GMMx", link = "logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                        Fractional logit regression model 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                  GMMx 
#> Number of observations:                                                     1000 
#> Standard errors:                                                          robust 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>           Estimate Std. Error z value Pr(>|z|)    
#> Constant  0.091323   0.015355   5.947 2.73e-09 ***
#> x1        0.460767   0.015385  29.950  < 2e-16 ***
#> var.endog 1.807916   0.009021 200.418  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 16:51:54 
#> -------------------------------------------------------------------------------- 
#> 

# Endogeneity, GMMz estimator (does not require reduced form for endog)
fracreghet(y = y_endog, x = X, z = Z, type = "GMMz", link = "logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                        Fractional logit regression model 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                  GMMz 
#> Number of observations:                                                     1000 
#> Standard errors:                                                          robust 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>           Estimate Std. Error z value Pr(>|z|)    
#> Constant   0.15277    0.02018   7.569 3.75e-14 ***
#> x1         0.47947    0.02051  23.377  < 2e-16 ***
#> var.endog  1.61252    0.01445 111.618  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 16:51:54 
#> -------------------------------------------------------------------------------- 
#> 

# Endogeneity, GMMxv estimator (assumes linear reduced form for var.endog)
fracreghet(y = y_endog, x = X, z = Z, var.endog = var.endog, type = "GMMxv", link = "logit")
#> Warning: NaNs produced
#> 
#> -------------------------------------------------------------------------------- 
#>                        Fractional logit regression model 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                 GMMxv 
#> Number of observations:                                                     1000 
#> Standard errors:                                                          robust 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>           Estimate Std. Error z value Pr(>|z|)    
#> Constant  -0.01262    0.01859  -0.679    0.497    
#> x1         0.48705    0.01884  25.853   <2e-16 ***
#> var.endog  1.59737    0.01339 119.332   <2e-16 ***
#> vhat       0.60263    0.01339  45.020   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                                  Reduced form: 
#> -------------------------------------------------------------------------------- 
#>             Estimate Std. Error z value Pr(>|z|)    
#> Z_INTERCEPT -0.02093    0.03082  -0.679    0.497    
#> Z_x1        -0.02149    0.03132  -0.686    0.493    
#> Z_z1         1.32751    0.02949  45.020   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 16:51:54 
#> -------------------------------------------------------------------------------- 
#> 

# Endogeneity, QMLxv control function approach
fracreghet(y = y_endog, x = X, z = Z, var.endog = var.endog, type = "QMLxv", link = "logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                        Fractional logit regression model 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                 QMLxv 
#> Number of observations:                                                     1000 
#> Standard errors:                                                          robust 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>           Estimate Std. Error z value Pr(>|z|)    
#> Constant  -0.01262    0.01859  -0.679    0.497    
#> x1         0.48705    0.01884  25.853   <2e-16 ***
#> var.endog  1.59737    0.01339 119.332   <2e-16 ***
#> vhat       0.60263    0.01339  45.020   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                                  Reduced form: 
#> -------------------------------------------------------------------------------- 
#>             Estimate Std. Error z value Pr(>|z|)    
#> Z_INTERCEPT -0.02093    0.03082  -0.679    0.497    
#> Z_x1        -0.02149    0.03132  -0.686    0.493    
#> Z_z1         1.32751    0.02949  45.020   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 16:51:54 
#> -------------------------------------------------------------------------------- 
#> 
```
