# Fitting Panel Data Fractional Response Models

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

`fracregpd` computes the GMM estimators proposed in Ramalho, Ramalho and
Coelho (2018) for panel data fractional response models with both
time-variant and time-invariant unobserved heterogeneity and endogeneous
covariates: GMMww, GMMc, GMMbgw, GMMpfe, GMMcre and GMMpre. In addition,
`fracregpd` also computes QMLcre, which was proposed by Papke and
Wooldridge (2008) and Wooldridge (2019). The Correlated Random Effects
(CRE) approach models the unobserved individual-specific effects as a
linear function of the time averages of covariates, allowing consistent
estimation even when unobserved effects are correlated with regressors.
For overidentified models, `fracregpd` calculates Hansen's J statistic.

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

Ramalho, E. A., Ramalho, J. J. S., & Coelho, L. M. S. (2018),
"Exponential Regression of Fractional-Response Fixed-Effects Models with
an Application to Firm Capital Structure", *Journal of Econometric
Methods*, 7(1), 20150019.

Wooldridge, J. M. (2019). Correlated random effects models with
unbalanced panels. *Journal of Econometrics*, 211(1), 137-150.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md),
for fitting standard cross-sectional fractional response models.  
[`fracreghet`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.md),
for fitting cross-sectional fractional response models with unobserved
heterogeneity.

## Examples

``` r
set.seed(123)
# Simulating Panel Data
N <- 100
T_periods <- 5
id <- rep(1:N, each = T_periods)
time <- rep(1:T_periods, times = N)
x_panel <- rnorm(N * T_periods)

# Unobserved individual effect (CRE)
c_i <- rep(rnorm(N), each = T_periods) 
y_panel <- exp(x_panel + c_i) / (1 + exp(x_panel + c_i))

X <- cbind(x_panel = x_panel)

# Endogenous variable and instrument simulation
z_panel <- rnorm(N * T_periods)
u_panel <- 0.5 * z_panel + rnorm(N * T_periods)
var_endog <- 0.8 * z_panel + u_panel
y_endog <- exp(x_panel + 1.2 * var_endog + c_i + u_panel) / 
             (1 + exp(x_panel + 1.2 * var_endog + c_i + u_panel))

X_endog <- cbind(x_panel = x_panel, var_endog = var_endog)
Z_inst <- cbind(x_panel = x_panel, z_panel = z_panel)

# Estimate a Correlated Random Effects (CRE) Model
fracregpd(id=id, time=time, y=y_panel, x=X, type="QMLcre", link="probit")
#> 
#> -------------------------------------------------------------------------------- 
#>                        Fractional probit regression model 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                QMLcre 
#> Exogeneity:                                                                 TRUE 
#> Use first lag of instruments:                                              FALSE 
#> Standard errors:                                                         cluster 
#> -------------------------------------------------------------------------------- 
#> Initial observations:                                                        500 
#> Estimation observations:                                                     500 
#> Initial cross-sectional units:                                               100 
#> Estimation cross-sectional units:                                            100 
#> Initial periods per unit (avg):                                                5 
#> Estimation periods per unit (avg):                                             5 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>                Estimate Std. Error z value Pr(>|z|)    
#> x_panel         0.52902    0.01143  46.300   <2e-16 ***
#> INTERCEPT_mean -0.01246    0.04901  -0.254    0.799    
#> x_panel_mean   -0.17409    0.13491  -1.290    0.197    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 17:25:28 
#> -------------------------------------------------------------------------------- 
#> 

# Exogeneity, no lags, no time dummies, clustered standard errors, GMMbgw estimator
fracregpd(id=id, time=time, y=y_panel, x=X, type="GMMbgw")
#> 
#> -------------------------------------------------------------------------------- 
#>                        Fractional logit regression model 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                GMMbgw 
#> Exogeneity:                                                                 TRUE 
#> Use first lag of instruments:                                              FALSE 
#> Standard errors:                                                         cluster 
#> -------------------------------------------------------------------------------- 
#> Initial observations:                                                        500 
#> Estimation observations:                                                     500 
#> Initial cross-sectional units:                                               100 
#> Estimation cross-sectional units:                                            100 
#> Initial periods per unit (avg):                                                5 
#> Estimation periods per unit (avg):                                             5 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>          Estimate Std. Error   z value Pr(>|z|)    
#> x_panel 1.000e+00  1.431e-16 6.986e+15   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 17:25:28 
#> -------------------------------------------------------------------------------- 
#> 

# Lagged covariates and instruments, robust standard errors, GMMww estimator
fracregpd(id=id, time=time, y=y_panel, x=X, lags=TRUE, type="GMMww", var.type="robust")
#> 
#> -------------------------------------------------------------------------------- 
#>                        Fractional logit regression model 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                 GMMww 
#> Exogeneity:                                                                 TRUE 
#> Use first lag of instruments:                                               TRUE 
#> Standard errors:                                                          robust 
#> -------------------------------------------------------------------------------- 
#> Initial observations:                                                        500 
#> Estimation observations:                                                     400 
#> Initial cross-sectional units:                                               100 
#> Estimation cross-sectional units:                                            100 
#> Initial periods per unit (avg):                                                5 
#> Estimation periods per unit (avg):                                             4 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>          Estimate Std. Error   z value Pr(>|z|)    
#> x_panel 1.000e+00  8.835e-17 1.132e+16   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 17:25:28 
#> -------------------------------------------------------------------------------- 
#> 

# Endogeneity, time dummies, GMMpfe estimator
fracregpd(id=id, time=time, y=y_endog, x=X_endog, z=Z_inst, 
          x.exogenous=FALSE, type="GMMpfe", tdummies=TRUE)
#> 
#> -------------------------------------------------------------------------------- 
#>                        Fractional logit regression model 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                GMMpfe 
#> Exogeneity:                                                                FALSE 
#> Use first lag of instruments:                                              FALSE 
#> Standard errors:                                                         cluster 
#> -------------------------------------------------------------------------------- 
#> Initial observations:                                                        500 
#> Estimation observations:                                                     500 
#> Initial cross-sectional units:                                               100 
#> Estimation cross-sectional units:                                            100 
#> Initial periods per unit (avg):                                                5 
#> Estimation periods per unit (avg):                                             5 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>           Estimate Std. Error z value Pr(>|z|)    
#> x_panel    0.99602    0.03117  31.957   <2e-16 ***
#> var_endog  1.57689    0.02524  62.474   <2e-16 ***
#> time.2    -0.15363    0.09752  -1.575   0.1152    
#> time.3    -0.04402    0.09147  -0.481   0.6303    
#> time.4    -0.10609    0.09379  -1.131   0.2580    
#> time.5    -0.19163    0.08917  -2.149   0.0316 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 17:25:29 
#> -------------------------------------------------------------------------------- 
#> 
```
