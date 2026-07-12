# Fitting Panel Data Fractional Response Regressions

`fracregpd` is used to fit panel data regression models when the
dependent variable has a bounded, fractional nature.

## Usage

``` r
fracregpd(
  id,
  time,
  y,
  x,
  z,
  var.endog,
  x.exogenous = T,
  lags,
  start,
  type,
  GMMww.cor = T,
  link = "logit",
  intercept = T,
  table = FALSE,
  variance = T,
  var.type = "cluster",
  tdummies = F,
  bootstrap = F,
  B = 200,
  offset = NULL,
  or = FALSE,
  level = 0.95,
  ...
)
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
  be optimised. Optional.

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
  `table = FALSERUE`.

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

- offset:

  an optional numeric vector containing an offset. It must be of the
  same dimension as the response variable. It specifies that the
  variable should be included in the model with its coefficient
  constrained to 1.

- or:

  a logical value indicating whether to report odds ratios. Only valid
  when the link function is `"logit"`. Defaults to `FALSE`.

- level:

  a numeric value between 0 and 1 indicating the confidence level for
  the confidence intervals. Defaults to `0.95`.

- ...:

  Arguments to pass to [nlminb](https://rdrr.io/r/stats/nlminb.html).

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

If `variance = TRUE` or `table = FALSERUE` and the algorithm converged
successfully, the previous list also contains the following elements:

- p.var:

  a named covariance matrix.

- var.type:

  covariance matrix type.

## Details

`fracregpd` computes the GMM estimators proposed in Ramalho, Ramalho and
Coelho (2018) for panel data fractional response models with both
time-variant and time-invariant unobserved heterogeneity and endogeneous
covariates: GMMww, GMMc, GMMbgw, GMMpfe, GMMcre and GMMpre. In addition,
`fracregpd` also computes QMLcre, which was proposed by Papke and
Wooldridge (2008) and Wooldridge (2019).

**Correlated Random Effects (CRE) - QMLcre:** In panel data, unobserved
individual-specific heterogeneity \\c_i\\ may be correlated with the
covariates \\x\_{it}\\. The CRE approach (Papke and Wooldridge, 2008)
models this dependence by projecting \\c_i\\ onto the time averages of
the strictly exogenous covariates \\\bar{x}\_i\\: \$\$c_i = \psi +
\bar{x}\_i \xi + a_i\$\$ where \\a_i\\ is an error term independent of
\\x_i\\. Assuming \\a_i \| x_i \sim N(0, \sigma_a^2)\\ and a probit
link, integrating out \\a_i\\ yields the "population-averaged" or scaled
conditional mean: \$\$E(y\_{it} \| x_i) = G(x\_{it} \beta_a + \psi_a +
\bar{x}\_i \xi_a)\$\$ where the parameters with subscript \\a\\ are
scaled by \\(1 + \sigma_a^2)^{-1/2}\\. This equation is estimated via
pooled Bernoulli QML.

**Generalised Method of Moments (GMM):** For models where strict
exogeneity fails or the link function is an exponential-type link,
Ramalho et al. (2018) propose GMM estimators based on the following
general moment conditions: \$\$E\[Z\_{it} (H(y\_{it}) -
\exp(x\_{it}\beta + c_i))\] = 0\$\$ where \\H(\cdot)\\ is a
transformation function and \\Z\_{it}\\ is a matrix of valid
instruments. Estimators such as GMMww, GMMc, and GMMbgw use different
transformations to eliminate the unobserved fixed effect \\c_i\\ before
applying GMM.

For overidentified models, `fracregpd` calculates Hansen's J statistic
to test the validity of the overidentifying restrictions.

## Odds Ratios

When `or=TRUE` and the fractional link function (`linkfrac` or `link`)
is `"logit"`, the model additionally computes odds ratios for the
coefficients. Odds Ratios are exponentiated coefficients. The
corresponding standard errors for the odds ratios are calculated using
the Delta method. The confidence intervals for the odds ratios are
calculated using the adjusted standard errors and the specified `level`
(defaulting to 95%). Odds ratios are particularly useful in fractional
logit models as they provide a direct multiplicative interpretation of
the independent variable on the odds of the fractional outcome.

## References

Papke, L. and Wooldridge, J.M. (2008), "Panel data methods for
fractional response variables with an application to test pass rates",
*Journal of Econometrics*, 145(1-2), 121-133.

Ramalho, E. A., Ramalho, J. J. S., & Coelho, L. M. S. (2018),
"Exponential Regression of Fractional-Response Fixed-Effects Models with
an Application to Firm Capital Structure", *Journal of Econometric
Methods*, 7(1), 20150019.

Wooldridge, J. M. (2019). Correlated random effects models with
unbalanced panels. *Journal of Econometrics*, 211(1), 137-150.

## See also

[`fracreg`](https://sulmanolieko.github.io/fracreg/reference/fracreg.md),
for fitting standard cross-sectional fractional response models.  
[`fracreghet`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.md),
for fitting cross-sectional fractional response models with unobserved
heterogeneity.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## Examples

``` r
### Empirical 401(k) Examples
data("fracreg_k401k")
y <- fracreg_k401k$prate
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)

# Artificial panel data structure for demonstration
N_emp <- nrow(X)
id_emp <- rep(1:(N_emp/2), each=2)
time_emp <- rep(1:2, times=N_emp/2)
mod <- fracregpd(id_emp, time_emp, y, X, type="QMLcre", link="probit")
summary(mod)
#> 
#> -------------------------------------------------------------------------------- 
#>            Fractional probit (correlated random effects)  regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                QMLcre 
#> Data type:                                                                 Panel 
#> Exogeneity:                                                                 TRUE 
#> Use first lag of instruments:                                              FALSE 
#> Standard errors:                                                         cluster 
#> -------------------------------------------------------------------------------- 
#> Number of observations:                                                     1534 
#> Number of groups:                                                            767 
#> Obs per group:                                                                 2 
#> Log pseudolikelihood:                                                  -554.0205 
#> Wald chi2(9):                                                          3191.0501 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>       Final (Correlated Random Effects) Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>                  Coefficient Cluster Std.Err.    z value [95% Conf. Interval]
#> mrate              3.915e-01        7.007e-02  5.588e+00  2.542e-01     0.529
#> age                1.386e-02        3.543e-03  3.911e+00  6.912e-03     0.021
#> totemp            -5.250e-06        2.649e-06 -1.982e+00 -1.044e-05     0.000
#> sole               2.378e-01        6.038e-02  3.938e+00  1.194e-01     0.356
#> (Intercept)_mean   6.235e-01        5.885e-02  1.060e+01  5.082e-01     0.739
#> mrate_mean         5.726e-02        6.343e-02  9.027e-01 -6.706e-02     0.182
#> age_mean           1.544e-03        4.233e-03  3.649e-01 -6.751e-03     0.010
#> totemp_mean        1.223e-06        2.946e-06  4.152e-01 -4.551e-06     0.000
#> sole_mean         -7.157e-02        8.424e-02 -8.495e-01 -2.367e-01     0.094
#>                  Pr(>|z|)    
#> mrate            2.30e-08 ***
#> age              9.20e-05 ***
#> totemp             0.0475 *  
#> sole             8.22e-05 ***
#> (Intercept)_mean  < 2e-16 ***
#> mrate_mean         0.3667    
#> age_mean           0.7152    
#> totemp_mean        0.6780    
#> sole_mean          0.3956    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-12 13:14:31 
#> -------------------------------------------------------------------------------- 
#> 

### Simulated Examples

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

# \donttest{
# Estimate a Correlated Random Effects (CRE) Model
mod <- fracregpd(id=id, time=time, y=y_panel, x=X, type="QMLcre", link="probit")
summary(mod)
#> 
#> -------------------------------------------------------------------------------- 
#>            Fractional probit (correlated random effects)  regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                QMLcre 
#> Data type:                                                                 Panel 
#> Exogeneity:                                                                 TRUE 
#> Use first lag of instruments:                                              FALSE 
#> Standard errors:                                                         cluster 
#> -------------------------------------------------------------------------------- 
#> Number of observations:                                                      500 
#> Number of groups:                                                            100 
#> Obs per group:                                                                 5 
#> Log pseudolikelihood:                                                  -313.6436 
#> Wald chi2(3):                                                          2156.6004 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>       Final (Correlated Random Effects) Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>                  Coefficient Cluster Std.Err.  z value [95% Conf. Interval]
#> x_panel              0.52902          0.01143 46.29965    0.50662     0.551
#> (Intercept)_mean    -0.01246          0.04901 -0.25433   -0.10851     0.084
#> x_panel_mean        -0.17409          0.13491 -1.29040   -0.43850     0.090
#>                  Pr(>|z|)    
#> x_panel            <2e-16 ***
#> (Intercept)_mean    0.799    
#> x_panel_mean        0.197    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-12 13:14:31 
#> -------------------------------------------------------------------------------- 
#> 

# Compute Partial Effects
pe_res <- fracregpd.pe(mod)
summary(pe_res)
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                             Average partial effects 
#> -------------------------------------------------------------------------------- 
#>                     Panel Data Fractional probit regression 
#> -------------------------------------------------------------------------------- 
#>         Estimate Std. Error z value Pr(>|z|)    
#> x_panel 0.189526   0.003544   53.48   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-12 13:14:31 
#> -------------------------------------------------------------------------------- 

# Exogeneity, no lags, no time dummies, clustered standard errors, GMMbgw estimator
mod <- fracregpd(id=id, time=time, y=y_panel, x=X, type="GMMbgw")
summary(mod)
#> 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit  regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                GMMbgw 
#> Data type:                                                                 Panel 
#> Exogeneity:                                                                 TRUE 
#> Use first lag of instruments:                                              FALSE 
#> Standard errors:                                                         cluster 
#> -------------------------------------------------------------------------------- 
#> Number of observations:                                                      500 
#> Number of groups:                                                            100 
#> Obs per group:                                                                 5 
#> Wald chi2(1):                                               4.87462988128311e+31 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                             Final GMM bgw estimates 
#> -------------------------------------------------------------------------------- 
#>         Coefficient Cluster Std.Err.   z value [95% Conf. Interval] Pr(>|z|)
#> x_panel   1.000e+00        1.432e-16 6.982e+15  1.000e+00         1   <2e-16
#>            
#> x_panel ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-12 13:14:31 
#> -------------------------------------------------------------------------------- 
#> 

# Estimate the GMMww estimator with odds ratios and 99% confidence intervals
mod <- fracregpd(id=id, time=time, y=y_panel, x=X, type="GMMww", or=TRUE, level=0.99)
summary(mod)
#> 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit  regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                 GMMww 
#> Data type:                                                                 Panel 
#> Exogeneity:                                                                 TRUE 
#> Use first lag of instruments:                                               TRUE 
#> Standard errors:                                                         cluster 
#> -------------------------------------------------------------------------------- 
#> Number of obs (initial):                                                     500 
#> Number of observations:                                                      400 
#> Number of groups (initial):                                                  100 
#> Number of groups:                                                            100 
#> Obs per group:                                                                 4 
#> Wald chi2(1):                                               8.09654990520078e+32 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                              Final GMM ww estimates 
#> -------------------------------------------------------------------------------- 
#>         Odds Ratio Cluster Std.Err.   z value [99% Conf. Interval] Pr(>|z|)    
#> x_panel  2.718e+00        2.597e-16 1.047e+16  2.718e+00     2.718   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-12 13:14:31 
#> -------------------------------------------------------------------------------- 
#> 

# Lagged covariates and instruments, robust standard errors, GMMww estimator
mod <- fracregpd(id=id, time=time, y=y_panel, x=X, lags=TRUE, type="GMMww", var.type="robust")
summary(mod)
#> 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit  regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                 GMMww 
#> Data type:                                                                 Panel 
#> Exogeneity:                                                                 TRUE 
#> Use first lag of instruments:                                               TRUE 
#> Standard errors:                                                          robust 
#> -------------------------------------------------------------------------------- 
#> Number of obs (initial):                                                     500 
#> Number of observations:                                                      400 
#> Number of groups (initial):                                                  100 
#> Number of groups:                                                            100 
#> Obs per group:                                                                 4 
#> Wald chi2(1):                                               1.27815198429371e+32 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                              Final GMM ww estimates 
#> -------------------------------------------------------------------------------- 
#>         Coefficient Robust Std.Err.   z value [95% Conf. Interval] Pr(>|z|)    
#> x_panel   1.000e+00       8.845e-17 1.131e+16  1.000e+00         1   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-12 13:14:31 
#> -------------------------------------------------------------------------------- 
#> 

# Endogeneity, time dummies, GMMpfe estimator
mod <- fracregpd(id=id, time=time, y=y_endog, x=X_endog, z=Z_inst,
                 x.exogenous=FALSE, type="GMMpfe", tdummies=TRUE)
summary(mod)
#> 
#> -------------------------------------------------------------------------------- 
#>               Fractional logit (pooled fixed effects)  regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                GMMpfe 
#> Data type:                                                                 Panel 
#> Exogeneity:                                                                FALSE 
#> Use first lag of instruments:                                              FALSE 
#> Standard errors:                                                         cluster 
#> -------------------------------------------------------------------------------- 
#> Number of observations:                                                      500 
#> Number of groups:                                                            100 
#> Obs per group:                                                                 5 
#> Wald chi2(6):                                                          7156.2911 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>       Final (Pooled Fixed Effects) Generalized Method of Moments estimates 
#> -------------------------------------------------------------------------------- 
#>           Coefficient Cluster Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
#> x_panel       0.99602          0.03117 31.95739    0.93493     1.057   <2e-16
#> var_endog     1.57689          0.02524 62.47420    1.52742     1.626   <2e-16
#> time.2       -0.15363          0.09752 -1.57541   -0.34476     0.038   0.1152
#> time.3       -0.04402          0.09147 -0.48128   -0.22331     0.135   0.6303
#> time.4       -0.10609          0.09379 -1.13110   -0.28992     0.078   0.2580
#> time.5       -0.19163          0.08917 -2.14916   -0.36640    -0.017   0.0316
#>              
#> x_panel   ***
#> var_endog ***
#> time.2       
#> time.3       
#> time.4       
#> time.5    *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-12 13:14:32 
#> -------------------------------------------------------------------------------- 
#> 
# }
```
