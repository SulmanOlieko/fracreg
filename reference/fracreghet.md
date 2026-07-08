# Fitting Fractional Response Regressions under Unobserved Heterogeneity

`fracreghet` is used to fit fractional response models under unobserved
heterogeneity, i.e. regression models for proportions, percentages or
fractions that suffer from neglected heterogeneity and/or endogeneity
issues.

## Usage

``` r
fracreghet(y, x, z = x, var.endog, start, type = "GMMx", link = "logit", 
           intercept = T, table = T, variance = T, 
           var.type = "robust", var.cluster, adjust = 0, offset = NULL, 
           or = FALSE, level = 0.95, ...)
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
  be optimised. Optional.

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

## Details

`fracreghet` computes the GMM estimators proposed in Ramalho and Ramalho
(2017) for fractional response models with unobserved heterogeneity:
GMMx, which allows for neglected heterogeneity but not for endogeneity;
GMMxv, which allows both issues and assumes a linear reduced form for
the endogeneous covariate (or for a transformation of it); and GMMz,
which also allows for both issues but does not require the assumption of
a reduced form for the endogenous covariate. In addition, `fracreghet`
also computes three linearised estimators (LINx, LINxv and LINz) that
have similar features to their GMM counterparts. It also provides a QML
estimator (QMLxv) that addresses endogeneity using a Control Function
(CF) approach, which includes the first-stage reduced-form residuals as
an additional regressor in the main fractional equation, providing a
Hausman-type test for endogeneity.

**Control Function (CF) Approach - QMLxv:** When a continuous regressor
\\y\_{2i}\\ is endogenous, the CF approach (Papke and Wooldridge, 2008;
Terza et al., 2008) uses a two-stage procedure. First, a linear reduced
form is estimated: \$\$y\_{2i} = z_i \pi + v_i\$\$ where \\z_i\\
includes all exogenous variables and external instruments. The residuals
\\\hat{v}\_i\\ are then included in the fractional response model:
\$\$E(y\_{1i} \| z_i, y\_{2i}, v_i) = G(x_i \beta + \gamma
\hat{v}\_i)\$\$ A test of \\H_0: \gamma = 0\\ serves as a robust
Hausman-type test for endogeneity.

**Generalised Method of Moments (GMM):** For estimators like GMMz, which
do not strictly require a linear reduced form, the estimation relies on
population orthogonality conditions between the instruments \\Z_i\\ and
the model residuals: \$\$E\[Z\_{i} (y_i - G(x_i \beta))\] = 0\$\$ or via
specific transformations of the dependent variable to eliminate
unobserved heterogeneity (Ramalho and Ramalho, 2017).

For overidentified models, `fracreghet` calculates Hansen's J statistic.
For `GMMx` and `LINx`, `fracreghet` stores the information needed to
implement the RESET test
([fracreghet.reset](https://sulmanolieko.github.io/fracreg/reference/fracreghet.reset.md)).
For all estimators, `fracreghet` stores the information needed to
calculate partial effects
([fracreghet.pe](https://sulmanolieko.github.io/fracreg/reference/fracreghet.pe.md)).

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

Papke, L. E. and Wooldridge, J. M. (2008), "Panel data methods for
fractional response variables with an application to test pass rates",
*Journal of Econometrics*, 145, 121-133.

Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of
nonlinear regression models with boundary outcomes and endogeneity, with
applications to nonnegative and fractional responses", *Econometric
Reviews*, 36(4), 397-420.

Terza, J. V., Basu, A., and Rathouz, P. J. (2008), "Two-stage residual
inclusion estimation: addressing endogeneity in health econometric
modeling", *Journal of Health Economics*, 27(3), 531-543.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreghet.reset`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.reset.md),
for the RESET test.  
[`fracreghet.pe`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.pe.md),
for computing partial effects.  
[`fracreg`](https://sulmanolieko.github.io/fracreg/reference/fracreg.md),
for fitting standard cross-sectional fractional response models.  
[`fracregpd`](https://sulmanolieko.github.io/fracreg/reference/fracregpd.md),
for fitting panel data fractional response models.

## Examples

``` r
### Empirical 401(k) Examples 
data("fracreg_k401k") 
y <- fracreg_k401k$prate 
X_het <- cbind(mrate = fracreg_k401k$mrate, ltotemp = fracreg_k401k$ltotemp)
 
# fracreghet estimators do not allow exact 1s or 0s
y_adj <- y
y_adj[y_adj == 1] <- 0.999

# Instrument mrate using age

Z_emp <- cbind(age = fracreg_k401k$age, ltotemp = fracreg_k401k$ltotemp) 
fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], type="QMLxv", link="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>         Fractional logit regression with heteroscedasticity/endogeneity 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                 QMLxv 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1534 
#> Standard errors:                                                          robust 
#> Wald chi2(6):                                                          1991.8748 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                   Final Quasi-Maximum Likelihood xv estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)    
#> Constant    -0.10561         0.75224 -0.14040   -1.57997     1.369    0.888    
#> mrate        3.72138         0.68952  5.39703    2.36994     5.073 6.78e-08 ***
#> ltotemp     -0.07009         0.05348 -1.31060   -0.17491     0.035    0.190    
#> vhat        -2.79515         0.70026 -3.99157   -4.16764    -1.423 6.56e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                                  Reduced form: 
#> -------------------------------------------------------------------------------- 
#>             Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
#> Z_INTERCEPT     0.95046         0.09550  9.95211    0.76327     1.138  < 2e-16
#> Z_age           0.01146         0.00231  4.95997    0.00693     0.016 7.05e-07
#> Z_ltotemp      -0.05534         0.01421 -3.89528   -0.08318    -0.027 9.81e-05
#>                
#> Z_INTERCEPT ***
#> Z_age       ***
#> Z_ltotemp   ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:25:12 
#> -------------------------------------------------------------------------------- 
#> 

# Compute the same QMLxv estimator reporting Odds Ratios with 90% confidence intervals
fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], type="QMLxv", 
           link="logit", or=TRUE, level=0.90) 
#> 
#> -------------------------------------------------------------------------------- 
#>         Fractional logit regression with heteroscedasticity/endogeneity 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                 QMLxv 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1534 
#> Standard errors:                                                          robust 
#> Wald chi2(6):                                                       2243425.9766 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                   Final Quasi-Maximum Likelihood xv estimates 
#> -------------------------------------------------------------------------------- 
#>          Odds Ratio Robust Std.Err.  z value [90% Conf. Interval] Pr(>|z|)    
#> Constant    0.89977         0.67684 -0.14040    0.26108     3.101    0.888    
#> mrate      41.32157        28.49225  5.39703   13.29273   128.452 6.78e-08 ***
#> ltotemp     0.93231         0.04986 -1.31060    0.85380     1.018    0.190    
#> vhat        0.06111         0.04279 -3.99157    0.01931     0.193 6.56e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                                  Reduced form: 
#> -------------------------------------------------------------------------------- 
#>             Odds Ratio Robust Std.Err.   z value [90% Conf. Interval] Pr(>|z|)
#> Z_INTERCEPT   2.586889        0.247056  9.952107   2.210829     3.027  < 2e-16
#> Z_age         1.011524        0.002337  4.959966   1.007688     1.015 7.05e-07
#> Z_ltotemp     0.946168        0.013441 -3.895285   0.924315     0.969 9.81e-05
#>                
#> Z_INTERCEPT ***
#> Z_age       ***
#> Z_ltotemp   ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:25:13 
#> -------------------------------------------------------------------------------- 
#> 
 
### Simulated Examples 

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
#>         Fractional logit regression with heteroscedasticity/endogeneity 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                  GMMx 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1000 
#> Standard errors:                                                          robust 
#> Wald chi2(2):                                                         42761.0276 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                               Final GMMx estimates 
#> -------------------------------------------------------------------------------- 
#>           Coefficient Robust Std.Err.   z value [95% Conf. Interval] Pr(>|z|)
#> Constant    9.132e-02       1.536e-02 5.947e+00  6.123e-02     0.121 2.73e-09
#> x1          4.608e-01       1.538e-02 2.995e+01  4.306e-01     0.491  < 2e-16
#> var.endog   1.808e+00       9.021e-03 2.004e+02  1.790e+00     1.826  < 2e-16
#>              
#> Constant  ***
#> x1        ***
#> var.endog ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:25:13 
#> -------------------------------------------------------------------------------- 
#> 

# Endogeneity, GMMz estimator (does not require reduced form for endog)
fracreghet(y = y_endog, x = X, z = Z, type = "GMMz", link = "logit")
#> 
#> -------------------------------------------------------------------------------- 
#>         Fractional logit regression with heteroscedasticity/endogeneity 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                  GMMz 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1000 
#> Standard errors:                                                          robust 
#> Wald chi2(2):                                                         14903.8562 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                               Final GMMz estimates 
#> -------------------------------------------------------------------------------- 
#>           Coefficient Robust Std.Err.   z value [95% Conf. Interval] Pr(>|z|)
#> Constant      0.15277         0.02018   7.56931    0.11321     0.192 3.75e-14
#> x1            0.47947         0.02051  23.37680    0.43927     0.520  < 2e-16
#> var.endog     1.61252         0.01445 111.61802    1.58421     1.641  < 2e-16
#>              
#> Constant  ***
#> x1        ***
#> var.endog ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:25:13 
#> -------------------------------------------------------------------------------- 
#> 

# Endogeneity, GMMxv estimator (assumes linear reduced form for var.endog)
fracreghet(y = y_endog, x = X, z = Z, var.endog = var.endog, type = "GMMxv", link = "logit")
#> Warning: NaNs produced
#> 
#> -------------------------------------------------------------------------------- 
#>         Fractional logit regression with heteroscedasticity/endogeneity 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                 GMMxv 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1000 
#> Standard errors:                                                          robust 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                              Final GMMxv estimates 
#> -------------------------------------------------------------------------------- 
#>           Coefficient Robust Std.Err.   z value [95% Conf. Interval] Pr(>|z|)
#> Constant     -0.01262         0.01859  -0.67857   -0.04905     0.024    0.497
#> x1            0.48705         0.01884  25.85287    0.45012     0.524   <2e-16
#> var.endog     1.59737         0.01339 119.33166    1.57113     1.624   <2e-16
#> vhat          0.60263         0.01339  45.01988    0.57640     0.629   <2e-16
#>              
#> Constant     
#> x1        ***
#> var.endog ***
#> vhat      ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                                  Reduced form: 
#> -------------------------------------------------------------------------------- 
#>             Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
#> Z_INTERCEPT    -0.02093         0.03082 -0.67933   -0.08133     0.039    0.497
#> Z_x1           -0.02149         0.03132 -0.68612   -0.08288     0.040    0.493
#> Z_z1            1.32751         0.02949 45.01988    1.26971     1.385   <2e-16
#>                
#> Z_INTERCEPT    
#> Z_x1           
#> Z_z1        ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:25:13 
#> -------------------------------------------------------------------------------- 
#> 

# Endogeneity, QMLxv control function approach
fracreghet(y = y_endog, x = X, z = Z, var.endog = var.endog, type = "QMLxv", link = "logit")
#> 
#> -------------------------------------------------------------------------------- 
#>         Fractional logit regression with heteroscedasticity/endogeneity 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                 QMLxv 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1000 
#> Standard errors:                                                          robust 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                   Final Quasi-Maximum Likelihood xv estimates 
#> -------------------------------------------------------------------------------- 
#>           Coefficient Robust Std.Err.   z value [95% Conf. Interval] Pr(>|z|)
#> Constant     -0.01262         0.01859  -0.67857   -0.04905     0.024    0.497
#> x1            0.48705         0.01884  25.85287    0.45012     0.524   <2e-16
#> var.endog     1.59737         0.01339 119.33166    1.57113     1.624   <2e-16
#> vhat          0.60263         0.01339  45.01988    0.57640     0.629   <2e-16
#>              
#> Constant     
#> x1        ***
#> var.endog ***
#> vhat      ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                                  Reduced form: 
#> -------------------------------------------------------------------------------- 
#>             Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
#> Z_INTERCEPT    -0.02093         0.03082 -0.67933   -0.08133     0.039    0.497
#> Z_x1           -0.02149         0.03132 -0.68612   -0.08288     0.040    0.493
#> Z_z1            1.32751         0.02949 45.01988    1.26971     1.385   <2e-16
#>                
#> Z_INTERCEPT    
#> Z_x1           
#> Z_z1        ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:25:13 
#> -------------------------------------------------------------------------------- 
#> 
```
