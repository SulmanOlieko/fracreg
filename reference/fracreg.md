# Fitting Fractional Response Regressions

`fracreg` is used to fit fractional response models, which are
appropriate for responses that are proportions, percentages, or
fractions restricted to the \[0, 1\] interval. It supports standard
one-part models, two-part hurdle models for modelling boundary values at
0 or 1, and three-part models for double inflation at both 0 and 1.

## Usage

``` r
fracreg(y, x, x2 = x, linkbin, linkfrac, type = "1P", inflation = 0, 
        intercept = TRUE, table = TRUE, variance = TRUE, var.type = "default", 
        var.eim = TRUE, var.cluster, dfc = FALSE, offset = NULL, 
        or = FALSE, level = 0.95, ...)
```

## Arguments

- y:

  a numeric vector containing the values of the response variable.

- x:

  a numeric matrix, with column names, containing the values of the
  covariates.

- x2:

  a numeric matrix, with column names, containing the values of the
  covariates in the fractional component of two-part models if option
  `type = "2P"` is defined. Defaults to `x`.

- linkbin:

  a description of the link function to use in the binary component of a
  two-part fractional response model, or a vector of two link functions
  for the two binary components of a three-part model (e.g.
  `c("logit", "probit")`). Available options: `logit`, `probit`,
  `cauchit`, `loglog`, `cloglog`.

- linkfrac:

  a description of the link function to use in standard fractional
  response models or in the fractional component of a two-part
  fractional response model. Available options: `logit`, `probit`,
  `cauchit`, `loglog`, `cloglog`.

- type:

  a description of the model to estimate: a standard one-part model
  (`1P`, the default), a two-part model (`2P`), the binary component of
  a two-part model (`2Pbin`), the fractional component of a two-part
  model (`2Pfrac`), or a three-part model (`3P`) for double boundary
  inflation.

- inflation:

  a numeric value indicating which of the extreme values of `0` (the
  default) or `1` is the relevant boundary value for defining two-part
  fractional response models.

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
  be calculated. Options are `standard` (recommended for models
  estimated by maximum likelihood, such as the binary component of
  two-part models), `robust` (recommended for models estimated by
  quasi-maximum likelihood, such as standard fractional response models
  or the fractional component of a two-part fractional response model),
  `cluster` (recommended in the case of panel data) and `default`
  (implements the `standard` or `robust` versions as appropriate).

- var.eim:

  a logical value indicating whether the expected information matrix
  should be used in the calculation of the variance. When false, the
  observation information matrix will be used. Defaults to `TRUE`.

- var.cluster:

  a numeric vector containing the values of the variable that specifies
  to which cluster each observation belongs.

- dfc:

  a logical value indicating whether a degrees of freedom correction
  should be applied to the covariance matrix. Defaults to `FALSE`.

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

  Arguments to pass to [glm](https://rdrr.io/r/stats/glm.html).

## Details

`fracreg` estimates one-part, two-part hurdle, and three-part
double-inflated fractional response models; see Ramalho, Ramalho and
Murteira (2011) and Fang and Ma (2013) for details on those models.

**One-Part Fractional Response Regressions (`type = "1P"`):** The
standard one-part model assumes that the conditional expectation of the
fractional response \\y_i \in \[0,1\]\\ is given by: \$\$E(y_i\|x_i) =
G(x_i \beta)\$\$ where \\G(\cdot)\\ is a known non-linear link function
mapping the linear predictor to the unit interval (e.g., logit, probit).
The parameters \\\beta\\ are estimated by maximising the Bernoulli-based
quasi-log-likelihood function: \$\$\ln L_i(\beta) = y_i \ln\[G(x_i
\beta)\] + (1 - y_i) \ln\[1 - G(x_i \beta)\]\$\$ This estimator requires
only the correct specification of the conditional mean to yield
consistent parameter estimates (Papke and Wooldridge, 1996).

**Two-Part Hurdle Models (`type = "2P"`):** When the data exhibits a
boundary mass (e.g., at \\y_i = 0\\), the two-part hurdle model handles
the boundary values separately from the interior fractional values. Let
\\y_i^\*\\ be a binary indicator such that \\y_i^\* = 1\\ if \\y_i \>
0\\ and \\y_i^\* = 0\\ otherwise. The probability of observing a
boundary value is modelled as: \$\$P(y_i = 0 \| x\_{1i}) = 1 - F(x\_{1i}
\gamma_1)\$\$ \$\$P(y_i \> 0 \| x\_{1i}) = F(x\_{1i} \gamma_1)\$\$ where
\\F(\cdot)\\ is a binary link function. Conditional on observing an
interior fractional value, the response is modelled as: \$\$E(y_i \|
x\_{2i}, y_i \> 0) = G(x\_{2i} \beta_2)\$\$ The unconditional mean of
the response is therefore: \$\$E(y_i\|x_i) = F(x\_{1i} \gamma_1) \times
G(x\_{2i} \beta_2)\$\$

**Three-Part Double Inflated Models (`type = "3P"`):** For data
containing boundary mass at both \\0\\ and \\1\\, the three-part model
estimates two separate binary mechanisms for each boundary and a
fractional component for the interior values \\(0, 1)\\, extending the
two-part logic to double inflation (Fang and Ma, 2013).

`fracreg` uses the standard [glm](https://rdrr.io/r/stats/glm.html)
command to perform the estimations. Therefore, `fracreg` is essentially
a convenience command, allowing estimation of several alternative
fractional response models using the same command. In addition,
`fracreg` provides an R-squared measure for all models (calculated as
the square of the correlation coefficient between the actual and fitted
values of the dependent variable), calculates the fitted values of the
dependent variable in two-part models and stores the information needed
to implement some very useful commands for fractional response models:
[fracreg.reset](https://sulmanolieko.github.io/fracreg/reference/fracreg.reset.md)
(RESET test),
[fracreg.ptest](https://sulmanolieko.github.io/fracreg/reference/fracreg.ptest.md)
(P test),
[fracreg.ggoff](https://sulmanolieko.github.io/fracreg/reference/fracreg.ggoff.md)
(GGOFF tests) and
[fracreg.pe](https://sulmanolieko.github.io/fracreg/reference/fracreg.pe.md)
(partial effects).

## Value

When `type = "1P" or "2Pfrac"`, `fracreg` returns a list with the
following elements:

- class:

  "fracreg".

- formula:

  the model formula.

- type:

  the name of the estimated model.

- link:

  the name of the specified link.

- method:

  estimation method. Currently, "QML" (quasi-maximum likelihood) for
  fractional components or models and"ML" (maximum likelihood) for the
  binary component of two-part models.

- p:

  a named vector of coefficients.

- yhat:

  the fitted mean values.

- xbhat:

  the fitted mean values of the linear predictor.

- converged:

  logical. Was the algorithm judged to have converged?

- x.names:

  a vector containing the names of the covariates.

If `variance = TRUE` or `table = TRUE`, the previous list also contains
the following elements:

- p.var:

  a named covariance matrix.

- var.type:

  covariance matrix type.

- var.eim:

  logical. Was the expected information matrix used in the computation
  of the covariance matrix?

- dfc:

  logical. Was a degrees of freedom correction used for the computation
  of the covariance matrix?

If `var.type = "cluster"`, the list also contains the following element:

- var.cluster:

  the variable that specifies to which cluster each observation belongs.

When `type = "2Pbin"`, `fracreg` returns a similar list with the
following additional element:

- LL:

  the value of the log-likelihood.

When `type = "2P"`, `fracreg` returns the previous lists, indexed by the
prefixes `resBIN` and `resFRAC`, and the following additional elements:

- class:

  "fracreg".

- type:

  "2P".

- ybase:

  a numeric vector containing the values of the response variable.

- x2base:

  a numeric matrix containing the values of the covariates.

- yhat2P:

  the overall fitted mean values.

- converged:

  logical. Were the algorithms judged to have converged in both parts of
  the model?

When `type = "3P"`, `fracreg` returns the previous lists, indexed by the
prefixes `resBIN0`, `resBIN1`, and `resFRAC`, and the following
additional elements:

- class:

  "fracreg".

- type:

  "3P".

- ybase:

  a numeric vector containing the values of the response variable.

- x2base:

  a numeric matrix containing the values of the covariates.

- yhat3P:

  the overall fitted mean values.

- converged:

  logical. Were the algorithms judged to have converged in all parts of
  the model?

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

Papke, L. E. and Wooldridge, J. M. (1996), "Econometric methods for
fractional response variables with an application to 401(k) plan
participation rates", *Journal of Applied Econometrics*, 11(6), 619-632.

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional response
models", *Journal of Economic Surveys*, 25(1), 19-68.

Fang, K., & Ma, S. (2013), "Three-part model for fractional response
variables with application to Chinese household health insurance
coverage", *Journal of Applied Statistics*, 40(5), 925-940.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg.reset`](https://sulmanolieko.github.io/fracreg/reference/fracreg.reset.md)
and
[`fracreg.ggoff`](https://sulmanolieko.github.io/fracreg/reference/fracreg.ggoff.md),
for specification tests.  
[`fracreg.ptest`](https://sulmanolieko.github.io/fracreg/reference/fracreg.ptest.md),
for non-nested hypothesis tests.  
[`fracreg.pe`](https://sulmanolieko.github.io/fracreg/reference/fracreg.pe.md),
for computing partial effects.  
`fracreghet`, for fitting cross-sectional fractional response models
with unobserved heterogeneity.  
`fracregpd`, for fitting panel data fractional response models.

## Examples

``` r
### Empirical 401(k) Examples
data("fracreg_k401k")
y <- fracreg_k401k$prate
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)

# 1P Model
fracreg(y, X, type="1P", linkfrac="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1534 
#> Pseudo R-squared:                                                        0.14667 
#> Log pseudolikelihood:                                                  -553.1626 
#> Wald chi2(4):                                                           147.3049 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                     Final Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient Robust Std.Err.    z value [95% Conf. Interval] Pr(>|z|)
#> Constant   9.316e-01       8.408e-02  1.108e+01  7.668e-01     1.096  < 2e-16
#> mrate      9.531e-01       1.371e-01  6.951e+00  6.843e-01     1.222 3.62e-12
#> age        2.791e-02       4.877e-03  5.723e+00  1.835e-02     0.037 1.05e-08
#> totemp    -8.182e-06       3.061e-06 -2.673e+00 -1.418e-05     0.000  0.00751
#> sole       3.405e-01       8.066e-02  4.222e+00  1.824e-01     0.499 2.43e-05
#>             
#> Constant ***
#> mrate    ***
#> age      ***
#> totemp   ** 
#> sole     ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:49:44 
#> -------------------------------------------------------------------------------- 
#> 

# 1P Model reporting odds ratios and 99% confidence intervals
fracreg(y, X, type="1P", linkfrac="logit", or=TRUE, level=0.99)
#> 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1534 
#> Pseudo R-squared:                                                        0.14667 
#> Log pseudolikelihood:                                                  -553.1626 
#> Wald chi2(4):                                                           147.3049 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                     Final Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Odds Ratio Robust Std.Err.    z value [99% Conf. Interval] Pr(>|z|)
#> Constant  2.539e+00       2.134e-01  1.108e+01  2.044e+00     3.153  < 2e-16
#> mrate     2.594e+00       3.556e-01  6.951e+00  1.822e+00     3.692 3.62e-12
#> age       1.028e+00       5.015e-03  5.723e+00  1.015e+00     1.041 1.05e-08
#> totemp    1.000e+00       3.061e-06 -2.673e+00  1.000e+00     1.000  0.00751
#> sole      1.406e+00       1.134e-01  4.222e+00  1.142e+00     1.730 2.43e-05
#>             
#> Constant ***
#> mrate    ***
#> age      ***
#> totemp   ** 
#> sole     ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:49:44 
#> -------------------------------------------------------------------------------- 
#> 

# 2P Model (modelling mass at 1)
fracreg(y, X, type="2P", inflation=1, linkbin="logit", linkfrac="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                         Part 1: Binary logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1534 
#> Pseudo R-squared:                                                         0.1485 
#> Log-likelihood:                                                        -938.1759 
#> Wald chi2(4):                                                           173.5169 
#> Prob > chi2:                                                              0.0000 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                        Final Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient EIM Std.Err.    z value [95% Conf. Interval] Pr(>|z|)    
#> Constant  -1.396e+00    1.270e-01 -1.099e+01 -1.645e+00    -1.147   <2e-16 ***
#> mrate      9.053e-01    9.699e-02  9.334e+00  7.152e-01     1.095   <2e-16 ***
#> age        1.156e-02    6.218e-03  1.858e+00 -6.312e-04     0.024   0.0631 .  
#> totemp    -1.418e-05    6.324e-06 -2.242e+00 -2.657e-05     0.000   0.0249 *  
#> sole       8.651e-01    1.131e-01  7.651e+00  6.435e-01     1.087    2e-14 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                       Part 2: Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                      852 
#> Pseudo R-squared:                                                        0.10004 
#> Log pseudolikelihood:                                                  -450.8391 
#> Wald chi2(4):                                                            65.4063 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                     Final Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient Robust Std.Err.    z value [95% Conf. Interval] Pr(>|z|)
#> Constant   7.460e-01       6.850e-02  1.089e+01  6.118e-01     0.880  < 2e-16
#> mrate      3.877e-01       9.725e-02  3.987e+00  1.971e-01     0.578 6.69e-05
#> age        2.562e-02       4.010e-03  6.390e+00  1.777e-02     0.033 1.66e-10
#> totemp    -4.061e-06       3.073e-06 -1.322e+00 -1.008e-05     0.000    0.186
#> sole      -1.510e-02       6.556e-02 -2.303e-01 -1.436e-01     0.113    0.818
#>             
#> Constant ***
#> mrate    ***
#> age      ***
#> totemp      
#> sole        
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>              Two-part regression - binary logit + fractional logit 
#> -------------------------------------------------------------------------------- 
#> Data type:                                                       Cross-sectional 
#> Pseudo R-squared:                                                        0.11243 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:49:44 
#> -------------------------------------------------------------------------------- 
#> 

# 3P Model (inject artificial 0s for demonstration)
y_3p <- y; y_3p[1:50] <- 0
fracreg(y_3p, X, type="3P", linkbin=c("logit","logit"), linkfrac="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                         Part 1: Binary logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1534 
#> Pseudo R-squared:                                                        0.00324 
#> Log-likelihood:                                                        -216.6222 
#> Wald chi2(4):                                                             3.7679 
#> Prob > chi2:                                                              0.4383 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                        Final Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient EIM Std.Err.    z value [95% Conf. Interval] Pr(>|z|)    
#> Constant   2.867e+00    3.157e-01  9.080e+00  2.248e+00     3.486   <2e-16 ***
#> mrate      1.147e-01    2.131e-01  5.381e-01 -3.030e-01     0.532    0.590    
#> age        8.375e-03    1.765e-02  4.745e-01 -2.622e-02     0.043    0.635    
#> totemp     1.036e-04    6.959e-05  1.489e+00 -3.275e-05     0.000    0.136    
#> sole       3.371e-01    2.983e-01  1.130e+00 -2.476e-01     0.922    0.259    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                         Part 2: Binary logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1484 
#> Pseudo R-squared:                                                        0.15199 
#> Log-likelihood:                                                        -903.5457 
#> Wald chi2(4):                                                            172.201 
#> Prob > chi2:                                                              0.0000 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                        Final Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient EIM Std.Err.    z value [95% Conf. Interval] Pr(>|z|)    
#> Constant  -1.435e+00    1.306e-01 -1.099e+01 -1.691e+00    -1.179  < 2e-16 ***
#> mrate      9.244e-01    9.877e-02  9.360e+00  7.309e-01     1.118  < 2e-16 ***
#> age        1.198e-02    6.339e-03  1.890e+00 -4.453e-04     0.024   0.0588 .  
#> totemp    -1.371e-05    6.317e-06 -2.170e+00 -2.609e-05     0.000   0.0300 *  
#> sole       8.852e-01    1.154e-01  7.674e+00  6.591e-01     1.111 1.67e-14 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                       Part 3: Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                      826 
#> Pseudo R-squared:                                                        0.09937 
#> Log pseudolikelihood:                                                  -437.3715 
#> Wald chi2(4):                                                            62.3685 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                     Final Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient Robust Std.Err.    z value [95% Conf. Interval] Pr(>|z|)
#> Constant   7.388e-01       7.039e-02  1.050e+01  6.008e-01     0.877  < 2e-16
#> mrate      3.960e-01       1.019e-01  3.885e+00  1.962e-01     0.596 0.000102
#> age        2.531e-02       4.068e-03  6.223e+00  1.734e-02     0.033 4.88e-10
#> totemp    -3.817e-06       3.088e-06 -1.236e+00 -9.869e-06     0.000 0.216472
#> sole      -4.483e-03       6.672e-02 -6.719e-02 -1.353e-01     0.126 0.946434
#>             
#> Constant ***
#> mrate    ***
#> age      ***
#> totemp      
#> sole        
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>      Three-part regression - binary logit , binary logit + fractional logit 
#> -------------------------------------------------------------------------------- 
#> Data type:                                                       Cross-sectional 
#> Pseudo R-squared:                                                        0.07934 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:49:44 
#> -------------------------------------------------------------------------------- 
#> 

### Simulated Examples

set.seed(123)
N <- 1000
x1 <- rnorm(N)
x2 <- runif(N)

# Generating a fractional dependent variable with inflation at 0 and 1
XB <- -0.5 + 0.8 * x1 + 1.2 * x2 + rnorm(N)
y_latent <- exp(XB) / (1 + exp(XB))

y <- y_latent
# Inflate at boundaries
y[y_latent < 0.2] <- 0
y[y_latent > 0.8] <- 1

X <- cbind(x1 = x1, x2 = x2)

# fracreg estimation of a logit fractional response model
fracreg(y, X, type="1P", linkfrac="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1000 
#> Pseudo R-squared:                                                         0.3903 
#> Log pseudolikelihood:                                                   -614.973 
#> Wald chi2(2):                                                           472.2453 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                     Final Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)    
#> Constant    -0.56969         0.06750 -8.43960   -0.70199    -0.437   <2e-16 ***
#> x1           0.78822         0.04015 19.63424    0.70954     0.867   <2e-16 ***
#> x2           1.35611         0.11899 11.39668    1.12289     1.589   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:49:44 
#> -------------------------------------------------------------------------------- 
#> 

# fracreg estimation of the binary logit component of the two-part fractional
# regression model with y=0 as the relevant boundary value
fracreg(y, X, type="2Pbin", inflation=0, linkbin="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                         Part 1: Binary logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1000 
#> Pseudo R-squared:                                                        0.14527 
#> Log-likelihood:                                                        -295.2068 
#> Wald chi2(2):                                                           102.5171 
#> Prob > chi2:                                                              0.0000 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                        Final Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient EIM Std.Err. z value [95% Conf. Interval] Pr(>|z|)    
#> Constant      1.4638       0.1933  7.5732     1.0849     1.843 3.64e-14 ***
#> x1            1.1857       0.1287  9.2126     0.9335     1.438  < 2e-16 ***
#> x2            2.2279       0.3932  5.6662     1.4572     2.999 1.46e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 

# fracreg estimation of the fractional component of the two-part fractional
# regression model with y=0 as the relevant boundary value and using a
# probit link function
fracreg(y, X, type="2Pfrac", inflation=0, linkfrac="probit")
#> 
#> -------------------------------------------------------------------------------- 
#>                       Part 2: Fractional probit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                      881 
#> Pseudo R-squared:                                                        0.32474 
#> Log pseudolikelihood:                                                  -555.3304 
#> Wald chi2(2):                                                           391.6944 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                     Final Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)    
#> Constant    -0.10334         0.03566 -2.89793   -0.17323    -0.033  0.00376 ** 
#> x1           0.37512         0.02129 17.62293    0.33340     0.417  < 2e-16 ***
#> x2           0.61326         0.06529  9.39237    0.48529     0.741  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 

# fracreg estimation of both components of a two-part fractional response model
# with y=0 as the relevant boundary value and using a cloglog binary link
# function and a logit fractional link function
fracreg(y, X, type="2P", inflation=0, linkbin="cloglog", linkfrac="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                        Part 1: Binary cloglog regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1000 
#> Pseudo R-squared:                                                        0.14837 
#> Log-likelihood:                                                        -291.7986 
#> Wald chi2(2):                                                           100.5683 
#> Prob > chi2:                                                              0.0000 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                        Final Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient EIM Std.Err. z value [95% Conf. Interval] Pr(>|z|)    
#> Constant     0.43428      0.08825 4.92113    0.26132     0.607 8.60e-07 ***
#> x1           0.55192      0.06009 9.18429    0.43414     0.670  < 2e-16 ***
#> x2           1.05359      0.17403 6.05409    0.71250     1.395 1.41e-09 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                       Part 2: Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                      881 
#> Pseudo R-squared:                                                        0.32412 
#> Log pseudolikelihood:                                                  -555.4465 
#> Wald chi2(2):                                                           368.0384 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                     Final Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)    
#> Constant    -0.17323         0.05786 -2.99384   -0.28664    -0.060  0.00275 ** 
#> x1           0.61205         0.03554 17.22078    0.54239     0.682  < 2e-16 ***
#> x2           1.00509         0.10703  9.39074    0.79531     1.215  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>             Two-part regression - binary cloglog + fractional logit 
#> -------------------------------------------------------------------------------- 
#> Data type:                                                       Cross-sectional 
#> Pseudo R-squared:                                                        0.38829 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:49:44 
#> -------------------------------------------------------------------------------- 
#> 

# Three-part double-inflated model (y has both 0s and 1s)
fracreg(y, X, type="3P", linkbin=c("logit","probit"), linkfrac="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                         Part 1: Binary logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1000 
#> Pseudo R-squared:                                                        0.14527 
#> Log-likelihood:                                                        -295.2068 
#> Wald chi2(2):                                                           102.5171 
#> Prob > chi2:                                                              0.0000 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                        Final Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient EIM Std.Err. z value [95% Conf. Interval] Pr(>|z|)    
#> Constant      1.4638       0.1933  7.5732     1.0849     1.843 3.64e-14 ***
#> x1            1.1857       0.1287  9.2126     0.9335     1.438  < 2e-16 ***
#> x2            2.2279       0.3932  5.6662     1.4572     2.999 1.46e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                         Part 2: Binary probit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                      881 
#> Pseudo R-squared:                                                        0.18634 
#> Log-likelihood:                                                        -348.8644 
#> Wald chi2(2):                                                           123.5731 
#> Prob > chi2:                                                              0.0000 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                        Final Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient EIM Std.Err.   z value [95% Conf. Interval] Pr(>|z|)    
#> Constant    -1.71537      0.12933 -13.26354   -1.96885    -1.462  < 2e-16 ***
#> x1           0.64810      0.06323  10.25024    0.52418     0.772  < 2e-16 ***
#> x2           1.07752      0.19174   5.61978    0.70172     1.453 1.91e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                       Part 3: Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                      715 
#> Pseudo R-squared:                                                        0.24292 
#> Log pseudolikelihood:                                                  -484.9654 
#> Wald chi2(2):                                                            261.436 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                     Final Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)    
#> Constant    -0.26763         0.04539 -5.89578   -0.35660    -0.179 3.73e-09 ***
#> x1           0.36198         0.02454 14.75144    0.31389     0.410  < 2e-16 ***
#> x2           0.58505         0.07969  7.34179    0.42886     0.741 2.11e-13 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>     Three-part regression - binary logit , binary probit + fractional logit 
#> -------------------------------------------------------------------------------- 
#> Data type:                                                       Cross-sectional 
#> Pseudo R-squared:                                                        0.38917 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-08 20:49:44 
#> -------------------------------------------------------------------------------- 
#> 
```
