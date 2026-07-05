# Fitting Fractional Regression Models

`fracreg` is used to fit fractional regression models, which are
appropriate for responses that are proportions, percentages, or
fractions restricted to the \[0, 1\] interval. It supports standard
one-part models, two-part hurdle models for modeling boundary values at
0 or 1, and three-part models for double inflation at both 0 and 1.

## Usage

``` r
fracreg(y, x, x2 = x, linkbin, linkfrac, type = "1P", inflation = 0, 
        intercept = TRUE, table = TRUE, variance = TRUE, var.type = "default", 
        var.eim = TRUE, var.cluster, dfc = FALSE, ...)
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
  two-part fractional regression model, or a vector of two link
  functions for the two binary components of a three-part model (e.g.
  `c("logit", "probit")`). Available options: `logit`, `probit`,
  `cauchit`, `loglog`, `cloglog`.

- linkfrac:

  a description of the link function to use in standard fractional
  regression models or in the fractional component of a two-part
  fractional regression model. Available options: `logit`, `probit`,
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
  fractional regression models.

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
  quasi-maximum likelihood, such as standard fractional regression
  models or the fractional component of a two-part fractional regression
  model), `cluster` (recommended in the case of panel data) and
  `default` (implements the `standard` or `robust` versions as
  appropriate).

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

- ...:

  Arguments to pass to [glm](https://rdrr.io/r/stats/glm.html).

## Details

`fracreg` estimates one-part, two-part hurdle, and three-part
double-inflated fractional regression models; see Ramalho, Ramalho and
Murteira (2011) and Fang and Ma (2013) for details on those models. The
one-part models and the fractional component of two- and three-part
models are estimated by Bernoulli-based quasi-maximum likelihood, while
the binary components of two- and three-part models are estimated by
maximum likelihood. `fracreg` uses the standard
[glm](https://rdrr.io/r/stats/glm.html) command to perform the
estimations. Therefore, `fracreg` is essentially a convenience command,
allowing estimation of several alternative fractional regression models
using the same command. In addition, `fracreg` provides an R-squared
measure for all models (calculated as the square of the correlation
coefficient between the actual and fitted values of the dependent
variable), calculates the fitted values of the dependent variable in
two-part models and stores the information needed to implement some very
useful commands for fractional regression models:
[fracreg.reset](https://SulmanOlieko.github.io/fracreg/reference/fracreg.reset.md)
(RESET test),
[fracreg.ptest](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ptest.md)
(P test),
[fracreg.ggoff](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ggoff.md)
(GGOFF tests) and
[fracreg.pe](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md)
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

## References

Papke, L. E. and Wooldridge, J. M. (1996), "Econometric methods for
fractional response variables with an application to 401(k) plan
participation rates", *Journal of Applied Econometrics*, 11(6), 619-632.

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional regression
models", *Journal of Economic Surveys*, 25(1), 19-68.

Fang, K., & Ma, S. (2013), "Three-part model for fractional response
variables with application to Chinese household health insurance
coverage", *Journal of Applied Statistics*, 40(5), 925-940.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg.reset`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.reset.md)
and
[`fracreg.ggoff`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ggoff.md),
for specification tests.  
[`fracreg.ptest`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ptest.md),
for non-nested hypothesis tests.  
[`fracreg.pe`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md),
for computing partial effects.  
`fracreghet`, for fitting cross-sectional fractional regression models
with unobserved heterogeneity.  
`fracregpd`, for fitting panel data fractional regression models.

## Examples

``` r
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

# fracreg estimation of a logit fractional regression model
fracreg(y, X, type="1P", linkfrac="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                Fractional regression model - logit specification 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Number of observations:                                                     1000 
#> Pseudo R-squared:                                                         0.3903 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>          Estimate Std. Error z value Pr(>|z|)    
#> Constant -0.56969    0.06750   -8.44   <2e-16 ***
#> x1        0.78822    0.04015   19.63   <2e-16 ***
#> x2        1.35611    0.11899   11.40   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 03:17:37 
#> -------------------------------------------------------------------------------- 
#> 

# fracreg estimation of the binary logit component of the two-part fractional
# regression model with y=0 as the relevant boundary value
fracreg(y, X, type="2Pbin", inflation=0, linkbin="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>            Binary component of a two-part model - logit specification 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Number of observations:                                                     1000 
#> Pseudo R-squared:                                                        0.14527 
#> Log-Likelihood:                                                        -295.2068 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>          Estimate Std. Error z value Pr(>|z|)    
#> Constant   1.4638     0.1933   7.573 3.64e-14 ***
#> x1         1.1857     0.1287   9.213  < 2e-16 ***
#> x2         2.2279     0.3932   5.666 1.46e-08 ***
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
#>         Fractional component of a two-part model - probit specification 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Number of observations:                                                      881 
#> Pseudo R-squared:                                                        0.32474 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>          Estimate Std. Error z value Pr(>|z|)    
#> Constant -0.10334    0.03566  -2.898  0.00376 ** 
#> x1        0.37512    0.02129  17.623  < 2e-16 ***
#> x2        0.61326    0.06529   9.392  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 

# fracreg estimation of both components of a two-part fractional regression model
# with y=0 as the relevant boundary value and using a cloglog binary link
# function and a logit fractional link function
fracreg(y, X, type="2P", inflation=0, linkbin="cloglog", linkfrac="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>           Binary component of a two-part model - cloglog specification 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Number of observations:                                                     1000 
#> Pseudo R-squared:                                                        0.14837 
#> Log-Likelihood:                                                        -291.7986 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>          Estimate Std. Error z value Pr(>|z|)    
#> Constant  0.43428    0.08825   4.921 8.60e-07 ***
#> x1        0.55192    0.06009   9.184  < 2e-16 ***
#> x2        1.05359    0.17403   6.054 1.41e-09 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>          Fractional component of a two-part model - logit specification 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Number of observations:                                                      881 
#> Pseudo R-squared:                                                        0.32412 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>          Estimate Std. Error z value Pr(>|z|)    
#> Constant -0.17323    0.05786  -2.994  0.00275 ** 
#> x1        0.61205    0.03554  17.221  < 2e-16 ***
#> x2        1.00509    0.10703   9.391  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                Two-part model - binary cloglog + fractional logit 
#> -------------------------------------------------------------------------------- 
#> Pseudo R-squared:                                                        0.38829 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 03:17:37 
#> -------------------------------------------------------------------------------- 
#> 

# Three-part double-inflated model (y has both 0s and 1s)
fracreg(y, X, type="3P", linkbin=c("logit","probit"), linkfrac="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>        First binary component of a three-part model - logit specification 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Number of observations:                                                     1000 
#> Pseudo R-squared:                                                        0.14527 
#> Log-Likelihood:                                                        -295.2068 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>          Estimate Std. Error z value Pr(>|z|)    
#> Constant   1.4638     0.1933   7.573 3.64e-14 ***
#> x1         1.1857     0.1287   9.213  < 2e-16 ***
#> x2         2.2279     0.3932   5.666 1.46e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>       Second binary component of a three-part model - probit specification 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                    ML 
#> Number of observations:                                                      881 
#> Pseudo R-squared:                                                        0.18634 
#> Log-Likelihood:                                                        -348.8644 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>          Estimate Std. Error z value Pr(>|z|)    
#> Constant -1.71537    0.12933  -13.26  < 2e-16 ***
#> x1        0.64810    0.06323   10.25  < 2e-16 ***
#> x2        1.07752    0.19174    5.62 1.91e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>         Fractional component of a three-part model - logit specification 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Number of observations:                                                      715 
#> Pseudo R-squared:                                                        0.24292 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                                 Final estimates 
#> -------------------------------------------------------------------------------- 
#>          Estimate Std. Error z value Pr(>|z|)    
#> Constant -0.26763    0.04539  -5.896 3.73e-09 ***
#> x1        0.36198    0.02454  14.751  < 2e-16 ***
#> x2        0.58505    0.07969   7.342 2.11e-13 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>        Three-part model - binary logit , binary probit + fractional logit 
#> -------------------------------------------------------------------------------- 
#> Pseudo R-squared:                                                        0.38917 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 03:17:37 
#> -------------------------------------------------------------------------------- 
#> 
```
