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
N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#fracreg estimation of a logit fractional regression model
fracreg(y,X,linkfrac="logit")
#> 
#> *** Fractional logit regression model ***
#> 
#>           Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT 0.086260   0.065448   1.318    0.188    
#> X1        0.836328   0.073667  11.353    0.000 ***
#> X2        0.862703   0.075733  11.391    0.000 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 250 
#> R-squared: 0.551 
#> 
#> 

#fracreg estimation of the binary logit component of the two-part fractional
#regression model with y=1 as the relevant boundary value
fracreg(y,X,linkbin="logit",type="2Pbin",inf=1)
#> 
#> *** Binary component of a two-part model - logit specification ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -3.191181   0.400029  -7.977    0.000 ***
#> X1         1.551070   0.293556   5.284    0.000 ***
#> X2         1.748795   0.331834   5.270    0.000 ***
#> 
#> Number of observations: 250 
#> R-squared: 0.389 
#> 
#> 

#fracreg estimation of the fractional component of the two-part fractional
#regression model with y=1 as the relevant boundary value and using a
#probit link function
fracreg(y,X,linkfrac="probit",type="2Pfrac",inf=1)
#> 
#> *** Fractional component of a two-part model - probit specification ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -0.056614   0.037488  -1.510    0.131    
#> X1         0.403554   0.040399   9.989    0.000 ***
#> X2         0.416087   0.044353   9.381    0.000 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 220 
#> R-squared: 0.434 
#> 
#> 

#fracreg estimation of both components of a two-part fractional regression model
#with y=1 as the relevant boundary value and using a cloglog binary link
#function and a logit fractional link function
fracreg(y,X,linkbin="cloglog",linkfrac="logit",type="2P",inf=1)
#> 
#> *** Binary component of a two-part model - cloglog specification ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -3.095296   0.349232  -8.863    0.000 ***
#> X1         1.247048   0.223632   5.576    0.000 ***
#> X2         1.447821   0.245933   5.887    0.000 ***
#> 
#> Number of observations: 250 
#> R-squared: 0.391 
#> 
#> 
#> 
#> *** Fractional component of a two-part model - logit specification ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -0.089035   0.061037  -1.459    0.145    
#> X1         0.662234   0.068483   9.670    0.000 ***
#> X2         0.681284   0.074740   9.115    0.000 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 220 
#> R-squared: 0.434 
#> 
#> 
#> 
#> *** Two-part model - binary cloglog + fractional logit  ***
#> 
#> R-squared: 0.293 
#> 
#> 

#Three-part model (y has both 0s and 1s)
y3p <- y
y3p[1:20] <- 0
y3p[21:40] <- 1
fracreg(y3p,X,linkbin=c("logit","probit"),linkfrac="logit",type="3P")
#> 
#> *** First binary component of a three-part model - logit specification ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT  2.469764   0.240642  10.263    0.000 ***
#> X1         0.218908   0.229740   0.953    0.341    
#> X2        -0.271501   0.245830  -1.104    0.269    
#> 
#> Number of observations: 250 
#> R-squared: 0.009 
#> 
#> 
#> 
#> *** Second binary component of a three-part model - probit specification ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -0.997159   0.115979  -8.598    0.000 ***
#> X1         0.601432   0.113444   5.302    0.000 ***
#> X2         0.646366   0.126233   5.120    0.000 ***
#> 
#> Number of observations: 230 
#> R-squared: 0.305 
#> 
#> 
#> 
#> *** Fractional component of a three-part model - logit specification ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -0.087891   0.070593  -1.245    0.213    
#> X1         0.638664   0.078939   8.091    0.000 ***
#> X2         0.651209   0.088764   7.336    0.000 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 183 
#> R-squared: 0.389 
#> 
#> 
#> 
#> *** Three-part model - binary logit , binary probit + fractional logit  ***
#> 
#> R-squared: 0.353 
#> 
#> 
```
