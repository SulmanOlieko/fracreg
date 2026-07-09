# RESET Test for Fractional Response Regressions

`fracreg.reset` is used to perform the Regression Equation Specification
Error Test (RESET) to check the functional form and specification of
fractional response models.

## Usage

``` r
fracreg.reset(object, lastpower.vec = 3, version = "LM", table = TRUE, ...)
```

## Arguments

- object:

  an object containing the results of an `fracreg` command.

- lastpower.vec:

  a numeric vector containing the maximum powers of the linear
  predictors to be used in RESET tests.

- version:

  a vector containing the test versions to use. Available options:
  `Wald`, `LM` (the default) and, only for the binary component of
  two-part models, `LR`. More than one option may be chosen.

- table:

  a logical value indicating whether a summary table with the test
  results should be printed.

- ...:

  Arguments to pass to [glm](https://rdrr.io/r/stats/glm.html), which is
  used to estimate the model under the alternative hypothesis when
  `version` is a vector containing `"Wald"` or `"LR"`.

## Details

`fracreg.reset` applies the RESET test statistic to fractional response
models estimated via `fracreg`. `fracreg.reset` may be used to test the
link specification of: (i) one-part fractional response models; (ii) the
binary components of two-part and three-part fractional response models;
and (iii) the fractional components of two-part and three-part
fractional response models.

**RESET Test Framework:** The Regression Equation Specification Error
Test (RESET) assesses whether the link function \\G(\cdot)\\ and the
linear index \\x\beta\\ are correctly specified. It tests the null
hypothesis \\H_0: \gamma = 0\\ in the augmented model: \$\$E(y\|x) =
G(x\beta + \sum\_{k=2}^P \gamma_k (x\hat{\beta})^k)\$\$ where \\P\\ is
the maximum power of the linear predictor (specified by `lastpower.vec`)
and \\\hat{\beta}\\ are the estimated parameters from the baseline
model.

When the `Wald` version is implemented, it is taken into account the
option that was chosen for computing standard errors in the model under
evaluation. For the `LM` version, a robust version is computed in cases
(i) and (iii) and a conventional version in case (ii). See Ramalho,
Ramalho and Murteira (2011) for details on the application of the RESET
test in the fractional response framework.

## Value

`fracreg.reset` returns a named vector with the test results.

## References

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional response
models", *Journal of Economic Surveys*, 25(1), 19-68.

Ramsey, J.B. (1969), "Tests for Specification Errors in Classical Linear
Least-Squares Regression Analysis", *Journal of the Royal Statistical
Society: Series B (Methodological)*, 31(2), 350-371.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg`](https://sulmanolieko.github.io/fracreg/reference/fracreg.md),
for fitting fractional response models.  
[`fracreg.ggoff`](https://sulmanolieko.github.io/fracreg/reference/fracreg.ggoff.md),
for asymptotically equivalent specification tests.  
[`fracreg.ptest`](https://sulmanolieko.github.io/fracreg/reference/fracreg.ptest.md),
for non-nested hypothesis tests.  
[`fracreg.pe`](https://sulmanolieko.github.io/fracreg/reference/fracreg.pe.md),
for computing partial effects.

## Examples

``` r
### Empirical 401(k) Examples
data("fracreg_k401k")
y <- fracreg_k401k$prate
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)

m <- fracreg(y, X, type="1P", linkfrac="logit")
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
#>                          Run Date: 2026-07-09 20:37:03 
#> -------------------------------------------------------------------------------- 
#> 
fracreg.reset(m)
#> 
#> -------------------------------------------------------------------------------- 
#>                                    RESET test 
#> -------------------------------------------------------------------------------- 
#> H0: Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#>       Statistic p-value   
#> LM(3)     10.29 0.00583 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-09 20:37:03 
#> -------------------------------------------------------------------------------- 

### Simulated Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Testing the logit specification of a standard fractional response model
#using LM and Wald versions of the RESET test, based on 1 or 2 fitted powers of
#the linear predictor
res <- fracreg(y,X,linkfrac="logit",table=FALSE)
fracreg.reset(res,2:3,c("Wald","LM"))
#> 
#> -------------------------------------------------------------------------------- 
#>                                    RESET test 
#> -------------------------------------------------------------------------------- 
#> H0: Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#>         Statistic p-value
#> LM(2)       0.013   0.909
#> Wald(2)     0.013   0.909
#> LM(3)       0.042   0.979
#> Wald(3)     0.045   0.978
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-09 20:37:03 
#> -------------------------------------------------------------------------------- 

#Testing the probit specification of the binary component of a two-part fractional
#regression model using LR-based RESET tests with quadratic and cubic fitted 
#powers of the linear predictor
res <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1,table=FALSE)
fracreg.reset(res,3,"LR")
#> 
#> -------------------------------------------------------------------------------- 
#>                                    RESET test 
#> -------------------------------------------------------------------------------- 
#> H0: Binary probit component of a two-part regression 
#> -------------------------------------------------------------------------------- 
#>       Statistic p-value
#> LR(3)     0.528   0.768
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-09 20:37:03 
#> -------------------------------------------------------------------------------- 
```
