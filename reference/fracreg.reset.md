# RESET Test for Fractional Response Regressions

`fracreg.reset` is used to perform the Regression Equation Specification
Error Test (RESET) to check the functional form and specification of
fractional response models.

## Usage

``` r
fracreg.reset(object, lastpower.vec = 3, version = "LM", table = FALSE, ...)
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

## Value

`fracreg.reset` returns a named vector with the test results.

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

## References

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional response
models", *Journal of Economic Surveys*, 25(1), 19-68.

Ramsey, J.B. (1969), "Tests for Specification Errors in Classical Linear
Least-Squares Regression Analysis", *Journal of the Royal Statistical
Society: Series B (Methodological)*, 31(2), 350-371.

## See also

[`fracreg`](https://sulmanolieko.github.io/fracreg/reference/fracreg.md),
for fitting fractional response models.  
[`fracreg.ggoff`](https://sulmanolieko.github.io/fracreg/reference/fracreg.ggoff.md),
for asymptotically equivalent specification tests.  
[`fracreg.ptest`](https://sulmanolieko.github.io/fracreg/reference/fracreg.ptest.md),
for non-nested hypothesis tests.  
[`fracreg.pe`](https://sulmanolieko.github.io/fracreg/reference/fracreg.pe.md),
for computing partial effects.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## Examples

``` r
### Empirical 401(k) Examples
data("fracreg_k401k")
y <- fracreg_k401k$prate
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)

m <- fracreg(y, X, type="1P", linkfrac="logit")
reset_res <- fracreg.reset(m)
summary(reset_res)
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
#>                          Run Date: 2026-07-12 14:57:29 
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
mod <- fracreg(y,X,linkfrac="logit")
reset_res <- fracreg.reset(mod,2:3,c("Wald","LM"))
summary(reset_res)
#> 
#> -------------------------------------------------------------------------------- 
#>                                    RESET test 
#> -------------------------------------------------------------------------------- 
#> H0: Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#>         Statistic p-value
#> LM(2)       0.134   0.715
#> Wald(2)     0.133   0.716
#> LM(3)       0.135   0.935
#> Wald(3)     0.133   0.936
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-12 14:57:29 
#> -------------------------------------------------------------------------------- 

#Testing the probit specification of the binary component of a two-part fractional
#regression model using LR-based RESET tests with quadratic and cubic fitted 
#powers of the linear predictor
mod <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1)
reset_res <- fracreg.reset(mod,3,"LR")
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
summary(reset_res)
#> 
#> -------------------------------------------------------------------------------- 
#>                                    RESET test 
#> -------------------------------------------------------------------------------- 
#> H0: Binary probit component of a two-part regression 
#> -------------------------------------------------------------------------------- 
#>       Statistic p-value
#> LR(3)     4.211   0.122
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-12 14:57:29 
#> -------------------------------------------------------------------------------- 
```
