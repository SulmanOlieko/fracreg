# RESET Test for Fractional Regression Models

`fracreg.reset` is used to perform the Regression Equation Specification
Error Test (RESET) to check the functional form and specification of
fractional regression models.

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

`fracreg.reset` applies the RESET test statistic to fractional
regression models estimated via `fracreg`. `fracreg.reset` may be used
to test the link specification of: (i) one-part fractional regression
models; (ii) the binary components of two-part and three-part fractional
regression models; and (iii) the fractional components of two-part and
three-part fractional regression models. When the `Wald` version is
implemented, it is taken into account the option that was chosen for
computing standard errors in the model under evaluation. For the `LM`
version, a robust version is computed in cases (i) and (iii) and a
conventional version in case (ii). See Ramalho, Ramalho and Murteira
(2011) for details on the application of the RESET test in the
fractional regression framework.

## Value

`fracreg.reset` returns a named vector with the test results.

## References

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional regression
models", *Journal of Economic Surveys*, 25(1), 19-68.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md),
for fitting fractional regression models.  
[`fracreg.ggoff`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ggoff.md),
for asymptotically equivalent specification tests.  
[`fracreg.ptest`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ptest.md),
for non-nested hypothesis tests.  
[`fracreg.pe`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md),
for computing partial effects.

## Examples

``` r
N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Testing the logit specification of a standard fractional regression model
#using LM and Wald versions of the RESET test, based on 1 or 2 fitted powers of
#the linear predictor
res <- fracreg(y,X,linkfrac="logit",table=FALSE)
fracreg.reset(res,2:3,c("Wald","LM"))
#> 
#> *** RESET test ***
#> 
#> H0:  Fractional logit model
#> 
#>  Version Statistic p-value  
#>    LM(2)     3.389   0.066 *
#>  Wald(2)     3.830   0.050 *
#>    LM(3)     4.222   0.121  
#>  Wald(3)     3.468   0.177  
#> 

#Testing the probit specification of the binary component of a two-part fractional
#regression model using LR-based RESET tests with quadratic and cubic fitted 
#powers of the linear predictor
res <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1,table=FALSE)
fracreg.reset(res,3,"LR")
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> 
#> *** RESET test ***
#> 
#> H0:  Binary probit component of a two-part model
#> 
#>  Version Statistic p-value 
#>    LR(3)     0.047   0.977 
#> 
```
