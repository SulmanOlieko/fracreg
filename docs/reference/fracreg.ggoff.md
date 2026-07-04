# GGOFF Tests for Fractional Regression Models

`fracreg.ggoff` is used to perform Generalized
Goodness-Of-Functional-Form (GGOFF) tests to check the adequacy of the
functional form and link specification of fractional regression models.

## Usage

``` r
fracreg.ggoff(object, version = "LM", table = TRUE, ...)
```

## Arguments

- object:

  an object containing the results of an `fracreg` command.

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

`fracreg.ggoff` applies the GGOFF, GOFF1 and GOOFF2 test statistics to
fractional regression models estimated via `fracreg`. `fracreg.ggoff`
may be used to test the link specification of: (i) one-part fractional
regression models; (ii) the binary component of two-part fractional
regression models; and (iii) the fractional component of two-part
fractional regression models. When the `Wald` version is implemented, it
is taken into account the option that was chosen for computing standard
errors in the model under evaluation. For the `LM` version, a robust
version is computed in cases (i) and (iii) and a conventional version in
case (ii). See Ramalho, Ramalho and Murteira (2014) for details on the
application of the GGOFF, GOFF1 and GOOFF2 tests in the fractional
regression framework.

## Value

`fracreg.ggoff` returns a named vector with the test results.

## References

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2014), "A generalized
goodness-of-functional form test for binary and fractional regression
models", *Manchester School*, 82(4), 488-507.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md),
for fitting fractional regression models.  
[`fracreg.reset`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.reset.md),
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
#using LM and Wald versions of the GGOFF test, based on 1 or 2 fitted powers of
#the linear predictor
res <- fracreg(y,X,linkfrac="logit",table=FALSE)
fracreg.ggoff(res,c("Wald","LM"))
#> 
#> *** GGOFF test ***
#> 
#> H0:  Fractional logit model
#> 
#>   Test Version Statistic p-value  
#>  GOFF1      LM     0.043   0.835  
#>  GOFF1    Wald     0.043   0.836  
#>  GOFF2      LM     0.100   0.751  
#>  GOFF2    Wald     0.101   0.751  
#>  GGOFF      LM     5.524   0.063 *
#>  GGOFF    Wald     4.973   0.083 *
#> 

#Testing the probit specification of the binary component of a two-part fractional
#regression model using a LR-based GGOFF test
res <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1,table=FALSE)
fracreg.ggoff(res,"LR")
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> 
#> *** GGOFF test ***
#> 
#> H0:  Binary probit component of a two-part model
#> 
#>   Test Version Statistic p-value 
#>  GOFF1      LR     0.163   0.686 
#>  GOFF2      LR     0.034   0.853 
#>  GGOFF      LR     4.009   0.135 
#> 
```
