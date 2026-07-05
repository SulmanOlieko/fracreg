# GGOFF Tests for Fractional Response Models

`fracreg.ggoff` is used to perform Generalized
Goodness-Of-Functional-Form (GGOFF) tests to check the adequacy of the
functional form and link specification of fractional response models.

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
fractional response models estimated via `fracreg`. `fracreg.ggoff` may
be used to test the link specification of: (i) one-part fractional
response models; (ii) the binary component of two-part fractional
response models; and (iii) the fractional component of two-part
fractional response models.

**GGOFF Test Framework:** The Generalized Goodness-of-Functional Form
(GGOFF) test evaluates the adequacy of the link function \\G(\cdot)\\.
It is based on augmenting the baseline model with specific directions of
departure. The auxiliary testing equation takes the form: \$\$E(y\|x) =
G\left(x\beta + \gamma_1 \frac{g'(x\hat{\beta})}{g(x\hat{\beta})} +
\gamma_2 x\hat{\beta} \right)\$\$ where \\g(\cdot)\\ and \\g'(\cdot)\\
are the first and second derivatives of \\G(\cdot)\\ evaluated at the
linear predictor \\x\hat{\beta}\\. The test checks \\H_0: \gamma_1 = 0,
\gamma_2 = 0\\. GOFF1 and GOFF2 are variants testing individual
components.

When the `Wald` version is implemented, it is taken into account the
option that was chosen for computing standard errors in the model under
evaluation. For the `LM` version, a robust version is computed in cases
(i) and (iii) and a conventional version in case (ii). See Ramalho,
Ramalho and Murteira (2014) for details on the application of the GGOFF,
GOFF1 and GOOFF2 tests in the fractional response framework.

## Value

`fracreg.ggoff` returns a named vector with the test results.

## References

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2014), "A generalized
goodness-of-functional form test for binary and fractional response
models", *Manchester School*, 82(4), 488-507.

Pregibon, D. (1980), "Goodness of Link Tests for Generalized Linear
Models", *Journal of the Royal Statistical Society: Series C (Applied
Statistics)*, 29(1), 15-24.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md),
for fitting fractional response models.  
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

#Testing the logit specification of a standard fractional response model
#using LM and Wald versions of the GGOFF test, based on 1 or 2 fitted powers of
#the linear predictor
res <- fracreg(y,X,linkfrac="logit",table=FALSE)
fracreg.ggoff(res,c("Wald","LM"))
#> 
#> -------------------------------------------------------------------------------- 
#>                                    GGOFF test 
#> -------------------------------------------------------------------------------- 
#> H0: Fractional logit model 
#> -------------------------------------------------------------------------------- 
#>              Statistic p-value  
#> GOFF1 - LM       0.060  0.8062  
#> GOFF1 - Wald     0.060  0.8064  
#> GOFF2 - LM       0.173  0.6777  
#> GOFF2 - Wald     0.171  0.6791  
#> GGOFF - LM       6.895  0.0318 *
#> GGOFF - Wald     5.989  0.0501 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 22:24:41 
#> -------------------------------------------------------------------------------- 

#Testing the probit specification of the binary component of a two-part fractional
#regression model using a LR-based GGOFF test
res <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1,table=FALSE)
fracreg.ggoff(res,"LR")
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> 
#> -------------------------------------------------------------------------------- 
#>                                    GGOFF test 
#> -------------------------------------------------------------------------------- 
#> H0: Binary probit component of a two-part model 
#> -------------------------------------------------------------------------------- 
#>            Statistic p-value
#> GOFF1 - LR     0.002   0.964
#> GOFF2 - LR     0.039   0.844
#> GGOFF - LR     0.717   0.699
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-05 22:24:41 
#> -------------------------------------------------------------------------------- 
```
