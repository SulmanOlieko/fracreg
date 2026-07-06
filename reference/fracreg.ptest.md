# P Test for Fractional Response Models

`fracreg.ptest` is used to perform the P test to evaluate the
specification of alternative, non-nested fractional response models by
testing against each other.

## Usage

``` r
fracreg.ptest(object1, object2, version = "Wald", table = TRUE)
```

## Arguments

- object1:

  an object containing the results of an `fracreg` command.

- object2:

  an object containing the results of another `fracreg` command.

- version:

  a vector containing the test versions to use. Available options:
  `Wald` (the default) and `LM`. Both options may be chosen at the same
  time and are computed in a robust way.

- table:

  a logical value indicating whether a summary table with the test
  results should be printed.

## Details

`fracreg.ptest` applies the P test statistic proposed by Davidson and
MacKinnon (1981) to fractional response models estimated via`fracreg`.
`fracreg.ptest` may be used to test against each other two alternative
specifications for the link function in: (i) one-part fractional
response models; (ii) the binary components of two-part and three-part
fractional response models; (iii) the fractional components of two-part
and three-part fractional response models; and (iv) two-part and
three-part fractional response models.

**P Test Framework:** The P test allows the comparison of non-nested
models (e.g., alternative link functions or non-nested regressors). Let
model 1 specify \\E_1(y\|x) = G(x\beta)\\ and model 2 specify
\\E_2(y\|x) = H(z\theta)\\. To test model 1 against model 2, the
baseline model is augmented with the difference between the fitted
values: \$\$E(y\|x) = G\left(x\beta + \gamma \left( \hat{y}\_{M2} -
\hat{y}\_{M1} \right)\right)\$\$ where \\\hat{y}\_{M1} =
G(x\hat{\beta})\\ and \\\hat{y}\_{M2} = H(z\hat{\theta})\\. The null
hypothesis that model 1 is correct is tested via \\H_0: \gamma = 0\\.

In addition, `fracreg.ptest` may be used to test one-part models against
two-part or three-part models and in cases where the link functions are
the same but the regressors are non-nested. See Ramalho, Ramalho and
Murteira (2011) for details on the application of the P test in the
fractional response framework.

## Value

`fracreg.reset` returns a named vector with the test results.

## References

Davidson, R. and J.G. MacKinnon (1981), "Several tests for model
specification on the presence of alternative hypotheses",
*Econometrica*, 49(3), 781-793.

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional response
models", *Journal of Economic Surveys*, 25(1), 19-68.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreg`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md),
for fitting fractional response models.  
[`fracreg.reset`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.reset.md)
and
[`fracreg.ggoff`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ggoff.md),
for specification tests.  
[`fracreg.pe`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md),
for computing partial effects.

## Examples

``` r
### Empirical 401(k) Examples
data("fracreg_k401k")
y <- fracreg_k401k$prate
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)

m1 <- fracreg(y, X, type="1P", linkfrac="logit")
#> 
#> -------------------------------------------------------------------------------- 
#>                              Fractional logit model 
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
#>                          Run Date: 2026-07-06 01:43:40 
#> -------------------------------------------------------------------------------- 
#> 
m2 <- fracreg(y, X, type="1P", linkfrac="probit")
#> 
#> -------------------------------------------------------------------------------- 
#>                             Fractional probit model 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                   QML 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1534 
#> Pseudo R-squared:                                                        0.14069 
#> Log pseudolikelihood:                                                  -554.2692 
#> Wald chi2(4):                                                           148.7649 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                          robust 
#> Small sample correction:                                                   FALSE 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                     Final Quasi-Maximum Likelihood estimates 
#> -------------------------------------------------------------------------------- 
#>          Coefficient Robust Std.Err.    z value [95% Conf. Interval] Pr(>|z|)
#> Constant   6.374e-01       4.517e-02  1.411e+01  5.489e-01     0.726  < 2e-16
#> mrate      4.190e-01       6.820e-02  6.143e+00  2.853e-01     0.553 8.07e-10
#> age        1.483e-02       2.536e-03  5.846e+00  9.855e-03     0.020 5.04e-09
#> totemp    -4.546e-06       1.732e-06 -2.625e+00 -7.940e-06     0.000  0.00866
#> sole       2.000e-01       4.311e-02  4.639e+00  1.155e-01     0.284 3.50e-06
#>             
#> Constant ***
#> mrate    ***
#> age      ***
#> totemp   ** 
#> sole     ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-06 01:43:40 
#> -------------------------------------------------------------------------------- 
#> 
fracreg.ptest(m1, m2)
#> 
#> -------------------------------------------------------------------------------- 
#>                                      P test 
#> -------------------------------------------------------------------------------- 
#> H0: Fractional logit model 
#> H1: Fractional probit model 
#> -------------------------------------------------------------------------------- 
#>      Statistic p-value
#> Wald    -1.483   0.138
#> -------------------------------------------------------------------------------- 
#> H0: Fractional probit model 
#> H1: Fractional logit model 
#> -------------------------------------------------------------------------------- 
#>      Statistic p-value   
#> Wald     2.754 0.00595 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-06 01:43:40 
#> -------------------------------------------------------------------------------- 

### Simulated Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Testing logit versus loglog specifications for standard fractional
#regression models using a LM version of the P test
res1 <- fracreg(y,X,linkfrac="logit",table=FALSE)
res2 <- fracreg(y,X,linkfrac="loglog",table=FALSE)
fracreg.ptest(res1,res2,"LM")
#> 
#> -------------------------------------------------------------------------------- 
#>                                      P test 
#> -------------------------------------------------------------------------------- 
#> H0: Fractional logit model 
#> H1: Fractional loglog model 
#> -------------------------------------------------------------------------------- 
#>    Statistic p-value  
#> LM     4.029  0.0447 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> H0: Fractional loglog model 
#> H1: Fractional logit model 
#> -------------------------------------------------------------------------------- 
#>    Statistic p-value   
#> LM     9.416 0.00215 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-06 01:43:40 
#> -------------------------------------------------------------------------------- 

#Testing a logit one-part fractional response model versus a binary logit +
#fractional probit two-part model using a Wald version of the P test
res1 <- fracreg(y,X,linkfrac="logit",table=FALSE)
res2 <- fracreg(y,X,linkbin="logit",linkfrac="probit",type="2P",inf=1,table=FALSE)
fracreg.ptest(res1,res2,"Wald")
#> 
#> -------------------------------------------------------------------------------- 
#>                                      P test 
#> -------------------------------------------------------------------------------- 
#> H0: Fractional logit model 
#> H1: Binary logit + Fractional probit two-part model 
#> -------------------------------------------------------------------------------- 
#>      Statistic p-value   
#> Wald    -2.659 0.00835 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#> H0: Binary logit + Fractional probit two-part model 
#> H1: Fractional logit model 
#> -------------------------------------------------------------------------------- 
#>      Statistic p-value    
#> Wald      12.6  <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-06 01:43:40 
#> -------------------------------------------------------------------------------- 
```
