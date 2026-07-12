# RESET Test for Fractional Response Regressions under Neglected Heterogeneity

`fracreghet.reset` is used to test the specification of fractional
response models estimated by GMMx or LINx.

## Usage

``` r
fracreghet.reset(
  object,
  lastpower.vec = 3,
  version = "Wald",
  table = FALSE,
  ...
)
```

## Arguments

- object:

  an object containing the results of an `fracreghet` command.

- lastpower.vec:

  a numeric vector containing the maximum powers of the linear
  predictors to be used in RESET tests.

- version:

  a vector containing the test versions to use. Available options:
  `Wald` (the default) and `LM` (only available for `GMMx`).

- table:

  a logical value indicating whether a summary table with the test
  results should be printed.

- ...:

  Arguments to pass to [nlminb](https://rdrr.io/r/stats/nlminb.html),
  which is used to estimate the model under the alternative hypothesis
  when `version` is equal to `"Wald"` and the null model was estimated
  by `GMMx`.

## Value

`fracreghet.reset` returns a named vector with the test results.

## Details

`fracreghet.reset` applies the RESET test statistic to fractional
response models estimated via `fracreghet` using the options `GMMx` or
`LINx`. `fracreghet.reset` may be used to test simultaneously the
validity of the link specification and the transformation applied to the
response variable by each estimator.

**RESET Test under Unobserved Heterogeneity:** The test is based on
augmenting the original model with powers of the linear predictor
\\x\hat{\beta}\\. For GMMx, it tests \\H_0: \gamma = 0\\ in the expanded
moment conditions: \$\$E\left\[Z_i \left(H(y_i) - \exp\left(x_i\beta +
\sum\_{k=2}^P \gamma_k
(x_i\hat{\beta})^k\right)E(e^{c_i})\right)\right\] = 0\$\$ This
simultaneously evaluates whether the mean function and the specific
heterogeneity transformation \\H(\cdot)\\ are correctly specified.

It is taken into account the option that was chosen for computing
standard errors in the model under evaluation. See Ramalho and Ramalho
(2017) for details.

## References

Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of
nonlinear regression models with boundary outcomes and endogeneity, with
applications to nonnegative and fractional responses", *Econometric
Reviews*, 36(4), 397-420.

Ramsey, J.B. (1969), "Tests for Specification Errors in Classical Linear
Least-Squares Regression Analysis", *Journal of the Royal Statistical
Society: Series B (Methodological)*, 31(2), 350-371.

## See also

[`fracreghet`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.md),
for fitting fractional response models under unobserved heterogeneity.  
[`fracreghet.pe`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.pe.md),
for computing partial effects.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

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
res_emp <- fracreghet(y_adj, X_het, type="GMMx", link="logit") 
#> 
#> -------------------------------------------------------------------------------- 
#>         Fractional logit regression with heteroscedasticity/endogeneity 
#> -------------------------------------------------------------------------------- 
#> Data type:                                                       Cross-sectional 
#> Estimator:                                                                  GMMx 
#> Convergence:                                                          Successful 
#> Number of observations:                                                     1534 
#> Wald chi2(2):                                                           153.0331 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                             HC0 
#> -------------------------------------------------------------------------------- 
#>                               Final GMMx estimates 
#> -------------------------------------------------------------------------------- 
#>             Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
#> (Intercept)     6.86160         0.15854 43.28078    6.55087     7.172  < 2e-16
#> mrate           0.39342         0.03450 11.40432    0.32581     0.461  < 2e-16
#> ltotemp        -0.16558         0.02516 -6.58147   -0.21489    -0.116 4.66e-11
#>                
#> (Intercept) ***
#> mrate       ***
#> ltotemp     ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-13 01:19:01 
#> -------------------------------------------------------------------------------- 
#> 
reset_res <- fracreghet.reset(res_emp)
summary(reset_res)
#> 
#> -------------------------------------------------------------------------------- 
#>                                    RESET test 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#> H0: Estimator: GMMx 
#> -------------------------------------------------------------------------------- 
#>         Statistic p-value    
#> Wald(3)     47.56 4.7e-11 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-13 01:19:01 
#> -------------------------------------------------------------------------------- 
 
### Simulated Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

Z <- cbind(rnorm(N),rnorm(N),rnorm(N))
dimnames(Z)[[2]] <- c("Z1","Z2","Z3")

y <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))

mod <- fracreghet(y,X,type="GMMx")
#> 
#> -------------------------------------------------------------------------------- 
#>         Fractional logit regression with heteroscedasticity/endogeneity 
#> -------------------------------------------------------------------------------- 
#> Data type:                                                       Cross-sectional 
#> Estimator:                                                                  GMMx 
#> Convergence:                                                          Successful 
#> Number of observations:                                                      250 
#> Wald chi2(2):                                                           345.2713 
#> Prob > chi2:                                                              0.0000 
#> Standard errors:                                                             HC0 
#> -------------------------------------------------------------------------------- 
#>                               Final GMMx estimates 
#> -------------------------------------------------------------------------------- 
#>             Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
#> (Intercept)     0.54440         0.08976  6.06483    0.36847     0.720 1.32e-09
#> X1              0.99639         0.08221 12.12082    0.83527     1.158  < 2e-16
#> X2              0.91896         0.09997  9.19210    0.72302     1.115  < 2e-16
#>                
#> (Intercept) ***
#> X1          ***
#> X2          ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-13 01:19:01 
#> -------------------------------------------------------------------------------- 
#> 

#LM and Wald versions of the RESET test, based on 1 or 2 fitted powers of xb
reset_res <- fracreghet.reset(mod,2:3,c("Wald","LM"))
summary(reset_res)
#> 
#> -------------------------------------------------------------------------------- 
#>                                    RESET test 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit regression 
#> -------------------------------------------------------------------------------- 
#> H0: Estimator: GMMx 
#> -------------------------------------------------------------------------------- 
#>         Statistic p-value  
#> LM(2)       0.236  0.6272  
#> Wald(2)     0.196  0.6577  
#> LM(3)       3.411  0.1817  
#> Wald(3)     6.050  0.0486 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-13 01:19:01 
#> -------------------------------------------------------------------------------- 
```
