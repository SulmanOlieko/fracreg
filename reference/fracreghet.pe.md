# Fractional Response Regressions under Unobserved Heterogeneity - Partial Effects

`fracreghet.pe` is used to compute average and/or conditional partial
effects in fractional response models under unobserved heterogeneity.

## Usage

``` r
fracreghet.pe(object, smearing = T, APE = T, CPE = F, at = NULL, 
              which.x = NULL, table = FALSE, variance = T)
```

## Arguments

- object:

  an object containing the results of an `fracreghet` command.

- smearing:

  a logical value indicating whether the smearing correction is to be
  applied

- APE:

  a logical value indicating whether average partial effects are to be
  computed.

- CPE:

  a logical value indicating whether conditional partial effects are to
  be computed.

- at:

  a numeric vector containing the covariates' values at which the
  conditional partial effects are to be computed or the strings `"mean"`
  (the default) or `"median"`, in which cases the covariates are
  evaluated at their mean or median values (or mode, in case of dummy
  variables), respectively.

- which.x:

  a vector containing the names of the covariates to which the partial
  effects are to be computed.

- table:

  a logical value indicating whether a summary table with the results
  should be printed.

- variance:

  a logical value indicating whether the variance of the estimated
  partial effects should be calculated. Defaults to `TRUE` whenever
  `table = TRUE`.

## Details

`fracreghet.pe` calculates partial effects for fractional response
models estimated via `fracreghet`. `fracreghet.pe` may be used to
compute average or conditional partial effects. These partial effects
may be conditional only on observables, using the smearing estimator, or
also on unobservables, setting the error term to zero.

**Partial Effects under Unobserved Heterogeneity:** When unobserved
heterogeneity or endogeneity is present, calculating partial effects
requires dealing with the unobserved error \\v_i\\. Let the conditional
mean be \\E(y\|x, v) = G(x\beta + \gamma v)\\. - **Conditional on
Observables (Smearing):** The unobserved heterogeneity is integrated out
over its empirical distribution. The average partial effect for a
continuous variable \\x_k\\ is computed as: \$\$PE_k(x) = \frac{1}{N}
\sum\_{i=1}^N g(x\beta + \gamma \hat{v}\_i) \beta_k\$\$ - **Conditional
on Unobservables (Error = 0):** The partial effect is evaluated for an
individual with the mean level of unobserved heterogeneity (\\v = 0\\):
\$\$PE_k(x) = g(x\beta) \beta_k\$\$

For discrete variables, the partial effects are calculated as the
discrete differences evaluated using either the smearing approach or
setting the error term to zero.

For calculating standard errors, it is taken into account the option
that was previously chosen for estimating the model. See Ramalho and
Ramalho (2017) for details on the computation of partial effects for
fractional response models under unobserved heterogeneity.

## Value

`fracreghet.pe` returns a list with the following element:

- PE.p:

  a named vector of partial effects.

If `variance = TRUE` or `table = TRUE`, the previous list also contains
the following element:

- PE.sd:

  a named vector of standard errors of the estimated partial effects.

When both average and conditional partial effects are requested, two
lists containing the previous elements are returned, indexed by the
prefixes `ape` and `cpe`.

## References

Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of
nonlinear regression models with boundary outcomes and endogeneity, with
applications to nonnegative and fractional responses", *Econometric
Reviews*, 36(4), 397-420.

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## See also

[`fracreghet`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.md),
for fitting fractional response models under unobserved heterogeneity.  
[`fracreghet.reset`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.reset.md),
for the RESET test.  

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
res_emp <- fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], 
                      type="QMLxv", link="logit") 
#> 
#> -------------------------------------------------------------------------------- 
#>         Fractional logit regression with heteroscedasticity/endogeneity 
#> -------------------------------------------------------------------------------- 
#> Estimator:                                                                 QMLxv 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                     1534 
#> Standard errors:                                                          robust 
#> Wald chi2(6):                                                          1991.8748 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                   Final Quasi-Maximum Likelihood xv estimates 
#> -------------------------------------------------------------------------------- 
#>             Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
#> (Intercept)    -0.10561         0.75224 -0.14040   -1.57997     1.369    0.888
#> mrate           3.72138         0.68952  5.39703    2.36994     5.073 6.78e-08
#> ltotemp        -0.07009         0.05348 -1.31060   -0.17491     0.035    0.190
#> vhat           -2.79515         0.70026 -3.99157   -4.16764    -1.423 6.56e-05
#>                
#> (Intercept)    
#> mrate       ***
#> ltotemp        
#> vhat        ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                                  Reduced form: 
#> -------------------------------------------------------------------------------- 
#>               Coefficient Robust Std.Err.  z value [95% Conf. Interval]
#> Z_(Intercept)     0.95046         0.09550  9.95211    0.76327     1.138
#> Z_age             0.01146         0.00231  4.95997    0.00693     0.016
#> Z_ltotemp        -0.05534         0.01421 -3.89528   -0.08318    -0.027
#>               Pr(>|z|)    
#> Z_(Intercept)  < 2e-16 ***
#> Z_age         7.05e-07 ***
#> Z_ltotemp     9.81e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-11 22:01:27 
#> -------------------------------------------------------------------------------- 
#> 
pe_res <- fracreghet.pe(res_emp, which.x="mrate")
summary(pe_res)
#> 
#> 
#> -------------------------------------------------------------------------------- 
#> Average partial effects (conditional only on observables, based on the smearing estimator) 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit regression 
#>                                 Estimator: QMLxv 
#> -------------------------------------------------------------------------------- 
#>       Estimate Std. Error z value Pr(>|z|)    
#> mrate  0.39498    0.09827   4.019 5.84e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-11 22:01:27 
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
#> Estimator:                                                                  GMMx 
#> Data type:                                                       Cross-sectional 
#> Number of observations:                                                      250 
#> Standard errors:                                                          robust 
#> Wald chi2(2):                                                           368.8403 
#> Prob > chi2:                                                              0.0000 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                               Final GMMx estimates 
#> -------------------------------------------------------------------------------- 
#>             Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
#> (Intercept)     0.43394         0.07745  5.60291    0.28214     0.586 2.11e-08
#> X1              0.98087         0.06730 14.57560    0.84897     1.113  < 2e-16
#> X2              0.88338         0.07066 12.50206    0.74489     1.022  < 2e-16
#>                
#> (Intercept) ***
#> X1          ***
#> X2          ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-11 22:01:27 
#> -------------------------------------------------------------------------------- 
#> 

#Smearing estimator of average partial effects for variable X1
pe_res <- fracreghet.pe(mod,which.x="X1")
summary(pe_res)
#> 
#> 
#> -------------------------------------------------------------------------------- 
#> Average partial effects (conditional only on observables, based on the smearing estimator) 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit regression 
#>                                 Estimator: GMMx 
#> -------------------------------------------------------------------------------- 
#>    Estimate Std. Error z value Pr(>|z|)    
#> X1 0.166732   0.009667   17.25   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-11 22:01:27 
#> -------------------------------------------------------------------------------- 

#Naive estimator of conditional partial effects for all covariates,
#which are evaluated at X1=1 and X2=-1
pe_res <- fracreghet.pe(mod,smearing=FALSE,APE=FALSE,CPE=TRUE,at=c(1,-1))
summary(pe_res)
#> 
#> 
#> -------------------------------------------------------------------------------- 
#> Conditional partial effects (conditional on both observables and unobservables, with error term = 0) 
#> -------------------------------------------------------------------------------- 
#>                           Fractional logit regression 
#>                                 Estimator: GMMx 
#> -------------------------------------------------------------------------------- 
#>    Estimate Std. Error z value Pr(>|z|)    
#> X1  0.22869    0.01337   17.11   <2e-16 ***
#> X2  0.20596    0.01454   14.17   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-11 22:01:27 
#> -------------------------------------------------------------------------------- 
#> 
#> Note: covariates evaluated at the following values:
#> 
#> X1 X2 
#>  1 -1 
```
