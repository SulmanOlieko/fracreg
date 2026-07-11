# Fractional Ridge Regression

`fracregridge` implements Fractional Ridge Regression (Rokem & Kay,
2020), which is an approach to regularized linear regression. Unlike
standard ridge regression where the penalty term \\\alpha\\ is chosen
directly, `fracregridge` allows you to specify the desired *fraction* of
the unregularized OLS coefficient vector length. The algorithm then
automatically determines the corresponding \\\alpha\\ penalties.

## Usage

``` r
fracregridge(y, x, fracs = seq(0.1, 1.0, by = 0.1), tol = 1e-10, intercept = TRUE, ...)
```

## Arguments

- y:

  A numeric vector or matrix of the dependent variable(s).

- x:

  A numeric matrix of the explanatory variables.

- fracs:

  A numeric vector indicating the desired fractions of the unregularized
  coefficient vector length. Default is `seq(0.1, 1.0, by=0.1)`. Values
  must be sorted in ascending order.

- tol:

  A numeric tolerance under which singular values of the `x` matrix are
  considered to be zero. Default is `1e-10`.

- intercept:

  logical. If `TRUE`, an intercept is included in the model. Default is
  `TRUE`.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class `fracregridge` containing:

- coef:

  The estimated ridge coefficients for each requested fraction.

- alphas:

  The corresponding \\\alpha\\ penalty values.

- fracs:

  The grid of fractions.

- call:

  The matched call.

## Details

Standard ridge regression minimizes the following objective function:
\$\$L = (y - X\beta)'(y - X\beta) + \alpha \beta'\beta\$\$

The penalty \\\alpha\\ shrinks the coefficients towards zero, reducing
the length of the coefficient vector \\\|\|\beta\|\|\_2\\. However,
choosing \\\alpha\\ can be unintuitive. `fracregridge` re-parameterizes
the problem so the user specifies `fracs`, the fraction of the
unregularized Ordinary Least Squares (OLS) vector length \\\gamma =
\frac{\|\|\beta(\alpha)\|\|\_2}{\|\|\beta(0)\|\|\_2}\\. The function
automatically determines the \\\alpha\\ values corresponding to these
fractions.

## References

Rokem, A., & Kay, K. (2020). Fractional ridge regression: a fast,
interpretable reparameterization of ridge regression. *GigaScience*,
9(12).

Rokem, A., and Kay, K., fracridge: Fractional Ridge Regression. Package
repository. \<https://github.com/nrdg/fracridge\>.

## See also

[`fracreg`](https://sulmanolieko.github.io/fracreg/reference/fracreg.md)

## Author

Sulman Olieko Owili \<oliekosulman@gmail.com\>

## Examples

``` r
# Empirical 401(k) Example
data("fracreg_k401k")
y_401k <- fracreg_k401k$prate
X_401k <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age,
                totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)

# Fit fractional ridge regression for the 401(k) participation rates
mod_401k <- fracregridge(y = y_401k, x = X_401k, fracs = seq(0.2, 1.0, by = 0.2))

# View full detailed summary
summary(mod_401k)
#> 
#> -------------------------------------------------------------------------------- 
#>                           Fractional Ridge Regression 
#> -------------------------------------------------------------------------------- 
#> Data type:                                                       Cross-sectional 
#> Convergence:                                                          Successful 
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.2 
#> -------------------------------------------------------------------------------- 
#> Number of observations:                                                     1534 
#> Pseudo R-squared:                                                        0.04571 
#> Degrees of freedom:                                                      1531.43 
#> -------------------------------------------------------------------------------- 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 1.218e-01  2.623e-03  46.434   <2e-16 ***
#> mrate       6.721e-02  3.418e-03  19.662   <2e-16 ***
#> age         3.471e-02  7.071e-04  49.080   <2e-16 ***
#> totemp      1.033e-06  9.158e-07   1.127     0.26    
#> sole        6.265e-02  2.701e-03  23.196   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.4 
#> -------------------------------------------------------------------------------- 
#> Number of observations:                                                     1534 
#> Pseudo R-squared:                                                        0.06873 
#> Degrees of freedom:                                                      1530.76 
#> -------------------------------------------------------------------------------- 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 2.741e-01  4.766e-03  57.513   <2e-16 ***
#> mrate       1.007e-01  5.152e-03  19.549   <2e-16 ***
#> age         2.461e-02  6.189e-04  39.769   <2e-16 ***
#> totemp      8.942e-07  6.957e-07   1.285    0.199    
#> sole        1.115e-01  4.878e-03  22.855   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.6 
#> -------------------------------------------------------------------------------- 
#> Number of observations:                                                     1534 
#> Pseudo R-squared:                                                        0.08697 
#> Degrees of freedom:                                                      1530.11 
#> -------------------------------------------------------------------------------- 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 4.439e-01  6.218e-03  71.390   <2e-16 ***
#> mrate       9.801e-02  5.457e-03  17.960   <2e-16 ***
#> age         1.586e-02  5.288e-04  29.996   <2e-16 ***
#> totemp      5.080e-07  5.229e-07   0.971    0.331    
#> sole        1.246e-01  6.223e-03  20.027   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.8 
#> -------------------------------------------------------------------------------- 
#> Number of observations:                                                     1534 
#> Pseudo R-squared:                                                        0.10095 
#> Degrees of freedom:                                                      1529.52 
#> -------------------------------------------------------------------------------- 
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  6.162e-01  7.193e-03  85.673   <2e-16 ***
#> mrate        7.665e-02  5.176e-03  14.808   <2e-16 ***
#> age          8.729e-03  4.567e-04  19.114   <2e-16 ***
#> totemp      -1.493e-07  4.081e-07  -0.366    0.714    
#> sole         9.781e-02  7.001e-03  13.970   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                             Target Fraction: frac_1 
#> -------------------------------------------------------------------------------- 
#> Number of observations:                                                     1534 
#> Pseudo R-squared:                                                        0.11402 
#> Degrees of freedom:                                                         1529 
#> -------------------------------------------------------------------------------- 
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  7.827e-01  8.711e-03  89.853  < 2e-16 ***
#> mrate        5.062e-02  5.259e-03   9.625  < 2e-16 ***
#> age          2.822e-03  4.487e-04   6.289 3.20e-10 ***
#> totemp      -9.814e-07  3.689e-07  -2.660  0.00782 ** 
#> sole         4.137e-02  8.275e-03   4.999 5.77e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-11 22:30:15 
#> -------------------------------------------------------------------------------- 
#> 

# Compute Average Partial Effects for Ridge
pe_401k <- fracregridge.pe(mod_401k)
#> 
#> Note: Fractional Ridge Regression is a linear model without a link function.
#> Therefore, the partial effects are mathematically identical to the coefficients themselves.
#> 
summary(pe_401k)
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                             Average partial effects 
#> -------------------------------------------------------------------------------- 
#>                           Fractional Ridge Regression 
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.2 
#> -------------------------------------------------------------------------------- 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 1.218e-01  2.623e-03  46.434   <2e-16 ***
#> mrate       6.721e-02  3.418e-03  19.662   <2e-16 ***
#> age         3.471e-02  7.071e-04  49.080   <2e-16 ***
#> totemp      1.033e-06  9.158e-07   1.127     0.26    
#> sole        6.265e-02  2.701e-03  23.196   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.4 
#> -------------------------------------------------------------------------------- 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 2.741e-01  4.766e-03  57.513   <2e-16 ***
#> mrate       1.007e-01  5.152e-03  19.549   <2e-16 ***
#> age         2.461e-02  6.189e-04  39.769   <2e-16 ***
#> totemp      8.942e-07  6.957e-07   1.285    0.199    
#> sole        1.115e-01  4.878e-03  22.855   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.6 
#> -------------------------------------------------------------------------------- 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 4.439e-01  6.218e-03  71.390   <2e-16 ***
#> mrate       9.801e-02  5.457e-03  17.960   <2e-16 ***
#> age         1.586e-02  5.288e-04  29.996   <2e-16 ***
#> totemp      5.080e-07  5.229e-07   0.971    0.331    
#> sole        1.246e-01  6.223e-03  20.027   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.8 
#> -------------------------------------------------------------------------------- 
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  6.162e-01  7.193e-03  85.673   <2e-16 ***
#> mrate        7.665e-02  5.176e-03  14.808   <2e-16 ***
#> age          8.729e-03  4.567e-04  19.114   <2e-16 ***
#> totemp      -1.493e-07  4.081e-07  -0.366    0.714    
#> sole         9.781e-02  7.001e-03  13.970   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                             Target Fraction: frac_1 
#> -------------------------------------------------------------------------------- 
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  7.827e-01  8.711e-03  89.853  < 2e-16 ***
#> mrate        5.062e-02  5.259e-03   9.625  < 2e-16 ***
#> age          2.822e-03  4.487e-04   6.289 3.20e-10 ***
#> totemp      -9.814e-07  3.689e-07  -2.660  0.00782 ** 
#> sole         4.137e-02  8.275e-03   4.999 5.77e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-11 22:30:15 
#> -------------------------------------------------------------------------------- 
#> 

# Simulated Data Example
set.seed(123)
n <- 100
p <- 10
y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

# Fit Fractional Ridge Regression
# We want the coefficients that correspond to 30\%, 50\%, and 80\% of the OLS length
mod_sim <- fracregridge(y, X, fracs = c(0.3, 0.5, 0.8))

# View brief summary
print(mod_sim)
#> 
#> Fractional Ridge Regression
#> 
#> Call:
#> fracregridge(y = y, x = X, fracs = c(0.3, 0.5, 0.8))
#> 
#> Ridge Coefficients at Target Fractions:
#>                 frac_0.3     frac_0.5    frac_0.8
#> (Intercept)  0.025969041  0.043431002  0.06995586
#> X1          -0.015190037 -0.025277903 -0.04006089
#> X2          -0.029634697 -0.050625055 -0.08391206
#> X3          -0.017972923 -0.035289788 -0.06765584
#> X4          -0.048194568 -0.081684776 -0.13262875
#> X5          -0.013807840 -0.021406027 -0.02901469
#> X6          -0.013152916 -0.022143313 -0.03593916
#> X7           0.048353290  0.078031858  0.11711702
#> X8          -0.006353284 -0.014086532 -0.03160760
#> X9           0.001567879  0.002529806  0.00567404
#> X10          0.020836124  0.031879315  0.04441221
#> 

# Compute Partial Effects
pe_sim <- fracregridge.pe(mod_sim)
#> 
#> Note: Fractional Ridge Regression is a linear model without a link function.
#> Therefore, the partial effects are mathematically identical to the coefficients themselves.
#> 
summary(pe_sim)
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                             Average partial effects 
#> -------------------------------------------------------------------------------- 
#>                           Fractional Ridge Regression 
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.3 
#> -------------------------------------------------------------------------------- 
#>              Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  0.025969   0.026424   0.983   0.3257  
#> X1          -0.015190   0.026179  -0.580   0.5618  
#> X2          -0.029635   0.026080  -1.136   0.2558  
#> X3          -0.017973   0.026498  -0.678   0.4976  
#> X4          -0.048195   0.026318  -1.831   0.0671 .
#> X5          -0.013808   0.025679  -0.538   0.5908  
#> X6          -0.013153   0.026949  -0.488   0.6255  
#> X7           0.048353   0.026532   1.822   0.0684 .
#> X8          -0.006353   0.027076  -0.235   0.8145  
#> X9           0.001568   0.026610   0.059   0.9530  
#> X10          0.020836   0.026763   0.779   0.4362  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.5 
#> -------------------------------------------------------------------------------- 
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  0.04343    0.04478   0.970   0.3322  
#> X1          -0.02528    0.04490  -0.563   0.5735  
#> X2          -0.05063    0.04501  -1.125   0.2607  
#> X3          -0.03529    0.04449  -0.793   0.4277  
#> X4          -0.08168    0.04471  -1.827   0.0677 .
#> X5          -0.02141    0.04462  -0.480   0.6314  
#> X6          -0.02214    0.04502  -0.492   0.6228  
#> X7           0.07803    0.04477   1.743   0.0814 .
#> X8          -0.01409    0.04497  -0.313   0.7541  
#> X9           0.00253    0.04489   0.056   0.9551  
#> X10          0.03188    0.04471   0.713   0.4758  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                            Target Fraction: frac_0.8 
#> -------------------------------------------------------------------------------- 
#>              Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  0.069956   0.074554   0.938   0.3481  
#> X1          -0.040061   0.075645  -0.530   0.5964  
#> X2          -0.083912   0.076139  -1.102   0.2704  
#> X3          -0.067656   0.073864  -0.916   0.3597  
#> X4          -0.132629   0.074935  -1.770   0.0767 .
#> X5          -0.029015   0.077330  -0.375   0.7075  
#> X6          -0.035939   0.072614  -0.495   0.6206  
#> X7           0.117117   0.074087   1.581   0.1139  
#> X8          -0.031608   0.072003  -0.439   0.6607  
#> X9           0.005674   0.073844   0.077   0.9388  
#> X10          0.044412   0.072992   0.608   0.5429  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-11 22:30:15 
#> -------------------------------------------------------------------------------- 
#> 
```
