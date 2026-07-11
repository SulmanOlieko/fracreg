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
#> Fractional Ridge Regression Summary
#> ========================================================================
#> Call:
#> fracregridge(y = y_401k, x = X_401k, fracs = seq(0.2, 1, by = 0.2))
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.2
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 1.2179e-01 2.6228e-03 46.4339   <2e-16 ***
#> mrate       6.7213e-02 3.4183e-03 19.6624   <2e-16 ***
#> age         3.4706e-02 7.0714e-04 49.0797   <2e-16 ***
#> totemp      1.0325e-06 9.1581e-07  1.1275   0.2595    
#> sole        6.2649e-02 2.7009e-03 23.1957   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.4
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 2.7411e-01 4.7660e-03 57.5130   <2e-16 ***
#> mrate       1.0072e-01 5.1522e-03 19.5485   <2e-16 ***
#> age         2.4613e-02 6.1890e-04 39.7691   <2e-16 ***
#> totemp      8.9421e-07 6.9571e-07  1.2853   0.1987    
#> sole        1.1149e-01 4.8779e-03 22.8554   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.6
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 4.4390e-01 6.2180e-03 71.3897   <2e-16 ***
#> mrate       9.8012e-02 5.4572e-03 17.9601   <2e-16 ***
#> age         1.5863e-02 5.2883e-04 29.9964   <2e-16 ***
#> totemp      5.0797e-07 5.2288e-07  0.9715   0.3313    
#> sole        1.2463e-01 6.2228e-03 20.0274   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.8
#>                Estimate  Std. Error z value Pr(>|z|)    
#> (Intercept)  6.1623e-01  7.1928e-03 85.6728   <2e-16 ***
#> mrate        7.6650e-02  5.1761e-03 14.8084   <2e-16 ***
#> age          8.7294e-03  4.5671e-04 19.1138   <2e-16 ***
#> totemp      -1.4935e-07  4.0815e-07 -0.3659   0.7144    
#> sole         9.7811e-02  7.0014e-03 13.9702   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------------
#> Table: frac_1
#>                Estimate  Std. Error z value  Pr(>|z|)    
#> (Intercept)  7.8274e-01  8.7113e-03 89.8533 < 2.2e-16 ***
#> mrate        5.0618e-02  5.2589e-03  9.6252 < 2.2e-16 ***
#> age          2.8218e-03  4.4870e-04  6.2888 3.199e-10 ***
#> totemp      -9.8137e-07  3.6894e-07 -2.6599  0.007815 ** 
#> sole         4.1367e-02  8.2753e-03  4.9989 5.766e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ========================================================================

# Compute Average Partial Effects for Ridge
pe_401k <- fracregridge.pe(mod_401k)
#> 
#> Note: Fractional Ridge Regression is a linear model without a link function.
#> Therefore, the partial effects are mathematically identical to the coefficients themselves.
#> 
summary(pe_401k)
#> 
#> Fractional Ridge Regression Summary
#> ========================================================================
#> Call:
#> fracregridge.pe(object = mod_401k)
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.2
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 1.2179e-01 2.6228e-03 46.4339   <2e-16 ***
#> mrate       6.7213e-02 3.4183e-03 19.6624   <2e-16 ***
#> age         3.4706e-02 7.0714e-04 49.0797   <2e-16 ***
#> totemp      1.0325e-06 9.1581e-07  1.1275   0.2595    
#> sole        6.2649e-02 2.7009e-03 23.1957   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.4
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 2.7411e-01 4.7660e-03 57.5130   <2e-16 ***
#> mrate       1.0072e-01 5.1522e-03 19.5485   <2e-16 ***
#> age         2.4613e-02 6.1890e-04 39.7691   <2e-16 ***
#> totemp      8.9421e-07 6.9571e-07  1.2853   0.1987    
#> sole        1.1149e-01 4.8779e-03 22.8554   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.6
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 4.4390e-01 6.2180e-03 71.3897   <2e-16 ***
#> mrate       9.8012e-02 5.4572e-03 17.9601   <2e-16 ***
#> age         1.5863e-02 5.2883e-04 29.9964   <2e-16 ***
#> totemp      5.0797e-07 5.2288e-07  0.9715   0.3313    
#> sole        1.2463e-01 6.2228e-03 20.0274   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.8
#>                Estimate  Std. Error z value Pr(>|z|)    
#> (Intercept)  6.1623e-01  7.1928e-03 85.6728   <2e-16 ***
#> mrate        7.6650e-02  5.1761e-03 14.8084   <2e-16 ***
#> age          8.7294e-03  4.5671e-04 19.1138   <2e-16 ***
#> totemp      -1.4935e-07  4.0815e-07 -0.3659   0.7144    
#> sole         9.7811e-02  7.0014e-03 13.9702   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------------
#> Table: frac_1
#>                Estimate  Std. Error z value  Pr(>|z|)    
#> (Intercept)  7.8274e-01  8.7113e-03 89.8533 < 2.2e-16 ***
#> mrate        5.0618e-02  5.2589e-03  9.6252 < 2.2e-16 ***
#> age          2.8218e-03  4.4870e-04  6.2888 3.199e-10 ***
#> totemp      -9.8137e-07  3.6894e-07 -2.6599  0.007815 ** 
#> sole         4.1367e-02  8.2753e-03  4.9989 5.766e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ========================================================================

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
#> Fractional Ridge Regression Summary
#> ========================================================================
#> Call:
#> fracregridge.pe(object = mod_sim)
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.3
#>               Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  0.0259690  0.0264237  0.9828  0.32571  
#> X1          -0.0151900  0.0261794 -0.5802  0.56176  
#> X2          -0.0296347  0.0260804 -1.1363  0.25584  
#> X3          -0.0179729  0.0264977 -0.6783  0.49759  
#> X4          -0.0481946  0.0263184 -1.8312  0.06707 .
#> X5          -0.0138078  0.0256786 -0.5377  0.59077  
#> X6          -0.0131529  0.0269490 -0.4881  0.62550  
#> X7           0.0483533  0.0265323  1.8224  0.06839 .
#> X8          -0.0063533  0.0270759 -0.2346  0.81448  
#> X9           0.0015679  0.0266103  0.0589  0.95302  
#> X10          0.0208361  0.0267628  0.7785  0.43625  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.5
#>               Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  0.0434310  0.0447842  0.9698  0.33215  
#> X1          -0.0252779  0.0449011 -0.5630  0.57346  
#> X2          -0.0506251  0.0450115 -1.1247  0.26071  
#> X3          -0.0352898  0.0444894 -0.7932  0.42765  
#> X4          -0.0816848  0.0447092 -1.8270  0.06770 .
#> X5          -0.0214060  0.0446157 -0.4798  0.63138  
#> X6          -0.0221433  0.0450223 -0.4918  0.62284  
#> X7           0.0780319  0.0447740  1.7428  0.08137 .
#> X8          -0.0140865  0.0449742 -0.3132  0.75412  
#> X9           0.0025298  0.0448853  0.0564  0.95505  
#> X10          0.0318793  0.0447049  0.7131  0.47578  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ------------------------------------------------------------------------
#> Table: frac_0.8
#>              Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  0.069956   0.074554  0.9383  0.34808  
#> X1          -0.040061   0.075645 -0.5296  0.59640  
#> X2          -0.083912   0.076139 -1.1021  0.27042  
#> X3          -0.067656   0.073864 -0.9160  0.35969  
#> X4          -0.132629   0.074935 -1.7699  0.07674 .
#> X5          -0.029015   0.077330 -0.3752  0.70751  
#> X6          -0.035939   0.072614 -0.4949  0.62065  
#> X7           0.117117   0.074087  1.5808  0.11392  
#> X8          -0.031608   0.072003 -0.4390  0.66068  
#> X9           0.005674   0.073844  0.0768  0.93875  
#> X10          0.044412   0.072992  0.6085  0.54289  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ========================================================================
```
