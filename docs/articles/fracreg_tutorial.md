# A Comprehensive Guide to Fractional Response Regressions with fracreg

## Introduction

The **`fracreg`** package is the most comprehensive tool available in R
for estimating, diagnosing, and interpreting fractional response models.
It supports:

- **One-part models (1P)**: A single process governs the entire
  distribution of $`y \in [0, 1]`$.
- **Two-part hurdle models (2P)**: Separately models the binary hurdle
  ($`y > 0`$) and the continuous fractional part ($`0 < y \le 1`$).
- **Three-part double-inflated models (3P)**: Handles inflation at
  **both** 0 and 1.
- **Panel data models** (`fracregpd`): CRE, QML, and GMM approaches for
  longitudinal data.
- **Endogeneity correction** (`fracreghet`): Control function and IV-GMM
  approaches.

This vignette walks you through comprehensive empirical examples (using
the 401(k) dataset) and simulated examples for each estimator.

``` r
library(fracreg)
```

## 1. Data Description and Preparation

We use the built-in `fracreg_k401k` dataset, which is the canonical
firm-level 401(k) plan participation data used in **Papke and Wooldridge
(1996)** (*“Econometric methods for fractional response variables with
an application to 401(k) plan participation rates”*, Journal of Applied
Econometrics). The dataset contains 1,534 observations of 401(k) plans.

The primary variables we use are: - `prate`: The plan participation rate
(the fraction of eligible employees who are active participants). It
strictly falls in the $`[0, 1]`$ interval and is our dependent variable
($`y`$). - `mrate`: The firm’s matching rate (the firm’s contribution
per \$1 of employee contribution). - `age`: The age of the 401(k)
plan. - `totemp`: Total number of employees at the firm. - `sole`: A
binary indicator equal to 1 if the 401(k) plan is the sole retirement
plan offered by the firm.

------------------------------------------------------------------------

## 2. Cross-Sectional Fractional Models (`fracreg`)

The core
[`fracreg()`](https://sulmanolieko.github.io/fracreg/reference/fracreg.md)
function is designed for univariate models where the dependent variable
is bounded between 0 and 1 inclusive ($`0 \le y \le 1`$). The examples
below demonstrate both empirical 401(k) plan participation data and
general simulated boundaries.

``` r
### Empirical 401(k) Examples 
data("fracreg_k401k") 
y <- fracreg_k401k$prate 
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age,  
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole) 
 
# 1P Model 
fracreg(y, X, type="1P", linkfrac="logit") 
```

Toggle to see the output

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
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:57 
    #> --------------------------------------------------------------------------------

``` r

# 1P Model reporting odds ratios and 99% confidence intervals
fracreg(y, X, type="1P", linkfrac="logit", or=TRUE, level=0.99)
```

Toggle to see the output

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
    #>          Odds Ratio Robust Std.Err.    z value [99% Conf. Interval] Pr(>|z|)
    #> Constant  2.539e+00       2.134e-01  1.108e+01  2.044e+00     3.153  < 2e-16
    #> mrate     2.594e+00       3.556e-01  6.951e+00  1.822e+00     3.692 3.62e-12
    #> age       1.028e+00       5.015e-03  5.723e+00  1.015e+00     1.041 1.05e-08
    #> totemp    1.000e+00       3.061e-06 -2.673e+00  1.000e+00     1.000  0.00751
    #> sole      1.406e+00       1.134e-01  4.222e+00  1.142e+00     1.730 2.43e-05
    #>             
    #> Constant ***
    #> mrate    ***
    #> age      ***
    #> totemp   ** 
    #> sole     ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:57 
    #> --------------------------------------------------------------------------------

``` r
 
# 2P Model (modelling mass at 1) 
fracreg(y, X, type="2P", inflation=1, linkbin="logit", linkfrac="logit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                         Part 1: Binary logit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1534 
    #> Pseudo R-squared:                                                         0.1485 
    #> Log-likelihood:                                                        -938.1759 
    #> Wald chi2(4):                                                           173.5169 
    #> Prob > chi2:                                                              0.0000 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                        Final Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient EIM Std.Err.    z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant  -1.396e+00    1.270e-01 -1.099e+01 -1.645e+00    -1.147   <2e-16 ***
    #> mrate      9.053e-01    9.699e-02  9.334e+00  7.152e-01     1.095   <2e-16 ***
    #> age        1.156e-02    6.218e-03  1.858e+00 -6.312e-04     0.024   0.0631 .  
    #> totemp    -1.418e-05    6.324e-06 -2.242e+00 -2.657e-05     0.000   0.0249 *  
    #> sole       8.651e-01    1.131e-01  7.651e+00  6.435e-01     1.087    2e-14 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                       Part 2: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                      852 
    #> Pseudo R-squared:                                                        0.10004 
    #> Log pseudolikelihood:                                                  -450.8391 
    #> Wald chi2(4):                                                            65.4063 
    #> Prob > chi2:                                                              0.0000 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                     Final Quasi-Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient Robust Std.Err.    z value [95% Conf. Interval] Pr(>|z|)
    #> Constant   7.460e-01       6.850e-02  1.089e+01  6.118e-01     0.880  < 2e-16
    #> mrate      3.877e-01       9.725e-02  3.987e+00  1.971e-01     0.578 6.69e-05
    #> age        2.562e-02       4.010e-03  6.390e+00  1.777e-02     0.033 1.66e-10
    #> totemp    -4.061e-06       3.073e-06 -1.322e+00 -1.008e-05     0.000    0.186
    #> sole      -1.510e-02       6.556e-02 -2.303e-01 -1.436e-01     0.113    0.818
    #>             
    #> Constant ***
    #> mrate    ***
    #> age      ***
    #> totemp      
    #> sole        
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>              Two-part regression - binary logit + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #> Data type:                                                       Cross-sectional 
    #> Pseudo R-squared:                                                        0.11243 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> --------------------------------------------------------------------------------

``` r
 
# 3P Model (inject artificial 0s for demonstration) 
y_3p <- y; y_3p[1:50] <- 0 
fracreg(y_3p, X, type="3P", linkbin=c("logit","logit"), linkfrac="logit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                         Part 1: Binary logit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1534 
    #> Pseudo R-squared:                                                        0.00324 
    #> Log-likelihood:                                                        -216.6222 
    #> Wald chi2(4):                                                             3.7679 
    #> Prob > chi2:                                                              0.4383 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                        Final Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient EIM Std.Err.    z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant   2.867e+00    3.157e-01  9.080e+00  2.248e+00     3.486   <2e-16 ***
    #> mrate      1.147e-01    2.131e-01  5.381e-01 -3.030e-01     0.532    0.590    
    #> age        8.375e-03    1.765e-02  4.745e-01 -2.622e-02     0.043    0.635    
    #> totemp     1.036e-04    6.959e-05  1.489e+00 -3.275e-05     0.000    0.136    
    #> sole       3.371e-01    2.983e-01  1.130e+00 -2.476e-01     0.922    0.259    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                         Part 2: Binary logit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1484 
    #> Pseudo R-squared:                                                        0.15199 
    #> Log-likelihood:                                                        -903.5457 
    #> Wald chi2(4):                                                            172.201 
    #> Prob > chi2:                                                              0.0000 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                        Final Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient EIM Std.Err.    z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant  -1.435e+00    1.306e-01 -1.099e+01 -1.691e+00    -1.179  < 2e-16 ***
    #> mrate      9.244e-01    9.877e-02  9.360e+00  7.309e-01     1.118  < 2e-16 ***
    #> age        1.198e-02    6.339e-03  1.890e+00 -4.453e-04     0.024   0.0588 .  
    #> totemp    -1.371e-05    6.317e-06 -2.170e+00 -2.609e-05     0.000   0.0300 *  
    #> sole       8.852e-01    1.154e-01  7.674e+00  6.591e-01     1.111 1.67e-14 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                       Part 3: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                      826 
    #> Pseudo R-squared:                                                        0.09937 
    #> Log pseudolikelihood:                                                  -437.3715 
    #> Wald chi2(4):                                                            62.3685 
    #> Prob > chi2:                                                              0.0000 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                     Final Quasi-Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient Robust Std.Err.    z value [95% Conf. Interval] Pr(>|z|)
    #> Constant   7.388e-01       7.039e-02  1.050e+01  6.008e-01     0.877  < 2e-16
    #> mrate      3.960e-01       1.019e-01  3.885e+00  1.962e-01     0.596 0.000102
    #> age        2.531e-02       4.068e-03  6.223e+00  1.734e-02     0.033 4.88e-10
    #> totemp    -3.817e-06       3.088e-06 -1.236e+00 -9.869e-06     0.000 0.216472
    #> sole      -4.483e-03       6.672e-02 -6.719e-02 -1.353e-01     0.126 0.946434
    #>             
    #> Constant ***
    #> mrate    ***
    #> age      ***
    #> totemp      
    #> sole        
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>      Three-part regression - binary logit , binary logit + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #> Data type:                                                       Cross-sectional 
    #> Pseudo R-squared:                                                        0.07934 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> --------------------------------------------------------------------------------

``` r
 
### Simulated Examples 
 
set.seed(123) 
N <- 1000 
x1 <- rnorm(N) 
x2 <- runif(N) 
 
# Generating a fractional dependent variable with inflation at 0 and 1 
XB <- -0.5 + 0.8 * x1 + 1.2 * x2 + rnorm(N) 
y_latent <- exp(XB) / (1 + exp(XB)) 
 
y <- y_latent 
# Inflate at boundaries 
y[y_latent < 0.2] <- 0 
y[y_latent > 0.8] <- 1 
 
X <- cbind(x1 = x1, x2 = x2) 
 
# fracreg estimation of a logit fractional response model 
fracreg(y, X, type="1P", linkfrac="logit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                           Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1000 
    #> Pseudo R-squared:                                                         0.3903 
    #> Log pseudolikelihood:                                                   -614.973 
    #> Wald chi2(2):                                                           472.2453 
    #> Prob > chi2:                                                              0.0000 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                     Final Quasi-Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant    -0.56969         0.06750 -8.43960   -0.70199    -0.437   <2e-16 ***
    #> x1           0.78822         0.04015 19.63424    0.70954     0.867   <2e-16 ***
    #> x2           1.35611         0.11899 11.39668    1.12289     1.589   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> --------------------------------------------------------------------------------

``` r
 
# fracreg estimation of the binary logit component of the two-part fractional 
# regression model with y=0 as the relevant boundary value 
fracreg(y, X, type="2Pbin", inflation=0, linkbin="logit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                         Part 1: Binary logit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1000 
    #> Pseudo R-squared:                                                        0.14527 
    #> Log-likelihood:                                                        -295.2068 
    #> Wald chi2(2):                                                           102.5171 
    #> Prob > chi2:                                                              0.0000 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                        Final Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient EIM Std.Err. z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant      1.4638       0.1933  7.5732     1.0849     1.843 3.64e-14 ***
    #> x1            1.1857       0.1287  9.2126     0.9335     1.438  < 2e-16 ***
    #> x2            2.2279       0.3932  5.6662     1.4572     2.999 1.46e-08 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> --------------------------------------------------------------------------------

``` r
 
# fracreg estimation of the fractional component of the two-part fractional 
# regression model with y=0 as the relevant boundary value and using a 
# probit link function 
fracreg(y, X, type="2Pfrac", inflation=0, linkfrac="probit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                       Part 2: Fractional probit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                      881 
    #> Pseudo R-squared:                                                        0.32474 
    #> Log pseudolikelihood:                                                  -555.3304 
    #> Wald chi2(2):                                                           391.6944 
    #> Prob > chi2:                                                              0.0000 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                     Final Quasi-Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant    -0.10334         0.03566 -2.89793   -0.17323    -0.033  0.00376 ** 
    #> x1           0.37512         0.02129 17.62293    0.33340     0.417  < 2e-16 ***
    #> x2           0.61326         0.06529  9.39237    0.48529     0.741  < 2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> --------------------------------------------------------------------------------

``` r
 
# fracreg estimation of both components of a two-part fractional response model 
# with y=0 as the relevant boundary value and using a cloglog binary link 
# function and a logit fractional link function 
fracreg(y, X, type="2P", inflation=0, linkbin="cloglog", linkfrac="logit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                        Part 1: Binary cloglog regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1000 
    #> Pseudo R-squared:                                                        0.14837 
    #> Log-likelihood:                                                        -291.7986 
    #> Wald chi2(2):                                                           100.5683 
    #> Prob > chi2:                                                              0.0000 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                        Final Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient EIM Std.Err. z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant     0.43428      0.08825 4.92113    0.26132     0.607 8.60e-07 ***
    #> x1           0.55192      0.06009 9.18429    0.43414     0.670  < 2e-16 ***
    #> x2           1.05359      0.17403 6.05409    0.71250     1.395 1.41e-09 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                       Part 2: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                      881 
    #> Pseudo R-squared:                                                        0.32412 
    #> Log pseudolikelihood:                                                  -555.4465 
    #> Wald chi2(2):                                                           368.0384 
    #> Prob > chi2:                                                              0.0000 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                     Final Quasi-Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant    -0.17323         0.05786 -2.99384   -0.28664    -0.060  0.00275 ** 
    #> x1           0.61205         0.03554 17.22078    0.54239     0.682  < 2e-16 ***
    #> x2           1.00509         0.10703  9.39074    0.79531     1.215  < 2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>             Two-part regression - binary cloglog + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #> Data type:                                                       Cross-sectional 
    #> Pseudo R-squared:                                                        0.38829 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> --------------------------------------------------------------------------------

``` r
 
# Three-part double-inflated model (y has both 0s and 1s) 
fracreg(y, X, type="3P", linkbin=c("logit","probit"), linkfrac="logit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                         Part 1: Binary logit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1000 
    #> Pseudo R-squared:                                                        0.14527 
    #> Log-likelihood:                                                        -295.2068 
    #> Wald chi2(2):                                                           102.5171 
    #> Prob > chi2:                                                              0.0000 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                        Final Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient EIM Std.Err. z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant      1.4638       0.1933  7.5732     1.0849     1.843 3.64e-14 ***
    #> x1            1.1857       0.1287  9.2126     0.9335     1.438  < 2e-16 ***
    #> x2            2.2279       0.3932  5.6662     1.4572     2.999 1.46e-08 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                         Part 2: Binary probit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                      881 
    #> Pseudo R-squared:                                                        0.18634 
    #> Log-likelihood:                                                        -348.8644 
    #> Wald chi2(2):                                                           123.5731 
    #> Prob > chi2:                                                              0.0000 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                        Final Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient EIM Std.Err.   z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant    -1.71537      0.12933 -13.26354   -1.96885    -1.462  < 2e-16 ***
    #> x1           0.64810      0.06323  10.25024    0.52418     0.772  < 2e-16 ***
    #> x2           1.07752      0.19174   5.61978    0.70172     1.453 1.91e-08 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                       Part 3: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                      715 
    #> Pseudo R-squared:                                                        0.24292 
    #> Log pseudolikelihood:                                                  -484.9654 
    #> Wald chi2(2):                                                            261.436 
    #> Prob > chi2:                                                              0.0000 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                     Final Quasi-Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant    -0.26763         0.04539 -5.89578   -0.35660    -0.179 3.73e-09 ***
    #> x1           0.36198         0.02454 14.75144    0.31389     0.410  < 2e-16 ***
    #> x2           0.58505         0.07969  7.34179    0.42886     0.741 2.11e-13 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>     Three-part regression - binary logit , binary probit + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #> Data type:                                                       Cross-sectional 
    #> Pseudo R-squared:                                                        0.38917 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> --------------------------------------------------------------------------------

### 1.1 Partial Effects (`fracreg.pe`)

Raw coefficients from fractional models are not directly interpretable
as marginal effects. Use
[`fracreg.pe()`](https://sulmanolieko.github.io/fracreg/reference/fracreg.pe.md)
to compute the Average Partial Effects (APE) and Conditional Partial
Effects (CPE) using the analytical delta method.

``` r
### Empirical 401(k) Examples 
data("fracreg_k401k") 
y <- fracreg_k401k$prate 
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age,  
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole) 
 
m <- fracreg(y, X, type="1P", linkfrac="logit") 
```

Toggle to see the output

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
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> --------------------------------------------------------------------------------

``` r
fracreg.pe(m) 
```

Toggle to see the output

    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                             Average partial effects 
    #> -------------------------------------------------------------------------------- 
    #>                           Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #>          Estimate Std. Error z value Pr(>|z|)    
    #> mrate   1.018e-01  1.456e-02   6.989 2.77e-12 ***
    #> age     2.980e-03  5.293e-04   5.630 1.80e-08 ***
    #> totemp -8.736e-07  3.279e-07  -2.665  0.00771 ** 
    #> sole    3.635e-02  8.515e-03   4.270 1.96e-05 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> --------------------------------------------------------------------------------

``` r
 
### Simulated Examples 
 
N <- 250 
u <- rnorm(N) 
 
X <- cbind(rnorm(N),rnorm(N)) 
dimnames(X)[[2]] <- c("X1","X2") 
 
ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u)) 
y <- rbeta(N,ym*20,20*(1-ym)) 
y[y > 0.9] <- 1 
 
#Computing average partial effects for a logit fractional response model 
res <- fracreg(y,X,linkfrac="logit",table=FALSE) 
fracreg.pe(res) 
```

Toggle to see the output

    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                             Average partial effects 
    #> -------------------------------------------------------------------------------- 
    #>                           Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #>    Estimate Std. Error z value Pr(>|z|)    
    #> X1 0.161653   0.009584   16.87   <2e-16 ***
    #> X2 0.165141   0.012853   12.85   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> --------------------------------------------------------------------------------

``` r
 
#Computing average partial effects for a binary logit + fractional probit 
#two-part model 
res <- fracreg(y,X,linkbin="logit",linkfrac="probit",type="2P",inf=1,table=FALSE) 
fracreg.pe(res) 
```

Toggle to see the output

    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                             Average partial effects 
    #> -------------------------------------------------------------------------------- 
    #>               Binary logit + Fractional probit two-part regression 
    #> -------------------------------------------------------------------------------- 
    #>    Estimate Std. Error z value Pr(>|z|)    
    #> X1  0.08479    0.01192   7.115 1.12e-12 ***
    #> X2  0.09068    0.01249   7.259 3.89e-13 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> --------------------------------------------------------------------------------

``` r
 
#Computing conditional partial effects for X2 in the logit component 
#of a two-part fractional response model, with the covariates evaluated 
#at their median values 
res <- fracreg(y,X,linkfrac="logit",type="2Pfrac",inf=1,table=FALSE) 
fracreg.pe(res,APE=FALSE,CPE=TRUE,at="median",which.x="X2") 
```

Toggle to see the output

    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                           Conditional partial effects 
    #> -------------------------------------------------------------------------------- 
    #>               Fractional logit component of a two-part regression 
    #> -------------------------------------------------------------------------------- 
    #>    Estimate Std. Error z value Pr(>|z|)    
    #> X2  0.17106    0.01895   9.026   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> -------------------------------------------------------------------------------- 
    #> 
    #> Note: covariates evaluated at median (or mode, for dummies) values

``` r
 
#Computing average partial effects for a three-part double-inflated model 
y3p <- y 
y3p[1:20] <- 0 
y3p[21:40] <- 1 
res3p <- fracreg(y3p,X,linkbin=c("logit","probit"),linkfrac="logit",type="3P",table=FALSE) 
fracreg.pe(res3p) 
```

Toggle to see the output

    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                             Average partial effects 
    #> -------------------------------------------------------------------------------- 
    #>     Three-part regression - binary logit , binary probit + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #>    Estimate Std. Error z value Pr(>|z|)    
    #> X1  0.15535    0.01440  10.790  < 2e-16 ***
    #> X2  0.11405    0.01854   6.153 7.62e-10 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:58 
    #> --------------------------------------------------------------------------------

------------------------------------------------------------------------

## 3. Hypothesis Testing and Specification Diagnostics

`fracreg` includes state-of-the-art specification tests to validate your
model’s functional form and link function assumptions.

### 3.1 Generalised Goodness-Of-Functional-Form (`fracreg.ggoff`)

The GGOFF test tests whether the chosen link function is adequate for
the data. A significant result suggests the link may be misspecified.

``` r
### Empirical 401(k) Examples 
data("fracreg_k401k") 
y <- fracreg_k401k$prate 
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age,  
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole) 
 
m <- fracreg(y, X, type="1P", linkfrac="logit") 
```

Toggle to see the output

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
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

``` r
fracreg.ggoff(m) 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    GGOFF test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #>            Statistic p-value   
    #> GOFF1 - LM     8.838 0.00295 **
    #> GOFF2 - LM     9.828 0.00172 **
    #> GGOFF - LM    10.351 0.00565 **
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

``` r
 
### Simulated Examples 
 
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
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    GGOFF test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #>              Statistic p-value
    #> GOFF1 - LM       1.256   0.262
    #> GOFF1 - Wald     1.274   0.259
    #> GOFF2 - LM       1.513   0.219
    #> GOFF2 - Wald     1.401   0.237
    #> GGOFF - LM       1.612   0.447
    #> GGOFF - Wald     1.336   0.513
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

``` r
 
#Testing the probit specification of the binary component of a two-part fractional 
#regression model using a LR-based GGOFF test 
res <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1,table=FALSE) 
fracreg.ggoff(res,"LR") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    GGOFF test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Binary probit component of a two-part regression 
    #> -------------------------------------------------------------------------------- 
    #>            Statistic p-value
    #> GOFF1 - LR     0.012   0.914
    #> GOFF2 - LR     0.017   0.895
    #> GGOFF - LR     0.024   0.988
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

### 3.2 RESET Test (`fracreg.reset`)

The RESET test detects general functional form misspecification by
testing whether powers of the fitted values have explanatory power.
Testing $`H_0: \gamma = 0`$ provides a robust diagnostic check.

``` r
### Empirical 401(k) Examples 
data("fracreg_k401k") 
y <- fracreg_k401k$prate 
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age,  
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole) 
 
m <- fracreg(y, X, type="1P", linkfrac="logit") 
```

Toggle to see the output

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
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

``` r
fracreg.reset(m) 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    RESET test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #>       Statistic p-value   
    #> LM(3)     10.29 0.00583 **
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

``` r
 
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
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    RESET test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #>         Statistic p-value  
    #> LM(2)       4.832  0.0279 *
    #> Wald(2)     5.816  0.0159 *
    #> LM(3)       4.866  0.0878 .
    #> Wald(3)     6.870  0.0322 *
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

``` r
 
#Testing the probit specification of the binary component of a two-part fractional 
#regression model using LR-based RESET tests with quadratic and cubic fitted  
#powers of the linear predictor 
res <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1,table=FALSE) 
fracreg.reset(res,3,"LR") 
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    RESET test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Binary probit component of a two-part regression 
    #> -------------------------------------------------------------------------------- 
    #>       Statistic p-value
    #> LR(3)     0.976   0.614
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

### 3.3 P-Test for Non-Nested Models (`fracreg.ptest`)

You can compare non-nested models (e.g., `logit` vs. `cloglog` link)
using the Davidson-MacKinnon (1981) P-test.

``` r
### Empirical 401(k) Examples 
data("fracreg_k401k") 
y <- fracreg_k401k$prate 
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age,  
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole) 
 
m1 <- fracreg(y, X, type="1P", linkfrac="logit") 
```

Toggle to see the output

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
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

``` r
m2 <- fracreg(y, X, type="1P", linkfrac="probit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                           Fractional probit regression 
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
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

``` r
fracreg.ptest(m1, m2) 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                      P test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit regression 
    #> H1: Fractional probit regression 
    #> -------------------------------------------------------------------------------- 
    #>      Statistic p-value
    #> Wald    -1.483   0.138
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional probit regression 
    #> H1: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #>      Statistic p-value   
    #> Wald     2.754 0.00595 **
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

``` r
 
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
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                      P test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit regression 
    #> H1: Fractional loglog regression 
    #> -------------------------------------------------------------------------------- 
    #>    Statistic p-value
    #> LM     0.239   0.625
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional loglog regression 
    #> H1: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #>    Statistic p-value   
    #> LM     8.502 0.00355 **
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

``` r
 
#Testing a logit one-part fractional response model versus a binary logit + 
#fractional probit two-part model using a Wald version of the P test 
res1 <- fracreg(y,X,linkfrac="logit",table=FALSE) 
res2 <- fracreg(y,X,linkbin="logit",linkfrac="probit",type="2P",inf=1,table=FALSE) 
fracreg.ptest(res1,res2,"Wald") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                      P test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit regression 
    #> H1: Binary logit + Fractional probit two-part regression 
    #> -------------------------------------------------------------------------------- 
    #>      Statistic p-value
    #> Wald     0.207   0.836
    #> -------------------------------------------------------------------------------- 
    #> H0: Binary logit + Fractional probit two-part regression 
    #> H1: Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #>      Statistic p-value    
    #> Wald     13.95  <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:49:59 
    #> --------------------------------------------------------------------------------

------------------------------------------------------------------------

## 4. Endogeneity & Heteroscedasticity (`fracreghet`)

When you suspect that one of your covariates is endogenous, or that the
variance is heteroscedastic,
[`fracreghet()`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.md)
provides instrumental variable correction. It natively supports IV
through a **Control Function (CF)** approach or **GMM** estimation.

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
fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], type="QMLxv", link="logit") 
```

Toggle to see the output

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
    #>          Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)    
    #> Constant    -0.10561         0.75224 -0.14040   -1.57997     1.369    0.888    
    #> mrate        3.72138         0.68952  5.39703    2.36994     5.073 6.78e-08 ***
    #> ltotemp     -0.07009         0.05348 -1.31060   -0.17491     0.035    0.190    
    #> vhat        -2.79515         0.70026 -3.99157   -4.16764    -1.423 6.56e-05 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #>                                  Reduced form: 
    #> -------------------------------------------------------------------------------- 
    #>             Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
    #> Z_INTERCEPT     0.95046         0.09550  9.95211    0.76327     1.138  < 2e-16
    #> Z_age           0.01146         0.00231  4.95997    0.00693     0.016 7.05e-07
    #> Z_ltotemp      -0.05534         0.01421 -3.89528   -0.08318    -0.027 9.81e-05
    #>                
    #> Z_INTERCEPT ***
    #> Z_age       ***
    #> Z_ltotemp   ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:00 
    #> --------------------------------------------------------------------------------

``` r

# Compute the same QMLxv estimator reporting Odds Ratios with 90% confidence intervals
fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], type="QMLxv", link="logit", or=TRUE, level=0.90)
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>         Fractional logit regression with heteroscedasticity/endogeneity 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                 QMLxv 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1534 
    #> Standard errors:                                                          robust 
    #> Wald chi2(6):                                                       2243425.7812 
    #> Prob > chi2:                                                              0.0000 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                   Final Quasi-Maximum Likelihood xv estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Odds Ratio Robust Std.Err.  z value [90% Conf. Interval] Pr(>|z|)    
    #> Constant    0.89977         0.67684 -0.14040    0.26108     3.101    0.888    
    #> mrate      41.32156        28.49224  5.39703   13.29272   128.452 6.78e-08 ***
    #> ltotemp     0.93231         0.04986 -1.31060    0.85380     1.018    0.190    
    #> vhat        0.06111         0.04279 -3.99157    0.01931     0.193 6.56e-05 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #>                                  Reduced form: 
    #> -------------------------------------------------------------------------------- 
    #>             Odds Ratio Robust Std.Err.   z value [90% Conf. Interval] Pr(>|z|)
    #> Z_INTERCEPT   2.586889        0.247056  9.952107   2.210829     3.027  < 2e-16
    #> Z_age         1.011524        0.002337  4.959966   1.007688     1.015 7.05e-07
    #> Z_ltotemp     0.946168        0.013441 -3.895285   0.924315     0.969 9.81e-05
    #>                
    #> Z_INTERCEPT ***
    #> Z_age       ***
    #> Z_ltotemp   ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:00 
    #> --------------------------------------------------------------------------------

``` r
 
### Simulated Examples 
 
set.seed(123) 
N <- 1000 
x1 <- rnorm(N) 
 
# Simulating an endogenous variable (var.endog) and an instrument (z1) 
z1 <- rnorm(N) 
u <- 0.5 * z1 + rnorm(N) 
var.endog <- 0.8 * z1 + u 
y_endog <- exp(0.5 * x1 + 1.2 * var.endog + u) / (1 + exp(0.5 * x1 + 1.2 * var.endog + u)) 
 
# Avoid exact 0 or 1 boundaries for some estimators 
y_endog[y_endog <= 0] <- 0.01 
y_endog[y_endog >= 1] <- 0.99 
 
X <- cbind(x1 = x1, var.endog = var.endog) 
Z <- cbind(x1 = x1, z1 = z1) 
 
# Exogeneity (assuming var.endog is exogenous for comparison), GMMx estimator 
fracreghet(y = y_endog, x = X, type = "GMMx", link = "logit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>         Fractional logit regression with heteroscedasticity/endogeneity 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                  GMMx 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1000 
    #> Standard errors:                                                          robust 
    #> Wald chi2(2):                                                         42761.0276 
    #> Prob > chi2:                                                              0.0000 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                               Final GMMx estimates 
    #> -------------------------------------------------------------------------------- 
    #>           Coefficient Robust Std.Err.   z value [95% Conf. Interval] Pr(>|z|)
    #> Constant    9.132e-02       1.536e-02 5.947e+00  6.123e-02     0.121 2.73e-09
    #> x1          4.608e-01       1.538e-02 2.995e+01  4.306e-01     0.491  < 2e-16
    #> var.endog   1.808e+00       9.021e-03 2.004e+02  1.790e+00     1.826  < 2e-16
    #>              
    #> Constant  ***
    #> x1        ***
    #> var.endog ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:00 
    #> --------------------------------------------------------------------------------

``` r
 
# Endogeneity, GMMz estimator (does not require reduced form for endog) 
fracreghet(y = y_endog, x = X, z = Z, type = "GMMz", link = "logit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>         Fractional logit regression with heteroscedasticity/endogeneity 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                  GMMz 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1000 
    #> Standard errors:                                                          robust 
    #> Wald chi2(2):                                                         14903.8562 
    #> Prob > chi2:                                                              0.0000 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                               Final GMMz estimates 
    #> -------------------------------------------------------------------------------- 
    #>           Coefficient Robust Std.Err.   z value [95% Conf. Interval] Pr(>|z|)
    #> Constant      0.15277         0.02018   7.56931    0.11321     0.192 3.75e-14
    #> x1            0.47947         0.02051  23.37680    0.43927     0.520  < 2e-16
    #> var.endog     1.61252         0.01445 111.61802    1.58421     1.641  < 2e-16
    #>              
    #> Constant  ***
    #> x1        ***
    #> var.endog ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:00 
    #> --------------------------------------------------------------------------------

``` r
 
# Endogeneity, GMMxv estimator (assumes linear reduced form for var.endog) 
fracreghet(y = y_endog, x = X, z = Z, var.endog = var.endog, type = "GMMxv", link = "logit") 
#> Warning in dgamma(y, 1/disp, scale = mu * disp, log = TRUE): NaNs produced
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>         Fractional logit regression with heteroscedasticity/endogeneity 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                 GMMxv 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1000 
    #> Standard errors:                                                          robust 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                              Final GMMxv estimates 
    #> -------------------------------------------------------------------------------- 
    #>           Coefficient Robust Std.Err.   z value [95% Conf. Interval] Pr(>|z|)
    #> Constant     -0.01262         0.01859  -0.67857   -0.04905     0.024    0.497
    #> x1            0.48705         0.01884  25.85287    0.45012     0.524   <2e-16
    #> var.endog     1.59737         0.01339 119.33166    1.57113     1.624   <2e-16
    #> vhat          0.60263         0.01339  45.01988    0.57640     0.629   <2e-16
    #>              
    #> Constant     
    #> x1        ***
    #> var.endog ***
    #> vhat      ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #>                                  Reduced form: 
    #> -------------------------------------------------------------------------------- 
    #>             Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
    #> Z_INTERCEPT    -0.02093         0.03082 -0.67933   -0.08133     0.039    0.497
    #> Z_x1           -0.02149         0.03132 -0.68612   -0.08288     0.040    0.493
    #> Z_z1            1.32751         0.02949 45.01988    1.26971     1.385   <2e-16
    #>                
    #> Z_INTERCEPT    
    #> Z_x1           
    #> Z_z1        ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:00 
    #> --------------------------------------------------------------------------------

``` r
 
# Endogeneity, QMLxv control function approach 
fracreghet(y = y_endog, x = X, z = Z, var.endog = var.endog, type = "QMLxv", link = "logit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>         Fractional logit regression with heteroscedasticity/endogeneity 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                 QMLxv 
    #> Data type:                                                       Cross-sectional 
    #> Number of observations:                                                     1000 
    #> Standard errors:                                                          robust 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                   Final Quasi-Maximum Likelihood xv estimates 
    #> -------------------------------------------------------------------------------- 
    #>           Coefficient Robust Std.Err.   z value [95% Conf. Interval] Pr(>|z|)
    #> Constant     -0.01262         0.01859  -0.67857   -0.04905     0.024    0.497
    #> x1            0.48705         0.01884  25.85287    0.45012     0.524   <2e-16
    #> var.endog     1.59737         0.01339 119.33166    1.57113     1.624   <2e-16
    #> vhat          0.60263         0.01339  45.01988    0.57640     0.629   <2e-16
    #>              
    #> Constant     
    #> x1        ***
    #> var.endog ***
    #> vhat      ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #>                                  Reduced form: 
    #> -------------------------------------------------------------------------------- 
    #>             Coefficient Robust Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
    #> Z_INTERCEPT    -0.02093         0.03082 -0.67933   -0.08133     0.039    0.497
    #> Z_x1           -0.02149         0.03132 -0.68612   -0.08288     0.040    0.493
    #> Z_z1            1.32751         0.02949 45.01988    1.26971     1.385   <2e-16
    #>                
    #> Z_INTERCEPT    
    #> Z_x1           
    #> Z_z1        ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:00 
    #> --------------------------------------------------------------------------------

### 4.1 Partial Effects for Endogenous Models (`fracreghet.pe`)

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
res_emp <- fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], type="QMLxv", link="logit", table=FALSE) 
fracreghet.pe(res_emp, which.x="mrate")
```

Toggle to see the output

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
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:01 
    #> --------------------------------------------------------------------------------

``` r

### Simulated Examples 
 
N <- 250 
u <- rnorm(N) 
 
X <- cbind(rnorm(N),rnorm(N)) 
dimnames(X)[[2]] <- c("X1","X2") 
 
Z <- cbind(rnorm(N),rnorm(N),rnorm(N)) 
dimnames(Z)[[2]] <- c("Z1","Z2","Z3") 
 
y <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u)) 
 
res <- fracreghet(y,X,type="GMMx",table=FALSE) 
 
#Smearing estimator of average partial effects for variable X1 
fracreghet.pe(res,which.x="X1") 
```

Toggle to see the output

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
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:01 
    #> --------------------------------------------------------------------------------

``` r
 
#Naive estimator of conditional partial effects for all covariates, 
#which are evaluated at X1=1 and X2=-1 
fracreghet.pe(res,smearing=FALSE,APE=FALSE,CPE=TRUE,at=c(1,-1)) 
```

Toggle to see the output

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
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:01 
    #> -------------------------------------------------------------------------------- 
    #> 
    #> Note: covariates evaluated at the following values:
    #> 
    #> X1 X2 
    #>  1 -1

### 4.2 RESET Test for Endogenous Models (`fracreghet.reset`)

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
res_emp <- fracreghet(y_adj, X_het, type="GMMx", link="logit", table=FALSE) 
fracreghet.reset(res_emp)
```

Toggle to see the output

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
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:01 
    #> --------------------------------------------------------------------------------

``` r

### Simulated Examples 
 
N <- 250 
u <- rnorm(N) 
 
X <- cbind(rnorm(N),rnorm(N)) 
dimnames(X)[[2]] <- c("X1","X2") 
 
Z <- cbind(rnorm(N),rnorm(N),rnorm(N)) 
dimnames(Z)[[2]] <- c("Z1","Z2","Z3") 
 
y <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u)) 
 
res <- fracreghet(y,X,type="GMMx",table=FALSE) 
 
#LM and Wald versions of the RESET test, based on 1 or 2 fitted powers of xb 
fracreghet.reset(res,2:3,c("Wald","LM")) 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    RESET test 
    #> -------------------------------------------------------------------------------- 
    #>                           Fractional logit regression 
    #> -------------------------------------------------------------------------------- 
    #> H0: Estimator: GMMx 
    #> -------------------------------------------------------------------------------- 
    #>         Statistic p-value  
    #> LM(2)       0.216  0.6424  
    #> Wald(2)     0.180  0.6711  
    #> LM(3)       3.481  0.1754  
    #> Wald(3)     6.273  0.0434 *
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:01 
    #> --------------------------------------------------------------------------------

------------------------------------------------------------------------

## 5. Panel Data Fractional Models (`fracregpd`)

For longitudinal or panel data where unobserved heterogeneity is a
concern,
[`fracregpd()`](https://sulmanolieko.github.io/fracreg/reference/fracregpd.md)
provides fixed-T panel estimators.

``` r
### Empirical 401(k) Examples 
data("fracreg_k401k") 
y <- fracreg_k401k$prate 
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age,  
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole) 
 
# Artificial panel data structure for demonstration 
N_emp <- nrow(X) 
id_emp <- rep(1:(N_emp/2), each=2) 
time_emp <- rep(1:2, times=N_emp/2) 
fracregpd(id_emp, time_emp, y, X, type="QMLcre", link="probit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>            Fractional probit (correlated random effects)  regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                QMLcre 
    #> Data type:                                                                 Panel 
    #> Exogeneity:                                                                 TRUE 
    #> Use first lag of instruments:                                              FALSE 
    #> Standard errors:                                                         cluster 
    #> -------------------------------------------------------------------------------- 
    #> Number of observations:                                                     1534 
    #> Number of groups:                                                            767 
    #> Obs per group:                                                                 2 
    #> Log pseudolikelihood:                                                  -554.0205 
    #> Wald chi2(9):                                                          3191.0501 
    #> Prob > chi2:                                                              0.0000 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>       Final (Correlated Random Effects) Quasi-Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>                Coefficient Cluster Std.Err.    z value [95% Conf. Interval]
    #> mrate            3.915e-01        7.007e-02  5.588e+00  2.542e-01     0.529
    #> age              1.386e-02        3.543e-03  3.911e+00  6.912e-03     0.021
    #> totemp          -5.250e-06        2.649e-06 -1.982e+00 -1.044e-05     0.000
    #> sole             2.378e-01        6.038e-02  3.938e+00  1.194e-01     0.356
    #> INTERCEPT_mean   6.235e-01        5.885e-02  1.060e+01  5.082e-01     0.739
    #> mrate_mean       5.726e-02        6.343e-02  9.027e-01 -6.706e-02     0.182
    #> age_mean         1.544e-03        4.233e-03  3.649e-01 -6.751e-03     0.010
    #> totemp_mean      1.223e-06        2.946e-06  4.152e-01 -4.551e-06     0.000
    #> sole_mean       -7.157e-02        8.424e-02 -8.495e-01 -2.367e-01     0.094
    #>                Pr(>|z|)    
    #> mrate          2.30e-08 ***
    #> age            9.20e-05 ***
    #> totemp           0.0475 *  
    #> sole           8.22e-05 ***
    #> INTERCEPT_mean  < 2e-16 ***
    #> mrate_mean       0.3667    
    #> age_mean         0.7152    
    #> totemp_mean      0.6780    
    #> sole_mean        0.3956    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:02 
    #> --------------------------------------------------------------------------------

``` r
 
### Simulated Examples 
 
set.seed(123) 
# Simulating Panel Data 
N <- 100 
T_periods <- 5 
id <- rep(1:N, each = T_periods) 
time <- rep(1:T_periods, times = N) 
x_panel <- rnorm(N * T_periods) 
 
# Unobserved individual effect (CRE) 
c_i <- rep(rnorm(N), each = T_periods)  
y_panel <- exp(x_panel + c_i) / (1 + exp(x_panel + c_i)) 
 
X <- cbind(x_panel = x_panel) 
 
# Endogenous variable and instrument simulation 
z_panel <- rnorm(N * T_periods) 
u_panel <- 0.5 * z_panel + rnorm(N * T_periods) 
var_endog <- 0.8 * z_panel + u_panel 
y_endog <- exp(x_panel + 1.2 * var_endog + c_i + u_panel) /  
             (1 + exp(x_panel + 1.2 * var_endog + c_i + u_panel)) 
 
X_endog <- cbind(x_panel = x_panel, var_endog = var_endog) 
Z_inst <- cbind(x_panel = x_panel, z_panel = z_panel) 
 
# Estimate a Correlated Random Effects (CRE) Model 
fracregpd(id=id, time=time, y=y_panel, x=X, type="QMLcre", link="probit") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>            Fractional probit (correlated random effects)  regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                QMLcre 
    #> Data type:                                                                 Panel 
    #> Exogeneity:                                                                 TRUE 
    #> Use first lag of instruments:                                              FALSE 
    #> Standard errors:                                                         cluster 
    #> -------------------------------------------------------------------------------- 
    #> Number of observations:                                                      500 
    #> Number of groups:                                                            100 
    #> Obs per group:                                                                 5 
    #> Log pseudolikelihood:                                                  -313.6436 
    #> Wald chi2(3):                                                          2156.6004 
    #> Prob > chi2:                                                              0.0000 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>       Final (Correlated Random Effects) Quasi-Maximum Likelihood estimates 
    #> -------------------------------------------------------------------------------- 
    #>                Coefficient Cluster Std.Err.  z value [95% Conf. Interval]
    #> x_panel            0.52902          0.01143 46.29965    0.50662     0.551
    #> INTERCEPT_mean    -0.01246          0.04901 -0.25433   -0.10851     0.084
    #> x_panel_mean      -0.17409          0.13491 -1.29040   -0.43850     0.090
    #>                Pr(>|z|)    
    #> x_panel          <2e-16 ***
    #> INTERCEPT_mean    0.799    
    #> x_panel_mean      0.197    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:02 
    #> --------------------------------------------------------------------------------

``` r
 
# Exogeneity, no lags, no time dummies, clustered standard errors, GMMbgw estimator 
fracregpd(id=id, time=time, y=y_panel, x=X, type="GMMbgw") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                           Fractional logit  regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                GMMbgw 
    #> Data type:                                                                 Panel 
    #> Exogeneity:                                                                 TRUE 
    #> Use first lag of instruments:                                              FALSE 
    #> Standard errors:                                                         cluster 
    #> -------------------------------------------------------------------------------- 
    #> Number of observations:                                                      500 
    #> Number of groups:                                                            100 
    #> Obs per group:                                                                 5 
    #> Wald chi2(1):                                               4.88095376809867e+31 
    #> Prob > chi2:                                                              0.0000 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                             Final GMM bgw estimates 
    #> -------------------------------------------------------------------------------- 
    #>         Coefficient Cluster Std.Err.   z value [95% Conf. Interval] Pr(>|z|)
    #> x_panel   1.000e+00        1.431e-16 6.986e+15  1.000e+00         1   <2e-16
    #>            
    #> x_panel ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:02 
    #> --------------------------------------------------------------------------------

``` r

# Estimate the GMMww estimator with odds ratios and 99% confidence intervals
fracregpd(id=id, time=time, y=y_panel, x=X, type="GMMww", or=TRUE, level=0.99)
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                           Fractional logit  regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                 GMMww 
    #> Data type:                                                                 Panel 
    #> Exogeneity:                                                                 TRUE 
    #> Use first lag of instruments:                                               TRUE 
    #> Standard errors:                                                         cluster 
    #> -------------------------------------------------------------------------------- 
    #> Number of obs (initial):                                                     500 
    #> Number of observations:                                                      400 
    #> Number of groups (initial):                                                  100 
    #> Number of groups:                                                            100 
    #> Obs per group:                                                                 4 
    #> Wald chi2(1):                                               8.12351101825493e+32 
    #> Prob > chi2:                                                              0.0000 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                              Final GMM ww estimates 
    #> -------------------------------------------------------------------------------- 
    #>         Odds Ratio Cluster Std.Err.   z value [99% Conf. Interval] Pr(>|z|)    
    #> x_panel  2.718e+00        2.592e-16 1.049e+16  2.718e+00     2.718   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:02 
    #> --------------------------------------------------------------------------------

``` r
 
# Lagged covariates and instruments, robust standard errors, GMMww estimator 
fracregpd(id=id, time=time, y=y_panel, x=X, lags=TRUE, type="GMMww", var.type="robust") 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                           Fractional logit  regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                 GMMww 
    #> Data type:                                                                 Panel 
    #> Exogeneity:                                                                 TRUE 
    #> Use first lag of instruments:                                               TRUE 
    #> Standard errors:                                                          robust 
    #> -------------------------------------------------------------------------------- 
    #> Number of obs (initial):                                                     500 
    #> Number of observations:                                                      400 
    #> Number of groups (initial):                                                  100 
    #> Number of groups:                                                            100 
    #> Obs per group:                                                                 4 
    #> Wald chi2(1):                                               1.28107689657633e+32 
    #> Prob > chi2:                                                              0.0000 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                              Final GMM ww estimates 
    #> -------------------------------------------------------------------------------- 
    #>         Coefficient Robust Std.Err.   z value [95% Conf. Interval] Pr(>|z|)    
    #> x_panel   1.000e+00       8.835e-17 1.132e+16  1.000e+00         1   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:02 
    #> --------------------------------------------------------------------------------

``` r
 
# Endogeneity, time dummies, GMMpfe estimator 
fracregpd(id=id, time=time, y=y_endog, x=X_endog, z=Z_inst,  
          x.exogenous=FALSE, type="GMMpfe", tdummies=TRUE) 
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>               Fractional logit (pooled fixed effects)  regression 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                GMMpfe 
    #> Data type:                                                                 Panel 
    #> Exogeneity:                                                                FALSE 
    #> Use first lag of instruments:                                              FALSE 
    #> Standard errors:                                                         cluster 
    #> -------------------------------------------------------------------------------- 
    #> Number of observations:                                                      500 
    #> Number of groups:                                                            100 
    #> Obs per group:                                                                 5 
    #> Wald chi2(6):                                                          7156.2911 
    #> Prob > chi2:                                                              0.0000 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>       Final (Pooled Fixed Effects) Generalized Method of Moments estimates 
    #> -------------------------------------------------------------------------------- 
    #>           Coefficient Cluster Std.Err.  z value [95% Conf. Interval] Pr(>|z|)
    #> x_panel       0.99602          0.03117 31.95739    0.93493     1.057   <2e-16
    #> var_endog     1.57689          0.02524 62.47420    1.52742     1.626   <2e-16
    #> time.2       -0.15363          0.09752 -1.57541   -0.34476     0.038   0.1152
    #> time.3       -0.04402          0.09147 -0.48128   -0.22331     0.135   0.6303
    #> time.4       -0.10609          0.09379 -1.13110   -0.28992     0.078   0.2580
    #> time.5       -0.19163          0.08917 -2.14916   -0.36640    -0.017   0.0316
    #>              
    #> x_panel   ***
    #> var_endog ***
    #> time.2       
    #> time.3       
    #> time.4       
    #> time.5    *  
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-06 17:50:03 
    #> --------------------------------------------------------------------------------

------------------------------------------------------------------------

## Conclusion

With `fracreg`, `fracreghet`, and `fracregpd`, you have a complete,
mathematically robust toolkit for modelling fractional responses bounded
between $`[0,1]`$, regardless of inflation, endogeneity, or unobserved
panel effects.

## Acknowledgements

This package builds upon, consolidates, and modernises the fractional
regression frameworks originally implemented in the `frm`, `frmhet`, and
`frmpd` R packages developed by Joaquim J.S. Ramalho. As those original
packages have been deprecated and removed from the active CRAN
repository, `fracreg` serves as an actively maintained successor,
ensuring these econometric tools remain available to the R community.

## References

- **Ramalho, J. J. S. (2022).** *frm: Fractional Regression Models*. R
  package. Formerly available on CRAN, currently archived.
- **Ramalho, J. J. S. (2023).** *frmhet: Fractional Regression Models
  under Heterogeneity*. R package. Formerly available on CRAN, currently
  archived.
- **Ramalho, J. J. S. (2023).** *frmpd: Fractional Regression Models for
  Panel Data*. R package. Formerly available on CRAN, currently
  archived.
- **Davidson, R. and MacKinnon, J.G. (1981).** “Several tests for model
  specification in the presence of alternative hypotheses”,
  *Econometrica*, 49(3), 781-793.
- **Fang, K., & Ma, S. (2013).** “Three-part model for fractional
  response variables with application to Chinese household health
  insurance coverage”, *Journal of Applied Statistics*, 40(5), 925-940.
- **Papke, L. E. and Wooldridge, J. M. (1996).** “Econometric methods
  for fractional response variables with an application to 401(k) plan
  participation rates”, *Journal of Applied Econometrics*, 11(6),
  619-632.
- **Papke, L. and Wooldridge, J.M. (2008).** “Panel data methods for
  fractional response variables with an application to test pass rates”,
  *Journal of Econometrics*, 145(1-2), 121-233.
- **Pregibon, D. (1980).** “Goodness of Link Tests for Generalized
  Linear Models”, *Journal of the Royal Statistical Society: Series C
  (Applied Statistics)*, 29(1), 15-24.
- **Ramalho, E. A., & Ramalho, J. J. S. (2017).** “Moment-based
  estimation of nonlinear regression models with boundary outcomes and
  endogeneity, with applications to nonnegative and fractional
  responses”, *Econometric Reviews*, 36(4), 397-420.
- **Ramalho, E.A., Ramalho, J.J.S. and Murteira, J.M.R. (2011).**
  “Alternative estimating and testing empirical strategies for
  fractional response models”, *Journal of Economic Surveys*, 25(1),
  19-68.
- **Ramalho, E.A., Ramalho, J.J.S. and Murteira, J.M.R. (2014).** “A
  generalized goodness-of-functional form test for binary and fractional
  response models”, *Manchester School*, 82(4), 488-507.
- **Ramsey, J.B. (1969).** “Tests for Specification Errors in Classical
  Linear Least-Squares Regression Analysis”, *Journal of the Royal
  Statistical Society: Series B (Methodological)*, 31(2), 350-371.

For more information, please visit the [package
website](https://sulmanolieko.github.io/fracreg/) or file an issue on
[GitHub](https://github.com/SulmanOlieko/fracreg/issues).

To cite this package in your research:

``` r
citation("fracreg")
```
