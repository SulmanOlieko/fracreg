
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fracreg: Fractional Response Models in R

<img src="man/figures/logo.png" align="right" height="139" alt="fracreg logo" />

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![](https://img.shields.io/badge/R-%3E=%203.5.0-blue.svg)](https://cran.r-project.org/)
[![](https://img.shields.io/github/last-commit/SulmanOlieko/fracreg.svg)](https://github.com/SulmanOlieko/fracreg/commits/main)
[![](https://img.shields.io/github/issues/SulmanOlieko/fracreg.svg)](https://github.com/SulmanOlieko/fracreg/issues)
[![](https://img.shields.io/badge/license-GPL--3-blue)](https://github.com/SulmanOlieko/fracreg)
[![CRAN Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/fracreg)](https://cran.r-project.org/package=fracreg)

An R package for fitting and testing specifications of fractional
regression models. It handles fractional, univariate proportions, and
bounded data ($0 \le y \le 1$), and provides estimators for one-part,
two-part (hurdle), and three-part (double-inflated) models.

## Features

`fracreg` is handles:

1.  **Univariate Fractional Models (`fracreg`)**: Fit standard 1-part
    models, hurdle 2-part models (mass at 0 or 1), and double-inflated
    3-part models.
2.  **Panel Data Models (`fracregpd`)**: Estimate fractional models with
    fixed-T longitudinal data, supporting Correlated Random Effects
    (CRE) to handle unobserved individual heterogeneity.
3.  **Endogeneity & Heteroscedasticity (`fracreghet`)**: Correct for
    endogenous covariates using Instrumental Variables (IV) via Control
    Function and GMM approaches.
4.  **Diagnostic Testing**: Includes generalized
    goodness-of-functional-form (`fracreg.ggoff`), RESET
    (`fracreg.reset`), and non-nested P-tests (`fracreg.ptest`).
5.  **Analytical Partial Effects**: Compute exact marginal effects for
    all model types (including 3-part composite effects) via the Delta
    method using `fracreg.pe()`.

------------------------------------------------------------------------

## Installation

You can install the released version of `fracreg` from CRAN with:

``` r
install.packages("fracreg")
library("fracreg")
```

Or install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("SulmanOlieko/fracreg")
```

------------------------------------------------------------------------

## Usage Examples

This guide walks you through empirical examples demonstrating each
estimator using the 401(k) dataset.

``` r
library(fracreg)
```

### 1. The Core Estimator (`fracreg`)

The main `fracreg()` function fits cross-sectional models.

#### 1.1 Data Description and Preparation

We use the built-in `fracreg_k401k` dataset, which is the canonical
firm-level 401(k) plan participation data used in **Papke and Wooldridge
(1996)** (*“Econometric methods for fractional response variables with
an application to 401(k) plan participation rates”*, Journal of Applied
Econometrics). The dataset contains 1,534 observations of 401(k) plans.

The primary variables we use are: - `prate`: The plan participation rate
(the fraction of eligible employees who are active participants). It
strictly falls in the $[0, 1]$ interval and is our dependent variable
($y$). - `mrate`: The firm’s matching rate (the firm’s contribution per
\$1 of employee contribution). - `age`: The age of the 401(k) plan. -
`totemp`: Total number of employees at the firm. - `sole`: A binary
indicator equal to 1 if the 401(k) plan is the sole retirement plan
offered by the firm.

``` r
data("fracreg_k401k")
y <- fracreg_k401k$prate
X <- cbind(
  mrate = fracreg_k401k$mrate, 
  age = fracreg_k401k$age, 
  totemp = fracreg_k401k$totemp, 
  sole = fracreg_k401k$sole
)
head(cbind(y, X))
```

<details>

<summary>

Toggle to see the output
</summary>

    #>          y mrate age totemp sole
    #> [1,] 0.261  0.21   8   8709    0
    #> [2,] 1.000  1.42   6    315    1
    #> [3,] 0.976  0.91  10    275    1
    #> [4,] 1.000  0.42   7    500    0
    #> [5,] 0.825  0.53  28    933    1
    #> [6,] 1.000  1.82   7    143    1

</details>

#### 1.2 The Standard 1-Part Model (1P)

The one-part model assumes a single mechanism governs the entire
distribution, estimating $E(y \mid x) = G(x\beta)$:

``` r
model_1p <- fracreg(
  y = y, 
  x = X, 
  type = "1P", 
  linkfrac = "logit"
)
```

<details>

<summary>

Toggle to see the output
</summary>

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                 Fractional response model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Number of observations:                                                     1534 
    #> Pseudo R-squared:                                                        0.14667 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>            Estimate Std. Error z value Pr(>|z|)    
    #> Constant  9.316e-01  8.408e-02  11.080  < 2e-16 ***
    #> mrate     9.531e-01  1.371e-01   6.951 3.62e-12 ***
    #> age       2.791e-02  4.877e-03   5.723 1.05e-08 ***
    #> totemp   -8.182e-06  3.061e-06  -2.673  0.00751 ** 
    #> sole      3.405e-01  8.066e-02   4.222 2.43e-05 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 23:18:43 
    #> --------------------------------------------------------------------------------

</details>

#### 1.3 Partial Effects (1P)

``` r
pe_1p <- fracreg.pe(model_1p)
```

<details>

<summary>

Toggle to see the output
</summary>

    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                             Average partial effects 
    #> -------------------------------------------------------------------------------- 
    #>                              Fractional logit model 
    #> -------------------------------------------------------------------------------- 
    #>          Estimate Std. Error z value Pr(>|z|)    
    #> mrate   1.018e-01  1.456e-02   6.989 2.77e-12 ***
    #> age     2.980e-03  5.293e-04   5.630 1.80e-08 ***
    #> totemp -8.736e-07  3.279e-07  -2.665  0.00771 ** 
    #> sole    3.635e-02  8.515e-03   4.270 1.96e-05 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 23:18:43 
    #> --------------------------------------------------------------------------------

</details>

#### 1.4 The Two-Part Hurdle Model (2P)

The two-part model separates the decision to participate (binary hurdle)
from the proportion chosen:
$E(y \mid x) = \Pr(y > 0 \mid x) \cdot E(y \mid y > 0, x)$.

``` r
# We use inflation=1 to indicate we are modelling the mass at 100% participation (1s)
model_2p <- fracreg(
  y = y, 
  x = X, 
  type = "2P", 
  inflation = 1,
  linkbin = "logit", 
  linkfrac = "logit"
)
```

<details>

<summary>

Toggle to see the output
</summary>

    #> 
    #> -------------------------------------------------------------------------------- 
    #>            Binary component of a two-part model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Number of observations:                                                     1534 
    #> Pseudo R-squared:                                                         0.1485 
    #> Log-Likelihood:                                                        -938.1759 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>            Estimate Std. Error z value Pr(>|z|)    
    #> Constant -1.396e+00  1.270e-01 -10.991   <2e-16 ***
    #> mrate     9.053e-01  9.699e-02   9.334   <2e-16 ***
    #> age       1.156e-02  6.218e-03   1.858   0.0631 .  
    #> totemp   -1.418e-05  6.324e-06  -2.242   0.0249 *  
    #> sole      8.651e-01  1.131e-01   7.651    2e-14 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>          Fractional component of a two-part model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Number of observations:                                                      852 
    #> Pseudo R-squared:                                                        0.10004 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>            Estimate Std. Error z value Pr(>|z|)    
    #> Constant  7.460e-01  6.850e-02  10.891  < 2e-16 ***
    #> mrate     3.877e-01  9.725e-02   3.987 6.69e-05 ***
    #> age       2.562e-02  4.010e-03   6.390 1.66e-10 ***
    #> totemp   -4.061e-06  3.073e-06  -1.322    0.186    
    #> sole     -1.510e-02  6.556e-02  -0.230    0.818    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                 Two-part model - binary logit + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #> Pseudo R-squared:                                                        0.11243 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 23:18:44 
    #> --------------------------------------------------------------------------------

</details>

``` r
pe_2p <- fracreg.pe(model_2p)
```

<details>

<summary>

Toggle to see the output
</summary>

    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                             Average partial effects 
    #> -------------------------------------------------------------------------------- 
    #>                  Binary logit + Fractional logit two-part model 
    #> -------------------------------------------------------------------------------- 
    #>          Estimate Std. Error z value Pr(>|z|)    
    #> mrate   1.769e-01  1.568e-02  11.284  < 2e-16 ***
    #> age     3.680e-03  1.061e-03   3.469 0.000523 ***
    #> totemp -2.633e-06  1.066e-06  -2.470 0.013507 *  
    #> sole    1.425e-01  1.805e-02   7.896 2.89e-15 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 23:18:44 
    #> --------------------------------------------------------------------------------

</details>

#### 1.5 The Three-Part Double-Inflated Model (3P)

The three-part model accommodates clustering at both $0$ and $1$ via two
separate binary hurdles and a fractional component for interior values.

``` r
# Since the k401k data has no exact 0s, we inject some artificial 0s just for demonstration
y_3p <- y
y_3p[1:50] <- 0

# Linkbin takes a vector of two links for the two binary hurdles
model_3p <- fracreg(
  y = y_3p, 
  x = X, 
  type = "3P", 
  linkbin = c("logit", "logit"), 
  linkfrac = "logit"
)
```

<details>

<summary>

Toggle to see the output
</summary>

    #> 
    #> -------------------------------------------------------------------------------- 
    #>        First binary component of a three-part model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Number of observations:                                                     1534 
    #> Pseudo R-squared:                                                        0.00324 
    #> Log-Likelihood:                                                        -216.6222 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>           Estimate Std. Error z value Pr(>|z|)    
    #> Constant 2.867e+00  3.157e-01   9.080   <2e-16 ***
    #> mrate    1.147e-01  2.131e-01   0.538    0.590    
    #> age      8.375e-03  1.765e-02   0.474    0.635    
    #> totemp   1.036e-04  6.959e-05   1.489    0.136    
    #> sole     3.371e-01  2.983e-01   1.130    0.259    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>       Second binary component of a three-part model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Number of observations:                                                     1484 
    #> Pseudo R-squared:                                                        0.15199 
    #> Log-Likelihood:                                                        -903.5457 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>            Estimate Std. Error z value Pr(>|z|)    
    #> Constant -1.435e+00  1.306e-01 -10.990  < 2e-16 ***
    #> mrate     9.244e-01  9.877e-02   9.360  < 2e-16 ***
    #> age       1.198e-02  6.339e-03   1.890   0.0588 .  
    #> totemp   -1.371e-05  6.317e-06  -2.170   0.0300 *  
    #> sole      8.852e-01  1.154e-01   7.674 1.67e-14 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>         Fractional component of a three-part model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Number of observations:                                                      826 
    #> Pseudo R-squared:                                                        0.09937 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>            Estimate Std. Error z value Pr(>|z|)    
    #> Constant  7.388e-01  7.039e-02  10.496  < 2e-16 ***
    #> mrate     3.960e-01  1.019e-01   3.885 0.000102 ***
    #> age       2.531e-02  4.068e-03   6.223 4.88e-10 ***
    #> totemp   -3.817e-06  3.088e-06  -1.236 0.216472    
    #> sole     -4.483e-03  6.672e-02  -0.067 0.946434    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>        Three-part model - binary logit , binary logit + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #> Pseudo R-squared:                                                        0.07934 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 23:18:44 
    #> --------------------------------------------------------------------------------

</details>

``` r
pe_3p <- fracreg.pe(model_3p)
```

<details>

<summary>

Toggle to see the output
</summary>

    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                             Average partial effects 
    #> -------------------------------------------------------------------------------- 
    #>        Three-part model - binary logit , binary logit + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #>         Estimate Std. Error z value Pr(>|z|)    
    #> mrate  8.110e-02  1.162e-02   6.977 3.02e-12 ***
    #> age    3.124e-03  6.764e-04   4.619 3.86e-06 ***
    #> totemp 1.858e-06  1.971e-06   0.942    0.346    
    #> sole   4.800e-02  1.116e-02   4.302 1.69e-05 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 23:18:44 
    #> --------------------------------------------------------------------------------

</details>

------------------------------------------------------------------------

### 2. Hypothesis Testing and Specification Diagnostics

#### 2.1 Generalized Goodness-Of-Functional-Form (GGOFF)

The GGOFF tests if the specified link function $G(\cdot)$ is appropriate
by testing $H_0: \gamma_1 = 0, \gamma_2 = 0$ in the augmented model
$E(y|x) = G\left(x\beta + \gamma_1 \frac{g'(x\hat{\beta})}{g(x\hat{\beta})} + \gamma_2 x\hat{\beta} \right)$.

``` r
# Perform the GGOFF test on the 1P model
fracreg.ggoff(model_1p, version = "LM")
```

<details>

<summary>

Toggle to see the output
</summary>

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    GGOFF test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit model 
    #> -------------------------------------------------------------------------------- 
    #>            Statistic p-value   
    #> GOFF1 - LM     8.838 0.00295 **
    #> GOFF2 - LM     9.828 0.00172 **
    #> GGOFF - LM    10.351 0.00565 **
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 23:18:44 
    #> --------------------------------------------------------------------------------

</details>

#### 2.2 RESET Test

The RESET test checks for general functional form misspecification by
augmenting the model with powers of the fitted index $(x\hat{\beta})^k$.

``` r
# Standard RESET test for functional form misspecification
fracreg.reset(model_1p)
```

<details>

<summary>

Toggle to see the output
</summary>

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    RESET test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit model 
    #> -------------------------------------------------------------------------------- 
    #>       Statistic p-value   
    #> LM(3)     10.29 0.00583 **
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 23:18:44 
    #> --------------------------------------------------------------------------------

</details>

------------------------------------------------------------------------

### 3. Endogeneity & Heteroscedasticity (`fracreghet`)

Use instrumental variables (`z`) to correct for endogenous covariates
via a control function approach
$E(y_{1i} \mid z_i, y_{2i}, v_i) = G(x_i \beta + \gamma \hat{v}_i)$, or
via IV-GMM.

``` r
# Simulating an endogenous variable (var.endog) and an instrument (z)
set.seed(42)
N <- length(y)
x1 <- rnorm(N)
z <- rnorm(N)
u <- 0.5 * z + rnorm(N)
var.endog <- 0.8 * z + u
y_endog <- exp(0.5 * x1 + 1.2 * var.endog + u) / (1 + exp(0.5 * x1 + 1.2 * var.endog + u))

# Estimate QML with control function for cross-sectional data
model_het <- fracreghet(
  y = y_endog, 
  x = cbind(x1, var.endog), 
  z = cbind(x1, z), 
  var.endog = var.endog,
  type = "QMLxv", 
  link = "logit"
)
```

<details>

<summary>

Toggle to see the output
</summary>

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                        Fractional logit regression model 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                 QMLxv 
    #> Number of observations:                                                     1534 
    #> Standard errors:                                                          robust 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>           Estimate Std. Error z value Pr(>|z|)    
    #> Constant  -0.01901    0.01643  -1.157    0.247    
    #> x1         0.49108    0.01622  30.273   <2e-16 ***
    #> var.endog  1.56341    0.01250 125.071   <2e-16 ***
    #> vhat       0.63659    0.01250  50.927   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #>                                  Reduced form: 
    #> -------------------------------------------------------------------------------- 
    #>             Estimate Std. Error z value Pr(>|z|)    
    #> Z_INTERCEPT -0.02986    0.02583  -1.156    0.248    
    #> Z_x1        -0.01401    0.02549  -0.550    0.583    
    #> Z_z          1.25670    0.02468  50.927   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 23:18:45 
    #> --------------------------------------------------------------------------------

</details>

------------------------------------------------------------------------

### 4. Panel Data Fractional Models (`fracregpd`)

For longitudinal or panel data where unobserved heterogeneity is a
concern, `fracregpd()` provides fixed-T panel estimators. It supports
Correlated Random Effects (CRE), where unobserved heterogeneity is
projected onto the time-averages of the covariates:
$c_i = \psi + \bar{x}_i \xi + a_i$.

``` r
# Simulating Panel Data
id <- rep(1:200, each = 5)
time <- rep(1:5, times = 200)
x_panel <- rnorm(1000)
# Unobserved individual effect (CRE)
c_i <- rep(rnorm(200), each = 5) 
y_panel <- exp(x_panel + c_i) / (1 + exp(x_panel + c_i))

# Estimate a Correlated Random Effects (CRE) Model
model_pd <- fracregpd(
  id = id, 
  time = time, 
  y = y_panel, 
  x = cbind(x_panel = x_panel), 
  type = "QMLcre", 
  link = "probit"
)
```

<details>

<summary>

Toggle to see the output
</summary>

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                        Fractional probit regression model 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                QMLcre 
    #> Exogeneity:                                                                 TRUE 
    #> Use first lag of instruments:                                              FALSE 
    #> Standard errors:                                                         cluster 
    #> -------------------------------------------------------------------------------- 
    #> Initial observations:                                                       1000 
    #> Estimation observations:                                                    1000 
    #> Initial cross-sectional units:                                               200 
    #> Estimation cross-sectional units:                                            200 
    #> Initial periods per unit (avg):                                                5 
    #> Estimation periods per unit (avg):                                             5 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>                Estimate Std. Error z value Pr(>|z|)    
    #> x_panel         0.50721    0.00864  58.705   <2e-16 ***
    #> INTERCEPT_mean -0.02762    0.03777  -0.731    0.465    
    #> x_panel_mean   -0.00904    0.08772  -0.103    0.918    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 23:18:45 
    #> --------------------------------------------------------------------------------

</details>

------------------------------------------------------------------------

## References

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
