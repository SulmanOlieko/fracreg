# fracreg: Fractional Regression Models in R

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![](https://img.shields.io/badge/R-%3E=%203.5.0-blue.svg)](https://cran.r-project.org/)
[![](https://img.shields.io/github/last-commit/SulmanOlieko/fracreg.svg)](https://github.com/SulmanOlieko/fracreg/commits/main)
[![](https://img.shields.io/github/issues/SulmanOlieko/fracreg.svg)](https://github.com/SulmanOlieko/fracreg/issues)
[![](https://img.shields.io/badge/license-GPL--3-blue)](https://github.com/SulmanOlieko/fracreg)
[![CRAN Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/fracreg)](https://cran.r-project.org/package=fracreg)

> **The Ultimate Toolkit for Fractional Data**

An advanced R package for the estimation, specification analysis, and
interpretation of fractional regression models. It handles univariate
proportions, participation rates, and bounded data ($`0 \le y \le 1`$),
providing state-of-the-art estimators for single-part, two-part
(hurdle), and three-part (double-inflated) models.

## Features

`fracreg` is fully equipped to handle complex data structures commonly
found in econometrics and biostatistics:

1.  **Univariate Fractional Models (`fracreg`)**: Fit standard 1-part
    models, hurdle 2-part models (mass at 0 or 1), and double-inflated
    3-part models.
2.  **Panel Data Models (`fracregpd`)**: Estimate fractional models with
    fixed-T longitudinal data, supporting Correlated Random Effects
    (CRE) to handle unobserved individual heterogeneity.
3.  **Endogeneity & Heteroskedasticity (`fracreghet`)**: Correct for
    endogenous covariates using Instrumental Variables (IV) via Control
    Function and GMM approaches.
4.  **Diagnostic Testing**: Includes generalized
    goodness-of-functional-form (`fracreg.ggoff`), RESET
    (`fracreg.reset`), and non-nested P-tests (`fracreg.ptest`).
5.  **Analytical Partial Effects**: Compute exact marginal effects for
    all model types (including 3-part composite effects) via the Delta
    method using
    [`fracreg.pe()`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md).

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

The following sections provide examples covering the three primary
estimation functions.

``` r
library(fracreg)
#> =========================================================
#>  fracreg: Fractional Regression Models
#> =========================================================
#>  For a comprehensive guide and tutorials, please visit:
#>  https://SulmanOlieko.github.io/fracreg
#> ---------------------------------------------------------
#>  To cite fracreg in publications, please type:
#>  citation("fracreg")
#> =========================================================
```

### 1. The Core Estimator (`fracreg`)

The main
[`fracreg()`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md)
function fits cross-sectional models.

#### 1.1 Simulating the Data

``` r
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

data <- data.frame(y = y, x1 = x1, x2 = x2)
head(data)
```

Toggle to see the output

``` R
#>           y          x1        x2
#> 1 0.0000000 -0.56047565 0.1596740
#> 2 0.3061926 -0.23017749 0.1445159
#> 3 0.5059710  1.55870831 0.1491804
#> 4 0.6901448  0.07050839 0.5144343
#> 5 0.7883783  0.12928774 0.4928273
#> 6 1.0000000  1.71506499 0.6163428
```

#### 1.2 The Standard 1-Part Model (1P)

``` r
model_1p <- fracreg(
  y = data$y, 
  x = cbind(x1 = data$x1, x2 = data$x2), 
  type = "1P", 
  linkfrac = "logit"
)
```

Toggle to see the output

``` R
#> 
#> *** Fractional logit regression model ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -0.569685   0.067502  -8.440    0.000 ***
#> x1         0.788220   0.040145  19.634    0.000 ***
#> x2         1.356110   0.118992  11.397    0.000 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 1000 
#> R-squared: 0.39
```

``` r
summary(model_1p)
```

Toggle to see the output

``` R
#>           Length Class   Mode     
#> class        1   -none-  character
#> formula      3   formula call     
#> type         1   -none-  character
#> link         1   -none-  character
#> method       1   -none-  character
#> p            3   -none-  numeric  
#> yhat      1000   -none-  numeric  
#> xbhat     1000   -none-  numeric  
#> converged    1   -none-  numeric  
#> x.names      3   -none-  character
#> p.var        9   -none-  numeric  
#> var.type     1   -none-  character
#> var.eim      1   -none-  logical  
#> dfc          1   -none-  logical
```

#### 1.3 Partial Effects (1P)

``` r
pe_1p <- fracreg.pe(model_1p)
```

Toggle to see the output

``` R
#> 
#> 
#> *** Average partial effects ***
#> 
#> Fractional logit model
#> 
#>    Estimate Std. Error t value Pr(>|t|)    
#> x1   0.1682     0.0068  24.694    0.000 ***
#> x2   0.2894     0.0239  12.132    0.000 ***
```

``` r
summary(pe_1p)
```

Toggle to see the output

``` R
#>       Length Class  Mode   
#> PE.p  2      -none- numeric
#> PE.sd 2      -none- numeric
```

#### 1.4 The Two-Part Hurdle Model (2P)

``` r
# We use inflation=0 to indicate we are modeling the mass at zero
model_2p <- fracreg(
  y = data$y, 
  x = cbind(x1 = data$x1, x2 = data$x2), 
  type = "2P", 
  inflation = 0,
  linkbin = "logit", 
  linkfrac = "logit"
)
```

Toggle to see the output

``` R
#> 
#> *** Binary component of a two-part model - logit specification ***
#> 
#>           Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT 1.463769   0.193283   7.573    0.000 ***
#> x1        1.185749   0.128710   9.213    0.000 ***
#> x2        2.227886   0.393191   5.666    0.000 ***
#> 
#> Number of observations: 1000 
#> R-squared: 0.145 
#> 
#> 
#> 
#> *** Fractional component of a two-part model - logit specification ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -0.173232   0.057863  -2.994    0.003 ***
#> x1         0.612050   0.035541  17.221    0.000 ***
#> x2         1.005089   0.107030   9.391    0.000 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 881 
#> R-squared: 0.324 
#> 
#> 
#> 
#> *** Two-part model - binary logit + fractional logit  ***
#> 
#> R-squared: 0.386
```

``` r
summary(model_2p)
```

Toggle to see the output

``` R
#>           Length Class  Mode     
#> resBIN      15   -none- list     
#> resFRAC     14   -none- list     
#> class        1   -none- character
#> type         1   -none- character
#> ybase     1000   -none- numeric  
#> x2base    3000   -none- numeric  
#> yhat2P    1000   -none- numeric  
#> converged    1   -none- numeric
```

``` r
pe_2p <- fracreg.pe(model_2p)
```

Toggle to see the output

``` R
#> 
#> 
#> *** Average partial effects ***
#> 
#> Binary logit + Fractional logit two-part model
#> 
#>    Estimate Std. Error t value Pr(>|t|)    
#> x1   0.1666     0.0067  24.804    0.000 ***
#> x2   0.2848     0.0245  11.612    0.000 ***
```

``` r
summary(pe_2p)
```

Toggle to see the output

``` R
#>       Length Class  Mode   
#> PE.p  2      -none- numeric
#> PE.sd 2      -none- numeric
```

#### 1.5 The Three-Part Double-Inflated Model (3P)

``` r
# Linkbin takes a vector of two links for the two binary hurdles
model_3p <- fracreg(
  y = data$y, 
  x = cbind(x1 = data$x1, x2 = data$x2), 
  type = "3P", 
  linkbin = c("logit", "logit"), 
  linkfrac = "logit"
)
```

Toggle to see the output

``` R
#> 
#> *** First binary component of a three-part model - logit specification ***
#> 
#>           Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT 1.463769   0.193283   7.573    0.000 ***
#> x1        1.185749   0.128710   9.213    0.000 ***
#> x2        2.227886   0.393191   5.666    0.000 ***
#> 
#> Number of observations: 1000 
#> R-squared: 0.145 
#> 
#> 
#> 
#> *** Second binary component of a three-part model - logit specification ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -2.953982   0.243671 -12.123    0.000 ***
#> x1         1.168175   0.116745  10.006    0.000 ***
#> x2         1.831493   0.345480   5.301    0.000 ***
#> 
#> Number of observations: 881 
#> R-squared: 0.188 
#> 
#> 
#> 
#> *** Fractional component of a three-part model - logit specification ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -0.267629   0.045393  -5.896    0.000 ***
#> x1         0.361985   0.024539  14.751    0.000 ***
#> x2         0.585047   0.079687   7.342    0.000 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 715 
#> R-squared: 0.243 
#> 
#> 
#> 
#> *** Three-part model - binary logit , binary logit + fractional logit  ***
#> 
#> R-squared: 0.389
```

``` r
summary(model_3p)
```

Toggle to see the output

``` R
#>           Length Class  Mode     
#> resBIN0     13   -none- list     
#> resBIN1     13   -none- list     
#> resFRAC     12   -none- list     
#> class        1   -none- character
#> type         1   -none- character
#> ybase     1000   -none- numeric  
#> x2base    3000   -none- numeric  
#> yhat3P    1000   -none- numeric  
#> converged    1   -none- numeric
```

``` r
pe_3p <- fracreg.pe(model_3p)
```

Toggle to see the output

``` R
#> 
#> 
#> *** Average partial effects ***
#> 
#> Three-part model - binary logit , binary logit + fractional logit
#> 
#>    Estimate Std. Error t value Pr(>|t|)    
#> x1   0.1663     0.0066  25.127    0.000 ***
#> x2   0.2789     0.0249  11.186    0.000 ***
```

``` r
summary(pe_3p)
```

Toggle to see the output

``` R
#>       Length Class  Mode   
#> PE.p  2      -none- numeric
#> PE.sd 2      -none- numeric
```

------------------------------------------------------------------------

### 2. Hypothesis Testing and Specification Diagnostics

#### 2.1 Generalized Goodness-Of-Functional-Form (GGOFF)

``` r
# Perform the GGOFF test on the 1P model
fracreg.ggoff(model_1p, version = "LM")
```

Toggle to see the output

``` R
#> 
#> *** GGOFF test ***
#> 
#> H0:  Fractional logit model
#> 
#>   Test Version Statistic p-value 
#>  GOFF1      LM     0.640   0.424 
#>  GOFF2      LM     0.485   0.486 
#>  GGOFF      LM     0.814   0.666
```

#### 2.2 RESET Test

``` r
# Standard RESET test for functional form misspecification
fracreg.reset(model_1p)
```

Toggle to see the output

``` R
#> 
#> *** RESET test ***
#> 
#> H0:  Fractional logit model
#> 
#>  Version Statistic p-value 
#>    LM(3)     0.852   0.653
```

------------------------------------------------------------------------

### 3. Endogeneity & Heteroskedasticity (`fracreghet`)

Use instrumental variables (`z`) to correct for endogenous covariates
via a control function.

``` r
# Simulating an endogenous variable (var.endog) and an instrument (z)
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

Toggle to see the output

``` R
#> 
#> *** Fractional logit regression model ***
#> *** Estimator: QMLxv
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -0.012887   0.018543  -0.695    0.487    
#> x1         0.521281   0.018894  27.589    0.000 ***
#> var.endog  1.611308   0.012726 126.614    0.000 ***
#> vhat       0.588692   0.012726  46.259    0.000 ***
#> 
#> Reduced form:
#>              Estimate Std. Error t value Pr(>|t|)    
#> Z_INTERCEPT -0.021891   0.031491  -0.695    0.487    
#> Z_x1         0.036150   0.032150   1.124    0.261    
#> Z_z          1.358944   0.029377  46.259    0.000 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 1000
```

``` r
summary(model_het)
```

Toggle to see the output

``` R
#>           Length Class   Mode     
#> class        1   -none-  character
#> formula      3   formula call     
#> type         1   -none-  character
#> link         1   -none-  character
#> adjust       1   -none-  numeric  
#> p            7   -none-  numeric  
#> Hy        1000   -none-  numeric  
#> xbhat     1000   -none-  numeric  
#> converged    1   -none-  logical  
#> x.names      3   -none-  character
#> p.var       49   -none-  numeric  
#> var.type     1   -none-  character
```

------------------------------------------------------------------------

### 4. Panel Data Fractional Models (`fracregpd`)

For longitudinal or panel data where unobserved heterogeneity is a
concern,
[`fracregpd()`](https://SulmanOlieko.github.io/fracreg/reference/fracregpd.md)
provides fixed-T panel estimators.

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

Toggle to see the output

``` R
#> 
#> *** Fractional probit regression model ***
#> *** Estimator: QMLcre
#> *** Exogeneity: TRUE
#> *** Use first lag of instruments: FALSE
#> 
#>                 Estimate Std. Error t value Pr(>|t|)    
#> x_panel         0.512705   0.008692  58.985    0.000 ***
#> INTERCEPT_mean  0.013508   0.037058   0.365    0.715    
#> x_panel_mean   -0.122045   0.081551  -1.497    0.135    
#> 
#> Note: cluster standard errors
#> 
#> Number of observations (initial): 1000 
#> Number of observations (for estimation): 1000 
#> Number of cross-sectional units (initial): 200 
#> Number of cross-sectional units (for estimation): 200 
#> Average number of time periods per cross-sectional unit (initial): 5 
#> Average number of time periods per cross-sectional unit (for estimation): 5
```

``` r
summary(model_pd)
```

Toggle to see the output

``` R
#>           Length Class  Mode     
#> type         1   -none- character
#> link         1   -none- character
#> Hy        1000   -none- numeric  
#> p            3   -none- numeric  
#> converged    1   -none- numeric  
#> p.var        9   -none- numeric  
#> var.type     1   -none- character
```

------------------------------------------------------------------------

## References

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), “Alternative
estimating and testing empirical strategies for fractional regression
models”, *Journal of Economic Surveys*, 25(1), 19-68.

Papke, L.E. and J.M. Wooldridge (1996), “Econometric methods for
fractional response variables with an application to 401(k) plan
participation rates”, *Journal of Applied Econometrics*, 11(6), 619-632.
