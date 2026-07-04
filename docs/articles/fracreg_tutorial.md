# A Comprehensive Guide to Fractional Regression Models with fracreg

## Introduction

The **`fracreg`** package is the most comprehensive tool available in R
for estimating, diagnosing, and interpreting fractional regression
models. It supports: - Univariate fractional models (1P) - Hurdle
two-part fractional models (2P) - Double-inflated three-part fractional
models (3P) - Fractional models with panel data (`fracregpd`) -
Fractional models with endogenous covariates or heteroskedasticity
(`fracreghet`)

This vignette will walk you through simulated examples to demonstrate
how to use each of these powerful estimators.

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

## 1. Cross-Sectional Fractional Models (`fracreg`)

The core
[`fracreg()`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md)
function is designed to handle univariate models where the dependent
variable is bounded between 0 and 1 inclusive ($`0 \le y \le 1`$). We
will simulate some data to demonstrate its usage.

### Simulating the Data

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
#>           y          x1        x2
#> 1 0.0000000 -0.56047565 0.1596740
#> 2 0.3061926 -0.23017749 0.1445159
#> 3 0.5059710  1.55870831 0.1491804
#> 4 0.6901448  0.07050839 0.5144343
#> 5 0.7883783  0.12928774 0.4928273
#> 6 1.0000000  1.71506499 0.6163428
```

### The Standard 1-Part Model (1P)

If you assume a single process governs the entire distribution of `y`,
you can use the standard one-part model.

``` r
model_1p <- fracreg(
  y = data$y, 
  x = cbind(x1 = data$x1, x2 = data$x2), 
  type = "1P", 
  linkfrac = "logit"
)
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
summary(model_1p)
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

### Partial Effects (1P)

To interpret the coefficients in terms of their marginal impact on
$`y`$, you can compute the partial effects using the
[`fracreg.pe()`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md)
command:

``` r
pe_1p <- fracreg.pe(model_1p)
#> 
#> 
#> *** Average partial effects ***
#> 
#> Fractional logit model
#> 
#>    Estimate Std. Error t value Pr(>|t|)    
#> x1   0.1682     0.0068  24.694    0.000 ***
#> x2   0.2894     0.0239  12.132    0.000 ***
summary(pe_1p)
#>       Length Class  Mode   
#> PE.p  2      -none- numeric
#> PE.sd 2      -none- numeric
```

### The Two-Part Hurdle Model (2P)

Often, fractional data (like participation rates or test scores) have a
large cluster of zeros. The two-part model estimates the binary “hurdle”
($`y>0`$) and the continuous fractional part ($`0 < y \le 1`$)
sequentially, then combines them.

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
summary(model_2p)
#>           Length Class  Mode     
#> resBIN      15   -none- list     
#> resFRAC     14   -none- list     
#> class        1   -none- character
#> type         1   -none- character
#> ybase     1000   -none- numeric  
#> x2base    3000   -none- numeric  
#> yhat2P    1000   -none- numeric  
#> converged    1   -none- numeric

# Calculate partial effects for the combined two-part model
pe_2p <- fracreg.pe(model_2p)
#> 
#> 
#> *** Average partial effects ***
#> 
#> Binary logit + Fractional logit two-part model
#> 
#>    Estimate Std. Error t value Pr(>|t|)    
#> x1   0.1666     0.0067  24.804    0.000 ***
#> x2   0.2848     0.0245  11.612    0.000 ***
summary(pe_2p)
#>       Length Class  Mode   
#> PE.p  2      -none- numeric
#> PE.sd 2      -none- numeric
```

### The Three-Part Double-Inflated Model (3P)

When your data has substantial clustering at **both** boundaries (0
and 1) simultaneously, the 3-part model decomposes the process into
three stages: 1. The probability of $`y > 0`$ 2. The probability of
$`y = 1`$ given $`y > 0`$ 3. The fractional outcome $`0 < y < 1`$

``` r
# Linkbin takes a vector of two links for the two binary hurdles
model_3p <- fracreg(
  y = data$y, 
  x = cbind(x1 = data$x1, x2 = data$x2), 
  type = "3P", 
  linkbin = c("logit", "logit"), 
  linkfrac = "logit"
)
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
summary(model_3p)
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

# The analytical delta method is natively supported for the 3P model!
pe_3p <- fracreg.pe(model_3p)
#> 
#> 
#> *** Average partial effects ***
#> 
#> Three-part model - binary logit , binary logit + fractional logit
#> 
#>    Estimate Std. Error t value Pr(>|t|)    
#> x1   0.1663     0.0066  25.127    0.000 ***
#> x2   0.2789     0.0249  11.186    0.000 ***
summary(pe_3p)
#>       Length Class  Mode   
#> PE.p  2      -none- numeric
#> PE.sd 2      -none- numeric
```

------------------------------------------------------------------------

## 2. Hypothesis Testing and Specification Diagnostics

`fracreg` includes state-of-the-art specification tests to check whether
your functional form and link functions are adequate.

### Generalized Goodness-Of-Functional-Form (GGOFF)

``` r
# Perform the GGOFF test on the 1P model
# Tests if the logit link is appropriate
fracreg.ggoff(model_1p, version = "LM")
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

### RESET Test

``` r
# Standard RESET test for functional form misspecification
fracreg.reset(model_1p)
#> 
#> *** RESET test ***
#> 
#> H0:  Fractional logit model
#> 
#>  Version Statistic p-value 
#>    LM(3)     0.852   0.653
```

### P-Test for Non-Nested Models

You can test non-nested models (e.g., comparing a `cloglog` link against
a `logit` link) using
[`fracreg.ptest()`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.ptest.md).

``` r
model_1p_clog <- fracreg(y = data$y, x = cbind(x1 = data$x1, x2 = data$x2), type = "1P", linkfrac = "cloglog")
#> 
#> *** Fractional cloglog regression model ***
#> 
#>            Estimate Std. Error t value Pr(>|t|)    
#> INTERCEPT -0.778070   0.048501 -16.043    0.000 ***
#> x1         0.522660   0.028127  18.582    0.000 ***
#> x2         0.885754   0.079842  11.094    0.000 ***
#> 
#> Note: robust standard errors
#> 
#> Number of observations: 1000 
#> R-squared: 0.388
# fracreg.ptest(model_1p, model_1p_clog) # Uncomment to run non-nested test
```

------------------------------------------------------------------------

## 3. Endogeneity & Heteroskedasticity (`fracreghet`)

When you suspect that one of your covariates is endogenous, or that the
variance is heteroskedastic, you can use `fracreghet`. It natively
supports instrumental variables (IV) through a Control Function approach
or GMM estimation.

### Control Function Approach

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

------------------------------------------------------------------------

## 4. Panel Data Fractional Models (`fracregpd`)

For longitudinal or panel data where unobserved heterogeneity is a
concern,
[`fracregpd()`](https://SulmanOlieko.github.io/fracreg/reference/fracregpd.md)
provides fixed-T panel estimators, including Correlated Random Effects
(CRE).

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

## Conclusion

With `fracreg`, `fracreghet`, and `fracregpd`, you have a complete,
mathematically robust toolkit for modeling fractional responses bounded
between \[0,1\], regardless of inflation, endogeneity, or unobserved
panel effects.
