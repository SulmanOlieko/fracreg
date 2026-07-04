# A Comprehensive Guide to Fractional Regression Models with fracreg

## Introduction

The **`fracreg`** package is the most comprehensive tool available in R
for estimating, diagnosing, and interpreting fractional regression
models. It supports:

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

This vignette walks you through simulated examples demonstrating each
estimator.

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

------------------------------------------------------------------------

## 1. Cross-Sectional Fractional Models (`fracreg`)

The core
[`fracreg()`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.md)
function is designed for univariate models where the dependent variable
is bounded between 0 and 1 inclusive ($`0 \le y \le 1`$).

### 1.1 Simulating the Data

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

    #>           y          x1        x2
    #> 1 0.0000000 -0.56047565 0.1596740
    #> 2 0.3061926 -0.23017749 0.1445159
    #> 3 0.5059710  1.55870831 0.1491804
    #> 4 0.6901448  0.07050839 0.5144343
    #> 5 0.7883783  0.12928774 0.4928273
    #> 6 1.0000000  1.71506499 0.6163428

### 1.2 The Standard 1-Part Model (1P)

If you assume a single process governs the entire distribution of `y`,
the standard one-part model is the natural choice. This corresponds to
the Papke and Wooldridge (1996) quasi-likelihood estimator.

``` r
model_1p <- fracreg(
  y = data$y, 
  x = cbind(x1 = data$x1, x2 = data$x2), 
  type = "1P", 
  linkfrac = "logit"
)
```

Toggle to see the output

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

``` r
summary(model_1p)
```

Toggle to see the output

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

### 1.3 Partial Effects (1P)

Raw coefficients from fractional models are not directly interpretable
as marginal effects. Use
[`fracreg.pe()`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md)
to compute the **Average Partial Effects (APE)** using the analytical
delta method:

``` r
pe_1p <- fracreg.pe(model_1p)
```

Toggle to see the output

    #> 
    #> 
    #> *** Average partial effects ***
    #> 
    #> Fractional logit model
    #> 
    #>    Estimate Std. Error t value Pr(>|t|)    
    #> x1   0.1682     0.0068  24.694    0.000 ***
    #> x2   0.2894     0.0239  12.132    0.000 ***

``` r
summary(pe_1p)
```

Toggle to see the output

    #>       Length Class  Mode   
    #> PE.p  2      -none- numeric
    #> PE.sd 2      -none- numeric

### 1.4 The Two-Part Hurdle Model (2P)

Often, fractional data (e.g., participation rates, insurance coverage,
expenditure shares) have a substantial mass of zeros. The two-part model
estimates:

1.  A **binary** component: $`\Pr(y > 0 \mid x)`$
2.  A **fractional** component: $`E(y \mid y > 0, x)`$

The overall conditional expectation is then
$`E(y \mid x) = \Pr(y > 0) \cdot E(y \mid y > 0)`$.

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

``` r
summary(model_2p)
```

Toggle to see the output

    #>           Length Class  Mode     
    #> resBIN      15   -none- list     
    #> resFRAC     14   -none- list     
    #> class        1   -none- character
    #> type         1   -none- character
    #> ybase     1000   -none- numeric  
    #> x2base    3000   -none- numeric  
    #> yhat2P    1000   -none- numeric  
    #> converged    1   -none- numeric

#### Partial Effects for the Two-Part Model

The
[`fracreg.pe()`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md)
command natively accounts for the product rule when computing APEs for
combined two-part models:

``` r
pe_2p <- fracreg.pe(model_2p)
```

Toggle to see the output

    #> 
    #> 
    #> *** Average partial effects ***
    #> 
    #> Binary logit + Fractional logit two-part model
    #> 
    #>    Estimate Std. Error t value Pr(>|t|)    
    #> x1   0.1666     0.0067  24.804    0.000 ***
    #> x2   0.2848     0.0245  11.612    0.000 ***

``` r
summary(pe_2p)
```

Toggle to see the output

    #>       Length Class  Mode   
    #> PE.p  2      -none- numeric
    #> PE.sd 2      -none- numeric

### 1.5 The Three-Part Double-Inflated Model (3P)

When data has substantial clustering at **both** boundaries (0 and 1)
simultaneously, the three-part model decomposes the process into:

1.  $`\Pr(y > 0 \mid x)`$ — first binary hurdle
2.  $`\Pr(y = 1 \mid y > 0, x)`$ — second binary hurdle
3.  $`E(y \mid 0 < y < 1, x)`$ — fractional component

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

``` r
summary(model_3p)
```

Toggle to see the output

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

#### Partial Effects for the Three-Part Model

The analytical delta method is natively supported for the 3P model,
correctly accounting for all three components:

``` r
pe_3p <- fracreg.pe(model_3p)
```

Toggle to see the output

    #> 
    #> 
    #> *** Average partial effects ***
    #> 
    #> Three-part model - binary logit , binary logit + fractional logit
    #> 
    #>    Estimate Std. Error t value Pr(>|t|)    
    #> x1   0.1663     0.0066  25.127    0.000 ***
    #> x2   0.2789     0.0249  11.186    0.000 ***

``` r
summary(pe_3p)
```

Toggle to see the output

    #>       Length Class  Mode   
    #> PE.p  2      -none- numeric
    #> PE.sd 2      -none- numeric

------------------------------------------------------------------------

## 2. Hypothesis Testing and Specification Diagnostics

`fracreg` includes state-of-the-art specification tests to validate your
model’s functional form and link function assumptions.

### 2.1 Generalized Goodness-Of-Functional-Form (GGOFF)

The GGOFF test (Ramalho and Ramalho, 2012) tests whether the chosen link
function is adequate for the data. A significant result suggests the
link may be misspecified.

``` r
# Perform the GGOFF test on the 1P model
# Tests if the logit link is appropriate
fracreg.ggoff(model_1p, version = "LM")
```

Toggle to see the output

    #> 
    #> *** GGOFF test ***
    #> 
    #> H0:  Fractional logit model
    #> 
    #>   Test Version Statistic p-value 
    #>  GOFF1      LM     0.640   0.424 
    #>  GOFF2      LM     0.485   0.486 
    #>  GGOFF      LM     0.814   0.666

### 2.2 RESET Test

The RESET test detects general functional form misspecification by
testing whether powers of the fitted values have explanatory power:

``` r
# Standard RESET test for functional form misspecification
fracreg.reset(model_1p)
```

Toggle to see the output

    #> 
    #> *** RESET test ***
    #> 
    #> H0:  Fractional logit model
    #> 
    #>  Version Statistic p-value 
    #>    LM(3)     0.852   0.653

### 2.3 P-Test for Non-Nested Models

You can compare non-nested models (e.g., `logit` vs. `cloglog` link)
using the Davidson-MacKinnon P-test:

``` r
model_1p_clog <- fracreg(
  y = data$y, 
  x = cbind(x1 = data$x1, x2 = data$x2), 
  type = "1P", 
  linkfrac = "cloglog"
)
```

Toggle to see the output

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

``` r
summary(model_1p_clog)
```

Toggle to see the output

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

``` r

# fracreg.ptest(model_1p, model_1p_clog) # Uncomment to run non-nested test
```

------------------------------------------------------------------------

## 3. Endogeneity & Heteroskedasticity (`fracreghet`)

When you suspect that one of your covariates is endogenous, or that the
variance is heteroskedastic,
[`fracreghet()`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.md)
provides instrumental variable correction. It natively supports IV
through a **Control Function (CF)** approach or **GMM** estimation.

### 3.1 Control Function Approach (QMLxv)

The control function approach adds the first-stage residual as an
additional regressor in the fractional model, providing a Hausman-type
test for endogeneity.

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

``` r
summary(model_het)
```

Toggle to see the output

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

------------------------------------------------------------------------

## 4. Panel Data Fractional Models (`fracregpd`)

For longitudinal or panel data where unobserved heterogeneity is a
concern,
[`fracregpd()`](https://SulmanOlieko.github.io/fracreg/reference/fracregpd.md)
provides fixed-T panel estimators. The Correlated Random Effects (CRE)
approach adds individual-specific time averages of covariates to the
model, allowing consistent estimation even when unobserved effects are
correlated with regressors.

### 4.1 Correlated Random Effects (QMLcre)

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

``` r
summary(model_pd)
```

Toggle to see the output

    #>           Length Class  Mode     
    #> type         1   -none- character
    #> link         1   -none- character
    #> Hy        1000   -none- numeric  
    #> p            3   -none- numeric  
    #> converged    1   -none- numeric  
    #> p.var        9   -none- numeric  
    #> var.type     1   -none- character

------------------------------------------------------------------------

## Conclusion

With `fracreg`, `fracreghet`, and `fracregpd`, you have a complete,
mathematically robust toolkit for modeling fractional responses bounded
between $`[0,1]`$, regardless of inflation, endogeneity, or unobserved
panel effects.

For more information, please visit the [package
website](https://SulmanOlieko.github.io/fracreg) or file an issue on
[GitHub](https://github.com/SulmanOlieko/fracreg/issues).

To cite this package in your research:

``` r
citation("fracreg")
```
