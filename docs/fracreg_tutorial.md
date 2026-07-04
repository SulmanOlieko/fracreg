# NA

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
#> Error in `library()`:
#> ! there is no package called 'fracreg'
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
#>           y          x1        x2
#> 1 0.0000000 -0.56047565 0.1596740
#> 2 0.3061926 -0.23017749 0.1445159
#> 3 0.5059710  1.55870831 0.1491804
#> 4 0.6901448  0.07050839 0.5144343
#> 5 0.7883783  0.12928774 0.4928273
#> 6 1.0000000  1.71506499 0.6163428
```

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
#> Error in `fracreg()`:
#> ! could not find function "fracreg"
summary(model_1p)
#> Error:
#> ! object 'model_1p' not found
```

### 1.3 Partial Effects (1P)

Raw coefficients from fractional models are not directly interpretable
as marginal effects. Use
[`fracreg.pe()`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md)
to compute the **Average Partial Effects (APE)** using the analytical
delta method:

``` r
pe_1p <- fracreg.pe(model_1p)
#> Error in `fracreg.pe()`:
#> ! could not find function "fracreg.pe"
summary(pe_1p)
#> Error:
#> ! object 'pe_1p' not found
```

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
#> Error in `fracreg()`:
#> ! could not find function "fracreg"
summary(model_2p)
#> Error:
#> ! object 'model_2p' not found
```

#### Partial Effects for the Two-Part Model

The
[`fracreg.pe()`](https://SulmanOlieko.github.io/fracreg/reference/fracreg.pe.md)
command natively accounts for the product rule when computing APEs for
combined two-part models:

``` r
pe_2p <- fracreg.pe(model_2p)
#> Error in `fracreg.pe()`:
#> ! could not find function "fracreg.pe"
summary(pe_2p)
#> Error:
#> ! object 'pe_2p' not found
```

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
#> Error in `fracreg()`:
#> ! could not find function "fracreg"
summary(model_3p)
#> Error:
#> ! object 'model_3p' not found
```

#### Partial Effects for the Three-Part Model

The analytical delta method is natively supported for the 3P model,
correctly accounting for all three components:

``` r
pe_3p <- fracreg.pe(model_3p)
#> Error in `fracreg.pe()`:
#> ! could not find function "fracreg.pe"
summary(pe_3p)
#> Error:
#> ! object 'pe_3p' not found
```

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
#> Error in `fracreg.ggoff()`:
#> ! could not find function "fracreg.ggoff"
```

### 2.2 RESET Test

The RESET test detects general functional form misspecification by
testing whether powers of the fitted values have explanatory power:

``` r
# Standard RESET test for functional form misspecification
fracreg.reset(model_1p)
#> Error in `fracreg.reset()`:
#> ! could not find function "fracreg.reset"
```

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
#> Error in `fracreg()`:
#> ! could not find function "fracreg"
summary(model_1p_clog)
#> Error:
#> ! object 'model_1p_clog' not found

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
#> Error in `fracreghet()`:
#> ! could not find function "fracreghet"
summary(model_het)
#> Error:
#> ! object 'model_het' not found
```

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
#> Error in `fracregpd()`:
#> ! could not find function "fracregpd"
summary(model_pd)
#> Error:
#> ! object 'model_pd' not found
```

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
