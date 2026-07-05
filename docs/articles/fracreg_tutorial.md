# A Comprehensive Guide to Fractional Response Models with fracreg

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

This vignette walks you through simulated examples demonstrating each
estimator.

``` r
library(fracreg)
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
the Papke and Wooldridge (1996) quasi-likelihood estimator, where the
conditional mean is specified as:
``` math
 E(y \mid x) = G(x\beta) 
```
where $`G(\cdot)`$ is a known non-linear link function (e.g., logit,
probit) satisfying $`0 \le G(\cdot) \le 1`$.

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
    #> -------------------------------------------------------------------------------- 
    #>                 Fractional response model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Number of observations:                                                     1000 
    #> Pseudo R-squared:                                                         0.3903 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Estimate Std. Error z value Pr(>|z|)    
    #> Constant -0.56969    0.06750   -8.44   <2e-16 ***
    #> x1        0.78822    0.04015   19.63   <2e-16 ***
    #> x2        1.35611    0.11899   11.40   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:49 
    #> --------------------------------------------------------------------------------

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
    #> -------------------------------------------------------------------------------- 
    #>                             Average partial effects 
    #> -------------------------------------------------------------------------------- 
    #>                              Fractional logit model 
    #> -------------------------------------------------------------------------------- 
    #>    Estimate Std. Error z value Pr(>|z|)    
    #> x1 0.168231   0.006813   24.69   <2e-16 ***
    #> x2 0.289436   0.023856   12.13   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:49 
    #> --------------------------------------------------------------------------------

### 1.4 The Two-Part Hurdle Model (2P)

Often, fractional data (e.g., participation rates, insurance coverage,
expenditure shares) have a substantial mass of zeros. The two-part model
estimates:

1.  A **binary** component:
    $`\Pr(y > 0 \mid x) = G_{BIN}(x\beta_{BIN})`$
2.  A **fractional** component:
    $`E(y \mid y > 0, x) = G_{FRAC}(x\beta_{FRAC})`$

The overall conditional expectation is then
$`E(y \mid x) = \Pr(y > 0 \mid x) \cdot E(y \mid y > 0, x)`$.

``` r
# We use inflation=0 to indicate we are modelling the mass at zero
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
    #> -------------------------------------------------------------------------------- 
    #>            Binary component of a two-part model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Number of observations:                                                     1000 
    #> Pseudo R-squared:                                                        0.14527 
    #> Log-Likelihood:                                                        -295.2068 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Estimate Std. Error z value Pr(>|z|)    
    #> Constant   1.4638     0.1933   7.573 3.64e-14 ***
    #> x1         1.1857     0.1287   9.213  < 2e-16 ***
    #> x2         2.2279     0.3932   5.666 1.46e-08 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>          Fractional component of a two-part model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Number of observations:                                                      881 
    #> Pseudo R-squared:                                                        0.32412 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Estimate Std. Error z value Pr(>|z|)    
    #> Constant -0.17323    0.05786  -2.994  0.00275 ** 
    #> x1        0.61205    0.03554  17.221  < 2e-16 ***
    #> x2        1.00509    0.10703   9.391  < 2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                 Two-part model - binary logit + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #> Pseudo R-squared:                                                        0.38621 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:49 
    #> --------------------------------------------------------------------------------

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
    #> -------------------------------------------------------------------------------- 
    #>                             Average partial effects 
    #> -------------------------------------------------------------------------------- 
    #>                  Binary logit + Fractional logit two-part model 
    #> -------------------------------------------------------------------------------- 
    #>    Estimate Std. Error z value Pr(>|z|)    
    #> x1 0.166600   0.006717   24.80   <2e-16 ***
    #> x2 0.284837   0.024530   11.61   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:49 
    #> --------------------------------------------------------------------------------

### 1.5 The Three-Part Double-Inflated Model (3P)

When data has substantial clustering at **both** boundaries (0 and 1)
simultaneously, the three-part model decomposes the process into:

1.  $`\Pr(y > 0 \mid x) = G_{BIN0}(x\beta_{BIN0})`$ — first binary
    hurdle
2.  $`\Pr(y = 1 \mid y > 0, x) = G_{BIN1}(x\beta_{BIN1})`$ — second
    binary hurdle
3.  $`E(y \mid 0 < y < 1, x) = G_{FRAC}(x\beta_{FRAC})`$ — fractional
    component

The overall conditional expectation integrates these three mechanisms
into a single framework.

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
    #> -------------------------------------------------------------------------------- 
    #>        First binary component of a three-part model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Number of observations:                                                     1000 
    #> Pseudo R-squared:                                                        0.14527 
    #> Log-Likelihood:                                                        -295.2068 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Estimate Std. Error z value Pr(>|z|)    
    #> Constant   1.4638     0.1933   7.573 3.64e-14 ***
    #> x1         1.1857     0.1287   9.213  < 2e-16 ***
    #> x2         2.2279     0.3932   5.666 1.46e-08 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>       Second binary component of a three-part model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                    ML 
    #> Number of observations:                                                      881 
    #> Pseudo R-squared:                                                        0.18752 
    #> Log-Likelihood:                                                        -348.7582 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Estimate Std. Error z value Pr(>|z|)    
    #> Constant  -2.9540     0.2437 -12.123  < 2e-16 ***
    #> x1         1.1682     0.1167  10.006  < 2e-16 ***
    #> x2         1.8315     0.3455   5.301 1.15e-07 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>         Fractional component of a three-part model - logit specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Number of observations:                                                      715 
    #> Pseudo R-squared:                                                        0.24292 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Estimate Std. Error z value Pr(>|z|)    
    #> Constant -0.26763    0.04539  -5.896 3.73e-09 ***
    #> x1        0.36198    0.02454  14.751  < 2e-16 ***
    #> x2        0.58505    0.07969   7.342 2.11e-13 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>        Three-part model - binary logit , binary logit + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #> Pseudo R-squared:                                                        0.38893 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:49 
    #> --------------------------------------------------------------------------------

#### Partial Effects for the Three-Part Model

The analytical delta method is natively supported for the 3P model,
correctly accounting for all three components:

``` r
pe_3p <- fracreg.pe(model_3p)
```

Toggle to see the output

    #> 
    #> 
    #> -------------------------------------------------------------------------------- 
    #>                             Average partial effects 
    #> -------------------------------------------------------------------------------- 
    #>        Three-part model - binary logit , binary logit + fractional logit 
    #> -------------------------------------------------------------------------------- 
    #>    Estimate Std. Error z value Pr(>|z|)    
    #> x1 0.166267   0.006617   25.13   <2e-16 ***
    #> x2 0.278858   0.024929   11.19   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:49 
    #> --------------------------------------------------------------------------------

------------------------------------------------------------------------

## 2. Hypothesis Testing and Specification Diagnostics

`fracreg` includes state-of-the-art specification tests to validate your
model’s functional form and link function assumptions.

### 2.1 Generalized Goodness-Of-Functional-Form (GGOFF)

The GGOFF test (Ramalho, Ramalho and Murteira, 2014) tests whether the
chosen link function is adequate for the data by testing
$`H_0: \gamma_1 = 0, \gamma_2 = 0`$ in the augmented model:
``` math
 E(y|x) = G\left(x\beta + \gamma_1 \frac{g'(x\hat{\beta})}{g(x\hat{\beta})} + \gamma_2 x\hat{\beta} \right) 
```
where $`g(\cdot)`$ and $`g'(\cdot)`$ are the first and second
derivatives of $`G(\cdot)`$. A significant result suggests the link may
be misspecified.

``` r
# Perform the GGOFF test on the 1P model
# Tests if the logit link is appropriate
fracreg.ggoff(model_1p, version = "LM")
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    GGOFF test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit model 
    #> -------------------------------------------------------------------------------- 
    #>            Statistic p-value
    #> GOFF1 - LM     0.640   0.424
    #> GOFF2 - LM     0.485   0.486
    #> GGOFF - LM     0.814   0.666
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:49 
    #> --------------------------------------------------------------------------------

### 2.2 RESET Test

The RESET test (Ramsey, 1969) detects general functional form
misspecification by testing whether powers of the fitted values have
explanatory power:
``` math
 E(y|x) = G\left(x\beta + \sum_{k=2}^P \gamma_k (x\hat{\beta})^k\right) 
```
Testing $`H_0: \gamma = 0`$ provides a robust diagnostic check.

``` r
# Standard RESET test for functional form misspecification
fracreg.reset(model_1p)
```

Toggle to see the output

    #> 
    #> -------------------------------------------------------------------------------- 
    #>                                    RESET test 
    #> -------------------------------------------------------------------------------- 
    #> H0: Fractional logit model 
    #> -------------------------------------------------------------------------------- 
    #>       Statistic p-value
    #> LM(3)     0.852   0.653
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:50 
    #> --------------------------------------------------------------------------------

### 2.3 P-Test for Non-Nested Models

You can compare non-nested models (e.g., `logit` vs. `cloglog` link)
using the Davidson-MacKinnon (1981) P-test. If model 1 specifies
$`E_1(y|x) = G(x\beta)`$ and model 2 specifies
$`E_2(y|x) = H(z\theta)`$, the baseline model is augmented with the
difference between the fitted values:
``` math
 E(y|x) = G\left(x\beta + \gamma \left( \hat{y}_{M2} - \hat{y}_{M1} \right)\right) 
```

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
    #> -------------------------------------------------------------------------------- 
    #>                Fractional response model - cloglog specification 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                   QML 
    #> Number of observations:                                                     1000 
    #> Pseudo R-squared:                                                        0.38847 
    #> Standard errors:                                                          robust 
    #> Small sample correction:                                                   FALSE 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>          Estimate Std. Error z value Pr(>|z|)    
    #> Constant -0.77807    0.04850  -16.04   <2e-16 ***
    #> x1        0.52266    0.02813   18.58   <2e-16 ***
    #> x2        0.88575    0.07984   11.09   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:50 
    #> --------------------------------------------------------------------------------

``` r

# fracreg.ptest(model_1p, model_1p_clog) # Uncomment to run non-nested test
```

------------------------------------------------------------------------

## 3. Endogeneity & Heteroscedasticity (`fracreghet`)

When you suspect that one of your covariates is endogenous, or that the
variance is heteroscedastic,
[`fracreghet()`](https://SulmanOlieko.github.io/fracreg/reference/fracreghet.md)
provides instrumental variable correction. It natively supports IV
through a **Control Function (CF)** approach or **GMM** estimation.

### 3.1 Control Function Approach (QMLxv)

The control function approach adds the first-stage residual as an
additional regressor in the fractional model, providing a Hausman-type
test for endogeneity. For an endogenous continuous regressor $`y_{2i}`$,
the first stage is:
``` math
 y_{2i} = z_i \pi + v_i 
```
The residuals $`\hat{v}_i`$ are then included in the second-stage
fractional response model:
``` math
 E(y_{1i} \mid z_i, y_{2i}, v_i) = G(x_i \beta + \gamma \hat{v}_i) 
```

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
    #> -------------------------------------------------------------------------------- 
    #>                        Fractional logit regression model 
    #> -------------------------------------------------------------------------------- 
    #> Estimator:                                                                 QMLxv 
    #> Number of observations:                                                     1000 
    #> Standard errors:                                                          robust 
    #> Convergence:                                                          Successful 
    #> -------------------------------------------------------------------------------- 
    #>                                 Final estimates 
    #> -------------------------------------------------------------------------------- 
    #>           Estimate Std. Error z value Pr(>|z|)    
    #> Constant  -0.01289    0.01854  -0.695    0.487    
    #> x1         0.52128    0.01889  27.589   <2e-16 ***
    #> var.endog  1.61131    0.01273 126.614   <2e-16 ***
    #> vhat       0.58869    0.01273  46.259   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #>                                  Reduced form: 
    #> -------------------------------------------------------------------------------- 
    #>             Estimate Std. Error z value Pr(>|z|)    
    #> Z_INTERCEPT -0.02189    0.03149  -0.695    0.487    
    #> Z_x1         0.03615    0.03215   1.124    0.261    
    #> Z_z          1.35894    0.02938  46.259   <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:50 
    #> --------------------------------------------------------------------------------

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

The Correlated Random Effects (CRE) approach models the dependence of
unobserved individual-specific heterogeneity $`c_i`$ on covariates by
projecting it onto the time averages of strictly exogenous covariates
$`\bar{x}_i`$:
``` math
 c_i = \psi + \bar{x}_i \xi + a_i 
```
Integrating out $`a_i`$ yields the “population-averaged” or scaled
conditional mean, which can be estimated via pooled Bernoulli QML:
``` math
 E(y_{it} \mid x_i) = G(x_{it} \beta_a + \psi_a + \bar{x}_i \xi_a) 
```

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
    #>                 Estimate Std. Error z value Pr(>|z|)    
    #> x_panel         0.512705   0.008692  58.985   <2e-16 ***
    #> INTERCEPT_mean  0.013508   0.037058   0.365    0.715    
    #> x_panel_mean   -0.122045   0.081551  -1.497    0.135    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> -------------------------------------------------------------------------------- 
    #>                          Run Date: 2026-07-05 22:24:50 
    #> --------------------------------------------------------------------------------

------------------------------------------------------------------------

## Conclusion

With `fracreg`, `fracreghet`, and `fracregpd`, you have a complete,
mathematically robust toolkit for modelling fractional responses bounded
between $`[0,1]`$, regardless of inflation, endogeneity, or unobserved
panel effects.

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

For more information, please visit the [package
website](https://SulmanOlieko.github.io/fracreg) or file an issue on
[GitHub](https://github.com/SulmanOlieko/fracreg/issues).

To cite this package in your research:

``` r
citation("fracreg")
```
