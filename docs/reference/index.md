# Package index

## Package Overview

General package description and help.

- [`fracreg-package`](https://sulmanolieko.github.io/fracreg/reference/fracreg-package.md)
  : Fractional Response Regressions
- [`fracreg_k401k`](https://sulmanolieko.github.io/fracreg/reference/fracreg_k401k.md)
  : 401(k) Plan Participation Data

## Base Fractional Estimators

Functions for cross-sectional, hurdle, and double-inflated models.

- [`fracreg()`](https://sulmanolieko.github.io/fracreg/reference/fracreg.md)
  : Fitting Fractional Response Regressions
- [`fracreg.pe()`](https://sulmanolieko.github.io/fracreg/reference/fracreg.pe.md)
  : Fractional Response Regressions - Partial Effects

## Panel Data & Endogeneity

Functions for panel data, unobserved heterogeneity, and endogenous
regressors.

- [`fracregpd()`](https://sulmanolieko.github.io/fracreg/reference/fracregpd.md)
  : Fitting Panel Data Fractional Response Regressions
- [`fracregpd.pe()`](https://sulmanolieko.github.io/fracreg/reference/fracregpd.pe.md)
  : Partial Effects for Fractional Panel Data Regression
- [`fracreghet()`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.md)
  : Fitting Fractional Response Regressions under Unobserved
  Heterogeneity
- [`fracreghet.pe()`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.pe.md)
  : Fractional Response Regressions under Unobserved Heterogeneity -
  Partial Effects

## Fractional Ridge Regression

Functions for applying ridge penalization to fractional regression
models.

- [`fracregridge()`](https://sulmanolieko.github.io/fracreg/reference/fracregridge.md)
  : Fractional Ridge Regression
- [`fracregridge.pe()`](https://sulmanolieko.github.io/fracreg/reference/fracregridge.pe.md)
  : Partial Effects for Fractional Ridge Regression

## Fractional Multinomial Logit

Functions for estimating fractional multinomial logit models.

- [`fracregmlogit()`](https://sulmanolieko.github.io/fracreg/reference/fracregmlogit.md)
  : Estimate Fractional Multinomial Logit Models
- [`fracregmlogit.pe()`](https://sulmanolieko.github.io/fracreg/reference/fracregmlogit.pe.md)
  : Fractional Multinomial Logit Average Partial Effects
- [`fracreg_spending`](https://sulmanolieko.github.io/fracreg/reference/fracreg_spending.md)
  : Government Spending Data
- [`print(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/summary.fracregmlogit.md)
  [`summary(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/summary.fracregmlogit.md)
  : Generate Summary Tables for fracregmlogit Objects
- [`plot(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/plot.fracregmlogit.md)
  : Plot Marginal or Discrete Effects of Willingness to Pay
- [`plot(`*`<fracregmlogit.pe>`*`)`](https://sulmanolieko.github.io/fracreg/reference/plot.fracregmlogit.pe.md)
  : Plot Marginal or Discrete Effects
- [`fitted(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/fitted.fracregmlogit.md)
  [`residuals(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/fitted.fracregmlogit.md)
  [`predict(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/fitted.fracregmlogit.md)
  : Extract Fitted Values, Residuals, and Predictions
- [`wtp()`](https://sulmanolieko.github.io/fracreg/reference/wtp.md) :
  "Willingness to Pay" for fracregmlogit models

## Specification Testing & Diagnostics

Tools for testing functional form, link specification, and non-nested
models.

- [`fracreg.ggoff()`](https://sulmanolieko.github.io/fracreg/reference/fracreg.ggoff.md)
  : GGOFF Tests for Fractional Response Regressions
- [`fracreg.reset()`](https://sulmanolieko.github.io/fracreg/reference/fracreg.reset.md)
  : RESET Test for Fractional Response Regressions
- [`fracreg.ptest()`](https://sulmanolieko.github.io/fracreg/reference/fracreg.ptest.md)
  : P Test for Fractional Response Regressions
- [`fracreghet.reset()`](https://sulmanolieko.github.io/fracreg/reference/fracreghet.reset.md)
  : RESET Test for Fractional Response Regressions under Neglected
  Heterogeneity

## Model Extractors & Utility Functions

S3 methods for extracting components from fitted models and data
cleaning.

- [`fracreg_clean_data()`](https://sulmanolieko.github.io/fracreg/reference/fracreg_clean_data.md)
  : Clean Data for Fractional Regression Models
- [`coef(`*`<fracreg>`*`)`](https://sulmanolieko.github.io/fracreg/reference/coef.fracreg.md)
  : Extract Model Coefficients for fracreg
- [`coef(`*`<fracreghet>`*`)`](https://sulmanolieko.github.io/fracreg/reference/coef.fracreghet.md)
  : Extract Model Coefficients for fracreghet
- [`coef(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/coef.fracregmlogit.md)
  : Extract Model Coefficients for fracregmlogit
- [`coef(`*`<fracregpd>`*`)`](https://sulmanolieko.github.io/fracreg/reference/coef.fracregpd.md)
  : Extract Model Coefficients for fracregpd
- [`coef(`*`<fracregridge>`*`)`](https://sulmanolieko.github.io/fracreg/reference/coef.fracregridge.md)
  : Extract Model Coefficients for fracregridge
- [`fitted(`*`<fracreg>`*`)`](https://sulmanolieko.github.io/fracreg/reference/fitted.fracreg.md)
  : Extract Fitted Values for fracreg
- [`fitted(`*`<fracreghet>`*`)`](https://sulmanolieko.github.io/fracreg/reference/fitted.fracreghet.md)
  : Extract Fitted Values for fracreghet
- [`fitted(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/fitted.fracregmlogit.md)
  [`residuals(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/fitted.fracregmlogit.md)
  [`predict(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/fitted.fracregmlogit.md)
  : Extract Fitted Values, Residuals, and Predictions
- [`fitted(`*`<fracregpd>`*`)`](https://sulmanolieko.github.io/fracreg/reference/fitted.fracregpd.md)
  : Extract Fitted Values for fracregpd
- [`fitted(`*`<fracregridge>`*`)`](https://sulmanolieko.github.io/fracreg/reference/fitted.fracregridge.md)
  : Extract Fitted Values for fracregridge
- [`residuals(`*`<fracreg>`*`)`](https://sulmanolieko.github.io/fracreg/reference/residuals.fracreg.md)
  : Extract Model Residuals for fracreg
- [`residuals(`*`<fracreghet>`*`)`](https://sulmanolieko.github.io/fracreg/reference/residuals.fracreghet.md)
  : Extract Model Residuals for fracreghet
- [`residuals(`*`<fracregpd>`*`)`](https://sulmanolieko.github.io/fracreg/reference/residuals.fracregpd.md)
  : Extract Model Residuals for fracregpd
- [`residuals(`*`<fracregridge>`*`)`](https://sulmanolieko.github.io/fracreg/reference/residuals.fracregridge.md)
  : Extract Model Residuals for fracregridge
- [`predict(`*`<fracreg>`*`)`](https://sulmanolieko.github.io/fracreg/reference/predict.fracreg.md)
  : Predict Method for fracreg
- [`predict(`*`<fracreghet>`*`)`](https://sulmanolieko.github.io/fracreg/reference/predict.fracreghet.md)
  : Predict Method for fracreghet
- [`predict(`*`<fracregpd>`*`)`](https://sulmanolieko.github.io/fracreg/reference/predict.fracregpd.md)
  : Predict Method for fracregpd
- [`predict(`*`<fracregridge>`*`)`](https://sulmanolieko.github.io/fracreg/reference/predict.fracregridge.md)
  : Predict Method for fracregridge
- [`vcov(`*`<fracreg>`*`)`](https://sulmanolieko.github.io/fracreg/reference/vcov.fracreg.md)
  : Extract Covariance Matrix for fracreg
- [`vcov(`*`<fracreghet>`*`)`](https://sulmanolieko.github.io/fracreg/reference/vcov.fracreghet.md)
  : Extract Covariance Matrix for fracreghet
- [`vcov(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/vcov.fracregmlogit.md)
  : Extract Covariance Matrix for fracregmlogit
- [`vcov(`*`<fracregpd>`*`)`](https://sulmanolieko.github.io/fracreg/reference/vcov.fracregpd.md)
  : Extract Covariance Matrix for fracregpd
- [`vcov(`*`<fracregridge>`*`)`](https://sulmanolieko.github.io/fracreg/reference/vcov.fracregridge.md)
  : Extract Covariance Matrix for fracregridge
- [`logLik(`*`<fracreg>`*`)`](https://sulmanolieko.github.io/fracreg/reference/logLik.fracreg.md)
  : Extract Log-Likelihood for fracreg
- [`logLik(`*`<fracreghet>`*`)`](https://sulmanolieko.github.io/fracreg/reference/logLik.fracreghet.md)
  : Extract Log-Likelihood for fracreghet
- [`logLik(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/logLik.fracregmlogit.md)
  : Extract Log-Likelihood for fracregmlogit
- [`logLik(`*`<fracregpd>`*`)`](https://sulmanolieko.github.io/fracreg/reference/logLik.fracregpd.md)
  : Extract Log-Likelihood for fracregpd
- [`logLik(`*`<fracregridge>`*`)`](https://sulmanolieko.github.io/fracreg/reference/logLik.fracregridge.md)
  : Extract Log-Likelihood for fracregridge
- [`nobs(`*`<fracreg>`*`)`](https://sulmanolieko.github.io/fracreg/reference/nobs.fracreg.md)
  : Extract the Number of Observations for fracreg
- [`nobs(`*`<fracreghet>`*`)`](https://sulmanolieko.github.io/fracreg/reference/nobs.fracreghet.md)
  : Extract the Number of Observations for fracreghet
- [`nobs(`*`<fracregmlogit>`*`)`](https://sulmanolieko.github.io/fracreg/reference/nobs.fracregmlogit.md)
  : Extract the Number of Observations for fracregmlogit
- [`nobs(`*`<fracregpd>`*`)`](https://sulmanolieko.github.io/fracreg/reference/nobs.fracregpd.md)
  : Extract the Number of Observations for fracregpd
- [`nobs(`*`<fracregridge>`*`)`](https://sulmanolieko.github.io/fracreg/reference/nobs.fracregridge.md)
  : Extract the Number of Observations for fracregridge
