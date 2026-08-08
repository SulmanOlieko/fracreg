# Fractional Multinomial Logit Average Partial Effects

Calculate average partial effects (APE) of independent variables from a
fractional multinomial logit model.

## Usage

``` r
fracregmlogit.pe(
  object,
  effect = c("marginal", "discrete"),
  marg.type = "atmean",
  se = TRUE,
  varlist = NULL,
  at = NULL,
  R = 1000
)
```

## Arguments

- object:

  A "fracregmlogit" object.

- effect:

  Can be "marginal", for marginal effects; or "discrete", for discrete
  changes from the min to the max.

- marg.type:

  Type of marginal or discrete effects to be computed. Default to
  "atmean", the effect at the mean of all covariates. Also takes
  "aveacr", the averaged effects across all observations. See details.

- se:

  Whether to calculate standard errors for those margins. Default to
  TRUE. See details.

- varlist:

  A string vector which provides the names of variables to calculate the
  marginal effect for. If missing, all variables except the (Intercept)
  will be calculated. Use "(Intercept)" if you wish to compute the
  marginal effect of the intercept.

- at:

  Specify values of the X-matrix at which the partial effect will be
  retrieved. Expects a vector input of length K-1. Only supported for
  `marg.type="atmean"`. See `predict.fracregmlogit(newdata)`.

- R:

  Number of times to sample for the Krinsky-Robb standard error. Default
  to 1000.

## Value

The function returns an object of class "fracregmlogit.pe". It contains
the following components:

`effects` A matrix of calculated effects.

`se` A matrix of standard errors corresponding to the effects. Shows up
if se=TRUE for the input parameter.

`ztable` A list of matrices containing effects, standard errors, z-stats
and p-values.

`R` Number of simulation times for Krinsky-Robb standard error
calculation. Null if se=FALSE.

`expl` String message explaining the effects calculated.

## Details

This module calculates the average partial effects (APEs) from a
fractional multinomial logit model. Partial effects are the counterpart
of the marginal effects in a linear model setting. In linear models,
usually the parameter estimate itself represents the marginal effect (if
the variable in question is continuous). In logit models, however, the
parameter estimate at hand is the effect on the log-ratio between the
choice variable and the baseline variable. This function is intended to
extract APEs from the coefficient estimates computed from the fractional
multinomial logit models.

This function allows for two types of partial effects: marginal effects,
and discrete effects. A marginal effect represents how a unit change in
one continuous variable x may influence the choice variable y. The
estimation of marginal effects is very straightforward. However, special
care is needed when averaging the marginal effect across observations to
acquire the APE. One approach is to use the estimate of the marginal
effect while setting other explanatory variables at the mean. We call
this the marginal effect at the mean (MEM), which corresponds to the
option `marg.type="atmean"`. Another approach is to take the average of
marginal effects for each individual. We call this the average marginal
effect (AME), which corresponds to the option `marg.type="aveacr"`.

The discrete effect represents how a discrete change in one specific x,
discrete or continuous, influences the choice variable y. This is more
useful for categorical variables, as calculating the "marginal effect"
makes little sense for them. In this function, we calculate the discrete
effect by changing the explanatory variable from its minimum to its
maximum. For a binary variable, this is just the difference between 0
and 1. Similar to the marginal effect case, we also have the discrete
effect at the mean (DEM), corresponding to `marg.type="atmean"` and the
average discrete effect (ADE), corresponding to `marg.type="aveacr"`.

Standard errors are provided for the effects by using the Krinsky-Robb
(KR) method. Krinsky-Robb is a simulation-based method that calculates
the empirical value of a function given a known distribution of its
variables. Here we provide Krinsky-Robb standard errors for MEM and DEM,
and the user can specify how many times of simulation `R` the
Krinsky-Robb algorithm should run.

The user can also specify a subset of explanatory variables when
calculating effects. This is done through specifying string vectors
containing the column names of the explanatory variables to `varlist`.
As the KR standard error can be computationally intensive, it is advised
to calculate it only for the variables of interest.

## See also

[`fracregmlogit`](https://sulmanolieko.github.io/fracreg/reference/fracregmlogit.md)
for the model estimation,
[`plot.fracregmlogit.pe`](https://sulmanolieko.github.io/fracreg/reference/plot.fracregmlogit.pe.md)
for plotting effects.

## Examples

``` r
data("fracreg_spending")
X = fracreg_spending[,2:5]
y = fracreg_spending[,6:11]
results1 = fracregmlogit(y, X)

# Calculate marginal effects at the mean (without standard errors for speed)
pe_marg = fracregmlogit.pe(results1, effect="marginal", se=FALSE)

# Calculate discrete effects for specific variables with standard errors
pe_disc = fracregmlogit.pe(results1, effect="discrete", 
                           varlist = colnames(results1$X)[c(1,3)], 
                           se=TRUE, R=50)
summary(pe_disc)
#> 
#> 
#> -------------------------------------------------------------------------------- 
#>                           Conditional partial effects 
#> -------------------------------------------------------------------------------- 
#>                     Fractional multinomial logit regression 
#> -------------------------------------------------------------------------------- 
#> 
#> Note: discrete effect at the mean, Krinsky-Robb standard error calculated 
#> 
#> -------------------------------------------------------------------------------- 
#>                                Choice: governing 
#> -------------------------------------------------------------------------------- 
#>                dy/dx Std. Error z value Pr(>|z|)    
#> (Intercept)  0.00000    0.00000     NaN      NaN    
#> popdens     -0.02383    0.00444  -5.368 7.98e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -------------------------------------------------------------------------------- 
#>                                  Choice: safety 
#> -------------------------------------------------------------------------------- 
#>                dy/dx Std. Error z value Pr(>|z|)    
#> (Intercept)  0.00000    0.00000     NaN      NaN    
#> popdens     -0.07750    0.01102  -7.031 2.06e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -------------------------------------------------------------------------------- 
#>                                Choice: education 
#> -------------------------------------------------------------------------------- 
#>                dy/dx Std. Error z value Pr(>|z|)  
#> (Intercept)  0.00000    0.00000     NaN      NaN  
#> popdens     -0.03498    0.01906  -1.835   0.0665 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -------------------------------------------------------------------------------- 
#>                                Choice: recreation 
#> -------------------------------------------------------------------------------- 
#>                dy/dx Std. Error z value Pr(>|z|)  
#> (Intercept)  0.00000    0.00000     NaN      NaN  
#> popdens     -0.02783    0.01308  -2.128   0.0333 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -------------------------------------------------------------------------------- 
#>                                  Choice: social 
#> -------------------------------------------------------------------------------- 
#>               dy/dx Std. Error z value Pr(>|z|)    
#> (Intercept) 0.00000    0.00000     NaN      NaN    
#> popdens     0.14039    0.03919   3.582  0.00034 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -------------------------------------------------------------------------------- 
#>                              Choice: urbanplanning 
#> -------------------------------------------------------------------------------- 
#>               dy/dx Std. Error z value Pr(>|z|)
#> (Intercept) 0.00000    0.00000     NaN      NaN
#> popdens     0.02376    0.05093   0.466    0.641
#> 
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-08-08 13:39:27 
#> -------------------------------------------------------------------------------- 
```
