# Generate Summary Tables for fracregmlogit Objects

Generate tables of coefficient estimates, partial effects, and
willingness to pay from fracregmlogit-type objects.

## Usage

``` r
# S3 method for class 'fracregmlogit.wtp'
summary(object, ...)
```

## Arguments

- object:

  an object with class "fracregmlogit", "fracregmlogit.pe", or
  "fracregmlogit.wtp".

- ...:

  Additional arguments passed to the printCoefmat function.

## Value

Returns the object invisibly.

## Details

This module provides summary methods for three fracregmlogit objects:
`fracregmlogit`, `fracregmlogit.pe` , and `fracregmlogit.wtp`.

For `fracregmlogit` objects, the summary prints the number of
observations, log pseudo-likelihood, baseline choice, and the
coefficient estimates with standard errors, z-statistics, and p-values
for each choice equation.

For `fracregmlogit.pe` objects, it displays the marginal or discrete
effects along with their computed standard errors (if Krinsky-Robb
sampling was performed) for each choice.

For `fracregmlogit.wtp` objects, it provides a table of the aggregated
willingness to pay along with its standard errors and test statistics.

## See also

[`fracregmlogit`](https://sulmanolieko.github.io/fracreg/reference/fracregmlogit.md),
[`fracregmlogit.pe`](https://sulmanolieko.github.io/fracreg/reference/fracregmlogit.pe.md)

## Examples

``` r
data("fracreg_spending")
X = fracreg_spending[,2:5]
y = fracreg_spending[,6:11]

# generate fracregmlogit summary
results1 = fracregmlogit(y, X)
#> [1] "Fractional logit model estimation completed. Time: 0.2 seconds"
summary(results1)
#> 
#> -------------------------------------------------------------------------------- 
#>                        Fractional multinomial logit model 
#> -------------------------------------------------------------------------------- 
#> Number of observations:                                                      392 
#> Baseline choice:                                                       governing 
#> Log pseudolikelihood:                                                  -673.1203 
#> Standard errors:                                                          robust 
#> 
#> -------------------------------------------------------------------------------- 
#>                                  Choice: safety 
#> -------------------------------------------------------------------------------- 
#>              Coefficient Std.Err. z value Pr(>|z|)    
#> houseval        -0.14001  0.03712  -3.772 0.000162 ***
#> popdens          0.01158  0.01875   0.618 0.536751    
#> noleft           0.08254  0.04572   1.805 0.071049 .  
#> minorityleft     0.18936  0.04450   4.255 2.09e-05 ***
#> constant         0.74898  0.06527  11.475  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -------------------------------------------------------------------------------- 
#>                                Choice: education 
#> -------------------------------------------------------------------------------- 
#>              Coefficient Std.Err. z value Pr(>|z|)    
#> houseval        -0.63715  0.10737  -5.934 2.96e-09 ***
#> popdens          0.09276  0.03044   3.048  0.00231 ** 
#> noleft          -0.36480  0.09164  -3.981 6.87e-05 ***
#> minorityleft     0.03874  0.09353   0.414  0.67875    
#> constant         1.21527  0.16536   7.349 1.99e-13 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -------------------------------------------------------------------------------- 
#>                                Choice: recreation 
#> -------------------------------------------------------------------------------- 
#>              Coefficient Std.Err. z value Pr(>|z|)    
#> houseval        -0.23088  0.03969  -5.817 5.98e-09 ***
#> popdens          0.07204  0.01577   4.569 4.89e-06 ***
#> noleft           0.01385  0.04307   0.322    0.748    
#> minorityleft     0.22266  0.04195   5.307 1.11e-07 ***
#> constant         0.42086  0.06640   6.339 2.32e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -------------------------------------------------------------------------------- 
#>                                  Choice: social 
#> -------------------------------------------------------------------------------- 
#>              Coefficient Std.Err. z value Pr(>|z|)    
#> houseval        -0.62082  0.06543  -9.488   <2e-16 ***
#> popdens          0.19818  0.01998   9.917   <2e-16 ***
#> noleft          -0.14671  0.05997  -2.446   0.0144 *  
#> minorityleft     0.13606  0.05900   2.306   0.0211 *  
#> constant         1.70671  0.11044  15.453   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -------------------------------------------------------------------------------- 
#>                              Choice: urbanplanning 
#> -------------------------------------------------------------------------------- 
#>              Coefficient Std.Err. z value Pr(>|z|)    
#> houseval        -0.17859  0.07388  -2.417  0.01564 *  
#> popdens          0.16048  0.03381   4.746 2.07e-06 ***
#> noleft           0.03022  0.08359   0.361  0.71773    
#> minorityleft     0.23444  0.07791   3.009  0.00262 ** 
#> constant         0.98183  0.12528   7.837 4.66e-15 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-11 22:01:39 
#> -------------------------------------------------------------------------------- 

# generate marginal effects summary
effects1 = fracregmlogit.pe(results1, effect="marginal", se=FALSE)
summary(effects1)
#> 
#> -------------------------------------------------------------------------------- 
#>             marginal effect at the mean, standard error not computed 
#> -------------------------------------------------------------------------------- 
#> Effects:
#>                  houseval      popdens       noleft minorityleft
#> governing      0.03034691 -0.010268477  0.005080197 -0.014758737
#> safety         0.03319717 -0.017802902  0.025010114  0.006114876
#> education     -0.03693917 -0.001846496 -0.036332555 -0.013701221
#> recreation     0.01074902 -0.004363436  0.008070371  0.007956663
#> social        -0.07220130  0.021568719 -0.022384532 -0.004810849
#> urbanplanning  0.03484736  0.012712591  0.020556405  0.019199268
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-07-11 22:01:39 
#> -------------------------------------------------------------------------------- 
```
