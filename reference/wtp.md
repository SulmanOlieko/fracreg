# "Willingness to Pay" for fracregmlogit models

Evaluate the "Willingness to Pay" given a set of arbitrary values for
outcome variables. Usually used for policy evaluations where the total
magnitude of marginal change matters.

## Usage

``` r
wtp(object, wtp.vec, varlist = NULL, indv.obs = FALSE)
```

## Arguments

- object:

  A `fracregmlogit.pe` object.

- wtp.vec:

  A numeric vector containing the arbitrary outcome values to be
  evaluated for each choice j.

- varlist:

  A string vector which provides the names of variables to calculate the
  wtp for. If missing, all variables in the object will be calculated.

- indv.obs:

  A logical value indicating whether to return individual observations.

## Value

A "fracregmlogit.wtp" object containing the estimates, standard error,
z-stats, and p-value.

## Details

This function calculates the aggregate effect of a variable on the
"willingness to pay" by linearly multiplying the average partial effect
with ex-ante (arbitrary) willingness to pay numbers associated with each
choice.

Suppose there are three choices A, B, C, each with a willingness to pay
(or cost, profit, budget), of 100, 200, and 300. The discrete effects of
variable X on A, B and C are 0.5, 0.5, and -1, with standard errors 0.2,
0.3 and 0.5. The aggregated discrete effect of X on the total
willingness to pay (or cost), is thus 100\*0.5 + 200\*0.5 + 300\*(-1) =
-150. The standard error can also be calculated to be 162.8, assuming
that the standard error is independent. A simple z-test is provided to
test whether the aggregate effect is different from zero.

Note that if the input `fracregmlogit.pe` object has no standard error
computation, then no standard error will be computed for the willingness
to pay.

## See also

[`fracregmlogit.pe`](https://sulmanolieko.github.io/fracreg/reference/fracregmlogit.pe.md),
[`plot.fracregmlogit`](https://sulmanolieko.github.io/fracreg/reference/plot.fracregmlogit.md)

## Examples

``` r
# \donttest{
data("fracreg_spending")
X = fracreg_spending[,2:5]
y = fracreg_spending[,6:11]
results1 = fracregmlogit(y, X)
pe = fracregmlogit.pe(results1)
# Assuming arbitrary WTP values for the 6 choices
wtp_est = wtp(pe, wtp.vec = c(1, 2, 3, 4, 5, 6), varlist = "houseval")
summary(wtp_est)
#> 
#> -------------------------------------------------------------------------------- 
#>                                Willingness to Pay 
#> -------------------------------------------------------------------------------- 
#>                     Fractional multinomial logit regression 
#> -------------------------------------------------------------------------------- 
#> 
#> Note: Krinsky-Robb standard error calculated
#>          Coefficient Std. Error z value Pr(>|z|)
#> houseval   -0.078693   0.120347 -0.6539   0.5132
#> 
#> -------------------------------------------------------------------------------- 
#>                          Run Date: 2026-08-08 15:50:34 
#> -------------------------------------------------------------------------------- 
#> 
# }
```
