# "Willingness to Pay" for fracregmlogit models

Calculates the willingness to pay for fractional multinomial logit
models.

## Usage

``` r
wtp(object, wtp.vec, varlist = NULL, indv.obs = FALSE)
```

## Arguments

- object:

  A fracregmlogit.pe object.

- wtp.vec:

  A numeric vector containing the outcome arbitrary values to be
  evaluated.

- varlist:

  A string or vector of strings of covariates to be evaluated. Default
  is all available explanatory variables.

- indv.obs:

  boolean. Provide WTP values for individual observations?

## Value

A matrix containing the estimates, standard error, z-stats, and p-value.

A "fracregmlogit.wtp" object.

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
