# Plot Marginal or Discrete Effects of Willingness to Pay

Plot marginal or discrete effects of willingness to pay, potentially
against another variable.

## Usage

``` r
# S3 method for class 'fracregmlogit'
plot(
  x,
  wtp.vec = NULL,
  varlist = NULL,
  against = NULL,
  mfrow = NULL,
  t = 500,
  effect = c("discrete", "marginal"),
  type = NULL,
  plot.show = TRUE,
  ...
)
```

## Arguments

- x:

  A "fracregmlogit" object.

- wtp.vec:

  A numeric vector for willingness to pay.

- varlist:

  A string vector which provides the names of variables to plot the
  effect for. If missing, all variables in the object will be plotted.

- against:

  A vector with the same length as the number of observations in the
  model, or the name of a variable. Serves as the x-axis in the plots.

- mfrow:

  A numeric vector with two elements. Specifies the number of rows and
  columns in a panel. Similar to par(mfrow=c()). Default to NULL, and
  the program will choose a square panel.

- t:

  Number of points to be used for smoothing.

- effect:

  The type of effect ("marginal" or "discrete").

- type:

  Plot type.

- plot.show:

  If TRUE, the plot will be created. Otherwise, the function returns raw
  data that can be used to create user-specified (custom) plots.

- ...:

  Additional arguments.

## Value

Panel plots of effects vs. chosen variables.

## Details

This function provides a visualisation tool for potentially
heterogeneous marginal and discrete effects of willingness to pay. The
function allows the user to plot marginal effects to detect any patterns
in the effects, in itself and against other variables. The plot also
allows visualisation of sub-groups in data, which can be very useful to
visualise categorical and dummy variables.

The function takes a `fracregmlogit` object, and internally calls `wtp`
and `fracregmlogit.pe` to compute the willingness to pay at different
data points.

Additional parameters include `varlist`, a vector of string variable
names to be plotted.

`against` allows a different variable to be chosen as the x-axis.
`against` can supply the column name of a variable in the original
dataset to plot against.

## See also

[`wtp`](https://sulmanolieko.github.io/fracreg/reference/wtp.md),
[`fracregmlogit.pe`](https://sulmanolieko.github.io/fracreg/reference/fracregmlogit.pe.md)

## Examples

``` r
data("fracreg_spending")
X = fracreg_spending[,2:5]
y = fracreg_spending[,6:11]
results1 = fracregmlogit(y, X)
#> [1] "Fractional logit model estimation completed. Time: 0.2 seconds"

# Define a willingness to pay vector
wtp.vec = c(1, 1, 1, 1, 1)

# Plot WTP for 'popdens'
plot(results1, wtp.vec=wtp.vec, varlist="popdens")
#> Error in wtp(fracregmlogit.pe(object, effect = effect, se = F, varlist = varlist,     at = newdata), wtp.vec): Wrong length of wtp.vec. Please check specification again.
```
