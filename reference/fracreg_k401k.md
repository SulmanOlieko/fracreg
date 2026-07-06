# 401(k) Plan Participation Data

A cross-sectional dataset on 401(k) plan participation rates and firm
characteristics, widely used in empirical applications of fractional
response models (e.g., Papke & Wooldridge, 1996). The data is derived
from the 'wooldridge' package.

## Usage

``` r
data("fracreg_k401k")
```

## Format

A data frame with 1,534 observations and 10 variables:

- prate:

  Participation rate: proportion of eligible employees participating in
  the 401(k) plan (0 to 1).

- mrate:

  Match rate: the firm's contribution matching rate per dollar.

- totpart:

  Total number of participants.

- totelg:

  Total number of eligible employees.

- age:

  Age of the 401(k) plan in years.

- totemp:

  Total number of firm employees.

- sole:

  Indicator variable: 1 if the 401(k) is the sole retirement plan
  offered.

- ltotemp:

  Natural log of total employees.

- age_sq:

  Square of plan age.

- mrate_sq:

  Square of match rate.

## Source

Papke, L. E., & Wooldridge, J. M. (1996). "Econometric Methods for
Fractional Response Variables with an Application to 401(k) Plan
Participation Rates." Journal of Applied Econometrics, 11(6), 619-632.

## Examples

``` r
data("fracreg_k401k")
summary(fracreg_k401k$prate)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0300  0.7802  0.9570  0.8736  1.0000  1.0000 
```
