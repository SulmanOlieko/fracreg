# Government Spending Data

Spending on different categories by Dutch cities in 2005. This dataset
is commonly used to demonstrate fractional multinomial logit models.

## Usage

``` r
data("fracreg_spending")
```

## Format

A data frame with 429 observations and 12 variables:

- muni:

  Name of municipality

- houseval:

  Average value of a house in 100,000 euros

- popdens:

  Population density in 1000s of persons per square km

- noleft:

  No left party in city government

- minorityleft:

  Minority left party in city government

- governing:

  Fraction of spending on governing

- safety:

  Fraction of spending on safety

- education:

  Fraction of spending on education

- recreation:

  Fraction of spending on recreation

- social:

  Fraction of spending on social services

- urbanplanning:

  Fraction of spending on urban planning

- tot:

  Total spending or population

## Source

\<http://fmwww.bc.edu/repec/bocode/c/citybudget.dta\>

## Examples

``` r
data("fracreg_spending")
head(fracreg_spending)
#>               muni houseval popdens noleft minorityleft  governing    safety
#> 1    's-Gravendeel     1.28   0.472      0            0 0.14832154 0.1838828
#> 2    's-Gravenhage     1.02   5.711      0            1 0.03985741 0.1016635
#> 3 's-Hertogenbosch     1.43   1.583      0            1         NA        NA
#> 4      Aa en Hunze     1.41   0.091      0            1 0.04114016 0.1750340
#> 5          Aalburg     1.59   0.241      1            0 0.09675697 0.2273455
#> 6         Aalsmeer     1.89   1.128      1            0 0.04417612 0.1001368
#>    education recreation    social urbanplanning      tot
#> 1 0.08417173 0.14274330 0.1531029     0.2877777   1.0039
#> 2 0.10785559 0.09804828 0.4179285     0.2346468 184.8314
#> 3         NA         NA        NA            NA   0.0000
#> 4 0.27075234 0.13644537 0.2220305     0.1545976   4.1152
#> 5 0.18297078 0.09675697 0.1705592     0.2256106   1.4986
#> 6 0.03706029 0.07046607 0.1161555     0.6320052   6.2115
```
