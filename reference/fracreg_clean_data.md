# Clean Data for Fractional Regression Models

An internal helper to gracefully drop missing values across an arbitrary
number of vectors and matrices, replicating the functionality of
na.action = na.omit for multi-array model inputs.

## Usage

``` r
fracreg_clean_data(..., na.action = stats::na.omit)
```

## Arguments

- ...:

  A variable number of vectors, matrices, or data frames.

- na.action:

  A function specifying how to handle missing values, default is
  [`stats::na.omit`](https://rdrr.io/r/stats/na.fail.html). If `NULL`,
  no action is taken.

## Value

A named list containing the subsets of the provided arrays without
missing values.
