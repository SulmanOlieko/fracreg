# Predict Method for fracreg

Predicts conditional mean values from a fitted `fracreg` model.

## Usage

``` r
# S3 method for class 'fracreg'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted model object of class `fracreg`.

- newdata:

  An optional data frame or matrix in which to look for variables with
  which to predict. If omitted, the fitted values are used.

- ...:

  Further arguments passed to or from other methods.

## Value

A numeric vector of predicted values.
