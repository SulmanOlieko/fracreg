# Extract Fitted Values, Residuals, and Predictions

Extract Fitted Values, Residuals, and Predictions

## Usage

``` r
# S3 method for class 'fracregmlogit'
fitted(object, ...)

# S3 method for class 'fracregmlogit'
residuals(object, ...)

# S3 method for class 'fracregmlogit'
predict(object, newdata = NULL, newbeta = NULL, ...)
```

## Arguments

- object:

  A "fracregmlogit" object.

- ...:

  Additional arguments.

- newdata:

  A new X matrix to perform model prediction. If NULL, defaults to the
  original dataset. X can be a vector with length k, or a matrix with k
  columns, where k is the number of explanatory variables in the
  original model.

- newbeta:

  A new augmented matrix of coefficients that can be used to predict
  outcome variables. Feeds into `object$coefficient`, which contains the
  baseline coefficient. Useful for constructing confidence intervals via
  simulation or bootstrapping.

## Value

An object of class `data.frame` containing numeric values where each
column corresponds to one of the choice alternatives in the response
variable matrix and each row corresponds to an observation.
Specifically:

- `fitted`: Returns the estimated fitted fractional response values
  (choice shares or predicted conditional probabilities).

- `residuals`: Returns the response residuals (the actual observed
  shares minus the estimated fitted shares).

- `predict`: Returns the predicted choice shares or conditional
  probabilities computed from the specified model object and `newdata`
  or `newbeta`.

## See also

[`fracregmlogit`](https://sulmanolieko.github.io/fracreg/reference/fracregmlogit.md)

## Examples

``` r
data("fracreg_spending")
df <- na.omit(fracreg_spending)
X = df[,2:5]
y = df[,6:11]
results1 = fracregmlogit(y, X)
#> [1] "Fractional logit model estimation completed. Time: 0.3 seconds"

# Extract fitted values
fit = fitted(results1)

# Extract residuals
res = residuals(results1)

# Predict using the first observation from the original dataset
pred = predict(results1, newdata = X[1,])
```
