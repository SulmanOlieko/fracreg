library(testthat)
library(fracreg)

test_that("fracreghet QMLxv estimator works", {
  set.seed(42)
  N <- 500
  x1 <- rnorm(N)
  z <- rnorm(N)
  u <- 0.5 * z + rnorm(N)
  var.endog <- 0.8 * z + u
  xb <- 0.5 * x1 + 1.2 * var.endog + u
  y_endog <- exp(xb) / (1 + exp(xb))
  
  X <- cbind(x1 = x1, var.endog = var.endog)
  Z <- cbind(x1 = x1, z = z)
  
  model_het <- fracreghet(
    y = y_endog, 
    x = X, 
    z = Z, 
    var.endog = var.endog,
    type = "QMLxv", 
    link = "logit",
    table = FALSE
  )
  
  expect_type(model_het, "list")
  expect_true(all(c("p", "p.var", "Hy", "xbhat") %in% names(model_het)))
})
