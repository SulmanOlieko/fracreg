library(testthat)
library(fracreg)

test_that("fracregpd QMLcre estimator works", {
  set.seed(42)
  N <- 100
  T <- 5
  id <- rep(1:N, each = T)
  time <- rep(1:T, times = N)
  x_panel <- rnorm(N * T)
  c_i <- rep(rnorm(N), each = T) 
  xb <- x_panel + c_i
  y_panel <- exp(xb) / (1 + exp(xb))
  
  X <- cbind(x_panel = x_panel)
  
  model_pd <- fracregpd(
    id = id, 
    time = time, 
    y = y_panel, 
    x = X, 
    type = "QMLcre", 
    link = "probit",
    table = FALSE
  )
  
  expect_type(model_pd, "list")
  expect_true(all(c("p", "p.var", "Hy") %in% names(model_pd)))
})
