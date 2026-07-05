library(testthat)
library(fracreg)

# Prepare data for testing
data("fracreg_k401k")
y <- fracreg_k401k$prate
X <- cbind(
  mrate = fracreg_k401k$mrate, 
  age = fracreg_k401k$age, 
  totemp = fracreg_k401k$totemp, 
  sole = fracreg_k401k$sole
)

test_that("fracreg.reset works", {
  model_1p <- fracreg(y = y, x = X, type = "1P", linkfrac = "logit")
  reset_res <- fracreg.reset(model_1p, table = FALSE)
  
  expect_type(reset_res, "double")
  expect_true(length(reset_res) > 0)
})

test_that("fracreg.ggoff works", {
  model_1p <- fracreg(y = y, x = X, type = "1P", linkfrac = "logit")
  ggoff_res <- fracreg.ggoff(model_1p, version = "LM", table = FALSE)
  
  expect_type(ggoff_res, "double")
  expect_true(length(ggoff_res) > 0)
})

test_that("fracreg.ptest works", {
  model_1p_logit <- fracreg(y = y, x = X, type = "1P", linkfrac = "logit")
  model_1p_clog <- fracreg(y = y, x = X, type = "1P", linkfrac = "cloglog")
  
  ptest_res <- fracreg.ptest(model_1p_logit, model_1p_clog, table = FALSE)
  
  expect_type(ptest_res, "double")
  expect_true(length(ptest_res) > 0)
})
