library(testthat)
library(fracreg)

# Prepare data for testing
data("fracreg_k401k")
y <- fracreg_k401k$prate
X <- as.matrix(cbind(
  mrate = fracreg_k401k$mrate, 
  age = fracreg_k401k$age, 
  totemp = fracreg_k401k$totemp, 
  sole = fracreg_k401k$sole
))

test_that("fracreg 1P estimator works", {
  model_1p <- fracreg(y = y, x = X, type = "1P", linkfrac = "logit", table = FALSE)
  
  expect_type(model_1p, "list")
  expect_true(all(c("p", "p.var", "yhat") %in% names(model_1p)))
  expect_false(any(is.na(model_1p$p)))
  expect_equal(length(model_1p$p), ncol(X) + 1)
})

test_that("fracreg 2P estimator works", {
  # Modify y slightly to ensure some 0s for the 2P test
  y_2p <- y
  y_2p[1:50] <- 0
  
  model_2p <- fracreg(y = y_2p, x = X, type = "2P", inflation = 0, linkbin = "logit", linkfrac = "logit", table = FALSE)
  
  expect_type(model_2p, "list")
  expect_true(all(c("resBIN", "resFRAC", "yhat2P") %in% names(model_2p)))
  expect_false(any(is.na(model_2p$resBIN$p)))
  expect_false(any(is.na(model_2p$resFRAC$p)))
})

test_that("fracreg 3P estimator works", {
  # Modify y slightly to ensure some 0s and 1s for the 3P test
  y_3p <- y
  y_3p[1:50] <- 0
  y_3p[51:100] <- 1
  
  model_3p <- fracreg(y = y_3p, x = X, type = "3P", linkbin = c("logit", "logit"), linkfrac = "logit", table = FALSE)
  
  expect_type(model_3p, "list")
  expect_true(all(c("resBIN0", "resBIN1", "resFRAC", "yhat3P") %in% names(model_3p)))
})

test_that("fracreg.pe works for 1P", {
  model_1p <- fracreg(y = y, x = X, type = "1P", linkfrac = "logit", table = FALSE)
  pe_1p <- fracreg.pe(model_1p, table = FALSE)
  
  expect_type(pe_1p, "list")
  expect_true(all(c("PE.p", "PE.sd") %in% names(pe_1p)))
  expect_equal(length(pe_1p$PE.p), ncol(X))
})
