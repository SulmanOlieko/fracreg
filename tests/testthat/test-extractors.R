library(testthat)
library(fracreg)

test_that("S3 Extractors return expected classes and dimensions across models", {
  N <- 100
  set.seed(42)
  X <- cbind(x1 = rnorm(N), x2 = runif(N))
  y <- exp(X[,1] - X[,2]) / (1 + exp(X[,1] - X[,2]))
  
  # 1. fracreg
  mod1 <- fracreg(y, X, type="1P", linkfrac="logit", table=FALSE)
  expect_equal(length(coef(mod1)), 3)
  expect_equal(length(fitted(mod1)), 100)
  expect_equal(length(residuals(mod1)), 100)
  expect_equal(nobs(mod1), 100)
  expect_true(is.numeric(logLik(mod1)))
  expect_equal(length(predict(mod1)), 100)
  
  # 2. fracregpd
  id <- rep(1:20, each=5)
  time <- rep(1:5, 20)
  mod2 <- fracregpd(id, time, y, X, type="GMMpfe", link="logit", table=FALSE)
  expect_equal(length(coef(mod2)), 2)
  expect_equal(length(fitted(mod2)), 100)
  expect_equal(length(residuals(mod2)), 100)
  expect_equal(nobs(mod2), 100)
  expect_true(is.na(logLik(mod2)))
  expect_equal(length(predict(mod2)), 100)
  
  # Bounded y to prevent intermediate glm(Gamma) optimization NaNs 
  y_bounded <- pmin(pmax(y, 0.1), 0.9)

  # 3. fracreghet
  mod3 <- fracreghet(y_bounded, X, type="GMMx", link="logit", table=FALSE)
  expect_equal(length(coef(mod3)), 3)
  expect_equal(length(fitted(mod3)), 100)
  expect_equal(length(residuals(mod3)), 100)
  expect_equal(nobs(mod3), 100)
  expect_true(is.na(logLik(mod3)))
  expect_equal(length(predict(mod3)), 100)
  
  # 4. fracregridge
  mod4 <- fracregridge(y, X, fracs=0.5)
  expect_true(is.array(coef(mod4)))
  expect_true(is.array(fitted(mod4)))
  expect_true(is.array(residuals(mod4)))
  expect_equal(nobs(mod4), 100)
  
  # 5. fracregmlogit
  y_ml <- cbind(y, 1-y)
  mod5 <- fracregmlogit(y_ml, X)
  expect_true(is.matrix(coef(mod5)))
  expect_s3_class(fitted(mod5), "data.frame")
  expect_s3_class(residuals(mod5), "data.frame")
  expect_s3_class(predict(mod5), "data.frame")
  expect_equal(nobs(mod5), 100)
  expect_true(is.numeric(logLik(mod5)))
})
