library(testthat)
library(fracreg)

test_that("global na.omit wrapper correctly drops NAs across all models", {
  N <- 100
  set.seed(42)
  x1 <- rnorm(N)
  x2 <- runif(N)
  y_latent <- exp(x1 - x2) / (1 + exp(x1 - x2))
  y <- y_latent
  
  y[c(5, 10)] <- NA
  x1[15] <- NA
  X <- cbind(x1, x2)
  
  # fracreg
  mod1 <- fracreg(y, X, type="1P", linkfrac="logit", table=FALSE)
  expect_equal(nobs(mod1), 97)
  
  # fracreghet
  y_bounded <- pmin(pmax(y, 0.1), 0.9)
  mod2 <- fracreghet(y_bounded, X, type="GMMx", link="logit", table=FALSE)
  expect_equal(nobs(mod2), 97)
  
  # fracregridge
  mod3 <- fracregridge(y, X, fracs=0.5)
  expect_equal(nobs(mod3), 97)
  
  # fracregmlogit
  y_ml <- cbind(y, 1-y)
  mod4 <- fracregmlogit(y_ml, X)
  expect_equal(nobs(mod4), 97)
  
  # fracregpd
  id <- rep(1:20, each=5)
  time <- rep(1:5, 20)
  mod5 <- fracregpd(id, time, y, X, type="GMMpfe", link="logit", table=FALSE)
  # For panel data, length is also 97 but nobs might not be fully reliable for pd, let's check nrow(mod5$xbase)
  expect_equal(nobs(mod5), 97)
})
