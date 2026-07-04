## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
devtools::load_all()


## ----data-sim-----------------------------------------------------------------
set.seed(123)
N <- 1000
x1 <- rnorm(N)
x2 <- runif(N)

# Generating a fractional dependent variable with inflation at 0 and 1
XB <- -0.5 + 0.8 * x1 + 1.2 * x2 + rnorm(N)
y_latent <- exp(XB) / (1 + exp(XB))

y <- y_latent
# Inflate at boundaries
y[y_latent < 0.2] <- 0
y[y_latent > 0.8] <- 1

data <- data.frame(y = y, x1 = x1, x2 = x2)
head(data)


## ----one-part-----------------------------------------------------------------
model_1p <- fracreg(
  y = data$y, 
  x = cbind(x1 = data$x1, x2 = data$x2), 
  type = "1P", 
  linkfrac = "logit"
)
summary(model_1p)


## ----pe-1p--------------------------------------------------------------------
pe_1p <- fracreg.pe(model_1p)
summary(pe_1p)


## ----two-part-----------------------------------------------------------------
# We use inflation=0 to indicate we are modeling the mass at zero
model_2p <- fracreg(
  y = data$y, 
  x = cbind(x1 = data$x1, x2 = data$x2), 
  type = "2P", 
  inflation = 0,
  linkbin = "logit", 
  linkfrac = "logit"
)
summary(model_2p)

# Calculate partial effects for the combined two-part model
pe_2p <- fracreg.pe(model_2p)
summary(pe_2p)


## ----three-part---------------------------------------------------------------
# Linkbin takes a vector of two links for the two binary hurdles
model_3p <- fracreg(
  y = data$y, 
  x = cbind(x1 = data$x1, x2 = data$x2), 
  type = "3P", 
  linkbin = c("logit", "logit"), 
  linkfrac = "logit"
)
summary(model_3p)

# The analytical delta method is natively supported for the 3P model!
pe_3p <- fracreg.pe(model_3p)
summary(pe_3p)


## ----ggoff--------------------------------------------------------------------
# Perform the GGOFF test on the 1P model
# Tests if the logit link is appropriate
fracreg.ggoff(model_1p, version = "LM")


## ----reset--------------------------------------------------------------------
# Standard RESET test for functional form misspecification
fracreg.reset(model_1p)


## ----ptest--------------------------------------------------------------------
model_1p_clog <- fracreg(y = data$y, x = cbind(x1 = data$x1, x2 = data$x2), type = "1P", linkfrac = "cloglog")
# fracreg.ptest(model_1p, model_1p_clog) # Uncomment to run non-nested test


## ----endogeneity--------------------------------------------------------------
# Simulating an endogenous variable (var.endog) and an instrument (z)
z <- rnorm(N)
u <- 0.5 * z + rnorm(N)
var.endog <- 0.8 * z + u
y_endog <- exp(0.5 * x1 + 1.2 * var.endog + u) / (1 + exp(0.5 * x1 + 1.2 * var.endog + u))

# Estimate QML with control function for cross-sectional data
model_het <- fracreghet(
  y = y_endog, 
  x = cbind(x1, var.endog), 
  z = cbind(x1, z), 
  var.endog = var.endog,
  type = "QMLxv", 
  link = "logit"
)


## ----panel-data---------------------------------------------------------------
# Simulating Panel Data
id <- rep(1:200, each = 5)
time <- rep(1:5, times = 200)
x_panel <- rnorm(1000)
# Unobserved individual effect (CRE)
c_i <- rep(rnorm(200), each = 5) 
y_panel <- exp(x_panel + c_i) / (1 + exp(x_panel + c_i))

# Estimate a Correlated Random Effects (CRE) Model
model_pd <- fracregpd(
  id = id, 
  time = time, 
  y = y_panel, 
  x = cbind(x_panel = x_panel), 
  type = "QMLcre", 
  link = "logit"
)

