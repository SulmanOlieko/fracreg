pkgname <- "fracreg"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "fracreg-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('fracreg')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("fracreg")
### * fracreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fracreg
### Title: Fitting Fractional Regression Models
### Aliases: fracreg

### ** Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#fracreg estimation of a logit fractional regression model
fracreg(y,X,linkfrac="logit")

#fracreg estimation of the binary logit component of the two-part fractional
#regression model with y=1 as the relevant boundary value
fracreg(y,X,linkbin="logit",type="2Pbin",inf=1)

#fracreg estimation of the fractional component of the two-part fractional
#regression model with y=1 as the relevant boundary value and using a
#probit link function
fracreg(y,X,linkfrac="probit",type="2Pfrac",inf=1)

#fracreg estimation of both components of a two-part fractional regression model
#with y=1 as the relevant boundary value and using a cloglog binary link
#function and a logit fractional link function
fracreg(y,X,linkbin="cloglog",linkfrac="logit",type="2P",inf=1)

#Three-part model (y has both 0s and 1s)
y3p <- y
y3p[1:20] <- 0
y3p[21:40] <- 1
fracreg(y3p,X,linkbin=c("logit","probit"),linkfrac="logit",type="3P")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fracreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fracreg.ggoff")
### * fracreg.ggoff

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fracreg.ggoff
### Title: GGOFF Tests for Fractional Regression Models
### Aliases: fracreg.ggoff

### ** Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Testing the logit specification of a standard fractional regression model
#using LM and Wald versions of the GGOFF test, based on 1 or 2 fitted powers of
#the linear predictor
res <- fracreg(y,X,linkfrac="logit",table=FALSE)
fracreg.ggoff(res,c("Wald","LM"))

#Testing the probit specification of the binary component of a two-part fractional
#regression model using a LR-based GGOFF test
res <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1,table=FALSE)
fracreg.ggoff(res,"LR")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fracreg.ggoff", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fracreg.pe")
### * fracreg.pe

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fracreg.pe
### Title: Fractional Regression Models - Partial Effects
### Aliases: fracreg.pe

### ** Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Computing average partial effects for a logit fractional regression model
res <- fracreg(y,X,linkfrac="logit",table=FALSE)
fracreg.pe(res)

#Computing average partial effects for a binary logit + fractional probit
#two-part model
res <- fracreg(y,X,linkbin="logit",linkfrac="probit",type="2P",inf=1,table=FALSE)
fracreg.pe(res)

#Computing conditional partial effects for X2 in the logit component
#of a two-part fractional regression model, with the covariates evaluated
#at their median values
res <- fracreg(y,X,linkfrac="logit",type="2Pfrac",inf=1,table=FALSE)
fracreg.pe(res,APE=FALSE,CPE=TRUE,at="median",which.x="X2")

#Computing average partial effects for a three-part double-inflated model
y3p <- y
y3p[1:20] <- 0
y3p[21:40] <- 1
res3p <- fracreg(y3p,X,linkbin=c("logit","probit"),linkfrac="logit",type="3P",table=FALSE)
fracreg.pe(res3p)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fracreg.pe", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fracreg.ptest")
### * fracreg.ptest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fracreg.ptest
### Title: P Test for Fractional Regression Models
### Aliases: fracreg.ptest

### ** Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Testing logit versus loglog specifications for standard fractional
#regression models using a LM version of the P test
res1 <- fracreg(y,X,linkfrac="logit",table=FALSE)
res2 <- fracreg(y,X,linkfrac="loglog",table=FALSE)
fracreg.ptest(res1,res2,"LM")

#Testing a logit one-part fractional regression model versus a binary logit +
#fractional probit two-part model using a Wald version of the P test
res1 <- fracreg(y,X,linkfrac="logit",table=FALSE)
res2 <- fracreg(y,X,linkbin="logit",linkfrac="probit",type="2P",inf=1,table=FALSE)
fracreg.ptest(res1,res2,"Wald")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fracreg.ptest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fracreg.reset")
### * fracreg.reset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fracreg.reset
### Title: RESET Test for Fractional Regression Models
### Aliases: fracreg.reset

### ** Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Testing the logit specification of a standard fractional regression model
#using LM and Wald versions of the RESET test, based on 1 or 2 fitted powers of
#the linear predictor
res <- fracreg(y,X,linkfrac="logit",table=FALSE)
fracreg.reset(res,2:3,c("Wald","LM"))

#Testing the probit specification of the binary component of a two-part fractional
#regression model using LR-based RESET tests with quadratic and cubic fitted 
#powers of the linear predictor
res <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1,table=FALSE)
fracreg.reset(res,3,"LR")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fracreg.reset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fracreghet")
### * fracreghet

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fracreghet
### Title: Fitting Fractional Regression Models under Unobserved
###   Heterogeneity
### Aliases: fracreghet

### ** Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

Z <- cbind(rnorm(N),rnorm(N),rnorm(N))
dimnames(Z)[[2]] <- c("Z1","Z2","Z3")

y <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))

#Exogeneity, GMMx estimator
fracreghet(y,X,type="GMMx")

#Endogeneity, GMMz estimator
fracreghet(y,X,Z,type="GMMz")

#Endogeneity, GMMxv estimator
fracreghet(y,X,Z,X[,1],type="GMMxv")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fracreghet", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fracreghet.pe")
### * fracreghet.pe

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fracreghet.pe
### Title: Fractional Regression Models under Unobserved Heterogeneity -
###   Partial Effects
### Aliases: fracreghet.pe

### ** Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

Z <- cbind(rnorm(N),rnorm(N),rnorm(N))
dimnames(Z)[[2]] <- c("Z1","Z2","Z3")

y <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))

res <- fracreghet(y,X,type="GMMx",table=FALSE)

#Smearing estimator of average partial effects for variable X1
fracreghet.pe(res,which.x="X1")

#Naive estimator of conditional partial effects for all covariates,
#which are evaluated at X1=1 and X2=-1
fracreghet.pe(res,smearing=FALSE,APE=FALSE,CPE=TRUE,at=c(1,-1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fracreghet.pe", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fracreghet.reset")
### * fracreghet.reset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fracreghet.reset
### Title: RESET Test for Fractional Regression Models under Neglected
###   Heterogeneity
### Aliases: fracreghet.reset

### ** Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

Z <- cbind(rnorm(N),rnorm(N),rnorm(N))
dimnames(Z)[[2]] <- c("Z1","Z2","Z3")

y <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))

res <- fracreghet(y,X,type="GMMx",table=FALSE)

#LM and Wald versions of the RESET test, based on 1 or 2 fitted powers of xb
fracreghet.reset(res,2:3,c("Wald","LM"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fracreghet.reset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fracregpd")
### * fracregpd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fracregpd
### Title: Fitting Panel Data Fractional Regression Models
### Aliases: fracregpd

### ** Examples

id <- rep(1:50,each=5)
time <- rep(1:5,50)
NT <- 250

XBu <- rnorm(NT)
y <- exp(XBu)/(1+exp(XBu))

X <- cbind(rnorm(NT),rnorm(NT))
dimnames(X)[[2]] <- c("X1","X2")

Z <- cbind(rnorm(NT),rnorm(NT),rnorm(NT))
dimnames(Z)[[2]] <- c("Z1","Z2","Z3")

# Exogeneity, no lags, no time dummies, clustered standard errors, GMMbgw estimator
fracregpd(id,time,y,X,type="GMMbgw")

# Lagged covariates and instruments, robust standard errors, GMMww estimator
fracregpd(id,time,y,X,lags=TRUE,type="GMMww",var.type="robust")

# Endogeneity, time dummies, GMMpfe estimator
fracregpd(id,time,y,X,Z,x.exogenous=FALSE,type="GMMpfe",tdummies=TRUE)

# Standard errors based on 100 bootstrap samples, QMLcre estimator (not run)
#fracregpd(id,time,y,X,Z,X[,1],x.exogenous=FALSE,type="QMLcre",link="probit",bootstrap=TRUE,B=100)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fracregpd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
