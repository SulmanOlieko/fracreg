\name{fracregpd}
\alias{fracregpd}

\title{
Fitting Panel Data Fractional Regression Models
}

\description{
\code{fracregpd} is used to fit panel data regression models when the dependent variable has a bounded, fractional nature.
}
\usage{
fracregpd(id, time, y, x, z, var.endog, x.exogenous = T, lags, start, type,
	GMMww.cor = T, link = "logit", intercept = T, table = T, variance = T,
	var.type = "cluster", tdummies = F, bootstrap = F, B = 200, ...)
}

\arguments{
  \item{id}{
a numeric vector identifying the cross-sectional units.
}
  \item{time}{
a numeric vector identifying the time periods in which the cross-sectional units were observed.
}
  \item{y}{
a numeric vector containing the values of the response variable.
}
  \item{x}{
a numeric matrix, with column names, containing the values of all covariates (exogenous and endogenous).
}
  \item{z}{
a numeric matrix, with column names, containing the values of all exogenous variables (covariates and external
instrumental variables). Only required in case of endogenous explanatory variables.
}
  \item{var.endog}{
a numeric vector containing the values of the endogenous covariate (or of some transformation of it), which will 
be used as dependent variable in the linear reduced form assumed for application of the \code{QMLcre} estimator. 
Only required for this estimator.
}
  \item{x.exogenous}{
a logical value indicating whether all explanatory variables are assumed to be exogenous or not.
}
  \item{lags}{
a logical value indicating whether the first lags of \code{x} or  \code{z} should be used as instruments for \code{x}. 
Defaults to \code{TRUE} for the GMMww and GMMc estimators and to \code{FALSE} for the remaining estimators. 
The \code{GMMcre} and \code{QMLcre} estimators do not admit lagged instruments.
}
  \item{start}{
a numeric vector containing the initial values for the parameters to be optimized. Optional.
}
  \item{type}{
a description of the estimator to compute: \code{GMMww}, \code{GMMc}, 
\code{GMMbgw}, \code{GMMpfe}, \code{GMMcre}, \code{GMMpre} or \code{QMLcre}.
}
  \item{GMMww.cor}{
a logical value indicating whether each explanatory variable should be transformed in deviations from its overall 
mean before computing the \code{GMMww} estimator.
}
  \item{link}{
a description of the link function to use. Available options for all GMM estimators: \code{logit} and \code{cloglog}. 
Only option for the \code{QMLcre} estimator: \code{probit}.
}
  \item{intercept}{
a logical value indicating whether the model should include a constant term or not. Only relevant for the 
\code{GMMpre} estimator.
}
  \item{table}{
a logical value indicating whether a summary table with the regression results should be printed.
}
  \item{variance}{
a logical value indicating whether the variance of the estimated parameters should be calculated. Defaults to 
\code{TRUE} whenever \code{table = TRUE}.
}
  \item{var.type}{
a description of the type of variance of the estimated parameters to be calculated. Options are \code{cluster}, the 
default, and \code{robust}. In overidentified models, it also affects the parameter estimates via the GMM weighting matrix.
}
  \item{tdummies}{
a logical value indicating whether time dummies should be included among the model explanatory variables.
}
  \item{bootstrap}{
a logical value indicating whether bootstrap should be used in the estimation of the parameter standard errors.
}
  \item{B}{
the number of bootstrap replications.
}
  \item{\dots}{
Arguments to pass to \link[stats]{nlminb}.
}
}

\details{
\code{fracregpd} computes the GMM estimators proposed in Ramalho and Ramalho (2016) 
for panel data fractional regression models with both time-variant and time-invariant
unobserved heterogeneity and endogeneous covariates: GMMww, GMMc, GMMbgw, GMMpfe, GMMcre
and GMMpre. In addition, \code{fracregpd} also computes QMLcre, which was proposed by Papke
and Wooldridge (2008) and Wooldridge (2010). For overidentified models, \code{fracregpd}
calculates Hansen's J statistic.
}

\value{
\code{fracregpd} returns a list with the following elements:
  \item{type}{the name of the estimator computed.
}
  \item{link}{the name of the specified link.
}
  \item{p}{a named vector of coefficients.
}
  \item{Hy}{the transformed values of the response variable when GMM estimators are computed or the 
   values of the response variable in the QML case.
}
  \item{converged}{logical. Was the algorithm judged to have converged?
}

In case of an overidentifying model, the following element is also returned:
  \item{J}{the result of Hansen's J test of overidentifying moment conditions.
}

If \code{variance = TRUE} or \code{table = TRUE} and the algorithm converged successfully, the previous list also contains the following elements:
  \item{p.var}{a named covariance matrix.
}
  \item{var.type}{covariance matrix type.
}
}

\references{
Papke, L. and Wooldridge, J.M. (2008), "Panel data methods for fractional response variables with an application to test pass rates", \emph{Journal of Econometrics}, 145(1-2), 121-233.

Ramalho, E.A. and J.J.S. Ramalho (2016), "Exponential regression of fractional-response fixed-efects models with an application to firm capital structure", mimeo.

Wooldridge, J.M. (2010), "Correlated random effects models with unbalanced panels", mimeo.
}

\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}

\seealso{
\code{\link{fracreg}}, for fitting standard cross-sectional fractional regression models.\cr
\code{\link{fracreghet}}, for fitting cross-sectional fractional regression models with unobserved heterogeneity.
}

\examples{
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

}
