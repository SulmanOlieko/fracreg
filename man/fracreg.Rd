\name{fracreg}
\alias{fracreg}

\title{
Fitting Fractional Regression Models
}

\description{
\code{fracreg} is used to fit fractional regression models, which are appropriate for responses that are proportions, percentages, or fractions restricted to the [0, 1] interval. It supports standard one-part models as well as two-part models for modeling boundary values.
}

\usage{
fracreg(y, x, x2 = x, linkbin, linkfrac, type = "1P", inflation = 0, intercept = TRUE,
	table = TRUE, variance = TRUE, var.type = "default", var.eim = TRUE, var.cluster,
	dfc = FALSE,...)
}

\arguments{
  \item{y}{
a numeric vector containing the values of the response variable.
}
  \item{x}{
a numeric matrix, with column names, containing the values of the covariates.
}
  \item{x2}{
a numeric matrix, with column names, containing the values of the covariates in the fractional component of two-part
models if option \code{type = "2P"} is defined. Defaults to \code{x}.
}
  \item{linkbin}{
a description of the link function to use in the binary component of a two-part fractional regression model. 
Available options: \code{logit}, \code{probit}, \code{cauchit}, \code{loglog}, \code{cloglog}.
}
  \item{linkfrac}{
a description of the link function to use in standard fractional regression models or in the fractional component 
of a two-part fractional regression model. Available options: \code{logit}, \code{probit}, \code{cauchit}, 
\code{loglog}, \code{cloglog}.
}
  \item{type}{
a description of the model to estimate: a standard one-part model (\code{1P}, the default), a two-part model 
(\code{2P}), the binary component of a two-part model (\code{2Pbin}), the fractional component of a two-part 
model (\code{2Pfrac}), or a three-part model (\code{3P}) for double boundary inflation.
}
  \item{inflation}{
a numeric value indicating which of the extreme values of \code{0} (the default) or \code{1} is the relevant boundary value for 
defining two-part fractional regression models.
}
  \item{intercept}{
a logical value indicating whether the model should include a constant term or not.
}
  \item{table}{
a logical value indicating whether a summary table with the regression results should be printed.
}
  \item{variance}{
a logical value indicating whether the variance of the estimated parameters should be calculated. Defaults to 
\code{TRUE} whenever \code{table = TRUE}.
}
  \item{var.type}{
a description of the type of variance of the estimated parameters to be calculated. Options are \code{standard} (recommended 
for models estimated by maximum likelihood, such as the binary component of two-part models), \code{robust} (recommended 
for models estimated by quasi-maximum likelihood, such as standard fractional regression models or the fractional 
component of a two-part fractional regression model), \code{cluster} (recommended in the case of panel data) and 
\code{default} (implements the \code{standard} or \code{robust} versions as appropriate).
}
  \item{var.eim}{
a logical value indicating whether the expected information matrix should be used in the calculation of the variance. When 
false, the observation information matrix will be used. Defaults to \code{TRUE}. 
}
  \item{var.cluster}{
a numeric vector containing the values of the variable that specifies to which cluster each observation belongs.
}
  \item{dfc}{
a logical value indicating whether a degrees of freedom correction should be applied to the covariance matrix. Defaults to 
\code{FALSE}.
}
  \item{\dots}{
Arguments to pass to \link[stats]{glm}.
}
}

\details{
\code{fracreg} estimates one- and two-part fractional regression models; see
Ramalho, Ramalho and Murteira (2011) for details on those models. The one-part
models and the fractional component of two-part models are estimated by
Bernoulli-based quasi-maximum likelihood, while the binary component of
two-part models is estimated by maximum likelihood. \code{fracreg} uses the standard
\link[stats]{glm} command to perform the estimations. Therefore,
\code{fracreg} is essentially a convenience command, allowing estimation of several
alternative fractional regression models using the same command. In addition,
\code{fracreg} provides an R-squared measure for all models (calculated as the square
of the correlation coefficient between the actual and fitted values of the
dependent variable), calculates the fitted values of the dependent variable in
two-part models and stores the information needed to implement some very useful
commands for fractional regression models: \link{fracreg.reset} (RESET test),
\link{fracreg.ptest} (P test), \link{fracreg.ggoff} (GGOFF tests) and \link{fracreg.pe} 
(partial effects).
}

\value{
When \code{type = "1P" or "2Pfrac"}, \code{fracreg} returns a list with the following elements:
  \item{class}{"fracreg".
}
  \item{formula}{the model formula.
}
  \item{type}{the name of the estimated model.
}
  \item{link}{the name of the specified link.
}
  \item{method}{estimation method. Currently, "QML" (quasi-maximum likelihood) for fractional components or models and
"ML" (maximum likelihood) for the binary component of two-part models.
}
  \item{p}{a named vector of coefficients.
}
  \item{yhat}{the fitted mean values.
}
  \item{xbhat}{the fitted mean values of the linear predictor.
}
  \item{converged}{logical. Was the algorithm judged to have converged?
}
  \item{x.names}{a vector containing the names of the covariates.
}

If \code{variance = TRUE} or \code{table = TRUE}, the previous list also contains the following elements:
  \item{p.var}{a named covariance matrix.
}
  \item{var.type}{covariance matrix type.
}
  \item{var.eim}{logical. Was the expected information matrix used in the computation of the covariance matrix?
}
  \item{dfc}{logical. Was a degrees of freedom correction used for the computation of the covariance matrix?
}

If \code{var.type = "cluster"}, the list also contains the following element:
  \item{var.cluster}{the variable that specifies to which cluster each observation belongs.
}

When \code{type = "2Pbin"}, \code{fracreg} returns a similar list with the following additional element:
  \item{LL}{the value of the log-likelihood.
}

When \code{type = "2P"}, \code{fracreg} returns the previous lists, indexed by the prefixes \code{resBIN} and 
\code{resFRAC}, and the following additional elements:
  \item{class}{"fracreg".
}
  \item{type}{"2P".
}
  \item{ybase}{a numeric vector containing the values of the response variable.
}
  \item{x2base}{a numeric matrix containing the values of the covariates.
}
  \item{yhat2P}{the overall fitted mean values.
}
  \item{converged}{logical. Were the algorithms judged to have converged in both parts of the model?
}
}

\references{
Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional regression models",
\emph{Journal of Economic Surveys}, 25(1), 19-68.
}

\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}

\seealso{
\code{\link{fracreg.reset}} and \code{\link{fracreg.ggoff}}, for specification tests.\cr
\code{\link{fracreg.ptest}}, for non-nested hypothesis tests.\cr
\code{\link{fracreg.pe}}, for computing partial effects.\cr
\code{fracreghet}, for fitting cross-sectional fractional regression models with unobserved heterogeneity.\cr
\code{fracregpd}, for fitting panel data fractional regression models.
}

\examples{
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

}

