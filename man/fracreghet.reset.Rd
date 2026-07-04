\name{fracreghet.reset}
\alias{fracreghet.reset}

\title{
RESET Test for Fractional Regression Models under Neglected Heterogeneity
}

\description{
\code{fracreghet.reset} is used to test the specification of fractional regression models estimated by GMMx or LINx.
}

\usage{
fracreghet.reset(object, lastpower.vec = 3, version = "Wald", table = T, ...)
}

\arguments{
  \item{object}{
an object containing the results of an \code{fracreghet} command.
}
  \item{lastpower.vec}{
a numeric vector containing the maximum powers of the linear predictors to be used in RESET tests.
}
  \item{version}{
a vector containing the test versions to use. Available options: \code{Wald} (the default) and \code{LM} (only 
available for \code{GMMx}).
}
  \item{table}{
a logical value indicating whether a summary table with the test results should be printed.
}
  \item{\dots}{
Arguments to pass to \link[stats]{nlminb}, which is used to estimate the model under the alternative hypothesis 
when \code{version} is equal to \code{"Wald"} and the null model was estimated by \code{GMMx}.
}
}

\details{
\code{fracreghet.reset} applies the RESET test statistic to fractional regression
models estimated via \code{fracreghet} using the options \code{GMMx} or \code{LINx}. 
\code{fracreghet.reset} may be used to test simultaneously the validity of the link
specification and the transformation applied to the response variable by each 
estimator.  It is taken into account the option that was chosen for computing 
standard errors in the model under evaluation. See Ramalho and Ramalho and 
Ramalho (2016) for details.
}

\value{
\code{fracreghet.reset} returns a named vector with the test results.
}

\references{
Ramalho, E.A. and J.J.S. Ramalho (2016), "Moment-based estimation of nonlinear 
regression models with boundary outcomes and endogeneity, with applications to 
nonnegative and fractional responses", \emph{Econometric Reviews}, forthcoming 
(DOI: 10.1080/07474938.2014.976531).
}

\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}

\seealso{
\code{\link{fracreghet}}, for fitting fractional regression models under unobserved heterogeneity.\cr
\code{\link{fracreghet.pe}}, for computing partial effects.
}

\examples{
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

}
