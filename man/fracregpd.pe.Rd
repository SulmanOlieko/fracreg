\name{fracregpd.pe}
\alias{fracregpd.pe}
\title{Partial Effects for Fractional Panel Data Regression}
\usage{
fracregpd.pe(object, APE = TRUE, CPE = FALSE, at = NULL, which.x = NULL, 
             variance = TRUE, table = FALSE, ...)
}
\arguments{
\item{object}{An object of class \code{fracregpd}.}

\item{APE}{logical. Compute Average Partial Effects?}

\item{CPE}{logical. Compute Conditional Partial Effects? (Not currently supported for panel data models).}

\item{at}{numeric vector. The values at which to evaluate the CPE.}

\item{which.x}{character vector. Variables for which to compute partial effects. By default, auxiliary CRE parameters (like \code{_mean} and \code{vhat}) are automatically excluded so that partial effects are only computed for the main structural parameters.}

\item{variance}{logical. Compute standard errors using the Delta method?}

\item{table}{logical. Print the resulting partial effects table?}

\item{...}{further arguments passed to or from other methods.}
}
\value{
An object of class \code{fracreg.pe} containing the standard coefficient tables with Average Partial Effects.
}
\description{
Computes Average Partial Effects (APEs) for fractional panel data models fit using \code{fracregpd}. In correlated random effects (CRE) models, the time-averages of the covariates are used merely to control for unobserved heterogeneity, and as such, they are automatically filtered out when computing the structural partial effects. Standard errors are computed using the Delta method and the bootstrapped parameter variance matrix.
}
\seealso{
\code{\link{fracregpd}}, \code{\link{fracreg.pe}}
}
\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}
