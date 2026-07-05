\name{fracreghet}
\alias{fracreghet}

\title{
Fitting Fractional Response Models under Unobserved Heterogeneity
}

\description{
\code{fracreghet} is used to fit fractional response models under unobserved heterogeneity, i.e. regression models for proportions, percentages or fractions that suffer from neglected heterogeneity and/or endogeneity issues.
}
\usage{
fracreghet(y, x, z = x, var.endog, start, type = "GMMx", link = "logit", 
           intercept = T, table = T, variance = T, 
           var.type = "robust", var.cluster, adjust = 0, ...)
}

\arguments{
  \item{y}{a numeric vector containing the values of the response variable.
}
  \item{x}{a numeric matrix, with column names, containing the values of all covariates (exogenous and endogenous).
}
  \item{z}{a numeric matrix, with column names, containing the values of all exogenous variables (covariates and instrumental 
variables). Defaults to \code{x}.
}
  \item{var.endog}{a numeric vector containing the values of the endogenous covariate (or of some transformation of it), which will be used as dependent variable in the linear reduced form assumed for application of xv-type estimators.
}
  \item{start}{a numeric vector containing the initial values for the parameters to be optimized. Optional.
}
  \item{type}{a description of the estimator to compute: \code{GMMx} (the default), \code{GMMxv}, \code{GMMz}, \code{LINx}, \code{LINxv}, \code{LINz} or \code{QMLxv}.
}
  \item{link}{a description of the link function to use. Available options for all estimators: \code{logit} and \code{cloglog}. Additional available options for QML and LIN estimators: \code{probit}, \code{cauchit} and \code{loglog}.
}
  \item{intercept}{a logical value indicating whether the model should include a constant term or not.
}
  \item{table}{a logical value indicating whether a summary table with the regression results should be printed.
}
  \item{variance}{a logical value indicating whether the variance of the estimated parameters should be calculated. Defaults to \code{TRUE} whenever \code{table = TRUE}.
}
  \item{var.type}{a description of the type of variance of the estimated parameters to be calculated. Options are \code{robust}, the default, and \code{cluster}.
}
  \item{var.cluster}{a numeric vector containing the values of the variable that specifies to which cluster each observation belongs.
}
  \item{adjust}{the numeric value to be added to the response variable in case of boundary observations when the LIN estimators are applied or the string \code{drop}, which implies that the boundary observations are dropped.
}
  \item{\dots}{Arguments to pass to \link[stats]{nlminb}.
}
}

\details{
\code{fracreghet} computes the GMM estimators proposed in Ramalho and Ramalho (2017) for fractional response models with unobserved heterogeneity: GMMx, which allows for neglected heterogeneity but not for endogeneity; GMMxv, which allows both issues and assumes a linear reduced form for the endogeneous covariate (or for a transformation of it); and GMMz, which also allows for both issues but does not require the assumption of a reduced form for the endogenous covariate. In addition, \code{fracreghet} also computes three linearized estimators (LINx, LINxv and LINz) that have similar features to their GMM counterparts. It also provides a QML estimator (QMLxv) that addresses endogeneity using a Control Function (CF) approach, which includes the first-stage reduced-form residuals as an additional regressor in the main fractional equation, providing a Hausman-type test for endogeneity. See Ramalho and Ramalho (2017) for details on each estimator. For overidentified models, \code{fracreghet} calculates Hansen's J statistic. For \code{GMMx} and \code{LINx}, \code{fracreghet} stores the information needed to implement the RESET test (\link{fracreghet.reset}). For all estimators, \code{fracreghet} stores the information needed to calculate partial effects (\link{fracreghet.pe}).
}

\value{\code{fracreghet} returns a list with the following elements:
  \item{class}{"fracreghet".
}
  \item{formula}{the model formula.
}
  \item{type}{the name of the estimator computed.
}
  \item{link}{the name of the specified link.
}
  \item{adjust}{The value or the type of the adjustment applied to LIN estimators.
}
  \item{p}{a named vector of coefficients.
}
  \item{Hy}{the transformed values of the response variable when GMM or LIN estimators are computed or the 
   values of the response variable in the QML case.
}
  \item{xbhat}{the fitted mean values of the linear predictor (for xv-type estimators, includes the term relative to the first-stage residual).
}
  \item{converged}{logical. Was the algorithm judged to have converged?
}
  \item{x.names}{a vector containing the names of the covariates.
}

In case of an overidentifying model, the following element is also returned:
  \item{J}{the result of Hansen's J test of overidentifying moment conditions.
}

If \code{variance = TRUE} or \code{table = TRUE} and the algorithm converged successfully, the previous list also contains the following elements:

  \item{p.var}{a named covariance matrix.
}
  \item{var.type}{covariance matrix type.
}

If \code{var.type = "cluster"}, the list also contains the following element:
  \item{var.cluster}{the variable that specifies to which cluster each observation belongs.
}
}

\references{
Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of nonlinear regression models with boundary outcomes and endogeneity, with applications to nonnegative and fractional responses", \emph{Econometric Reviews}, 36(4), 397-420.
}

\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}

\seealso{
\code{\link{fracreghet.reset}}, for the RESET test.\cr
\code{\link{fracreghet.pe}}, for computing partial effects.\cr
\code{\link{fracreg}}, for fitting standard cross-sectional fractional response models.\cr
\code{\link{fracregpd}}, for fitting panel data fractional response models.
}

\examples{
set.seed(123)
N <- 1000
x1 <- rnorm(N)

# Simulating an endogenous variable (var.endog) and an instrument (z1)
z1 <- rnorm(N)
u <- 0.5 * z1 + rnorm(N)
var.endog <- 0.8 * z1 + u
y_endog <- exp(0.5 * x1 + 1.2 * var.endog + u) / (1 + exp(0.5 * x1 + 1.2 * var.endog + u))

# Avoid exact 0 or 1 boundaries for some estimators
y_endog[y_endog <= 0] <- 0.01
y_endog[y_endog >= 1] <- 0.99

X <- cbind(x1 = x1, var.endog = var.endog)
Z <- cbind(x1 = x1, z1 = z1)

# Exogeneity (assuming var.endog is exogenous for comparison), GMMx estimator
fracreghet(y = y_endog, x = X, type = "GMMx", link = "logit")

# Endogeneity, GMMz estimator (does not require reduced form for endog)
fracreghet(y = y_endog, x = X, z = Z, type = "GMMz", link = "logit")

# Endogeneity, GMMxv estimator (assumes linear reduced form for var.endog)
fracreghet(y = y_endog, x = X, z = Z, var.endog = var.endog, type = "GMMxv", link = "logit")

# Endogeneity, QMLxv control function approach
fracreghet(y = y_endog, x = X, z = Z, var.endog = var.endog, type = "QMLxv", link = "logit")
}
