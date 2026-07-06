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
           var.type = "robust", var.cluster, adjust = 0, offset = NULL, ...)
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
  \item{start}{a numeric vector containing the initial values for the parameters to be optimised. Optional.
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
  \item{offset}{an optional numeric vector containing an offset. It must be of the same dimension as the response variable. It specifies that the variable should be included in the model with its coefficient constrained to 1.
}
  \item{\dots}{Arguments to pass to \link[stats]{nlminb}.
}
}

\details{
\code{fracreghet} computes the GMM estimators proposed in Ramalho and Ramalho (2017) for fractional response models with unobserved heterogeneity: GMMx, which allows for neglected heterogeneity but not for endogeneity; GMMxv, which allows both issues and assumes a linear reduced form for the endogeneous covariate (or for a transformation of it); and GMMz, which also allows for both issues but does not require the assumption of a reduced form for the endogenous covariate. In addition, \code{fracreghet} also computes three linearised estimators (LINx, LINxv and LINz) that have similar features to their GMM counterparts. It also provides a QML estimator (QMLxv) that addresses endogeneity using a Control Function (CF) approach, which includes the first-stage reduced-form residuals as an additional regressor in the main fractional equation, providing a Hausman-type test for endogeneity. 

\strong{Control Function (CF) Approach - QMLxv:}
When a continuous regressor \eqn{y_{2i}} is endogenous, the CF approach (Papke and Wooldridge, 2008; Terza et al., 2008) uses a two-stage procedure. First, a linear reduced form is estimated:
\deqn{y_{2i} = z_i \pi + v_i}
where \eqn{z_i} includes all exogenous variables and external instruments. The residuals \eqn{\hat{v}_i} are then included in the fractional response model:
\deqn{E(y_{1i} | z_i, y_{2i}, v_i) = G(x_i \beta + \gamma \hat{v}_i)}
A test of \eqn{H_0: \gamma = 0} serves as a robust Hausman-type test for endogeneity.

\strong{Generalised Method of Moments (GMM):}
For estimators like GMMz, which do not strictly require a linear reduced form, the estimation relies on population orthogonality conditions between the instruments \eqn{Z_i} and the model residuals:
\deqn{E[Z_{i} (y_i - G(x_i \beta))] = 0}
or via specific transformations of the dependent variable to eliminate unobserved heterogeneity (Ramalho and Ramalho, 2017).

For overidentified models, \code{fracreghet} calculates Hansen's J statistic. For \code{GMMx} and \code{LINx}, \code{fracreghet} stores the information needed to implement the RESET test (\link{fracreghet.reset}). For all estimators, \code{fracreghet} stores the information needed to calculate partial effects (\link{fracreghet.pe}).
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
Papke, L. E. and Wooldridge, J. M. (2008), "Panel data methods for fractional response variables with an application to test pass rates", \emph{Journal of Econometrics}, 145, 121-133.

Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of nonlinear regression models with boundary outcomes and endogeneity, with applications to nonnegative and fractional responses", \emph{Econometric Reviews}, 36(4), 397-420.

Terza, J. V., Basu, A., and Rathouz, P. J. (2008), "Two-stage residual inclusion estimation: addressing endogeneity in health econometric modeling", \emph{Journal of Health Economics}, 27(3), 531-543.
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
### Empirical 401(k) Examples 
data("fracreg_k401k") 
y <- fracreg_k401k$prate 
X_het <- cbind(mrate = fracreg_k401k$mrate, ltotemp = fracreg_k401k$ltotemp)
 
# fracreghet estimators do not allow exact 1s or 0s
y_adj <- y
y_adj[y_adj == 1] <- 0.999

# Instrument mrate using age

Z_emp <- cbind(age = fracreg_k401k$age, ltotemp = fracreg_k401k$ltotemp) 
fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], type="QMLxv", link="logit") 
 
### Simulated Examples 

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
