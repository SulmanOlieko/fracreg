\name{fracregpd}
\alias{fracregpd}

\title{
Fitting Panel Data Fractional Response Models
}

\description{
\code{fracregpd} is used to fit panel data regression models when the dependent variable has a bounded, fractional nature.
}
\usage{fracregpd(id, time, y, x, z, var.endog, x.exogenous = T, lags, start, type,
	GMMww.cor = T, link = "logit", intercept = T, table = T, variance = T,
	var.type = "cluster", tdummies = F, bootstrap = F, B = 200, ...)
}

\arguments{
  \item{id}{a numeric vector identifying the cross-sectional units.
}
  \item{time}{a numeric vector identifying the time periods in which the cross-sectional units were observed.
}
  \item{y}{a numeric vector containing the values of the response variable.
}
  \item{x}{a numeric matrix, with column names, containing the values of all covariates (exogenous and endogenous).
}
  \item{z}{a numeric matrix, with column names, containing the values of all exogenous variables (covariates and external
instrumental variables). Only required in case of endogenous explanatory variables.
}
  \item{var.endog}{a numeric vector containing the values of the endogenous covariate (or of some transformation of it), which will be used as dependent variable in the linear reduced form assumed for application of the \code{QMLcre} estimator. Only required for this estimator.
}
  \item{x.exogenous}{a logical value indicating whether all explanatory variables are assumed to be exogenous or not.
}
  \item{lags}{a logical value indicating whether the first lags of \code{x} or  \code{z} should be used as instruments for \code{x}. Defaults to \code{TRUE} for the GMMww and GMMc estimators and to \code{FALSE} for the remaining estimators. The \code{GMMcre} and \code{QMLcre} estimators do not admit lagged instruments.
}
  \item{start}{a numeric vector containing the initial values for the parameters to be optimized. Optional.
}
  \item{type}{a description of the estimator to compute: \code{GMMww}, \code{GMMc}, \code{GMMbgw}, \code{GMMpfe}, \code{GMMcre}, \code{GMMpre} or \code{QMLcre}.
}
  \item{GMMww.cor}{a logical value indicating whether each explanatory variable should be transformed in deviations from its overall mean before computing the \code{GMMww} estimator.
}
  \item{link}{a description of the link function to use. Available options for all GMM estimators: \code{logit} and \code{cloglog}. Only option for the \code{QMLcre} estimator: \code{probit}.
}
  \item{intercept}{a logical value indicating whether the model should include a constant term or not. Only relevant for the 
\code{GMMpre} estimator.
}
  \item{table}{a logical value indicating whether a summary table with the regression results should be printed.
}
  \item{variance}{a logical value indicating whether the variance of the estimated parameters should be calculated. Defaults to 
\code{TRUE} whenever \code{table = TRUE}.
}
  \item{var.type}{a description of the type of variance of the estimated parameters to be calculated. Options are \code{cluster}, the default, and \code{robust}. In overidentified models, it also affects the parameter estimates via the GMM weighting matrix.
}
  \item{tdummies}{a logical value indicating whether time dummies should be included among the model explanatory variables.
}
  \item{bootstrap}{a logical value indicating whether bootstrap should be used in the estimation of the parameter standard errors.
}
  \item{B}{the number of bootstrap replications.
}
  \item{\dots}{Arguments to pass to \link[stats]{nlminb}.
}
}

\details{
\code{fracregpd} computes the GMM estimators proposed in Ramalho, Ramalho and Coelho (2018) for panel data fractional response models with both time-variant and time-invariant unobserved heterogeneity and endogeneous covariates: GMMww, GMMc, GMMbgw, GMMpfe, GMMcre and GMMpre. In addition, \code{fracregpd} also computes QMLcre, which was proposed by Papke and Wooldridge (2008) and Wooldridge (2019). 

\strong{Correlated Random Effects (CRE) - QMLcre:}
In panel data, unobserved individual-specific heterogeneity \eqn{c_i} may be correlated with the covariates \eqn{x_{it}}. The CRE approach (Papke and Wooldridge, 2008) models this dependence by projecting \eqn{c_i} onto the time averages of the strictly exogenous covariates \eqn{\bar{x}_i}:
\deqn{c_i = \psi + \bar{x}_i \xi + a_i}
where \eqn{a_i} is an error term independent of \eqn{x_i}. Assuming \eqn{a_i | x_i \sim N(0, \sigma_a^2)} and a probit link, integrating out \eqn{a_i} yields the "population-averaged" or scaled conditional mean:
\deqn{E(y_{it} | x_i) = G(x_{it} \beta_a + \psi_a + \bar{x}_i \xi_a)}
where the parameters with subscript \eqn{a} are scaled by \eqn{(1 + \sigma_a^2)^{-1/2}}. This equation is estimated via pooled Bernoulli QML.

\strong{Generalized Method of Moments (GMM):}
For models where strict exogeneity fails or the link function is an exponential-type link, Ramalho et al. (2018) propose GMM estimators based on the following general moment conditions:
\deqn{E[Z_{it} (H(y_{it}) - \exp(x_{it}\beta + c_i))] = 0}
where \eqn{H(\cdot)} is a transformation function and \eqn{Z_{it}} is a matrix of valid instruments. Estimators such as GMMww, GMMc, and GMMbgw use different transformations to eliminate the unobserved fixed effect \eqn{c_i} before applying GMM.

For overidentified models, \code{fracregpd} calculates Hansen's J statistic to test the validity of the overidentifying restrictions.
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

Ramalho, E. A., Ramalho, J. J. S., & Coelho, L. M. S. (2018), "Exponential Regression of Fractional-Response Fixed-Effects Models with an Application to Firm Capital Structure", \emph{Journal of Econometric Methods}, 7(1), 20150019.

Wooldridge, J. M. (2019). Correlated random effects models with unbalanced panels. \emph{Journal of Econometrics}, 211(1), 137-150.
}

\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}

\seealso{
\code{\link{fracreg}}, for fitting standard cross-sectional fractional response models.\cr
\code{\link{fracreghet}}, for fitting cross-sectional fractional response models with unobserved heterogeneity.
}

\examples{
set.seed(123)
# Simulating Panel Data
N <- 100
T_periods <- 5
id <- rep(1:N, each = T_periods)
time <- rep(1:T_periods, times = N)
x_panel <- rnorm(N * T_periods)

# Unobserved individual effect (CRE)
c_i <- rep(rnorm(N), each = T_periods) 
y_panel <- exp(x_panel + c_i) / (1 + exp(x_panel + c_i))

X <- cbind(x_panel = x_panel)

# Endogenous variable and instrument simulation
z_panel <- rnorm(N * T_periods)
u_panel <- 0.5 * z_panel + rnorm(N * T_periods)
var_endog <- 0.8 * z_panel + u_panel
y_endog <- exp(x_panel + 1.2 * var_endog + c_i + u_panel) / 
             (1 + exp(x_panel + 1.2 * var_endog + c_i + u_panel))

X_endog <- cbind(x_panel = x_panel, var_endog = var_endog)
Z_inst <- cbind(x_panel = x_panel, z_panel = z_panel)

# Estimate a Correlated Random Effects (CRE) Model
fracregpd(id=id, time=time, y=y_panel, x=X, type="QMLcre", link="probit")

# Exogeneity, no lags, no time dummies, clustered standard errors, GMMbgw estimator
fracregpd(id=id, time=time, y=y_panel, x=X, type="GMMbgw")

# Lagged covariates and instruments, robust standard errors, GMMww estimator
fracregpd(id=id, time=time, y=y_panel, x=X, lags=TRUE, type="GMMww", var.type="robust")

# Endogeneity, time dummies, GMMpfe estimator
fracregpd(id=id, time=time, y=y_endog, x=X_endog, z=Z_inst, 
          x.exogenous=FALSE, type="GMMpfe", tdummies=TRUE)
}
