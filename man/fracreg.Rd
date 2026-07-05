\name{fracreg}
\alias{fracreg}

\title{
Fitting Fractional Response Models
}

\description{
\code{fracreg} is used to fit fractional response models, which are appropriate for responses that are proportions, percentages, or fractions restricted to the [0, 1] interval. It supports standard one-part models, two-part hurdle models for modelling boundary values at 0 or 1, and three-part models for double inflation at both 0 and 1.
}

\usage{
fracreg(y, x, x2 = x, linkbin, linkfrac, type = "1P", inflation = 0, 
        intercept = TRUE, table = TRUE, variance = TRUE, var.type = "default", 
        var.eim = TRUE, var.cluster, dfc = FALSE, ...)
}

\arguments{
  \item{y}{a numeric vector containing the values of the response variable.
}
  \item{x}{a numeric matrix, with column names, containing the values of the covariates.
}
  \item{x2}{a numeric matrix, with column names, containing the values of the covariates in the fractional component of two-part
models if option \code{type = "2P"} is defined. Defaults to \code{x}.
}
  \item{linkbin}{a description of the link function to use in the binary component of a two-part fractional response model, or a vector of two link functions for the two binary components of a three-part model (e.g. \code{c("logit", "probit")}). Available options: \code{logit}, \code{probit}, \code{cauchit}, \code{loglog}, \code{cloglog}.
}
  \item{linkfrac}{a description of the link function to use in standard fractional response models or in the fractional component of a two-part fractional response model. Available options: \code{logit}, \code{probit}, \code{cauchit}, \code{loglog}, \code{cloglog}.
}
  \item{type}{a description of the model to estimate: a standard one-part model (\code{1P}, the default), a two-part model (\code{2P}), the binary component of a two-part model (\code{2Pbin}), the fractional component of a two-part model (\code{2Pfrac}), or a three-part model (\code{3P}) for double boundary inflation.
}
  \item{inflation}{a numeric value indicating which of the extreme values of \code{0} (the default) or \code{1} is the relevant boundary value for defining two-part fractional response models.
}
  \item{intercept}{a logical value indicating whether the model should include a constant term or not.
}
  \item{table}{a logical value indicating whether a summary table with the regression results should be printed.
}
  \item{variance}{a logical value indicating whether the variance of the estimated parameters should be calculated. Defaults to \code{TRUE} whenever \code{table = TRUE}.
}
  \item{var.type}{a description of the type of variance of the estimated parameters to be calculated. Options are \code{standard} (recommended for models estimated by maximum likelihood, such as the binary component of two-part models), \code{robust} (recommended for models estimated by quasi-maximum likelihood, such as standard fractional response models or the fractional component of a two-part fractional response model), \code{cluster} (recommended in the case of panel data) and \code{default} (implements the \code{standard} or \code{robust} versions as appropriate).
}
  \item{var.eim}{a logical value indicating whether the expected information matrix should be used in the calculation of the variance. When false, the observation information matrix will be used. Defaults to \code{TRUE}. 
}
  \item{var.cluster}{a numeric vector containing the values of the variable that specifies to which cluster each observation belongs.
}
  \item{dfc}{a logical value indicating whether a degrees of freedom correction should be applied to the covariance matrix. Defaults to \code{FALSE}.
}
  \item{\dots}{Arguments to pass to \link[stats]{glm}.
}
}

\details{
\code{fracreg} estimates one-part, two-part hurdle, and three-part double-inflated fractional response models; see Ramalho, Ramalho and Murteira (2011) and Fang and Ma (2013) for details on those models. 

\strong{One-Part Fractional Response Models (\code{type = "1P"}):}
The standard one-part model assumes that the conditional expectation of the fractional response \eqn{y_i \in [0,1]} is given by:
\deqn{E(y_i|x_i) = G(x_i \beta)}
where \eqn{G(\cdot)} is a known non-linear link function mapping the linear predictor to the unit interval (e.g., logit, probit). The parameters \eqn{\beta} are estimated by maximising the Bernoulli-based quasi-log-likelihood function:
\deqn{\ln L_i(\beta) = y_i \ln[G(x_i \beta)] + (1 - y_i) \ln[1 - G(x_i \beta)]}
This estimator requires only the correct specification of the conditional mean to yield consistent parameter estimates (Papke and Wooldridge, 1996).

\strong{Two-Part Hurdle Models (\code{type = "2P"}):}
When the data exhibits a boundary mass (e.g., at \eqn{y_i = 0}), the two-part hurdle model handles the boundary values separately from the interior fractional values. Let \eqn{y_i^*} be a binary indicator such that \eqn{y_i^* = 1} if \eqn{y_i > 0} and \eqn{y_i^* = 0} otherwise. The probability of observing a boundary value is modelled as:
\deqn{P(y_i = 0 | x_{1i}) = 1 - F(x_{1i} \gamma_1)}
\deqn{P(y_i > 0 | x_{1i}) = F(x_{1i} \gamma_1)}
where \eqn{F(\cdot)} is a binary link function. Conditional on observing an interior fractional value, the response is modelled as:
\deqn{E(y_i | x_{2i}, y_i > 0) = G(x_{2i} \beta_2)}
The unconditional mean of the response is therefore:
\deqn{E(y_i|x_i) = F(x_{1i} \gamma_1) \times G(x_{2i} \beta_2)}

\strong{Three-Part Double Inflated Models (\code{type = "3P"}):}
For data containing boundary mass at both \eqn{0} and \eqn{1}, the three-part model estimates two separate binary mechanisms for each boundary and a fractional component for the interior values \eqn{(0, 1)}, extending the two-part logic to double inflation (Fang and Ma, 2013).

\code{fracreg} uses the standard \link[stats]{glm} command to perform the estimations. Therefore, \code{fracreg} is essentially a convenience command, allowing estimation of several alternative fractional response models using the same command. In addition, \code{fracreg} provides an R-squared measure for all models (calculated as the square of the correlation coefficient between the actual and fitted values of the dependent variable), calculates the fitted values of the dependent variable in two-part models and stores the information needed to implement some very useful commands for fractional response models: \link{fracreg.reset} (RESET test), \link{fracreg.ptest} (P test), \link{fracreg.ggoff} (GGOFF tests) and \link{fracreg.pe} (partial effects).
}

\value{When \code{type = "1P" or "2Pfrac"}, \code{fracreg} returns a list with the following elements:
  \item{class}{"fracreg".
}
  \item{formula}{the model formula.
}
  \item{type}{the name of the estimated model.
}
  \item{link}{the name of the specified link.
}
  \item{method}{estimation method. Currently, "QML" (quasi-maximum likelihood) for fractional components or models and"ML" (maximum likelihood) for the binary component of two-part models.
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

When \code{type = "2P"}, \code{fracreg} returns the previous lists, indexed by the prefixes \code{resBIN} and \code{resFRAC}, and the following additional elements:
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

When \code{type = "3P"}, \code{fracreg} returns the previous lists, indexed by the prefixes \code{resBIN0}, \code{resBIN1}, and \code{resFRAC}, and the following additional elements:
  \item{class}{"fracreg".
}
  \item{type}{"3P".
}
  \item{ybase}{a numeric vector containing the values of the response variable.
}
  \item{x2base}{a numeric matrix containing the values of the covariates.
}
  \item{yhat3P}{the overall fitted mean values.
}
  \item{converged}{logical. Were the algorithms judged to have converged in all parts of the model?
}
}

\references{
Papke, L. E. and Wooldridge, J. M. (1996), "Econometric methods for fractional response variables with an application to 401(k) plan participation rates", \emph{Journal of Applied Econometrics}, 11(6), 619-632.

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional response models",
\emph{Journal of Economic Surveys}, 25(1), 19-68.

Fang, K., & Ma, S. (2013), "Three-part model for fractional response variables with application to Chinese household health insurance coverage", \emph{Journal of Applied Statistics}, 40(5), 925-940.
}

\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}

\seealso{
\code{\link{fracreg.reset}} and \code{\link{fracreg.ggoff}}, for specification tests.\cr
\code{\link{fracreg.ptest}}, for non-nested hypothesis tests.\cr
\code{\link{fracreg.pe}}, for computing partial effects.\cr
\code{fracreghet}, for fitting cross-sectional fractional response models with unobserved heterogeneity.\cr
\code{fracregpd}, for fitting panel data fractional response models.
}

\examples{
### Empirical 401(k) Examples
data("fracreg_k401k")
y <- fracreg_k401k$prate
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)

# 1P Model
fracreg(y, X, type="1P", linkfrac="logit")

# 2P Model (modelling mass at 1)
fracreg(y, X, type="2P", inflation=1, linkbin="logit", linkfrac="logit")

# 3P Model (inject artificial 0s for demonstration)
y_3p <- y; y_3p[1:50] <- 0
fracreg(y_3p, X, type="3P", linkbin=c("logit","logit"), linkfrac="logit")

### Simulated Examples

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

X <- cbind(x1 = x1, x2 = x2)

# fracreg estimation of a logit fractional response model
fracreg(y, X, type="1P", linkfrac="logit")

# fracreg estimation of the binary logit component of the two-part fractional
# regression model with y=0 as the relevant boundary value
fracreg(y, X, type="2Pbin", inflation=0, linkbin="logit")

# fracreg estimation of the fractional component of the two-part fractional
# regression model with y=0 as the relevant boundary value and using a
# probit link function
fracreg(y, X, type="2Pfrac", inflation=0, linkfrac="probit")

# fracreg estimation of both components of a two-part fractional response model
# with y=0 as the relevant boundary value and using a cloglog binary link
# function and a logit fractional link function
fracreg(y, X, type="2P", inflation=0, linkbin="cloglog", linkfrac="logit")

# Three-part double-inflated model (y has both 0s and 1s)
fracreg(y, X, type="3P", linkbin=c("logit","probit"), linkfrac="logit")
}

