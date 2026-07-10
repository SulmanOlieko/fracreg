\name{fracreg.pe}
\alias{fracreg.pe}

\title{
Fractional Response Regressions - Partial Effects
}

\description{
\code{fracreg.pe} is used to compute average and/or conditional partial effects in fractional response models.
}

\usage{
fracreg.pe(object, APE = TRUE, CPE = FALSE, at = NULL, which.x = NULL, 
           variance = TRUE, table = FALSE)
}

\arguments{
  \item{object}{an object containing the results of an \code{fracreg} command.
}
  \item{APE}{a logical value indicating whether average partial effects are to be computed.
}
  \item{CPE}{a logical value indicating whether conditional partial effects are to be computed.
}
  \item{at}{a numeric vector containing the covariates' values at which the conditional partial effects are to be computed or 
the strings \code{"mean"} (the default) or \code{"median"}, in which cases the covariates are evaluated at their 
mean or median values (or mode, in case of dummy variables), respectively.
}
  \item{which.x}{a vector containing the names of the covariates to which the partial effects are to be computed.
}
  \item{variance}{a logical value indicating whether the variance of the estimated partial effects should be calculated. Defaults to 
\code{TRUE} whenever \code{table = TRUE}.
}
  \item{table}{a logical value indicating whether a summary table with the results should be printed.
}
}

\details{
\code{fracreg.pe} calculates partial effects for fractional response models estimated via \code{fracreg}. \code{fracreg.pe} may be used to compute average or conditional partial effects for: (i) one-part fractional response models; (ii) the binary components of two-part and three-part fractional response models; (iii) the fractional components of two-part and three-part fractional response models; and (iv) two-part and three-part fractional response models overall. 

\strong{Partial Effects for Continuous Variables:}
For a continuous covariate \eqn{x_k}, the partial effect on the conditional mean \eqn{E(y|x) = G(x\beta)} is the first derivative with respect to \eqn{x_k}:
\deqn{PE_k(x) = \frac{\partial E(y|x)}{\partial x_k} = g(x\beta)\beta_k}
where \eqn{g(\cdot)} is the probability density function corresponding to the link function \eqn{G(\cdot)}.

\strong{Partial Effects for Discrete Variables:}
For a discrete or dummy covariate \eqn{x_k}, the partial effect is calculated as the discrete difference in the expected value when \eqn{x_k} changes from 0 to 1, holding all other variables \eqn{x_{-k}} constant:
\deqn{PE_k(x) = G(x_{-k}\beta_{-k} + \beta_k) - G(x_{-k}\beta_{-k})}

\strong{Average vs. Conditional Partial Effects:}
- \strong{Average Partial Effects (APE):} Evaluated for each observation \eqn{i} in the sample and then averaged:
\deqn{APE_k = \frac{1}{N} \sum_{i=1}^N PE_k(x_i)}
- \strong{Conditional Partial Effects (CPE):} Evaluated at a specific vector of covariate values \eqn{x^*} (e.g., the sample mean or median):
\deqn{CPE_k = PE_k(x^*)}

For calculating standard errors, it is taken into account the option that was previously chosen for estimating the model. See Ramalho, Ramalho and Murteira (2011) and Fang and Ma (2013) for details on the computation of partial effects in the fractional response framework.
}

\value{
\code{fracreg.pe} returns a list with the following element:
  \item{PE.p}{a named vector of partial effects.
}

If \code{variance = TRUE} or \code{table = TRUE}, the previous list also contains the following element:
  \item{PE.sd}{a named vector of standard errors of the estimated partial effects.
}

When both average and conditional partial effects are requested, two lists containing the previous elements are returned, indexed by the prefixes \code{ape} and \code{cpe}.
}

\references{
Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional response models",
\emph{Journal of Economic Surveys}, 25(1), 19-68.

Fang, K., & Ma, S. (2013), "Three-part model for fractional response variables with application to Chinese household health insurance coverage", \emph{Journal of Applied Statistics}, 40(5), 925-940.
}

\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}

\seealso{
\code{\link{fracreg}}, for fitting fractional response models.\cr
\code{\link{fracreg.reset}} and \code{\link{fracreg.ggoff}}, for specification tests.\cr
\code{\link{fracreg.ptest}}, for non-nested hypothesis tests.
}

\examples{
### Empirical 401(k) Examples
data("fracreg_k401k")
y <- fracreg_k401k$prate
X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
           totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)

m <- fracreg(y, X, type="1P", linkfrac="logit")
pe_res <- fracreg.pe(m)
summary(pe_res)

### Simulated Examples

N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Computing average partial effects for a logit fractional response model
mod <- fracreg(y,X,linkfrac="logit")
pe_res <- fracreg.pe(mod)
summary(pe_res)

#Computing average partial effects for a binary logit + fractional probit
#two-part model
mod <- fracreg(y,X,linkbin="logit",linkfrac="probit",type="2P",inf=1)
pe_res <- fracreg.pe(mod)
summary(pe_res)

#Computing conditional partial effects for X2 in the logit component
#of a two-part fractional response model, with the covariates evaluated
#at their median values
mod <- fracreg(y,X,linkfrac="logit",type="2Pfrac",inf=1)
pe_res <- fracreg.pe(mod,APE=FALSE,CPE=TRUE,at="median",which.x="X2")
summary(pe_res)

#Computing average partial effects for a three-part double-inflated model
y3p <- y
y3p[1:20] <- 0
y3p[21:40] <- 1
res3p <- fracreg(y3p,X,linkbin=c("logit","probit"),linkfrac="logit",type="3P")
pe_res <- fracreg.pe(res3p)
summary(pe_res)
}
