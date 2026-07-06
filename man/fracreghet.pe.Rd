\name{fracreghet.pe}
\alias{fracreghet.pe}

\title{
Fractional Response Models under Unobserved Heterogeneity - Partial Effects
}
\description{
\code{fracreghet.pe} is used to compute average and/or conditional partial effects in fractional response models 
under unobserved heterogeneity.
}
\usage{
fracreghet.pe(object, smearing = T, APE = T, CPE = F, at = NULL, 
              which.x = NULL, table = T, variance = T)
}

\arguments{
  \item{object}{an object containing the results of an \code{fracreghet} command.
}
  \item{smearing}{a logical value indicating whether the smearing correction is to be applied
}
  \item{APE}{a logical value indicating whether average partial effects are to be computed.
}
  \item{CPE}{a logical value indicating whether conditional partial effects are to be computed.
}
  \item{at}{a numeric vector containing the covariates' values at which the conditional partial effects are to be computed or the strings \code{"mean"} (the default) or \code{"median"}, in which cases the covariates are evaluated at their mean or median values (or mode, in case of dummy variables), respectively.
}
  \item{which.x}{a vector containing the names of the covariates to which the partial effects are to be computed.
}
  \item{table}{a logical value indicating whether a summary table with the results should be printed.
}
  \item{variance}{a logical value indicating whether the variance of the estimated partial effects should be calculated. Defaults to \code{TRUE} whenever \code{table = TRUE}.
}
}

\details{
\code{fracreghet.pe} calculates partial effects for fractional response models estimated via \code{fracreghet}. \code{fracreghet.pe} may be used to compute average or conditional partial effects. These partial effects may be conditional only on observables, using the smearing estimator, or also on unobservables, setting the error term to zero.

\strong{Partial Effects under Unobserved Heterogeneity:}
When unobserved heterogeneity or endogeneity is present, calculating partial effects requires dealing with the unobserved error \eqn{v_i}. Let the conditional mean be \eqn{E(y|x, v) = G(x\beta + \gamma v)}.
- \strong{Conditional on Observables (Smearing):} The unobserved heterogeneity is integrated out over its empirical distribution. The average partial effect for a continuous variable \eqn{x_k} is computed as:
\deqn{PE_k(x) = \frac{1}{N} \sum_{i=1}^N g(x\beta + \gamma \hat{v}_i) \beta_k}
- \strong{Conditional on Unobservables (Error = 0):} The partial effect is evaluated for an individual with the mean level of unobserved heterogeneity (\eqn{v = 0}):
\deqn{PE_k(x) = g(x\beta) \beta_k}

For discrete variables, the partial effects are calculated as the discrete differences evaluated using either the smearing approach or setting the error term to zero. 

For calculating standard errors, it is taken into account the option that was previously chosen for estimating the model. See Ramalho and Ramalho (2017) for details on the computation of partial effects for fractional response models under unobserved heterogeneity.
}

\value{
\code{fracreghet.pe} returns a list with the following element:
  \item{PE.p}{a named vector of partial effects.
}

If \code{variance = TRUE} or \code{table = TRUE}, the previous list also contains the following element:
  \item{PE.sd}{a named vector of standard errors of the estimated partial effects.
}

When both average and conditional partial effects are requested, two lists containing the previous elements are returned, indexed by the prefixes \code{ape} and \code{cpe}.
}

\references{
Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of nonlinear regression models with boundary outcomes and endogeneity, with applications to nonnegative and fractional responses", \emph{Econometric Reviews}, 36(4), 397-420.
}

\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}

\seealso{
\code{\link{fracreghet}}, for fitting fractional response models under unobserved heterogeneity.\cr
\code{\link{fracreghet.reset}}, for the RESET test.\cr
}

\examples{
### Empirical 401(k) Examples 
data("fracreg_k401k") 
y <- fracreg_k401k$prate 
X_het <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age)
 
# fracreghet estimators do not allow exact 1s or 0s
y_adj <- y
y_adj[y_adj == 1] <- 0.999

# Artificial instrumental variable for demonstration 
set.seed(42)
Z_emp <- cbind(X_het, z = fracreg_k401k$mrate * rnorm(length(y))) 
res_emp <- fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], 
                      type="QMLxv", link="logit", table=FALSE) 
fracreghet.pe(res_emp, which.x="mrate")
 
### Simulated Examples

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
}
