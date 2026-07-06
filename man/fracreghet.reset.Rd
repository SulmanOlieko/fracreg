\name{fracreghet.reset}
\alias{fracreghet.reset}

\title{
RESET Test for Fractional Response Models under Neglected Heterogeneity
}

\description{
\code{fracreghet.reset} is used to test the specification of fractional response models estimated by GMMx or LINx.
}

\usage{fracreghet.reset(object, lastpower.vec = 3, version = "Wald", table = T, ...)
}

\arguments{
  \item{object}{an object containing the results of an \code{fracreghet} command.
}
  \item{lastpower.vec}{a numeric vector containing the maximum powers of the linear predictors to be used in RESET tests.
}
  \item{version}{a vector containing the test versions to use. Available options: \code{Wald} (the default) and \code{LM} (only 
available for \code{GMMx}).
}
  \item{table}{a logical value indicating whether a summary table with the test results should be printed.
}
  \item{\dots}{Arguments to pass to \link[stats]{nlminb}, which is used to estimate the model under the alternative hypothesis when \code{version} is equal to \code{"Wald"} and the null model was estimated by \code{GMMx}.
}
}

\details{
\code{fracreghet.reset} applies the RESET test statistic to fractional response
models estimated via \code{fracreghet} using the options \code{GMMx} or \code{LINx}. \code{fracreghet.reset} may be used to test simultaneously the validity of the link specification and the transformation applied to the response variable by each estimator.  

\strong{RESET Test under Unobserved Heterogeneity:}
The test is based on augmenting the original model with powers of the linear predictor \eqn{x\hat{\beta}}. For GMMx, it tests \eqn{H_0: \gamma = 0} in the expanded moment conditions:
\deqn{E\left[Z_i \left(H(y_i) - \exp\left(x_i\beta + \sum_{k=2}^P \gamma_k (x_i\hat{\beta})^k\right)E(e^{c_i})\right)\right] = 0}
This simultaneously evaluates whether the mean function and the specific heterogeneity transformation \eqn{H(\cdot)} are correctly specified.

It is taken into account the option that was chosen for computing standard errors in the model under evaluation. See Ramalho and Ramalho (2017) for details.
}

\value{
\code{fracreghet.reset} returns a named vector with the test results.
}

\references{
Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of nonlinear regression models with boundary outcomes and endogeneity, with applications to nonnegative and fractional responses", \emph{Econometric Reviews}, 36(4), 397-420.

Ramsey, J.B. (1969), "Tests for Specification Errors in Classical Linear Least-Squares Regression Analysis", \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, 31(2), 350-371.
}

\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}

\seealso{
\code{\link{fracreghet}}, for fitting fractional response models under unobserved heterogeneity.\cr
\code{\link{fracreghet.pe}}, for computing partial effects.
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
res_emp <- fracreghet(y_adj, X_het, type="GMMx", link="logit", table=FALSE) 
fracreghet.reset(res_emp)
 
### Simulated Examples

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
