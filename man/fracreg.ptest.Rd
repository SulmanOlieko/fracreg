\name{fracreg.ptest}
\alias{fracreg.ptest}

\title{
P Test for Fractional Response Models
}

\description{
\code{fracreg.ptest} is used to perform the P test to evaluate the specification of alternative, non-nested fractional response models by testing against each other.
}

\usage{fracreg.ptest(object1, object2, version = "Wald", table = TRUE)
}

\arguments{
  \item{object1}{an object containing the results of an \code{fracreg} command.
}
  \item{object2}{an object containing the results of another \code{fracreg} command.
}
  \item{version}{a vector containing the test versions to use. Available options: \code{Wald} (the default) and \code{LM}. 
Both options may be chosen at the same time and are computed in a robust way.
}
  \item{table}{a logical value indicating whether a summary table with the test results should be printed.
}
}

\details{
\code{fracreg.ptest} applies the P test statistic proposed by Davidson and
MacKinnon (1981) to fractional response models estimated via
\code{fracreg}. \code{fracreg.ptest} may be used to test against each other two alternative
specifications for the link function  in: (i) one-part fractional response models;
(ii) the binary components of two-part and three-part fractional response models; (iii) the
fractional components of two-part and three-part fractional response models; and (iv) two-part
and three-part fractional response models. In addition, \code{fracreg.ptest} may be used to test
one-part models against two-part or three-part models and in cases where the link functions are the
same but the regressors are non-nested. See Ramalho, Ramalho and Murteira (2011)
for details on the application of the P test in the fractional response framework.
}

\value{
\code{fracreg.reset} returns a named vector with the test results.
}

\references{
Davidson, R. and J.G. MacKinnon (1981), "Several tests for model specification
on the presence of alternative hypotheses", \emph{Econometrica}, 49(3), 781-793.

Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
estimating and testing empirical strategies for fractional response models",
\emph{Journal of Economic Surveys}, 25(1), 19-68.
}

\author{
Sulman Olieko Owili <oliekosulman@gmail.com>
}

\seealso{
\code{\link{fracreg}}, for fitting fractional response models.\cr
\code{\link{fracreg.reset}} and \code{\link{fracreg.ggoff}}, for specification tests.\cr
\code{\link{fracreg.pe}}, for computing partial effects.
}

\examples{
N <- 250
u <- rnorm(N)

X <- cbind(rnorm(N),rnorm(N))
dimnames(X)[[2]] <- c("X1","X2")

ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
y <- rbeta(N,ym*20,20*(1-ym))
y[y > 0.9] <- 1

#Testing logit versus loglog specifications for standard fractional
#regression models using a LM version of the P test
res1 <- fracreg(y,X,linkfrac="logit",table=FALSE)
res2 <- fracreg(y,X,linkfrac="loglog",table=FALSE)
fracreg.ptest(res1,res2,"LM")

#Testing a logit one-part fractional response model versus a binary logit +
#fractional probit two-part model using a Wald version of the P test
res1 <- fracreg(y,X,linkfrac="logit",table=FALSE)
res2 <- fracreg(y,X,linkbin="logit",linkfrac="probit",type="2P",inf=1,table=FALSE)
fracreg.ptest(res1,res2,"Wald")
}
