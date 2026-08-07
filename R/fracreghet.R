fracreghet.links <- function(link) 
{
	switch(link,
		logit = {
			linkfun <- function(mu) qlogis(mu)
			linkinv <- function(eta) plogis(eta)
			mu.eta <- function(eta) dlogis(eta)
			H1 <- function(dep) dep/(1-dep)
			G2 <- function(eta) exp(eta)
			g2 <- function(eta) exp(eta)
			H2 <- function(uu) log(uu)
			gd <- function(eta) exp(eta)*(1-exp(eta))/((1+exp(eta))^3)
			valideta <- function(eta) TRUE
		},

		probit = {
			linkfun <- function(mu) qnorm(mu)
			linkinv <- function(eta) {
				thresh <- -qnorm(.Machine$double.eps)
				eta <- pmin(pmax(eta, -thresh), thresh)
	 			pnorm(eta)
			}
			mu.eta <- function(eta) pmax(dnorm(eta),.Machine$double.eps)
			H1 <- NULL
			G2 <- NULL
			g2 <- NULL
			H2 <- NULL
			gd <- function(eta) -eta*dnorm(eta)
			valideta <- function(eta) TRUE
		},

		cauchit = {
		  	linkfun <- function(mu) qcauchy(mu)
			linkinv <- function(eta) {
				thresh <- -qcauchy(.Machine$double.eps)
				eta <- pmin(pmax(eta, -thresh), thresh)
				pcauchy(eta)
			}
			mu.eta <- function(eta) pmax(dcauchy(eta),.Machine$double.eps)
			H1 <- NULL
			G2 <- NULL
			g2 <- NULL
			H2 <- NULL
			gd <- function(eta) -2*eta/(pi*(eta^2+1)^2)
			valideta <- function(eta) TRUE
		},

		cloglog = {
			linkfun <- function(mu) log(-log(1 - mu))
			linkinv <- function(eta) pmax(pmin(-expm1(-exp(eta)),1-.Machine$double.eps),.Machine$double.eps)
			mu.eta <- function(eta) {
				eta <- pmin(eta, 700)
				pmax(exp(eta) * exp(-exp(eta)),.Machine$double.eps)
			}
			H1 <- function(dep) -log(1-dep)
			G2 <- function(eta) exp(eta)
			g2 <- function(eta) -exp(-eta)
			H2 <- function(uu) log(uu)
			gd <- function(eta) (exp(-exp(eta))*exp(eta))*(1-exp(eta))
			valideta <- function(eta) TRUE
		},

		loglog = {
			linkfun <- function(mu) -log(-log(mu))
			linkinv <- function(eta) exp(-exp(-eta))
			mu.eta <- function(eta) exp(-exp(-eta)-eta)
			H1 <- NULL
			G2 <- NULL
			g2 <- NULL
			H2 <- NULL
			gd <- function(eta) (exp(-exp(-eta))*exp(-eta))*(exp(-eta)-1)
			valideta <- function(eta) TRUE
		}
	)

	structure(list(linkfun=linkfun,linkinv=linkinv,mu.eta=mu.eta,H1=H1,G2=G2,g2=g2,H2=H2,gd=gd,valideta=valideta,name=link),class="link-glm")
}

fracreghet.gi <- function(type,z,Hy,XB,link)
{
	if(type!="QMLxv")
	{
		G2 <- fracreghet.links(link)$G2(XB)
		u <- Hy/G2-1
	}
	if(type=="QMLxv")
	{
		yhat <- fracreghet.links(link)$linkinv(XB)
		g <- fracreghet.links(link)$mu.eta(XB)
		u <- Hy-yhat
		u <- g*u/(yhat*(1-yhat))
	}

	gi <- t(z*u)

	return(list(gi=gi,u=u))
}

fracreghet.Gn <- function(type,x,z,z.in,u,bet,p)
{
	N <- length(u)

	if(all(type!=c("GMMxv","LINxv","QMLxv"))) Gn <- -(1/N)*t(z*bet)%*%x
	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		k <- ncol(x)
		
		G11 <- -(1/N)*t(z*bet)%*%x

		G12 <- (1/N)*t(x*(bet*p[k]))%*%z.in
		G12[k,] <- G12[k,]-(1/N)*apply(u*z.in,2,sum)

		G21 <- matrix(0,nrow=ncol(z.in),ncol=k)
		G22 <- -(1/N)*t(z.in)%*%z.in

		Gn <- rbind(cbind(G11,G12),cbind(G21,G22))
	}

	return(Gn)
}

fracreghet.est <- function(type,x,z,link,start,Hy,variance,var.type,var.cluster,gixv,vhat,offset=NULL,...)
{
	if(any(type==c("GMMxv","QMLxv")))
	{
		z.in <- z
		z <- x
	}

	GMM.est <- TRUE
	if(any(type==c("GMMx","GMMxv")) & !any(Hy==0))
	{
		results <- tryCatch(glm(Hy ~ x-1,family=Gamma(link=log),maxit=100,offset=offset),error=function(e) return(NULL))
		if(any(is.null(results))) converged <- FALSE
		else
		{
			p <- results$coefficients
			XB <- results$linear.predictors
			converged <- results$converged*(1-results$boundary)
		}

		if(converged==TRUE) GMM.est <- FALSE
	}
	if(GMM.est==TRUE)
	{
		GMMn <- function(p)
		{
			XB <- as.vector(x%*%p)
			if (!is.null(offset)) XB <- XB + offset
			gi <- fracreghet.gi(type,z,Hy,XB,link)$gi

			gn <- as.matrix(apply(gi,1,mean))
			Qn <- t(gn)%*%S%*%gn

			return(Qn)
		}

		S <- diag(ncol(z))
		results <- nlminb(start=start,objective=GMMn,...)
		p <- results$par
		XB <- as.vector(x%*%p)
		if (!is.null(offset)) XB <- XB + offset

		if(type=="GMMz" & ncol(z)>ncol(x))
		{
			fi.inv <- fracreghet.var(type,p,XB,x,z,link,Hy,var.type,var.cluster,TRUE,gixv,vhat)$fi.inv
			if(!is.character(fi.inv))
			{
				S <- fi.inv
				results <- nlminb(start=start,objective=GMMn,...)
				p <- results$par
				XB <- as.vector(x%*%p)
				if (!is.null(offset)) XB <- XB + offset
			}

			Qn <- results$objective
		}

		converged <- ifelse(results$convergence==0,TRUE,FALSE)
	}

	ret.list <- list(p=p,XB=XB,converged=converged)
	if(type=="GMMz" & ncol(z)>ncol(x)) ret.list[["Qn"]] <- Qn
	if (any(type == c("QMLxv", "QMLz"))) {
		yhat <- fracreghet.links(link)$linkinv(XB)
		eps <- 1e-16
		LL <- sum(ifelse(Hy > 0, Hy * log(pmax(yhat, eps)), 0) + ifelse(Hy < 1, (1-Hy) * log(pmax(1-yhat, eps)), 0))
		ret.list[["LL"]] <- LL
	}

	if(variance==FALSE | converged==FALSE) return(ret.list)

	if(any(type==c("GMMxv","QMLxv"))) z <- z.in
	p.var <- fracreghet.var(type,p,XB,x,z,link,Hy,var.type,var.cluster,FALSE,gixv,vhat)$p.var
	ret.list[["p.var"]] <- p.var

	return(ret.list)
}

fracreghet.var <- function(type,p,XB,x,z,link,Hy,var.type,id,step.one,gixv,vhat)
{
	N <- length(Hy)

	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		z.in <- z
		z <- x
	}
	else z.in <- z

	if(any(type==c("LINx","LINxv","LINz")))
	{
		u <- Hy-XB
		gi <- t(z*u)
		bet <- rep(1,N)
	}
	if(any(type==c("GMMx","GMMz","GMMxv","QMLxv")))
	{
		results <- fracreghet.gi(type,z,Hy,XB,link)
		gi <- results$gi
		u <- results$u

		if(type!="QMLxv")
		{
			G2 <- fracreghet.links(link)$G2(XB)
			g2 <- fracreghet.links(link)$g2(XB)
			bet <- Hy*g2/(G2^2)
		}
		if(type=="QMLxv")
		{
			yhat <- fracreghet.links(link)$linkinv(XB)
			g <- fracreghet.links(link)$mu.eta(XB)
			bet <- (g^2)/(yhat*(1-yhat))
		}
	}

	if(any(type==c("GMMxv","LINxv","QMLxv"))) gi <- rbind(gi,gixv)

	if(var.type=="robust") fi <- (1/N)*gi%*%t(gi)

	if(var.type=="cluster")
	{
		fi <- 0

		for(j in unique(id))
		{
			Zi <- matrix(z[id==j,],ncol=ncol(z))
			ui <- u[id==j]
			zu <- t(Zi)%*%ui

			if(any(type==c("GMMxv","LINxv","QMLxv")))
			{
				Zi <- matrix(z.in[id==j,],ncol=ncol(z.in))
				vi <- vhat[id==j]
				zu2 <- t(Zi)%*%vi
				zu <- rbind(zu,zu2)
			}

			fi <- fi+zu%*%t(zu)
		}

		fi <- fi/N
	}

	fi.inv <- tryCatch(solve(fi),error=function(e) NaN)
	if(any(is.nan(fi.inv))) fi.inv <- "singular"

	if(step.one==TRUE) return(list(fi.inv=fi.inv))

	Gn <- fracreghet.Gn(type,x,z,z.in,u,bet,p)

	if(is.numeric(fi.inv))
	{
		sigma <- N*t(Gn)%*%fi.inv%*%Gn
		p.var <- tryCatch(solve(sigma),error=function(e) NaN)
		if(any(is.nan(p.var))) p.var <- "singular"
	}
	else p.var <- "singular"

	if(is.character(p.var))
	{
		if(ncol(z)>ncol(x)) p.var <- "singular"
		else
		{
			Gn.inv <- tryCatch(solve(Gn),error=function(e) NaN)
			if(any(is.nan(Gn.inv))) p.var <- "singular"
			else p.var <- (1/N)*Gn.inv%*%fi%*%t(Gn.inv)
		}
	}

	ret.list <- list(p.var=p.var)

	return(ret.list)
}

#' @title Fitting Fractional Response Regressions under Unobserved Heterogeneity
#'
#' @description
#' \code{fracreghet} is used to fit fractional response models under unobserved heterogeneity, i.e. regression models for proportions, percentages or fractions that suffer from neglected heterogeneity and/or endogeneity issues.
#' @param y a numeric vector containing the values of the response variable.
#' @param x a numeric matrix, with column names, containing the values of all covariates (exogenous and endogenous).
#' @param z a numeric matrix, with column names, containing the values of all exogenous variables (covariates and instrumental  variables). Defaults to \code{x}.
#' @param var.endog a numeric vector containing the values of the endogenous covariate (or of some transformation of it), which will be used as dependent variable in the linear reduced form assumed for application of xv-type estimators.
#' @param start a numeric vector containing the initial values for the parameters to be optimised. Optional.
#' @param type a description of the estimator to compute: \code{GMMx} (the default), \code{GMMxv}, \code{GMMz}, \code{LINx}, \code{LINxv}, \code{LINz} or \code{QMLxv}.
#' @param link a description of the link function to use. Available options for all estimators: \code{logit} and \code{cloglog}. Additional available options for QML and LIN estimators: \code{probit}, \code{cauchit} and \code{loglog}.
#' @param intercept a logical value indicating whether the model should include a constant term or not.
#' @param table a logical value indicating whether a summary table with the regression results should be printed.
#' @param variance a logical value indicating whether the variance of the estimated parameters should be calculated. Defaults to \code{TRUE} whenever \code{table = TRUE}.
#' @param var.type a description of the type of variance of the estimated parameters to be calculated. Options are \code{robust}, the default, and \code{cluster}.
#' @param var.cluster a numeric vector containing the values of the variable that specifies to which cluster each observation belongs.
#' @param adjust the numeric value to be added to the response variable in case of boundary observations when the LIN estimators are applied or the string \code{drop}, which implies that the boundary observations are dropped.
#' @param offset an optional numeric vector containing an offset. It must be of the same dimension as the response variable. It specifies that the variable should be included in the model with its coefficient constrained to 1.
#' @param or a logical value indicating whether to report odds ratios. Only valid when the link function is \code{"logit"}. Defaults to \code{FALSE}.
#' @param level a numeric value between 0 and 1 indicating the confidence level for the confidence intervals. Defaults to \code{0.95}.
#' @param na.action A function specifying how to handle missing values, default is \code{stats::na.omit}. If \code{NULL}, no action is taken.
#' @param \dots Arguments to pass to \link[stats]{nlminb}.
#'
#' @details
#' \code{fracreghet} computes the GMM estimators proposed in Ramalho and Ramalho (2017) for fractional response models with unobserved heterogeneity: GMMx, which allows for neglected heterogeneity but not for endogeneity; GMMxv, which allows both issues and assumes a linear reduced form for the endogeneous covariate (or for a transformation of it); and GMMz, which also allows for both issues but does not require the assumption of a reduced form for the endogenous covariate. In addition, \code{fracreghet} also computes three linearised estimators (LINx, LINxv and LINz) that have similar features to their GMM counterparts. It also provides a QML estimator (QMLxv) that addresses endogeneity using a Control Function (CF) approach, which includes the first-stage reduced-form residuals as an additional regressor in the main fractional equation, providing a Hausman-type test for endogeneity. 
#' 
#' \strong{Control Function (CF) Approach - QMLxv:}
#' When a continuous regressor \eqn{y_{2i}} is endogenous, the CF approach (Papke and Wooldridge, 2008; Terza et al., 2008) uses a two-stage procedure. First, a linear reduced form is estimated:
#' \deqn{y_{2i} = z_i \pi + v_i}
#' where \eqn{z_i} includes all exogenous variables and external instruments. The residuals \eqn{\hat{v}_i} are then included in the fractional response model:
#' \deqn{E(y_{1i} | z_i, y_{2i}, v_i) = G(x_i \beta + \gamma \hat{v}_i)}
#' A test of \eqn{H_0: \gamma = 0} serves as a robust Hausman-type test for endogeneity.
#' 
#' \strong{Generalised Method of Moments (GMM):}
#' For estimators like GMMz, which do not strictly require a linear reduced form, the estimation relies on population orthogonality conditions between the instruments \eqn{Z_i} and the model residuals:
#' \deqn{E[Z_{i} (y_i - G(x_i \beta))] = 0}
#' or via specific transformations of the dependent variable to eliminate unobserved heterogeneity (Ramalho and Ramalho, 2017).
#' 
#' For overidentified models, \code{fracreghet} calculates Hansen's J statistic. For \code{GMMx} and \code{LINx}, \code{fracreghet} stores the information needed to implement the RESET test (\link{fracreghet.reset}). For all estimators, \code{fracreghet} stores the information needed to calculate partial effects (\link{fracreghet.pe}).
#'
#' @return
#' \code{fracreghet} returns a list with the following elements:
#'   \item{class}{"fracreghet".
#' }
#'   \item{formula}{the model formula.
#' }
#'   \item{type}{the name of the estimator computed.
#' }
#'   \item{link}{the name of the specified link.
#' }
#'   \item{adjust}{The value or the type of the adjustment applied to LIN estimators.
#' }
#'   \item{p}{a named vector of coefficients.
#' }
#'   \item{Hy}{the transformed values of the response variable when GMM or LIN estimators are computed or the 
#'    values of the response variable in the QML case.
#' }
#'   \item{xbhat}{the fitted mean values of the linear predictor (for xv-type estimators, includes the term relative to the first-stage residual).
#' }
#'   \item{converged}{logical. Was the algorithm judged to have converged?
#' }
#'   \item{x.names}{a vector containing the names of the covariates.
#' }
#' 
#' In case of an overidentifying model, the following element is also returned:
#'   \item{J}{the result of Hansen's J test of overidentifying moment conditions.
#' }
#' 
#' If \code{variance = TRUE} or \code{table = TRUE} and the algorithm converged successfully, the previous list also contains the following elements:
#' 
#'   \item{p.var}{a named covariance matrix.
#' }
#'   \item{var.type}{covariance matrix type.
#' }
#' 
#' If \code{var.type = "cluster"}, the list also contains the following element:
#'   \item{var.cluster}{the variable that specifies to which cluster each observation belongs.
#' }
#'
#' @section Odds Ratios:
#' When \code{or=TRUE} and the fractional link function (\code{linkfrac} or \code{link}) is \code{"logit"}, the model additionally computes odds ratios for the coefficients. 
#' Odds Ratios are exponentiated coefficients.
#' The corresponding standard errors for the odds ratios are calculated using the Delta method.
#' The confidence intervals for the odds ratios are calculated using the adjusted standard errors and the specified \code{level} (defaulting to 95\%).
#' Odds ratios are particularly useful in fractional logit models as they provide a direct multiplicative interpretation of the independent variable on the odds of the fractional outcome.
#'
#' @references
#' Papke, L. E. and Wooldridge, J. M. (2008), "Panel data methods for fractional response variables with an application to test pass rates", \emph{Journal of Econometrics}, 145, 121-133.
#' 
#' Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of nonlinear regression models with boundary outcomes and endogeneity, with applications to nonnegative and fractional responses", \emph{Econometric Reviews}, 36(4), 397-420.
#' 
#' Terza, J. V., Basu, A., and Rathouz, P. J. (2008), "Two-stage residual inclusion estimation: addressing endogeneity in health econometric modeling", \emph{Journal of Health Economics}, 27(3), 531-543.
#'
#' @author Sulman Olieko Owili <oliekosulman@gmail.com>
#'
#' @seealso
#' \code{\link{fracreghet.reset}}, for the RESET test.\cr
#' \code{\link{fracreghet.pe}}, for computing partial effects.\cr
#' \code{\link{fracreg}}, for fitting standard cross-sectional fractional response models.\cr
#' \code{\link{fracregpd}}, for fitting panel data fractional response models.
#'
#' @examples
#' ### Empirical 401(k) Examples 
#' data("fracreg_k401k") 
#' y <- fracreg_k401k$prate 
#' X_het <- cbind(mrate = fracreg_k401k$mrate, ltotemp = fracreg_k401k$ltotemp)
#'  
#' # fracreghet estimators do not allow exact 1s or 0s
#' y_adj <- y
#' y_adj[y_adj == 1] <- 0.999
#' 
#' # Instrument mrate using age
#' 
#' Z_emp <- cbind(age = fracreg_k401k$age, ltotemp = fracreg_k401k$ltotemp) 
#' mod <- fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], type="QMLxv", link="logit")
#' summary(mod)
#' 
#' # Compute the same QMLxv estimator reporting Odds Ratios with 90% confidence intervals
#' mod <- fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], type="QMLxv", 
#'            link="logit", or=TRUE, level=0.90) 
#'  
#' ### Simulated Examples 
#' 
#' set.seed(123)
#' N <- 1000
#' x1 <- rnorm(N)
#' 
#' # Simulating an endogenous variable (var.endog) and an instrument (z1)
#' z1 <- rnorm(N)
#' u <- 0.5 * z1 + rnorm(N)
#' var.endog <- 0.8 * z1 + u
#' y_endog <- exp(0.5 * x1 + 1.2 * var.endog + u) / (1 + exp(0.5 * x1 + 1.2 * var.endog + u))
#' 
#' # Avoid exact 0 or 1 boundaries for some estimators
#' y_endog[y_endog <= 0] <- 0.01
#' y_endog[y_endog >= 1] <- 0.99
#' 
#' X <- cbind(x1 = x1, var.endog = var.endog)
#' Z <- cbind(x1 = x1, z1 = z1)
#' 
#' # Exogeneity (assuming var.endog is exogenous for comparison), GMMx estimator
#' mod <- fracreghet(y = y_endog, x = X, type = "GMMx", link = "logit")
#' summary(mod)
#' 
#' # Endogeneity, GMMz estimator (does not require reduced form for endog)
#' mod <- fracreghet(y = y_endog, x = X, z = Z, type = "GMMz", link = "logit")
#' summary(mod)
#' 
#' # Endogeneity, GMMxv estimator (assumes linear reduced form for var.endog)
#' mod <- fracreghet(y = y_endog, x = X, z = Z, var.endog = var.endog, type = "GMMxv", link = "logit")
#' summary(mod)
#' 
#' # Endogeneity, QMLxv control function approach
#' mod <- fracreghet(y = y_endog, x = X, z = Z, var.endog = var.endog, type = "QMLxv", link = "logit")
#' summary(mod)
#' @export
fracreghet <- function(y,x,z=x,var.endog,start,type="GMMx",link="logit",intercept=TRUE,table=FALSE,variance=TRUE,var.type="robust",var.cluster,adjust=0,offset=NULL,or=FALSE,level=0.95,na.action=stats::na.omit,...)
{
	cl <- match.call()
	LL <- NULL
	
	if (!missing(y) && !missing(x)) {
	    args_to_clean <- list(y=y, x=x, offset=offset)
	    if (!missing(z)) args_to_clean$z <- z
	    if (!missing(var.endog)) args_to_clean$var.endog <- var.endog
	    if (!missing(var.cluster)) args_to_clean$var.cluster <- var.cluster
	    
	    args_to_clean$na.action <- na.action
	    cleaned <- do.call(fracreg_clean_data, args_to_clean)
	    
	    y <- cleaned$y
	    x <- cleaned$x
	    offset <- cleaned$offset
	    if (!missing(z)) z <- cleaned$z
	    if (!missing(var.endog)) var.endog <- cleaned$var.endog
	    if (!missing(var.cluster)) var.cluster <- cleaned$var.cluster
	}
	
	### 1. Error and warning messages

	if(missing(y)) stop("dependent variable is missing")
	if(missing(x)) stop("explanatory variables are missing")
	if(any(y>1) | any(y<0)) stop("The dependent variable has values outside the unit interval")
	if((any(y==1) | any(y==0)) & adjust==0 & any(type==c("LINx","LINz","LINxv"))) stop("0/1 values for the response variable: LIN estimators require adjustment")
	if(all(link!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(link)," - link not recognised")
	if(all(type!=c("GMMx","GMMz","GMMxv","LINx","LINz","LINxv","QMLxv"))) stop(sQuote(type)," - type not recognised")
	if(any(type==c("GMMx","GMMz","GMMxv")) & all(link!=c("logit","cloglog"))) stop("type and link not compatible")
	if(any(type==c("GMMx","GMMz","GMMxv")) & any(y==1)) stop("estimator does not allow y = 1")
	if(!is.numeric(adjust) & adjust!="drop") stop("adjust not defined properly")
	if(any(type==c("GMMxv","LINxv","QMLxv")) & missing(var.endog)) stop(sQuote(type)," requires var.endog to be specified")
	if(all(type!=c("GMMxv","LINxv","QMLxv")) & !missing(var.endog)) stop("var.endog should not be specified for this estimator")
	if(table==TRUE & variance==FALSE)
	{
		variance <- TRUE
		warning("option variance changed from FALSE to TRUE, as required by table=TRUE")
	}
	if(all(var.type!=c("robust","cluster"))) stop(sQuote(var.type)," - var.type not recognised")
	if(var.type=="cluster" & missing(var.cluster)) stop("option cluster for covariance matrix but no var.cluster supplied")
	if(var.type=="robust" & !missing(var.cluster)) stop("option robust for covariance matrix but var.cluster supplied")

	if(!is.logical(intercept)) stop("non-logical value assigned to option intercept")
	if(!is.logical(table)) stop("non-logical value assigned to option table")
	if(!is.logical(variance)) stop("non-logical value assigned to option variance")

	### 2. Data and variables preparation

	y <- as.vector(y)

	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.data.frame(z)) z <- as.matrix(z)

	if(!is.matrix(x)) stop("x is not a matrix")
	if(!is.matrix(z)) stop("z is not a matrix")

	x.names <- dimnames(x)[[2]]
	z.names <- dimnames(z)[[2]]

	if(is.null(x.names)) stop("x has no column names")
	if(is.null(z.names)) stop("z has no column names")

	if(intercept==TRUE)
	{
		x <- cbind(1,x)
		z <- cbind(1,z)
		x.names <- c("(Intercept)",x.names)
		z.names <- c("(Intercept)",z.names)
	}

	if(length(x.names)!=length(unique(x.names))) stop("some covariate names in x are identical")
	if(length(z.names)!=length(unique(z.names))) stop("some instrument names in z are identical")
	if(identical(x.names,z.names) & any(type==c("GMMz","GMMxv","LINz","LINxv","QMLxv"))) stop("instruments and covariates are identical")
	if(!identical(x.names,z.names) & any(type==c("GMMx","LINx"))) stop("instruments should not be specified for this estimator")
	if(length(x.names)>length(z.names)) stop("number of instruments not enough")

	if(adjust=="drop")
	{
		x <- as.matrix(x[y>0 & y<1,])
		z <- as.matrix(z[y>0 & y<1,])
		if(var.type=="cluster") var.cluster <- var.cluster[y>0 & y<1]
		if(!is.null(offset)) offset <- offset[y>0 & y<1]
		if(!missing(var.endog)) var.endog <- var.endog[y>0 & y<1]
		y <- y[y>0 & y<1]
	}

	if(any(length(y)!=c(nrow(x),nrow(z)))) stop("the number of observations for y, x and/or z is different")
	if(var.type=="cluster")
	{
		if(length(y)!=length(var.cluster)) stop("var.cluster does not have the appropriate dimension")
	}
	if(!missing(var.endog))
	{
		if(length(y)!=length(var.endog)) stop("var.endog does not have the appropriate dimension")
	}
	if(!is.null(offset) && length(y)!=length(offset)) stop("offset does not have the appropriate dimension")

	class <- "fracreghet"

	### 3. Estimation

	if(any(type==c("GMMx","LINx"))) z <- x

	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		results <- lm(var.endog ~ 0+z)
		PIhat <- results$coefficients
		vhat <- results$residuals

		gixv <- t(z*vhat)
		x <- cbind(x,vhat)
	}
	else
	{
		gixv <- NA
		vhat <- NA
	}

	k <- ncol(x)
	kz <- ncol(z)
	N <- nrow(x)
	J <- NA
	dfJ <- kz-k

	if(any(type==c("GMMx","GMMz","GMMxv","QMLxv")))
	{
		if(any(type==c("GMMx","GMMz","GMMxv"))) Hy <- fracreghet.links(link)$H1(y)
		if(type=="QMLxv") Hy <- y

		if(missing(start)) start <- rep(0,k)
		if(length(start)!=k) stop("start is not of the same dimension as the covariate vector (including vhat, in case of GMMxv or QMLxv)")
		
		results <- fracreghet.est(type,x,z,link,start,Hy,variance,var.type,var.cluster,gixv,vhat,offset=offset,...)
		p <- results$p
		converged <- results$converged

		XB <- results$XB
		if(any(type=="GMMz") & dfJ>0)
		{
			Qn <- results$Qn
			J <- N*Qn
		}
	}

	if(any(type==c("LINx","LINz","LINxv")))
	{
		if(is.numeric(adjust))
		{
			y <- y+adjust
			if(any(y>=1) | any(y<=0)) stop("After adjustment, the dependent variable has values outside or at the boundaries of the unit interval")
		}

		Hy <- fracreghet.links(link)$linkfun(y)

		if(any(type==c("LINx","LINxv")))
		{
			results <- lm(Hy ~ 0+x, offset=offset)
			p <- results$coefficients

			XB <- results$fitted.values
		}
		if(type==c("LINz"))
		{
			Hy_adj <- Hy
			if (!is.null(offset)) Hy_adj <- Hy - offset
			XZ <- t(x)%*%z
			p <- solve(XZ%*%t(XZ))%*%(XZ%*%(t(z)%*%Hy_adj))
			if(kz>k)
			{
				u <- as.vector(Hy_adj-x%*%p)

				gi <- t(z*u)
				fi <- (1/N)*gi%*%t(gi)
				fi.inv <- solve(fi)
 				p <- as.vector(solve(XZ%*%fi.inv%*%t(XZ))%*%(XZ%*%fi.inv%*%(t(z)%*%Hy_adj)))
			}

			XB <- as.vector(x%*%p)
			if(!is.null(offset)) XB <- XB + offset

			if(kz>k) J <- N*t(Hy-XB)%*%z%*%fi.inv%*%t(z)%*%(Hy-XB)
		}

		converged <- TRUE
		if(variance==TRUE) results <- fracreghet.var(type,p,XB,x,z,link,Hy,var.type,var.cluster,FALSE,gixv,vhat)
	}

	if(variance==TRUE & converged==TRUE) p.var <- results$p.var
	else p.var <- "singular"

	x.names.in <- x.names

	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		p <- c(p,PIhat)
		x.names <- c(x.names,"vhat",paste("Z",z.names,sep="_"))
	}
	names(p) <- x.names

	table.info <- list(p=p,p.var=p.var,x.names=x.names,type=type,link=link,converged=converged,N=N,var.type=var.type,adjust=adjust,k=k,J=J,dfJ=dfJ,LL=LL,or=or,level=level)

	if(table==TRUE) do.call(fracreghet.table, table.info)

	formula <- y ~ x - 1

	res <- list(call=cl, class=class,formula=formula,type=type,link=link,adjust=adjust,p=p,Hy=Hy,xbhat=as.vector(XB),converged=converged,x.names=x.names.in,table.info=table.info)
	if(!is.null(offset)) res[["offset"]] <- offset
	if(any(type==c("GMMz","LINz")) & kz>k) res[["J"]] <- J

	if(variance==TRUE & converged==TRUE)
	{
		if(is.character(p.var)) p.var <- matrix(NA,nrow=length(p),ncol=length(p))
		dimnames(p.var) <- list(x.names,x.names)
		res[["p.var"]] <- p.var
		res[["var.type"]] <- var.type
		if(var.type=="cluster") res[["var.cluster"]] <- var.cluster
	}

	### 4. Return results

	class(res) <- "fracreghet"
	return(invisible(res))
}
