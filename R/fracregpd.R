fracregpd.links <- function(link) 
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
		}
	)

	structure(list(linkfun=linkfun,linkinv=linkinv,mu.eta=mu.eta,H1=H1,G2=G2,g2=g2,H2=H2,gd=gd,valideta=valideta,name=link),class="link-glm")
}

fracregpd.gi <- function(type,id,Ti,Hy,z,XB,link,at,at1)
{
	if(any(type==c("GMMc","GMMww","GMMbgw")) | (type=="GMMpfe" & !any(Hy==0)))
	{
		XBexp <- exp(XB)

		if(type=="GMMc") u <- (XBexp[at1]/XBexp[at])*Hy[at]-Hy[at1]
		if(type=="GMMww")
		{
			u <- Hy/XBexp
			u <- u[at]-u[at1]
		}
		if(type=="GMMbgw")
		{
			Hy.m <- rep(as.vector(by(Hy[at],id,mean)),times=Ti)
			XBexp.m <- rep(as.vector(by(XBexp[at],id,mean)),times=Ti)
			u <- Hy[at]-(XBexp[at]/XBexp.m)*Hy.m
		}
		if(type=="GMMpfe")
		{
			HyXBexp.m <- rep(as.vector(by(Hy[at]/XBexp[at],id,mean)),times=Ti)
			u <- Hy[at]/(XBexp[at]*HyXBexp.m)-1
		}
	}
	if(any(type==c("GMMpre","GMMcre")) | (type=="GMMpfe" & any(Hy==0)))
	{
		G2 <- fracregpd.links(link)$G2(XB)
		u <- Hy[at]/G2[at]-1
	}
	if(type=="QMLcre")
	{
		yhat <- fracregpd.links(link)$linkinv(XB)
		g <- fracregpd.links(link)$mu.eta(XB)
		u <- (Hy-yhat)*g/(yhat*(1-yhat))
	}

	gi <- t(z*u)

	return(list(gi=gi,u=u))
}

fracregpd.Gn <- function(type,x.exogenous,id,Ti,Hy,x,z,XB,link,at,at1,NT,k,p,z.in,u)
{
	if(any(type==c("GMMc","GMMww","GMMbgw")) | (type=="GMMpfe" & !any(Hy==0)))
	{
		XBexp <- exp(XB)

		if(type=="GMMc") Gn <- t(z*((XBexp[at1]/XBexp[at])*Hy[at]))%*%(x[at1,]-x[at,])
		if(type=="GMMww")
		{
			uu <- -Hy/XBexp
			Gn <- t(z)%*%(uu[at]*x[at,]-uu[at1]*x[at1,])
		}
		if(type=="GMMbgw")
		{
			Hy.m <- rep(as.vector(by(Hy[at],id,mean)),times=Ti)
			XBexp.m <- rep(as.vector(by(XBexp[at],id,mean)),times=Ti)

			XBexpX.m <- matrix(NA,nrow=NT,ncol=k)
			for(j in 1:k) XBexpX.m[,j] <- rep(as.vector(by(XBexp[at]*x[at,j],id,mean)),times=Ti)

			Gn <- -t(z*(XBexp[at]*Hy.m/XBexp.m))%*%x[at,]+t(z)%*%((XBexp[at]*Hy.m/((XBexp.m)^2))*XBexpX.m)
		}
		if(type=="GMMpfe")
		{
			HyXBexp.m <- rep(as.vector(by(Hy[at]/XBexp[at],id,mean)),times=Ti)

			HyXBexpX.m <- matrix(NA,nrow=NT,ncol=k)
			for(j in 1:k) HyXBexpX.m[,j] <- rep(as.vector(by((Hy[at]/XBexp[at])*x[at,j],id,mean)),times=Ti)
			Gn <- -t(z*(Hy[at]/(XBexp[at]*HyXBexp.m)))%*%x[at,]+t(z)%*%(Hy[at]/(XBexp[at]*(HyXBexp.m^2))*HyXBexpX.m)
		}
	}
	if(any(type==c("GMMpre","GMMcre")) | (type=="GMMpfe" & any(Hy==0)))
	{
		G2 <- fracregpd.links(link)$G2(XB)
		g2 <- fracregpd.links(link)$g2(XB)
		bet <- Hy*g2/(G2^2)
		Gn <- -t(z*bet[at])%*%x[at,]
	}
	if(type=="QMLcre")
	{
		yhat <- fracregpd.links(link)$linkinv(XB)
		g <- fracregpd.links(link)$mu.eta(XB)
		bet <- g^2/(yhat*(1-yhat))
		Gn <- -t(z*bet)%*%x

		if(x.exogenous==F)
		{
			G12 <- t(x*(bet*p[k]))%*%z.in
			G12[k,] <- G12[k,]-apply(u*z.in,2,sum)

			G21 <- matrix(0,nrow=ncol(z.in),ncol=k)
			G22 <- -t(z.in)%*%z.in

			Gn <- rbind(cbind(Gn,G12),cbind(G21,G22))
		}
	}

	Gn <- Gn/NT

	return(Gn)
}

fracregpd.est <- function(type,x.exogenous,lags,id,Ti,Hy,x,z,link,var.type,start,at,at1,variance,NT,k,kz,gixv,vhat,bootstrap,offset=NULL,...)
{
	if(type=="QMLcre" & x.exogenous==F)
	{
		z.in <- z
		z <- x
		kz <- ncol(z)
	}

	GMM.est <- T
	if(x.exogenous==T & ((any(type==c("GMMpre","GMMcre")) & !any(Hy==0) & lags==F) | type=="QMLcre"))
	{
		if(type=="QMLcre") results <- tryCatch(glm(Hy ~ x-1,family=quasibinomial(link=fracregpd.links(link)),maxit=100,offset=offset),error=function(e) return(NULL))
		if(any(type==c("GMMpre","GMMcre"))) results <- tryCatch(glm(Hy ~ x-1,family=Gamma(link=log),maxit=100,offset=offset),error=function(e) return(NULL))

		if(any(is.null(results))) converged <- F
		else
		{
			p <- results$coefficients
			XB <- results$linear.predictors
			converged <- results$converged*(1-results$boundary)
		}

		if(converged==T) GMM.est <- F
	}
	if(GMM.est==T)
	{
		GMMn <- function(p)
		{
			XB <- as.vector(x%*%p)
			if(!is.null(offset)) XB <- XB + offset
			gi <- fracregpd.gi(type,id,Ti,Hy,z,XB,link,at,at1)$gi
			gn <- as.matrix(apply(gi,1,mean))
			Qn <- t(gn)%*%S%*%gn

			return(Qn)
		}

		S <- diag(kz)
		results <- nlminb(start=start,objective=GMMn,...)
		p <- results$par
		XB <- as.vector(x%*%p)
		if(!is.null(offset)) XB <- XB + offset

		if(type!="QMLcre" & kz>k)
		{
			fi.inv <- fracregpd.var(type,x.exogenous,var.type,id,Ti,Hy,x,z,XB,link,at,at1,NT,k,kz,T,gixv,vhat,p)$fi.inv
			if(!is.character(fi.inv))
			{
				S <- fi.inv
				results <- nlminb(start=start,objective=GMMn,...)
				p <- results$par
				XB <- as.vector(x%*%p)
				if(!is.null(offset)) XB <- XB + offset
			}

			Qn <- results$objective
		}

		converged <- ifelse(results$convergence==0,T,F)
	}

	ret.list <- list(p=p,converged=converged)
	if(type!="QMLcre" & kz>k) ret.list[["Qn"]] <- Qn
	if (type == "QMLcre") {
		yhat <- fracregpd.links(link)$linkinv(XB)
		eps <- 1e-16
		LL <- sum(ifelse(Hy > 0, Hy * log(pmax(yhat, eps)), 0) + ifelse(Hy < 1, (1-Hy) * log(pmax(1-yhat, eps)), 0))
		ret.list[["LL"]] <- LL
	}
	
	ret.list[["xbhat"]] <- XB

	if(variance==F | converged==F | bootstrap==T) return(ret.list)

	if(type=="QMLcre" & x.exogenous==F) z <- z.in
	p.var <- fracregpd.var(type,x.exogenous,var.type,id,Ti,Hy,x,z,XB,link,at,at1,NT,k,kz,F,gixv,vhat,p)$p.var
	ret.list[["p.var"]] <- p.var

	return(ret.list)
}

fracregpd.var <- function(type,x.exogenous,var.type,id,Ti,Hy,x,z,XB,link,at,at1,NT,k,kz,step.one,gixv,vhat,p)
{
	if(type=="QMLcre" & x.exogenous==F)
	{
		z.in <- z
		z <- x
		kz <- ncol(z)
	}
	else z.in <- z

	results <- fracregpd.gi(type,id,Ti,Hy,z,XB,link,at,at1)
	u <- results$u

	if(var.type=="robust")
	{
		gi <- results$gi
		if(type=="QMLcre" & x.exogenous==F) gi <- rbind(gi,gixv)
		fi <- (1/NT)*gi%*%t(gi)
	}

	if(var.type=="cluster")
	{
		fi <- 0

		for(j in unique(id))
		{
			Zi <- matrix(z[id==j,],ncol=kz)
			ui <- u[id==j]
			zu <- t(Zi)%*%ui

			if(type=="QMLcre" & x.exogenous==F)
			{
				Zi <- matrix(z.in[id==j,],ncol=ncol(z.in))
				vi <- vhat[id==j]
				zu2 <- t(Zi)%*%vi
				zu <- rbind(zu,zu2)
			}

			fi <- fi+zu%*%t(zu)
		}

		fi <- fi/NT
	}

	fi.inv <- tryCatch(solve(fi),error=function(e) NaN)
	if(any(is.nan(fi.inv))) fi.inv <- "singular"

	if(step.one==T) return(list(fi.inv=fi.inv))

	Gn <- fracregpd.Gn(type,x.exogenous,id,Ti,Hy,x,z,XB,link,at,at1,NT,k,p,z.in,u)

	if(is.numeric(fi.inv))
	{
		sigma <- NT*t(Gn)%*%fi.inv%*%Gn
		p.var <- tryCatch(solve(sigma),error=function(e) NaN)
		if(any(is.nan(p.var))) p.var <- "singular"
	}
	else p.var <- "singular"

	if(is.character(p.var))
	{
		if(k!=kz) p.var <- "singular"
		else
		{
			Gn.inv <- tryCatch(solve(Gn),error=function(e) NaN)
			if(any(is.nan(Gn.inv))) p.var <- "singular"
			else p.var <- (1/NT)*Gn.inv%*%fi%*%t(Gn.inv)
		}
	}

	ret.list <- list(p.var=p.var)

	return(ret.list)
}

#' @title Fitting Panel Data Fractional Response Regressions
#'
#' @description
#' \code{fracregpd} is used to fit panel data regression models when the dependent variable has a bounded, fractional nature.
#' @param id a numeric vector identifying the cross-sectional units.
#' @param time a numeric vector identifying the time periods in which the cross-sectional units were observed.
#' @param y a numeric vector containing the values of the response variable.
#' @param x a numeric matrix, with column names, containing the values of all covariates (exogenous and endogenous).
#' @param z a numeric matrix, with column names, containing the values of all exogenous variables (covariates and external instrumental variables). Only required in case of endogenous explanatory variables.
#' @param var.endog a numeric vector containing the values of the endogenous covariate (or of some transformation of it), which will be used as dependent variable in the linear reduced form assumed for application of the \code{QMLcre} estimator. Only required for this estimator.
#' @param x.exogenous a logical value indicating whether all explanatory variables are assumed to be exogenous or not.
#' @param lags a logical value indicating whether the first lags of \code{x} or  \code{z} should be used as instruments for \code{x}. Defaults to \code{TRUE} for the GMMww and GMMc estimators and to \code{FALSE} for the remaining estimators. The \code{GMMcre} and \code{QMLcre} estimators do not admit lagged instruments.
#' @param start a numeric vector containing the initial values for the parameters to be optimised. Optional.
#' @param type a description of the estimator to compute: \code{GMMww}, \code{GMMc}, \code{GMMbgw}, \code{GMMpfe}, \code{GMMcre}, \code{GMMpre} or \code{QMLcre}.
#' @param GMMww.cor a logical value indicating whether each explanatory variable should be transformed in deviations from its overall mean before computing the \code{GMMww} estimator.
#' @param link a description of the link function to use. Available options for all GMM estimators: \code{logit} and \code{cloglog}. Only option for the \code{QMLcre} estimator: \code{probit}.
#' @param intercept a logical value indicating whether the model should include a constant term or not. Only relevant for the  \code{GMMpre} estimator.
#' @param table a logical value indicating whether a summary table with the regression results should be printed.
#' @param variance a logical value indicating whether the variance of the estimated parameters should be calculated. Defaults to  \code{TRUE} whenever \code{table = FALSERUE}.
#' @param var.type a description of the type of variance of the estimated parameters to be calculated. Options are \code{cluster}, the default, and \code{robust}. In overidentified models, it also affects the parameter estimates via the GMM weighting matrix.
#' @param tdummies a logical value indicating whether time dummies should be included among the model explanatory variables.
#' @param bootstrap a logical value indicating whether bootstrap should be used in the estimation of the parameter standard errors.
#' @param B the number of bootstrap replications.
#' @param offset an optional numeric vector containing an offset. It must be of the same dimension as the response variable. It specifies that the variable should be included in the model with its coefficient constrained to 1.
#' @param or a logical value indicating whether to report odds ratios. Only valid when the link function is \code{"logit"}. Defaults to \code{FALSE}.
#' @param level a numeric value between 0 and 1 indicating the confidence level for the confidence intervals. Defaults to \code{0.95}.
#' @param na.action A function specifying how to handle missing values, default is \code{stats::na.omit}. If \code{NULL}, no action is taken.
#' @param \dots Arguments to pass to \link[stats]{nlminb}.
#'
#' @details
#' \code{fracregpd} computes the GMM estimators proposed in Ramalho, Ramalho and Coelho (2018) for panel data fractional response models with both time-variant and time-invariant unobserved heterogeneity and endogeneous covariates: GMMww, GMMc, GMMbgw, GMMpfe, GMMcre and GMMpre. In addition, \code{fracregpd} also computes QMLcre, which was proposed by Papke and Wooldridge (2008) and Wooldridge (2019). 
#' 
#' \strong{Correlated Random Effects (CRE) - QMLcre:}
#' In panel data, unobserved individual-specific heterogeneity \eqn{c_i} may be correlated with the covariates \eqn{x_{it}}. The CRE approach (Papke and Wooldridge, 2008) models this dependence by projecting \eqn{c_i} onto the time averages of the strictly exogenous covariates \eqn{\bar{x}_i}:
#' \deqn{c_i = \psi + \bar{x}_i \xi + a_i}
#' where \eqn{a_i} is an error term independent of \eqn{x_i}. Assuming \eqn{a_i | x_i \sim N(0, \sigma_a^2)} and a probit link, integrating out \eqn{a_i} yields the "population-averaged" or scaled conditional mean:
#' \deqn{E(y_{it} | x_i) = G(x_{it} \beta_a + \psi_a + \bar{x}_i \xi_a)}
#' where the parameters with subscript \eqn{a} are scaled by \eqn{(1 + \sigma_a^2)^{-1/2}}. This equation is estimated via pooled Bernoulli QML.
#' 
#' \strong{Generalised Method of Moments (GMM):}
#' For models where strict exogeneity fails or the link function is an exponential-type link, Ramalho et al. (2018) propose GMM estimators based on the following general moment conditions:
#' \deqn{E[Z_{it} (H(y_{it}) - \exp(x_{it}\beta + c_i))] = 0}
#' where \eqn{H(\cdot)} is a transformation function and \eqn{Z_{it}} is a matrix of valid instruments. Estimators such as GMMww, GMMc, and GMMbgw use different transformations to eliminate the unobserved fixed effect \eqn{c_i} before applying GMM.
#' 
#' For overidentified models, \code{fracregpd} calculates Hansen's J statistic to test the validity of the overidentifying restrictions.
#'
#' @return
#' \code{fracregpd} returns a list with the following elements:
#'   \item{type}{the name of the estimator computed.
#' }
#'   \item{link}{the name of the specified link.
#' }
#'   \item{p}{a named vector of coefficients.
#' }
#'   \item{Hy}{the transformed values of the response variable when GMM estimators are computed or the 
#'    values of the response variable in the QML case.
#' }
#'   \item{converged}{logical. Was the algorithm judged to have converged?
#' }
#' 
#' In case of an overidentifying model, the following element is also returned:
#'   \item{J}{the result of Hansen's J test of overidentifying moment conditions.
#' }
#' 
#' If \code{variance = TRUE} or \code{table = FALSERUE} and the algorithm converged successfully, the previous list also contains the following elements:
#'   \item{p.var}{a named covariance matrix.
#' }
#'   \item{var.type}{covariance matrix type.
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
#' Papke, L. and Wooldridge, J.M. (2008), "Panel data methods for fractional response variables with an application to test pass rates", \emph{Journal of Econometrics}, 145(1-2), 121-133.
#' 
#' Ramalho, E. A., Ramalho, J. J. S., & Coelho, L. M. S. (2018), "Exponential Regression of Fractional-Response Fixed-Effects Models with an Application to Firm Capital Structure", \emph{Journal of Econometric Methods}, 7(1), 20150019.
#' 
#' Wooldridge, J. M. (2019). Correlated random effects models with unbalanced panels. \emph{Journal of Econometrics}, 211(1), 137-150.
#'
#' @author Sulman Olieko Owili <oliekosulman@gmail.com>
#'
#' @seealso
#' \code{\link{fracreg}}, for fitting standard cross-sectional fractional response models.\cr
#' \code{\link{fracreghet}}, for fitting cross-sectional fractional response models with unobserved heterogeneity.
#'
#' @examples
#' ### Empirical 401(k) Examples
#' data("fracreg_k401k")
#' y <- fracreg_k401k$prate
#' X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
#'            totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)
#' 
#' # Artificial panel data structure for demonstration
#' N_emp <- nrow(X)
#' id_emp <- rep(1:(N_emp/2), each=2)
#' time_emp <- rep(1:2, times=N_emp/2)
#' mod <- fracregpd(id_emp, time_emp, y, X, type="QMLcre", link="probit")
#' summary(mod)
#' 
#' ### Simulated Examples
#' 
#' set.seed(123)
#' # Simulating Panel Data
#' N <- 100
#' T_periods <- 5
#' id <- rep(1:N, each = T_periods)
#' time <- rep(1:T_periods, times = N)
#' x_panel <- rnorm(N * T_periods)
#' 
#' # Unobserved individual effect (CRE)
#' c_i <- rep(rnorm(N), each = T_periods) 
#' y_panel <- exp(x_panel + c_i) / (1 + exp(x_panel + c_i))
#' 
#' X <- cbind(x_panel = x_panel)
#' 
#' # Endogenous variable and instrument simulation
#' z_panel <- rnorm(N * T_periods)
#' u_panel <- 0.5 * z_panel + rnorm(N * T_periods)
#' var_endog <- 0.8 * z_panel + u_panel
#' y_endog <- exp(x_panel + 1.2 * var_endog + c_i + u_panel) / 
#'              (1 + exp(x_panel + 1.2 * var_endog + c_i + u_panel))
#' 
#' X_endog <- cbind(x_panel = x_panel, var_endog = var_endog)
#' Z_inst <- cbind(x_panel = x_panel, z_panel = z_panel)
#' 
#' \donttest{
#' # Estimate a Correlated Random Effects (CRE) Model
#' mod <- fracregpd(id=id, time=time, y=y_panel, x=X, type="QMLcre", link="probit")
#' summary(mod)
#' 
#' # Compute Partial Effects
#' pe_res <- fracregpd.pe(mod)
#' summary(pe_res)
#' 
#' # Exogeneity, no lags, no time dummies, clustered standard errors, GMMbgw estimator
#' mod <- fracregpd(id=id, time=time, y=y_panel, x=X, type="GMMbgw")
#' summary(mod)
#' 
#' # Estimate the GMMww estimator with odds ratios and 99% confidence intervals
#' mod <- fracregpd(id=id, time=time, y=y_panel, x=X, type="GMMww", or=TRUE, level=0.99)
#' summary(mod)
#' 
#' # Lagged covariates and instruments, robust standard errors, GMMww estimator
#' mod <- fracregpd(id=id, time=time, y=y_panel, x=X, lags=TRUE, type="GMMww", var.type="robust")
#' summary(mod)
#' 
#' # Endogeneity, time dummies, GMMpfe estimator
#' mod <- fracregpd(id=id, time=time, y=y_endog, x=X_endog, z=Z_inst,
#'                  x.exogenous=FALSE, type="GMMpfe", tdummies=TRUE)
#' summary(mod)
#' }
#' @export
fracregpd <- function(id,time,y,x,z,var.endog,x.exogenous=T,lags,start,type,GMMww.cor=T,link="logit",intercept=T,table=FALSE,variance=T,var.type="cluster",tdummies=F,bootstrap=F,B=200,offset=NULL,or=FALSE,level=0.95,na.action=stats::na.omit,...)
{
	cl <- match.call()
	
	if (!missing(id) && !missing(time) && !missing(y) && !missing(x)) {
	    args_to_clean <- list(id=id, time=time, y=y, x=x, offset=offset)
	    if (!missing(z)) args_to_clean$z <- z
	    if (!missing(var.endog)) args_to_clean$var.endog <- var.endog
	    
	    args_to_clean$na.action <- na.action
	    cleaned <- do.call(fracreg_clean_data, args_to_clean)
	    
	    id <- cleaned$id
	    time <- cleaned$time
	    y <- cleaned$y
	    x <- cleaned$x
	    offset <- cleaned$offset
	    if (!missing(z)) z <- cleaned$z
	    if (!missing(var.endog)) var.endog <- cleaned$var.endog
	}
	
	### 1. Error and warning messages

	if(missing(id)) stop("variable id is missing")
	if(missing(time)) stop("variable time is missing")
	if(missing(y)) stop("dependent variable is missing")
	if(missing(x)) stop("explanatory variables are missing")

	if(missing(type)) stop("type must be specified")
	if(any(y>1) | any(y<0)) stop("The dependent variable has values outside the unit interval")
	if(!is.logical(x.exogenous)) stop("x.exogenous must be logical")
	if(!is.logical(GMMww.cor)) stop("GMMww.cor must be logical")
	if(!is.logical(intercept)) stop("intercept must be logical")
	if(!is.logical(table)) stop("table must be logical")
	if(!is.logical(variance)) stop("variance must be logical")
	if(!is.logical(tdummies)) stop("tdummies must be logical")
	if(all(link!=c("logit","probit","cloglog"))) stop(sQuote(link)," - link not recognised")
	if(all(type!=c("GMMc","GMMww","GMMbgw","GMMpre","GMMpfe","GMMcre","QMLcre"))) stop(sQuote(type)," - type not recognised")
	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpre","GMMpfe","GMMcre")) & all(link!=c("logit","cloglog"))) stop("type and link not compatible")
	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpre","GMMpfe","GMMcre")) & any(y==1)) stop("estimator does not allow y = 1")
	if(type=="QMLcre" & link!="probit") stop("type and link not compatible")
	if(!missing(z) & x.exogenous==T) stop("z must not be specified in case of exogeneity")
	if(missing(z) & x.exogenous==F) stop("z needs to be specified in case of endogeneity")
	if(type=="QMLcre" & x.exogenous==F & missing(var.endog)) stop("QMLcre under endogeneity requires var.endog to be specified")
	if((type!="QMLcre" | x.exogenous==T) & !missing(var.endog)) stop("var.endog should not be specified for this estimator")
	if(table==T & variance==F)
	{
		variance <- T
		warning("option variance changed from F to T, as required by table=T")
	}
	if(all(var.type!=c("robust","cluster"))) stop(sQuote(var.type)," - var.type not recognised")

	if(any(unlist(by(time,id,duplicated)))) stop("variable time has duplicated entries for the same id")

	### 2. Data and variables preparation - general

	id <- as.vector(id)
	time <- as.vector(time)
	y <- as.vector(y)
	if(missing(z)) z <- x

	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.data.frame(z)) z <- as.matrix(z)

	if(!is.matrix(x)) stop("x is not a matrix")
	if(!is.matrix(z)) stop("z is not a matrix")

	if(any(is.na(id))) stop("id has missing values")
	if(any(is.na(time))) stop("time has missing values")
	if(any(is.na(y))) stop("y has missing values")
	if(any(is.na(x))) stop("x has missing values")
	if(any(is.na(z))) stop("z has missing values")

	x.names <- dimnames(x)[[2]]
	z.names <- dimnames(z)[[2]]

	if(is.null(x.names)) stop("x has no column names")
	if(is.null(z.names)) stop("z has no column names")

	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpfe","GMMcre","QMLcre"))) intercept <- F
	if(intercept==T)
	{
		x <- cbind(1,x)
		z <- cbind(1,z)
		x.names <- c("(Intercept)",x.names)
		z.names <- c("(Intercept)",z.names)
	}

	if(length(x.names)!=length(unique(x.names))) stop("some covariate names in x are identical")
	if(length(z.names)!=length(unique(z.names))) stop("some instrument names in z are identical")
	if(identical(x.names,z.names) & x.exogenous==F) stop("instruments and covariates are identical")
	if(!identical(x.names,z.names) & x.exogenous==T) stop("instruments should not be specified for this estimator")
	if(length(x.names)>length(z.names)) stop("number of instruments not enough")

	if(any(length(y)!=c(nrow(x),nrow(z)))) stop("the number of observations for y, x and/or z is different")
	if(!missing(var.endog))
	{
		if(length(y)!=length(var.endog)) stop("var.endog does not have the appropriate dimension")
	}
	if(!is.null(offset) && length(y)!=length(offset)) stop("offset does not have the appropriate dimension")

	### 3. Data and variables preparation - panel / estimator specifics

	if(missing(lags))
	{
		if(any(type==c("GMMc","GMMww"))) lags <- T
		if(any(type==c("GMMbgw","GMMpre","GMMpfe","GMMcre","QMLcre"))) lags <- F
	}
	if(!is.logical(lags)) stop("lags must be logical")
	if(any(type==c("GMMcre","QMLcre")) & lags==T) stop("GMMcre/QMLcre cannot be used with lagged instruments")

	N.ini <- length(unique(id))
	Ti.ini <- as.vector(by(id,id,length))
	NT.ini <- length(id)

	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpfe","GMMcre","QMLcre")))
	{
		for(j in 1:ncol(x))
		{
			xm <- rep(as.vector(by(x[,j],id,mean)),times=Ti.ini)
			x.dif <- x[,j]-xm

			if(all(x.dif==0)) stop("Time-invariant covariates not allowed")

			if(any(type==c("GMMcre","QMLcre")) & x.exogenous==T)
			{
				if(j==1) x.m <- cbind(xm)
				if(j>1) x.m <- cbind(x.m,xm)
			}
		}
		for(j in 1:ncol(z))
		{
			zm <- rep(as.vector(by(z[,j],id,mean)),times=Ti.ini)
			z.dif <- z[,j]-zm

			if(all(z.dif==0)) stop("Time-invariant instruments not allowed")

			if(any(type==c("GMMcre","QMLcre")) & x.exogenous==F)
			{
				if(j==1) z.m <- cbind(zm)
				if(j>1) z.m <- cbind(z.m,zm)
			}
		}
	}

	T.ord <- rep(NA,NT.ini)

	if(tdummies==T)
	{
		if(lags==T) out <- 2
		else out <- 1

		td <- matrix(0,nrow=NT.ini,ncol=length(unique(time))-out)
	}

	a <- 0
	for(j in sort(unique(time)))
	{
		a <- a+1
		T.ord[time==j] <- a
		if(tdummies==T)
		{
			if(a>out)
			{
				td[time==j,a-out] <- 1
				if(a==out+1) td.names <- paste("time",j,sep=".")
				if(a>(out+1)) td.names <- c(td.names,paste("time",j,sep="."))
			}
		}
	}

	if(type=="GMMpfe" & any(y==0))
	{
		D <- matrix(0,nrow=NT.ini,ncol=N.ini)
		a <- 0
		for(j in sort(unique(id)))
		{
			a <- a+1
			D[id==j,a] <- 1
			if(a==1) D.names <- paste("D",j,sep=".")
			if(a>1) D.names <- c(D.names,paste("D",j,sep="."))
		}

	}

	ord <- order(id,time)
	id <- id[ord]
	T.ord <- T.ord[ord]
	y <- y[ord]
	x <- cbind(x[ord,])
	if(tdummies==T) td <- cbind(td[ord,])
	if(type=="GMMpfe" & any(y==0)) D <- D[ord,]
	z <- cbind(z[ord,])
	if(!is.null(offset)) offset <- offset[ord]

	if(any(type==c("GMMc","GMMww")) | lags==T)
	{
		at <- rep(F,NT.ini)
		at1 <- rep(F,NT.ini)

		if(any(Ti.ini!=max(Ti.ini)))
		{
			keep <- rep(F,NT.ini)

			for(j in 1:NT.ini)
			{
				if(j==1 & (T.ord[2]-T.ord[1])==1 & id[1]==id[2]) at1[1] <- T
				if(j>1 & j<NT.ini)
				{
					if((T.ord[j]-T.ord[j-1])==1 & id[j]==id[j-1]) at[j] <- T
					if((T.ord[j+1]-T.ord[j])==1 & id[j]==id[j+1]) at1[j] <- T
				}
				if(j==NT.ini & (T.ord[NT.ini]-T.ord[NT.ini-1])==1 & id[NT.ini]==id[NT.ini-1]) at[NT.ini] <- T
			}

			keep <- at==T | at1==T
			id <- id[keep==T]
			y <- y[keep==T]
			x <- cbind(x[keep==T,])
			if(tdummies==T) td <- cbind(td[keep==T,])
			if(type=="GMMpfe" & any(y==0)) D <- D[keep==T,]
			z <- cbind(z[keep==T,])
			at <- at[keep==T]
			at1 <- at1[keep==T]
			if(!is.null(offset)) offset <- offset[keep==T]
		}
		else
		{
			at[T.ord!=1] <- T
			at1[T.ord!=max(Ti.ini)] <- T			
		}

		z.full <- z
		id.full <- id

		if(lags==T)
		{
			z <- cbind(z[at1,])
			id <- id[at1]
		}
		if(lags==F)
		{
			z <- cbind(z[at,])
			id <- id[at]
		}
	}
	else
	{
		at <- rep(T,NT.ini)
		at1 <- rep(T,NT.ini)

		z.full <- z
		id.full <- id
	}

	if(tdummies==T)
	{
		x <- cbind(x,td)
		z <- cbind(z,td[at,])
		z.full <- cbind(z.full,td)

		if(any(type==c("GMMcre","QMLcre")))
		{
			xa.names <- x.names
			za.names <- z.names
		}

		x.names <- c(x.names,td.names)
		z.names <- c(z.names,td.names)
	}
	else
	{
		xa.names <- x.names
		za.names <- z.names
	}

	if(type=="GMMpfe" & any(y==0))
	{
		x <- cbind(x,D)
		z <- cbind(z,D[at,])
		z.full <- cbind(z.full,D)
		x.names <- c(x.names,D.names)
		z.names <- c(z.names,D.names)
	}

	if(any(type==c("GMMcre","QMLcre")))
	{
		if(any(Ti.ini!=max(Ti.ini)))
		{
			Ti.unique <- unique(Ti.ini)
			Ti.unique <- Ti.unique[Ti.unique!=1]
			Ti.unique <- sort(Ti.unique)

			T.dum <- matrix(NA,nrow=NT.ini,ncol=length(Ti.unique))
			a <- 0

			Ti.ini.exp <- rep(Ti.ini,Ti.ini)

			for(j in Ti.unique)
			{
				a <- a+1
				T.dum[,a] <- Ti.ini.exp==j
			}

			if(x.exogenous==T)
			{
				x.m.aug <- T.dum*x.m[,1]
				if(ncol(x.m)>1) for(j in 2:ncol(x.m)) x.m.aug <- cbind(x.m.aug,T.dum*x.m[,j])
				x.m <- cbind(T.dum,x.m.aug)
				x.m.names <- c(paste("(Intercept)",Ti.unique,sep="_"),paste(rep(paste(xa.names,"mean",sep="_"),each=length(Ti.unique)),rep(Ti.unique,length(xa.names)),sep="_"))
				x.m <- cbind(x.m[Ti.ini.exp!=1,])
			}
			if(x.exogenous==F)
			{
				z.m.aug <- T.dum*z.m[,1]
				if(ncol(z.m)>1) for(j in 2:ncol(z.m)) z.m.aug <- cbind(z.m.aug,T.dum*z.m[,j])
				z.m <- cbind(T.dum,z.m.aug)
				z.m.names <- c(paste("(Intercept)",Ti.unique,sep="_"),paste(rep(paste(za.names,"mean",sep="_"),each=length(Ti.unique)),rep(Ti.unique,length(za.names)),sep="_"))
				z.m <- cbind(z.m[Ti.ini.exp!=1,])
			}

			id <- id[Ti.ini.exp!=1]
			y <- y[Ti.ini.exp!=1]
			x <- cbind(x[Ti.ini.exp!=1,])
			z <- cbind(z[Ti.ini.exp!=1,])
			at <- at[Ti.ini.exp!=1]
			if(type=="QMLcre" & x.exogenous==F) var.endog <- var.endog[Ti.ini.exp!=1]
			if(!is.null(offset)) offset <- offset[Ti.ini.exp!=1]
		}
		else
		{
			if(x.exogenous==T)
			{
				x.m <- cbind(1,x.m)
				x.m.names <- c("(Intercept)_mean",paste(xa.names,"mean",sep="_"))
			}
			if(x.exogenous==F)	
			{
				z.m <- cbind(1,z.m)
				z.m.names <- c("(Intercept)_mean",paste(za.names,"mean",sep="_"))
			}
		}

		if(x.exogenous==F)
		{
			x <- cbind(x,z.m)
			z <- cbind(z,z.m)
			x.names <- c(x.names,z.m.names)

			if(type=="QMLcre")
			{
				results <- lm(var.endog ~ 0+z)
				PIhat <- results$coefficients
				vhat <- results$residuals

				gixv <- t(z*vhat)

				x <- cbind(x,vhat)
				x1.names <- c(x.names,"vhat")
				x2.names <- paste("Z",c(z.names,z.m.names),sep="_")
				x.names <- c(x1.names,x2.names)
			}
		}
		else
		{
			x <- cbind(x,x.m)
			z <- x
			x.names <- c(x.names,x.m.names)

			gixv <- NA
			vhat <- NA
		}

		z.full <- z
	}
	else
	{
		gixv <- NA
		vhat <- NA
	}

	N <- length(unique(id))
	Ti <- as.vector(by(id,id,length))
	NT <- sum(at)

	k <- ncol(x)
	kz <- ncol(z)
	J <- NA
	dfJ <- kz-k

	### 4. Estimation

	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpre","GMMpfe","GMMcre"))) Hy <- fracregpd.links(link)$H1(y)
	if(type=="QMLcre") Hy <- y

	if(type=="GMMww" & GMMww.cor==T)
	{
		x <- x-matrix(apply(x,2,mean),nrow=nrow(x),ncol=k,byrow=T)
		z <- z-matrix(apply(z.full,2,mean),nrow=nrow(z),ncol=kz,byrow=T)
	}

	if(missing(start)) start <- rep(0,k)
	if(length(start)!=k) stop("start is not of the same dimension as the covariate vector (includes all auxiliary parameters)")

	results <- fracregpd.est(type,x.exogenous,lags,id,Ti,Hy,x,z,link,var.type,start,at,at1,variance,NT,k,kz,gixv,vhat,bootstrap,offset=offset,...)
	p <- results$p
	converged <- results$converged
	LL <- results$LL

	if(type!="QMLcre" & dfJ>0)
	{
		Qn <- results$Qn
		J <- NT*Qn
	}

	if(variance==T & converged==T)
	{
 		if(bootstrap==F) p.var <- results$p.var
		else
		{
			Ti.full <- as.vector(by(id.full,id.full,length))

			pboot <- matrix(NA,nrow=B,ncol=length(x.names))

			for(j in 1:B)
			{
				index.id <- sample(unique(id),N,replace=T)
				a <- 1
				aa <- 1

				for(jj in unique(id))
				{
					n.id <- sum(index.id==jj)
					if(n.id>=1)
					{
						bb.1 <- (a!=1)*sum(Ti.full[1:(a-1)])+1
						bb.2 <- sum(Ti.full[1:a])

						refa <- rep(bb.1:bb.2,n.id)
						refe <- rep(aa:(aa+n.id-1),each=Ti.full[a])

						if(aa==1)
						{
							ref <- refa
							id.B <- refe
						}
						else
						{
							ref <- c(ref,refa)
							id.B <- c(id.B,refe)
						}

						aa <- aa+n.id
					}

					a <- a+1
				}

				Hy.B <- Hy[ref]
				x.B <- cbind(x[ref,])
				z.B <- cbind(z.full[ref,])
				if(type=="QMLcre" & x.exogenous==F) var.endog.B <- var.endog[ref]
				at.B <- at[ref]
				at1.B <- at1[ref]
				if(!is.null(offset)) offset.B <- offset[ref] else offset.B <- NULL

				if(any(type==c("GMMc","GMMww")) | lags==T)
				{
					if(lags==T)
					{
						z.B <- cbind(z.B[at1.B,])
						id.B <- id.B[at1.B]
					}
					if(lags==F)
					{
						z.B <- cbind(z.B[at.B,])
						id.B <- id.B[at.B]
					}
				}

				Ti.B <- as.vector(by(id.B,id.B,length))
				NT.B <- sum(at.B)

				if(type=="QMLcre" & x.exogenous==F)
				{
					results <- lm(var.endog.B ~ 0+z.B)
					PIres <- results$coefficients
				}

				results <- fracregpd.est(type,x.exogenous,lags,id.B,Ti.B,Hy.B,x.B,z.B,link,var.type,start,at.B,at1.B,variance,NT.B,k,kz,NA,NA,bootstrap,offset=offset.B,...)
				if(results$converged==T)
				{
					if(table==T) cat("1")

					pres <- results$p
					if(type=="QMLcre" & x.exogenous==F) pres <- c(pres,PIres)
					pboot[j,] <- pres
				}
				else if(table==T) cat("0")

				if(any(j==seq(50,100000,50)) & table==T) cat("\n"
)			}

			p.var <- matrix(NA,nrow=length(x.names),ncol=length(x.names))
			diag(p.var) <- apply(pboot,2,var,na.rm=T)
		}
	}
	else p.var <- "singular"

	if(type=="QMLcre" & x.exogenous==F) p <- c(p,PIhat)

	table.info <- list(p=p,p.var=p.var,x.names=x.names,x.exogenous=x.exogenous,lags=lags,type=type,link=link,converged=converged,N.ini=N.ini,N=N,NT.ini=NT.ini,NT=NT,J=J,dfJ=dfJ,k=k,var.type=var.type,bootstrap=bootstrap,LL=LL,or=or,level=level)
	if(table==T) do.call(fracregpd.table, table.info)

	names(p) <- x.names
	res <- list(call=cl, type=type,link=link,Hy=Hy,p=p,converged=converged,table.info=table.info,x=x,xbhat=results$xbhat)
	if(!is.null(offset)) res[["offset"]] <- offset
	if(dfJ>0) res[["J"]] <- J

	if(variance==T & converged==T)
	{ 
		if(is.character(p.var)) p.var <- matrix(NA,nrow=length(p),ncol=length(p))
		dimnames(p.var) <- list(x.names,x.names)
		res[["p.var"]] <- p.var
		res[["var.type"]] <- var.type
	}

	### 5. Return results

	class(res) <- "fracregpd"
	return(invisible(res))
}


