#' @title Fractional Response Regressions under Unobserved Heterogeneity - Partial Effects
#'
#' @description
#' \code{fracreghet.pe} is used to compute average and/or conditional partial effects in fractional response models 
#' under unobserved heterogeneity.
#' @param object an object containing the results of an \code{fracreghet} command.
#' @param smearing a logical value indicating whether the smearing correction is to be applied
#' @param APE a logical value indicating whether average partial effects are to be computed.
#' @param CPE a logical value indicating whether conditional partial effects are to be computed.
#' @param at a numeric vector containing the covariates' values at which the conditional partial effects are to be computed or the strings \code{"mean"} (the default) or \code{"median"}, in which cases the covariates are evaluated at their mean or median values (or mode, in case of dummy variables), respectively.
#' @param which.x a vector containing the names of the covariates to which the partial effects are to be computed.
#' @param table a logical value indicating whether a summary table with the results should be printed.
#' @param variance a logical value indicating whether the variance of the estimated partial effects should be calculated. Defaults to \code{TRUE} whenever \code{table = TRUE}.
#'
#' @details
#' \code{fracreghet.pe} calculates partial effects for fractional response models estimated via \code{fracreghet}. \code{fracreghet.pe} may be used to compute average or conditional partial effects. These partial effects may be conditional only on observables, using the smearing estimator, or also on unobservables, setting the error term to zero.
#' 
#' \strong{Partial Effects under Unobserved Heterogeneity:}
#' When unobserved heterogeneity or endogeneity is present, calculating partial effects requires dealing with the unobserved error \eqn{v_i}. Let the conditional mean be \eqn{E(y|x, v) = G(x\beta + \gamma v)}.
#' - \strong{Conditional on Observables (Smearing):} The unobserved heterogeneity is integrated out over its empirical distribution. The average partial effect for a continuous variable \eqn{x_k} is computed as:
#' \deqn{PE_k(x) = \frac{1}{N} \sum_{i=1}^N g(x\beta + \gamma \hat{v}_i) \beta_k}
#' - \strong{Conditional on Unobservables (Error = 0):} The partial effect is evaluated for an individual with the mean level of unobserved heterogeneity (\eqn{v = 0}):
#' \deqn{PE_k(x) = g(x\beta) \beta_k}
#' 
#' For discrete variables, the partial effects are calculated as the discrete differences evaluated using either the smearing approach or setting the error term to zero. 
#' 
#' For calculating standard errors, it is taken into account the option that was previously chosen for estimating the model. See Ramalho and Ramalho (2017) for details on the computation of partial effects for fractional response models under unobserved heterogeneity.
#'
#' @return
#' \code{fracreghet.pe} returns a list with the following element:
#'   \item{PE.p}{a named vector of partial effects.
#' }
#' 
#' If \code{variance = TRUE} or \code{table = TRUE}, the previous list also contains the following element:
#'   \item{PE.sd}{a named vector of standard errors of the estimated partial effects.
#' }
#' 
#' When both average and conditional partial effects are requested, two lists containing the previous elements are returned, indexed by the prefixes \code{ape} and \code{cpe}.
#'
#' @references
#' Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of nonlinear regression models with boundary outcomes and endogeneity, with applications to nonnegative and fractional responses", \emph{Econometric Reviews}, 36(4), 397-420.
#'
#' @author Sulman Olieko Owili <oliekosulman@gmail.com>
#'
#' @seealso
#' \code{\link{fracreghet}}, for fitting fractional response models under unobserved heterogeneity.\cr
#' \code{\link{fracreghet.reset}}, for the RESET test.\cr
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
#' res_emp <- fracreghet(y_adj, X_het, Z_emp, var.endog = X_het[, "mrate"], 
#'                       type="QMLxv", link="logit") 
#' pe_res <- fracreghet.pe(res_emp, which.x="mrate")
#' summary(pe_res)
#'  
#' ### Simulated Examples
#' 
#' N <- 250
#' u <- rnorm(N)
#' 
#' X <- cbind(rnorm(N),rnorm(N))
#' dimnames(X)[[2]] <- c("X1","X2")
#' 
#' Z <- cbind(rnorm(N),rnorm(N),rnorm(N))
#' dimnames(Z)[[2]] <- c("Z1","Z2","Z3")
#' 
#' y <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
#' 
#' mod <- fracreghet(y,X,type="GMMx")
#' 
#' #Smearing estimator of average partial effects for variable X1
#' pe_res <- fracreghet.pe(mod,which.x="X1")
#' summary(pe_res)
#' 
#' #Naive estimator of conditional partial effects for all covariates,
#' #which are evaluated at X1=1 and X2=-1
#' pe_res <- fracreghet.pe(mod,smearing=FALSE,APE=FALSE,CPE=TRUE,at=c(1,-1))
#' summary(pe_res)
#' @export
fracreghet.pe <- function(object,smearing=T,APE=T,CPE=F,at=NULL,which.x=NULL,table=FALSE,variance=T)
{
	### 1. Error and warning messages

	if(missing(object)) stop("object is missing")
	if(is.null(object$class)) stop("object is not the output of an fracreghet command")
	if(object$class!="fracreghet") stop("object is not the output of an fracreghet command")

	if(all(c(APE,CPE)==F)) stop("You must specify at least one option: APE and/or CPE")
	if(CPE==F & !is.null(at)) stop("option at is only required for cpe")

	if(object$converged==0) stop("object is not the output of a successful (converged) fracreghet command")

	if(table==T & variance==F)
	{
		variance <- T
		warning("option variance changed from F to T, as required by table=T")
	}

	### 2. Recovering definitions and estimates and other definitions

	type <- object$type
	link <- object$link
	adjust <- object$adjust
	p <- object$p
	p.var <- object$p.var
	x <- model.matrix(object$formula)
	x.names <- object$x.names
	Hy <- object$Hy

	if(any(x.names=="(Intercept)")) xvar.names <- x.names[-1]
	else xvar.names <- x.names

	k <- length(xvar.names)
	npar <- ncol(x)
	n <- nrow(x)

	if(any(type==c("QMLxv"))) pv <- p[npar]
	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		npar <- npar-1
		p <- p[1:npar]
		p.var <- p.var[1:npar,1:npar]
	}

	if(all(type!=c("GMMxv","LINxv","QMLxv"))) xbhat <- object$xbhat
	if(any(type==c("GMMxv","LINxv","QMLxv"))) xbhat <- as.vector(x[,-(npar+1)]%*%as.matrix(p))

	if(is.null(which.x)) which.x <- xvar.names
	xw.names <- unique(c(xvar.names,which.x))
	if(!identical(xvar.names,xw.names)) stop("option which not appropriately defined")

	### 3. Average partial effects

	if(smearing==T) title1 <- "(conditional only on observables, based on the smearing estimator)"
	else title1 <- "(conditional on both observables and unobservables, with error term = 0)"

	title2 <- paste("Fractional",link,"regression")
	title3 <- paste("Estimator:",type)
	title <- c(title1,title2,title3)
	if(adjust!=0)
	{
		if(is.numeric(adjust)) title4 <- paste("Adjustment:",adjust,"added to all observations")
		else title4 <- "Adjustment: all boundary observations dropped"
		title <- c(title,title4)
	}

	if(APE==T)
	{
		PE.type <- "APE"

		if(any(x.names=="(Intercept)")) p.pe <- matrix(rep(p[-1],each=n),ncol=k)
 		else p.pe <- matrix(rep(p,each=n),ncol=k)
		dimnames(p.pe) <- list(NULL,xvar.names)

		if(smearing==F)
		{
			g <- fracreghet.links(link)$mu.eta(xbhat)
			what <- NA
		}
		if(smearing==T)
		{
			if(type=="QMLxv") what <- x[,npar+1]*pv

			if(any(type==c("LINx","LINxv","LINz"))) what <- Hy-xbhat
			if(any(type==c("GMMx","GMMxv","GMMz")))
			{
				G2 <- fracreghet.links(link)$G2(xbhat)
				uhat.star <- Hy/G2-1
				what <- fracreghet.links(link)$H2(uhat.star+1)
			}

			g <- rep(NA,n)
			for(j in 1:n) g[j] <- mean(fracreghet.links(link)$mu.eta(xbhat[j]+what))
		}

		PE.p <- as.matrix(p.pe[,which.x])*g
		PE.p <- apply(PE.p,2,mean)

		resAPE <- list(PE.p=PE.p)

		if(variance==T)
		{
			PE.sd <- fracreghet.pe.var(x,npar,which.x,x.names,xvar.names,type,p,xbhat,g,link,p.var,smearing,what,n)
			resAPE[["PE.sd"]] <- PE.sd
		}

		table.info.APE <- list(PE.p=PE.p,PE.sd=if(variance) PE.sd else NA,PE.type=PE.type,which.x=which.x,xvar.names=xvar.names,title=title,adjust=adjust,at=at)
		if(table==T) do.call(fracreghet.pe.table, table.info.APE)
		resAPE[["table.info"]] <- table.info.APE
	}

	### 4. Conditional partial effects

	if(CPE==T)
	{
		PE.type <- "CPE"

		if(is.null(at)) at <- "mean"

		if(length(at)==1)
		{
			if((!any(at==c("mean","median")) & k!=1)) stop("at not appropriately specified")

			if(any(at==c("mean","median")))
			{
				if(at=="mean") xm <- apply(x,2,mean)
				if(at=="median") xm <- apply(x,2,median)

				xdum <- apply(x,2,function(a) all(a %in% c(0,1)))

				if(any(type==c("GMMxv","LINxv","QMLxv")))
				{
					xm <- xm[-npar]
					xdum <- xdum[-npar]
				}
				xm[xdum==T] <- round(xm,0)[xdum==T]
			}
			else
			{
				if(is.numeric(at))
				{
					if(any(x.names=="(Intercept)")) xm <- c(1,at)
					else xm <- at
				}
				else stop("at not appropriately specified")
			}
		}
		else
		{
			if(length(at)!=k) stop("at not appropriately specified")
			else
			{
				if(any(x.names=="(Intercept)")) xm <- c(1,at)
				else xm <- at
			}
		}

		xmbhat <- as.vector(xm%*%p)

		if(smearing==F) g <- fracreghet.links(link)$mu.eta(xmbhat)

		if(smearing==T)
		{
			if(type=="QMLxv") what <- x[,npar]*pv
			if(any(type==c("LINx","LINxv","LINz"))) what <- Hy-xbhat
			if(any(type==c("GMMx","GMMxv","GMMz")))
			{
				G2 <- fracreghet.links(link)$G2(xbhat)
				uhat.star <- Hy/G2-1
				what <- fracreghet.links(link)$H2(uhat.star+1)
			}

			g <- mean(fracreghet.links(link)$mu.eta(xmbhat+what))
		}

		PE.p <- p[which.x]*g

		resCPE <- list(PE.p=PE.p)

		if(variance==T)
		{
			PE.sd <- fracreghet.pe.var(x,npar,which.x,x.names,xvar.names,type,p,xbhat,g,link,p.var,smearing,what,n)
			resCPE[["PE.sd"]] <- PE.sd
		}

		table.info.CPE <- list(PE.p=PE.p,PE.sd=if(variance) PE.sd else NA,PE.type=PE.type,which.x=which.x,xvar.names=xvar.names,title=title,adjust=adjust,at=at)
		if(table==T) do.call(fracreghet.pe.table, table.info.CPE)
		resCPE[["table.info"]] <- table.info.CPE
	}

	### 5. Return results

	if(APE==T & CPE==T) res <- list(ape=resAPE,cpe=resCPE)
	else if(APE==T & CPE==F) res <- resAPE
	else if(APE==F & CPE==T) res <- resCPE
	
	class(res) <- "fracreghet.pe"
	return(invisible(res))
}

fracreghet.pe.var <- function(x,npar,which.x,x.names,xvar.names,type,p,xbhat,g,link,p.var,smearing,what,N)
{
	if(smearing==F) gd <- fracreghet.links(link)$gd(xbhat)
	if(smearing==T) 
	{
		gd <- rep(NA,N)
		for(j in 1:N) gd[j] <- mean(fracreghet.links(link)$gd(xbhat[j]+what))
	}

	PE.sd <- matrix(NA,nrow=npar,ncol=npar)

	for(i in 1:npar)
	{
		for(j in 1:npar) PE.sd[i,j] <- mean(p[i]*gd*x[,j]+(i==j)*g)
	}

	PE.sd <- PE.sd%*%p.var%*%t(PE.sd)

	PE.sd <- diag(PE.sd)^0.5
	if(any(x.names=="(Intercept)")) PE.sd <- PE.sd[-1]

	names(PE.sd) <- xvar.names
	PE.sd <- PE.sd[which.x]

	return(PE.sd)
}

