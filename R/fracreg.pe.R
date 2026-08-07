#' @title Fractional Response Regressions - Partial Effects
#'
#' @description
#' \code{fracreg.pe} is used to compute average and/or conditional partial effects in fractional response models.
#' @param object an object containing the results of an \code{fracreg} command.
#' @param APE a logical value indicating whether average partial effects are to be computed.
#' @param CPE a logical value indicating whether conditional partial effects are to be computed.
#' @param at a numeric vector containing the covariates' values at which the conditional partial effects are to be computed or  the strings \code{"mean"} (the default) or \code{"median"}, in which cases the covariates are evaluated at their  mean or median values (or mode, in case of dummy variables), respectively.
#' @param which.x a vector containing the names of the covariates to which the partial effects are to be computed.
#' @param variance a logical value indicating whether the variance of the estimated partial effects should be calculated. Defaults to  \code{TRUE} whenever \code{table = TRUE}.
#' @param table a logical value indicating whether a summary table with the results should be printed.
#'
#' @details
#' \code{fracreg.pe} calculates partial effects for fractional response models estimated via \code{fracreg}. \code{fracreg.pe} may be used to compute average or conditional partial effects for: (i) one-part fractional response models; (ii) the binary components of two-part and three-part fractional response models; (iii) the fractional components of two-part and three-part fractional response models; and (iv) two-part and three-part fractional response models overall. 
#' 
#' \strong{Partial Effects for Continuous Variables:}
#' For a continuous covariate \eqn{x_k}, the partial effect on the conditional mean \eqn{E(y|x) = G(x\beta)} is the first derivative with respect to \eqn{x_k}:
#' \deqn{PE_k(x) = \frac{\partial E(y|x)}{\partial x_k} = g(x\beta)\beta_k}
#' where \eqn{g(\cdot)} is the probability density function corresponding to the link function \eqn{G(\cdot)}.
#' 
#' \strong{Partial Effects for Discrete Variables:}
#' For a discrete or dummy covariate \eqn{x_k}, the partial effect is calculated as the discrete difference in the expected value when \eqn{x_k} changes from 0 to 1, holding all other variables \eqn{x_{-k}} constant:
#' \deqn{PE_k(x) = G(x_{-k}\beta_{-k} + \beta_k) - G(x_{-k}\beta_{-k})}
#' 
#' \strong{Average vs. Conditional Partial Effects:}
#' - \strong{Average Partial Effects (APE):} Evaluated for each observation \eqn{i} in the sample and then averaged:
#' \deqn{APE_k = \frac{1}{N} \sum_{i=1}^N PE_k(x_i)}
#' - \strong{Conditional Partial Effects (CPE):} Evaluated at a specific vector of covariate values \eqn{x^*} (e.g., the sample mean or median):
#' \deqn{CPE_k = PE_k(x^*)}
#' 
#' For calculating standard errors, it is taken into account the option that was previously chosen for estimating the model. See Ramalho, Ramalho and Murteira (2011) and Fang and Ma (2013) for details on the computation of partial effects in the fractional response framework.
#'
#' @return
#' \code{fracreg.pe} returns a list with the following element:
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
#' Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
#' estimating and testing empirical strategies for fractional response models",
#' \emph{Journal of Economic Surveys}, 25(1), 19-68.
#' 
#' Fang, K., & Ma, S. (2013), "Three-part model for fractional response variables with application to Chinese household health insurance coverage", \emph{Journal of Applied Statistics}, 40(5), 925-940.
#'
#' @author Sulman Olieko Owili <oliekosulman@gmail.com>
#'
#' @seealso
#' \code{\link{fracreg}}, for fitting fractional response models.\cr
#' \code{\link{fracreg.reset}} and \code{\link{fracreg.ggoff}}, for specification tests.\cr
#' \code{\link{fracreg.ptest}}, for non-nested hypothesis tests.
#'
#' @examples
#' ### Empirical 401(k) Examples
#' data("fracreg_k401k")
#' y <- fracreg_k401k$prate
#' X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
#'            totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)
#' 
#' m <- fracreg(y, X, type="1P", linkfrac="logit")
#' pe_res <- fracreg.pe(m)
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
#' ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
#' y <- rbeta(N,ym*20,20*(1-ym))
#' y[y > 0.9] <- 1
#' 
#' #Computing average partial effects for a logit fractional response model
#' mod <- fracreg(y,X,linkfrac="logit")
#' pe_res <- fracreg.pe(mod)
#' summary(pe_res)
#' 
#' #Computing average partial effects for a binary logit + fractional probit
#' #two-part model
#' mod <- fracreg(y,X,linkbin="logit",linkfrac="probit",type="2P",inf=1)
#' pe_res <- fracreg.pe(mod)
#' summary(pe_res)
#' 
#' #Computing conditional partial effects for X2 in the logit component
#' #of a two-part fractional response model, with the covariates evaluated
#' #at their median values
#' mod <- fracreg(y,X,linkfrac="logit",type="2Pfrac",inf=1)
#' pe_res <- fracreg.pe(mod,APE=FALSE,CPE=TRUE,at="median",which.x="X2")
#' summary(pe_res)
#' 
#' #Computing average partial effects for a three-part double-inflated model
#' y3p <- y
#' y3p[1:20] <- 0
#' y3p[21:40] <- 1
#' res3p <- fracreg(y3p,X,linkbin=c("logit","probit"),linkfrac="logit",type="3P")
#' pe_res <- fracreg.pe(res3p)
#' summary(pe_res)
#' @export
fracreg.pe <- function(object,APE=TRUE,CPE=FALSE,at=NULL,which.x=NULL,variance=TRUE,table=FALSE)
{
	### 1. Error and warning messages

	if(missing(object)) stop("object is missing")
	if(is.null(object$class)) stop("object is not the output of an fracreg command")
	if(object$class!="fracreg") stop("object is not the output of an fracreg command")

	if(!is.logical(APE)) stop("non-logical value assigned to option APE")
	if(!is.logical(CPE)) stop("non-logical value assigned to option CPE")
	if(!is.logical(variance)) stop("non-logical value assigned to option variance")
	if(!is.logical(table)) stop("non-logical value assigned to option table")

	if(all(c(APE,CPE)==FALSE)) stop("You must specify at least one option: APE and/or CPE")
	if(CPE==FALSE & !is.null(at)) stop("option at is only required for CPE")

	if(object$converged==0) stop("object is not the output of a successful (converged) fracreg command")

	if(object$type!="2P" & object$type!="3P")
	{
		if(is.null(object$p.var)) stop("fracreg command was run with variance = FALSE")
	}
	else if (object$type=="2P")
	{
		if(is.null(object$resBIN$p.var)) stop("fracreg command was run with variance = FALSE")
	}
	else if (object$type=="3P")
	{
		if(is.null(object$resBIN0$p.var)) stop("fracreg command was run with variance = FALSE")
	}

	if(table==TRUE & variance==FALSE)
	{
		variance <- TRUE
		warning("option variance changed from FALSE to TRUE, as required by table=TRUE")
	}

	### 2. Recovering definitions and estimates and other definitions

	type <- object$type

	if(type!="2P" & type!="3P")
	{
		linka <- object$link
		pa <- object$p
		pa.var <- object$p.var
		x <- model.matrix(object$formula)
		x.names <- object$x.names

		if(type=="1P") title <- paste("Fractional",linka,"regression")
		if(type=="2Pbin") title <- paste("Binary",linka,"component of a two-part regression")
		if(type=="2Pfrac") title <- paste("Fractional",linka,"component of a two-part regression")
		if(type=="3Pbin0") title <- paste("Binary",linka,"component 1 of a three-part regression")
		if(type=="3Pbin1") title <- paste("Binary",linka,"component 2 of a three-part regression")
		if(type=="3Pfrac") title <- paste("Fractional",linka,"component of a three-part regression")
	}
	if(type=="3P")
	{
		linka <- object$resBIN0$link
		pa <- object$resBIN0$p
		pa.var <- object$resBIN0$p.var
		xa <- model.matrix(object$resBIN0$formula)
		xa.names <- object$resBIN0$x.names

		linkb <- object$resBIN1$link
		pb <- object$resBIN1$p
		pb.var <- object$resBIN1$p.var

		linkc <- object$resFRAC$link
		pc <- object$resFRAC$p
		pc.var <- object$resFRAC$p.var
		xc <- object$x2base
		xc.names <- object$resFRAC$x.names

		if(!identical(xa.names,xc.names)) stop("currently fracreg.pe requires all components of 3P models to use the same covariates")
		x <- xa
		x.names <- xa.names

		title <- paste("Three-part regression - binary",linka,", binary",linkb,"+ fractional",linkc)
	}
	if(type=="2P")
	{
		linka <- object$resBIN$link
		pa <- object$resBIN$p
		pa.var <- object$resBIN$p.var
		xa <- model.matrix(object$resBIN$formula)
		xa.names <- object$resBIN$x.names

		linkb <- object$resFRAC$link
		pb <- object$resFRAC$p
		pb.var <- object$resFRAC$p.var
		xb <- object$x2base
		xb.names <- object$resFRAC$x.names

		if(!identical(xa.names,xb.names)) stop("currently fracreg.pe requires both components of two-part models to use the same covariates")
		x <- xa
		x.names <- xa.names

		title <- paste("Binary",linka,"+ Fractional",linkb,"two-part regression")
	}

	if(any(x.names=="(Intercept)")) xvar.names <- x.names[-1]
	else xvar.names <- x.names

	k <- length(xvar.names)
	npar <- ncol(x)
	n <- nrow(x)

	if(is.null(which.x)) which.x <- xvar.names
	xw.names <- unique(c(xvar.names,which.x))
	if(!identical(xvar.names,xw.names)) stop("option which not appropriately defined")

	### 3. Average partial effects

	if(APE==TRUE)
	{
		PE.type <- "APE"

		if(any(x.names=="(Intercept)")) p.pe <- matrix(rep(pa[-1],each=n),ncol=k)
 		else p.pe <- matrix(rep(pa,each=n),ncol=k)
		dimnames(p.pe) <- list(NULL,xvar.names)

		if(type!="2P" & type!="3P") xbhata <- object$xbhat
		if(type=="2P") xbhata <- object$resBIN$xbhat
		if(type=="3P") xbhata <- object$resBIN0$xbhat

		ga <- fracreg.links(linka)$mu.eta(xbhata)
		PEa.p <- as.matrix(p.pe[,which.x])*ga

		if(type!="2P" & type!="3P") PE.p <- apply(PEa.p,2,mean)

		if(type=="3P")
		{
			yhata <- object$resBIN0$yhat
			if(any(x.names=="(Intercept)")) p.pe2 <- matrix(rep(pb[-1],each=n),ncol=k) else p.pe2 <- matrix(rep(pb,each=n),ncol=k)
			if(any(x.names=="(Intercept)")) p.pe3 <- matrix(rep(pc[-1],each=n),ncol=k) else p.pe3 <- matrix(rep(pc,each=n),ncol=k)
			dimnames(p.pe2) <- list(NULL,xvar.names)
			dimnames(p.pe3) <- list(NULL,xvar.names)
			
			xbhatb <- as.vector(x%*%pb)
			gb <- fracreg.links(linkb)$mu.eta(xbhatb)
			yhatb <- fracreg.links(linkb)$linkinv(xbhatb)
			PEb.p <- as.matrix(p.pe2[,which.x])*gb
			
			xbhatc <- as.vector(x%*%pc)
			gc <- fracreg.links(linkc)$mu.eta(xbhatc)
			yhatc <- fracreg.links(linkc)$linkinv(xbhatc)
			PEc.p <- as.matrix(p.pe3[,which.x])*gc
			
			H <- yhatb + (1-yhatb)*yhatc
			PE.p <- apply( PEa.p*H + PEb.p*yhata*(1-yhatc) + PEc.p*yhata*(1-yhatb), 2, mean )
		}
		else if(type=="2P")
		{
			yhata <- object$resBIN$yhat
			PEa.p <- as.matrix(p.pe[,which.x])*ga

			if(any(x.names=="(Intercept)")) p.pe <- matrix(rep(pb[-1],each=n),ncol=k)
 			else p.pe <- matrix(rep(pb,each=n),ncol=k)
			dimnames(p.pe) <- list(NULL,xvar.names)

			xbhatb <- as.vector(x%*%pb)
			gb <- fracreg.links(linkb)$mu.eta(xbhatb)
			yhatb <- fracreg.links(linkb)$linkinv(xbhatb)
			PEb.p <- as.matrix(p.pe[,which.x])*gb

			PE.p <- apply(PEb.p*yhata+PEa.p*yhatb,2,mean)
		}

		resAPE <- list(PE.p=PE.p)

		if(variance==TRUE)
		{
			if(type!="2P" & type!="3P") PE.sd <- fracreg.pe.var(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var)
			if(type=="2P") PE.sd <- fracreg.pe.var(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var,pb,xbhatb,gb,linkb,pb.var,yhata,yhatb)
			if(type=="3P") PE.sd <- fracreg.pe.var(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var,pb,xbhatb,gb,linkb,pb.var,yhata,yhatb,pc,xbhatc,gc,linkc,pc.var,yhatc)
			resAPE[["PE.sd"]] <- PE.sd
		}

		table.info.APE <- list(PE.p=PE.p,PE.sd=if(variance) PE.sd else NA,PE.type=PE.type,which.x=which.x,xvar.names=xvar.names,title=title,at=at)
		if(table==TRUE) do.call(fracreg.pe.table, table.info.APE)
		resAPE[["table.info"]] <- table.info.APE
	}

	### 4. Conditional partial effects

	if(CPE==TRUE)
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
				xm[xdum==TRUE] <- round(xm,0)[xdum==TRUE]
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

		xbhata <- as.vector(xm%*%pa)
		ga <- fracreg.links(linka)$mu.eta(xbhata)
		PEa.p <- pa[which.x]*ga
		PE.p <- PEa.p

		if(type=="3P")
		{
			yhata <- fracreg.links(linka)$linkinv(xbhata)
			
			xbhatb <- as.vector(xm%*%pb)
			gb <- fracreg.links(linkb)$mu.eta(xbhatb)
			yhatb <- fracreg.links(linkb)$linkinv(xbhatb)
			PEb.p <- pb[which.x]*gb
			
			xbhatc <- as.vector(xm%*%pc)
			gc <- fracreg.links(linkc)$mu.eta(xbhatc)
			yhatc <- fracreg.links(linkc)$linkinv(xbhatc)
			PEc.p <- pc[which.x]*gc
			
			H <- yhatb + (1-yhatb)*yhatc
			PE.p <- PEa.p*H + PEb.p*yhata*(1-yhatc) + PEc.p*yhata*(1-yhatb)
		}
		else if(type=="2P")
		{
			yhata <- fracreg.links(linka)$linkinv(xbhata)

			xbhatb <- as.vector(xm%*%pb)
			gb <- fracreg.links(linkb)$mu.eta(xbhatb)
			yhatb <- fracreg.links(linkb)$linkinv(xbhatb)
			PEb.p <- pb[which.x]*gb

			PE.p <- PEb.p*yhata+PEa.p*yhatb
		}

		resCPE <- list(PE.p=PE.p)

		if(variance==TRUE)
		{
			if(type!="2P" & type!="3P") PE.sd <- fracreg.pe.var(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var)
			if(type=="2P") PE.sd <- fracreg.pe.var(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var,pb,xbhatb,gb,linkb,pb.var,yhata,yhatb)
			if(type=="3P") PE.sd <- fracreg.pe.var(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var,pb,xbhatb,gb,linkb,pb.var,yhata,yhatb,pc,xbhatc,gc,linkc,pc.var,yhatc)
			resCPE[["PE.sd"]] <- PE.sd
		}

		table.info.CPE <- list(PE.p=PE.p,PE.sd=if(variance) PE.sd else NA,PE.type=PE.type,which.x=which.x,xvar.names=xvar.names,title=title,at=at)
		if(table==TRUE) do.call(fracreg.pe.table, table.info.CPE)
		resCPE[["table.info"]] <- table.info.CPE
	}

	### 5. Return results

	if(APE==TRUE & CPE==TRUE) res <- list(ape=resAPE,cpe=resCPE)
	else if(APE==TRUE & CPE==FALSE) res <- resAPE
	else if(APE==FALSE & CPE==TRUE) res <- resCPE
	
	res[["link"]] <- object$link
	res[["inf"]] <- object$inf
	
	if(type=="2P")
	{
		res[["resBIN"]] <- fracreg.pe(object$resBIN, table=FALSE, APE=APE, CPE=CPE, at=at, which.x=which.x, variance=variance)
		res[["resFRAC"]] <- fracreg.pe(object$resFRAC, table=FALSE, APE=APE, CPE=CPE, at=at, which.x=which.x, variance=variance)
	}
	if(type=="3P")
	{
		res[["resBIN0"]] <- fracreg.pe(object$resBIN0, table=FALSE, APE=APE, CPE=CPE, at=at, which.x=which.x, variance=variance)
		res[["resBIN1"]] <- fracreg.pe(object$resBIN1, table=FALSE, APE=APE, CPE=CPE, at=at, which.x=which.x, variance=variance)
		res[["resFRAC"]] <- fracreg.pe(object$resFRAC, table=FALSE, APE=APE, CPE=CPE, at=at, which.x=which.x, variance=variance)
	}
	
	class(res) <- "fracreg.pe"
	return(invisible(res))
}
