fracreg.pe <- function(object,APE=TRUE,CPE=FALSE,at=NULL,which.x=NULL,variance=TRUE,table=TRUE)
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
		if(is.null(object$p.var)) stop("fracreg command was run with variance = F")
	}
	else if (object$type=="2P")
	{
		if(is.null(object$resBIN$p.var)) stop("fracreg command was run with variance = F")
	}
	else if (object$type=="3P")
	{
		if(is.null(object$resBIN0$p.var)) stop("fracreg command was run with variance = F")
	}

	if(table==TRUE & variance==FALSE)
	{
		variance <- TRUE
		warning("option variance changed from F to TRUE, as required by table=T")
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

	if(any(x.names=="INTERCEPT")) xvar.names <- x.names[-1]
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

		if(any(x.names=="INTERCEPT")) p.pe <- matrix(rep(pa[-1],each=n),ncol=k)
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
			if(any(x.names=="INTERCEPT")) p.pe2 <- matrix(rep(pb[-1],each=n),ncol=k) else p.pe2 <- matrix(rep(pb,each=n),ncol=k)
			if(any(x.names=="INTERCEPT")) p.pe3 <- matrix(rep(pc[-1],each=n),ncol=k) else p.pe3 <- matrix(rep(pc,each=n),ncol=k)
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

			if(any(x.names=="INTERCEPT")) p.pe <- matrix(rep(pb[-1],each=n),ncol=k)
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

		if(table==TRUE) fracreg.pe.table(PE.p,PE.sd,PE.type,which.x,xvar.names,title,at)
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
					if(any(x.names=="INTERCEPT")) xm <- c(1,at)
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
				if(any(x.names=="INTERCEPT")) xm <- c(1,at)
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

		if(table==TRUE) fracreg.pe.table(PE.p,PE.sd,PE.type,which.x,xvar.names,title,at)
	}

	### 5. Return results

	if(APE==TRUE & CPE==TRUE) return(invisible(list(ape=resAPE,cpe=resCPE)))
	if(APE==TRUE & CPE==FALSE) return(invisible(resAPE))
	if(APE==FALSE & CPE==TRUE) return(invisible(resCPE))
}
