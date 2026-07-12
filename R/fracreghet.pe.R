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

