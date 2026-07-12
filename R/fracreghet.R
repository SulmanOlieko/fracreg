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

	GMM.est <- T
	if(any(type==c("GMMx","GMMxv")) & !any(Hy==0))
	{
		results <- tryCatch(glm(Hy ~ x-1,family=Gamma(link=log),maxit=100,offset=offset),error=function(e) return(NULL))
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
			fi.inv <- fracreghet.var(type,p,XB,x,z,link,Hy,var.type,var.cluster,T,gixv,vhat)$fi.inv
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

		converged <- ifelse(results$convergence==0,T,F)
	}

	ret.list <- list(p=p,XB=XB,converged=converged)
	if(type=="GMMz" & ncol(z)>ncol(x)) ret.list[["Qn"]] <- Qn
	if (any(type == c("QMLxv", "QMLz"))) {
		yhat <- fracreghet.links(link)$linkinv(XB)
		eps <- 1e-16
		LL <- sum(ifelse(Hy > 0, Hy * log(pmax(yhat, eps)), 0) + ifelse(Hy < 1, (1-Hy) * log(pmax(1-yhat, eps)), 0))
		ret.list[["LL"]] <- LL
	}

	if(variance==F | converged==F) return(ret.list)

	if(any(type==c("GMMxv","QMLxv"))) z <- z.in
	p.var <- fracreghet.var(type,p,XB,x,z,link,Hy,var.type,var.cluster,F,gixv,vhat)$p.var
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

	if(step.one==T) return(list(fi.inv=fi.inv))

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

fracreghet <- function(y,x,z=x,var.endog,start,type="GMMx",link="logit",intercept=T,table=T,variance=T,var.type="robust",var.cluster,adjust=0,offset=NULL,or=FALSE,level=0.95,...)
{
	cl <- match.call()
	LL <- NULL

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
	if(table==T & variance==F)
	{
		variance <- T
		warning("option variance changed from F to T, as required by table=T")
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

	if(intercept==T)
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

		converged <- T
		if(variance==T) results <- fracreghet.var(type,p,XB,x,z,link,Hy,var.type,var.cluster,F,gixv,vhat)
	}

	if(variance==T & converged==T) p.var <- results$p.var
	else p.var <- "singular"

	x.names.in <- x.names

	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		p <- c(p,PIhat)
		x.names <- c(x.names,"vhat",paste("Z",z.names,sep="_"))
	}
	names(p) <- x.names

	if(table==T) fracreghet.table(p,p.var,x.names,type,link,converged,N,var.type,adjust,k,J,dfJ,LL=LL,or=or,level=level)

	formula <- y ~ x - 1

	res <- list(class=class,formula=formula,type=type,link=link,adjust=adjust,p=p,Hy=Hy,xbhat=as.vector(XB),converged=converged,x.names=x.names.in)
	if(!is.null(offset)) res[["offset"]] <- offset
	if(any(type==c("GMMz","LINz")) & kz>k) res[["J"]] <- J

	if(variance==T & converged==T)
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
