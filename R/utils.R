fracreg.links <- function(link) 
{
	switch(link,
		logit = {
			linkfun <- function(mu) qlogis(mu)
			linkinv <- function(eta) plogis(eta)
			mu.eta <- function(eta) dlogis(eta)
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
			gd <- function(eta) (exp(-exp(eta))*exp(eta))*(1-exp(eta))
			valideta <- function(eta) TRUE
		},

		loglog = {
			linkfun <- function(mu) -log(-log(mu))
			linkinv <- function(eta) exp(-exp(-eta))
			mu.eta <- function(eta) exp(-exp(-eta)-eta)
			gd <- function(eta) (exp(-exp(-eta))*exp(-eta))*(exp(-eta)-1)
			valideta <- function(eta) TRUE
		},

		stop(sQuote(link), " - link not recognised")
	)

	structure(list(linkfun=linkfun,linkinv=linkinv,mu.eta=mu.eta,gd=gd,valideta=valideta,name=link),class="link-glm")
}

fracreg.est <- function(y,x,link,method,variance,var.type,var.eim,var.cluster,dfc,...)
{
	if(method=="ML") results <- glm(y ~ x-1,family=binomial(link=fracreg.links(link)),...)
	if(method=="QML") results <- glm(y ~ x-1,family=quasibinomial(link=fracreg.links(link)),...)
	p <- results$coefficients
	xbhat <- results$linear.predictors
	yhat <- results$fitted.values
	converged <- results$converged*(1-results$boundary)
	if(method=="ML") LL <- logLik(results)

	ret.list <- list(p=p,yhat=yhat,xbhat=xbhat,converged=converged)
	if(method=="ML") ret.list[["LL"]] <- LL

	if(variance==FALSE) return(ret.list)

	p.var <- fracreg.var(y,x,yhat,xbhat,link,var.type,var.eim,var.cluster,dfc)$var
	ret.list[["p.var"]] <- p.var

	return(ret.list)
}

fracreg.var <- function(y,x,yhat,xbhat,link,var.type,var.eim,var.cluster,dfc)
{
	n <- nrow(x)
	uhat <- y-yhat
	g <- fracreg.links(link)$mu.eta(xbhat)
	gd <- fracreg.links(link)$gd(xbhat)

	A <- 0
	B <- 0

	for(jj in 1:n)
	{
		xx <- x[jj,]%*%t(x[jj,])
		A1 <- xx*(g[jj]^2)/(yhat[jj]*(1-yhat[jj]))
		if(var.eim==TRUE) A <- A+A1
		if(var.eim==FALSE)
		{
			div <- (yhat[jj]*(1-yhat[jj]))
			A2 <- -xx*uhat[jj]*gd[jj]/div
			A3 <- xx*uhat[jj]*(g[jj]^2)/(div*yhat[jj])
			A4 <- -xx*uhat[jj]*(g[jj]^2)/(div*(1-yhat[jj]))

			A <- A+A1+A2+A3+A4
		}
		if(var.type=="robust") B <- B+xx*(uhat[jj]^2)*(g[jj]^2)/((yhat[jj]*(1-yhat[jj]))^2)
	}

	if(var.type=="cluster")
	{
		id <- var.cluster
		id.uni <- unique(id)

		for(j in id.uni)
		{
			Xi <- matrix(x[id==j,],ncol=ncol(x))
			yhati <- yhat[id==j]
			gi <- g[id==j]
			ui <- uhat[id==j]

			ugGi <- ui*gi/(yhati*(1-yhati))

			B <- B+t(Xi)%*%ugGi%*%t(ugGi)%*%Xi
		}
	}

	if(dfc==TRUE)
	{
		if(any(var.type==c("standard","robust"))) df <- n/(n-1)
		if(var.type=="cluster") df <- length(id.uni)/(length(id.uni)-1)
	}
	else df <- 1

	A.inv <- solve(A)

	if(var.type=="standard") var <- df*A.inv
	if(var.type!="standard") var <- df*A.inv%*%B%*%A.inv


	return(list(var=var))
}

fracreg.table <- function(y,yhat,p,p.var,x.names,type,link,converged,var.type)
{
	if(converged==TRUE)
	{
		R2 <- cor(y,yhat)^2

		if(type!="2P" & type!="3P")
		{
			n <- length(y)

			p.sd <- diag(p.var)^0.5
			z.ratio <- p/p.sd
			p.value <- 2*(1-pnorm(abs(z.ratio)))

			stars <- rep("",length(p))
			stars[p.value<=0.01] <- "***"
			stars[p.value>0.01 & p.value<=0.05] <- "**"
			stars[p.value>0.05 & p.value<=0.1] <- "*"

			p <- formatC(p,digits=6,format="f")
			p.sd <- formatC(p.sd,digits=6,format="f")
			z.ratio <- formatC(z.ratio,digits=3,format="f")
			p.value <- formatC(p.value,digits=3,format="f")
			stars <- format(stars,justify="left")

			results <- data.frame(cbind(p,p.sd,z.ratio,p.value,stars),row.names=NULL)

			namcol <- c("Estimate","Std. Error","t value","Pr(>|t|)","")
			dimnames(results) <- list(x.names,namcol)

			cat("\n")
			if(type=="2Pbin") cat("*** Binary component of a two-part model -",link,"specification ***")
			if(type=="1P") cat("*** Fractional",link,"regression model ***")
			if(type=="2Pfrac") cat("*** Fractional component of a two-part model -",link,"specification ***")
			if(type=="3Pbin0") cat("*** First binary component of a three-part model -",link,"specification ***")
			if(type=="3Pbin1") cat("*** Second binary component of a three-part model -",link,"specification ***")
			if(type=="3Pfrac") cat("*** Fractional component of a three-part model -",link,"specification ***")
			cat("\n\n")

			print(results)
			cat("\n")
			if(var.type!="standard")
			{
				cat("Note:",var.type,"standard errors")
				cat("\n\n")
			}
			cat("Number of observations:",n,"\n")
			cat("R-squared:",round(R2,3),"\n")
			cat("\n")
		}
		else
		{
			cat("\n")
		if(type=="2P") cat("*** Two-part model - binary",link[1],"+ fractional",link[2]," ***")
			if(type=="3P") cat("*** Three-part model - binary",link[1],", binary",link[2],"+ fractional",link[3]," ***")
			cat("\n\n")
			cat("R-squared:",round(R2,3),"\n")
			cat("\n")
		}
	}
	else cat("ALGORITHM DID NOT CONVERGE OR STOPPED AT A BOUNDARY VALUE")
	cat("\n")
}

fracreg.lm <- function(y,yhat,gx,gz,type)
{
	u <- y-yhat
	w <- as.vector(sqrt(yhat*(1-yhat)))

	uw <- u/w
	gxw <- gx/w

	gzw <- gz/w

	if(type=="2Pbin")
	{
		gxzw <- cbind(gxw,gzw)

		res <- lm(uw ~ gxzw-1)

		LM <- t(res$fitted.values)%*%uw

	}
	else
	{
		res <- lm(gzw ~ gxw-1)
		ur <- as.matrix(res$residuals)

		n <- length(u)
		ones <- rep(1,n)

		uwr <- as.matrix(uw*ur)
		res <- lm(ones ~ uwr-1)
		uu <- res$residuals
		RSS <- sum(uu^2)

		LM <- n-RSS
	}

	return(list(LM=LM))
}

fracreg.tests.table <- function(test.which,S,Sp,ver,title1,title2=NA,n.ver=NA,test.ggoff=NA)
{
	stars <- rep("",length(S))
	stars[Sp<=0.01] <- "***"
	stars[Sp>0.01 & Sp<=0.05] <- "**"
	stars[Sp>0.05 & Sp<=0.1] <- "*"

	if(all(is.na(S))) stop("NO WALD/LR TEST HAS ACHIEVED CONVERGENCE. USE OTHER TEST VERSIONS INSTEAD")

	S <- formatC(S,digits=3,format="f")
	Sp <- formatC(Sp,digits=3,format="f")
	stars <- format(stars,justify="left")

	if(test.which!="GGOFF")
	{
		results <- data.frame(cbind(ver,S,Sp,stars))
		namcol <- c("Version","Statistic","p-value","")	}
	else
	{
		results <- data.frame(cbind(test.ggoff,ver,S,Sp,stars))
		namcol <- c("Test","Version","Statistic","p-value","")
	}

	results <- results[-1,]
	dimnames(results) <- list(seq_len(nrow(results)),namcol)

	cat("\n")
	cat("***",test.which,"test ***")
	cat("\n\n")
	cat("H0: ",title1)
	cat("\n")
	if(test.which!="P")
	{
		cat("\n")
		print(results,row.names=FALSE)
	}
	else
	{
		cat("H1: ",title2)
		cat("\n\n")
		print(results[1:n.ver,],row.names=FALSE)
		cat("\n")
		cat("H0: ",title2)
		cat("\n")
		cat("H1: ",title1)
		cat("\n\n")
		print(results[(n.ver+1):(2*n.ver),],row.names=FALSE)
	}
	cat("\n")
}

fracreg.pe.var <- function(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var,pb=NA,xbhatb=NA,gb=NA,linkb=NA,pb.var=NA,yhata=NA,yhatb=NA,pc=NA,xbhatc=NA,gc=NA,linkc=NA,pc.var=NA,yhatc=NA)
{
	gda <- fracreg.links(linka)$gd(xbhata)
	PE1.sd <- matrix(NA,nrow=npar,ncol=npar)

	if(type!="2P")
	{
		for(i in 1:npar)
		{
			for(j in 1:npar) PE1.sd[i,j] <- mean(pa[i]*gda*x[,j]+(i==j)*ga)
		}

		PE.sd <- PE1.sd%*%pa.var%*%t(PE1.sd)
	}
	if(type=="3P")
	{
		gdb <- fracreg.links(linkb)$gd(xbhatb)
		gdc <- fracreg.links(linkc)$gd(xbhatc)
		PE2.sd <- matrix(NA,nrow=npar,ncol=npar)
		PE3.sd <- matrix(NA,nrow=npar,ncol=npar)

		for(i in 1:npar)
		{
			for(j in 1:npar)
			{
				H <- yhatb + (1-yhatb)*yhatc
				PE1.sd[i,j] <- mean( (i==j)*ga*H + pa[i]*H*gda*x[,j] + pb[i]*gb*ga*(1-yhatc)*x[,j] + pc[i]*gc*ga*(1-yhatb)*x[,j] )
				PE2.sd[i,j] <- mean( pa[i]*ga*gb*(1-yhatc)*x[,j] + (i==j)*gb*yhata*(1-yhatc) + pb[i]*gdb*yhata*(1-yhatc)*x[,j] - pc[i]*gc*yhata*gb*x[,j] )
				PE3.sd[i,j] <- mean( pa[i]*ga*gc*(1-yhatb)*x[,j] - pb[i]*gb*yhata*gc*x[,j] + (i==j)*gc*yhata*(1-yhatb) + pc[i]*gdc*yhata*(1-yhatb)*x[,j] )
			}
		}

		PE.sd <- PE1.sd%*%pa.var%*%t(PE1.sd) + PE2.sd%*%pb.var%*%t(PE2.sd) + PE3.sd%*%pc.var%*%t(PE3.sd)
	}
	if(type=="2P")
	{
		gdb <- fracreg.links(linkb)$gd(xbhatb)
		PE2.sd <- matrix(NA,nrow=npar,ncol=npar)

		for(i in 1:npar)
		{
			for(j in 1:npar)
			{
				PE1.sd[i,j] <- mean(pb[i]*gb*ga*x[,j]+(i==j)*ga*yhatb+pa[i]*gda*yhatb*x[,j])
				PE2.sd[i,j] <- mean(pa[i]*ga*gb*x[,j]+(i==j)*gb*yhata+pb[i]*gdb*yhata*x[,j])
			}
		}

		PE.sd <- PE1.sd%*%pa.var%*%t(PE1.sd)+PE2.sd%*%pb.var%*%t(PE2.sd)
	}

	PE.sd <- diag(PE.sd)^0.5
	if(any(x.names=="INTERCEPT")) PE.sd <- PE.sd[-1]

	names(PE.sd) <- xvar.names
	PE.sd <- PE.sd[which.x]

	return(PE.sd)
}

fracreg.pe.table <- function(PE.p,PE.sd,PE.type,which.x,xvar.names,title,at)
{
	z.ratio <- PE.p/PE.sd
	p.value <- 2*(1-pnorm(abs(z.ratio)))

	stars <- rep("",length(PE.p))
	stars[p.value<=0.01] <- "***"
	stars[p.value>0.01 & p.value<=0.05] <- "**"
	stars[p.value>0.05 & p.value<=0.1] <- "*"

	PE.p <- formatC(PE.p,digits=4,format="f")
	PE.sd <- formatC(PE.sd,digits=4,format="f")
	z.ratio <- formatC(z.ratio,digits=3,format="f")
	p.value <- formatC(p.value,digits=3,format="f")
	stars <- format(stars,justify="left")

	results <- data.frame(cbind(PE.p,PE.sd,z.ratio,p.value,stars),row.names=NULL)

	namcol <- c("Estimate","Std. Error","t value","Pr(>|t|)","")
	dimnames(results) <- list(which.x,namcol)

	cat("\n\n")
	if(PE.type=="APE") cat("*** Average partial effects ***")
	if(PE.type=="CPE") cat("*** Conditional partial effects ***")
	cat("\n\n")
	cat(title)
	cat("\n\n")
	print(results)
	cat("\n")
	if(PE.type=="CPE")
	{
		cat("------------------")

		if(length(at)==1)
		{
			if(any(at==c("mean","median"))) cat("\nNote: covariates evaluated at",at,"(or mode, for dummies) values\n")
			else
			{
				names(at) <- xvar.names
				cat("\nNote: covariates evaluated at the following values:\n\n")
				print(at)
			}
		}
		else
		{
			names(at) <- xvar.names
			cat("\nNote: covariates evaluated at the following values:\n\n")
			print(at)
		}
	}
}
