fracreg.reset <- function(object,lastpower.vec=3,version="LM",table=TRUE,...)
{
	### 1. Error and warning messages

	if(missing(object)) stop("object is missing")
	if(is.null(object$class)) stop("object is not the output of an fracreg command")
	if(object$class!="fracreg") stop("object is not the output of an fracreg command")
	if(object$type=="2P" | object$type=="3P") stop("The RESET test is not applicable to two-part or three-part models")
	if(object$converged==0) stop("object is not the output of a successful (converged) fracreg command")
	if(object$method!="ML" & any(version=="LR")) stop("LR tests require ML estimation")
	if(any(lastpower.vec<2)) stop(sQuote(lastpower.vec)," - lastpower.vec contains elements lower than 2")
	if(all(version!="LM") & all(version!="Wald") & all(version!="LR")) stop("test version not correctly specified")
	if(!is.logical(table)) stop("non-logical value assigned to option table")

	### 2. Recovering definitions and estimates

	mf <- model.frame(object$formula)
	y <- model.response(mf)
	x <- model.matrix(object$formula)

	yhat <- object$yhat
	xbhat <- object$xbhat
	method <- object$method
	link <- object$link
	type <- object$type

	if(any(version=="Wald"))
	{
		if(is.null(object$var.type)) stop("Wald test required but fracreg command was run with variance = F")
		var.type <- object$var.type
		var.eim <- object$var.eim
		dfc <- object$dfc
		if(var.type=="cluster") var.cluster <- object$var.cluster
	}
	if(any(version=="LR")) LL0 <- object$LL
	if(type=="1P") title <- paste("Fractional",link,"regression")
	if(type=="2Pbin") title <- paste("Binary",link,"component of a two-part regression")
	if(type=="2Pfrac") title <- paste("Fractional",link,"component of a two-part regression")

	### 3. Tests

	lastpower.vec <- round(lastpower.vec,0)

	g <- fracreg.links(link)$mu.eta(xbhat)
	gx <- g*x

	z.all <- as.matrix(xbhat^2)
	if(max(lastpower.vec)>2) for(i in 3:max(lastpower.vec)) z.all <- cbind(z.all,xbhat^i)

	ver <- NA
	S <- NA
	Sp <- NA

	for(j in lastpower.vec)
	{
		df <- j-1
		z <- z.all[,1:(j-1)]
		gz <- g*z

		if(any(version=="LM"))
		{
			name <- paste("LM(",j,")",sep="")
			ver <- c(ver,name)

			results <- fracreg.lm(y,yhat,gx,gz,type)

			Sj <- results$LM
			S <- c(S,Sj)
			Sp <- c(Sp,1-pchisq(Sj,df))
		}
		if(any(version=="LR") | any(version=="Wald"))
		{
			if(any(version=="Wald")) results <- fracreg.est(y,cbind(x,z),link,method,TRUE,var.type,var.eim,var.cluster,dfc,...)
			else results <- fracreg.est(y,cbind(x,z),link,method,FALSE,...)

			if(any(version=="LR"))
			{
				name <- paste("LR(",j,")",sep="")
				ver <- c(ver,name)

				if(results$converged==TRUE)
				{
					LL1 <- results$LL

					Sj <- 2*(LL1-LL0)
					S <- c(S,Sj)
					Sp <- c(Sp,1-pchisq(Sj,df))
				}
				else
				{
					S <- c(S,NA)
					Sp <- c(Sp,NA)
				}
			}
			if(any(version=="Wald"))
			{
				name <- paste("Wald(",j,")",sep="")
				ver <- c(ver,name)

				if(results$converged==TRUE)
				{
					p1 <- results$p
					p.var1 <- results$p.var

					p.n <- length(p1)
					p1 <- p1[(p.n-df+1):p.n]
					p.var1 <- p.var1[(p.n-df+1):p.n,(p.n-df+1):p.n]

					Sj <- t(p1)%*%solve(p.var1)%*%p1
					S <- c(S,Sj)
					Sp <- c(Sp,1-pchisq(Sj,df))
				}
				else
				{
					S <- c(S,NA)
					Sp <- c(Sp,NA)
				}
			}
		}
	}

	if(table==TRUE) fracreg.tests.table("RESET",S,Sp,ver,title)

	### 4. Return results

	statistics <- S[-1]
	names(statistics) <- ver[-1]

	return(invisible(statistics))
}
