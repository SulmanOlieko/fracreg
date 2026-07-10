fracreg.ggoff <- function(object,version="LM",table=FALSE,...)
{
	### 1. Error and warning messages

	if(missing(object)) stop("object is missing")
	if(is.null(object$class)) stop("object is not the output of an fracreg command")
	if(object$class!="fracreg") stop("object is not the output of an fracreg command")
	if(object$type=="2P" | object$type=="3P") stop("GGOFF tests are not applicable to two-part or three-part models")
	if(object$converged==0) stop("object is not the output of a successful (converged) fracreg command")
	if(object$method!="ML" & any(version=="LR")) stop("LR tests require ML estimation")
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

	ver <- NA
	S <- NA
	Sp <- NA
	test <- NA
	tests <- c("GOFF1","GOFF2","GGOFF")

	g <- fracreg.links(link)$mu.eta(xbhat)
	gx <- g*x
	if(link!="loglog") z1 <- yhat*log(yhat)/g
	if(link!="cloglog") z2 <- (1-yhat)*log(1-yhat)/g

	for(j in 1:3)
	{
		if((j==1 & link!="loglog") | (j==2 & link!="cloglog") | j==3)
		{
			if(j==1) z <- z1
			if(j==2) z <- z2
			if(j==3)
			{
				if(all(link!=c("loglog","cloglog"))) z <- cbind(z1,z2)
				if(link=="loglog") z <- z2
				if(link=="cloglog") z <- z1
			}
			z <- as.matrix(z)

			df <- ncol(z)

			if(any(version=="LM"))
			{
				test <- c(test,tests[j])
				ver <- c(ver,"LM")

				gz <- g*z
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
					test <- c(test,tests[j])
					ver <- c(ver,"LR")

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
					test <- c(test,tests[j])
					ver <- c(ver,"Wald")

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
	}

	table.info <- list(test.which="GGOFF",S=S,Sp=Sp,ver=ver,title1=title,test.ggoff=test)
	if(table==TRUE) do.call(fracreg.tests.table, table.info)

	### 4. Return results

	statistics <- S[-1]
	names(statistics) <- paste(test[-1],ver[-1],sep="-")

	class(statistics) <- c("fracreg.ggoff", "numeric")
	attr(statistics, "table.info") <- table.info

	return(invisible(statistics))
}
