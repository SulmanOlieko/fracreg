fracreghet.reset <- function(object,lastpower.vec=3,version="Wald",table=FALSE,...)
{
	### 1. Error and warning messages

	if(missing(object)) stop("object is missing")
	if(is.null(object$class)) stop("object is not the output of an fracreghet command")
	if(object$class!="fracreghet") stop("object is not the output of an fracreghet command")
	if(object$converged==0) stop("object is not the output of a successful (converged) fracreghet command")
	if(any(lastpower.vec<2)) stop(sQuote(lastpower.vec)," - lastpower.vec contains elements lower than 2")
	if(all(object$type!=c("GMMx","LINx"))) stop("fracreghet.reset is only implemented for GMMx and LINx estimators")
	if(all(version!="LM") & all(version!="Wald")) stop("test version not correctly specified")
	if(any(version=="LM") & object$type=="LINx") stop("LM version not implemented for LINx; choose Wald version")
	if(is.null(object$var.type)) stop("fracreghet command was run with variance = F")

	### 2. Recovering definitions and estimates

	mf <- model.frame(object$formula)
	y <- model.response(mf)
	x <- model.matrix(object$formula)

	type <- object$type
	xbhat <- object$xbhat
	Hy <- object$Hy
	link <- object$link
	adjust <- object$adjust

	if(any(version=="Wald") | type=="GMMx")
	{
		var.type <- object$var.type
		if(var.type=="cluster") var.cluster <- object$var.cluster
	}

	title1 <- paste("Fractional",link,"regression")
	title2 <- paste("Estimator:",type)
	if(adjust!=0)
	{
		if(is.numeric(adjust)) title3 <- paste("(adjustment:",adjust,"added to all observations)")
		else title3 <- "(adjustment: all boundary observations dropped)"
		title2 <- paste(title2,title3,sep=" ")
	}

	### 3. Test

	lastpower.vec <- round(lastpower.vec,0)

	N <- length(Hy)
	g <- fracreghet.links(link)$mu.eta(xbhat)
	gx <- g*x

	xx.all <- as.matrix(xbhat^2)
	if(max(lastpower.vec)>2) for(i in 3:max(lastpower.vec)) xx.all <- cbind(xx.all,xbhat^i)

	ver <- NA
	S <- NA
	Sp <- NA

	for(m in lastpower.vec)
	{
		df <- m-1
		xx <- xx.all[,1:(m-1)]

		X <- cbind(x,xx)
		k <- ncol(X)

		if(any(version=="LM"))
		{
			name <- paste("LM(",m,")",sep="")
			ver <- c(ver,name)

			if(type=="GMMx")
			{
				gi <- fracreghet.gi(type,X,Hy,xbhat,link)$gi
				gn <- as.matrix(apply(gi,1,mean))

				if(var.type=="robust") fi <- (1/N)*gi%*%t(gi)
				if(var.type=="cluster")
				{
					fi <- 0
					id <- var.cluster
					u <- Hy-xbhat

					for(j in unique(id))
					{
						Zi <- matrix(X[id==j,],ncol=k)
						ui <- u[id==j]
						zu <- t(Zi)%*%ui

						fi <- fi+zu%*%t(zu)
					}

					fi <- fi/N
				}

				fi.inv <- tryCatch(solve(fi),error=function(e) NaN)
				if(any(is.nan(fi.inv))) fi.inv <- "singular"

				if(!is.character(fi.inv))
				{
					Sj <- N*t(gn)%*%fi.inv%*%gn
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

		if(any(version=="Wald"))
		{
			name <- paste("Wald(",m,")",sep="")
			ver <- c(ver,name)

			if(type=="GMMx")
			{
				results <- fracreghet.est(type,X,X,link,c(object$p,rep(0,df)),Hy,T,var.type,var.cluster,NA,NA,...)
				converged <- results$converged
			}
			if(type=="LINx")
			{
				results <- lm(Hy ~ 0+X)
				converged <- T
			}

			if(converged==T)
			{
				if(type=="GMMx")
				{
					p1 <- results$p
					p.var1 <- results$p.var
				}
				if(type=="LINx")
				{ 
					p1 <- results$coefficients
					XB <- results$fitted.values
					p.var1 <- fracreghet.var(type,p1,XB,X,X,link,Hy,var.type,var.cluster,F,NA,NA)$p.var
				}

				if(!is.character(p.var1))
				{
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
			else
			{
				S <- c(S,NA)
				Sp <- c(Sp,NA)
			}
		}
	}

	table.info <- list(test.which="RESET",S=S,Sp=Sp,ver=ver,title1=title1,title2=title2)
	if(any(!is.na(S)) & table==T) do.call(fracreghet.tests.table, table.info)
	if(all(is.na(S))) warning("RESET test could not be computed; either algorithm did not converge (Wald version) or covariance matrix is singular (Wald/LM versions)")

	### 4. Return results

	statistics <- S[-1]
	names(statistics) <- ver[-1]
	
	class(statistics) <- c("fracreghet.reset", "numeric")
	attr(statistics, "table.info") <- table.info

	return(invisible(statistics))
}
