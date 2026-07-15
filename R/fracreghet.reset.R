#' @title RESET Test for Fractional Response Regressions under Neglected Heterogeneity
#'
#' @description
#' \code{fracreghet.reset} is used to test the specification of fractional response models estimated by GMMx or LINx.
#' @param object an object containing the results of an \code{fracreghet} command.
#' @param lastpower.vec a numeric vector containing the maximum powers of the linear predictors to be used in RESET tests.
#' @param version a vector containing the test versions to use. Available options: \code{Wald} (the default) and \code{LM} (only  available for \code{GMMx}).
#' @param table a logical value indicating whether a summary table with the test results should be printed.
#' @param \dots Arguments to pass to \link[stats]{nlminb}, which is used to estimate the model under the alternative hypothesis when \code{version} is equal to \code{"Wald"} and the null model was estimated by \code{GMMx}.
#'
#' @details
#' \code{fracreghet.reset} applies the RESET test statistic to fractional response
#' models estimated via \code{fracreghet} using the options \code{GMMx} or \code{LINx}. \code{fracreghet.reset} may be used to test simultaneously the validity of the link specification and the transformation applied to the response variable by each estimator.  
#' 
#' \strong{RESET Test under Unobserved Heterogeneity:}
#' The test is based on augmenting the original model with powers of the linear predictor \eqn{x\hat{\beta}}. For GMMx, it tests \eqn{H_0: \gamma = 0} in the expanded moment conditions:
#' \deqn{E\left[Z_i \left(H(y_i) - \exp\left(x_i\beta + \sum_{k=2}^P \gamma_k (x_i\hat{\beta})^k\right)E(e^{c_i})\right)\right] = 0}
#' This simultaneously evaluates whether the mean function and the specific heterogeneity transformation \eqn{H(\cdot)} are correctly specified.
#' 
#' It is taken into account the option that was chosen for computing standard errors in the model under evaluation. See Ramalho and Ramalho (2017) for details.
#'
#' @return
#' \code{fracreghet.reset} returns a named vector with the test results.
#'
#' @references
#' Ramalho, E. A., & Ramalho, J. J. S. (2017), "Moment-based estimation of nonlinear regression models with boundary outcomes and endogeneity, with applications to nonnegative and fractional responses", \emph{Econometric Reviews}, 36(4), 397-420.
#' 
#' Ramsey, J.B. (1969), "Tests for Specification Errors in Classical Linear Least-Squares Regression Analysis", \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, 31(2), 350-371.
#'
#' @author Sulman Olieko Owili <oliekosulman@gmail.com>
#'
#' @seealso
#' \code{\link{fracreghet}}, for fitting fractional response models under unobserved heterogeneity.\cr
#' \code{\link{fracreghet.pe}}, for computing partial effects.
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
#' res_emp <- fracreghet(y_adj, X_het, type="GMMx", link="logit") 
#' reset_res <- fracreghet.reset(res_emp)
#' summary(reset_res)
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
#' #LM and Wald versions of the RESET test, based on 1 or 2 fitted powers of xb
#' reset_res <- fracreghet.reset(mod,2:3,c("Wald","LM"))
#' summary(reset_res)
#' @export
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
	if(is.null(object$var.type)) stop("fracreghet command was run with variance = FALSE")

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
				results <- fracreghet.est(type,X,X,link,c(object$p,rep(0,df)),Hy,TRUE,var.type,var.cluster,NA,NA,...)
				converged <- results$converged
			}
			if(type=="LINx")
			{
				results <- lm(Hy ~ 0+X)
				converged <- TRUE
			}

			if(converged==TRUE)
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
					p.var1 <- fracreghet.var(type,p1,XB,X,X,link,Hy,var.type,var.cluster,FALSE,NA,NA)$p.var
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
	if(any(!is.na(S)) & table==TRUE) do.call(fracreghet.tests.table, table.info)
	if(all(is.na(S))) warning("RESET test could not be computed; either algorithm did not converge (Wald version) or covariance matrix is singular (Wald/LM versions)")

	### 4. Return results

	statistics <- S[-1]
	names(statistics) <- ver[-1]
	
	class(statistics) <- c("fracreghet.reset", "numeric")
	attr(statistics, "table.info") <- table.info

	return(invisible(statistics))
}
