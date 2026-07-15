#' @title RESET Test for Fractional Response Regressions
#'
#' @description
#' \code{fracreg.reset} is used to perform the Regression Equation Specification Error Test (RESET) to check the functional form and specification of fractional response models.
#' @param object an object containing the results of an \code{fracreg} command.
#' @param lastpower.vec a numeric vector containing the maximum powers of the linear predictors to be used in RESET tests.
#' @param version a vector containing the test versions to use. Available options: \code{Wald}, \code{LM} (the default) and, only for the binary component of two-part models, \code{LR}. More than one option may be chosen.
#' @param table a logical value indicating whether a summary table with the test results should be printed.
#' @param \dots Arguments to pass to \link[stats]{glm}, which is used to estimate the model under the alternative hypothesis when \code{version} is a vector containing \code{"Wald"} or \code{"LR"}.
#'
#' @details
#' \code{fracreg.reset} applies the RESET test statistic to fractional response
#' models estimated via \code{fracreg}. \code{fracreg.reset} may be used to test the link specification of: (i) one-part fractional response models; (ii) the binary
#' components of two-part and three-part fractional response models; and (iii) the fractional components of two-part and three-part fractional response models. 
#' 
#' \strong{RESET Test Framework:}
#' The Regression Equation Specification Error Test (RESET) assesses whether the link function \eqn{G(\cdot)} and the linear index \eqn{x\beta} are correctly specified. It tests the null hypothesis \eqn{H_0: \gamma = 0} in the augmented model:
#' \deqn{E(y|x) = G(x\beta + \sum_{k=2}^P \gamma_k (x\hat{\beta})^k)}
#' where \eqn{P} is the maximum power of the linear predictor (specified by \code{lastpower.vec}) and \eqn{\hat{\beta}} are the estimated parameters from the baseline model.
#' 
#' When the \code{Wald} version is implemented, it is taken into account the option that was chosen for computing standard errors in the model under evaluation. For the \code{LM} version, a robust version is computed in cases (i) and (iii) and a conventional version in case (ii). See Ramalho, Ramalho and Murteira (2011) for details on the application of the RESET test in the fractional response framework.
#'
#' @return
#' \code{fracreg.reset} returns a named vector with the test results.
#'
#' @references
#' Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
#' estimating and testing empirical strategies for fractional response models",
#' \emph{Journal of Economic Surveys}, 25(1), 19-68.
#' 
#' Ramsey, J.B. (1969), "Tests for Specification Errors in Classical Linear Least-Squares Regression Analysis", \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, 31(2), 350-371.
#'
#' @author Sulman Olieko Owili <oliekosulman@gmail.com>
#'
#' @seealso
#' \code{\link{fracreg}}, for fitting fractional response models.\cr
#' \code{\link{fracreg.ggoff}}, for asymptotically equivalent specification tests.\cr
#' \code{\link{fracreg.ptest}}, for non-nested hypothesis tests.\cr
#' \code{\link{fracreg.pe}}, for computing partial effects.
#'
#' @examples
#' ### Empirical 401(k) Examples
#' data("fracreg_k401k")
#' y <- fracreg_k401k$prate
#' X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
#'            totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)
#' 
#' m <- fracreg(y, X, type="1P", linkfrac="logit")
#' reset_res <- fracreg.reset(m)
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
#' ym <- exp(X[,1]+X[,2]+u)/(1+exp(X[,1]+X[,2]+u))
#' y <- rbeta(N,ym*20,20*(1-ym))
#' y[y > 0.9] <- 1
#' 
#' #Testing the logit specification of a standard fractional response model
#' #using LM and Wald versions of the RESET test, based on 1 or 2 fitted powers of
#' #the linear predictor
#' mod <- fracreg(y,X,linkfrac="logit")
#' reset_res <- fracreg.reset(mod,2:3,c("Wald","LM"))
#' summary(reset_res)
#' 
#' #Testing the probit specification of the binary component of a two-part fractional
#' #regression model using LR-based RESET tests with quadratic and cubic fitted 
#' #powers of the linear predictor
#' mod <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1)
#' reset_res <- fracreg.reset(mod,3,"LR")
#' summary(reset_res)
#' @export
fracreg.reset <- function(object,lastpower.vec=3,version="LM",table=FALSE,...)
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
		if(is.null(object$var.type)) stop("Wald test required but fracreg command was run with variance = FALSE")
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

	table.info <- list(test.which="RESET",S=S,Sp=Sp,ver=ver,title1=title)
	if(table==TRUE) do.call(fracreg.tests.table, table.info)

	### 4. Return results

	statistics <- S[-1]
	names(statistics) <- ver[-1]
	
	class(statistics) <- c("fracreg.reset", "numeric")
	attr(statistics, "table.info") <- table.info

	return(invisible(statistics))
}
