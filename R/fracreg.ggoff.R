#' @title GGOFF Tests for Fractional Response Regressions
#'
#' @description
#' \code{fracreg.ggoff} is used to perform Generalised Goodness-Of-Functional-Form (GGOFF) tests to check the adequacy of the functional form and link specification of fractional response models.
#' @param object an object containing the results of an \code{fracreg} command.
#' @param version a vector containing the test versions to use. Available options: \code{Wald}, \code{LM} (the default) and, only for the binary component of two-part models, \code{LR}. More than one option may be chosen.
#' @param table a logical value indicating whether a summary table with the test results should be printed.
#' @param \dots Arguments to pass to \link[stats]{glm}, which is used to estimate the model under the alternative hypothesis when \code{version} is a vector containing \code{"Wald"} or \code{"LR"}.
#'
#' @details
#' \code{fracreg.ggoff} applies the GGOFF, GOFF1 and GOOFF2 test statistics to fractional response
#' models estimated via \code{fracreg}. \code{fracreg.ggoff} may be used to test the link
#' specification of: (i) one-part fractional response models; (ii) the binary
#' component of two-part fractional response models; and (iii) the fractional
#' component of two-part fractional response models.
#' 
#' \strong{GGOFF Test Framework:}
#' The Generalised Goodness-of-Functional Form (GGOFF) test evaluates the adequacy of the link function \eqn{G(\cdot)}. It is based on augmenting the baseline model with specific directions of departure. The auxiliary testing equation takes the form:
#' \deqn{E(y|x) = G\left(x\beta + \gamma_1 \frac{g'(x\hat{\beta})}{g(x\hat{\beta})} + \gamma_2 x\hat{\beta} \right)}
#' where \eqn{g(\cdot)} and \eqn{g'(\cdot)} are the first and second derivatives of \eqn{G(\cdot)} evaluated at the linear predictor \eqn{x\hat{\beta}}. The test checks \eqn{H_0: \gamma_1 = 0, \gamma_2 = 0}. GOFF1 and GOFF2 are variants testing individual components.
#' 
#' When the \code{Wald} version is implemented, it is taken into account the option that was chosen for computing standard errors in the model under evaluation. For the \code{LM} version, a robust version is computed in cases (i) and (iii) and a conventional version in case (ii). See Ramalho, Ramalho and Murteira (2014) for details on the application of the GGOFF, GOFF1 and GOOFF2 tests in the fractional response framework.
#'
#' @return
#' \code{fracreg.ggoff} returns a named vector with the test results.
#'
#' @references
#' Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2014), "A generalized goodness-of-functional form test for binary and fractional response models", \emph{Manchester School}, 82(4), 488-507.
#' 
#' Pregibon, D. (1980), "Goodness of Link Tests for Generalized Linear Models", \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, 29(1), 15-24.
#'
#' @author Sulman Olieko Owili <oliekosulman@gmail.com>
#'
#' @seealso
#' \code{\link{fracreg}}, for fitting fractional response models.\cr
#' \code{\link{fracreg.reset}}, for asymptotically equivalent specification tests.\cr
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
#' ggoff_res <- fracreg.ggoff(m)
#' summary(ggoff_res)
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
#' #using LM and Wald versions of the GGOFF test, based on 1 or 2 fitted powers of
#' #the linear predictor
#' mod <- fracreg(y,X,linkfrac="logit")
#' ggoff_res <- fracreg.ggoff(mod,c("Wald","LM"))
#' summary(ggoff_res)
#' 
#' #Testing the probit specification of the binary component of a two-part fractional
#' #regression model using a LR-based GGOFF test
#' mod <- fracreg(y,X,linkbin="probit",type="2Pbin",inf=1)
#' ggoff_res <- fracreg.ggoff(mod,"LR")
#' summary(ggoff_res)
#' @export
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
