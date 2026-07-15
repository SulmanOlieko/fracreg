#' @title Fitting Fractional Response Regressions
#'
#' @description
#' \code{fracreg} is used to fit fractional response models, which are appropriate for responses that are proportions, percentages, or fractions restricted to the [0, 1] interval. It supports standard one-part models, two-part hurdle models for modelling boundary values at 0 or 1, and three-part models for double inflation at both 0 and 1.
#' @param y a numeric vector containing the values of the response variable.
#' @param x a numeric matrix, with column names, containing the values of the covariates.
#' @param x2 a numeric matrix, with column names, containing the values of the covariates in the fractional component of two-part models if option \code{type = "2P"} is defined. Defaults to \code{x}.
#' @param linkbin a description of the link function to use in the binary component of a two-part fractional response model, or a vector of two link functions for the two binary components of a three-part model (e.g. \code{c("logit", "probit")}). Available options: \code{logit}, \code{probit}, \code{cauchit}, \code{loglog}, \code{cloglog}.
#' @param linkfrac a description of the link function to use in standard fractional response models or in the fractional component of a two-part fractional response model. Available options: \code{logit}, \code{probit}, \code{cauchit}, \code{loglog}, \code{cloglog}.
#' @param type a description of the model to estimate: a standard one-part model (\code{1P}, the default), a two-part model (\code{2P}), the binary component of a two-part model (\code{2Pbin}), the fractional component of a two-part model (\code{2Pfrac}), or a three-part model (\code{3P}) for double boundary inflation.
#' @param inflation a numeric value indicating which of the extreme values of \code{0} (the default) or \code{1} is the relevant boundary value for defining two-part fractional response models.
#' @param intercept a logical value indicating whether the model should include a constant term or not.
#' @param table a logical value indicating whether a summary table with the regression results should be printed.
#' @param variance a logical value indicating whether the variance of the estimated parameters should be calculated. Defaults to \code{TRUE} whenever \code{table = TRUE}.
#' @param var.type a description of the type of variance of the estimated parameters to be calculated. Options are \code{standard} (recommended for models estimated by maximum likelihood, such as the binary component of two-part models), \code{robust} (recommended for models estimated by quasi-maximum likelihood, such as standard fractional response models or the fractional component of a two-part fractional response model), \code{cluster} (recommended in the case of panel data) and \code{default} (implements the \code{standard} or \code{robust} versions as appropriate).
#' @param var.eim a logical value indicating whether the expected information matrix should be used in the calculation of the variance. When false, the observation information matrix will be used. Defaults to \code{TRUE}.
#' @param var.cluster a numeric vector containing the values of the variable that specifies to which cluster each observation belongs.
#' @param dfc a logical value indicating whether a degrees of freedom correction should be applied to the covariance matrix. Defaults to \code{FALSE}.
#' @param offset an optional numeric vector containing an offset. It must be of the same dimension as the response variable. It specifies that the variable should be included in the model with its coefficient constrained to 1.
#' @param or a logical value indicating whether to report odds ratios. Only valid when the link function is \code{"logit"}. Defaults to \code{FALSE}.
#' @param level a numeric value between 0 and 1 indicating the confidence level for the confidence intervals. Defaults to \code{0.95}.
#' @param na.action A function specifying how to handle missing values, default is \code{stats::na.omit}. If \code{NULL}, no action is taken.
#' @param \dots Arguments to pass to \link[stats]{glm}.
#'
#' @details
#' \code{fracreg} estimates one-part, two-part hurdle, and three-part double-inflated fractional response models; see Ramalho, Ramalho and Murteira (2011) and Fang and Ma (2013) for details on those models. 
#' 
#' \strong{One-Part Fractional Response Regressions (\code{type = "1P"}):}
#' The standard one-part model assumes that the conditional expectation of the fractional response \eqn{y_i \in [0,1]} is given by:
#' \deqn{E(y_i|x_i) = G(x_i \beta)}
#' where \eqn{G(\cdot)} is a known non-linear link function mapping the linear predictor to the unit interval (e.g., logit, probit). The parameters \eqn{\beta} are estimated by maximising the Bernoulli-based quasi-log-likelihood function:
#' \deqn{\ln L_i(\beta) = y_i \ln[G(x_i \beta)] + (1 - y_i) \ln[1 - G(x_i \beta)]}
#' This estimator requires only the correct specification of the conditional mean to yield consistent parameter estimates (Papke and Wooldridge, 1996).
#' 
#' \strong{Two-Part Hurdle Models (\code{type = "2P"}):}
#' When the data exhibits a boundary mass (e.g., at \eqn{y_i = 0}), the two-part hurdle model handles the boundary values separately from the interior fractional values. Let \eqn{y_i^*} be a binary indicator such that \eqn{y_i^* = 1} if \eqn{y_i > 0} and \eqn{y_i^* = 0} otherwise. The probability of observing a boundary value is modelled as:
#' \deqn{P(y_i = 0 | x_{1i}) = 1 - F(x_{1i} \gamma_1)}
#' \deqn{P(y_i > 0 | x_{1i}) = F(x_{1i} \gamma_1)}
#' where \eqn{F(\cdot)} is a binary link function. Conditional on observing an interior fractional value, the response is modelled as:
#' \deqn{E(y_i | x_{2i}, y_i > 0) = G(x_{2i} \beta_2)}
#' The unconditional mean of the response is therefore:
#' \deqn{E(y_i|x_i) = F(x_{1i} \gamma_1) \times G(x_{2i} \beta_2)}
#' 
#' \strong{Three-Part Double Inflated Models (\code{type = "3P"}):}
#' For data containing boundary mass at both \eqn{0} and \eqn{1}, the three-part model estimates two separate binary mechanisms for each boundary and a fractional component for the interior values \eqn{(0, 1)}, extending the two-part logic to double inflation (Fang and Ma, 2013).
#' 
#' \code{fracreg} uses the standard \link[stats]{glm} command to perform the estimations. Therefore, \code{fracreg} is essentially a convenience command, allowing estimation of several alternative fractional response models using the same command. In addition, \code{fracreg} provides an R-squared measure for all models (calculated as the square of the correlation coefficient between the actual and fitted values of the dependent variable), calculates the fitted values of the dependent variable in two-part models and stores the information needed to implement some very useful commands for fractional response models: \link{fracreg.reset} (RESET test), \link{fracreg.ptest} (P test), \link{fracreg.ggoff} (GGOFF tests) and \link{fracreg.pe} (partial effects).
#'
#' @return
#' When \code{type = "1P" or "2Pfrac"}, \code{fracreg} returns a list with the following elements:
#'   \item{class}{"fracreg".
#' }
#'   \item{formula}{the model formula.
#' }
#'   \item{type}{the name of the estimated model.
#' }
#'   \item{link}{the name of the specified link.
#' }
#'   \item{method}{estimation method. Currently, "QML" (quasi-maximum likelihood) for fractional components or models and"ML" (maximum likelihood) for the binary component of two-part models.
#' }
#'   \item{p}{a named vector of coefficients.
#' }
#'   \item{yhat}{the fitted mean values.
#' }
#'   \item{xbhat}{the fitted mean values of the linear predictor.
#' }
#'   \item{converged}{logical. Was the algorithm judged to have converged?
#' }
#'   \item{x.names}{a vector containing the names of the covariates.
#' }
#' 
#' If \code{variance = TRUE} or \code{table = TRUE}, the previous list also contains the following elements:
#'   \item{p.var}{a named covariance matrix.
#' }
#'   \item{var.type}{covariance matrix type.
#' }
#'   \item{var.eim}{logical. Was the expected information matrix used in the computation of the covariance matrix?
#' }
#'   \item{dfc}{logical. Was a degrees of freedom correction used for the computation of the covariance matrix?
#' }
#' 
#' If \code{var.type = "cluster"}, the list also contains the following element:
#'   \item{var.cluster}{the variable that specifies to which cluster each observation belongs.
#' }
#' 
#' When \code{type = "2Pbin"}, \code{fracreg} returns a similar list with the following additional element:
#'   \item{LL}{the value of the log-likelihood.
#' }
#' 
#' When \code{type = "2P"}, \code{fracreg} returns the previous lists, indexed by the prefixes \code{resBIN} and \code{resFRAC}, and the following additional elements:
#'   \item{class}{"fracreg".
#' }
#'   \item{type}{"2P".
#' }
#'   \item{ybase}{a numeric vector containing the values of the response variable.
#' }
#'   \item{x2base}{a numeric matrix containing the values of the covariates.
#' }
#'   \item{yhat2P}{the overall fitted mean values.
#' }
#'   \item{converged}{logical. Were the algorithms judged to have converged in both parts of the model?
#' }
#' 
#' When \code{type = "3P"}, \code{fracreg} returns the previous lists, indexed by the prefixes \code{resBIN0}, \code{resBIN1}, and \code{resFRAC}, and the following additional elements:
#'   \item{class}{"fracreg".
#' }
#'   \item{type}{"3P".
#' }
#'   \item{ybase}{a numeric vector containing the values of the response variable.
#' }
#'   \item{x2base}{a numeric matrix containing the values of the covariates.
#' }
#'   \item{yhat3P}{the overall fitted mean values.
#' }
#'   \item{converged}{logical. Were the algorithms judged to have converged in all parts of the model?
#' }
#'
#' @section Odds Ratios:
#' When \code{or=TRUE} and the fractional link function (\code{linkfrac} or \code{link}) is \code{"logit"}, the model additionally computes odds ratios for the coefficients. 
#' Odds Ratios are exponentiated coefficients.
#' The corresponding standard errors for the odds ratios are calculated using the Delta method.
#' The confidence intervals for the odds ratios are calculated using the adjusted standard errors and the specified \code{level} (defaulting to 95\%).
#' Odds ratios are particularly useful in fractional logit models as they provide a direct multiplicative interpretation of the independent variable on the odds of the fractional outcome.
#'
#' @references
#' Papke, L. E. and Wooldridge, J. M. (1996), "Econometric methods for fractional response variables with an application to 401(k) plan participation rates", \emph{Journal of Applied Econometrics}, 11(6), 619-632.
#' 
#' Ramalho, E.A., J.J.S. Ramalho and J.M.R. Murteira (2011), "Alternative
#' estimating and testing empirical strategies for fractional response models",
#' \emph{Journal of Economic Surveys}, 25(1), 19-68.
#' 
#' Fang, K., & Ma, S. (2013), "Three-part model for fractional response variables with application to Chinese household health insurance coverage", \emph{Journal of Applied Statistics}, 40(5), 925-940.
#'
#' @author Sulman Olieko Owili <oliekosulman@gmail.com>
#'
#' @seealso
#' \code{\link{fracreg.reset}} and \code{\link{fracreg.ggoff}}, for specification tests.\cr
#' \code{\link{fracreg.ptest}}, for non-nested hypothesis tests.\cr
#' \code{\link{fracreg.pe}}, for computing partial effects.\cr
#' \code{fracreghet}, for fitting cross-sectional fractional response models with unobserved heterogeneity.\cr
#' \code{fracregpd}, for fitting panel data fractional response models.
#'
#' @examples
#' ### Empirical 401(k) Examples
#' data("fracreg_k401k")
#' y <- fracreg_k401k$prate
#' X <- cbind(mrate = fracreg_k401k$mrate, age = fracreg_k401k$age, 
#'            totemp = fracreg_k401k$totemp, sole = fracreg_k401k$sole)
#' 
#' # 1P Model
#' mod <- fracreg(y, X, type="1P", linkfrac="logit")
#' summary(mod)
#' 
#' # 1P Model reporting odds ratios and 99% confidence intervals
#' mod <- fracreg(y, X, type="1P", linkfrac="logit", or=TRUE, level=0.99)
#' summary(mod)
#' 
#' # 2P Model (modelling mass at 1)
#' mod <- fracreg(y, X, type="2P", inflation=1, linkbin="logit", linkfrac="logit")
#' summary(mod)
#' 
#' # 3P Model (inject artificial 0s for demonstration)
#' y_3p <- y; y_3p[1:50] <- 0
#' mod <- fracreg(y_3p, X, type="3P", linkbin=c("logit","logit"), linkfrac="logit")
#' summary(mod)
#' 
#' ### Simulated Examples
#' 
#' set.seed(123)
#' N <- 1000
#' x1 <- rnorm(N)
#' x2 <- runif(N)
#' 
#' # Generating a fractional dependent variable with inflation at 0 and 1
#' XB <- -0.5 + 0.8 * x1 + 1.2 * x2 + rnorm(N)
#' y_latent <- exp(XB) / (1 + exp(XB))
#' 
#' y <- y_latent
#' # Inflate at boundaries
#' y[y_latent < 0.2] <- 0
#' y[y_latent > 0.8] <- 1
#' 
#' X <- cbind(x1 = x1, x2 = x2)
#' 
#' # fracreg estimation of a logit fractional response model
#' mod <- fracreg(y, X, type="1P", linkfrac="logit")
#' summary(mod)
#' 
#' # fracreg estimation of the binary logit component of the two-part fractional
#' # regression model with y=0 as the relevant boundary value
#' mod <- fracreg(y, X, type="2Pbin", inflation=0, linkbin="logit")
#' summary(mod)
#' 
#' # fracreg estimation of the fractional component of the two-part fractional
#' # regression model with y=0 as the relevant boundary value and using a
#' # probit link function
#' mod <- fracreg(y, X, type="2Pfrac", inflation=0, linkfrac="probit")
#' summary(mod)
#' 
#' # fracreg estimation of both components of a two-part fractional response model
#' # with y=0 as the relevant boundary value and using a cloglog binary link
#' # function and a logit fractional link function
#' mod <- fracreg(y, X, type="2P", inflation=0, linkbin="cloglog", linkfrac="logit")
#' summary(mod)
#' 
#' # Three-part double-inflated model (y has both 0s and 1s)
#' mod <- fracreg(y, X, type="3P", linkbin=c("logit","probit"), linkfrac="logit")
#' summary(mod)
#' @export
fracreg <- function(y,x,x2=x,linkbin,linkfrac,type="1P",inflation=0,intercept=TRUE,table=FALSE,variance=TRUE,var.type="default",var.eim=TRUE,var.cluster,dfc=FALSE,offset=NULL,or=FALSE,level=0.95,na.action=stats::na.omit,...)
{
	cl <- match.call()
	
	if (!missing(y) && !missing(x)) {
	    if (!missing(var.cluster)) {
	        cleaned <- fracreg_clean_data(y=y, x=x, x2=x2, var.cluster=var.cluster, offset=offset, na.action=na.action)
	        var.cluster <- cleaned$var.cluster
	    } else {
	        cleaned <- fracreg_clean_data(y=y, x=x, x2=x2, offset=offset, na.action=na.action)
	    }
	    y <- cleaned$y
	    x <- cleaned$x
	    x2 <- cleaned$x2
	    offset <- cleaned$offset
	}
	
	### 1. Error and warning messages

	if(missing(y)) stop("dependent variable is missing")
	if(missing(x)) stop("explanatory variables are missing")

	if(all(type!=c("1P","2Pbin","2Pfrac","2P","3P"))) stop(sQuote(type)," - type not recognised")
	if(any(y>1) | any(y<0)) stop("The dependent variable has values outside the unit interval")
	if(all(inflation!=c(0,1))) stop(inflation," - value not recognised for inflation")
	if(length(inflation)>1) stop(inflation," - only a single value allowed for inflation")
	if(type!="1P" & !any(y==inflation)) stop("The dependent variable has no ",sQuote(inflation)," values")

	if(type=="2Pbin")
	{
		if(missing(linkbin)) stop("linkbin is missing")
		if(!missing(linkfrac)) warning(sQuote(type)," and",sQuote(linkfrac)," - type does not use linkfrac")

		if(all(linkbin[1]!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkbin[1])," - linkbin not recognised")
	}
	if(type=="3P")
	{
		if(missing(linkbin)) stop("linkbin is missing")
		if(missing(linkfrac)) stop("linkfrac is missing")
		if(length(linkbin) == 1) linkbin <- c(linkbin, linkbin)
		if(all(linkbin[1]!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkbin[1])," - linkbin not recognised")
		if(all(linkbin[2]!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkbin[2])," - linkbin not recognised")
		if(all(linkfrac!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkfrac)," - linkfrac not recognised")
	}
	if(type=="2P")
	{
		if(missing(linkbin)) stop("linkbin is missing")
		if(missing(linkfrac)) stop("linkfrac is missing")
		if(all(linkbin[1]!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkbin[1])," - linkbin not recognised")
		if(all(linkfrac!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkfrac)," - linkfrac not recognised")
	}
	if(any(type==c("1P","2Pfrac")))
	{
		if(missing(linkfrac)) stop("linkfrac is missing")
		if(!missing(linkbin)) warning(sQuote(type)," and",sQuote(linkbin)," - type does not use linkbin")

		if(all(linkfrac!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkfrac)," - linkfrac not recognised")
	}

	if(table==TRUE & variance==FALSE)
	{
		variance <- TRUE
		warning("option variance changed from FALSE to TRUE, as required by table=TRUE")
	}

	if(all(var.type!=c("standard","robust","cluster","default"))) stop(sQuote(var.type)," - var.type not recognised")
	if(var.type=="cluster" & missing(var.cluster)) stop("option cluster for covariance matrix but no var.cluster supplied")
	if(missing(var.cluster)) var.cluster <- NULL

	if(!is.logical(intercept)) stop("non-logical value assigned to option intercept")
	if(!is.logical(table)) stop("non-logical value assigned to option table")
	if(!is.logical(variance)) stop("non-logical value assigned to option variance")
	if(!is.logical(var.eim)) stop("non-logical value assigned to option var.eim")
	if(!is.logical(dfc)) stop("non-logical value assigned to option dfc")

	### 2. Data and variables preparation

	if(any(type==c("1P","2Pfrac"))) x2 <- x

	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.data.frame(x2)) x2 <- as.matrix(x2)

	if(!is.matrix(x)) stop("x is not a matrix")
	if(!is.matrix(x2)) stop("x2 is not a matrix")

	x.names <- dimnames(x)[[2]]
	x2.names <- dimnames(x2)[[2]]

	if(is.null(x.names)) stop("x has no column names")
	if(is.null(x2.names)) stop("x2 has no column names")

	if(intercept==TRUE)
	{
		x <- cbind(1,x)
		x2 <- cbind(1,x2)

		x.names <- c("(Intercept)",x.names)
		x2.names <- c("(Intercept)",x2.names)
	}

	if(length(x.names)!=length(unique(x.names))) stop("some covariate names in x are identical")
	if(length(x2.names)!=length(unique(x2.names))) stop("some covariate names in x2 are identical")

	if(length(y)!=nrow(x)) stop("the number of observations for y and x are different")
	if(!is.null(offset) && length(y)!=length(offset)) stop("offset does not have the appropriate dimension")
	if(var.type=="cluster")
	{
		if(length(y)!=length(var.cluster)) stop("var.cluster does not have the appropriate dimension")
	}

	if(any(type==c("2Pbin","2P")))
	{
		if(inflation==0) yb <- y>0
		if(inflation==1) yb <- y==1
	}

	if(any(type==c("2Pfrac","2P")))
	{
		var.cluf <- NULL
		offsetf <- NULL
		if(inflation==0)
		{
			yf <- y[y>0]
			x2f <- x2[y>0,]
			if(var.type=="cluster") var.cluf <- var.cluster[y>0]
			if(!is.null(offset)) offsetf <- offset[y>0]
		}
		if(inflation==1)
		{
			yf <- y[y<1]
			x2f <- x2[y<1,]
			if(var.type=="cluster") var.cluf <- var.cluster[y<1]
			if(!is.null(offset)) offsetf <- offset[y<1]
		}

		if(length(yf)!=nrow(x2f)) stop("the number of observations for y and x2 are different")
	}

	if(type=="3P")
	{
		var.clu1 <- NULL
		var.cluf <- NULL
		offset1 <- NULL
		offsetf <- NULL
		yb0 <- y>0
		y1_subset <- y[y>0]
		x_subset <- x[y>0,]
		yb1 <- y1_subset==1
		
		yf <- y[y>0 & y<1]
		x2f <- x2[y>0 & y<1,]
		
		if(var.type=="cluster") {
			var.clu1 <- var.cluster[y>0]
			var.cluf <- var.cluster[y>0 & y<1]
		}
		if(!is.null(offset)) {
			offset1 <- offset[y>0]
			offsetf <- offset[y>0 & y<1]
		}
	}

	### 3. Estimation

	class <- "fracreg"

	if(any(type==c("2Pbin","2P")))
	{
		if(var.type=="default") var.ty <- "standard"
		else var.ty <- var.type

		method <- "ML"
		results <- fracreg.est(yb,x,linkbin,method,variance,var.ty,var.eim,var.cluster,dfc,offset=offset,...)
		p <- results$p
		if(variance==TRUE) p.var <- results$p.var
		yhat1 <- results$yhat
		xbhat <- results$xbhat
		converged1 <- results$converged
		LL <- results$LL

		table.info <- list(y=yb,yhat=yhat1,p=p,p.var=if(variance) p.var else NA,x.names=x.names,type="2Pbin",link=linkbin,converged=converged1,var.type=var.ty,LL=LL,method=method,var.cluster=if(var.type=="cluster") var.cluster else NULL,dfc=dfc,or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info)

		formula <- yb ~ x - 1
		names(p) <- x.names

		resBIN <- list(call=cl, class=class,formula=formula,type=type,link=linkbin,method=method,p=p,yhat=yhat1,xbhat=xbhat,converged=converged1,LL=LL,x.names=x.names,table.info=table.info)
		if(variance==TRUE)
		{ 
			dimnames(p.var) <- list(x.names,x.names)
			resBIN[["p.var"]] <- p.var
			resBIN[["var.type"]] <- var.ty
			resBIN[["var.eim"]] <- var.eim
			resBIN[["dfc"]] <- dfc
			if(var.type=="cluster") resBIN[["var.cluster"]] <- var.cluster
		}
		if(!is.null(offset)) resBIN[["offset"]] <- offset
	}

	if(any(type==c("1P","2Pfrac","2P")))
	{
		var.clu <- NULL
		off.curr <- NULL
		if(type=="1P")
		{
			yy <- y
			xx2 <- x2
			ty <- "1P"
			if(var.type=="cluster") var.clu <- var.cluster
			if(!is.null(offset)) off.curr <- offset
		}
		else
		{
			yy <- yf
			xx2 <- x2f
			ty <- "2Pfrac"
			if(var.type=="cluster") var.clu <- var.cluf
			if(!is.null(offset)) off.curr <- offsetf
		}

		if(var.type=="default") var.ty <- "robust"
		else var.ty <- var.type

		method <- "QML"
		results <- fracreg.est(yy,xx2,linkfrac,method,variance,var.ty,var.eim,var.clu,dfc,offset=off.curr,...)
		p <- results$p
		if(variance==TRUE) p.var <- results$p.var
		yhat <- results$yhat
		xbhat <- results$xbhat
		converged2 <- results$converged

		table.info <- list(y=yy,yhat=yhat,p=p,p.var=if(variance) p.var else NA,x.names=x2.names,type=ty,link=linkfrac,converged=converged2,var.type=var.ty,LL=results$LL,method=method,var.cluster=var.clu,dfc=dfc,or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info)

		formula <- yy ~ xx2 - 1
		names(p) <- x2.names

		resFRAC <- list(call=cl, class=class,formula=formula,type=type,link=linkfrac,method=method,p=p,yhat=yhat,xbhat=xbhat,converged=converged2,x.names=x2.names,table.info=table.info)
		if(variance==TRUE)
		{ 
			dimnames(p.var) <- list(x2.names,x2.names)
			resFRAC[["p.var"]] <- p.var
			resFRAC[["var.type"]] <- var.ty
			resFRAC[["var.eim"]] <- var.eim
			resFRAC[["dfc"]] <- dfc
			if(var.type=="cluster") resFRAC[["var.cluster"]] <- var.clu
		}
		if(!is.null(off.curr)) resFRAC[["offset"]] <- off.curr
	}

	if(type=="3P")
	{
		if(var.type=="default") var.ty <- "standard"
		else var.ty <- var.type

		# Component 1: y > 0
		method <- "ML"
		results0 <- fracreg.est(yb0,x,linkbin[1],method,variance,var.ty,var.eim,var.cluster,dfc,offset=offset,...)
		p0 <- results0$p
		if(variance==TRUE) p.var0 <- results0$p.var
		yhat1 <- results0$yhat
		converged1 <- results0$converged
		table.info0 <- list(y=yb0,yhat=yhat1,p=p0,p.var=if(variance) p.var0 else NA,x.names=x.names,type="3Pbin0",link=linkbin[1],converged=converged1,var.type=var.ty,LL=results0$LL,method=method,var.cluster=if(var.type=="cluster") var.cluster else NULL,dfc=dfc,or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info0)
		names(p0) <- x.names
		resBIN0 <- list(call=cl, class=class,formula=yb0 ~ x - 1,type="3Pbin0",link=linkbin[1],method=method,p=p0,yhat=yhat1,xbhat=results0$xbhat,converged=converged1,LL=results0$LL,x.names=x.names,table.info=table.info0)
		if(variance==TRUE) { dimnames(p.var0) <- list(x.names,x.names); resBIN0[["p.var"]] <- p.var0; resBIN0[["var.type"]] <- var.ty }
		if(!is.null(offset)) resBIN0[["offset"]] <- offset

		# Component 2: y = 1 | y > 0
		if(var.type=="cluster") var.clu1 <- var.cluster[y>0] else var.clu1 <- NULL
		results1 <- fracreg.est(yb1,x_subset,linkbin[2],method,variance,var.ty,var.eim,var.clu1,dfc,offset=offset1,...)
		p1 <- results1$p
		if(variance==TRUE) p.var1 <- results1$p.var
		converged2 <- results1$converged
		table.info1 <- list(y=yb1,yhat=results1$yhat,p=p1,p.var=if(variance) p.var1 else NA,x.names=x.names,type="3Pbin1",link=linkbin[2],converged=converged2,var.type=var.ty,LL=results1$LL,method=method,var.cluster=var.clu1,dfc=dfc,or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info1)
		names(p1) <- x.names
		resBIN1 <- list(call=cl, class=class,formula=yb1 ~ x_subset - 1,type="3Pbin1",link=linkbin[2],method=method,p=p1,yhat=results1$yhat,xbhat=results1$xbhat,converged=converged2,LL=results1$LL,x.names=x.names,table.info=table.info1)
		if(variance==TRUE) { dimnames(p.var1) <- list(x.names,x.names); resBIN1[["p.var"]] <- p.var1; resBIN1[["var.type"]] <- var.ty }
		if(!is.null(offset1)) resBIN1[["offset"]] <- offset1

		# Component 3: 0 < y < 1
		if(var.type=="default") var.tyF <- "robust" else var.tyF <- var.type
		resultsF <- fracreg.est(yf,x2f,linkfrac,"QML",variance,var.tyF,var.eim,var.cluf,dfc,offset=offsetf,...)
		pF <- resultsF$p
		if(variance==TRUE) p.varF <- resultsF$p.var
		converged3 <- resultsF$converged
		table.infoF <- list(y=yf,yhat=resultsF$yhat,p=pF,p.var=if(variance) p.varF else NA,x.names=x2.names,type="3Pfrac",link=linkfrac,converged=converged3,var.type=var.tyF,LL=resultsF$LL,method="QML",var.cluster=var.cluf,dfc=dfc,or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.infoF)
		names(pF) <- x2.names
		resFRAC <- list(call=cl, class=class,formula=yf ~ x2f - 1,type="3Pfrac",link=linkfrac,method="QML",p=pF,yhat=resultsF$yhat,xbhat=resultsF$xbhat,converged=converged3,x.names=x2.names,table.info=table.infoF)
		if(variance==TRUE) { dimnames(p.varF) <- list(x2.names,x2.names); resFRAC[["p.var"]] <- p.varF; resFRAC[["var.type"]] <- var.tyF }
		if(!is.null(offsetf)) resFRAC[["offset"]] <- offsetf
	}

	if(type=="2P")
	{
		off2Pfrac <- if(!is.null(offset)) offset else 0
		yhat2 <- fracreg.links(linkfrac)$linkinv(x2%*%p + off2Pfrac)
		yhat <- yhat1*yhat2

		converged <- converged1*converged2

		table.info <- list(y=y,yhat=yhat,p=NA,p.var=NA,x.names=NA,type=type,link=c(linkbin,linkfrac),converged=converged,var.type="standard",or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info)

		ybase <- y
		x2base <- x2
	}

	if(type=="3P")
	{
		# Overall Expected Value
		off3Pbin2 <- if(!is.null(offset)) offset else 0
		off3Pfrac <- if(!is.null(offset)) offset else 0

		F1 <- yhat1
		F2 <- fracreg.links(linkbin[2])$linkinv(x%*%p1 + off3Pbin2)
		G <- fracreg.links(linkfrac)$linkinv(x2%*%pF + off3Pfrac)
		yhat <- F1 * (F2 + (1 - F2) * G)

		converged <- converged1*converged2*converged3

		table.info <- list(y=y,yhat=yhat,p=NA,p.var=NA,x.names=NA,type=type,link=c(linkbin,linkfrac),converged=converged,var.type="standard",or=or,level=level)
		if(table==TRUE) do.call(fracreg.table, table.info)

		ybase <- y
		x2base <- x2
	}

	### 4. Return results

	if(type=="2Pbin") {
		class(resBIN) <- "fracreg"
		return(invisible(resBIN))
	}
	if(any(type==c("1P","2Pfrac"))) {
		class(resFRAC) <- "fracreg"
		return(invisible(resFRAC))
	}
	if(type=="2P") {
		res <- list(call=cl, resBIN=resBIN,resFRAC=resFRAC,class=class,type=type,ybase=ybase,x2base=x2base,yhat2P=yhat,converged=converged,table.info=table.info)
		class(res) <- "fracreg"
		return(invisible(res))
	}
	if(type=="3P") {
		res <- list(call=cl, resBIN0=resBIN0,resBIN1=resBIN1,resFRAC=resFRAC,class=class,type=type,ybase=ybase,x2base=x2base,yhat3P=yhat,converged=converged,table.info=table.info)
		class(res) <- "fracreg"
		return(invisible(res))
	}
}
