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

.fracreg.sep <- function(width = 80) paste0(rep("-", width), collapse = "")

.fracreg.center <- function(text, width = 80) {
    spaces <- max(0, floor((width - nchar(text)) / 2))
    paste0(paste0(rep(" ", spaces), collapse = ""), text)
}

.fracreg.cat.right <- function(label, value, width = 80) {
    val_str <- as.character(value)
    spaces <- max(1, width - nchar(label) - nchar(val_str) - 2)
    cat(label, paste0(rep(" ", spaces), collapse=""), val_str, "\n")
}

fracreg.est <- function(y,x,link,method,variance,var.type,var.eim,var.cluster,dfc,...)
{
	if(method=="ML") results <- glm(y ~ x-1,family=binomial(link=fracreg.links(link)),...)
	if(method=="QML") results <- glm(y ~ x-1,family=quasibinomial(link=fracreg.links(link)),...)
	p <- results$coefficients
	xbhat <- results$linear.predictors
	yhat <- results$fitted.values
	converged <- results$converged*(1-results$boundary)
	if(method=="ML") LL <- as.numeric(logLik(results))
	if(method=="QML") {
		eps <- 1e-16
		LL <- sum(ifelse(y > 0, y * log(pmax(yhat, eps)), 0) + ifelse(y < 1, (1-y) * log(pmax(1-yhat, eps)), 0))
	}

	ret.list <- list(p=p,yhat=yhat,xbhat=xbhat,converged=converged)
	ret.list[["LL"]] <- LL

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

fracreg.pe.var <- function(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var,pb=NA,xbhatb=NA,gb=NA,linkb=NA,pb.var=NA,yhata=NA,yhatb=NA,pc=NA,xbhatc=NA,gc=NA,linkc=NA,pc.var=NA,yhatc=NA)
{
	gda <- fracreg.links(linka)$gd(xbhata)
	PE1.sd <- matrix(NA,nrow=npar,ncol=npar)

	if(type!="2P" & type!="3P")
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
	if(any(x.names=="(Intercept)")) PE.sd <- PE.sd[-1]

	names(PE.sd) <- xvar.names
	PE.sd <- PE.sd[which.x]

	return(PE.sd)
}

#' @export
summary.fracreg <- function(object, ...) {
    if(!is.null(object$resBIN0) && !is.null(object$resBIN0$table.info)) do.call(fracreg.table, object$resBIN0$table.info)
    if(!is.null(object$resBIN1) && !is.null(object$resBIN1$table.info)) do.call(fracreg.table, object$resBIN1$table.info)
    if(!is.null(object$resBIN) && !is.null(object$resBIN$table.info)) do.call(fracreg.table, object$resBIN$table.info)
    if(!is.null(object$resFRAC) && !is.null(object$resFRAC$table.info)) do.call(fracreg.table, object$resFRAC$table.info)
    if(!is.null(object$table.info)) do.call(fracreg.table, object$table.info)
    invisible(object)
}

#' @export
print.fracreg <- function(x, ...) {
        title <- "Fractional response regression"
    if(!is.null(x$type)) {
        if(x$type == "1P") title <- paste("Fractional", x$link, "regression")
        if(x$type == "2Pbin") title <- paste("Part 1: Binary", x$link, "regression")
        if(x$type == "2Pfrac") title <- paste("Part 2: Fractional", x$link, "regression")
        if(x$type == "3Pbin0") title <- paste("Part 1: Binary", x$link, "regression")
        if(x$type == "3Pbin1") title <- paste("Part 2: Binary", x$link, "regression")
        if(x$type == "3Pfrac") title <- paste("Part 3: Fractional", x$link, "regression")
        if(x$type == "2P") title <- paste("Two-part fractional regression: binary", x$link[1], "+ fractional", x$link[2])
        if(x$type == "3P") title <- paste("Three-part fractional regression: binary", x$link[1], ", binary", x$link[2], "+ fractional", x$link[3])
    }
    cat(paste0("\n", title, "\n"))
    if(!is.null(x$call)) {
        cat("\nCall:\n")
        print(x$call)
    }
    
    if(!is.null(x$resBIN0)) {
        cat("\nComponent 1 (Binary 0):\nCoefficients:\n")
        print.default(x$resBIN0$p)
    }
    if(!is.null(x$resBIN1)) {
        cat("\nComponent 2 (Binary 1):\nCoefficients:\n")
        print.default(x$resBIN1$p)
    }
    if(!is.null(x$resBIN)) {
        cat("\nBinary Component:\nCoefficients:\n")
        print.default(x$resBIN$p)
    }
    if(!is.null(x$resFRAC)) {
        cat("\nFractional Component:\nCoefficients:\n")
        print.default(x$resFRAC$p)
    }
    if(is.null(x$resBIN0) && is.null(x$resBIN1) && is.null(x$resBIN) && is.null(x$resFRAC) && !is.null(x$p)) {
        cat("\nCoefficients:\n")
        print.default(x$p)
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.fracregpd <- function(object, ...) {
    if(!is.null(object$table.info)) do.call(fracregpd.table, object$table.info)
    invisible(object)
}

#' @export
print.fracregpd <- function(x, ...) {
        model_desc <- "panel data"
    if (!is.null(x$type) && x$type %in% c("QMLcre", "GMMcre")) model_desc <- "correlated random effects"
    link_str <- if(!is.null(x$link)) x$link else ""
    title <- paste("Fractional", link_str, model_desc, "regression")
    cat(paste0("\n", title, "\n"))
    if(!is.null(x$call)) {
        cat("\nCall:\n")
        print(x$call)
    }
    if(!is.null(x$p)) {
        cat("\nCoefficients:\n")
        print.default(x$p)
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreghet <- function(object, ...) {
    if(!is.null(object$table.info)) do.call(fracreghet.table, object$table.info)
    invisible(object)
}

#' @export
print.fracreghet <- function(x, ...) {
        link_str <- if(!is.null(x$link)) x$link else ""
    title <- paste("Fractional", link_str, "regression with heteroscedasticity/endogeneity")
    cat(paste0("\n", title, "\n"))
    if(!is.null(x$call)) {
        cat("\nCall:\n")
        print(x$call)
    }
    if(!is.null(x$p)) {
        cat("\nCoefficients:\n")
        print.default(x$p)
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreg.pe <- function(object, ...) {
    if(!is.null(object$table.info)) {
        do.call(fracreg.pe.table, object$table.info)
    } else if(!is.null(object$ape) || !is.null(object$cpe)) {
        if(!is.null(object$ape$table.info)) do.call(fracreg.pe.table, object$ape$table.info)
        if(!is.null(object$cpe$table.info)) do.call(fracreg.pe.table, object$cpe$table.info)
    }
    invisible(object)
}

#' @export
print.fracreg.pe <- function(x, ...) {
    cat("\nPartial Effects for Fractional Regression\n")
    if(!is.null(x$ape)) {
        cat("\nAverage Partial Effects (APE):\n")
        print.default(x$ape$table.info$PE.p)
    }
    if(!is.null(x$cpe)) {
        cat("\nConditional Partial Effects (CPE):\n")
        print.default(x$cpe$table.info$PE.p)
    }
    if(!is.null(x$table.info)) {
        if(x$table.info$PE.type == "APE") cat("\nAverage Partial Effects (APE):\n")
        else cat("\nConditional Partial Effects (CPE):\n")
        print.default(x$table.info$PE.p)
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.fracregpd.pe <- function(object, ...) {
    if(!is.null(object$table.info)) {
        do.call(fracreg.pe.table, object$table.info)
    } else if(!is.null(object$ape) || !is.null(object$cpe)) {
        if(!is.null(object$ape$table.info)) do.call(fracreg.pe.table, object$ape$table.info)
        if(!is.null(object$cpe$table.info)) do.call(fracreg.pe.table, object$cpe$table.info)
    }
    invisible(object)
}

#' @export
print.fracregpd.pe <- function(x, ...) {
    cat("\nPartial Effects for Fractional Panel Data Regression\n")
    if(!is.null(x$ape)) {
        cat("\nAverage Partial Effects (APE):\n")
        print.default(x$ape$table.info$PE.p)
    }
    if(!is.null(x$cpe)) {
        cat("\nConditional Partial Effects (CPE):\n")
        print.default(x$cpe$table.info$PE.p)
    }
    if(!is.null(x$table.info)) {
        if(x$table.info$PE.type == "APE") cat("\nAverage Partial Effects (APE):\n")
        else cat("\nConditional Partial Effects (CPE):\n")
        print.default(x$table.info$PE.p)
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreghet.pe <- function(object, ...) {
    if(!is.null(object$table.info)) {
        do.call(fracreghet.pe.table, object$table.info)
    } else if(!is.null(object$ape) || !is.null(object$cpe)) {
        if(!is.null(object$ape$table.info)) do.call(fracreghet.pe.table, object$ape$table.info)
        if(!is.null(object$cpe$table.info)) do.call(fracreghet.pe.table, object$cpe$table.info)
    }
    invisible(object)
}

#' @export
print.fracreghet.pe <- function(x, ...) {
    cat("\nPartial Effects for Fractional Regression with Heteroskedasticity\n")
    if(!is.null(x$ape)) {
        cat("\nAverage Partial Effects (APE):\n")
        print.default(x$ape$table.info$PE.p)
    }
    if(!is.null(x$cpe)) {
        cat("\nConditional Partial Effects (CPE):\n")
        print.default(x$cpe$table.info$PE.p)
    }
    if(!is.null(x$table.info)) {
        if(x$table.info$PE.type == "APE") cat("\nAverage Partial Effects (APE):\n")
        else cat("\nConditional Partial Effects (CPE):\n")
        print.default(x$table.info$PE.p)
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreg.reset <- function(object, ...) {
    ti <- attr(object, "table.info")
    if(!is.null(ti)) do.call(fracreg.tests.table, ti)
    invisible(object)
}

#' @export
print.fracreg.reset <- function(x, ...) {
    cat("\nRESET Test for Fractional Regression\n\n")
    print.default(c(x))
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreghet.reset <- function(object, ...) {
    ti <- attr(object, "table.info")
    if(!is.null(ti)) do.call(fracreghet.tests.table, ti)
    invisible(object)
}

#' @export
print.fracreghet.reset <- function(x, ...) {
    cat("\nRESET Test for Fractional Regression with Heteroskedasticity\n\n")
    print.default(c(x))
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreg.ggoff <- function(object, ...) {
    ti <- attr(object, "table.info")
    if(!is.null(ti)) do.call(fracreg.tests.table, ti)
    invisible(object)
}

#' @export
print.fracreg.ggoff <- function(x, ...) {
    cat("\nGOFF Test for Fractional Regression\n\n")
    print.default(c(x))
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreg.ptest <- function(object, ...) {
    ti <- attr(object, "table.info")
    if(!is.null(ti)) do.call(fracreg.tests.table, ti)
    invisible(object)
}

#' @export
print.fracreg.ptest <- function(x, ...) {
    cat("\nP-Test for Fractional Regression\n\n")
    print.default(c(x))
    cat("\n")
    invisible(x)
}

#' @export
print.fracregridge <- function(x, ...) {
    cat("\nFractional Ridge Regression\n")
    if(!is.null(x$call)) {
        cat("\nCall:\n")
        print(x$call)
    }
    cat("\nRidge Coefficients at Target Fractions:\n")
    print(x$coef)
    cat("\n")
    invisible(x)
}

#' @export
summary.fracregridge <- function(object, ...) {
    cat("\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center("Fractional Ridge Regression"), "\n")
    cat(.fracreg.sep(), "\n")
    
    .fracreg.cat.right("Data type:", "Cross-sectional")
    .fracreg.cat.right("Convergence:", "Successful")
    cat(.fracreg.sep(), "\n")
    
    for (name in names(object$table.info)) {
        cat(.fracreg.center(paste("Target Fraction:", name)), "\n")
        cat(.fracreg.sep(), "\n")
        if(!is.null(object$stats.info)) {
            stats <- object$stats.info[[name]]
            .fracreg.cat.right("Number of observations:", stats$n_obs)
            .fracreg.cat.right("Pseudo R-squared:", round(stats$R2, 5))
            .fracreg.cat.right("Degrees of freedom:", round(stats$n_obs - stats$df_alpha, 2))
            cat(.fracreg.sep(), "\n")
        }
        suppressWarnings(stats::printCoefmat(object$table.info[[name]], P.values=TRUE, has.Pvalue=TRUE, digits=4, signif.legend=TRUE))
        cat(.fracreg.sep(), "\n")
    }
    
    cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
    cat(.fracreg.sep(), "\n\n")
    invisible(object)
}

#' @export
summary.fracregridge.pe <- function(object, ...) {
    cat("\n\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center("Average partial effects"), "\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center("Fractional Ridge Regression"), "\n")
    cat(.fracreg.sep(), "\n")
    
    for (name in names(object$table.info)) {
        cat(.fracreg.center(paste("Target Fraction:", name)), "\n")
        cat(.fracreg.sep(), "\n")
        suppressWarnings(stats::printCoefmat(object$table.info[[name]], P.values=TRUE, has.Pvalue=TRUE, digits=4, signif.legend=TRUE))
        cat(.fracreg.sep(), "\n")
    }
    cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
    cat(.fracreg.sep(), "\n\n")
    invisible(object)
}

#' @export
print.fracregridge.pe <- function(x, ...) {
    cat("\nPartial Effects for Fractional Ridge Regression\n")
    if(!is.null(x$call)) {
        cat("\nCall:\n")
        print(x$call)
    }
    cat("\n")
    invisible(x)
}

#' @export
print.fracregmlogit <- function(x, ...) {
    cat("\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center("Fractional multinomial logit model"), "\n")
    cat(.fracreg.sep(), "\n")
    if(!is.null(x$call)) {
        cat("\nCall:\n")
        print(x$call)
    }
    cat("\nCoefficients:\n")
    print.default(x$coefficient)
    cat("\n")
    invisible(x)
}

#' Generate Summary Tables for fracregmlogit Objects
#' 
#' Generate tables of coefficient estimates, partial effects, and willingness to pay from
#' fracregmlogit-type objects. 
#' 
#' @name summary.fracregmlogit
#' @aliases summary.fracregmlogit.pe
#' @aliases summary.fracregmlogit.wtp
#' 
#' @param object an object with class "fracregmlogit", "fracregmlogit.pe", or "fracregmlogit.wtp".
#' @param ... Additional arguments passed to the printCoefmat function.
#' @return Returns the object invisibly.
#' 
#' @details This module provides summary methods for three fracregmlogit objects: \code{fracregmlogit}, \code{fracregmlogit.pe}
#' , and \code{fracregmlogit.wtp}. 
#' 
#' For \code{fracregmlogit} objects, the summary prints the number of observations, log pseudo-likelihood,
#' baseline choice, and the coefficient estimates with standard errors, z-statistics, and p-values
#' for each choice equation. 
#' 
#' For \code{fracregmlogit.pe} objects, it displays the marginal or discrete effects 
#' along with their computed standard errors (if Krinsky-Robb sampling was performed) for each choice.
#' 
#' For \code{fracregmlogit.wtp} objects, it provides a table of the aggregated willingness to pay 
#' along with its standard errors and test statistics.
#' @seealso \code{\link{fracregmlogit}}, \code{\link{fracregmlogit.pe}}
#' @examples
#' data("fracreg_spending")
#' X = fracreg_spending[,2:5]
#' y = fracreg_spending[,6:11]
#' 
#' # generate fracregmlogit summary
#' results1 = fracregmlogit(y, X)
#' summary(results1)
#' 
#' # generate marginal effects summary
#' effects1 = fracregmlogit.pe(results1, effect="marginal", se=FALSE)
#' summary(effects1)
#' 
#' @rdname summary.fracregmlogit
#' @exportS3Method summary fracregmlogit
summary.fracregmlogit <- function(object, ...) {
    cat("\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center("Fractional multinomial logit model"), "\n")
    cat(.fracreg.sep(), "\n")
    
    .fracreg.cat.right("Number of observations:", object$count["Obs"])
    .fracreg.cat.right("Baseline choice:", object$baseline)
    if (!is.null(object$likelihood)) {
        .fracreg.cat.right("Log pseudolikelihood:", round(object$likelihood, 4))
    }
    if (object$reps > 0) {
        .fracreg.cat.right("Standard errors:", "bootstrap")
        .fracreg.cat.right("Bootstrap reps:", object$reps)
    } else {
        .fracreg.cat.right("Standard errors:", "robust")
    }
    cat("\n")

    for(i in 1:length(object$estimates)){
        cat(.fracreg.sep(), "\n")
        cat(.fracreg.center(paste("Choice:", names(object$estimates)[i])), "\n")
        cat(.fracreg.sep(), "\n")
        
        res <- object$estimates[[i]]
        colnames(res) <- c("Coefficient", "Std.Err.", "z value", "Pr(>|z|)")
        rownames(res)[rownames(res) == "(Intercept)"] <- "(Intercept)"
        stats::printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
        cat("\n")
    }
    
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
    cat(.fracreg.sep(), "\n")

    invisible(object)
}

#' @export
print.fracregmlogit.pe <- function(x, ...) {
    cat("\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center(x$expl), "\n")
    cat(.fracreg.sep(), "\n")
    cat("\nEffects:\n")
    print.default(x$effects)
    cat("\n")
    invisible(x)
}

#' @export
#' @exportS3Method summary fracregmlogit.pe
summary.fracregmlogit.pe <- function(object, ...) {
    cat("\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center(object$expl), "\n")
    cat(.fracreg.sep(), "\n")
    
    if(!is.null(object$ztable)){
        for(i in 1:length(object$ztable)){
            cat("\n")
            cat(.fracreg.sep(), "\n")
            cat(.fracreg.center(paste("Variable:", names(object$ztable)[i])), "\n")
            cat(.fracreg.sep(), "\n")
            
            res <- object$ztable[[i]]
            colnames(res) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
            stats::printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
        }
        cat("\n")
    } else {
        cat("Effects:\n")
        print.default(object$effects)
    }
    
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
    cat(.fracreg.sep(), "\n")

    invisible(object)
}
############
# generate willingness to pay tables
############

#' @export
summary.fracregmlogit.wtp = function(object, ...) {
  if(!inherits(object, "fracregmlogit.wtp")) stop("Expect an fracregmlogit.wtp object. Wrong object type given.")
  if (is.null(dim(object$wtp))) return(object$wtp)
  if(is.null(colnames(object$wtp)) || colnames(object$wtp)[1]!="estimate") return(object$wtp) # no need to summary.
  
  cat("\nFractional Multinomial Logit Model - Willingness to Pay\n")
  stats::printCoefmat(object$wtp, digits = max(3, getOption("digits") - 2), 
                      signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")
  invisible(object)
}
