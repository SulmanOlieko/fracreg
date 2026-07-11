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

fracreg.table <- function(y,yhat,p,p.var,x.names,type,link,converged,var.type,R2=NA,LL=NULL,method=NULL,var.cluster=NULL,dfc=NULL,or=FALSE,level=0.95)
{
	if(converged==TRUE)
	{
		R2 <- cor(y,yhat)^2

		if(type!="2P" & type!="3P")
		{
			n <- length(y)
			
			Wald <- NA
			p_wald <- NA
			df_wald <- NA
			if(!is.character(p.var)) {
				if(length(p) > 1 && x.names[1] == "(Intercept)") {
					p_idx <- 2:length(p)
					W <- tryCatch(t(p[p_idx]) %*% solve(p.var[p_idx, p_idx]) %*% p[p_idx], error = function(e) NA)
					if (!is.na(W)) {
						Wald <- as.numeric(W)
						df_wald <- length(p_idx)
						p_wald <- 1 - pchisq(Wald, df = df_wald)
					}
				} else if(length(p) > 0 && x.names[1] != "(Intercept)") {
					W <- tryCatch(t(p) %*% solve(p.var) %*% p, error = function(e) NA)
					if (!is.na(W)) {
						Wald <- as.numeric(W)
						df_wald <- length(p)
						p_wald <- 1 - pchisq(Wald, df = df_wald)
					}
				}
			}

			p.sd <- diag(p.var)^0.5
			z.ratio <- p/p.sd
			p.value <- 2*pnorm(abs(z.ratio), lower.tail=FALSE)
			qval <- qnorm((1 + level) / 2)
			ci.low <- p - qval * p.sd
			ci.high <- p + qval * p.sd

			if(or && link == "logit") {
				p <- exp(p)
				p.sd <- p.sd * p
				ci.low <- exp(ci.low)
				ci.high <- exp(ci.high)
			}
			
			conf_str <- paste0("[", round(level * 100), "% Conf.")
			res <- cbind(Coefficient = p, `Std.Err.` = p.sd, `z value` = z.ratio, `ci_low` = ci.low, `ci_high` = ci.high, `Pr(>|z|)` = p.value)
			colnames(res)[4] <- conf_str
			colnames(res)[5] <- "Interval]"

			if(or && link == "logit") {
				colnames(res)[1] <- "Odds Ratio"
			}

			se_type <- var.type
			if(se_type == "standard") se_type <- "EIM"
			colnames(res)[2] <- paste0(tools::toTitleCase(se_type), " Std.Err.")
			rownames(res) <- x.names
			rownames(res)[rownames(res) == "(Intercept)"] <- "(Intercept)"

			cat("\n")
			cat(.fracreg.sep(), "\n")
			if(type=="2Pbin") cat(.fracreg.center(paste("Part 1: Binary", link, "regression")), "\n")
			if(type=="1P") cat(.fracreg.center(paste("Fractional", link, "regression")), "\n")
			if(type=="2Pfrac") cat(.fracreg.center(paste("Part 2: Fractional", link, "regression")), "\n")
			if(type=="3Pbin0") cat(.fracreg.center(paste("Part 1: Binary", link, "regression")), "\n")
			if(type=="3Pbin1") cat(.fracreg.center(paste("Part 2: Binary", link, "regression")), "\n")
			if(type=="3Pfrac") cat(.fracreg.center(paste("Part 3: Fractional", link, "regression")), "\n")
			cat(.fracreg.sep(), "\n")

			if(!is.null(method)) .fracreg.cat.right("Estimator:", method)
			.fracreg.cat.right("Data type:", "Cross-sectional")
			.fracreg.cat.right("Number of observations:", n)
			if(!is.null(var.cluster) && var.type=="cluster") .fracreg.cat.right("Number of clusters:", length(unique(var.cluster)))
			.fracreg.cat.right("Pseudo R-squared:", round(R2, 5))
			if(!is.null(LL)) {
				if(!is.null(method) && method == "ML") {
					.fracreg.cat.right("Log-likelihood:", round(LL, 4))
				} else {
					.fracreg.cat.right("Log pseudolikelihood:", round(LL, 4))
				}
			}
			if(!is.na(Wald)) {
				.fracreg.cat.right(paste0("Wald chi2(", df_wald, "):"), round(Wald, 4))
				.fracreg.cat.right("Prob > chi2:", sprintf("%.4f", p_wald))
			}
			if(var.type!="standard") .fracreg.cat.right("Standard errors:", var.type)
			if(is.null(dfc) || !dfc) {
				.fracreg.cat.right("Small sample correction:", "FALSE")
			} else {
				.fracreg.cat.right("Small sample correction:", "TRUE")
			}
			.fracreg.cat.right("Convergence:", "Successful")
			
			cat(.fracreg.sep(), "\n")
			method_full <- method
			if(!is.null(method)) {
				if(method == "QML") method_full <- "Quasi-Maximum Likelihood"
				else if(method == "ML") method_full <- "Maximum Likelihood"
			}
			cat(.fracreg.center(if(!is.null(method_full)) paste("Final", method_full, "estimates") else "Final estimates"), "\n")
			cat(.fracreg.sep(), "\n")

			stats::printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
			
			cat(.fracreg.sep(), "\n")
			if(type == "1P") cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
			if(type == "1P") cat(.fracreg.sep(), "\n")
		}
		else
		{
			cat("\n")
			cat(.fracreg.sep(), "\n")
			if(type=="2P") cat(.fracreg.center(paste("Two-part fractional regression: binary", link[1], "+ fractional", link[2])) , "\n")
			if(type=="3P") cat(.fracreg.center(paste("Three-part fractional regression: binary", link[1], ", binary", link[2], "+ fractional", link[3])) , "\n")
			cat(.fracreg.sep(), "\n")
			.fracreg.cat.right("Data type:", "Cross-sectional")
			.fracreg.cat.right("Pseudo R-squared:", round(R2, 5))
			.fracreg.cat.right("Convergence:", "Successful")
			cat(.fracreg.sep(), "\n")
			cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
			cat(.fracreg.sep(), "\n")
		}
	}
	else {
		cat(.fracreg.sep(), "\n")
		.fracreg.cat.right("Convergence:", "FAILED OR STOPPED AT BOUNDARY")
		cat(.fracreg.sep(), "\n")
		if(type %in% c("1P", "2P", "3P")) cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
		if(type %in% c("1P", "2P", "3P")) cat(.fracreg.sep(), "\n")
	}
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

	res <- cbind(`Statistic` = S[-1], `p-value` = Sp[-1])
	if(test.which != "GGOFF") {
		rownames(res) <- ver[-1]
	} else {
		rownames(res) <- paste(test.ggoff[-1], ver[-1], sep=" - ")
	}

	cat("\n")
	cat(.fracreg.sep(), "\n")
	cat(.fracreg.center(paste(test.which, "test")), "\n")
	cat(.fracreg.sep(), "\n")
	
	if(test.which != "P")
	{
		cat("H0:", title1, "\n")
		cat(.fracreg.sep(), "\n")
		printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
	}
	else
	{
		cat("H0:", title1, "\n")
		cat("H1:", title2, "\n")
		cat(.fracreg.sep(), "\n")
		stats::printCoefmat(res[1:n.ver, , drop = FALSE], P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
		
		cat(.fracreg.sep(), "\n")
		cat("H0:", title2, "\n")
		cat("H1:", title1, "\n")
		cat(.fracreg.sep(), "\n")
		stats::printCoefmat(res[(n.ver+1):(2*n.ver), , drop = FALSE], P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
	}
	cat(.fracreg.sep(), "\n")
	cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
	cat(.fracreg.sep(), "\n")
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

fracreg.pe.table <- function(PE.p,PE.sd,PE.type,which.x,xvar.names,title,at)
{
	z.ratio <- PE.p/PE.sd
	p.value <- 2*(1-pnorm(abs(z.ratio)))

	res <- cbind(Estimate = PE.p, `Std. Error` = PE.sd, `z value` = z.ratio, `Pr(>|z|)` = p.value)
	rownames(res) <- which.x
	rownames(res)[rownames(res) == "(Intercept)"] <- "(Intercept)"

	cat("\n\n")
	cat(.fracreg.sep(), "\n")
	if(PE.type=="APE") cat(.fracreg.center("Average partial effects"), "\n")
	if(PE.type=="CPE") cat(.fracreg.center("Conditional partial effects"), "\n")
	cat(.fracreg.sep(), "\n")
	
	cat(.fracreg.center(title), "\n")
	cat(.fracreg.sep(), "\n")
	
	stats::printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
	
	cat(.fracreg.sep(), "\n")
	cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
	cat(.fracreg.sep(), "\n")
	
	if(PE.type=="CPE")
	{
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
