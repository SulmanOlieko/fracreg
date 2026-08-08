#' @title Clean Data for Fractional Regression Models
#' @description An internal helper to gracefully drop missing values across an arbitrary number of vectors and matrices, replicating the functionality of na.action = na.omit for multi-array model inputs.
#' @param ... A variable number of vectors, matrices, or data frames.
#' @param na.action A function specifying how to handle missing values, default is \code{stats::na.omit}. If \code{NULL}, no action is taken.
#' @return A named list containing the subsets of the provided arrays without missing values.
fracreg_clean_data <- function(..., na.action = stats::na.omit) {
  if (is.null(na.action)) return(list(...))
  
  args <- list(...)
  non_null <- !sapply(args, is.null)
  
  if (sum(non_null) > 0) {
    # 1. Determine the target number of observations (N)
    # We look for the first matrix, data.frame, or vector that likely represents the dataset.
    N_target <- 0
    for (i in seq_along(args)) {
        if (!is.null(args[[i]])) {
            if (is.matrix(args[[i]]) || is.data.frame(args[[i]])) {
                N_target <- max(N_target, nrow(args[[i]]))
            } else if (is.vector(args[[i]]) && !is.character(args[[i]])) {
                N_target <- max(N_target, length(args[[i]]))
            }
        }
    }
    
    # If no data structure of length > 1 is found, return as-is
    if (N_target <= 1) return(args)
    
    # 2. Extract ONLY arrays/vectors that match the target N exactly
    valid_idx <- sapply(args, function(x) {
        if (is.null(x)) return(FALSE)
        if (is.matrix(x) || is.data.frame(x)) return(nrow(x) == N_target)
        if (is.vector(x)) return(length(x) == N_target)
        return(FALSE)
    })
    
    if (sum(valid_idx) > 0) {
      valid_args <- args[valid_idx]
      valid_rows <- do.call(stats::complete.cases, valid_args)
      
      if (!all(valid_rows)) {
        for (i in seq_along(args)) {
          if (valid_idx[i]) {
            if (is.matrix(args[[i]]) || is.data.frame(args[[i]])) {
              args[[i]] <- args[[i]][valid_rows, , drop = FALSE]
            } else {
              args[[i]] <- args[[i]][valid_rows]
            }
          }
        }
      }
    }
  }
  return(args)
}

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
    ans <- list()
    ans$call <- object$call
    ans$type <- object$type
    ans$link <- object$link
    ans$p <- object$p
    ans$table.info <- object$table.info
    
    if(!is.null(object$resBIN0) && !is.null(object$resBIN0$table.info)) ans$resBIN0 <- do.call(summary_fracreg_table, object$resBIN0$table.info)
    if(!is.null(object$resBIN1) && !is.null(object$resBIN1$table.info)) ans$resBIN1 <- do.call(summary_fracreg_table, object$resBIN1$table.info)
    if(!is.null(object$resBIN) && !is.null(object$resBIN$table.info)) ans$resBIN <- do.call(summary_fracreg_table, object$resBIN$table.info)
    if(!is.null(object$resFRAC) && !is.null(object$resFRAC$table.info)) ans$resFRAC <- do.call(summary_fracreg_table, object$resFRAC$table.info)
    if(!is.null(object$table.info)) ans$main <- do.call(summary_fracreg_table, object$table.info)
    
    if (ans$type %in% c("1P", "2Pbin", "2Pfrac")) {
         if(!is.null(ans$main) && !is.null(ans$main$coefficients)) ans$coefficients <- ans$main$coefficients
    } else {
         coef_list <- list()
         if(!is.null(ans$resBIN0) && !is.null(ans$resBIN0$coefficients)) coef_list$BIN0 <- ans$resBIN0$coefficients
         if(!is.null(ans$resBIN1) && !is.null(ans$resBIN1$coefficients)) coef_list$BIN1 <- ans$resBIN1$coefficients
         if(!is.null(ans$resBIN) && !is.null(ans$resBIN$coefficients)) coef_list$BIN <- ans$resBIN$coefficients
         if(!is.null(ans$resFRAC) && !is.null(ans$resFRAC$coefficients)) coef_list$FRAC <- ans$resFRAC$coefficients
         ans$coefficients <- coef_list
    }
    class(ans) <- "summary.fracreg"
    return(ans)
}

#' @export
print.summary.fracreg <- function(x, ...) {
    if(!is.null(x$resBIN0)) print_fracreg_table(x$resBIN0)
    if(!is.null(x$resBIN1)) print_fracreg_table(x$resBIN1)
    if(!is.null(x$resBIN)) print_fracreg_table(x$resBIN)
    if(!is.null(x$resFRAC)) print_fracreg_table(x$resFRAC)
    if(!is.null(x$main)) print_fracreg_table(x$main)
    invisible(x)
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
        if(x$type == "2P") title <- paste0("Two-part fractional regression: binary ", x$link[1], " + fractional ", x$link[2])
        if(x$type == "3P") title <- paste0("Three-part fractional regression: binary ", x$link[1], ", binary ", x$link[2], " + fractional ", x$link[3])
    }
    cat(paste0("\n", title, "\n"))
    cat(paste0("\nCall:  "))
    if(!is.null(x$call)) {
        print(x$call)
    }
    cat("\nCoefficients:\n")
    if(!is.null(x$resBIN0)) {
        cat(paste0("Binary ", x$link[1], " regression (Pr(y > 0)):\n"))
        print.default(x$resBIN0$p)
        cat("\n")
    }
    if(!is.null(x$resBIN1)) {
        cat(paste0("Binary ", x$link[2], " regression (Pr(y = 1 | y > 0)):\n"))
        print.default(x$resBIN1$p)
        cat("\n")
    }
    if(!is.null(x$resBIN)) {
        cond_str <- if (!is.null(x$inf) && x$inf == 1) "Pr(y = 1)" else "Pr(y > 0)"
        cat(paste0("Binary ", x$link[1], " regression (", cond_str, "):\n"))
        print.default(x$resBIN$p)
        cat("\n")
    }
    if(!is.null(x$resFRAC)) {
        link_idx <- if (x$type == "3P") 3 else 2
        cat(paste0("Fractional ", x$link[link_idx], " regression (Fractional Component):\n"))
        print.default(x$resFRAC$p)
        cat("\n")
    }
    if(is.null(x$resBIN0) && is.null(x$resBIN1) && is.null(x$resBIN) && is.null(x$resFRAC) && !is.null(x$p)) {
        print.default(x$p)
        cat("\n")
    }
    
    if (!is.null(x$table.info) && !is.null(x$table.info$y)) {
        n <- length(x$table.info$y)
        k <- length(x$p)
        if(is.null(k)) k <- 0
        df.residual <- max(0, n - k)
        cat(paste0("Degrees of Freedom: ", n, " Total; ", df.residual, " Residual\n"))
    }
    if (!is.null(x$table.info) && !is.null(x$table.info$LL)) {
        cat(paste0("Log-likelihood/Log-pseudolikelihood: ", round(x$table.info$LL, 4), "\n"))
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.fracregpd <- function(object, ...) {
    ans <- list()
    ans$call <- object$call
    ans$type <- object$type
    ans$link <- object$link
    ans$p <- object$p
    ans$table.info <- object$table.info
    if(!is.null(object$table.info)) {
        ans$main <- do.call(summary_fracregpd_table, object$table.info)
        ans$coefficients <- ans$main$coefficients
    }
    class(ans) <- "summary.fracregpd"
    return(ans)
}

#' @export
print.summary.fracregpd <- function(x, ...) {
    if(!is.null(x$main)) print_fracregpd_table(x$main)
    invisible(x)
}

#' @export
print.fracregpd <- function(x, ...) {
    model_desc <- "panel data"
    if (!is.null(x$type) && x$type %in% c("QMLcre", "GMMcre")) model_desc <- "correlated random effects"
    link_str <- if(!is.null(x$link)) x$link else ""
    title <- paste("Fractional", link_str, model_desc, "regression")
    cat(paste0("\n", title, "\n"))
    cat(paste0("\nCall:  "))
    if(!is.null(x$call)) {
        print(x$call)
    }
    cat("\nCoefficients:\n")
    if(!is.null(x$p)) {
        print.default(x$p)
        cat("\n")
    }
    if (!is.null(x$table.info) && !is.null(x$table.info$N)) {
        n <- x$table.info$NT
        k <- length(x$p)
        df.residual <- max(0, n - k)
        cat(paste0("Degrees of Freedom: ", n, " Total; ", df.residual, " Residual\n"))
    }
    if (!is.null(x$table.info) && !is.null(x$table.info$LL)) {
        cat(paste0("Log-likelihood/Log-pseudolikelihood: ", round(x$table.info$LL, 4), "\n"))
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreghet <- function(object, ...) {
    ans <- list()
    ans$call <- object$call
    ans$type <- object$type
    ans$link <- object$link
    ans$p <- object$p
    ans$table.info <- object$table.info
    if(!is.null(object$table.info)) {
        ans$main <- do.call(summary_fracreghet_table, object$table.info)
        ans$coefficients <- ans$main$coefficients
    }
    class(ans) <- "summary.fracreghet"
    return(ans)
}

#' @export
print.summary.fracreghet <- function(x, ...) {
    if(!is.null(x$main)) print_fracreghet_table(x$main)
    invisible(x)
}

#' @export
print.fracreghet <- function(x, ...) {
    link_str <- if(!is.null(x$link)) x$link else ""
    title <- paste("Fractional", link_str, "regression with heteroscedasticity/endogeneity")
    cat(paste0("\n", title, "\n"))
    cat(paste0("\nCall:  "))
    if(!is.null(x$call)) {
        print(x$call)
    }
    cat("\nCoefficients:\n")
    if(!is.null(x$p)) {
        print.default(x$p)
        cat("\n")
    }
    if (!is.null(x$table.info) && !is.null(x$table.info$N)) {
        n <- x$table.info$N
        k <- length(x$p)
        df.residual <- max(0, n - k)
        cat(paste0("Degrees of Freedom: ", n, " Total; ", df.residual, " Residual\n"))
    }
    if (!is.null(x$table.info) && !is.null(x$table.info$LL)) {
        cat(paste0("Log-likelihood/Log-pseudolikelihood: ", round(x$table.info$LL, 4), "\n"))
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreg.pe <- function(object, ...) {
    ans <- list()
    ans$link <- object$link
    ans$inf <- object$inf
    if(!is.null(object$table.info)) {
        ans$main <- do.call(summary_fracreg_pe_table, object$table.info)
        ans$coefficients <- list(overall = ans$main$coefficients)
    } else if(!is.null(object$ape) || !is.null(object$cpe)) {
        ans$ape <- if(!is.null(object$ape$table.info)) do.call(summary_fracreg_pe_table, object$ape$table.info) else NULL
        ans$cpe <- if(!is.null(object$cpe$table.info)) do.call(summary_fracreg_pe_table, object$cpe$table.info) else NULL
        
        coef_list <- list()
        if(!is.null(ans$ape)) coef_list$ape_overall <- ans$ape$coefficients
        if(!is.null(ans$cpe)) coef_list$cpe_overall <- ans$cpe$coefficients
        ans$coefficients <- coef_list
    }

    # Extract summaries for components if they exist
    for (comp in c("resBIN0", "resBIN1", "resBIN", "resFRAC")) {
        if (!is.null(object[[comp]])) {
            comp_sum <- summary(object[[comp]])
            ans[[comp]] <- comp_sum
            
            clean_comp <- sub("^res", "", comp)
            if (!is.null(comp_sum$coefficients)) {
                if (is.list(comp_sum$coefficients)) {
                    for (nm in names(comp_sum$coefficients)) {
                        ans$coefficients[[paste0(clean_comp, "_", nm)]] <- comp_sum$coefficients[[nm]]
                    }
                } else {
                    ans$coefficients[[clean_comp]] <- comp_sum$coefficients
                }
            }
        }
    }

    # Reorder coefficients so overall blocks are at the bottom for tidy()
    overall_names <- intersect(names(ans$coefficients), c("overall", "ape_overall", "cpe_overall"))
    other_names <- setdiff(names(ans$coefficients), overall_names)
    ans$coefficients <- ans$coefficients[c(other_names, overall_names)]

    # Unlist if only one 'overall' block so standard logic remains backward compatible
    if (length(ans$coefficients) == 1 && names(ans$coefficients)[1] %in% c("overall", "ape_overall", "cpe_overall")) {
        ans$coefficients <- ans$coefficients[[1]]
    }

    class(ans) <- "summary.fracreg.pe"
    return(ans)
}

#' @export
print.summary.fracreg.pe <- function(x, subcomponent = FALSE, ...) {
    
    has_subcomponents <- !is.null(x$resBIN) || !is.null(x$resBIN0) || !is.null(x$resBIN1) || !is.null(x$resFRAC)
    
    # Print the banner at the top if it's the main call AND we have subcomponents
    if(!subcomponent && has_subcomponents) {
        pe_obj <- if(!is.null(x$main)) x$main else if(!is.null(x$ape)) x$ape else if(!is.null(x$cpe)) x$cpe
        if(!is.null(pe_obj)) {
            cat("\n\n")
            cat(.fracreg.sep(), "\n")
            if(pe_obj$PE.type=="APE") cat(.fracreg.center("Average partial effects"), "\n")
            if(pe_obj$PE.type=="CPE") cat(.fracreg.center("Conditional partial effects"), "\n")
            cat(.fracreg.sep(), "\n")
        }
    }

    if(!is.null(x$resBIN0)) {
        cat(paste0("\nBinary ", x$link[1], " regression (Pr(y > 0)):\n"))
        print.summary.fracreg.pe(x$resBIN0, subcomponent = TRUE)
    }
    if(!is.null(x$resBIN1)) {
        cat(paste0("\nBinary ", x$link[2], " regression (Pr(y = 1 | y > 0)):\n"))
        print.summary.fracreg.pe(x$resBIN1, subcomponent = TRUE)
    }
    if(!is.null(x$resBIN)) {
        cond_str <- if (!is.null(x$inf) && x$inf == 1) "Pr(y = 1)" else "Pr(y > 0)"
        cat(paste0("\nBinary ", x$link[1], " regression (", cond_str, "):\n"))
        print.summary.fracreg.pe(x$resBIN, subcomponent = TRUE)
    }
    if(!is.null(x$resFRAC)) {
        link_idx <- if (is.null(x$resBIN0)) 2 else 3
        cat(paste0("\nFractional ", x$link[link_idx], " regression (Fractional Component):\n"))
        print.summary.fracreg.pe(x$resFRAC, subcomponent = TRUE)
    }
    if(!is.null(x$main) || !is.null(x$ape) || !is.null(x$cpe)) {
        # Ensure only necessary spacing
    }
    if(!is.null(x$main)) print_fracreg_pe_table(x$main, subcomponent, print.banner = !has_subcomponents)
    if(!is.null(x$ape)) print_fracreg_pe_table(x$ape, subcomponent, print.banner = !has_subcomponents)
    if(!is.null(x$cpe)) print_fracreg_pe_table(x$cpe, subcomponent, print.banner = !has_subcomponents)
    
    invisible(x)
}

#' @export
print.fracreg.pe <- function(x, ...) {
    print_pe_object <- function(obj) {
        if(!is.null(obj$table.info) && !is.null(obj$table.info$title)) {
            title <- obj$table.info$title
        } else {
            title <- "Fractional response regression"
        }
        
        pe_type_name <- if(!is.null(obj$table.info) && !is.null(obj$table.info$PE.type)) obj$table.info$PE.type else "Average"
        if(pe_type_name == "APE") pe_type_name <- "Average partial effects"
        if(pe_type_name == "CPE") pe_type_name <- "Conditional partial effects"
        
        cat(paste0("\n", title, "\n"))
        
        cat(paste0("\n", pe_type_name, ":\n"))
        if(!is.null(obj$resBIN0)) {
            cat(paste0("Binary ", obj$link[1], " regression (Pr(y > 0)):\n"))
            print.default(if(!is.null(obj$resBIN0$table.info)) obj$resBIN0$table.info$PE.p else obj$resBIN0$PE.p)
            cat("\n")
        }
        if(!is.null(obj$resBIN1)) {
            cat(paste0("Binary ", obj$link[2], " regression (Pr(y = 1 | y > 0)):\n"))
            print.default(if(!is.null(obj$resBIN1$table.info)) obj$resBIN1$table.info$PE.p else obj$resBIN1$PE.p)
            cat("\n")
        }
        if(!is.null(obj$resBIN)) {
            cond_str <- if (!is.null(obj$inf) && obj$inf == 1) "Pr(y = 1)" else "Pr(y > 0)"
            cat(paste0("Binary ", obj$link[1], " regression (", cond_str, "):\n"))
            print.default(if(!is.null(obj$resBIN$table.info)) obj$resBIN$table.info$PE.p else obj$resBIN$PE.p)
            cat("\n")
        }
        if(!is.null(obj$resFRAC)) {
            link_idx <- if (is.null(obj$resBIN0)) 2 else 3
            cat(paste0("Fractional ", obj$link[link_idx], " regression (Fractional Component):\n"))
            print.default(if(!is.null(obj$resFRAC$table.info)) obj$resFRAC$table.info$PE.p else obj$resFRAC$PE.p)
            cat("\n")
        }
        if(is.null(obj$resBIN0) && is.null(obj$resBIN1) && is.null(obj$resBIN) && is.null(obj$resFRAC)) {
            if(!is.null(obj$table.info)) print.default(obj$table.info$PE.p)
            else print.default(obj$PE.p)
            cat("\n")
        } else {
            cat("Overall:\n")
            if(!is.null(obj$table.info)) print.default(obj$table.info$PE.p)
            else print.default(obj$PE.p)
            cat("\n")
        }
    }
    
    if(!is.null(x$ape)) print_pe_object(x$ape)
    if(!is.null(x$cpe)) print_pe_object(x$cpe)
    if(is.null(x$ape) && is.null(x$cpe)) print_pe_object(x)
    
    invisible(x)
}


#' @export
summary.fracreghet.pe <- function(object, ...) {
    ans <- list()
    if(!is.null(object$table.info)) {
        ans$main <- do.call(summary_fracreghet_pe_table, object$table.info)
        ans$coefficients <- ans$main$coefficients
    } else if(!is.null(object$ape) || !is.null(object$cpe)) {
        ans$ape <- if(!is.null(object$ape$table.info)) do.call(summary_fracreghet_pe_table, object$ape$table.info) else NULL
        ans$cpe <- if(!is.null(object$cpe$table.info)) do.call(summary_fracreghet_pe_table, object$cpe$table.info) else NULL
        coef_list <- list()
        if(!is.null(ans$ape)) coef_list$ape <- ans$ape$coefficients
        if(!is.null(ans$cpe)) coef_list$cpe <- ans$cpe$coefficients
        ans$coefficients <- coef_list
    }
    class(ans) <- "summary.fracreghet.pe"
    return(ans)
}

#' @export
print.summary.fracreghet.pe <- function(x, ...) {
    if(!is.null(x$main)) print_fracreghet_pe_table(x$main)
    if(!is.null(x$ape)) print_fracreghet_pe_table(x$ape)
    if(!is.null(x$cpe)) print_fracreghet_pe_table(x$cpe)
    invisible(x)
}

#' @export
print.fracreghet.pe <- function(x, ...) {
    print_pe_object <- function(obj) {
        if(!is.null(obj$table.info) && !is.null(obj$table.info$title)) {
            title_vec <- obj$table.info$title
        } else {
            title_vec <- c("", "Fractional regression with heteroscedasticity/endogeneity")
        }
        
        pe_type_name <- if(!is.null(obj$table.info) && !is.null(obj$table.info$PE.type)) obj$table.info$PE.type else "Average"
        if(pe_type_name == "APE") pe_type_name <- "Average partial effects"
        if(pe_type_name == "CPE") pe_type_name <- "Conditional partial effects"
        
        # Format the main title using the components of title_vec
        title_str <- paste(title_vec[-1], collapse = "\n")
        
        # Append title_vec[1] (e.g. "(conditional on observables...)") to the PE type
        if (nzchar(title_vec[1])) {
            pe_type_full <- paste(pe_type_name, title_vec[1])
        } else {
            pe_type_full <- pe_type_name
        }
        
        cat(paste0("\n", title_str, "\n"))
        cat(paste0("\n", pe_type_full, ":\n"))
        
        if(!is.null(obj$table.info)) print.default(obj$table.info$PE.p)
        else print.default(obj$PE.p)
        cat("\n")
    }
    
    if(!is.null(x$ape)) print_pe_object(x$ape)
    if(!is.null(x$cpe)) print_pe_object(x$cpe)
    if(is.null(x$ape) && is.null(x$cpe)) print_pe_object(x)
    
    invisible(x)
}

#' @export
summary.fracreg.reset <- function(object, ...) {
    ans <- list()
    ti <- attr(object, "table.info")
    if(!is.null(ti)) {
        ans$main <- do.call(summary_fracreg_tests_table, ti)
        ans$coefficients <- ans$main$coefficients
    }
    class(ans) <- "summary.fracreg.reset"
    return(ans)
}

#' @export
print.summary.fracreg.reset <- function(x, ...) {
    if(!is.null(x$main)) print_fracreg_tests_table(x$main)
    invisible(x)
}

#' @export
print.fracreg.reset <- function(x, ...) {
    ti <- attr(x, "table.info")
    if(!is.null(ti) && !is.null(ti$title1)) {
        cat(paste0("\nRESET Test for ", ti$title1, "\n\n"))
    } else {
        cat("\nRESET Test for Fractional Regression\n\n")
    }
    print.default(c(x))
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreghet.reset <- function(object, ...) {
    ans <- list()
    ti <- attr(object, "table.info")
    if(!is.null(ti)) {
        ans$main <- do.call(summary_fracreghet_tests_table, ti)
        ans$coefficients <- ans$main$coefficients
    }
    class(ans) <- "summary.fracreghet.reset"
    return(ans)
}

#' @export
print.summary.fracreghet.reset <- function(x, ...) {
    if(!is.null(x$main)) print_fracreghet_tests_table(x$main)
    invisible(x)
}

#' @export
print.fracreghet.reset <- function(x, ...) {
    ti <- attr(x, "table.info")
    if(!is.null(ti) && !is.null(ti$title1)) {
        cat(paste0("\nRESET Test for ", ti$title1, " with Heteroskedasticity\n\n"))
    } else {
        cat("\nRESET Test for Fractional Regression with Heteroskedasticity\n\n")
    }
    print.default(c(x))
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreg.ggoff <- function(object, ...) {
    ans <- list()
    ti <- attr(object, "table.info")
    if(!is.null(ti)) {
        ans$main <- do.call(summary_fracreg_tests_table, ti)
        ans$coefficients <- ans$main$coefficients
    }
    class(ans) <- "summary.fracreg.ggoff"
    return(ans)
}

#' @export
print.summary.fracreg.ggoff <- function(x, ...) {
    if(!is.null(x$main)) print_fracreg_tests_table(x$main)
    invisible(x)
}

#' @export
print.fracreg.ggoff <- function(x, ...) {
    ti <- attr(x, "table.info")
    if(!is.null(ti) && !is.null(ti$title1)) {
        cat(paste0("\nGOFF Test for ", ti$title1, "\n\n"))
    } else {
        cat("\nGOFF Test for Fractional Regression\n\n")
    }
    print.default(c(x))
    cat("\n")
    invisible(x)
}

#' @export
summary.fracreg.ptest <- function(object, ...) {
    ans <- list()
    ti <- attr(object, "table.info")
    if(!is.null(ti)) {
        ans$main <- do.call(summary_fracreg_tests_table, ti)
        ans$coefficients <- ans$main$coefficients
    }
    class(ans) <- "summary.fracreg.ptest"
    return(ans)
}

#' @export
print.summary.fracreg.ptest <- function(x, ...) {
    if(!is.null(x$main)) print_fracreg_tests_table(x$main)
    invisible(x)
}

#' @export
print.fracreg.ptest <- function(x, ...) {
    ti <- attr(x, "table.info")
    if(!is.null(ti) && !is.null(ti$title1)) {
        cat(paste0("\nP-Test for ", ti$title1, "\n\n"))
    } else {
        cat("\nP-Test for Fractional Regression\n\n")
    }
    print.default(c(x))
    cat("\n")
    invisible(x)
}

#' @export
print.fracregridge <- function(x, ...) {
    cat("\nFractional Ridge Regression\n")
    cat(paste0("\nCall:  "))
    if(!is.null(x$call)) {
        print(x$call)
    }
    cat("\nCoefficients:\n")
    print(x$coef)
    cat("\n")
    if (!is.null(x$y) && !is.null(x$coef)) {
        n <- nrow(x$y)
        k <- nrow(x$coef)
        df.residual <- max(0, n - k)
        cat(paste0("Degrees of Freedom: ", n, " Total; ", df.residual, " Residual\n"))
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.fracregridge <- function(object, ...) {
    ans <- list()
    ans$table.info <- object$table.info
    ans$stats.info <- object$stats.info
    ans$coefficients <- object$table.info
    class(ans) <- "summary.fracregridge"
    return(ans)
}

#' @export
print.summary.fracregridge <- function(x, ...) {
    cat("\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center("Fractional ridge regression"), "\n")
    cat(.fracreg.sep(), "\n")
    
    .fracreg.cat.right("Data type:", "Cross-sectional")
    .fracreg.cat.right("Convergence:", "Successful")
    .fracreg.cat.right("Standard errors:", "homoskedastic")
    cat(.fracreg.sep(), "\n")
    
    for (name in names(x$table.info)) {
        cat(.fracreg.center(paste("Target Fraction:", sub(".*frac_", "", name))), "\n")
        cat(.fracreg.sep(), "\n")
        if(!is.null(x$stats.info)) {
            stats <- x$stats.info[[name]]
            .fracreg.cat.right("Number of observations:", stats$n_obs)
            .fracreg.cat.right("Pseudo R-squared:", round(stats$R2, 5))
            .fracreg.cat.right("Degrees of freedom:", round(stats$n_obs - stats$df_alpha, 2))
            if (!is.null(stats$W) && !is.na(stats$W)) {
                .fracreg.cat.right(paste0("Wald chi2(", stats$df_W, "):"), round(stats$W, 4))
                .fracreg.cat.right("Prob > chi2:", sprintf("%.4f", stats$p_W))
            }
            cat(.fracreg.sep(), "\n")
        }
        suppressWarnings(stats::printCoefmat(x$table.info[[name]], P.values=TRUE, has.Pvalue=TRUE, digits=4, signif.legend=TRUE))
        cat(.fracreg.sep(), "\n")
    }
    
    cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
    cat(.fracreg.sep(), "\n\n")
    invisible(x)
}

#' @export
summary.fracregridge.pe <- function(object, ...) {
    ans <- list()
    ans$table.info <- object$table.info
    ans$coefficients <- object$table.info
    class(ans) <- "summary.fracregridge.pe"
    return(ans)
}

#' @export
print.summary.fracregridge.pe <- function(x, ...) {
    cat("\n\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center("Average partial effects"), "\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center("Fractional ridge regression"), "\n")
    cat(.fracreg.sep(), "\n")
    
    cat("\nNote: Fractional ridge regression is a linear model without a link function.\n")
    cat("Therefore, the partial effects are mathematically identical to the coefficients themselves.\n\n")
    
    for (name in names(x$table.info)) {
        cat(.fracreg.center(paste("Target Fraction:", sub(".*frac_", "", name))), "\n")
        cat(.fracreg.sep(), "\n")
        res <- x$table.info[[name]]
        colnames(res)[1] <- "dy/dx"
        suppressWarnings(stats::printCoefmat(res, P.values=TRUE, has.Pvalue=TRUE, digits=4, signif.legend=TRUE))
        cat(.fracreg.sep(), "\n")
    }
    cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
    cat(.fracreg.sep(), "\n\n")
    invisible(x)
}

#' @export
print.fracregridge.pe <- function(x, ...) {
    cat("\nFractional ridge regression\n")
    cat("\nAverage partial effects:\n")
    
    cat("\nNote: Fractional ridge regression is a linear model without a link function.\n")
    cat("Therefore, the partial effects are mathematically identical to the coefficients themselves.\n\n")
    
    for (name in names(x$table.info)) {
        cat(paste0("Target Fraction: ", sub(".*frac_", "", name), "\n"))
        res <- x$table.info[[name]]
        colnames(res)[1] <- "dy/dx"
        print.default(res[,"dy/dx", drop=FALSE])
        cat("\n")
    }
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
#' @param x an object of class "fracregmlogit", "fracregmlogit.pe", or "fracregmlogit.wtp" (used by print methods).
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
#' @export
print.fracregmlogit <- function(x, ...) {
    cat("\nFractional Multinomial Logit Regression\n")
    cat(paste0("\nCall:  "))
    if(!is.null(x$call)) {
        print(x$call)
    }
    cat("\nCoefficients:\n")
    print.default(t(x$coefficient[-1, , drop = FALSE]))
    cat("\n")
    if (!is.null(x$count)) {
        n <- unname(x$count["Obs"])
        k <- unname(x$count["Explanatories"]) + 1
        df.residual <- max(0, n - k)
        cat(paste0("Degrees of Freedom: ", n, " Total; ", df.residual, " Residual\n"))
    }
    if (!is.null(x$likelihood)) {
        cat(paste0("Log-likelihood/Log-pseudolikelihood: ", round(x$likelihood, 4), "\n"))
    }
    cat("\n")
    invisible(x)
}

#' @rdname summary.fracregmlogit
#' @exportS3Method summary fracregmlogit
summary.fracregmlogit <- function(object, ...) {
    ans <- list()
    ans$count <- object$count
    ans$likelihood <- object$likelihood
    ans$pseudo_R2 <- object$pseudo_R2
    ans$baseline <- object$baseline
    ans$reps <- object$reps
    ans$estimates <- object$estimates
    ans$wald <- object$wald
    
    # Store coefficients for broom::tidy
    coef_list <- list()
    for(i in 1:length(object$estimates)){
        res <- object$estimates[[i]]
        se_type <- if (object$reps > 0) "Bootstrap" else "Robust"
        colnames(res) <- c("Coefficient", paste(se_type, "Std.Err."), "z value", "Pr(>|z|)")
        rownames(res)[rownames(res) == "(Intercept)"] <- "(Intercept)"
        coef_list[[names(object$estimates)[i]]] <- res
    }
    ans$coefficients <- coef_list
    class(ans) <- "summary.fracregmlogit"
    return(ans)
}

#' @export
print.summary.fracregmlogit <- function(x, ...) {
    cat("\n")
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center("Fractional multinomial logit model"), "\n")
    cat(.fracreg.sep(), "\n")
    
    .fracreg.cat.right("Data type:", "Cross-sectional")
    .fracreg.cat.right("Convergence:", "Successful")
    .fracreg.cat.right("Number of observations:", x$count["Obs"])
    if (!is.null(x$likelihood)) {
        .fracreg.cat.right("Log pseudolikelihood:", round(x$likelihood, 4))
    }
    if (!is.null(x$pseudo_R2)) {
        .fracreg.cat.right("Pseudo R-squared:", round(x$pseudo_R2, 5))
    }
    .fracreg.cat.right("Baseline choice:", x$baseline)
    
    if (x$reps > 0) {
        .fracreg.cat.right("Standard errors:", "bootstrap")
        .fracreg.cat.right("Bootstrap reps:", x$reps)
    } else {
        .fracreg.cat.right("Standard errors:", "HC0")
    }
    cat("\n")

    for(i in 1:length(x$coefficients)){
        cat(.fracreg.sep(), "\n")
        cat(.fracreg.center(paste("Choice:", names(x$coefficients)[i])), "\n")
        cat(.fracreg.sep(), "\n")
        
        if (!is.null(x$wald[[i]]) && !is.na(x$wald[[i]]$W)) {
            .fracreg.cat.right(paste0("Wald chi2(", x$wald[[i]]$df, "):"), round(x$wald[[i]]$W, 4))
            .fracreg.cat.right("Prob > chi2:", sprintf("%.4f", x$wald[[i]]$p))
            cat(.fracreg.sep(), "\n")
        }
        
        res <- x$coefficients[[i]]
        stats::printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
        cat("\n")
    }
    
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
    cat(.fracreg.sep(), "\n")

    invisible(x)
}

#' @export
print.fracregmlogit.pe <- function(x, ...) {
    cat("\nFractional multinomial logit regression\n")
    cat("\nPartial effects:\n")
    cat("\nNote:", x$expl, "\n\n")
    print.default(x$effects)
    cat("\n")
    invisible(x)
}

#' @export
#' @exportS3Method summary fracregmlogit.pe
summary.fracregmlogit.pe <- function(object, ...) {
    ans <- list()
    ans$expl <- object$expl
    ans$ztable <- object$ztable
    ans$effects <- object$effects
    
    coef_list <- list()
    if(!is.null(object$ztable)){
        for(i in 1:length(object$ztable)){
            res <- object$ztable[[i]]
            colnames(res) <- c("dy/dx", "Std. Error", "z value", "Pr(>|z|)")
            coef_list[[names(object$ztable)[i]]] <- res
        }
    } else {
        for(i in 1:ncol(object$effects)){
            res <- matrix(object$effects[,i], ncol=1)
            rownames(res) <- rownames(object$effects)
            colnames(res) <- c("dy/dx")
            coef_list[[colnames(object$effects)[i]]] <- res
        }
    }
    ans$coefficients <- coef_list
    class(ans) <- "summary.fracregmlogit.pe"
    return(ans)
}

#' @export
print.summary.fracregmlogit.pe <- function(x, ...) {
    cat("\n\n")
    cat(.fracreg.sep(), "\n")
    if(grepl("average across observations", x$expl)) {
        cat(.fracreg.center("Average partial effects"), "\n")
    } else {
        cat(.fracreg.center("Conditional partial effects"), "\n")
    }
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center("Fractional multinomial logit regression"), "\n")
    cat(.fracreg.sep(), "\n")
    cat("\nNote:", x$expl, "\n")
    
    if(!is.null(x$ztable)){
        for(i in 1:length(x$coefficients)){
            cat("\n")
            cat(.fracreg.sep(), "\n")
            cat(.fracreg.center(paste("Choice:", names(x$coefficients)[i])), "\n")
            cat(.fracreg.sep(), "\n")
            
            res <- x$coefficients[[i]]
            stats::printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
        }
        cat("\n")
    } else {
        cat("Effects:\n")
        print.default(x$effects)
    }
    
    cat(.fracreg.sep(), "\n")
    cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
    cat(.fracreg.sep(), "\n")

    invisible(x)
}
############
# generate willingness to pay tables
############

#' @export
summary.fracregmlogit.wtp = function(object, ...) {
  if(!inherits(object, "fracregmlogit.wtp")) stop("Expect an fracregmlogit.wtp object. Wrong object type given.")
  ans <- list()
  ans$wtp <- object$wtp
  ans$coefficients <- object$wtp
  class(ans) <- "summary.fracregmlogit.wtp"
  return(ans)
}

#' @export
print.fracregmlogit.wtp = function(x, ...) {
    cat("\nFractional multinomial logit regression\n")
    cat("\nWillingness to Pay:\n\n")
    print.default(x$wtp)
    cat("\n")
    invisible(x)
}

#' @export
print.summary.fracregmlogit.wtp = function(x, ...) {
  if (is.null(dim(x$wtp))) { print(x$wtp); return(invisible(x)) }
  if(is.null(colnames(x$wtp)) || colnames(x$wtp)[1]!="Coefficient") { print(x$wtp); return(invisible(x)) } 
  
  cat("\n")
  cat(.fracreg.sep(), "\n")
  cat(.fracreg.center("Willingness to Pay"), "\n")
  cat(.fracreg.sep(), "\n")
  cat(.fracreg.center("Fractional multinomial logit regression"), "\n")
  cat(.fracreg.sep(), "\n\n")
  
  cat("Note: Krinsky-Robb standard error calculated\n")
  stats::printCoefmat(x$wtp, digits = max(3, getOption("digits") - 2), 
                      signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE)
  
  cat("\n")
  cat(.fracreg.sep(), "\n")
  cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
  cat(.fracreg.sep(), "\n\n")
  
  invisible(x)
}

#' @exportS3Method broom::tidy
tidy.fracreg <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
    res <- summary(x)$coefficients
    
    if(is.list(res) && !is.data.frame(res)) {
        out <- do.call(rbind, lapply(names(res), function(nm) {
            d <- as.data.frame(res[[nm]])
            comp <- sub("^frac_", "", nm)
            colnames(d)[colnames(d) %in% c("Coefficient", "dy/dx", "Odds Ratio", "Estimate")] <- "estimate"
            colnames(d)[grepl("Std\\. ?Err", colnames(d), ignore.case = TRUE)] <- "std.error"
            colnames(d)[colnames(d) == "z value"] <- "statistic"
            colnames(d)[colnames(d) == "Pr(>|z|)"] <- "p.value"
            
            data.frame(
                component = ifelse(comp == "BIN", "BIN0", comp),
                term = rownames(d),
                estimate = if (!is.null(d$estimate)) as.numeric(d$estimate) else rep(NA_real_, nrow(d)),
                std.error = if (!is.null(d$std.error)) as.numeric(d$std.error) else rep(NA_real_, nrow(d)),
                statistic = if (!is.null(d$statistic)) as.numeric(d$statistic) else rep(NA_real_, nrow(d)),
                p.value = if (!is.null(d$p.value)) as.numeric(d$p.value) else rep(NA_real_, nrow(d)),
                stringsAsFactors = FALSE
            )
        }))
    } else {
        d <- as.data.frame(res)
        colnames(d)[colnames(d) %in% c("Coefficient", "dy/dx", "Odds Ratio", "Estimate")] <- "estimate"
        colnames(d)[grepl("Std\\. ?Err", colnames(d), ignore.case = TRUE)] <- "std.error"
        colnames(d)[colnames(d) == "z value"] <- "statistic"
        colnames(d)[colnames(d) == "Pr(>|z|)"] <- "p.value"
        
        out <- data.frame(
            term = rownames(d),
            estimate = if (!is.null(d$estimate)) as.numeric(d$estimate) else rep(NA_real_, nrow(d)),
            std.error = if (!is.null(d$std.error)) as.numeric(d$std.error) else rep(NA_real_, nrow(d)),
            statistic = if (!is.null(d$statistic)) as.numeric(d$statistic) else rep(NA_real_, nrow(d)),
            p.value = if (!is.null(d$p.value)) as.numeric(d$p.value) else rep(NA_real_, nrow(d)),
            stringsAsFactors = FALSE
        )
    }
    
    rownames(out) <- NULL
    
    if (conf.int) {
        alpha <- 1 - conf.level
        z_crit <- qnorm(1 - alpha / 2)
        out$conf.low <- out$estimate - z_crit * out$std.error
        out$conf.high <- out$estimate + z_crit * out$std.error
    }
    
    if (requireNamespace("tibble", quietly = TRUE)) {
        out <- tibble::as_tibble(out)
    }
    
    out
}

tidy_fracreg_test <- function(x, ...) {
    res <- summary(x)$coefficients
    d <- as.data.frame(res)
    
    out <- data.frame(
        term = rownames(d),
        statistic = if (!is.null(d[["Statistic"]])) as.numeric(d[["Statistic"]]) else rep(NA_real_, nrow(d)),
        p.value = if (!is.null(d[["p-value"]])) as.numeric(d[["p-value"]]) else rep(NA_real_, nrow(d)),
        stringsAsFactors = FALSE
    )
    
    rownames(out) <- NULL
    
    if (requireNamespace("tibble", quietly = TRUE)) {
        out <- tibble::as_tibble(out)
    }
    
    out
}


#' @exportS3Method broom::glance
glance.fracreg <- function(x, ...) {
    out <- data.frame(
        nobs = nobs(x),
        logLik = as.numeric(logLik(x)),
        stringsAsFactors = FALSE
    )
    if (requireNamespace("tibble", quietly = TRUE)) {
        out <- tibble::as_tibble(out)
    }
    out
}

#' @exportS3Method stats::formula
formula.fracreg <- function(x, ...) {
    # Extract original call
    cl <- x$call
    
    # Try to extract variable names
    x_names <- cl$xnames
    
    # If xnames is not stored, attempt to evaluate x to get its colnames
    if (is.null(x_names)) {
        x_val <- tryCatch(eval(cl$x, envir = parent.frame()), error = function(e) NULL)
        if (!is.null(x_val) && inherits(x_val, "fracreg")) {
            x_val <- tryCatch(eval(cl$x, envir = globalenv()), error = function(e) NULL)
        } else if (is.null(x_val)) {
            x_val <- tryCatch(eval(cl$x, envir = globalenv()), error = function(e) NULL)
        }
        if (is.matrix(x_val) || is.data.frame(x_val)) {
            x_names <- colnames(x_val)
        }
    }
    
    # Construct formula string
    if (is.null(x_names)) {
        x_vars <- "1"
    } else {
        # Remove intercept if it's there
        x_vars <- x_names[x_names != "(Intercept)"]
        if (length(x_vars) == 0) x_vars <- "1"
    }
    
    f <- stats::as.formula(paste("y ~", paste(x_vars, collapse = " + ")))
    # Set the formula environment to parent.frame() to avoid capturing
    # the local execution environment where `x` is the model object.
    environment(f) <- parent.frame()
    return(f)
}

#' @exportS3Method stats::model.frame
model.frame.fracreg <- function(formula, ...) {
    object <- formula
    if (!is.null(object$call$formula)) {
        return(stats::model.frame(stats::formula(object), data = eval(object$call$data, envir = parent.frame())))
    } else if (!is.null(object$x) && !is.null(object$y)) {
        y_name <- "y"
        if (is.character(object$call$y) || is.symbol(object$call$y)) y_name <- as.character(object$call$y)
        
        df <- data.frame(y = as.vector(object$y))
        colnames(df)[1] <- y_name
        df <- cbind(df, as.data.frame(object$x))
        return(df)
    } else {
        env <- parent.frame()
        x_val <- tryCatch(eval(object$call$x, envir = env), error = function(e) NULL)
        if (!is.null(x_val) && inherits(x_val, "fracreg")) {
            x_val <- tryCatch(eval(object$call$x, envir = globalenv()), error = function(e) NULL)
        } else if (is.null(x_val)) {
            x_val <- tryCatch(eval(object$call$x, envir = globalenv()), error = function(e) NULL)
        }
        
        y_val <- tryCatch(eval(object$call$y, envir = env), error = function(e) NULL)
        if (!is.null(y_val) && inherits(y_val, "fracreg")) {
            y_val <- tryCatch(eval(object$call$y, envir = globalenv()), error = function(e) NULL)
        } else if (is.null(y_val)) {
            y_val <- tryCatch(eval(object$call$y, envir = globalenv()), error = function(e) NULL)
        }
        if (is.null(x_val) || is.null(y_val)) {
            stop("Unable to construct model.frame: original 'x' or 'y' variables are not accessible from the current environment.")
        }
        
        y_name <- "y"
        if (is.character(object$call$y) || is.symbol(object$call$y)) y_name <- as.character(object$call$y)
        else if (!is.null(colnames(y_val))) y_name <- colnames(y_val)[1]
        
        df <- data.frame(y = as.vector(y_val))
        colnames(df)[1] <- y_name
        
        if (is.matrix(x_val) || is.data.frame(x_val)) {
             df <- cbind(df, as.data.frame(x_val))
        } else {
             df$x <- as.vector(x_val)
             if (is.character(object$call$x) || is.symbol(object$call$x)) colnames(df)[2] <- as.character(object$call$x)
        }
        
        f <- stats::formula(object)
        attr(df, "terms") <- stats::terms(f)
        class(df) <- "data.frame"
        return(df)
    }
}

#' @exportS3Method stats::model.matrix
model.matrix.fracreg <- function(object, ...) {
    env <- parent.frame()
    
    # Evaluate variables in the caller's environment, falling back to globalenv()
    x_val <- tryCatch(eval(object$call$x, envir = env), error = function(e) NULL)
    if (!is.null(x_val) && inherits(x_val, "fracreg")) {
        x_val <- tryCatch(eval(object$call$x, envir = globalenv()), error = function(e) NULL)
    } else if (is.null(x_val)) {
        x_val <- tryCatch(eval(object$call$x, envir = globalenv()), error = function(e) NULL)
    }
    
    if (is.null(x_val)) {
        stop("Unable to construct model.matrix: original 'x' variable is not accessible from the current environment.")
    }
    
    mat <- as.matrix(x_val)
    
    # Ensure intercept is prepended if needed
    coef_names <- names(coef(object))
    if (is.null(coef_names) && object$type %in% c("2P", "3P")) {
        coef_names <- unlist(lapply(summary(object)$coefficients, rownames))
    }
    
    if ("(Intercept)" %in% coef_names) {
        mat <- cbind("(Intercept)" = 1, mat)
    }
    
    return(mat)
}

#' @exportS3Method stats::family
family.fracreg <- function(object, ...) {
    link_name <- object$call$linkfrac
    if (is.null(link_name)) link_name <- "logit"
    else {
        link_name <- tryCatch(eval(link_name, envir = parent.frame()), error = function(e) "logit")
        if (!is.character(link_name)) {
            link_name <- as.character(object$call$linkfrac)
        }
    }
    
    if (link_name %in% c("logit", "probit", "cauchit", "cloglog")) {
        return(stats::quasibinomial(link = link_name))
    }
    return(stats::gaussian())
}

#' @exportS3Method stats::family
family.fracregpd <- family.fracreg
#' @exportS3Method stats::family
family.fracreghet <- family.fracreg
#' @exportS3Method stats::family
family.fracregridge <- family.fracreg
#' @exportS3Method stats::family
family.fracregmlogit <- family.fracreg

#' @exportS3Method stats::terms
terms.fracreg <- function(x, ...) {
    if (!is.null(x$call$formula)) {
        return(stats::terms(stats::formula(x), data = eval(x$call$data, envir = parent.frame())))
    } else {
        mf <- stats::model.frame(x)
        return(attr(mf, "terms"))
    }
}

#' @exportS3Method stats::terms
terms.fracreghet <- terms.fracreg
#' @exportS3Method stats::terms
terms.fracregpd <- terms.fracreg
#' @exportS3Method stats::terms
terms.fracregridge <- terms.fracreg
#' @exportS3Method stats::terms
terms.fracregmlogit <- terms.fracreg

#' @exportS3Method broom::tidy
tidy.fracreghet <- tidy.fracreg
#' @exportS3Method broom::tidy
tidy.fracregpd <- tidy.fracreg
#' @exportS3Method broom::tidy
tidy.fracregridge <- tidy.fracreg
#' @exportS3Method broom::tidy
tidy.fracregmlogit <- tidy.fracreg

#' @exportS3Method broom::tidy
tidy.fracreg.pe <- tidy.fracreg
#' @exportS3Method broom::tidy
tidy.fracreghet.pe <- tidy.fracreg
#' @exportS3Method broom::tidy
tidy.fracregpd.pe <- tidy.fracreg
#' @exportS3Method broom::tidy
tidy.fracregridge.pe <- tidy.fracreg
#' @exportS3Method broom::tidy
tidy.fracregmlogit.pe <- tidy.fracreg

#' @exportS3Method broom::tidy
tidy.fracregmlogit.wtp <- tidy.fracreg

#' @exportS3Method broom::tidy
tidy.fracreg.reset <- tidy_fracreg_test
#' @exportS3Method broom::tidy
tidy.fracreghet.reset <- tidy_fracreg_test
#' @exportS3Method broom::tidy
tidy.fracreg.ggoff <- tidy_fracreg_test
#' @exportS3Method broom::tidy
tidy.fracreg.ptest <- tidy_fracreg_test

#' @exportS3Method broom::glance
glance.fracreghet <- glance.fracreg
#' @exportS3Method broom::glance
glance.fracregpd <- glance.fracreg
#' @exportS3Method broom::glance
glance.fracregridge <- glance.fracreg
#' @exportS3Method broom::glance
glance.fracregmlogit <- glance.fracreg

#' @exportS3Method stats::formula
formula.fracreghet <- formula.fracreg
#' @exportS3Method stats::formula
formula.fracregpd <- formula.fracreg
#' @exportS3Method stats::formula
formula.fracregridge <- formula.fracreg
#' @exportS3Method stats::formula
formula.fracregmlogit <- formula.fracreg

#' @exportS3Method stats::model.frame
model.frame.fracreghet <- model.frame.fracreg
#' @exportS3Method stats::model.frame
model.frame.fracregpd <- model.frame.fracreg
#' @exportS3Method stats::model.frame
model.frame.fracregridge <- model.frame.fracreg
#' @exportS3Method stats::model.frame
model.frame.fracregmlogit <- model.frame.fracreg

#' @exportS3Method stats::model.matrix
model.matrix.fracreghet <- model.matrix.fracreg
#' @exportS3Method stats::model.matrix
model.matrix.fracregpd <- model.matrix.fracreg
#' @exportS3Method stats::model.matrix
model.matrix.fracregridge <- model.matrix.fracreg
#' @exportS3Method stats::model.matrix
model.matrix.fracregmlogit <- model.matrix.fracreg

#' @exportS3Method gtsummary::tbl_regression
tbl_regression.fracreg <- function(x, ...) {
  if (!requireNamespace("gtsummary", quietly = TRUE)) {
    stop("Package 'gtsummary' is required for tbl_regression.")
  }
  
  NextMethod()
}
