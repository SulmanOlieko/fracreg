fracreghet.table <- function(p,p.var,x.names,type,link,converged,N,var.type,adjust,k,J,dfJ,LL=NULL,or=FALSE,level=0.95)
{
	if(converged==TRUE)
	{
		stars <- rep("",length(p))

		if(!is.character(p.var))
		{
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
			
			se_label <- var.type
			if(se_label == "standard") se_label <- "EIM"
			colnames(res)[2] <- paste0(tools::toTitleCase(se_label), " Std.Err.")
		}
		else
		{
			conf_str <- paste0("[", round(level * 100), "% Conf.")
			res <- cbind(Coefficient = p, `Std.Err.` = rep(NA,length(p)), `z value` = rep(NA,length(p)), `ci_low` = rep(NA,length(p)), `ci_high` = rep(NA,length(p)), `Pr(>|z|)` = rep(NA,length(p)))
			colnames(res)[4] <- conf_str
			colnames(res)[5] <- "Interval]"

			if(or && link == "logit") {
				res[,"Coefficient"] <- exp(p)
				colnames(res)[1] <- "Odds Ratio"
			}
			se_label <- var.type
			if(se_label == "standard") se_label <- "EIM"
			colnames(res)[2] <- paste0(tools::toTitleCase(se_label), " Std.Err.")
		}
		rownames(res) <- x.names
		rownames(res)[rownames(res) == "(Intercept)"] <- "(Intercept)"

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

		cat("\n")
		cat(.fracreg.sep(), "\n")
		cat(.fracreg.center(paste("Fractional", link, "regression with heteroscedasticity/endogeneity")), "\n")
		cat(.fracreg.sep(), "\n")
		
		.fracreg.cat.right("Data type:", "Cross-sectional")
		.fracreg.cat.right("Estimator:", type)
		.fracreg.cat.right("Convergence:", "Successful")
		.fracreg.cat.right("Number of observations:", N)
		if(!is.null(LL)) .fracreg.cat.right("Log pseudolikelihood:", round(LL, 4))
		
		if(any(type==c("GMMz","LINz")) & dfJ>0)
		{
			p.value <- 1-pchisq(J,dfJ)
			.fracreg.cat.right("J test (p-value):", paste0(round(J, 4), " (", round(p.value, 4), ")"))
		}
		if(!is.na(Wald)) {
			.fracreg.cat.right(paste0("Wald chi2(", df_wald, "):"), round(Wald, 4))
			.fracreg.cat.right("Prob > chi2:", sprintf("%.4f", p_wald))
		}
		
		if(adjust!=0)
		{
			if(is.numeric(adjust)) .fracreg.cat.right("Adjustment:", paste(adjust,"added to all observations"))
			else .fracreg.cat.right("Adjustment:", "all boundary observations dropped")
		}
		
		if(var.type == "robust") {
			.fracreg.cat.right("Standard errors:", "HC0")
		} else if(var.type == "cluster") {
			.fracreg.cat.right("Standard errors:", "CRVE")
		} else {
			.fracreg.cat.right("Standard errors:", var.type)
		}
		
		cat(.fracreg.sep(), "\n")
		type_full <- type
		type_full <- sub("^QML", "Quasi-Maximum Likelihood ", type_full)
		type_full <- sub("^ML", "Maximum Likelihood ", type_full)
		type_full <- trimws(type_full)
		cat(.fracreg.center(paste("Final", type_full, "estimates")), "\n")
		cat(.fracreg.sep(), "\n")
		
		if(all(type!=c("GMMxv","LINxv","QMLxv"))) printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE, na.print = ".")
		if(any(type==c("GMMxv","LINxv","QMLxv")))
		{
			printCoefmat(res[1:k, , drop = FALSE], P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE, na.print = ".")
			cat("\n")
			cat(.fracreg.center("Reduced form:"), "\n")
			cat(.fracreg.sep(), "\n")
			printCoefmat(res[-(1:k), , drop = FALSE], P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE, na.print = ".")
		}
		cat(.fracreg.sep(), "\n")
		cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
		cat(.fracreg.sep(), "\n")
	}
	else {
		cat(.fracreg.sep(), "\n")
		.fracreg.cat.right("Convergence:", "FAILED")
		cat(.fracreg.sep(), "\n")
		cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
		cat(.fracreg.sep(), "\n")
	}
	cat("\n")
}

