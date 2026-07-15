fracregpd.table <- function(p,p.var,x.names,x.exogenous,lags,type,link,converged,N.ini,N,NT.ini,NT,J,dfJ,k,var.type,bootstrap,LL=NULL,or=FALSE,level=0.95)
{
	if(converged==TRUE)
	{
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
			res_table <- cbind(Coefficient = p, `Std.Err.` = p.sd, `z value` = z.ratio, `ci_low` = ci.low, `ci_high` = ci.high, `Pr(>|z|)` = p.value)
			colnames(res_table)[4] <- conf_str
			colnames(res_table)[5] <- "Interval]"
			
			if(or && link == "logit") {
				colnames(res_table)[1] <- "Odds Ratio"
			}
			
			se_type <- if(bootstrap) "bootstrap" else var.type
			if(se_type == "standard") se_type <- "EIM"
			colnames(res_table)[2] <- paste0(tools::toTitleCase(se_type), " Std.Err.")
		}
		else
		{
			conf_str <- paste0("[", round(level * 100), "% Conf.")
			res_table <- cbind(Coefficient = p, `Std.Err.` = NA, `z value` = NA, `ci_low` = NA, `ci_high` = NA, `Pr(>|z|)` = NA)
			colnames(res_table)[4] <- conf_str
			colnames(res_table)[5] <- "Interval]"

			if(or && link == "logit") {
				res_table[,"Coefficient"] <- exp(p)
				colnames(res_table)[1] <- "Odds Ratio"
			}
			se_type <- if(bootstrap) "bootstrap" else var.type
			if(se_type == "standard") se_type <- "EIM"
			colnames(res_table)[2] <- paste0(tools::toTitleCase(se_type), " Std.Err.")
		}
		rownames(res_table) <- x.names
		rownames(res_table)[rownames(res_table) == "(Intercept)"] <- "(Intercept)"

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
		model_desc <- switch(type,
			"QMLcre" = "(correlated random effects) ",
			"QMLfe" = "(fixed effects) ",
			"GMMcre" = "(correlated random effects) ",
			"GMMpre" = "(pooled random effects) ",
			"GMMpfe" = "(pooled fixed effects) ",
			"")
		cat(.fracreg.center(paste("Fractional", link, model_desc, "regression")), "\n")
		cat(.fracreg.sep(), "\n")
		.fracreg.cat.right("Data type:", "Panel")
		.fracreg.cat.right("Estimator:", type)
		.fracreg.cat.right("Convergence:", "Successful")
		
		if(N.ini!=N | NT.ini!=NT) .fracreg.cat.right("Number of obs (initial):", NT.ini)
		.fracreg.cat.right("Number of observations:", NT)
		if(N.ini!=N | NT.ini!=NT) .fracreg.cat.right("Number of groups (initial):", N.ini)
		.fracreg.cat.right("Number of groups:", N)
		.fracreg.cat.right("Obs per group:", NT/N)
		if(!is.null(LL)) .fracreg.cat.right("Log pseudolikelihood:", round(LL, 4))
		
		if(type!="QMLcre" & dfJ>0)
		{
			p.value.J <- 1-pchisq(J,dfJ)
			.fracreg.cat.right("J test (p-value):", paste0(round(J,4), " (", round(p.value.J,4), ")"))
		}
		if(!is.na(Wald)) {
			.fracreg.cat.right(paste0("Wald chi2(", df_wald, "):"), round(Wald, 4))
			.fracreg.cat.right("Prob > chi2:", sprintf("%.4f", p_wald))
		}
		
		.fracreg.cat.right("Exogeneity:", x.exogenous)
		.fracreg.cat.right("Use first lag of instruments:", lags)
		
		if(bootstrap==FALSE) {
			if(var.type == "robust") {
				.fracreg.cat.right("Standard errors:", "HC0")
			} else if(var.type == "cluster") {
				.fracreg.cat.right("Standard errors:", "CRVE")
			} else {
				.fracreg.cat.right("Standard errors:", var.type)
			}
		}
		if(bootstrap==TRUE) .fracreg.cat.right("Standard errors:", "bootstrap")
		cat(.fracreg.sep(), "\n")

		type_full <- switch(type,
			"QMLcre" = "(Correlated Random Effects) Quasi-Maximum Likelihood",
			"QMLfe" = "(Fixed Effects) Quasi-Maximum Likelihood",
			"GMMcre" = "(Correlated Random Effects) Generalized Method of Moments",
			"GMMpre" = "(Pooled Random Effects) Generalized Method of Moments",
			"GMMpfe" = "(Pooled Fixed Effects) Generalized Method of Moments",
			type)
		if (type_full == type) {
			if (grepl("^GMM", type)) type_full <- paste("GMM", sub("^GMM", "", type))
			if (grepl("^QML", type)) type_full <- paste("Quasi-Maximum Likelihood", sub("^QML", "", type))
		}
		cat(.fracreg.center(paste("Final", type_full, "estimates")), "\n")
		cat(.fracreg.sep(), "\n")

		if(type!="QMLcre" | x.exogenous==TRUE) {
			printCoefmat(res_table, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
		}
		else
		{
			printCoefmat(res_table[1:k, , drop=FALSE], P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
			cat("\n")
			cat(.fracreg.center("Reduced form:"), "\n")
			cat(.fracreg.sep(), "\n")
			printCoefmat(res_table[-(1:k), , drop=FALSE], P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
		}

		cat(.fracreg.sep(), "\n")
		cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
		cat(.fracreg.sep(), "\n")
	}
	else
	{
		cat("\n")
		cat(.fracreg.sep(), "\n")
		.fracreg.cat.right("Convergence:", "FAILED OR STOPPED AT BOUNDARY")
		cat(.fracreg.sep(), "\n")
		cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
		cat(.fracreg.sep(), "\n")
	}
	cat("\n")
}

