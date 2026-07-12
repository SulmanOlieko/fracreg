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

