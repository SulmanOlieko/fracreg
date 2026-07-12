fracreg.pe.table <- function(PE.p,PE.sd,PE.type,which.x,xvar.names,title,at)
{
	z.ratio <- PE.p/PE.sd
	p.value <- 2*(1-pnorm(abs(z.ratio)))

	res <- cbind(`dy/dx` = PE.p, `Std. Error` = PE.sd, `z value` = z.ratio, `Pr(>|z|)` = p.value)
	rownames(res) <- which.x
	rownames(res)[rownames(res) == "(Intercept)"] <- "(Intercept)"

	cat("\n\n")
	cat(.fracreg.sep(), "\n")
	if(PE.type=="APE") cat(.fracreg.center("Average partial effects"), "\n")
	if(PE.type=="CPE") cat(.fracreg.center("Conditional partial effects"), "\n")
	cat(.fracreg.sep(), "\n")
	
	cat(.fracreg.center(title), "\n")
	cat(.fracreg.sep(), "\n")
	
	if(!is.na(PE.sd[1])) cat("\nNote: Standard errors computed using the Delta method\n")
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

