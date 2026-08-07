summary_fracreghet_pe_table <- function(PE.p,PE.sd,PE.type,which.x,xvar.names,title,adjust,at)
{
	z.ratio <- PE.p/PE.sd
	p.value <- 2*(1-pnorm(abs(z.ratio)))

	res <- cbind(`dy/dx` = PE.p, `Std. Error` = PE.sd, `z value` = z.ratio, `Pr(>|z|)` = p.value)
	rownames(res) <- which.x
	rownames(res)[rownames(res) == "(Intercept)"] <- "(Intercept)"

    ans <- list(
        coefficients = res,
        PE.type = PE.type,
        title = title,
        adjust = adjust,
        at = at,
        which.x = which.x,
        xvar.names = xvar.names,
        PE.sd = PE.sd
    )
    return(ans)
}

print_fracreghet_pe_table <- function(x)
{
    res <- x$coefficients
    PE.type <- x$PE.type
    title <- x$title
    adjust <- x$adjust
    at <- x$at
    xvar.names <- x$xvar.names
    PE.sd <- x$PE.sd

	cat("\n")
	cat(.fracreg.sep(), "\n")
	if(PE.type=="APE") cat(.fracreg.center("Average partial effects"), "\n")
	if(PE.type=="CPE") cat(.fracreg.center("Conditional partial effects"), "\n")
	cat(.fracreg.sep(), "\n")
	
	if(nzchar(title[1])) cat(.fracreg.center(title[1]), "\n")
	cat(.fracreg.center(title[2]), "\n")
	cat(.fracreg.center(title[3]), "\n")
	
	if(adjust!=0) cat(.fracreg.center(title[4]), "\n")
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
	cat("\nNote: Partial effects for the reduced form are omitted because it is a linear model where they strictly equal the estimated coefficients.\n")
}

fracreghet.pe.table <- function(PE.p,PE.sd,PE.type,which.x,xvar.names,title,adjust,at)
{
    ans <- summary_fracreghet_pe_table(PE.p,PE.sd,PE.type,which.x,xvar.names,title,adjust,at)
    print_fracreghet_pe_table(ans)
}

