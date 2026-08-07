summary_fracreg_tests_table <- function(test.which,S,Sp,ver,title1,title2=NA,n.ver=NA,test.ggoff=NA)
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

    ans <- list(
        coefficients = res,
        test.which = test.which,
        title1 = title1,
        title2 = title2,
        n.ver = n.ver
    )
    return(ans)
}

print_fracreg_tests_table <- function(x)
{
    test.which <- x$test.which
    title1 <- x$title1
    title2 <- x$title2
    n.ver <- x$n.ver
    res <- x$coefficients

	cat("\n")
	cat(.fracreg.sep(), "\n")
	cat(.fracreg.center(paste(test.which, "test")), "\n")
	cat(.fracreg.sep(), "\n")
	
	if(test.which != "P")
	{
		cat("H0:", title1, "\n")
		cat(.fracreg.sep(), "\n")
		stats::printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
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

fracreg.tests.table <- function(test.which,S,Sp,ver,title1,title2=NA,n.ver=NA,test.ggoff=NA)
{
    ans <- summary_fracreg_tests_table(test.which,S,Sp,ver,title1,title2,n.ver,test.ggoff)
    print_fracreg_tests_table(ans)
}

