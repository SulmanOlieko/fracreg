fracreghet.tests.table <- function(test.which,S,Sp,ver,title1,title2)
{
	res <- cbind(`Statistic` = S[-1], `p-value` = Sp[-1])
	rownames(res) <- ver[-1]

	cat("\n")
	cat(.fracreg.sep(), "\n")
	cat(.fracreg.center(paste(test.which, "test")), "\n")
	cat(.fracreg.sep(), "\n")
	cat(.fracreg.center(title1), "\n")
	cat(.fracreg.sep(), "\n")
	
	cat("H0:", title2, "\n")
	cat(.fracreg.sep(), "\n")
	printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
	
	cat(.fracreg.sep(), "\n")
	cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
	cat(.fracreg.sep(), "\n")
}

