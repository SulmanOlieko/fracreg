summary_fracreghet_tests_table <- function(test.which,S,Sp,ver,title1,title2)
{
	res <- cbind(`Statistic` = S[-1], `p-value` = Sp[-1])
	rownames(res) <- ver[-1]

    ans <- list(
        coefficients = res,
        test.which = test.which,
        title1 = title1,
        title2 = title2
    )
    return(ans)
}

print_fracreghet_tests_table <- function(x)
{
    res <- x$coefficients
    test.which <- x$test.which
    title1 <- x$title1
    title2 <- x$title2

	cat("\n")
	cat(.fracreg.sep(), "\n")
	cat(.fracreg.center(paste(test.which, "test")), "\n")
	cat(.fracreg.sep(), "\n")
	cat(.fracreg.center(title1), "\n")
	cat(.fracreg.sep(), "\n")
	
	cat("H0:", title2, "\n")
	cat(.fracreg.sep(), "\n")
	stats::printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE, digits = 4, signif.legend = TRUE)
	
	cat(.fracreg.sep(), "\n")
	cat(.fracreg.center(paste("Run Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), "\n")
	cat(.fracreg.sep(), "\n")
}

fracreghet.tests.table <- function(test.which,S,Sp,ver,title1,title2)
{
    ans <- summary_fracreghet_tests_table(test.which,S,Sp,ver,title1,title2)
    print_fracreghet_tests_table(ans)
}

