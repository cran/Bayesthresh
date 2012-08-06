Bayes.factor <- function(model1, model2)
{
	if(!inherits(model1, "Bayesthresh"))
		stop("Model1 does not belong a class Bayesthresh")
	if(!inherits(model2, "Bayesthresh"))
		stop("Model2 does not belong a class Bayesthresh")
	veros1 <- 1/((1/model1$inter$ef.iter)*(as.numeric(model1$lik)))
	Postmean1 <- veros1[1]
	veros2 <- 1/((1/model2$inter$ef.iter)*(as.numeric(model2$lik)))
	Postmean2 <- veros2[1]
	bf <- Postmean1/Postmean2
	result <- data.frame(bf)
	rownames(result) <- c("model1/model2")
	colnames(result) <- c("Bayes factor")
	cat("\nBayes factor for comparison two models\n")
	cat("\n \n")
	cat("Model 1:", deparse(model1$formula), fill=TRUE )
	cat("Model 2:", deparse(model2$formula), "\n")
	cat("\n \n")
	print(result)
}

