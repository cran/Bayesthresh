\name{plot.random.effects}
\alias{plot.random.effects}

\title{Plot random effects of model}

\description{Plot the estimates of random effects with confidence intervals}

\usage{
	\method{plot}{random.effects}(x, by = NULL, level = 0.95, xlab = NULL, ylab = NULL, 
				      xlim = NULL, ylim = NULL, pch=19, 
				      col.points='red', col.seg=gray(0.5), ...)
}

\arguments{
  \item{x}{\code{"random.effects"} object}
   \item{by}{defines the element to be plotted}
   \item{level}{the confidence level required}
   \item{xlab}{a title for the x axis}
   \item{ylab}{a title for the y axis}
   \item{ylim}{the x limits of the plot}
   \item{xlim}{the y limits of the plot}
   \item{pch}{either an integer specifying a symbol or a single character 
	   to be used as the default in plotting points}
   \item{col.points}{a specification for the default plotting color points}
   \item{col.seg}{a specification for the default plotting color segments}
   \item{...}{Arguments to be passed to method, such as \code{\link{plot}}}
}

\examples{

# Not run
data(sensory)

Consumer <- factor(sensory$consumer)
Sacarose <- factor(sensory$sacarose)

#### Model 
# Not run
dex1 <- Bayesthresh(flavor ~ (1|Consumer) + Sacarose, 
																			burn = 10, jump = 2, ef.iter = 100, data=sensory) 

rand <- random.effects(dex1)
plot(rand, by="Consumer")

}
\keyword{plot.random.effects}
\keyword{random.effects}

