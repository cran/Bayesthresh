\name{random.effects}
\alias{random.effects}

\title{Extract the random effects of the model}
\description{
  Extract the random effects of the model of the class Bayesthres
}

\usage{random.effects(object, HPDinterval=FALSE, prob=0.95)}

\arguments{
  \item{object}{is an object of class \code{"Bayesthres"}}
		\item{HPDinterval}{If HPDinterval=TRUE, return the credibility interval for fixed effects model. 
						By default is FALSE. See package \code{coda}}
		\item{prob}{A numeric scalar in the interval (0,1) giving the target probability content of the intervals. 
						The nominal probability content of the intervals is the multiple of 1/nrow(obj) nearest to prob}
}


\examples{
# Not run
data(sensory)

Consumer <- factor(sensory$consumer)
Sacarose <- factor(sensory$sacarose)

#### Model 
# Not run
dex1 <- Bayesthresh(cor ~ (1|Consumer) + Sacarose, burn = 0, Write=TRUE,
																			jump = 1, ef.iter = 10, data=sensory) 
random.effects(dex1)
random.effects(dex1, HPDinterval=TRUE)
}
\keyword{models}
\keyword{methods}
