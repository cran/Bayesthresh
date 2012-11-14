\name{coef.Bayesthresh}
\alias{coef.Bayesthresh}

\title{Coefficients for the fixed effects model}

\description{Extract the fixed effects of model. Coefficients correspond to columns of space spanned by the 
				columns of design matrix for fixed effects}

\usage{
	\method{coef}{Bayesthresh}(object, HPDinterval=FALSE, prob=0.95,...)
}

\arguments{
  \item{object}{is an object of class \code{"Bayesthresh"}}
		\item{HPDinterval}{If HPDinterval=TRUE, return the credibility interval for fixed effects model. 
						By default is FALSE. See package coda}
		\item{prob}{A numeric scalar in the interval (0,1) giving the target probability content of the intervals. 
						The nominal probability content of the intervals is the multiple of 1/nrow(obj) nearest to prob}
  \item{...}{no usage}
}

\examples{

# Not run
data(sensory)

Consumer <- factor(sensory$consumer)
Sacarose <- factor(sensory$sacarose)

#### Model 
# Not run
dex1 <- Bayesthresh(cor ~ (1|Consumer) + Sacarose, burn = 0,
																			jump = 1, ef.iter = 10, data=sensory, Write=TRUE) 

coef(dex1)
coef(dex1, HPDinterval=TRUE)

# Not run
dex2 <- Bayesthresh(cor ~ (1|Consumer) + Sacarose, burn = 0,
																			jump = 1, ef.iter = 10, data=sensory) 
coef(dex2)


}
