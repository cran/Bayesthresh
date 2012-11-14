\name{MCMCsample}
\alias{MCMCsample}

\title{MCMC sample}
\description{Returns the chains of the MCMC process}

\usage{MCMCsample(object)}

\arguments{
  \item{object}{is an object of class \code{"Bayesthres"}}
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
MCMCsample(dex1)
}
