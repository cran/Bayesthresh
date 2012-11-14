\name{compVar}
\alias{compVar}

\title{Variance component estimates}

\description{Summaries for variance components for the random effects in the model}

\usage{compVar(object)}

\arguments{
  \item{object}{object class \code{"Bayesthresh"}  }
}

\examples{

### Not run
data(sensory)

Consumer <- factor(sensory$consumer)
Sacarose <- factor(sensory$sacarose)

#### Model 
dex1 <- Bayesthresh(cor ~ (1|Consumer) + Sacarose, burn = 0, 
																			jump = 1, ef.iter = 10, data=sensory) 
compVar(dex1)

}
\keyword{models}
\keyword{variance components}
