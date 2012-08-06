\name{random.effects}
\alias{random.effects}

\title{Extract the random effects of the model}
\description{
  Extract the random effects of the model of the class Bayesthres
}

\usage{random.effects(object)}

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
dex1 <- Bayesthresh(flavor ~ (1|Consumer) + Sacarose, burn = 0, 
																			jump = 1, ef.iter = 10, data=sensory) 
random.effects(dex1)
}
\keyword{models}
\keyword{methods}
