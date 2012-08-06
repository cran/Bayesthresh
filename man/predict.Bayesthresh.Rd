\name{predict.Bayesthresh}
\alias{predict.Bayesthresh}

\title{Extract the predict values}

\description{Predict values based on the Bayesthres model object}

\usage{\method{predict}{Bayesthresh}(object, ...)}

\arguments{
  \item{object}{ Object of class \code{"Bayesthresh"} }
  \item{...}{No usage}
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

predict(dex1)

}
\keyword{models}
\keyword{methods}
