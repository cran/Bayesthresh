\name{coef.Bayesthresh}
\alias{coef.Bayesthresh}

\title{Coefficients for the fixed effects model}

\description{Extract coefficents the fixed effects for and object of the class Bayesthresh. Coefficients 
				correspond to columns of space spanned by the columns of design matrix for fixed effects}

\usage{\method{coef}{Bayesthresh}(object,...)}

\arguments{
  \item{object}{is an object of class \code{"Bayesthresh"}}
  \item{...}{no usage}
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

coef(dex1)

}
