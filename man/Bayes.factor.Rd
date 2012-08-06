\name{Bayes.factor}
\alias{Bayes.factor}

\title{Bayes factor of the two models}

\description{Compares two models by evaluating their Bayes factor}

\usage{Bayes.factor(model1, model2)}

\arguments{
  \item{model1}{object of the class model \code{"Bayesthresh"}} 
  \item{model2}{object of the class model \code{"Bayesthresh"}}
}

\details{
				At each step during the Markov chains, the marginal likelihood for a model is evaluated, 
				conditioning on actual values for the parameters in that step. Bayes factor is then estimated 
				by the ratios of the arithmetic means of marginal likelihoods from both models. Details of the 
				implementation can be found in Sorensen and Gianola (2004). For a discussion of the possible
				interpretation of Bayes factors, see Jeffreys(1961)
	}

\references{

	SORENSEN, D.; GIANOLA, D. Likelihood, bayesian and MCMC methods in
	quantitative genetics. United States of America: Springer, 2004. 740 p.

	JEFFREYS, H. Theory of probability. Oxford: Claredon Press, 1961. 470 p.
}

\examples{
data(sensory)

Consumer <- factor(sensory$consumer)
Sacarose <- factor(sensory$sacarose)

# Not run

#### Model 1

# Model with Gaussian link

dex1 <- Bayesthresh(flavor ~ (1|Consumer) + Sacarose, burn = 0, jump = 1,
																		 ef.iter = 10, data=sensory) 
summary(dex1)

#### Model 2

# Model with t-Student link

dex2 <- Bayesthresh(flavor ~ (1|Consumer) + Sacarose, burn = 0, jump = 1, 
																			ef.iter = 10, algor=list(algorithm="NC", link="t"), 
																			data=sensory) 
summary(dex2)

Bayes.factor(dex1,dex2)

}
