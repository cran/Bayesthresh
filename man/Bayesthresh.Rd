\name{Bayesthresh}
\alias{Bayesthresh}

\title{Bayesian thresholds mixed-effects models for categorial data}

\description{
  This library fits a linear mixed model for ordinal categorical responses using
		Bayesian inference via Monte Carlo Markov Chains. Default is Nandran & Chen algorithm
		using Gaussian link function and saving just the summaries of the chains.
  Among the options, package allow for two other options of algorithms, for using
		Student's "t" link function and for saving the full chains.
}

\usage{
       Bayesthresh(formula, data, subset, na.action, A=NULL,
																		algor = list(algorithm="NC", link="Gaussian"),
							           Write, priors = list(ru = 10, su = 2, dre = 20, dse = 5), 
                  burn = 100, jump = 10, ef.iter = 4000, model=TRUE)
}
\arguments{
  \item{formula}{a two-sided linear formula object describing the fixed-effects part of the model, 
						with the response on the left of a \code{~} operator and the terms, separated by \code{+} 
						operators, on the right.  The vertical bar character \code{"|"} separates an expression for
						a model matrix and a grouping factor.}
  \item{data}{an optional data frame containing the variables named in \code{formula}.  By default 
						the variables are taken from the environment from which \code{Bayesthres} is called.}
  \item{subset, na.action}{further model specification arguments as in \code{\link[stats]{lm}}; see 
						there for details.}
  \item{A}{Matrix of variance-covariance of random effects.}
  \item{algor}{is a list that contains the name of the algorithm to be used. By default the algorithm 
						is the NC with function link Gaussian}
  \item{Write}{the Write is a function that by default is FALSE. If TRUE, the function saves the 
						iterations of the sampling process in the file output.txt}
  \item{priors}{priors is a list that contains the parameters of the priors used to estimate the
					 variance components of random effects}
  \item{burn, jump, ef.iter}{are of the arguments of iteraction. By default the burn, jump and ef.iter
    (effective iteractions) are 100, 5 e 1000 respectively}
  \item{model}{logical scalar.  If \code{FALSE} the model frame in slot \code{frame} is truncated to zero rows.}
}
\details{
  \describe{
    \item{subset}{an optional expression indicating the subset of the rows of
      \code{data} that should be used in the fit. This can be a logical
      vector, or a numeric vector indicating which observation numbers are
      to be included, or a  character  vector of the row names to be
      included.  All observations are included by default.}
    \item{na.action}{a function that indicates what should happen when the
      data contain \code{NA}s.  The default action (\code{na.fail}) prints
      an error message and terminate if there are any incomplete
      observations.}
    \item{algor}{the are three options algorithms, AC, MC and NC, with link function Gaussian and 
	    t-Student distribution. The object algort, by default is \code{list(algorithm="NC", link="Gaussian")}.}
    \item{Write}{if Write=TRUE, the chain of iteractions is saved in the file output.txt. 
	    The convergence process can be analyzed by the library \code{coda}}
    \item{priors}{The object defines the priors for the variance components of the AC and MC algorithms. 
								For the NC algorithm can be also defined a prior of residual variance. Objects ru (shape parameter)
								and su (scale parameter) are the parameters of inverse gamma for the variance components. The NC
								algorithm allows to change parameters of the residual variance. dre (shape parameter) and dse
								(scale parameter) define the prior of the residual variance. By default, algorithms AC and MC have
								a residual variance equal to 1.}
    }
  }

\examples{

# Not run
data(sensory)

Consumer <- factor(sensory$consumer)
Sacarose <- factor(sensory$sacarose)

#### Model 
# Not run
dex1 <- Bayesthresh(flavor ~ (1|Consumer) + Sacarose, burn = 0, jump = 1, ef.iter = 10, data=sensory) 
summary(dex1)

}
