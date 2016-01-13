\encoding{latin1}
\name{robVcov}
\alias{robVcov}
\title{
Robust Variance Calculation
}
\description{
\code{robVcov} calculates the asymptotic variance for Z-estimators.
}
\usage{
robVcov(u, d.u, id=NULL)
}

\arguments{
  \item{u }{
  A n x q matrix of the estimating equations evaluated at the estimated
  model parameters,
  where n is the number of observations and q is the number of estimating equations.
  }
  \item{d.u }{
  The mean of the jacobian of \code{u} evaluated at the solution to the
  estimating equations, with rows corresponding to the estimating
  equations and columns corresponding to the model parameters. The
  number of model parameters is assumed to equal the number of
  estimating equations such that \code{d.u} is a q x q square matrix.
  }
  \item{id }{
  A factor with levels corresponding the clusters in the data. Default
  is \code{NULL} in which case all observations are considered to be independent.
}
}
\details{
 For robust variance estimation, see van der Vaart (2000).

  For clustered data, the rows in \code{u} are added clusterwise
  resulting in a cluster robust estimate of the variance.
}
\value{
	The estimated covariance matrix.
}
\author{
  Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}
\references{
  van der Vaart, A.W. (2000), \emph{Asymptotic Statistics},
  Cambridge University Press, pp. 52--53.
}
\keyword{cluster}

