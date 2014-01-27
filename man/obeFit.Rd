\encoding{latin1}
\name{obeFit}
\alias{obeFit}
\title{
Outcome Nuisance Model Based estimation
}
\description{
	Outcome nuisance based estimation given a \code{drgeeData} object.
}
\usage{
obeFit(object, ...)
}

\arguments{
  \item{object}{
	A \code{drgeeData} object.
}
  \item{\dots}{
  	Further arguments to be passed to \code{geeFit}
}
}
\details{
  \code{obeFit} performs outcome nuisance model based estimation
  (as described in the man page for \code{\link{drgee}})) given a
  \code{drgeeData} object.

  The model is specified by variables \code{y, a, x, ax, v}
  and as a character variable \code{olink} which are assumed to be
  found in the input argument \code{object} as a \code{drgeeData} object.

  Robust variance is calculated by \code{robVcov}.

  \code{obeFit} is meant to be called by \code{drgee}.

}
\value{
	A list containing:
	\item{coefficients }{The estimated parameters in the main model.}
	\item{vcov }{The estimated covariance matrix.}
}
\note{
	When \code{object$olink} is \code{"logit"}, the result is
	asymptotically equivalent with logistic regression.
}

\author{
Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}
\seealso{
\code{\link{drgeeData}}, \code{\link{drgee}} and \code{\link{robVcov}}.
}
\keyword{models}
\keyword{regression}
