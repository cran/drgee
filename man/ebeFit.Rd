\encoding{latin1}
\name{ebeFit}
\alias{ebeFit}
\title{
Exposure nuisance model based estimation
}
\description{
Exposure nuisance model based estimation given a \code{drgeeData} object.
}
\usage{
ebeFit(object, rootFinder=findRoots, ...)
}
\arguments{
  \item{object}{
	A \code{drgeeData} object.
}
  \item{rootFinder}{
	A function used to solve a system of nonlinear equations. Default is \code{findRoots}.
}
  \item{\dots}{
  	Further arguments to be passed to \code{drgeeFit} or to the function \code{rootFinder}.
}
}
\details{
  \code{ebeFit} performs exposure  nuisance model based estimation
  (as described in the man page for \code{\link{drgee}}) given a
  \code{drgeeData} object.

  The models are specified by variables \code{y, a, x, yx, ax, z}
  and as character variables \code{olink} and \code{elink}
  which are assumed to be found in the input
  argument \code{object} as a \code{drgeeData} object.

  The \code{rootFinder} argument is only needed when \code{object$olink} is
  \code{"log"}. The \code{rootFinder} argument is not needed when
  \code{object$olink} is \code{"identity"}. Any user supplied
  function to \code{rootFinder} is expected to have the same input
  and output as \code{findRoots}.

  Robust variance is calculated by \code{robVcov}.

  \code{ebeFit} is meant to be called by \code{drgee}.

}
\value{
	A list containing:
	\item{coefficients }{The estimated parameters in the main model.}
	\item{vcov }{The estimated covariance matrix.}
	\item{optim.object }{The optimization object returned from the
	  function \code{rootFinder} function. Is \code{NULL} when
	  \code{object$olink} is \code{"logit"}.}
}
\note{
	When \code{object$olink} is \code{"logit"}, the result is
	equivalent with logistic regression where the outcome and
	exposure have switched places in the model and robust variance
	is used.
}
\author{
Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}
\seealso{
\code{\link{drgeeData}},
\code{\link{drgee}}, \code{\link{findRoots}} and \code{\link{robVcov}}.
}
\keyword{models}
\keyword{regression}
