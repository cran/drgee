\encoding{latin1}
\name{drFit}
\alias{drFit}
\title{
Doubly Robust Generalized Estimating Equations
}
\description{
Doubly robust estimation given a \code{drgeeData} object.
}
\usage{
drFit(object, rootFinder=findRoots, ...)
}
\arguments{
  \item{object }{
	A \code{drgeeData} object.
}
  \item{rootFinder }{
	A function used to solve a system of non linear
	equations. Default is \code{findRoots}.
}
  \item{\dots }{
  	Further arguments to be passed to \code{\link{geeFit}} and \code{rootFinder}.
}
}
\details{
  \code{drFit} performs doubly robust estimation (as described in
  the man page for \code{\link{drgee}}) given a \code{drgeeData} object.

  The models are specified by variables \code{y, a, x, yx, ax, v,
    z} and as character variables \code{olink} and \code{elink}
  which are assumed to be found in the input
  argument \code{object} as a \code{drgeeData} object.

  The \code{rootFinder} argument is not needed when
  \code{object$olink} is \code{"identity"}. Any user supplied
  function to the argument \code{rootFinder} is expected to have
  the same input and output as \code{findRoots}.

  Robust variance is calculated by \code{robVcov}.

  \code{drFit} is meant to be called by \code{drgee}.
}
\value{
  A list containing:

	\item{coefficients }{The estimated parameters in the main model.}
	\item{vcov }{The estimated covariance matrix.}
	\item{optim.object }{The optimization object returned from the
	  function \code{rootFinder}, if called. }
}
\author{
Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}
\seealso{
\code{\link{drgee}}, \code{\link{drgeeData}}, \code{\link{findRoots}}
and \code{\link{robVcov}}.
}
\keyword{models}
\keyword{regression}
\keyword{robust}
