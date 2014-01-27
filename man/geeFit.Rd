\encoding{latin1}
\name{geeFit}
\alias{geeFit}
\title{
Generalized Estimating Equation Solver
}
\description{
\code{geeFit} is basically a wrapper around \code{glm.fit} but also
computes the residuals and their derivatives with respect to the estimated coefficients.
}
\usage{
geeFit(y, x, link=c("identity","log","logit"), ...)
}

\arguments{
  \item{y }{
  The outcome observations.
}
  \item{x }{
  The design matrix.
}
  \item{link }{
  A character string naming the link function. Have to be \code{"identity"},
  \code{"log"} or \code{"logit"}.
}
  \item{\dots}{
  Further arguments to be passed to \code{glm.fit}
}
}
\details{
\code{geeFit} estimates coefficients \eqn{\theta}{theta} using an estimating equation
\eqn{W\{R-k^{-1}(\theta^T W)\}=0}{W(R-k^(-1)(theta^T W))=0} where \eqn{k} is the identity,
log or logit link function. \code{glm.fit} is used to estimate
the coefficients. The \code{family} object used in the \code{glm.fit} call depends
on the argument \code{link}. It will be \code{gaussian()}, \code{quasipoisson()} or
\code{binomial()} when \code{link} is \code{"identity"}, \code{"log"} and \code{"logit"}
respectively.

\code{geeFit} is intended to be used as help function in \code{drFit}, \code{ebeFit}
and \code{obeFit}.
}
\value{
  A list with the following elements:
  \item{coefficients }{A named vector of the estimated coefficients. }
  \item{fitted.values }{A vector of the predicted outcomes evaluated with the
  estimated coefficients.}
  \item{res }{A vector of the residuals, i.e. the difference between the
    the observed outcome \code{y} and \code{fitted.values}.}
  \item{d.res }{The derivatives of the residuals with respect to the model
    parameters evaluated with the estimated coefficients. A matrix with
  rows corresponding to observations and columns corresponding the the
  estimated coefficients.}
}
\author{
  Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}
\seealso{
\code{\link{glm.fit}}
}
\keyword{models}
\keyword{regression}
