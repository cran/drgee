\encoding{latin1}
\name{gee}
\alias{gee}
\alias{coef.gee}
\alias{vcov.gee}
\alias{print.gee}
\alias{summary.gee}
\alias{print.summary.gee}
\title{Generalized Estimating Equations}
\description{
  \code{gee} performs estimation of parameters in a restricted mean
  model using standard GEEs with independent working correlation matrix. For
  clustered data, cluster-robust standard errors are calculated. When
  \code{cond=TRUE}, cluster-specific
  intercepts are assumed.
}
\usage{
gee(formula, link = c("identity", "log", "logit"), data, cond = FALSE,
clusterid, rootFinder = findRoots, ...)
}

\arguments{
  \item{formula}{
    An expression or formula representing the expected outcome
    conditional on covariates.
  }

\item{link}{
  A character string naming the link function to use. Has to be
  \code{"identity"}, \code{"log"} or \code{"logit"}. Default is \code{"identity"}.
}

  \item{data}{
   A data frame or environment containing the variables appearing in
   \code{formula}. If missing, the variables are expected to be found
   in the environment of the \code{formula} argument.
 }

\item{cond}{
  A logical value indicating whether cluster-specific intercepts should
  be assumed. Requires a \code{clusterid} argument.
}

\item{clusterid}{
  A cluster-defining variable or a character string naming a
  cluster-defining variable in the \code{data} argument. If it is not
  found in the \code{data} argument, it will be searched for in the
  calling frame. If missing, each observation will be considered to be
  a separate cluster. This argument is required when \code{cond = TRUE}.
  }

\item{rootFinder}{
  A function to solve a system of non linear equations. Default
  is \code{findRoots}.
}

\item{\dots}{
  Further arguments to be passed to the function \code{rootFinder}.
}
}
\details{
  Estimates parameters in a regression model, defined by
  \code{formula}. When \code{cond=FALSE}, the estimated coefficients are
  identical to those obtained with \code{glm}, but since no
  distributional assumptions are made, a robust variance is
  calculated. When \code{cond=TRUE} and \code{link} is \code{"identity"}
  or \code{"log"}, the coefficients will be calculated
  using conditional estimating equations as described in Goetgeluk and
  Vansteelandt (2008). When \code{cond=TRUE} and \code{link="logit"},
  the coefficients will be calculated by conditional logistic regression
  (with robust standard errors).
}
\value{
  \code{gee} return an object of class \code{gee} containing:

  \item{coefficients }{Estimates of the parameters.}
  \item{vcov }{Robust variance of the estimates.}
  \item{call }{The matched call.}
  \item{gee.data }{The \code{drgee.data} object used in the calculations.}
  \item{optim.object }{An estimation object returned from the function specified
    in the \code{rootFinder}, if this function is called.}
  \item{res }{The residuals from the estimating equations.}
  \item{d.res }{The derivative of the residuals from the estimating equations.}
  \item{naive.var }{The non-robust estimate of covariance of the
    estimates, if available. }
  
  The class methods \code{coef} and \code{vcov} can be used to extract
  the estimated parameters and their covariance matrix from a
  \code{gee} object. \code{summary.drgee} produces a summary of the
  calculations.
}
\seealso{
  \code{glm}
}
\author{
  Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}
\keyword{models}
\keyword{regression}
\keyword{cluster}

\references{
  Goetgeluk S., & Vansteelandt S. (2008), Conditional generalized
  estimating equations for the analysis of clustered and longitudinal
  data. \emph{Biometrics}, \bold{64}(3), pp. 772--780.
}
